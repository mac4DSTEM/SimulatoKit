//
//  Detectors.swift
//  Role: the detectors (abTEM `detectors.py`), applied to the full unshifted pattern P(k) = |FFT₂(ψ)|² of one exit
//        wave (a probe: no 1/N; a plane wave: ψ / N first, as `Waves.diffractionPattern`). Angles in mrad on the
//        `fftfreq` grid with Δα = λ / L · 10³.
//

import Foundation
import SimulatoCore

/// The full, unshifted diffraction intensity of one exit wave, the input every detector reads.
public struct DiffractionIntensity: Sendable {
    public let grid: Grid
    public let beam: ElectronBeam
    /// |FFT₂(ψ)|² in FFT order, `[ix * ny + iy]`.
    public let values: [Double]

    /// An intensity given directly, |FFT₂(ψ)|² in FFT order (for example to test a custom `Detector`).
    /// - Precondition: `values` holds nx × ny values.
    public init(grid: Grid, beam: ElectronBeam, values: [Double]) {
        precondition(values.count == grid.count, "DiffractionIntensity: \(grid.nx) x \(grid.ny) values required")
        self.grid = grid
        self.beam = beam
        self.values = values
    }

    /// The intensity of an exit wave (abTEM's `diffraction_patterns(max_angle="full", fftshift=False)`).
    public init(_ waves: Waves) throws {
        var f = waves.field
        if waves.normalization == .values {
            let n = Double(waves.grid.count)
            for i in f.real.indices { f.real[i] /= n; f.imaginary[i] /= n }
        }
        try FFT2D(nx: waves.grid.nx, ny: waves.grid.ny).forward(&f)
        self.init(grid: waves.grid, beam: waves.beam, values: f.squaredMagnitude())
    }

    /// Angular sampling per axis [mrad per pixel]: (1 / extent) · λ · 10³.
    public var angularSampling: (x: Double, y: Double) {
        (grid.reciprocalSamplingX * beam.wavelength * 1e3, grid.reciprocalSamplingY * beam.wavelength * 1e3)
    }

    /// The angles [mrad] of the pixels of an m-point axis at angular sampling d (NumPy's `fftfreq(m, 1 / d / m)`).
    static func angles(_ m: Int, _ d: Double) -> [Double] { Grid.fftfreq(m, 1 / d / Double(m)) }

    /// The antialias cutoff angles [mrad] (abTEM's `cutoff_angles`: ⌊n'/2⌋ Δα, n' the cutoff gpts with the grid's parity).
    public var cutoffAngles: (x: Double, y: Double) {
        let g = Self.cutoffGpts(grid)
        let a = angularSampling
        return (Double(g.nx / 2) * a.x, Double(g.ny / 2) * a.y)
    }

    /// Cutoff gpts with the grid's own parity (abTEM's `_antialias_cutoff_gpts`).
    static func cutoffGpts(_ grid: Grid) -> (nx: Int, ny: Int) {
        let kcut = 2.0 / 3.0 / max(grid.samplingX, grid.samplingY)
        func size(_ n: Int, _ extent: Double) -> Int {
            var m = Int(floor(roundToDecimals(kcut * extent, 7)))
            if m % 2 != n % 2 { m += 1 }
            return m
        }
        return (size(grid.nx, Double(grid.nx) * grid.samplingX), size(grid.ny, Double(grid.ny) * grid.samplingY))
    }

    static func roundToDecimals(_ x: Double, _ decimals: Int) -> Double {
        let f = pow(10.0, Double(decimals))
        return (x * f).rounded(.toNearestOrEven) / f
    }

    /// The pattern cropped to m_x × m_y (abTEM's `fft_crop`) and `fftshift`ed.
    func cropped(_ mx: Int, _ my: Int) -> [Double] {
        let sx = Waves.cropIndices(n: grid.nx, m: mx), sy = Waves.cropIndices(n: grid.ny, m: my)
        var out = [Double](repeating: 0, count: sx.count * sy.count)
        for (cx, ix) in sx.enumerated() {
            let ox = (cx + sx.count / 2) % sx.count
            for (cy, iy) in sy.enumerated() {
                out[ox * sy.count + (cy + sy.count / 2) % sy.count] = values[ix * grid.ny + iy]
            }
        }
        return out
    }
}

/// A detector: one exit wave's intensity in, values out (row-major in `shape`).
public protocol Detector: Sendable {
    /// The shape of one position's values ([] for one number).
    func shape(for grid: Grid, beam: ElectronBeam) throws -> [Int]
    func detect(_ intensity: DiffractionIntensity) throws -> [Double]
}

/// The annular detector (abTEM's `AnnularDetector`, `integrate_radial`): Σ P over inner² ≤ α² < outer².
public struct AnnularDetector: Detector, Hashable {
    public let inner: Double
    /// `nil`: ⌊min(cutoff angles)⌋.
    public let outer: Double?

    public init(inner: Double = 0, outer: Double? = nil) throws {
        guard inner.isFinite, inner >= 0 else { throw SimulatoError.invalidDetector("inner must be finite and ≥ 0") }
        if let o = outer, !(o.isFinite && o >= inner) {
            throw SimulatoError.invalidDetector("outer must be finite and ≥ inner, got \(o)")
        }
        self.inner = inner
        self.outer = outer
    }

    public func shape(for grid: Grid, beam: ElectronBeam) throws -> [Int] { [] }

    /// The outer angle used on this intensity, refused beyond the pattern's maximum angle (abTEM's check).
    public func outerAngle(_ p: DiffractionIntensity) throws -> Double {
        let a = p.angularSampling
        let maximum = min(Double(p.grid.nx / 2) * a.x, Double(p.grid.ny / 2) * a.y)
        guard let o = outer else { let c = p.cutoffAngles; return floor(min(c.x, c.y)) }
        // abTEM: np.isclose(maximum, outer, atol=1e-5), whose default rtol 1e-5 applies to outer (measurements.py 2933).
        if o > maximum && abs(o - maximum) > 1e-5 + 1e-5 * abs(o) {
            throw SimulatoError.invalidDetector("outer \(o) mrad exceeds the maximum simulated angle \(maximum) mrad")
        }
        return o
    }

    public func detect(_ p: DiffractionIntensity) throws -> [Double] {
        let o = try outerAngle(p)
        guard inner <= o else {
            throw SimulatoError.invalidDetector("inner \(inner) mrad exceeds the outer angle \(o) mrad")
        }
        let a = p.angularSampling
        let ax = DiffractionIntensity.angles(p.grid.nx, a.x), ay = DiffractionIntensity.angles(p.grid.ny, a.y)
        let i2 = inner * inner, o2 = o * o
        var sum = 0.0
        for ix in 0..<p.grid.nx {
            for iy in 0..<p.grid.ny {
                let a2 = ax[ix] * ax[ix] + ay[iy] * ay[iy]
                if a2 >= i2 && a2 < o2 { sum += p.values[ix * p.grid.ny + iy] }
            }
        }
        return [sum]
    }
}

/// The pixelated detector (abTEM's `PixelatedDetector`, `diffraction_patterns(max_angle, parity="same")`): the pattern
/// cropped and `fftshift`ed.
public struct PixelatedDetector: Detector, Hashable {
    public enum MaxAngle: Sendable, Hashable { case cutoff, valid, full }
    public let maxAngle: MaxAngle

    public init(maxAngle: MaxAngle = .valid) { self.maxAngle = maxAngle }

    func size(_ grid: Grid) -> (nx: Int, ny: Int) {
        switch maxAngle {
        case .full: return (grid.nx, grid.ny)
        case .cutoff: return DiffractionIntensity.cutoffGpts(grid)
        case .valid:
            let c = DiffractionIntensity.cutoffGpts(grid)
            func v(_ m: Int, _ n: Int) -> Int {
                var k = Int(floor(DiffractionIntensity.roundToDecimals(Double(m) / 2.0.squareRoot(), 7)))
                if k % 2 != n % 2 { k += 1 }
                return k
            }
            return (v(c.nx, grid.nx), v(c.ny, grid.ny))
        }
    }

    public func shape(for grid: Grid, beam: ElectronBeam) throws -> [Int] { let s = size(grid); return [s.nx, s.ny] }

    public func detect(_ p: DiffractionIntensity) throws -> [Double] {
        let s = size(p.grid)
        return p.cropped(s.nx, s.ny)
    }
}

/// The centre of mass of the pattern (abTEM's `DiffractionPatterns.center_of_mass`) over a pixelated detector's region:
/// Σ P · α along x and along y, not divided by Σ P (as abTEM, whose patterns of a normalised probe sum to about 1),
/// with the pixel at index m/2 of the shifted axis at zero angle. Two values: x then y.
public struct CenterOfMassDetector: Detector, Hashable {
    public enum Units: Sendable, Hashable {
        /// Scattering angle [mrad].
        case milliradians
        /// Spatial frequency [1/Å].
        case inverseAngstroms
    }

    public let region: PixelatedDetector
    public let units: Units

    public init(maxAngle: PixelatedDetector.MaxAngle = .valid, units: Units = .inverseAngstroms) {
        region = PixelatedDetector(maxAngle: maxAngle)
        self.units = units
    }

    public func shape(for grid: Grid, beam: ElectronBeam) throws -> [Int] { [2] }

    public func detect(_ p: DiffractionIntensity) throws -> [Double] {
        let s = region.size(p.grid)
        let pattern = p.cropped(s.nx, s.ny)
        let a = p.angularSampling
        let perAngle = units == .milliradians ? 1.0 : 1.0 / (p.beam.wavelength * 1e3)
        var cx = 0.0, cy = 0.0
        for ix in 0..<s.nx {
            let ax = Double(ix - s.nx / 2) * a.x * perAngle
            for iy in 0..<s.ny {
                let value = pattern[ix * s.ny + iy]
                cx += value * ax
                cy += value * Double(iy - s.ny / 2) * a.y * perAngle
            }
        }
        return [cx, cy]
    }
}

/// The segmented detector (abTEM's `SegmentedDetector`, `polar_binning`): n_r × n_φ bins between inner and outer, the
/// azimuth measured from `rotation`; summed in float64 (abTEM sums in float32, so the kit's sums are the closer to exact).
public struct SegmentedDetector: Detector, Hashable {
    public let radialBins: Int
    public let azimuthalBins: Int
    public let inner: Double
    public let outer: Double
    /// [rad]
    public let rotation: Double

    public init(radialBins: Int, azimuthalBins: Int, inner: Double, outer: Double, rotation: Double = 0) throws {
        guard radialBins >= 1, azimuthalBins >= 1 else { throw SimulatoError.invalidDetector("bins must be at least 1") }
        guard inner.isFinite, outer.isFinite, inner >= 0, outer > inner, rotation.isFinite else {
            throw SimulatoError.invalidDetector("need 0 ≤ inner < outer, finite, and a finite rotation")
        }
        self.radialBins = radialBins
        self.azimuthalBins = azimuthalBins
        self.inner = inner
        self.outer = outer
        self.rotation = rotation
    }

    public func shape(for grid: Grid, beam: ElectronBeam) throws -> [Int] { [radialBins, azimuthalBins] }

    public func detect(_ p: DiffractionIntensity) throws -> [Double] {
        let a = p.angularSampling
        func m(_ n: Int, _ d: Double) -> Int {
            var k = Int(2 * (outer / d).rounded(.up)) + 1
            if k % 2 != n % 2 { k += 1 }
            return min(k, n)
        }
        let mx = m(p.grid.nx, a.x), my = m(p.grid.ny, a.y)
        // Bins on the cropped grid in FFT order; the cropped pattern in FFT order is the full pattern at the crop's
        // source indices (shifting both, as abTEM does, changes nothing).
        let sx = Waves.cropIndices(n: p.grid.nx, m: mx), sy = Waves.cropIndices(n: p.grid.ny, m: my)
        let ax = DiffractionIntensity.angles(mx, a.x), ay = DiffractionIntensity.angles(my, a.y)
        var out = [Double](repeating: 0, count: radialBins * azimuthalBins)
        let twoPi = 2 * Double.pi
        for cx in 0..<mx {
            for cy in 0..<my {
                let alpha = (ax[cx] * ax[cx] + ay[cy] * ay[cy]).squareRoot()
                guard alpha >= inner && alpha < outer else { continue }
                var phi = (atan2(ay[cy], ax[cx]) - rotation).truncatingRemainder(dividingBy: twoPi)
                if phi < 0 { phi += twoPi }
                let r = Int((Double(radialBins) * (alpha - inner) / (outer - inner)).rounded(.down))
                // α just below outer can round to r = n_r; abTEM's bin label then exceeds n_r n_φ − 1 and the pixel is
                // dropped (`label_to_index`); so is it here (it used to index past the end).
                guard r < radialBins else { continue }
                let f = (Double(azimuthalBins) * (phi / twoPi)).rounded(.down)
                let az = Int(min(max(f, 0), Double(azimuthalBins - 1)))
                out[az + r * azimuthalBins] += p.values[sx[cx] * p.grid.ny + sy[cy]]
            }
        }
        return out
    }
}

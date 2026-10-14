//
//  Waves.swift
//  Role: a wave function on a grid with its beam: the plane wave, the intensity |ψ|² and the diffraction pattern at the
//        antialias cutoff (registration sections 3.7 and 3.10; abTEM `waves.py`).
//

/// A wave function in real space, with the grid, the beam and abTEM's normalisation convention.
public struct Waves: Sendable, Hashable {
    /// How the wave is normalised, which decides the diffraction pattern's scale (abTEM's metadata `normalization`).
    public enum Normalization: Sendable, Hashable {
        /// ψ = 1 per pixel at the entrance (a plane wave): the pattern is computed from ψ / (nx ny), so it sums to ~1.
        case values
        /// Σ_k |Ψ(k)|² = 1 at the entrance (a probe): the pattern is |FFT(ψ)|² with no division.
        case reciprocalSpace
    }

    /// How far out the diffraction pattern reaches.
    public enum MaxAngle: Sendable, Hashable {
        /// Cropped to the antialias cutoff, odd-sized, the zero frequency at the centre (abTEM's default).
        case cutoff
        /// The whole grid, the zero frequency at index (nx / 2, ny / 2).
        case full
    }

    public var field: ComplexField
    public let grid: Grid
    public let beam: ElectronBeam
    public let normalization: Normalization

    public init(field: ComplexField, grid: Grid, beam: ElectronBeam, normalization: Normalization) {
        precondition(field.nx == grid.nx && field.ny == grid.ny, "Waves: the field must match the grid")
        self.field = field
        self.grid = grid
        self.beam = beam
        self.normalization = normalization
    }

    /// The plane wave ψ = 1 at every pixel (abTEM's `PlaneWave(normalize=False)`).
    public static func planeWave(beam: ElectronBeam, grid: Grid) -> Waves {
        Waves(field: ComplexField(nx: grid.nx, ny: grid.ny,
                                  real: [Double](repeating: 1, count: grid.count),
                                  imaginary: [Double](repeating: 0, count: grid.count)),
              grid: grid, beam: beam, normalization: .values)
    }

    /// |ψ|² per pixel, sampled as the grid.
    public func intensity() -> RealField {
        RealField(nx: grid.nx, ny: grid.ny, values: field.squaredMagnitude(),
                  samplingX: grid.samplingX, samplingY: grid.samplingY)
    }

    /// Σ |ψ|² over the grid.
    public var totalIntensity: Double { field.squaredMagnitude().reduce(0, +) }

    /// The diffraction pattern |Ψ(k)|² (abTEM's `diffraction_patterns(max_angle=..., fftshift=True)`), sampled in 1/Å.
    public func diffractionPattern(maxAngle: MaxAngle = .cutoff) throws -> RealField {
        var f = field
        if normalization == .values {
            let n = Double(grid.count)
            for i in f.real.indices {
                f.real[i] /= n
                f.imaginary[i] /= n
            }
        }
        try FFT2D(nx: grid.nx, ny: grid.ny).forward(&f)
        let size: (nx: Int, ny: Int) = maxAngle == .cutoff ? Antialias.cutoffPatternSize(on: grid) : (grid.nx, grid.ny)
        let sx = Self.cropIndices(n: grid.nx, m: size.nx)
        let sy = Self.cropIndices(n: grid.ny, m: size.ny)
        var out = [Double](repeating: 0, count: size.nx * size.ny)
        for (cx, ix) in sx.enumerated() {
            let ox = (cx + size.nx / 2) % size.nx
            for (cy, iy) in sy.enumerated() {
                let oy = (cy + size.ny / 2) % size.ny
                let re = f.real[ix * grid.ny + iy], im = f.imaginary[ix * grid.ny + iy]
                out[ox * size.ny + oy] = re * re + im * im
            }
        }
        return RealField(nx: size.nx, ny: size.ny, values: out,
                         samplingX: grid.reciprocalSamplingX, samplingY: grid.reciprocalSamplingY)
    }

    /// The source index of each kept frequency, in FFT order (abTEM's `_fft_interpolation_masks_1d` for m ≤ n):
    /// the first ⌊m/2⌋ + 1 (m odd) or m/2 (m even), then the last ⌊m/2⌋.
    static func cropIndices(n: Int, m: Int) -> [Int] {
        precondition(m <= n || m == n + 1, "crop larger than the grid")
        if m >= n { return Array(0..<n) }
        let low = m % 2 == 1 ? m / 2 + 1 : m / 2
        let high = m / 2
        return Array(0..<low) + Array((n - high)..<n)
    }
}

//
//  AzimuthalAverage.swift
//  Role: abTEM 1.0.10's `DiffractionPatterns.azimuthal_average` (measurements.py, `_azimuthal_average`) on a centred
//        diffraction pattern: rings at radii 0, Δ, 2Δ, … below a maximum angle, each the mean of the pixels a step or a
//        Gaussian weight selects around that radius. The angles follow abTEM bit for bit where it matters: the axis
//        coordinates are a float32 `linspace` between the pattern limits [mrad], the radius is float32, the ring test is
//        float64, so a pixel on a ring edge falls on the side abTEM puts it. The sums are float64; abTEM sums a float32
//        pattern in float32, so on float32 data the two differ at about 1e-7 relative.
//

import Foundation
import SimulatoCore

/// How `azimuthalAverage` weights the pixels around each ring (abTEM's `weighting_function`).
public enum AzimuthalWeighting: Sendable, Hashable {
    /// The pixels with |α − ring| < width, each weight 1 (abTEM "step").
    case step
    /// exp(−(α − ring)² / (width² / 2)) over every pixel (abTEM "gaussian").
    case gaussian
}

extension RealField {
    /// The azimuthal average of a centred diffraction pattern (abTEM `DiffractionPatterns.azimuthal_average`): for each
    /// ring radius 0, Δ, 2Δ, … below `maxAngle`, Σ P · w / Σ w with w the weighting around that radius, and 0 where
    /// Σ w is 0. The pattern is centred as `Waves.diffractionPattern` and `PixelatedDetector` return it: sampling in 1/Å,
    /// the zero frequency at index n / 2 of each axis (abTEM `fftshift=True`).
    /// - Parameters:
    ///   - beam: the beam that sets the angles (α = k λ 10³ mrad).
    ///   - maxAngle: the last ring lies below it [mrad]; `nil` is abTEM's default, the largest |lowest angle| of the two
    ///     axes (the edge of the axis that reaches further, not the corner).
    ///   - radialSampling: the ring spacing Δ, in units of the finer angular sampling of the two axes (default 1).
    ///   - weighting: `.step` (default) or `.gaussian`.
    ///   - width: the step half-width or the Gaussian width, in units of the finer angular sampling (default 1).
    /// - Returns: one reciprocal-space profile, sampling Δ / (λ 10³) [1/Å].
    /// - Throws: `invalidDetector` unless the radial sampling, the width and a given maximum angle are finite and
    ///   positive (abTEM returns zeros, NaN or an empty profile), or when the rings would number a million or more.
    public func azimuthalAverage(beam: ElectronBeam, maxAngle: Double? = nil, radialSampling: Double = 1,
                                 weighting: AzimuthalWeighting = .step, width: Double = 1) throws
        -> ReciprocalSpaceLineProfiles {
        guard radialSampling.isFinite, radialSampling > 0 else {
            throw SimulatoError.invalidDetector("the radial sampling must be finite and positive, got \(radialSampling)")
        }
        guard width.isFinite, width > 0 else {
            throw SimulatoError.invalidDetector("the width must be finite and positive, got \(width)")
        }
        if let m = maxAngle, !(m.isFinite && m > 0) {
            throw SimulatoError.invalidDetector("the maximum angle must be finite and positive, got \(m)")
        }
        let wavelength = beam.wavelength
        // abTEM `limits`, then `angular_limits` (k × λ × 10³, in that order), then `angular_sampling`.
        func limits(_ n: Int, _ s: Double) -> (lower: Double, upper: Double) {
            let l = n % 2 == 1 ? (Double(-((n - 1) / 2)) * s, Double((n - 1) / 2) * s)
                               : (Double(-(n / 2)) * s, Double(n / 2 - 1) * s)
            return (l.0 * wavelength * 1e3, l.1 * wavelength * 1e3)
        }
        let lx = limits(nx, samplingX), ly = limits(ny, samplingY)
        let minAngular = Swift.min(samplingX * wavelength * 1e3, samplingY * wavelength * 1e3)
        // abTEM: -min(min(angular_limits)), the lower limit of the lexicographically smaller (lower, upper) pair.
        let maximum = maxAngle ?? -Swift.min(lx.lower, ly.lower)
        let step = radialSampling * minAngular
        let w = width * minAngular

        // NumPy `arange(0, maximum, step)`: ceil(maximum / step) rings at i × step.
        let ringsDouble = (maximum / step).rounded(.up)
        guard ringsDouble < 1e6 else {
            throw SimulatoError.invalidDetector("a maximum angle of \(maximum) mrad at \(step) mrad gives \(ringsDouble) rings")
        }
        let rings = Swift.max(Int(ringsDouble), 0)

        // NumPy `linspace(lower, upper, n, dtype=float32)`: i × step + lower in float64, the last point set to upper,
        // then rounded to float32.
        func linspace32(_ l: (lower: Double, upper: Double), _ n: Int) -> [Float] {
            guard n > 1 else { return [Float(l.lower)] }
            let d = (l.upper - l.lower) / Double(n - 1)
            var out = (0..<n).map { Float(Double($0) * d + l.lower) }
            out[n - 1] = Float(l.upper)
            return out
        }
        let ax = linspace32(lx, nx), ay = linspace32(ly, ny)
        var radius = [Double](repeating: 0, count: nx * ny)
        for ix in 0..<nx {
            for iy in 0..<ny {
                let r: Float = (ax[ix] * ax[ix] + ay[iy] * ay[iy]).squareRoot()
                radius[ix * ny + iy] = Double(r)
            }
        }

        var out = [Double](repeating: 0, count: rings)
        let gaussianScale = w * w / 2
        for i in 0..<rings {
            let center = Double(i) * step
            var weight = 0.0, sum = 0.0
            switch weighting {
            case .step:
                for k in radius.indices where abs(radius[k] - center) < w {
                    weight += 1
                    sum += values[k]
                }
            case .gaussian:
                for k in radius.indices {
                    let d = radius[k] - center
                    let m = exp(-(d * d) / gaussianScale)
                    weight += m
                    sum += values[k] * m
                }
            }
            out[i] = weight > 0 ? sum / weight : 0
        }
        return ReciprocalSpaceLineProfiles(count: rings, values: out, sampling: step / (wavelength * 1e3))
    }
}

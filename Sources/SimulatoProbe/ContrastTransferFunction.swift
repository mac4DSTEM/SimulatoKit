//
//  ContrastTransferFunction.swift
//  Role: slice 3 (registration-slice-3.md section 3): abTEM's CTF — the aberration phase exp(−iχ), the spatial and
//        temporal envelopes of the quasi-coherent approximation, the aperture — and `apply_ctf` for HRTEM images; the
//        closed forms Scherzer defocus, point resolution and Nyquist sampling.
//

import Foundation
import SimulatoCore

extension Aberrations {
    /// abTEM's ∂χ/∂α at (α, φ) for wavelength λ (`SpatialEnvelope`), in rad per rad.
    public func dChiDAlpha(alpha a: Double, phi p: Double, wavelength: Double) -> Double {
        let t1 = (C12 * cos(2 * (p - phi12)) + C10) * a
        let t2 = (C23 * cos(3 * (p - phi23)) + C21 * cos(1 * (p - phi21))) * pow(a, 2)
        let t3 = (C34 * cos(4 * (p - phi34)) + C32 * cos(2 * (p - phi32)) + C30) * pow(a, 3)
        let t4 = (C45 * cos(5 * (p - phi45)) + C43 * cos(3 * (p - phi43)) + C41 * cos(1 * (p - phi41))) * pow(a, 4)
        let t5 = (C56 * cos(6 * (p - phi56)) + C54 * cos(4 * (p - phi54)) + C52 * cos(2 * (p - phi52)) + C50) * pow(a, 5)
        return 2 * Double.pi / wavelength * (t1 + t2 + t3 + t4 + t5)
    }

    /// abTEM's ∂χ/∂φ term (the azimuthal derivative divided by α, as abTEM writes it).
    public func dChiDPhi(alpha a: Double, phi p: Double, wavelength: Double) -> Double {
        let t1 = 1 / 2.0 * (2.0 * C12 * sin(2.0 * (p - phi12))) * a
        let t2 = 1 / 3.0 * (3.0 * C23 * sin(3.0 * (p - phi23)) + 1.0 * C21 * sin(1.0 * (p - phi21))) * pow(a, 2)
        let t3 = 1 / 4.0 * (4.0 * C34 * sin(4.0 * (p - phi34)) + 2.0 * C32 * sin(2.0 * (p - phi32))) * pow(a, 3)
        let t4 = 1 / 5.0 * (5.0 * C45 * sin(5.0 * (p - phi45)) + 3.0 * C43 * sin(3.0 * (p - phi43))
                            + 1.0 * C41 * sin(1.0 * (p - phi41))) * pow(a, 4)
        let t5 = (1 / 6.0) * (6.0 * C56 * sin(6.0 * (p - phi56)) + 4.0 * C54 * sin(4.0 * (p - phi54))
                              + 2.0 * C52 * sin(2.0 * (p - phi52))) * pow(a, 5)
        return -2 * Double.pi / wavelength * (t1 + t2 + t3 + t4 + t5)
    }
}

/// The contrast transfer function of an objective lens (abTEM's `CTF`).
public struct ContrastTransferFunction: Sendable, Hashable {
    public var aberrations: Aberrations
    /// The objective aperture's semiangle [mrad]; `.infinity` for none (abTEM's default).
    public var semiangleCutoff: Double
    public var softAperture: Bool
    /// The 1/e half-width of the focal spread Δ [Å] (temporal coherence); 0 for none.
    public var focalSpread: Double
    /// The angular spread β [mrad] (spatial coherence); 0 for none.
    public var angularSpread: Double

    public init(aberrations: Aberrations = Aberrations(), semiangleCutoff: Double = .infinity, softAperture: Bool = true,
                focalSpread: Double = 0, angularSpread: Double = 0) throws {
        guard aberrations.isFinite else { throw SimulatoError.invalidAberration("every coefficient must be finite") }
        guard semiangleCutoff > 0, !semiangleCutoff.isNaN else { throw SimulatoError.invalidAperture(semiangleCutoff) }
        guard focalSpread.isFinite, focalSpread >= 0 else {
            throw SimulatoError.invalidAberration("focal spread must be finite and ≥ 0 Å, got \(focalSpread)")
        }
        guard angularSpread.isFinite, angularSpread >= 0 else {
            throw SimulatoError.invalidAberration("angular spread must be finite and ≥ 0 mrad, got \(angularSpread)")
        }
        self.aberrations = aberrations
        self.semiangleCutoff = semiangleCutoff
        self.softAperture = softAperture
        self.focalSpread = focalSpread
        self.angularSpread = angularSpread
    }

    /// The temporal envelope exp(−((½π/λ) Δ α²)²) at α [rad].
    public static func temporalEnvelope(alpha: Double, wavelength: Double, focalSpread: Double) -> Double {
        let x = 0.5 * Double.pi / wavelength * focalSpread * alpha * alpha
        return exp(-(x * x))
    }

    /// The spatial envelope exp(−sign(β)(β/2)²((∂χ/∂α)² + (∂χ/∂φ)²)) at (α, φ), β in mrad.
    public static func spatialEnvelope(alpha: Double, phi: Double, wavelength: Double, aberrations: Aberrations,
                                       angularSpread: Double) -> Double {
        let b = angularSpread / 1e3
        let dk = aberrations.dChiDAlpha(alpha: alpha, phi: phi, wavelength: wavelength)
        let dp = aberrations.dChiDPhi(alpha: alpha, phi: phi, wavelength: wavelength)
        let sign: Double = b > 0 ? 1 : (b < 0 ? -1 : 0)
        return exp(-sign * (b / 2) * (b / 2) * (dk * dk + dp * dp))
    }

    /// The kernel H on the grid (FFT order): exp(−iχ) × E_s (β ≠ 0) × E_t (Δ ≠ 0) × aperture (finite cutoff).
    public func kernel(on grid: Grid, wavelength: Double) throws -> ComplexField {
        let aperture = semiangleCutoff.isFinite
            ? try Aperture(semiangleCutoff: semiangleCutoff, soft: softAperture).values(on: grid, wavelength: wavelength)
            : nil
        let kx = grid.frequenciesX, ky = grid.frequenciesY
        var h = ComplexField(nx: grid.nx, ny: grid.ny)
        for ix in 0..<grid.nx {
            for iy in 0..<grid.ny {
                let i = ix * grid.ny + iy
                let alpha = (kx[ix] * kx[ix] + ky[iy] * ky[iy]).squareRoot() * wavelength
                let phi = atan2(ky[iy], kx[ix])
                let chi = aberrations.chi(alpha: alpha, phi: phi, wavelength: wavelength)
                var re = cos(-chi), im = sin(-chi)
                var factor = 1.0
                if angularSpread != 0 {
                    factor *= Self.spatialEnvelope(alpha: alpha, phi: phi, wavelength: wavelength, aberrations: aberrations,
                                                   angularSpread: angularSpread)
                }
                if focalSpread != 0 {
                    factor *= Self.temporalEnvelope(alpha: alpha, wavelength: wavelength, focalSpread: focalSpread)
                }
                if let a = aperture { factor *= a[i] }
                re *= factor
                im *= factor
                h.real[i] = re
                h.imaginary[i] = im
            }
        }
        return h
    }

    /// Scherzer defocus [Å]: sign(Cs) √(1.5 |Cs| λ) (abTEM's convention: defocus = −C10).
    public static func scherzerDefocus(Cs: Double, energy: Double) throws -> Double {
        let lambda = try ElectronBeam(energy: energy).wavelength
        let s: Double = Cs > 0 ? 1 : (Cs < 0 ? -1 : 0)
        return s * (3 / 2 * abs(Cs) * lambda).squareRoot()
    }

    /// Point resolution [Å]: (λ³ |Cs| / 6)^{1/4}.
    public static func pointResolution(Cs: Double, energy: Double) throws -> Double {
        let lambda = try ElectronBeam(energy: energy).wavelength
        return pow(pow(lambda, 3) * abs(Cs) / 6, 1.0 / 4)
    }

    /// Nyquist sampling [Å] for a semiangle α_c [mrad]: 1 / (4 α_c / λ · 10⁻³).
    public static func nyquistSampling(semiangleCutoff: Double, energy: Double) throws -> Double {
        // abTEM raises ZeroDivisionError at 0; a semiangle must be finite and positive (slice-3 refuter).
        guard semiangleCutoff.isFinite, semiangleCutoff > 0 else { throw SimulatoError.invalidAperture(semiangleCutoff) }
        let lambda = try ElectronBeam(energy: energy).wavelength
        return 1 / (4 * semiangleCutoff / lambda * 1e-3)
    }
}

extension Waves {
    /// The waves with a CTF applied (abTEM's `apply_ctf`): IFFT(FFT(ψ) · H). The intensity of the result is the HRTEM image.
    public func applying(_ ctf: ContrastTransferFunction) throws -> Waves {
        let h = try ctf.kernel(on: grid, wavelength: beam.wavelength)
        let fft = try FFT2D(nx: grid.nx, ny: grid.ny)
        var f = field
        fft.forward(&f)
        for i in f.real.indices {
            let a = f.real[i], b = f.imaginary[i], c = h.real[i], d = h.imaginary[i]
            f.real[i] = a * c - b * d
            f.imaginary[i] = a * d + b * c
        }
        fft.inverse(&f)
        return Waves(field: f, grid: grid, beam: beam, normalization: normalization)
    }
}

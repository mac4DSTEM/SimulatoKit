//
//  Probe.swift
//  Role: the convergent probe (registration section 3.8; abTEM `Probe._calculate_array`, `soft_aperture`,
//        `fft_shift_kernel`, `_WavesNormalization`): Ψ(k) = shift · aperture · exp(−iχ), normalised so Σ_k |Ψ(k)|² = 1,
//        then ψ = IFFT(Ψ).
//

import Foundation
import SimulatoCore

/// The objective aperture: a disk of semiangle α_c with abTEM's one-pixel soft edge (or a hard edge).
public struct Aperture: Sendable, Hashable {
    /// Semiangle cutoff [mrad]; `.infinity` for no aperture.
    public let semiangleCutoff: Double
    /// abTEM's default soft edge (True).
    public let soft: Bool

    public init(semiangleCutoff: Double, soft: Bool = true) throws {
        guard semiangleCutoff > 0, !semiangleCutoff.isNaN else { throw SimulatoError.invalidAperture(semiangleCutoff) }
        self.semiangleCutoff = semiangleCutoff
        self.soft = soft
    }

    /// The aperture on the grid's reciprocal space for wavelength λ.
    public func values(on grid: Grid, wavelength: Double) -> [Double] {
        if semiangleCutoff == .infinity { return [Double](repeating: 1, count: grid.count) }
        let kx = grid.frequenciesX, ky = grid.frequenciesY
        let cutoff = semiangleCutoff * 1e-3
        let dax = grid.reciprocalSamplingX * 1e3 * wavelength * 1e-3
        let day = grid.reciprocalSamplingY * 1e3 * wavelength * 1e-3
        var out = [Double](repeating: 0, count: grid.count)
        for ix in 0..<grid.nx {
            for iy in 0..<grid.ny {
                let i = ix * grid.ny + iy
                let alpha = (kx[ix] * kx[ix] + ky[iy] * ky[iy]).squareRoot() * wavelength
                if soft {
                    if ix == 0 && iy == 0 { out[i] = 1; continue }
                    let phi = atan2(ky[iy], kx[ix])
                    let a = cos(phi) * dax, b = sin(phi) * day
                    let d = (a * a + b * b).squareRoot()
                    out[i] = min(max((cutoff - alpha) / d + 0.5, 0), 1)
                } else {
                    out[i] = alpha <= cutoff ? 1 : 0
                }
            }
        }
        return out
    }
}

/// A convergent electron probe.
public struct Probe: Sendable, Hashable {
    public let beam: ElectronBeam
    public let grid: Grid
    public let aperture: Aperture
    public let aberrations: Aberrations

    public init(beam: ElectronBeam, grid: Grid, aperture: Aperture, aberrations: Aberrations = Aberrations()) throws {
        guard aberrations.isFinite else { throw SimulatoError.invalidAberration("every coefficient must be finite") }
        self.beam = beam
        self.grid = grid
        self.aperture = aperture
        self.aberrations = aberrations
    }

    /// The probe centred at `position` (x, y) [Å], in real space, Σ_k |Ψ(k)|² = 1 (so Σ_r |ψ|² = 1 / (nx ny)).
    public func build(at position: SIMD2<Double>) throws -> Waves {
        guard position.x.isFinite, position.y.isFinite else { throw SimulatoError.invalidPosition(index: 0) }
        let lambda = beam.wavelength
        let kx = grid.frequenciesX, ky = grid.frequenciesY
        let ap = aperture.values(on: grid, wavelength: lambda)
        let ux = position.x / grid.samplingX, uy = position.y / grid.samplingY
        let px = Grid.fftfreq(grid.nx, 1.0), py = Grid.fftfreq(grid.ny, 1.0)
        var psi = ComplexField(nx: grid.nx, ny: grid.ny)
        for ix in 0..<grid.nx {
            let sx = -2 * Double.pi * px[ix] * ux
            let (cx, snx) = (cos(sx), sin(sx))
            for iy in 0..<grid.ny {
                let i = ix * grid.ny + iy
                let sy = -2 * Double.pi * py[iy] * uy
                let (cy, sny) = (cos(sy), sin(sy))
                // shift = e^{i sx} · e^{i sy}
                var re = cx * cy - snx * sny
                var im = cx * sny + snx * cy
                re *= ap[i]
                im *= ap[i]
                let alpha = (kx[ix] * kx[ix] + ky[iy] * ky[iy]).squareRoot() * lambda
                let phi = atan2(ky[iy], kx[ix])
                let chi = aberrations.chi(alpha: alpha, phi: phi, wavelength: lambda)
                let (c, s) = (cos(-chi), sin(-chi))
                psi.real[i] = re * c - im * s
                psi.imaginary[i] = re * s + im * c
            }
        }
        var norm = 0.0
        for i in psi.real.indices { norm += psi.real[i] * psi.real[i] + psi.imaginary[i] * psi.imaginary[i] }
        let f = norm.squareRoot()
        for i in psi.real.indices {
            psi.real[i] /= f
            psi.imaginary[i] /= f
        }
        try FFT2D(nx: grid.nx, ny: grid.ny).inverse(&psi)
        return Waves(field: psi, grid: grid, beam: beam, normalization: .reciprocalSpace)
    }

    /// Ψ₀ = aperture · exp(−iχ) / ‖·‖ in reciprocal space, unshifted (registration-performance.md 1.6): `build(at:)`
    /// without the shift ramp, whose modulus is 1, so its norm does not depend on the position.
    package func centredSpectrum() -> ComplexField {
        let lambda = beam.wavelength
        let kx = grid.frequenciesX, ky = grid.frequenciesY
        let ap = aperture.values(on: grid, wavelength: lambda)
        var psi = ComplexField(nx: grid.nx, ny: grid.ny)
        var norm = 0.0
        for ix in 0..<grid.nx {
            for iy in 0..<grid.ny {
                let i = ix * grid.ny + iy
                let alpha = (kx[ix] * kx[ix] + ky[iy] * ky[iy]).squareRoot() * lambda
                let phi = atan2(ky[iy], kx[ix])
                let chi = aberrations.chi(alpha: alpha, phi: phi, wavelength: lambda)
                psi.real[i] = ap[i] * cos(-chi)
                psi.imaginary[i] = ap[i] * sin(-chi)
                norm += psi.real[i] * psi.real[i] + psi.imaginary[i] * psi.imaginary[i]
            }
        }
        let f = norm.squareRoot()
        for i in psi.real.indices {
            psi.real[i] /= f
            psi.imaginary[i] /= f
        }
        return psi
    }

    /// The shift to `position` as two ramps, r_x[ix] = exp(−2πi p_x[ix] u_x) (an nx × 1 field) and r_y[iy] likewise
    /// (1 × ny), computed in float64 as `build(at:)` computes them; the shift is r_x[ix] · r_y[iy].
    package func shiftRamps(at position: SIMD2<Double>) throws -> (x: ComplexField, y: ComplexField) {
        guard position.x.isFinite, position.y.isFinite else { throw SimulatoError.invalidPosition(index: 0) }
        let ux = position.x / grid.samplingX, uy = position.y / grid.samplingY
        let px = Grid.fftfreq(grid.nx, 1.0), py = Grid.fftfreq(grid.ny, 1.0)
        var rx = ComplexField(nx: grid.nx, ny: 1), ry = ComplexField(nx: 1, ny: grid.ny)
        for ix in 0..<grid.nx {
            let sx = -2 * Double.pi * px[ix] * ux
            rx.real[ix] = cos(sx)
            rx.imaginary[ix] = sin(sx)
        }
        for iy in 0..<grid.ny {
            let sy = -2 * Double.pi * py[iy] * uy
            ry.real[iy] = cos(sy)
            ry.imaginary[iy] = sin(sy)
        }
        return (rx, ry)
    }
}

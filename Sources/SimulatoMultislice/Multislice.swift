//
//  Multislice.swift
//  Role: the multislice loop (abTEM `multislice.py`): per slice, transmit by the
//        bandlimited t = IFFT(FFT(exp(iσV)) · A), then propagate by P = exp(−iπλΔz k²) · A. One API, two backends:
//        `.cpu` (float64, Accelerate, the reference) and `.metal` (float32, MPSGraph), the tolerance between them
//        in `MetalAccuracy.registered`.
//

import Foundation
import SimulatoCore

/// Which slices' exit waves to return.
public enum ExitPlanes: Sendable, Hashable {
    /// Only the wave after the last slice.
    case last
    /// abTEM's integer `exit_planes = m`: the incident wave, then after slices m, 2m, ..., and the last.
    case every(Int)
    /// After the given zero-based slice indices (ascending), the last slice always included; −1 is the incident wave.
    case afterSlices([Int])

    /// The plane indices (−1 = incident) for n slices.
    public func planes(sliceCount n: Int) throws -> [Int] {
        switch self {
        case .last:
            return [n - 1]
        case .every(let m):
            guard m >= 1 else { throw SimulatoError.invalidExitPlanes("every(\(m)): the step must be at least 1") }
            if m >= n { return [n - 1] }
            var p = Array(stride(from: m - 1, to: n, by: m))
            if p.last != n - 1 { p.append(n - 1) }
            return [-1] + p
        case .afterSlices(let list):
            guard list == list.sorted(), Set(list).count == list.count,
                  list.allSatisfy({ $0 >= -1 && $0 < n }) else {
                throw SimulatoError.invalidExitPlanes("afterSlices must be ascending, distinct, within -1 ..< \(n)")
            }
            return list.last == n - 1 ? list : list + [n - 1]
        }
    }
}

/// The compute backend. `.cpu` is the float64 reference; `.metal` is float32 on the GPU, within the registered tolerance.
public enum Backend: Sendable, Hashable {
    case cpu
    case metal
}

/// The multislice algorithm (Fourier space, propagator order 1, transmit then propagate).
public enum Multislice {
    /// Propagates `incident` through `potential` and returns the waves at the requested exit planes, in plane order.
    /// `options` reach `.metal` (renormalisation, cache); `.cpu` runs one wave on one thread whatever they say.
    public static func run(_ incident: Waves, through potential: ProjectedPotential, exitPlanes: ExitPlanes = .last,
                           backend: Backend = .cpu, options: ExecutionOptions = .standard) throws -> [Waves] {
        guard incident.grid.matches(potential.grid) else { throw SimulatoError.gridMismatch }
        let planes = try exitPlanes.planes(sliceCount: potential.slices.count)
        switch backend {
        case .cpu:
            let transmissions = try transmissionFunctions(potential, beam: incident.beam)
            return try runCPU(incident, potential: potential, transmissions: transmissions, planes: planes)
        case .metal:
            let engine = try MetalEngine.require()
            let prepared = try MetalMultislice.prepare(potential, beam: incident.beam, tilt: incident.tilt, engine: engine)
            return try MetalMultislice.run(incident, prepared: prepared, planes: planes, options: options)
        }
    }

    /// The same for a `Potential`: on `.cpu` built first (`build()`), on `.metal` built slice by slice on the GPU.
    public static func run(_ incident: Waves, through potential: Potential, exitPlanes: ExitPlanes = .last,
                           backend: Backend = .cpu, options: ExecutionOptions = .standard) throws -> [Waves] {
        guard incident.grid.matches(potential.grid) else { throw SimulatoError.gridMismatch }
        switch backend {
        case .cpu:
            return try run(incident, through: potential.build(), exitPlanes: exitPlanes, backend: .cpu, options: options)
        case .metal:
            let planes = try exitPlanes.planes(sliceCount: potential.sliceCount)
            let engine = try MetalEngine.require()
            let prepared = try MetalMultislice.prepare(potential, beam: incident.beam, tilt: incident.tilt, engine: engine)
            return try MetalMultislice.run(incident, prepared: prepared, planes: planes, options: options)
        }
    }

    /// The bandlimited transmission function of each slice, in float64 (both backends start from these).
    public static func transmissionFunctions(_ potential: ProjectedPotential, beam: ElectronBeam) throws -> [ComplexField] {
        let grid = potential.grid
        let fft = try FFT2D(nx: grid.nx, ny: grid.ny)
        let aperture = Antialias.aperture(on: grid)
        let sigma = beam.interactionParameter
        return potential.slices.map { slice in
            var t = ComplexField(nx: grid.nx, ny: grid.ny)
            for i in slice.values.indices {
                let phase = sigma * slice.values[i]
                t.real[i] = cos(phase)
                t.imaginary[i] = sin(phase)
            }
            fft.forward(&t)
            for i in aperture.indices {
                t.real[i] *= aperture[i]
                t.imaginary[i] *= aperture[i]
            }
            fft.inverse(&t)
            return t
        }
    }

    /// The Fresnel propagator over Δz, bandlimited: exp(−iπλΔz k²) · A (abTEM `_fresnel_propagator_array`, order 1).
    /// With a beam `tilt` [mrad] it carries e^{−2πi (kx tan(tx/1000) + ky tan(ty/1000)) Δz} as well (abTEM
    /// `_apply_tilt_to_fresnel_propagator_array`), applied after the aperture.
    public static func propagator(grid: Grid, beam: ElectronBeam, thickness dz: Double,
                                  tilt: SIMD2<Double> = .zero) -> ComplexField {
        let lambda = beam.wavelength
        let kx = grid.frequenciesX, ky = grid.frequenciesY
        let aperture = Antialias.aperture(on: grid)
        var p = ComplexField(nx: grid.nx, ny: grid.ny)
        for ix in 0..<grid.nx {
            for iy in 0..<grid.ny {
                let i = ix * grid.ny + iy
                let k2 = kx[ix] * kx[ix] + ky[iy] * ky[iy]
                let x = -k2 * Double.pi * dz * lambda
                var re = cos(x) * aperture[i], im = sin(x) * aperture[i]
                if tilt != .zero {
                    let t = -2 * Double.pi * dz * (kx[ix] * tan(tilt.x / 1e3) + ky[iy] * tan(tilt.y / 1e3))
                    let (c, s) = (cos(t), sin(t))
                    (re, im) = (re * c - im * s, re * s + im * c)
                }
                p.real[i] = re
                p.imaginary[i] = im
            }
        }
        return p
    }

    package static func runCPU(_ incident: Waves, potential: ProjectedPotential, transmissions: [ComplexField],
                       planes: [Int]) throws -> [Waves] {
        let grid = potential.grid
        let fft = try FFT2D(nx: grid.nx, ny: grid.ny)
        var propagators: [Double: ComplexField] = [:]
        var psi = incident.field
        var out: [Waves] = []
        if planes.first == -1 { out.append(incident) }
        let wanted = Set(planes)
        for (j, slice) in potential.slices.enumerated() {
            let t = transmissions[j]
            multiply(&psi, by: t)
            fft.forward(&psi)
            if propagators[slice.thickness] == nil {
                propagators[slice.thickness] = propagator(grid: grid, beam: incident.beam, thickness: slice.thickness,
                                                          tilt: incident.tilt)
            }
            multiply(&psi, by: propagators[slice.thickness]!)
            fft.inverse(&psi)
            if wanted.contains(j) {
                out.append(Waves(field: psi, grid: grid, beam: incident.beam, normalization: incident.normalization,
                                 tilt: incident.tilt))
            }
        }
        return out
    }

    /// ψ ← ψ · t, (a + ib)(c + id) = (ac − bd) + i(ad + bc).
    static func multiply(_ psi: inout ComplexField, by t: ComplexField) {
        for i in psi.real.indices {
            let a = psi.real[i], b = psi.imaginary[i], c = t.real[i], d = t.imaginary[i]
            psi.real[i] = a * c - b * d
            psi.imaginary[i] = a * d + b * c
        }
    }
}

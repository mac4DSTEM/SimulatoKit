//
//  MultisliceChecks.swift
//  Role: SimulatoMultislice from outside: the Lobato parameterisation, slicing, the projected potential, the multislice
//        loop with its exit planes and backends, and frozen phonons with their generator.
//

import Foundation
import SimulatoCore
import SimulatoMultislice

/// One Al atom at the centre of a 4 × 4 × 4 Å box.
func singleAtomStructure(z: Int = 13) throws -> AtomicStructure {
    try AtomicStructure(box: SIMD3(4, 4, 4), atoms: [.init(atomicNumber: z, position: SIMD3(2, 2, 2))])
}

func checkMultislice(_ r: inout Report) {
    r.section("Multislice: Lobato parameterisation") { r in
        r.check("supports Z = 1, 13, 103; not 0, 104, -1",
                [1, 13, 103].allSatisfy { LobatoParameterization.supports(atomicNumber: $0) }
                && ![0, 104, -1].contains { LobatoParameterization.supports(atomicNumber: $0) })
        let p = try LobatoParameterization.parameters(atomicNumber: 13)
        r.check("parameters(atomicNumber: 13): five a and five b", p.a.count == 5 && p.b.count == 5 && p.b.allSatisfy { $0 > 0 })
        let s = try LobatoParameterization.scaledParameters(atomicNumber: 13)
        // B = 2π / √b; A = π² a / b^1.5 / κ
        var ok = true
        for i in 0..<5 {
            ok = ok && abs(s.B[i] - 2 * Double.pi / p.b[i].squareRoot()) < 1e-12 * s.B[i]
            ok = ok && abs(s.A[i] - Double.pi * Double.pi * p.a[i] / pow(p.b[i], 1.5) / PhysicalConstants.kappa) <= 1e-12 * abs(s.A[i])
        }
        r.check("scaledParameters: B = 2π/√b, A = π² a / b^1.5 / κ", ok)
        let exact = LobatoParameterization()
        let pi: LobatoParameterization.PiConstant = .exact
        r.equal("default piConstant is .exact", exact.piConstant, pi)
        // f(0) = 8π Σ [A/B/B² + A B/B⁴] = 16π Σ A / B³
        var f0 = 0.0
        for i in 0..<5 { f0 += s.A[i] / pow(s.B[i], 3) }
        f0 *= 16 * Double.pi
        r.nearRelative("projectedScatteringFactor(k2: 0) = 16π Σ A/B³", try exact.projectedScatteringFactor(atomicNumber: 13, k2: 0),
                       f0, 1e-12)
        let ks = [0.0, 0.1, 0.5, 2.0]
        let arr = try exact.projectedScatteringFactor(atomicNumber: 13, k2: ks)
        r.check("projectedScatteringFactor: array overload equals the scalar one, f decreases with k²",
                try ks.indices.allSatisfy { arr[$0] == (try exact.projectedScatteringFactor(atomicNumber: 13, k2: ks[$0])) }
                && zip(arr, arr.dropFirst()).allSatisfy { $0 > $1 })
        var f32 = LobatoParameterization(piConstant: .abTEMFloat32)
        let rel = abs(try f32.projectedScatteringFactor(atomicNumber: 13, k2: 0.3) - exact.projectedScatteringFactor(atomicNumber: 13, k2: 0.3))
            / (try exact.projectedScatteringFactor(atomicNumber: 13, k2: 0.3))
        r.check("abTEMFloat32 moves f by about 5e-8 relative (between 1e-9 and 1e-6)", rel > 1e-9 && rel < 1e-6, "got \(rel)")
        f32.piConstant = .exact
        r.check("piConstant is a settable property; LobatoParameterization is Hashable", f32 == exact)
        r.equal("float32Pi = 3.1415927410125732", LobatoParameterization.float32Pi, 3.1415927410125732)
        r.equal("float32PiSquared = 9.869604110717773", LobatoParameterization.float32PiSquared, 9.869604110717773)
        r.refuses("parameters(atomicNumber: 0)", .unsupportedElement(0)) { try LobatoParameterization.parameters(atomicNumber: 0) }
        r.refuses("scaledParameters(atomicNumber: 104)", .unsupportedElement(104)) { try LobatoParameterization.scaledParameters(atomicNumber: 104) }
        r.refuses("projectedScatteringFactor(atomicNumber: 200)", .unsupportedElement(200)) {
            try exact.projectedScatteringFactor(atomicNumber: 200, k2: 0)
        }
    }

    r.section("Multislice: slicing") { r in
        r.equal("Slicing.maximumSliceCount", Slicing.maximumSliceCount, 1_000_000)
        // ⌈10 / 3⌉ = 4 slices of 2.5
        r.equal("thicknesses(depth: 10, requested: 3) = 4 × 2.5", try Slicing.thicknesses(depth: 10, requested: 3), [2.5, 2.5, 2.5, 2.5])
        r.equal("thicknesses(depth: 10, requested: 20) = one slice", try Slicing.thicknesses(depth: 10, requested: 20), [10])
        r.refuses("thicknesses requested 0", .invalidSliceThickness(0)) { try Slicing.thicknesses(depth: 10, requested: 0) }
        r.refusesCase("thicknesses requested NaN", "invalidSliceThickness") { try Slicing.thicknesses(depth: 10, requested: .nan) }
        r.refusesCase("thicknesses requested 1e-12 (more than a million slices)", "invalidSliceThickness") {
            try Slicing.thicknesses(depth: 10, requested: 1e-12)
        }
        r.equal("validated keeps a list that sums to the depth", try Slicing.validated([4, 6], depth: 10), [4, 6])
        r.refusesCase("validated: sum 9 against depth 10", "sliceThicknessesDoNotSumToBox") { try Slicing.validated([4, 5], depth: 10) }
        r.refuses("validated: empty list", .invalidSliceThickness(0)) { try Slicing.validated([], depth: 10) }
        r.refusesCase("validated: a negative slice", "invalidSliceThickness") { try Slicing.validated([12, -2], depth: 10) }
        // fractional coordinate mod 1: x = −1 in a 4 Å box → 3; x = 5 → 1; z within 1e-10 of the top face → 0.
        let box = try AtomicStructure(box: SIMD3(4, 4, 4), atoms: [
            .init(atomicNumber: 6, position: SIMD3(-1, 5, 4 - 1e-11)), .init(atomicNumber: 6, position: SIMD3(4, 0, 1)),
        ])
        let w = Slicing.wrapped(box)
        r.check("wrapped: −1 → 3, 5 → 1, top face → 0, x = box → 0",
                w[0].position == SIMD3(3, 1, 0) && w[1].position == SIMD3(0, 0, 1), "got \(w.map(\.position))")
        // edges 2.5, 5, 7.5, 10 (minus 1e-12): z = 2.5 is already in slice 1; z past the last edge → last slice
        r.equal("sliceIndices at 0, 2.4, 2.5, 9.99, 10", Slicing.sliceIndices(z: [0, 2.4, 2.5, 9.99, 10], thicknesses: [2.5, 2.5, 2.5, 2.5]),
                [0, 0, 1, 3, 3])
    }

    r.section("Multislice: projected potential") { r in
        let s = try singleAtomStructure()
        let v = try ProjectedPotential(structure: s, nx: 32, ny: 32, sliceThickness: 2)
        r.check("ProjectedPotential: grid 32², 2 slices of 2 Å, thickness 4",
                v.grid.nx == 32 && v.grid.extentX == 4 && v.slices.count == 2 && v.slices.allSatisfy { $0.thickness == 2 }
                && v.thickness == 4)
        // z = 2 lies on the edge 2 − 1e-12, so the atom is in slice 1.
        r.equal("sliceOfAtom = [1]", v.sliceOfAtom, [1])
        r.check("slice 0 (no atoms) is zero, slice 1 is not", v.slices[0].values.allSatisfy { $0 == 0 }
                && v.slices[1].values.contains { $0 > 0 })
        // Σ_r V = V̂(k = 0) = f(0) / (Δx Δy): the delta weights sum to 1.
        let f0 = try LobatoParameterization().projectedScatteringFactor(atomicNumber: 13, k2: 0)
        let dxdy = 4.0 / 32 * 4.0 / 32
        r.nearRelative("Σ V over a slice = f(0) / (Δx Δy)", v.slices[1].values.reduce(0, +), f0 / dxdy, 1e-9)
        // positions wrap periodically: x = −2 ≡ x = 2 (the box is 4 Å)
        let wrapped = try ProjectedPotential(structure: AtomicStructure(box: s.box, atoms: [.init(atomicNumber: 13, position: SIMD3(-2, 6, 2))]),
                                             nx: 32, ny: 32, sliceThickness: 2)
        r.check("an atom outside the box is wrapped (periodic)", maxDifference(wrapped.slices[1].values, v.slices[1].values) < 1e-9)
        let explicit = try ProjectedPotential(structure: s, nx: 32, ny: 32, sliceThicknesses: [1, 3])
        r.check("init(sliceThicknesses: [1, 3]): two slices, atom z = 2 in slice 1", explicit.slices.map(\.thickness) == [1, 3]
                && explicit.sliceOfAtom == [1])
        let f32 = try ProjectedPotential(structure: s, nx: 32, ny: 32, sliceThickness: 2,
                                         parameterization: LobatoParameterization(piConstant: .abTEMFloat32))
        let d = maxDifference(f32.slices[1].values, v.slices[1].values)
        r.check("parameterization: .abTEMFloat32 differs from .exact by under 1e-6 relative", d > 0
                && d < 1e-6 * (v.slices[1].values.max() ?? 1), "got \(d)")
        let rebuilt = ProjectedPotential(grid: v.grid, slices: v.slices)
        r.check("init(grid:slices:) from another potential's slices; sliceOfAtom empty", rebuilt.slices == v.slices
                && rebuilt.sliceOfAtom.isEmpty && rebuilt.thickness == 4)
        r.check("ProjectedPotential and PotentialSlice are Hashable", v == v && v != rebuilt && Set(v.slices).count == 2)
        r.check("PotentialSlice stores thickness and values", v.slices[1].thickness == 2 && v.slices[1].values.count == 1024)
        r.refuses("unsupported element Z = 104", .unsupportedElement(104)) {
            try ProjectedPotential(structure: singleAtomStructure(z: 104), nx: 32, ny: 32, sliceThickness: 2)
        }
        r.refusesCase("grid length 7", "unsupportedGridLength") { try ProjectedPotential(structure: s, nx: 7, ny: 32, sliceThickness: 2) }
        r.refusesCase("slice thicknesses not summing to the box", "sliceThicknessesDoNotSumToBox") {
            try ProjectedPotential(structure: s, nx: 32, ny: 32, sliceThicknesses: [1, 1])
        }
        r.refuses("slice thickness 0", .invalidSliceThickness(0)) { try ProjectedPotential(structure: s, nx: 32, ny: 32, sliceThickness: 0) }
    }

    r.section("Multislice: exit planes") { r in
        r.equal(".last, 5 slices", try ExitPlanes.last.planes(sliceCount: 5), [4])
        // the incident wave (−1), then after slices 2, 4 (zero-based 1, 3), and the last
        r.equal(".every(2), 5 slices", try ExitPlanes.every(2).planes(sliceCount: 5), [-1, 1, 3, 4])
        r.equal(".every(1), 3 slices", try ExitPlanes.every(1).planes(sliceCount: 3), [-1, 0, 1, 2])
        r.equal(".every(5), 5 slices (m ≥ n: the last only)", try ExitPlanes.every(5).planes(sliceCount: 5), [4])
        r.equal(".afterSlices([0, 2]), 5 slices (the last appended)", try ExitPlanes.afterSlices([0, 2]).planes(sliceCount: 5), [0, 2, 4])
        r.equal(".afterSlices([-1, 4]), 5 slices", try ExitPlanes.afterSlices([-1, 4]).planes(sliceCount: 5), [-1, 4])
        r.refusesCase(".every(0)", "invalidExitPlanes") { try ExitPlanes.every(0).planes(sliceCount: 5) }
        r.refusesCase(".afterSlices([2, 1]) descending", "invalidExitPlanes") { try ExitPlanes.afterSlices([2, 1]).planes(sliceCount: 5) }
        r.refusesCase(".afterSlices([1, 1]) duplicate", "invalidExitPlanes") { try ExitPlanes.afterSlices([1, 1]).planes(sliceCount: 5) }
        r.refusesCase(".afterSlices([5]) out of range", "invalidExitPlanes") { try ExitPlanes.afterSlices([5]).planes(sliceCount: 5) }
        r.check("ExitPlanes is Hashable", ExitPlanes.every(2) == .every(2) && ExitPlanes.every(2) != .last)
    }

    r.section("Multislice: run") { r in
        let beam = try ElectronBeam(energy: 200_000)
        let empty = try AtomicStructure(box: SIMD3(4, 4, 6), atoms: [])
        let vacuum = try ProjectedPotential(structure: empty, nx: 32, ny: 32, sliceThickness: 2)
        r.equal("vacuum potential has 3 slices", vacuum.slices.count, 3)
        let plane = Waves.planeWave(beam: beam, grid: vacuum.grid)
        // t = IFFT(FFT(1) · A) = 1 and P(k = 0) = 1: a plane wave in vacuum stays 1.
        let out = try Multislice.run(plane, through: vacuum)
        r.equal("run(.last): one exit wave", out.count, 1)
        r.check("vacuum: a plane wave stays 1 (|ψ − 1| < 1e-12)",
                maxDifference(out[0].field.real, plane.field.real) < 1e-12 && out[0].field.maxMagnitude < 1 + 1e-12
                && maxDifference(out[0].field.imaginary, plane.field.imaginary) < 1e-12)
        r.check("run keeps grid, beam and normalisation", out[0].grid == plane.grid && out[0].beam == beam
                && out[0].normalization == .values)
        let every = try Multislice.run(plane, through: vacuum, exitPlanes: .every(1), backend: .cpu)
        r.equal("run(.every(1)): incident + 3 exit waves", every.count, 4)
        r.check("run(.every(1)): the first wave is the incident one", every[0] == plane)
        let some = try Multislice.run(plane, through: vacuum, exitPlanes: .afterSlices([0]))
        r.equal("run(.afterSlices([0])): slice 0 and the last", some.count, 2)
        let t = try Multislice.transmissionFunctions(vacuum, beam: beam)
        r.check("transmissionFunctions: one per slice, 1 + 0i in vacuum", t.count == 3
                && t.allSatisfy { maxDifference($0.real, plane.field.real) < 1e-12 && $0.maxMagnitude < 1 + 1e-12 })
        // P = exp(−iπλΔz k²) A: at dz = 0 it is A; at (1, 0): k = 1/4 Å⁻¹ (extent 4)
        let g = vacuum.grid
        let p0 = Multislice.propagator(grid: g, beam: beam, thickness: 0)
        let aperture = Antialias.aperture(on: g)
        r.check("propagator(thickness: 0) = the aperture, real", p0.real == aperture && p0.imaginary.allSatisfy { $0 == 0 })
        let p = Multislice.propagator(grid: g, beam: beam, thickness: 2)
        let k2 = 1.0 / 16
        let chi = -Double.pi * 2 * lambda200 * k2
        r.check("propagator at (1, 0) = exp(−iπλΔz k²) (inside the aperture)",
                abs(p.real[1 * 32] - cos(chi)) < 1e-14 && abs(p.imaginary[1 * 32] - sin(chi)) < 1e-14)
        // a real atom: the CPU path conserves the intensity to the band limit, and a batch equals single runs bit for bit
        let atom = try ProjectedPotential(structure: singleAtomStructure(), nx: 32, ny: 32, sliceThickness: 2)
        let one = try Multislice.run(Waves.planeWave(beam: beam, grid: atom.grid), through: atom)[0]
        r.check("an Al atom: the exit wave has total intensity ≤ the incident one (band limit) and ≥ 0.9 of it",
                one.totalIntensity <= 1024 * (1 + 1e-9) && one.totalIntensity > 0.9 * 1024, "got \(one.totalIntensity)")
        r.check("an Al atom changes the phase (the wave is not 1)", maxDifference(one.field.imaginary, [Double](repeating: 0, count: 1024)) > 1e-4)
        let batch = try Multislice.runBatch([Waves.planeWave(beam: beam, grid: atom.grid), one], through: atom)
        let second = try Multislice.run(one, through: atom)[0]
        r.check("runBatch(.cpu) = run, bit for bit", batch.count == 2 && batch[0] == one && batch[1] == second)
        r.equal("runBatch([]) = []", try Multislice.runBatch([], through: atom).count, 0)
        r.equal("maximumCPUThreads = 8", Multislice.maximumCPUThreads, 8)
        r.check("Backend cases are Hashable", Backend.cpu != Backend.metal && Set([Backend.cpu, .cpu]).count == 1)
        let otherGrid = try Grid(nx: 32, ny: 32, extentX: 5, extentY: 4)
        r.refuses("run on another grid", .gridMismatch) { try Multislice.run(Waves.planeWave(beam: beam, grid: otherGrid), through: vacuum) }
        r.refuses("runBatch on another grid", .gridMismatch) {
            try Multislice.runBatch([Waves.planeWave(beam: beam, grid: otherGrid)], through: vacuum)
        }
        r.refusesCase("run with .every(0)", "invalidExitPlanes") { try Multislice.run(plane, through: vacuum, exitPlanes: .every(0)) }
    }

    r.section("Multislice: Metal") { r in
        guard Multislice.metalIsAvailable else {
            r.skip(".metal backend", "Multislice.metalIsAvailable is false (no Metal device on this machine)")
            r.refusesCase(".metal without a device is refused", "metalUnavailable") {
                let beam = try ElectronBeam(energy: 200_000)
                let v = try ProjectedPotential(structure: AtomicStructure(box: SIMD3(4, 4, 4), atoms: []), nx: 32, ny: 32, sliceThickness: 2)
                return try Multislice.run(Waves.planeWave(beam: beam, grid: v.grid), through: v, backend: .metal)
            }
            return
        }
        let beam = try ElectronBeam(energy: 200_000)
        let atom = try ProjectedPotential(structure: singleAtomStructure(), nx: 32, ny: 32, sliceThickness: 2)
        let plane = Waves.planeWave(beam: beam, grid: atom.grid)
        let cpu = try Multislice.run(plane, through: atom, backend: .cpu)[0]
        let gpu = try Multislice.run(plane, through: atom, backend: .metal)[0]
        r.check(".metal agrees with .cpu to float32 (1e-4 absolute on ψ)", maxDifference(cpu.field.real, gpu.field.real) < 1e-4
                && maxDifference(cpu.field.imaginary, gpu.field.imaginary) < 1e-4)
        let vac = try ProjectedPotential(structure: AtomicStructure(box: SIMD3(4, 4, 4), atoms: []), nx: 32, ny: 32, sliceThickness: 2)
        let v = try Multislice.run(Waves.planeWave(beam: beam, grid: vac.grid), through: vac, backend: .metal)[0]
        r.check(".metal: a plane wave in vacuum stays 1 (float32)", maxDifference(v.field.real, plane.field.real) < 1e-5)
        let gb = try Multislice.runBatch([plane, cpu], through: atom, backend: .metal)
        let cb = try Multislice.runBatch([plane, cpu], through: atom, backend: .cpu)
        r.check(".metal runBatch agrees with .cpu runBatch (1e-4)", gb.count == 2
                && zip(gb, cb).allSatisfy { maxDifference($0.0.field.real, $0.1.field.real) < 1e-4 })
    }

    r.section("Multislice: frozen phonons") { r in
        let s = try AtomicStructure(box: SIMD3(4, 4, 4), atoms: [
            .init(atomicNumber: 13, position: SIMD3(1, 1, 1)), .init(atomicNumber: 8, position: SIMD3(3, 2, 3)),
        ])
        r.check("Directions: raw values 1, 2, 4 and .all", FrozenPhonons.Directions.x.rawValue == 1
                && FrozenPhonons.Directions.y.rawValue == 2 && FrozenPhonons.Directions.z.rawValue == 4
                && FrozenPhonons.Directions.all == [.x, .y, .z] && FrozenPhonons.Directions(rawValue: 3) == [.x, .y])
        // σ = 0 leaves every atom still: the ensemble's potentials equal the static one.
        let still = try FrozenPhonons(structure: s, sigmas: [13: 0, 8: 0], configurations: 2, seed: 7)
        r.check("σ = 0: count 2, every displacement zero", still.count == 2 && still.displacements.count == 2
                && still.displacements.allSatisfy { $0.allSatisfy { $0 == .zero } })
        r.check("σ = 0: configuration(c) is the structure", try still.configuration(1) == s && still.structure == s)
        let static_ = try ProjectedPotential(structure: s, nx: 32, ny: 32, sliceThickness: 2)
        let ensemble = try still.potentials(nx: 32, ny: 32, sliceThickness: 2)
        r.check("σ = 0: FrozenPhonons.potentials equal the static potential",
                ensemble.count == 2 && ensemble.allSatisfy { $0 == static_ })
        let ensemble32 = try still.potentials(nx: 32, ny: 32, sliceThickness: 2,
                                              parameterization: LobatoParameterization(piConstant: .abTEMFloat32))
        r.check("potentials(parameterization:) is passed on", ensemble32[0] != static_)
        // the same seed, the same ensemble; another seed, another one; consecutive seeds do not share configurations
        let a = try FrozenPhonons(structure: s, sigmas: [13: 0.1, 8: 0.1], configurations: 3, seed: 5)
        let b = try FrozenPhonons(structure: s, sigmas: [13: 0.1, 8: 0.1], configurations: 3, seed: 5)
        let c = try FrozenPhonons(structure: s, sigmas: [13: 0.1, 8: 0.1], configurations: 3, seed: 6)
        r.check("same seed → same displacements (Hashable equal); another seed differs", a == b && a != c)
        r.check("configurationSeed(5, 1) ≠ configurationSeed(6, 0), and is deterministic",
                FrozenPhonons.configurationSeed(5, 1) != FrozenPhonons.configurationSeed(6, 0)
                && FrozenPhonons.configurationSeed(5, 1) == FrozenPhonons.configurationSeed(5, 1)
                && FrozenPhonons.configurationSeed(5, 0) != FrozenPhonons.configurationSeed(5, 1))
        let xOnly = try FrozenPhonons(structure: s, sigmas: [13: 0.1, 8: 0.1], configurations: 2, seed: 5, directions: .x)
        r.check("directions: .x leaves y and z at zero", xOnly.displacements.allSatisfy { $0.allSatisfy { $0.y == 0 && $0.z == 0 } }
                && xOnly.displacements[0].contains { $0.x != 0 })
        // the configuration is the structure plus the displacement, not wrapped
        let moved = try a.configuration(2)
        r.check("configuration(c) = position + displacement", zip(moved.atoms, a.displacements[2]).enumerated().allSatisfy {
            $0.element.0.position == s.atoms[$0.offset].position + $0.element.1 })
        // statistics: x² averages σ² over many configurations of one atom
        let one = try AtomicStructure(box: SIMD3(4, 4, 4), atoms: [.init(atomicNumber: 13, position: SIMD3(2, 2, 2))])
        let many = try FrozenPhonons(structure: one, sigmas: [13: 0.2], configurations: 4000, seed: 11)
        let xs = many.displacements.map { $0[0].x }
        let mean = xs.reduce(0, +) / 4000, variance = xs.map { $0 * $0 }.reduce(0, +) / 4000
        r.check("4000 draws with σ = 0.2: mean 0 (±0.02), variance 0.04 (±10%)", abs(mean) < 0.02 && abs(variance - 0.04) < 0.004,
                "mean \(mean), variance \(variance)")
        let given = try FrozenPhonons(structure: one, displacements: [[SIMD3(0.1, 0, 0)], [SIMD3(0, 0.2, 0)]])
        let givenMoved = try given.configuration(1).atoms[0].position
        r.check("init(displacements:) keeps them; configuration applies them", given.count == 2 && givenMoved == SIMD3(2, 2.2, 2))
        r.refusesCase("no displacements", "invalidFrozenPhonons") { try FrozenPhonons(structure: one, displacements: []) }
        r.refusesCase("displacements for the wrong atom count", "invalidFrozenPhonons") {
            try FrozenPhonons(structure: one, displacements: [[.zero, .zero]])
        }
        r.refusesCase("a non-finite displacement", "invalidFrozenPhonons") {
            try FrozenPhonons(structure: one, displacements: [[SIMD3(.nan, 0, 0)]])
        }
        r.refusesCase("configurations: 0", "invalidFrozenPhonons") { try FrozenPhonons(structure: one, sigmas: [13: 0.1], configurations: 0, seed: 1) }
        r.refusesCase("σ missing for an element", "invalidFrozenPhonons") { try FrozenPhonons(structure: s, sigmas: [13: 0.1], configurations: 1, seed: 1) }
        r.refusesCase("σ negative", "invalidFrozenPhonons") { try FrozenPhonons(structure: one, sigmas: [13: -0.1], configurations: 1, seed: 1) }
        r.refusesCase("σ NaN", "invalidFrozenPhonons") { try FrozenPhonons(structure: one, sigmas: [13: .nan], configurations: 1, seed: 1) }
    }

    r.section("Multislice: generator") { r in
        var g = GaussianGenerator(seed: 42)
        var h = GaussianGenerator(seed: 42)
        var other = GaussianGenerator(seed: 43)
        let first = g.nextUInt64()
        r.check("GaussianGenerator: same seed, same stream; another seed differs", first == h.nextUInt64()
                && first != other.nextUInt64() && g.nextUInt64() == h.nextUInt64())
        var sum = 0.0, sumSquares = 0.0
        let n = 20_000
        for _ in 0..<n { let z = g.nextNormal(); sum += z; sumSquares += z * z }
        r.check("20 000 normals: mean 0 (±0.03), variance 1 (±0.05)", abs(sum / Double(n)) < 0.03 && abs(sumSquares / Double(n) - 1) < 0.05,
                "mean \(sum / Double(n)), variance \(sumSquares / Double(n))")
        r.check("nextNormal is finite", h.nextNormal().isFinite)
    }

    r.section("Multislice: namespace") { r in
        r.check("SimulatoMultislice namespace type exists", String(describing: SimulatoMultislice.self) == "SimulatoMultislice")
    }

    r.section("Multislice: the lazy Potential") { r in
        let beam = try ElectronBeam(energy: 200_000)
        let structure = try singleAtomStructure().repeated(2, 2, 2)
        let lazy = try Potential(structure: structure, nx: 64, ny: 64, sliceThickness: 2)
        let uneven = try Potential(structure: structure, nx: 64, ny: 64, sliceThicknesses: [1, 3, 4])
        let eager = try ProjectedPotential(structure: structure, nx: 64, ny: 64, sliceThickness: 2)
        r.check("build() equals ProjectedPotential(structure:...) exactly", try lazy.build() == eager)
        r.check("sliceCount, sliceThicknesses, thickness, grid", lazy.sliceCount == 4 && lazy.sliceThicknesses == [2, 2, 2, 2]
                && lazy.thickness == 8 && lazy.grid == eager.grid && uneven.sliceCount == 3)
        r.check("atoms wrapped, one slice each, atoms(inSlice:) partitions them",
                lazy.atoms.count == 8 && lazy.sliceOfAtom == eager.sliceOfAtom
                && (0..<4).map { lazy.atoms(inSlice: $0).count }.reduce(0, +) == 8)
        r.check("parameterization kept", lazy.parameterization == LobatoParameterization())
        let wave = Waves.planeWave(beam: beam, grid: lazy.grid)
        let a = try Multislice.run(wave, through: lazy, exitPlanes: .every(2))
        let b = try Multislice.run(wave, through: eager, exitPlanes: .every(2))
        r.check("Multislice.run through a Potential equals the run through its build() on .cpu", a == b)
        let batch = try Multislice.runBatch([wave, wave], through: lazy)
        r.check("Multislice.runBatch through a Potential", batch.count == 2 && batch[0] == b.last!)
        let fp = try FrozenPhonons(structure: structure, sigmas: [13: 0], configurations: 2, seed: 1)
        let lazies = try fp.lazyPotentials(nx: 64, ny: 64, sliceThickness: 2)
        r.check("lazyPotentials: σ = 0 configurations equal the static Potential", lazies.count == 2 && lazies.allSatisfy { $0 == lazy })
        let slice = PotentialSlice(thickness: 2, values: eager.slices[0].values)
        r.check("PotentialSlice(thickness:values:) and ProjectedPotential(grid:slices:)",
                ProjectedPotential(grid: eager.grid, slices: [slice]).slices == [eager.slices[0]])
        if Multislice.metalIsAvailable {
            let gpu = try Multislice.run(wave, through: lazy, backend: .metal)[0]
            var worst = 0.0
            for i in gpu.field.real.indices {
                worst = max(worst, hypot(gpu.field.real[i] - b.last!.field.real[i], gpu.field.imaginary[i] - b.last!.field.imaginary[i]))
            }
            r.check("a Potential on .metal (built on the GPU) within 1e-4 of the CPU", worst <= 1e-4 * b.last!.field.maxMagnitude,
                    "worst \(worst)")
        } else {
            r.skip("a Potential on .metal", "no Metal device")
        }
        r.refusesCase("slice thicknesses that do not sum to the box", "sliceThicknessesDoNotSumToBox") {
            try Potential(structure: structure, nx: 64, ny: 64, sliceThicknesses: [1, 1])
        }
    }
}

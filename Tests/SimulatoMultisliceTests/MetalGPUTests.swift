//
//  MetalGPUTests.swift
//  Role: registration-performance.md sections 3 and 4: the Metal path on the GPU end to end (the potential built on the
//        device for a `Potential`, the transmissions built on the device for both forms, the probes on the device)
//        against the CPU path, and its truth. `Potential.build()` against `ProjectedPotential` runs without a device.
//

import Testing
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe

enum GPUCases {
    static func potential(_ id: String, pi: LobatoParameterization.PiConstant = .exact) throws -> Potential {
        let c = Pins.case(id)
        let g = c["gpts"] as! [Int]
        return try Potential(structure: Pins.structure(id), nx: g[0], ny: g[1],
                             sliceThickness: c["slice_thickness_requested"] as! Double,
                             parameterization: LobatoParameterization(piConstant: pi))
    }

    static func aluminium(cells: (Int, Int, Int)) throws -> AtomicStructure {
        let cell = try AtomicStructure(box: SIMD3(4.05, 4.05, 4.05), atoms: [
            .init(atomicNumber: 13, position: SIMD3(0, 0, 0)), .init(atomicNumber: 13, position: SIMD3(0, 2.025, 2.025)),
            .init(atomicNumber: 13, position: SIMD3(2.025, 0, 2.025)), .init(atomicNumber: 13, position: SIMD3(2.025, 2.025, 0)),
        ])
        return cell.repeated(cells.0, cells.1, cells.2)
    }

    struct Worst {
        var wave = 0.0, intensity = 0.0, pattern = 0.0, total = 0.0, totalPerSlice = 0.0
    }

    /// Slice 1's bounds (section 4) on every plane, the total-intensity bound of amendment 2.
    static func compare(_ gpu: [Waves], _ cpu: [Waves], planes: [Int], label: String) throws -> Worst {
        var worst = Worst()
        #expect(gpu.count == cpu.count, "\(label): plane count")
        for ((a, b), plane) in zip(zip(gpu, cpu), planes) {
            let w = relMax(a.field, b.field)
            let i = relMax(a.intensity().values, b.intensity().values)
            let p = relMax(try a.diffractionPattern().values, try b.diffractionPattern().values)
            let t = relDiff(a.totalIntensity, b.totalIntensity)
            worst.wave = max(worst.wave, w)
            worst.intensity = max(worst.intensity, i)
            worst.pattern = max(worst.pattern, p)
            worst.total = max(worst.total, t)
            if plane >= 0 { worst.totalPerSlice = max(worst.totalPerSlice, t / Double(plane + 1)) }
            let bound = MetalAccuracy.registered
            #expect(w <= bound.wave, "\(label): exit wave \(w)")
            #expect(i <= bound.intensity, "\(label): intensity \(i)")
            #expect(p <= bound.pattern, "\(label): diffraction pattern \(p)")
            #expect(t <= bound.totalIntensityPerSlice * Double(plane + 1), "\(label): total intensity \(t)")
        }
        return worst
    }
}

@Suite("Potential: the lazy potential (registration-performance.md 1.1)")
struct PotentialTests {
    @Test(arguments: ["A1", "A2", "B1", "D1"])
    func buildEqualsTheEagerPotentialBitForBit(id: String) throws {
        for pi in [LobatoParameterization.PiConstant.exact, .abTEMFloat32] {
            let lazy = try GPUCases.potential(id, pi: pi)
            let eager = try Pins.potential(id, pi: pi)
            #expect(try lazy.build() == eager, "\(id)")
            #expect(lazy.sliceCount == eager.slices.count && lazy.sliceOfAtom == eager.sliceOfAtom)
        }
    }
}

@Suite("Metal on the GPU end to end against the CPU path (registration-performance.md)",
       .enabled(if: Multislice.metalIsAvailable,
                "no Metal device (MTLCreateSystemDefaultDevice() returned nil): the Metal path is not tested here"))
struct MetalGPUTests {
    /// Section 3, the GPU's V against the CPU's (1e−5 × max V), and truth 1 (Σ V ΔxΔy = Σ f(0), 1e−5 relative).
    @Test(arguments: ["A1", "A2", "B1", "D1"])
    func gpuPotentialMatchesTheCPUPotential(id: String) throws {
        let lazy = try GPUCases.potential(id)
        let cpu = try lazy.build()
        let gpu = try Multislice.gpuPotentialSlices(lazy)
        #expect(gpu.count == cpu.slices.count)
        var worst = 0.0
        for (j, (g, c)) in zip(gpu, cpu.slices).enumerated() {
            let d = relMax(g, c.values)
            worst = max(worst, d)
            #expect(d <= 1e-5, "\(id) slice \(j): \(d)")
            let atoms = lazy.atoms(inSlice: j)
            guard !atoms.isEmpty else { continue }
            let f0 = try atoms.map { try lazy.parameterization.projectedScatteringFactor(atomicNumber: $0.atomicNumber, k2: 0) }
                .reduce(0, +)
            let integral = g.reduce(0, +) * lazy.grid.samplingX * lazy.grid.samplingY
            #expect(relDiff(integral, f0) <= 1e-5, "\(id) slice \(j): k = 0 identity \(relDiff(integral, f0))")
        }
        print("GPU V vs CPU V, \(id): \(worst)")
    }

    @Test(arguments: CaseBeam.all)
    func potentialOnMetalMatchesTheCPUPath(r: CaseBeam) throws {
        let c = Pins.case(r.id)
        let lazy = try GPUCases.potential(r.potentialCase)
        let incident = try Pins.incident(r.id, r.label, grid: lazy.grid)
        let planes = ExitPlanes.afterSlices(c["exit_planes"] as! [Int])
        let cpu = try Multislice.run(incident, through: lazy, exitPlanes: planes, backend: .cpu)
        let gpu = try Multislice.run(incident, through: lazy, exitPlanes: planes, backend: .metal)
        let w = try GPUCases.compare(gpu, cpu, planes: planes.planes(sliceCount: lazy.sliceCount), label: "\(r.id)/\(r.label)")
        print("GPU potential vs CPU, \(r.id)/\(r.label): wave \(w.wave), intensity \(w.intensity), pattern \(w.pattern), per slice \(w.totalPerSlice)")
    }

    /// Section 3, new: 1024² × 100 slices, the potential built on the GPU, against the CPU path.
    @Test func potentialOnMetalMatchesTheCPUPathAt1024Squared() throws {
        let lazy = try Potential(structure: GPUCases.aluminium(cells: (12, 12, 50)), nx: 1024, ny: 1024, sliceThickness: 2.025)
        #expect(lazy.sliceCount == 100)
        let beam = try ElectronBeam(energy: 300e3)
        let incident = Waves.planeWave(beam: beam, grid: lazy.grid)
        let planes = ExitPlanes.every(50)
        let cpu = try Multislice.run(incident, through: lazy, exitPlanes: planes, backend: .cpu)
        let gpu = try Multislice.run(incident, through: lazy, exitPlanes: planes, backend: .metal)
        let w = try GPUCases.compare(gpu, cpu, planes: planes.planes(sliceCount: 100), label: "Al 1024²")
        print("GPU potential vs CPU, 1024² x 100 slices: wave \(w.wave), intensity \(w.intensity), pattern \(w.pattern), per slice \(w.totalPerSlice)")
    }

    /// Truth 2: vacuum stays 1; the GPU's probe has Σ_k |Ψ(k)|² = 1 and equals `Probe.build` within 1e−5 × max.
    @Test func vacuumAndProbesOnTheGPU() throws {
        let empty = try AtomicStructure(box: SIMD3(16.2, 16.2, 20), atoms: [])
        let lazy = try Potential(structure: empty, nx: 256, ny: 256, sliceThickness: 2)
        let beam = try ElectronBeam(energy: 200e3)
        let exit = try Multislice.run(.planeWave(beam: beam, grid: lazy.grid), through: lazy, backend: .metal)[0]
        var worst = 0.0
        for i in exit.field.real.indices {
            worst = max(worst, abs(exit.field.real[i] - 1), abs(exit.field.imaginary[i]))
        }
        #expect(worst <= 1e-6, "vacuum plane wave \(worst)")
        var ab = Aberrations()
        ab.defocus = -50
        ab.Cs = 1e4
        let probe = try Probe(beam: beam, grid: lazy.grid, aperture: Aperture(semiangleCutoff: 20), aberrations: ab)
        let positions = [SIMD2(8.37, 7.91), SIMD2(0, 0), SIMD2(16.1, 3.3)]
        let gpu = try Multislice.gpuProbes(spectrum: probe.centredSpectrum(),
                                           ramps: positions.map { try probe.shiftRamps(at: $0) }, grid: lazy.grid)
        for (g, r) in zip(gpu, positions) {
            let cpu = try probe.build(at: r).field
            let d = relMax(g, cpu)
            let norm = g.squaredMagnitude().reduce(0, +) * Double(lazy.grid.count)
            #expect(d <= 1e-5, "probe at \(r): \(d)")
            #expect(abs(norm - 1) <= 1e-6, "probe at \(r): Σ|Ψ(k)|² − 1 = \(norm - 1)")
            print("GPU probe vs CPU at \(r): \(d), Σ|Ψ(k)|² − 1 = \(norm - 1)")
        }
    }

    /// Truth 3: two runs give bit-equal waves (the host-merged scatter is deterministic).
    @Test func runsAreDeterministic() throws {
        let lazy = try GPUCases.potential("D1")
        let incident = try Pins.incident("D1", "probe", grid: lazy.grid)
        let a = try Multislice.run(incident, through: lazy, backend: .metal)[0]
        let b = try Multislice.run(incident, through: lazy, backend: .metal)[0]
        #expect(a.field == b.field)
        let stacked = try AtomicStructure(box: SIMD3(8, 8, 8), atoms: (0..<8).map {
            .init(atomicNumber: 79, position: SIMD3(4.01, 3.99, Double($0)))
        })
        let column = try Potential(structure: stacked, nx: 128, ny: 128, sliceThickness: 8)
        let w = Waves.planeWave(beam: incident.beam, grid: column.grid)
        #expect(try Multislice.run(w, through: column, backend: .metal)[0].field
                == Multislice.run(w, through: column, backend: .metal)[0].field)
        let cpu = try Multislice.run(w, through: column, backend: .cpu)[0]
        let gpu = try Multislice.run(w, through: column, backend: .metal)[0]
        #expect(relMax(gpu.field, cpu.field) <= 1e-4, "a column of eight gold atoms in one slice (the merge)")
    }

    /// Amendment 1: the per-step renormalisation's effect on Σ|ψ|², both potential forms, 512² × 100 slices.
    @Test func renormalizationTightensTheTotalIntensity() throws {
        let lazy = try Potential(structure: GPUCases.aluminium(cells: (8, 8, 50)), nx: 512, ny: 512, sliceThickness: 2.025)
        let eager = try lazy.build()
        let beam = try ElectronBeam(energy: 300e3)
        var ab = Aberrations()
        ab.defocus = 30
        let probe = try Probe(beam: beam, grid: lazy.grid, aperture: Aperture(semiangleCutoff: 25), aberrations: ab)
        let planes = ExitPlanes.every(25)
        let indices = try planes.planes(sliceCount: 100)
        for (label, incident) in [("plane", Waves.planeWave(beam: beam, grid: lazy.grid)),
                                  ("probe", try probe.build(at: SIMD2(15.3, 16.9)))] {
            let cpu = try Multislice.run(incident, through: eager, exitPlanes: planes, backend: .cpu)
            for renormalize in [false, true] {
                var gpu: [Waves] = []
                let clock = ContinuousClock()
                _ = try Multislice.runMetal(incident, through: lazy, exitPlanes: .last, renormalize: renormalize)  // warm
                let time = try clock.measure {
                    gpu = try Multislice.runMetal(incident, through: lazy, exitPlanes: planes, renormalize: renormalize)
                }
                print("512² x 100 slices, \(label), renormalize \(renormalize): \(time) on the GPU (potential included)")
                let a = try GPUCases.compare(gpu, cpu, planes: indices, label: "\(label) atoms r=\(renormalize)")
                let b = try GPUCases.compare(try Multislice.runMetal(incident, through: eager, exitPlanes: planes,
                                                                     renormalize: renormalize),
                                             cpu, planes: indices, label: "\(label) values r=\(renormalize)")
                if renormalize {
                    // Amendment 1: the measured drift with the renormalisation, guarded at 3e−7 per slice (worst
                    // measured 1.4e−7; without it 9.3e−7, inside the accepted 1.5e−6).
                    #expect(a.totalPerSlice <= 3e-7 && b.totalPerSlice <= 3e-7, "\(label): renormalised drift")
                }
                print("Σ|ψ|² per slice, 512² x 100, \(label), renormalize \(renormalize): atoms \(a.totalPerSlice), values \(b.totalPerSlice); wave \(a.wave) / \(b.wave), pattern \(a.pattern) / \(b.pattern)")
            }
        }
    }

    /// Uneven slices (every distinct thickness its own propagator), both forms, against the CPU path.
    @Test func unevenSlicesMatchTheCPUPath() throws {
        let structure = try Pins.structure("D1")
        let thicknesses = [1.0, 2.5, 0.7, 4.0, 3.325, 2.0, 6.0]
        let lazy = try Potential(structure: structure, nx: 120, ny: 160, sliceThicknesses: thicknesses)
        let eager = try lazy.build()
        let incident = try Pins.incident("D1", "probe", grid: lazy.grid)
        let planes = ExitPlanes.every(2)
        let indices = try planes.planes(sliceCount: thicknesses.count)
        let cpu = try Multislice.run(incident, through: eager, exitPlanes: planes, backend: .cpu)
        let a = try GPUCases.compare(Multislice.run(incident, through: lazy, exitPlanes: planes, backend: .metal), cpu,
                                     planes: indices, label: "uneven, atoms")
        let b = try GPUCases.compare(Multislice.run(incident, through: eager, exitPlanes: planes, backend: .metal), cpu,
                                     planes: indices, label: "uneven, values")
        print("uneven slices: wave \(a.wave) / \(b.wave), pattern \(a.pattern) / \(b.pattern)")
    }

    /// Gate B finding 4: a potential with no slice returns the incident wave on `.metal`, as on `.cpu` (it trapped).
    @Test func noSliceReturnsTheIncidentWave() throws {
        let grid = try Grid(nx: 64, ny: 64, extentX: 8, extentY: 8)
        let empty = ProjectedPotential(grid: grid, slices: [])
        let wave = Waves.planeWave(beam: try ElectronBeam(energy: 200e3), grid: grid)
        let cpu = try Multislice.run(wave, through: empty, backend: .cpu)
        let gpu = try Multislice.run(wave, through: empty, backend: .metal)
        #expect(gpu.count == cpu.count && gpu.first?.field == wave.field)
        let batch = try Multislice.runBatch([wave, wave], through: empty, backend: .metal)
        #expect(batch.count == 2 && batch.allSatisfy { $0.field == wave.field })
    }

    /// Gate B finding 5: a wave of amplitude 1e8 stays finite and within the bounds (the renormalisation's products
    /// overflowed float32 and gave NaN).
    @Test func largeAmplitudesStayFinite() throws {
        let lazy = try GPUCases.potential("B1")
        let beam = try Pins.beam("B1")
        let one = Waves.planeWave(beam: beam, grid: lazy.grid)
        var field = one.field
        field.real = field.real.map { $0 * 1e8 }
        let big = Waves(field: field, grid: lazy.grid, beam: beam, normalization: .values)
        let gpu = try Multislice.run(big, through: lazy, backend: .metal)[0]
        let cpu = try Multislice.run(big, through: lazy, backend: .cpu)[0]
        #expect(gpu.field.real.allSatisfy(\.isFinite) && gpu.field.imaginary.allSatisfy(\.isFinite), "finite")
        #expect(relMax(gpu.field, cpu.field) <= 1e-4, "amplitude 1e8: \(relMax(gpu.field, cpu.field))")
    }

    /// Gate B finding 6: the compiled executables are bounded (each held device memory; twelve batch sizes held 2.9 GB),
    /// and a remainder batch reuses a power-of-two shape; results do not depend on the batch a wave ran in.
    @Test func compiledExecutablesAreBounded() throws {
        let lazy = try Potential(structure: GPUCases.aluminium(cells: (2, 2, 2)), nx: 64, ny: 64, sliceThickness: 2.025)
        let wave = Waves.planeWave(beam: try ElectronBeam(energy: 300e3), grid: lazy.grid)
        let single = try Multislice.runBatch([wave], through: lazy, backend: .metal)[0]
        for n in 1...12 {
            let out = try Multislice.runBatch([Waves](repeating: wave, count: n), through: lazy, backend: .metal)
            #expect(out.count == n && relMax(out[n - 1].field, single.field) <= 1e-6, "batch of \(n)")
        }
        #expect(Multislice.metalCachedExecutables <= 12, "\(Multislice.metalCachedExecutables) executables held")
    }
}

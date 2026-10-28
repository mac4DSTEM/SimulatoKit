//
//  GeneralityTests.swift
//  Role: every policy that became a parameter, off its default. The defaults themselves are every other test.
//

import Testing
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe
import SimulatoKitTestSupport

@Suite("Policies as parameters")
struct GeneralityTests {
    static func batch(_ n: Int) throws -> (Potential, [Waves]) {
        let lazy = try Potential(structure: aluminiumStructure(cells: (2, 2, 3)), nx: 64, ny: 64, sliceThickness: 2.025)
        let probe = try Probe(beam: ElectronBeam(energy: 200e3), grid: lazy.grid, aperture: Aperture(semiangleCutoff: 20))
        return (lazy, try (0..<n).map { try probe.build(at: SIMD2(0.2 * Double($0), 0.13 * Double($0))) })
    }

    @Test func theDefaultsAreTheRegisteredValues() {
        let o = ExecutionOptions.standard
        #expect(o.maximumCPUThreads == 8 && Multislice.maximumCPUThreads == 8)
        #expect(o.metalBatchBytes == 512 << 20 && o.scanChunkBytes == 256 << 20)
        #expect(o.metalRenormalization && o.metalExecutableCacheLimit == 12 && o.metalBatchShapes == .powersOfTwo)
        let a = MetalAccuracy.registered
        #expect(a.wave == 1e-4 && a.intensity == 2e-4 && a.pattern == 2e-4 && a.detector == 2e-4 && a.detectorFloor == 1e-5)
        #expect(a.detectorBound(imageMaximum: 1) == 2e-4 && a.detectorBound(imageMaximum: 1e-7) == 2e-9)
        #expect(a.totalIntensityPerSlice == 1.5e-6 && a.sliceLimit == 100 && a.gridLimit == 1024)
        #expect(Slicing.maximumSliceCount == 1_000_000)
        #expect(ExecutionOptions(maximumCPUThreads: 0, metalBatchBytes: 0, metalExecutableCacheLimit: 0)
                == ExecutionOptions(maximumCPUThreads: 1, metalBatchBytes: 1, metalExecutableCacheLimit: 1), "clamped to 1")
    }

    @Test func oneCPUThreadEqualsTheDefaultBitForBit() throws {
        let (lazy, incidents) = try Self.batch(11)
        let eager = try lazy.build()
        let a = try Multislice.runBatch(incidents, through: eager)
        let b = try Multislice.runBatch(incidents, through: eager, options: ExecutionOptions(maximumCPUThreads: 1))
        #expect(a.map(\.field) == b.map(\.field))
    }

    @Test func theMetalBatchPlanFollowsItsOptions() throws {
        let grid = try Grid(nx: 64, ny: 64, extentX: 8, extentY: 8)
        let perWave = grid.count * 24
        #expect(Multislice.metalBatchPlan(count: 11, grid: grid) == [[0, 11, 16]])
        #expect(Multislice.metalBatchPlan(count: 11, grid: grid, options: ExecutionOptions(metalBatchShapes: .exact))
                == [[0, 11, 11]])
        #expect(Multislice.metalBatchPlan(count: 11, grid: grid, options: ExecutionOptions(metalBatchBytes: 4 * perWave))
                == [[0, 4, 4], [4, 4, 4], [8, 3, 4]])
        #expect(Multislice.metalBatchPlan(count: 11, grid: grid,
                                          options: ExecutionOptions(metalBatchBytes: 4 * perWave, metalBatchShapes: .exact))
                == [[0, 4, 4], [4, 4, 4], [8, 3, 3]])
        #expect(Multislice.metalBatchPlan(count: 3, grid: grid, options: ExecutionOptions(metalBatchBytes: 1))
                == [[0, 1, 1], [1, 1, 1], [2, 1, 1]])
    }

    @Test(.enabled(if: Multislice.metalIsAvailable, "no Metal device"))
    func metalBatchingAndShapesChangeRoundingOnly() throws {
        let (lazy, incidents) = try Self.batch(11)
        let reference = try Multislice.runBatch(incidents, through: lazy, backend: .metal)
        // The default run itself (11 waves in a padded batch of 16) returns every wave, in order.
        #expect(reference.count == incidents.count)
        let cpu = try Multislice.runBatch(incidents, through: lazy)
        for (x, y) in zip(reference, cpu) { #expect(relMax(x.field, y.field) <= MetalAccuracy.registered.wave) }
        let perWave = lazy.grid.count * 24
        for options in [ExecutionOptions(metalBatchBytes: perWave), ExecutionOptions(metalBatchBytes: 3 * perWave),
                        ExecutionOptions(metalBatchBytes: 3 * perWave, metalBatchShapes: .exact),
                        ExecutionOptions(metalBatchShapes: .exact)] {
            let out = try Multislice.runBatch(incidents, through: lazy, backend: .metal, options: options)
            #expect(out.count == incidents.count)
            for (x, y) in zip(out, reference) { #expect(relMax(x.field, y.field) <= 1e-6, "\(options)") }
        }
    }

    @Test(.enabled(if: Multislice.metalIsAvailable, "no Metal device"))
    func theExecutableCacheFollowsItsLimit() throws {
        let (lazy, incidents) = try Self.batch(7)
        let held = try Multislice.withMetalExclusive {
            for n in 1...7 {
                _ = try Multislice.runBatch(Array(incidents.prefix(n)), through: lazy, backend: .metal,
                                            options: ExecutionOptions(metalExecutableCacheLimit: 2, metalBatchShapes: .exact))
            }
            return Multislice.metalCachedExecutables
        }
        #expect(held <= 2, "\(held) held")
    }

    @Test(.enabled(if: Multislice.metalIsAvailable, "no Metal device"))
    func renormalizationOffDriftsMore() throws {
        let lazy = try Potential(structure: aluminiumStructure(cells: (4, 4, 50)), nx: 256, ny: 256, sliceThickness: 2.025)
        let probe = try Probe(beam: ElectronBeam(energy: 300e3), grid: lazy.grid, aperture: Aperture(semiangleCutoff: 25))
        let wave = try probe.build(at: SIMD2(4.1, 3.7))
        let cpu = try Multislice.run(wave, through: lazy, backend: .cpu)[0].totalIntensity
        let on = try Multislice.run(wave, through: lazy, backend: .metal)[0].totalIntensity
        let off = try Multislice.run(wave, through: lazy, backend: .metal,
                                     options: ExecutionOptions(metalRenormalization: false))[0].totalIntensity
        #expect(abs(off - cpu) > 2 * abs(on - cpu), "off \(abs(off - cpu) / cpu), on \(abs(on - cpu) / cpu)")
    }

    @Test func theGridToleranceIsAParameter() throws {
        let a = try Grid(nx: 64, ny: 64, extentX: 10, extentY: 10)
        let b = try Grid(nx: 64, ny: 64, extentX: 10 * (1 + 1e-10), extentY: 10)
        #expect(!a.matches(b) && a.matches(b, relativeTolerance: 1e-9) && a.matches(a, relativeTolerance: 0))
    }

    @Test func theSliceCeilingIsAParameter() throws {
        let s = try aluminiumStructure(cells: (1, 1, 1))
        #expect(throws: SimulatoError.invalidSliceThickness(0.3)) {
            try Potential(structure: s, nx: 32, ny: 32, sliceThickness: 0.3, maximumSliceCount: 10)
        }
        #expect(try Potential(structure: s, nx: 32, ny: 32, sliceThickness: 0.3).sliceCount == 14)
        #expect(throws: SimulatoError.invalidSliceThickness(0)) {
            try Potential(structure: s, nx: 32, ny: 32, sliceThicknesses: [Double](repeating: 4.05 / 11, count: 11),
                          maximumSliceCount: 10)
        }
        #expect(throws: SimulatoError.invalidSliceThickness(1)) { try Slicing.thicknesses(depth: 5, requested: 1, maximumCount: 4) }
    }

    /// Counts set through the public vars are clamped too, and the engine survives them;
    /// the cache limit holds on hits; a huge ceiling does not trap; odd grid tolerances keep a grid matching itself.
    @Test func edgeValuesAreSafe() throws {
        var o = ExecutionOptions.standard
        o.metalExecutableCacheLimit = -1
        o.maximumCPUThreads = Int.min
        o.metalBatchBytes = 0
        o.scanChunkBytes = -5
        #expect(o == ExecutionOptions(maximumCPUThreads: 1, metalBatchBytes: 1, scanChunkBytes: 1, metalExecutableCacheLimit: 1))
        #expect(try Slicing.thicknesses(depth: 4, requested: 1, maximumCount: Int.max).count == 4)
        #expect(throws: SimulatoError.invalidSliceThickness(1)) {
            try Slicing.thicknesses(depth: 9.3e18, requested: 1, maximumCount: Int.max)
        }
        let g = try Grid(nx: 64, ny: 64, extentX: 10, extentY: 10)
        #expect(g.matches(g, relativeTolerance: -1) && g.matches(g, relativeTolerance: .nan) && g.matches(g, relativeTolerance: 0))
        let far = try Grid(nx: 64, ny: 64, extentX: 20, extentY: 10)
        #expect(!g.matches(far, relativeTolerance: .nan) && g.matches(far, relativeTolerance: .infinity))
        let s = try aluminiumStructure(cells: (1, 1, 1))
        #expect(throws: SimulatoError.invalidSliceThickness(0.3)) {
            try ProjectedPotential(structure: s, nx: 32, ny: 32, sliceThickness: 0.3, maximumSliceCount: 10)
        }
        guard Multislice.metalIsAvailable else { return }
        let (lazy, incidents) = try Self.batch(5)
        let held = try Multislice.withMetalExclusive {
            let gpu = try Multislice.runBatch(incidents, through: lazy, backend: .metal, options: o)
            #expect(gpu.count == 5)
            _ = try Multislice.runBatch(incidents, through: lazy, backend: .metal)  // hits and misses at the default limit
            _ = try Multislice.runBatch(incidents, through: lazy, backend: .metal, options: o)  // a hit at limit 1
            return Multislice.metalCachedExecutables
        }
        #expect(held == 1, "\(held) held under limit 1")
    }
}

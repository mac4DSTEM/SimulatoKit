//
//  STEMPinsTests.swift
//  Role: registration-slice-2.md section 4: the STEM driver and every detector against abTEM (S1 with three frozen-phonon
//        configurations, S2 static, S3 vacuum, S4 single atoms), CPU path; and the Metal path against the CPU path.
//

import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoMultislice
import SimulatoProbe

@Suite("STEM against abTEM")
struct STEMPinsTests {
    static func run(_ id: String, pi: LobatoParameterization.PiConstant, backend: Backend = .cpu)
        throws -> [String: ScanMeasurement] {
        let potentials = try Pins2.potentials(id, pi: pi)
        let probe = try Pins2.probe(id, grid: potentials[0].grid)
        let detectors = try Pins2.detectors()
        let m = try STEM.scan(probe: probe, potentials: potentials, scan: Pins2.scan(id), detectors: detectors.map(\.1),
                              backend: backend)
        return Dictionary(uniqueKeysWithValues: zip(detectors.map(\.0), m))
    }

    @Test(arguments: ["S1", "S2", "S3"])
    func scanPositionsEqualAbTEM(id: String) throws {
        let p = try Pins2.scan(id).positions.flatMap { [$0.x, $0.y] }
        #expect(p == Pins2.doubles("\(id).scan_positions"))
    }

    @Test(arguments: ["S1", "S2", "S3"], [LobatoParameterization.PiConstant.abTEMFloat32, .exact])
    func detectorsEqualAbTEM(id: String, pi: LobatoParameterization.PiConstant) throws {
        let measured = try Self.run(id, pi: pi)
        let tight = pi == .abTEMFloat32
        let entries = Pins2.case(id)["detectors"] as! [String: [String: Any]]
        for (name, entry) in entries {
            let m = measured[name]!
            let scale = entry["scale"] as! Double
            switch entry["kind"] as! String {
            case "annular":
                #expect(scaledMax(m.values, Pins2.doubles("\(id).\(name)"), scale: scale) <= (tight ? 1e-10 : 1e-6),
                        "\(id) \(name)")
            case "segmented":
                #expect(scaledMax(m.values, Pins2.doubles("\(id).\(name).float64"), scale: scale) <= (tight ? 1e-10 : 1e-6),
                        "\(id) \(name) against the float64 sum")
                #expect(scaledMax(m.values, Pins2.doubles("\(id).\(name)"), scale: scale) <= 1e-6,
                        "\(id) \(name) against abTEM's float32 sum")
            default:
                let sums = (entry["sums"] as! [[Double]]).flatMap { $0 }
                let mine = (0..<(m.scanShape[0] * m.scanShape[1])).map { q in
                    m[q / m.scanShape[1], q % m.scanShape[1]].reduce(0, +)
                }
                #expect(scaledMax(mine, sums, scale: sums.map(abs).max()!) <= (tight ? 1e-10 : 1e-6), "\(id) \(name) sums")
                for pos in entry["stored_positions"] as! [[Int]] {
                    let pinned = Pins2.doubles("\(id).\(name).position.\(pos[0]).\(pos[1])")
                    #expect(relMax2(m[pos[0], pos[1]], pinned) <= (tight ? 1e-10 : 1e-6), "\(id) \(name) at \(pos)")
                }
            }
        }
    }

    @Test(arguments: [("C", 6), ("Si", 14), ("Cu", 29), ("Au", 79)])
    func singleAtomHAADFEqualsAbTEM(symbol: String, z: Int) throws {
        let s4 = Pins2.case("S4")
        let pinned = (s4["HAADF"] as! [String: Double])[symbol]!
        #expect(abs(try Self.singleAtomHAADF(z) - pinned) / pinned <= 1e-10, "\(symbol)")
    }

    static func singleAtomHAADF(_ z: Int, pi: LobatoParameterization.PiConstant = .abTEMFloat32) throws -> Double {
        let s = try AtomicStructure(box: SIMD3(8, 8, 2), atoms: [.init(atomicNumber: z, position: SIMD3(4, 4, 1))])
        let p = try ProjectedPotential(structure: s, nx: 128, ny: 128, sliceThickness: 2,
                                       parameterization: LobatoParameterization(piConstant: pi))
        let probe = try Probe(beam: ElectronBeam(energy: 200e3), grid: p.grid, aperture: Aperture(semiangleCutoff: 25))
        let exit = try Multislice.run(probe.build(at: SIMD2(4, 4)), through: p)[0]
        return try AnnularDetector(inner: 70, outer: 130).detect(DiffractionIntensity(exit))[0]
    }

    /// The parallel CPU batch equals one sequential run per position, bit for bit, and its time is printed.
    @Test func parallelCPUBatchEqualsSequentialRuns() throws {
        let potentials = try Pins2.potentials("S1", pi: .exact)
        let probe = try Pins2.probe("S1", grid: potentials[0].grid)
        let incidents = try Pins2.scan("S1").positions.map { try probe.build(at: $0) }
        let clock = ContinuousClock()
        var batch: [Waves] = []
        let parallel = try clock.measure { batch = try Multislice.runBatch(incidents, through: potentials[0]) }
        var single: [Waves] = []
        let sequential = try clock.measure {
            single = try incidents.map { try Multislice.run($0, through: potentials[0])[0] }
        }
        #expect(batch.count == single.count)
        for (a, b) in zip(batch, single) { #expect(a.field == b.field, "bit for bit") }
        print("CPU batch of \(incidents.count) positions: parallel \(parallel), sequential \(sequential)")
    }

    @Test(.enabled(if: Multislice.metalIsAvailable, "no Metal device: the batched Metal scan is not tested here"),
          arguments: ["S1", "S3"])
    func metalScanMatchesTheCPUScan(id: String) throws {
        let cpu = try Self.run(id, pi: .exact, backend: .cpu)
        let gpu = try Self.run(id, pi: .exact, backend: .metal)
        let entries = Pins2.case(id)["detectors"] as! [String: [String: Any]]
        for (name, entry) in entries {
            let scale = entry["scale"] as! Double
            let a = gpu[name]!.values, b = cpu[name]!.values
            let s = (entry["kind"] as! String) == "pixelated" ? max(b.map(abs).max()!, 1e-300) : scale
            let d = scaledMax(a, b, scale: s)
            print("Metal vs CPU scan, \(id) \(name): \(d)")
            #expect(d <= 2e-4, "\(id) \(name): Metal against CPU \(d)")
        }
    }

    /// registration-performance.md 1.6 and truth 3: the scan through unbuilt potentials (each slice built on the GPU,
    /// probes and intensities on the GPU) against the CPU scan, and two runs bit-equal.
    @Test(.enabled(if: Multislice.metalIsAvailable, "no Metal device: the GPU-built scan is not tested here"),
          arguments: ["S1", "S3"])
    func metalScanThroughUnbuiltPotentialsMatchesTheCPUScan(id: String) throws {
        let c = Pins2.case(id)
        let g = c["gpts"] as! [Int]
        let lazy = try Pins2.frozenPhonons(id).lazyPotentials(nx: g[0], ny: g[1],
                                                              sliceThickness: c["slice_thickness_requested"] as! Double)
        let probe = try Pins2.probe(id, grid: lazy[0].grid)
        let detectors = try Pins2.detectors()
        let scan = try Pins2.scan(id)
        let gpu = try STEM.scan(probe: probe, potentials: lazy, scan: scan, detectors: detectors.map(\.1), backend: .metal)
        let again = try STEM.scan(probe: probe, potentials: lazy, scan: scan, detectors: detectors.map(\.1), backend: .metal)
        let cpu = try Self.run(id, pi: .exact, backend: .cpu)
        let entries = c["detectors"] as! [String: [String: Any]]
        for ((name, _), (m, m2)) in zip(detectors, zip(gpu, again)) {
            #expect(m.values == m2.values, "\(id) \(name): deterministic")
            let entry = entries[name]!
            let b = cpu[name]!.values
            let s = (entry["kind"] as! String) == "pixelated" ? max(b.map(abs).max()!, 1e-300) : entry["scale"] as! Double
            let d = scaledMax(m.values, b, scale: s)
            print("GPU-built scan vs CPU, \(id) \(name): \(d)")
            #expect(d <= 2e-4, "\(id) \(name): \(d)")
        }
    }

    /// registration-generality.md: a CPU scan in chunks of one position equals the default bit for bit.
    @Test func theScanChunkIsAParameter() throws {
        let potentials = try Pins2.potentials("S1", pi: .exact)
        let probe = try Pins2.probe("S1", grid: potentials[0].grid)
        let scan = try Pins2.scan("S1")
        let d = try Pins2.detectors().map(\.1)
        let a = try STEM.scan(probe: probe, potentials: potentials, scan: scan, detectors: d)
        let b = try STEM.scan(probe: probe, potentials: potentials, scan: scan, detectors: d,
                              options: ExecutionOptions(scanChunkBytes: 1))
        #expect(zip(a, b).allSatisfy { $0.values == $1.values })
    }

    /// Performance amendment 3, item 6, accepted 2026-10-15 (amendment 4): the Gate B refuter's weak signal, one H atom
    /// seen by an 85–100 mrad annulus, whose image maximum is 9.5e−7 of the probe's intensity; Metal against CPU is about
    /// 2.4e−10 absolute, beyond 2e−4 × the image maximum (1.9e−10) and within the floored bound (2e−9).
    @Test(.enabled(if: Multislice.metalIsAvailable, "no Metal device: the floored detector bound is not tested here"))
    func weakHighAngleSignalIsWithinTheFlooredBound() throws {
        let s = try AtomicStructure(box: SIMD3(16, 16, 40), atoms: [.init(atomicNumber: 1, position: SIMD3(8, 8, 2))])
        let p = try ProjectedPotential(structure: s, nx: 256, ny: 256, sliceThickness: 0.4)
        let probe = try Probe(beam: ElectronBeam(energy: 300e3), grid: p.grid, aperture: Aperture(semiangleCutoff: 20))
        let scan = try GridScan(start: SIMD2(7.7, 7.7), end: SIMD2(8.3, 8.3), nx: 3, ny: 3)
        let d: [any Detector] = [try AnnularDetector(inner: 85, outer: 100)]
        let cpu = try STEM.scan(probe: probe, potentials: [p], scan: scan, detectors: d)[0].values
        let gpu = try STEM.scan(probe: probe, potentials: [p], scan: scan, detectors: d, backend: .metal)[0].values
        let maximum = cpu.map(abs).max()!
        let error = zip(gpu, cpu).map { abs($0 - $1) }.max()!
        let bound = MetalAccuracy.registered.detectorBound(imageMaximum: maximum, incidentIntensity: 1)
        print("weak signal: image max \(maximum), Metal error \(error), relative \(error / maximum), floored bound \(bound)")
        #expect(maximum < 1e-5, "the case is below the floor")
        #expect(error <= bound)
    }
}

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
        throws -> [String: Measurement] {
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
}

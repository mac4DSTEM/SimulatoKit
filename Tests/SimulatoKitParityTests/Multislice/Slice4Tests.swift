//
//  Slice4Tests.swift
//  Role: any grid size, `Potential(sampling:)` and beam tilt against abTEM 1.0.10 (tools/slice-4/pins.py). abTEM states
//        no tolerance for these paths, so the bound is assert_allclose's float64 default, rtol 1e-7 of the wave's
//        maximum; the Metal checks use the kit's registered float32 bound.
//

import SimulatoKitTestSupport
import CryptoKit
import Foundation
import Testing
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe

enum Slice4Pins {
    static let directory = URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent().appendingPathComponent("Fixtures")

    struct Loaded: @unchecked Sendable {
        let json: [String: Any]
        let blob: Data
    }

    static let loaded: Loaded = {
        let json = try! JSONSerialization.jsonObject(
            with: Data(contentsOf: directory.appendingPathComponent("slice-4-pins.json"))) as! [String: Any]
        let arrays = json["arrays"] as! [String: Any]
        let blob = try! Data(contentsOf: directory.appendingPathComponent(arrays["file"] as! String))
        let digest = SHA256.hash(data: blob).map { String(format: "%02x", $0) }.joined()
        precondition(digest == arrays["sha256"] as! String, "slice-4-arrays.bin does not match its sha256")
        return Loaded(json: json, blob: blob)
    }()

    static var json: [String: Any] { loaded.json }

    static func complex(_ name: String) -> ComplexField {
        let index = (json["arrays"] as! [String: Any])["index"] as! [String: Any]
        let entry = index[name] as! [String: Any]
        let offset = entry["offset"] as! Int, count = (entry["count"] as! Int) * 2
        let shape = entry["shape"] as! [Int]
        let d: [Double] = loaded.blob.withUnsafeBytes { raw in
            (0..<count).map { Double(bitPattern: UInt64(littleEndian: raw.loadUnaligned(fromByteOffset: offset + 8 * $0, as: UInt64.self))) }
        }
        return ComplexField(nx: shape[0], ny: shape[1], real: stride(from: 0, to: count, by: 2).map { d[$0] },
                            imaginary: stride(from: 1, to: count, by: 2).map { d[$0] })
    }

    static func vec(_ any: Any) -> [Double] { (any as! [NSNumber]).map(\.doubleValue) }

    static var structure: AtomicStructure {
        let b = vec(json["box"]!)
        let atoms = (json["atoms"] as! [[Any]]).map { row -> AtomicStructure.Atom in
            let p = vec(row[1])
            return .init(atomicNumber: (row[0] as! NSNumber).intValue, position: SIMD3(p[0], p[1], p[2]))
        }
        return try! AtomicStructure(box: SIMD3(b[0], b[1], b[2]), atoms: atoms)
    }

    static var beam: ElectronBeam { try! ElectronBeam(energy: (json["energy"] as! NSNumber).doubleValue) }
    static var sliceThickness: Double { (json["sliceThickness"] as! NSNumber).doubleValue }
    static var gpts: (Int, Int) { let g = json["gpts"] as! [Int]; return (g[0], g[1]) }
    static func tilt(_ run: String) -> SIMD2<Double> {
        let t = vec(((json["runs"] as! [String: Any])[run] as! [String: Any])["tilt"]!)
        return SIMD2(t[0], t[1])
    }
}

@Suite("Slice 4: any grid size, sampling, tilt (abTEM 1.0.10)")
struct Slice4Tests {
    static let rtol = 1e-7  // numpy.testing.assert_allclose's default for float64

    @Test func samplingGivesAbTEMsGrid() throws {
        for case let g as [String: Any] in Slice4Pins.json["grids"] as! [Any] {
            let s = g["sampling"]!
            let (sx, sy) = (s as? [NSNumber]).map { ($0[0].doubleValue, $0[1].doubleValue) }
                ?? ((s as! NSNumber).doubleValue, (s as! NSNumber).doubleValue)
            let want = g["gpts"] as! [Int]
            let potential = try Potential(structure: Slice4Pins.structure, sampling: (x: sx, y: sy),
                                          sliceThickness: Slice4Pins.sliceThickness)
            #expect([potential.grid.nx, potential.grid.ny] == want, "sampling (\(sx), \(sy))")
        }
        let scalar = try Potential(structure: Slice4Pins.structure, sampling: 0.1, sliceThickness: 2.025)
        #expect(scalar.grid.nx == 81 && scalar.grid.ny == 81)
    }

    @Test func samplingRefusesBadValues() {
        for bad in [0.0, -0.1, .nan, .infinity] {
            #expect(throws: SimulatoError.self) { _ = try Potential(structure: Slice4Pins.structure, sampling: bad, sliceThickness: 2.025) }
        }
    }

    static func potential() throws -> Potential {
        let (nx, ny) = Slice4Pins.gpts
        return try Potential(structure: Slice4Pins.structure, nx: nx, ny: ny, sliceThickness: Slice4Pins.sliceThickness)
    }

    static func incident(_ run: String, _ potential: Potential) throws -> Waves {
        let tilt = Slice4Pins.tilt(run)
        if run.hasPrefix("plane") { return Waves.planeWave(beam: Slice4Pins.beam, grid: potential.grid, tilt: tilt) }
        let probe = try Probe(beam: Slice4Pins.beam, grid: potential.grid, aperture: Aperture(semiangleCutoff: 20), tilt: tilt)
        let p = Slice4Pins.vec(Slice4Pins.json["probePosition"]!)
        return try probe.build(at: SIMD2(p[0], p[1]))
    }

    @Test(arguments: ["plane0", "planeT", "probe0", "probeT"])
    func exitWavesOnAnOddGridMatchAbTEM(run: String) throws {
        let potential = try Self.potential()
        #expect(potential.grid.nx == 90 && potential.grid.ny == 75, "90 x 75 has no vDSP length on the y axis")
        let exit = try Multislice.run(Self.incident(run, potential), through: potential)[0]
        #expect(exit.tilt == Slice4Pins.tilt(run), "the exit wave keeps the tilt")
        #expect(relMax(exit.field, Slice4Pins.complex(run)) <= Self.rtol, "\(run)")
    }

    @Test func tiltChangesTheExitWave() throws {
        let potential = try Self.potential()
        let a = try Multislice.run(Self.incident("planeT", potential), through: potential)[0].field
        let b = try Multislice.run(Self.incident("plane0", potential), through: potential)[0].field
        #expect(relMax(a, b) > 1e-3, "a tilt of 1.5 mrad moves the exit wave")
    }

    @Test(.enabled(if: Multislice.metalIsAvailable, "no Metal device"), arguments: ["plane0", "planeT", "probeT"])
    func metalMatchesCPUOnAnOddGridWithTilt(run: String) throws {
        let potential = try Self.potential()
        let incident = try Self.incident(run, potential)
        let cpu = try Multislice.run(incident, through: potential)[0]
        let gpu = try Multislice.run(incident, through: potential, backend: .metal)[0]
        #expect(gpu.tilt == cpu.tilt)
        #expect(relMax(gpu.field, cpu.field) <= MetalAccuracy.registered.wave, "\(run)")
    }

    @Test(.enabled(if: Multislice.metalIsAvailable, "no Metal device"))
    func aMetalBatchRefusesMixedTilts() throws {
        let potential = try Self.potential()
        let waves = [try Self.incident("planeT", potential), try Self.incident("plane0", potential)]
        #expect(throws: SimulatoError.self) { _ = try Multislice.runBatch(waves, through: potential, backend: .metal) }
        let cpu = try Multislice.runBatch(waves, through: potential)
        #expect(cpu.map(\.tilt) == waves.map(\.tilt))
    }

    @Test func probeRefusesANonFiniteTilt() throws {
        let grid = try Grid(nx: 32, ny: 32, extentX: 5, extentY: 5)
        #expect(throws: SimulatoError.self) {
            _ = try Probe(beam: Slice4Pins.beam, grid: grid, aperture: Aperture(semiangleCutoff: 20), tilt: SIMD2(.nan, 0))
        }
    }
}

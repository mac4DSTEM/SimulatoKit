//
//  WavesDetectorTests.swift
//  Role: WavesDetector over a scan against abTEM 1.0.10 (tools/waves/pins.py): the exit waves of a 3 x 2 scan with and
//        without `gpts`, at assert_allclose's float64 defaults (rtol 1e-7, atol 1e-12 of the wave maximum); `.metal`
//        against `.cpu` at MetalAccuracy.registered.wave; the typed refusals.
//

import SimulatoKitTestSupport
import CryptoKit
import Foundation
import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoMultislice
import SimulatoProbe

enum WavesPins {
    static let directory = URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent().appendingPathComponent("Fixtures")

    struct Loaded: @unchecked Sendable {
        let json: [String: Any]
        let blob: Data
    }

    static let loaded: Loaded = {
        let json = try! JSONSerialization.jsonObject(
            with: Data(contentsOf: directory.appendingPathComponent("waves-pins.json"))) as! [String: Any]
        let arrays = json["arrays"] as! [String: Any]
        let blob = try! Data(contentsOf: directory.appendingPathComponent(arrays["file"] as! String))
        let digest = SHA256.hash(data: blob).map { String(format: "%02x", $0) }.joined()
        precondition(digest == arrays["sha256"] as! String, "waves-arrays.bin does not match its sha256")
        return Loaded(json: json, blob: blob)
    }()

    static var json: [String: Any] { loaded.json }

    /// Interleaved (re, im) doubles of a complex array.
    static func doubles(_ name: String) -> [Double] {
        let index = (json["arrays"] as! [String: Any])["index"] as! [String: Any]
        let entry = index[name] as! [String: Any]
        let offset = entry["offset"] as! Int
        let count = (entry["count"] as! Int) * 2
        return loaded.blob.withUnsafeBytes { raw in
            (0..<count).map { Double(bitPattern: UInt64(littleEndian: raw.loadUnaligned(fromByteOffset: offset + 8 * $0, as: UInt64.self))) }
        }
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

    static func scan() throws -> GridScan {
        let s = json["scan"] as! [String: Any]
        let a = vec(s["start"]!), b = vec(s["end"]!), g = s["gpts"] as! [Int]
        return try GridScan(start: SIMD2(a[0], a[1]), end: SIMD2(b[0], b[1]), nx: g[0], ny: g[1])
    }

    static func setup() throws -> (probe: Probe, potential: Potential) {
        let n = json["gpts"] as! Int
        let potential = try Potential(structure: structure, nx: n, ny: n, sliceThickness: (json["sliceThickness"] as! NSNumber).doubleValue,
                                    parameterization: LobatoParameterization(piConstant: .abTEMFloat32))
        let beam = try ElectronBeam(energy: (json["energy"] as! NSNumber).doubleValue)
        let probe = try Probe(beam: beam, grid: potential.grid,
                              aperture: Aperture(semiangleCutoff: (json["semiangle"] as! NSNumber).doubleValue))
        return (probe, potential)
    }

    /// assert_allclose: |a − b| ≤ atol + rtol |b| per number, rtol 1e-7, atol 1e-12 of the maximum of |b|; returns the
    /// largest |a − b| / (atol + rtol |b|) (≤ 1 passes) and the largest |a − b| / max|b|.
    static func allclose(_ a: ScanWaves, _ b: [Double]) -> (ratio: Double, relMax: Double) {
        var maxAbs = 0.0
        for k in 0..<(b.count / 2) { maxAbs = max(maxAbs, hypot(b[2 * k], b[2 * k + 1])) }
        let atol = 1e-12 * maxAbs
        var ratio = 0.0, worst = 0.0
        for k in 0..<a.real.count {
            for (x, y) in [(a.real[k], b[2 * k]), (a.imaginary[k], b[2 * k + 1])] {
                ratio = max(ratio, abs(x - y) / (atol + 1e-7 * abs(y)))
                worst = max(worst, abs(x - y) / maxAbs)
            }
        }
        return (ratio, worst)
    }
}

@Suite("WavesDetector against abTEM 1.0.10")
struct WavesDetectorTests {
    static var runs: [[String: Any]] { WavesPins.json["runs"] as! [[String: Any]] }

    @Test(arguments: [0, 1, 2, 3, 4, 5])
    func exitWavesMatchAbTEM(run: Int) throws {
        let entry = Self.runs[run]
        let (probe, potential) = try WavesPins.setup()
        let detector: WavesDetector
        if let g = entry["gpts"] as? [Int] { detector = try WavesDetector(rows: g[0], columns: g[1]) } else { detector = WavesDetector() }
        let out = try STEM.scan(probe: probe, potentials: [potential], scan: WavesPins.scan(), detector: detector)
        let shape = entry["shape"] as! [Int]
        #expect([out.scanShape[0], out.scanShape[1], out.rows, out.columns] == shape)
        let c = WavesPins.allclose(out, WavesPins.doubles(entry["name"] as! String))
        print("WavesDetector \(entry["name"]!): max |a-b|/(atol+rtol|b|) = \(c.ratio), max|a-b|/max|b| = \(c.relMax)")
        #expect(c.ratio <= 1, "\(entry["name"]!): outside assert_allclose rtol 1e-7, atol 1e-12 of the maximum")
    }

    @Test func projectedPotentialPathEqualsThePotentialPath() throws {
        let (probe, potential) = try WavesPins.setup()
        let detector = try WavesDetector(rows: 33, columns: 40)
        let a = try STEM.scan(probe: probe, potentials: [potential], scan: WavesPins.scan(), detector: detector)
        let b = try STEM.scan(probe: probe, potentials: [potential.build()], scan: WavesPins.scan(), detector: detector)
        #expect(a.real == b.real && a.imaginary == b.imaginary)
        let w = a[2, 1]
        #expect(w.grid.nx == 33 && w.grid.ny == 40)
    }

    @Test func memoryRefusalIsTyped() throws {
        let (probe, potential) = try WavesPins.setup()
        let scan = try WavesPins.scan()
        // 6 positions x 48 x 48 x 16 bytes = 221 184 bytes.
        let need = WavesDetector.requiredBytes(positions: 6, rows: 48, columns: 48)
        #expect(need == 221_184)
        #expect(throws: WavesDetectorError.tooLarge(positions: 6, rows: 48, columns: 48, required: need, limit: need - 1)) {
            try STEM.scan(probe: probe, potentials: [potential], scan: scan, detector: WavesDetector(), memoryLimit: need - 1)
        }
        // The resampled size counts, not the simulation grid.
        let small = try WavesDetector(rows: 10, columns: 11)
        let ok = try STEM.scan(probe: probe, potentials: [potential], scan: scan, detector: small, memoryLimit: 6 * 10 * 11 * 16)
        #expect(ok.rows == 10 && ok.columns == 11)
        #expect(throws: WavesDetectorError.self) {
            try STEM.scan(probe: probe, potentials: [potential.build()], scan: scan, detector: small, memoryLimit: 6 * 10 * 11 * 16 - 1)
        }
        // Default: half of the physical memory; a huge scan is refused before any work.
        let huge = try GridScan(start: SIMD2(0, 0), end: SIMD2(4, 4), nx: 100_000, ny: 100_000)
        #expect(throws: WavesDetectorError.self) {
            try STEM.scan(probe: probe, potentials: [potential], scan: huge, detector: WavesDetector())
        }
        #expect(WavesDetector.requiredBytes(positions: .max, rows: 2, columns: 2) == .max)
        #expect("\(WavesDetectorError.tooLarge(positions: 1, rows: 2, columns: 3, required: 96, limit: 5))".contains("96 bytes"))
    }

    @Test func frozenPhononsAreRefused() throws {
        let (probe, potential) = try WavesPins.setup()
        #expect(throws: WavesDetectorError.severalPotentials(2)) {
            try STEM.scan(probe: probe, potentials: [potential, potential], scan: WavesPins.scan(), detector: WavesDetector())
        }
        #expect(throws: WavesDetectorError.severalPotentials(0)) {
            try STEM.scan(probe: probe, potentials: [Potential](), scan: WavesPins.scan(), detector: WavesDetector())
        }
        #expect(throws: (any Error).self) { try WavesDetector(rows: 1, columns: 8) }
    }

    @Test(.enabled(if: Multislice.metalIsAvailable, "no Metal device: the Metal path is not tested here"),
          arguments: [0, 2, 3])
    func metalMatchesCPU(run: Int) throws {
        let entry = Self.runs[run]
        let (probe, potential) = try WavesPins.setup()
        let detector: WavesDetector
        if let g = entry["gpts"] as? [Int] { detector = try WavesDetector(rows: g[0], columns: g[1]) } else { detector = WavesDetector() }
        let cpu = try STEM.scan(probe: probe, potentials: [potential], scan: WavesPins.scan(), detector: detector)
        let gpu = try STEM.scan(probe: probe, potentials: [potential], scan: WavesPins.scan(), detector: detector, backend: .metal)
        let d = relMax(gpu.real + gpu.imaginary, cpu.real + cpu.imaginary)
        // The complex bound of max |difference| / max |wave|: the real/imaginary max is within a factor of 2 of it.
        var worst = 0.0, top = 0.0
        for k in cpu.real.indices {
            worst = max(worst, hypot(gpu.real[k] - cpu.real[k], gpu.imaginary[k] - cpu.imaginary[k]))
            top = max(top, hypot(cpu.real[k], cpu.imaginary[k]))
        }
        print("WavesDetector \(entry["name"]!) Metal vs CPU: \(worst / top) (registered \(MetalAccuracy.registered.wave)); parts \(d)")
        #expect(worst / top <= MetalAccuracy.registered.wave)
    }
}

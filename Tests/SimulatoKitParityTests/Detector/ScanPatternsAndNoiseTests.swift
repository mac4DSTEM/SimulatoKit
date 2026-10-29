//
//  ScanPatternsAndNoiseTests.swift
//  Role: GridScan(sampling:), the centre-of-mass detector, the Gaussian source size and Poisson noise against abTEM 1.0.10
//        (tools/pins/scan-patterns/pins.py). The centre of mass is abTEM's complex64, so its bound is 1e-5 of the largest value;
//        the Gaussian filter is float64 (assert_allclose's default, rtol 1e-7). Poisson noise is statistical.
//

import SimulatoKitTestSupport
import CryptoKit
import Foundation
import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoMultislice
import SimulatoProbe

enum Slice4ScanPins {
    static let directory = URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent().appendingPathComponent("Fixtures")

    struct Loaded: @unchecked Sendable {
        let json: [String: Any]
        let blob: Data
    }

    static let loaded: Loaded = {
        let json = try! JSONSerialization.jsonObject(
            with: Data(contentsOf: directory.appendingPathComponent("slice-4-scan-pins.json"))) as! [String: Any]
        let arrays = json["arrays"] as! [String: Any]
        let blob = try! Data(contentsOf: directory.appendingPathComponent(arrays["file"] as! String))
        let digest = SHA256.hash(data: blob).map { String(format: "%02x", $0) }.joined()
        precondition(digest == arrays["sha256"] as! String, "slice-4-scan-arrays.bin does not match its sha256")
        return Loaded(json: json, blob: blob)
    }()

    static var json: [String: Any] { loaded.json }

    /// Doubles of an array; a complex array comes as interleaved (re, im) pairs.
    static func doubles(_ name: String) -> [Double] {
        let index = (json["arrays"] as! [String: Any])["index"] as! [String: Any]
        let entry = index[name] as! [String: Any]
        let offset = entry["offset"] as! Int
        let count = (entry["count"] as! Int) * ((entry["kind"] as! String) == "c16" ? 2 : 1)
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

    /// The scan the pins ran, with the three detectors of interest, on `.cpu`.
    static let measurements: [ScanMeasurement] = {
        let n = json["gpts"] as! Int
        let potential = try! Potential(structure: structure, nx: n, ny: n, sliceThickness: (json["sliceThickness"] as! NSNumber).doubleValue)
        let beam = try! ElectronBeam(energy: (json["energy"] as! NSNumber).doubleValue)
        let probe = try! Probe(beam: beam, grid: potential.grid,
                               aperture: Aperture(semiangleCutoff: (json["semiangle"] as! NSNumber).doubleValue))
        return try! STEM.scan(probe: probe, potentials: [potential], scan: scan(),
                              detectors: [PixelatedDetector(maxAngle: .cutoff),
                                          CenterOfMassDetector(maxAngle: .cutoff, units: .milliradians),
                                          CenterOfMassDetector(maxAngle: .cutoff, units: .inverseAngstroms)])
    }()

    static func relMax(_ a: [Double], _ b: [Double]) -> Double {
        precondition(a.count == b.count, "\(a.count) vs \(b.count) values")
        var d = 0.0, m = 0.0
        for i in a.indices { d = max(d, abs(a[i] - b[i])); m = max(m, abs(b[i])) }
        return m > 0 ? d / m : d
    }
}

@Suite("Scan sampling, centre of mass, source size, noise (abTEM 1.0.10)")
struct ScanPatternsAndNoiseTests {
    @Test func gridScanSamplingGivesAbTEMsPoints() throws {
        for case let s as [String: Any] in Slice4ScanPins.json["scans"] as! [Any] {
            let a = Slice4ScanPins.vec(s["start"]!), b = Slice4ScanPins.vec(s["end"]!)
            let (sx, sy) = (s["sampling"] as? [NSNumber]).map { ($0[0].doubleValue, $0[1].doubleValue) }
                ?? ((s["sampling"] as! NSNumber).doubleValue, (s["sampling"] as! NSNumber).doubleValue)
            let scan = try GridScan(start: SIMD2(a[0], a[1]), end: SIMD2(b[0], b[1]), samplingX: sx, samplingY: sy)
            #expect([scan.nx, scan.ny] == s["gpts"] as! [Int], "sampling (\(sx), \(sy))")
        }
        let withEnd = try GridScan(start: SIMD2(0, 0), end: SIMD2(10, 10), samplingX: 1, samplingY: 1, endpoint: true)
        #expect(withEnd.nx == 11 && withEnd.ny == 11, "abTEM adds one point with endpoint")
        #expect(throws: SimulatoError.self) {
            _ = try GridScan(start: SIMD2(2, 0), end: SIMD2(1, 1), samplingX: 1, samplingY: 1)
        }
        #expect(throws: SimulatoError.self) {
            _ = try GridScan(start: SIMD2(0, 0), end: SIMD2(1, 1), samplingX: 0, samplingY: 0.1)
        }
    }

    @Test func scanPatternsMatchAbTEM() {
        let dp = Slice4ScanPins.measurements[0]
        let shape = Slice4ScanPins.json["dpShape"] as! [Int]
        #expect(dp.scanShape == [6, 5])
        #expect(dp.valueShape == Array(shape.suffix(2)))
        #expect(Slice4ScanPins.relMax(dp.values, Slice4ScanPins.doubles("dp")) <= 1e-6, "the default Lobato parameterisation (.exact) is registered at 1e-6 of abTEM float32 π")
    }

    @Test func centreOfMassMatchesAbTEM() {
        for (index, name) in [(1, "com_mrad"), (2, "com_inv")] {
            let kit = Slice4ScanPins.measurements[index]
            #expect(kit.valueShape == [2])
            // abTEM's complex values are interleaved (re, im) = (x, y), the kit's pairs (x, y): the same order.
            let want = Slice4ScanPins.doubles(name)
            #expect(kit.values.count == want.count)
            #expect(Slice4ScanPins.relMax(kit.values, want) <= 1e-5, "\(name)")
        }
    }

    @Test func gaussianSourceSizeMatchesScipyViaAbTEM() throws {
        let dp = Slice4ScanPins.measurements[0]
        let sampling = Slice4ScanPins.vec(Slice4ScanPins.json["scanSampling"]!)
        let sigma = Slice4ScanPins.vec(Slice4ScanPins.json["sigma"]!)
        let a = try dp.gaussianSourceSize(sigmaX: sigma[0], sigmaY: sigma[1], samplingX: sampling[0], samplingY: sampling[1])
        #expect(Slice4ScanPins.relMax(a.values, Slice4ScanPins.doubles("dp_gauss")) <= 1e-7)
        // A σ larger than the scan wraps around it (σ / d = 2.4 pixels on a 5-point axis).
        let big = (Slice4ScanPins.json["sigmaBig"] as! NSNumber).doubleValue
        let b = try dp.gaussianSourceSize(sigma: big, samplingX: sampling[0], samplingY: sampling[1])
        #expect(Slice4ScanPins.relMax(b.values, Slice4ScanPins.doubles("dp_gauss_big")) <= 1e-7)
        // No source size leaves the scan alone, and the total is conserved.
        let zero = try dp.gaussianSourceSize(sigma: 0, samplingX: sampling[0], samplingY: sampling[1])
        #expect(zero.values == dp.values)
        let total: Double = dp.values.reduce(0, +)
        #expect(abs(a.values.reduce(0, +) - total) <= 1e-9 * total)
        #expect(throws: SimulatoError.self) { _ = try dp.gaussianSourceSize(sigma: -1, samplingX: 1, samplingY: 1) }
    }

    @Test func poissonNoiseHasThePoissonMeanAndVariance() throws {
        let n = 20000
        for mean in [0.5, 4.0, 9.9, 10.0, 37.0, 800.0] {
            let m = ScanMeasurement(scanShape: [n, 1], valueShape: [], values: [Double](repeating: mean / 3.0, count: n))
            let noisy = try m.poissonNoise(dose: 3.0, seed: 11)
            let mu = noisy.values.reduce(0, +) / Double(n)
            let variance = noisy.values.reduce(0) { $0 + ($1 - mu) * ($1 - mu) } / Double(n - 1)
            let meanBound: Double = 5 * (mean / Double(n)).squareRoot()
            let varianceBound: Double = 6 * mean * (2.0 / Double(n)).squareRoot() + 0.02 * mean
            #expect(abs(mu - mean) <= meanBound, "mean \(mean): got \(mu)")
            #expect(abs(variance - mean) <= varianceBound, "variance \(mean): got \(variance)")
            #expect(noisy.values.allSatisfy { $0 >= 0 && $0 == $0.rounded() })
        }
    }

    @Test func poissonNoiseIsSeededClippedAndRefusesBadDoses() throws {
        let m = ScanMeasurement(scanShape: [2, 2], valueShape: [], values: [1, -3, 0, 5])
        let a = try m.poissonNoise(dose: 2, seed: 7), b = try m.poissonNoise(dose: 2, seed: 7)
        let other = try m.poissonNoise(dose: 2, seed: 8)
        #expect(a.values == b.values)
        #expect(a.values != other.values)
        #expect(a.values[1] == 0 && a.values[2] == 0, "negative and zero values stay 0 electrons")
        #expect(try m.poissonNoise(dose: 0, seed: 1).values == [0, 0, 0, 0])
        for bad in [-1.0, .nan, .infinity] { #expect(throws: SimulatoError.self) { _ = try m.poissonNoise(dose: bad, seed: 1) } }
        let perArea = try m.poissonNoise(dosePerArea: 1, samplingX: 0.5, samplingY: 4, seed: 7)
        #expect(perArea.values == a.values, "dose per pattern = area per pixel (2 Å²) x 1")
    }
}

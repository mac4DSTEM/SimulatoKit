//
//  LinescanPins.swift
//  Role: the pins written by tools/pins/linescan/pins.py (abTEM 1.0.10, float64): the JSON and its float64 sidecar (checked
//        against the sha256 the generator recorded), and the comparison abTEM's own tests use, `assert_allclose`.
//

import CryptoKit
import Foundation
import SimulatoCore

enum LinescanPins {
    static let directory = URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent().appendingPathComponent("Fixtures")

    struct Loaded: @unchecked Sendable {  // immutable after load
        let json: [String: Any]
        let blob: Data
    }

    static let loaded: Loaded = {
        let json = try! JSONSerialization.jsonObject(
            with: Data(contentsOf: directory.appendingPathComponent("linescan-pins.json"))) as! [String: Any]
        let arrays = json["arrays"] as! [String: Any]
        let blob = try! Data(contentsOf: directory.appendingPathComponent(arrays["file"] as! String))
        let digest = SHA256.hash(data: blob).map { String(format: "%02x", $0) }.joined()
        precondition(digest == arrays["sha256"] as! String, "linescan-arrays.bin does not match its sha256 in the pins")
        return Loaded(json: json, blob: blob)
    }()

    static var json: [String: Any] { loaded.json }

    static func doubles(_ name: String) -> [Double] {
        let index = (json["arrays"] as! [String: Any])["index"] as! [String: Any]
        guard let entry = index[name] as? [String: Any] else { fatalError("no pinned array \(name)") }
        let offset = entry["offset"] as! Int
        let count = entry["count"] as! Int
        return loaded.blob.withUnsafeBytes { raw in
            (0..<count).map { Double(bitPattern: UInt64(littleEndian: raw.loadUnaligned(fromByteOffset: offset + 8 * $0, as: UInt64.self))) }
        }
    }

    static func vec(_ any: Any) -> [Double] { (any as! [NSNumber]).map(\.doubleValue) }
    static func num(_ any: Any?) -> Double { (any as! NSNumber).doubleValue }
    static func list(_ key: String) -> [[String: Any]] { json[key] as! [[String: Any]] }

    static var structure: AtomicStructure {
        let b = vec(json["box"]!)
        let atoms = (json["atoms"] as! [[Any]]).map { row -> AtomicStructure.Atom in
            let p = vec(row[1])
            return .init(atomicNumber: (row[0] as! NSNumber).intValue, position: SIMD3(p[0], p[1], p[2]))
        }
        return try! AtomicStructure(box: SIMD3(b[0], b[1], b[2]), atoms: atoms)
    }
}

/// `numpy.testing.assert_allclose(actual, desired, rtol, atol)`: the largest |a − d| / (atol + rtol |d|) (passes at 1 or
/// less) and the largest |a − d| / max |d|, which is the number reported as the measured agreement.
struct Closeness {
    let ratio: Double
    let relativeToMax: Double
    var passes: Bool { ratio <= 1 }
}

func allClose(_ a: [Double], _ d: [Double], rtol: Double = 1e-7, atolOfMax: Double = 1e-12) -> Closeness {
    guard a.count == d.count else { return Closeness(ratio: .infinity, relativeToMax: .infinity) }
    let m = d.map(abs).max() ?? 0
    let atol = atolOfMax * m
    var ratio = 0.0, diff = 0.0
    for i in a.indices {
        let e = abs(a[i] - d[i])
        diff = max(diff, e)
        ratio = max(ratio, e / (atol + rtol * abs(d[i])))
    }
    return Closeness(ratio: ratio, relativeToMax: m > 0 ? diff / m : diff)
}

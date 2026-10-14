//
//  Slice2Support.swift
//  Role: the slice-2 pins (tools/slice-2/pins.py; registration-slice-2.md section 4): the JSON and its sidecar (checked
//        against the recorded sha256), the cases' structures, displacements, probes and detectors.
//

import CryptoKit
import Foundation
import SimulatoCore
import SimulatoDetector
import SimulatoMultislice
import SimulatoProbe

enum Pins2 {
    static let directory = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures")

    struct Loaded: @unchecked Sendable {  // immutable after load
        let json: [String: Any]
        let blob: Data
    }

    static let loaded: Loaded = {
        let json = try! JSONSerialization.jsonObject(
            with: Data(contentsOf: directory.appendingPathComponent("slice-2-pins.json"))) as! [String: Any]
        let arrays = json["arrays"] as! [String: Any]
        let blob = try! Data(contentsOf: directory.appendingPathComponent(arrays["file"] as! String))
        let digest = SHA256.hash(data: blob).map { String(format: "%02x", $0) }.joined()
        precondition(digest == arrays["sha256"] as! String, "slice-2-arrays.bin does not match its sha256 in the pins")
        return Loaded(json: json, blob: blob)
    }()

    static var json: [String: Any] { loaded.json }

    static func doubles(_ name: String) -> [Double] {
        let index = (json["arrays"] as! [String: Any])["index"] as! [String: Any]
        guard let entry = index[name] as? [String: Any] else { fatalError("no pinned array \(name)") }
        let offset = entry["offset"] as! Int
        let count = (entry["count"] as! Int) * ((entry["kind"] as! String) == "c16" ? 2 : 1)
        return loaded.blob.withUnsafeBytes { raw in
            (0..<count).map { i in
                Double(bitPattern: UInt64(littleEndian: raw.loadUnaligned(fromByteOffset: offset + 8 * i, as: UInt64.self)))
            }
        }
    }

    static func `case`(_ id: String) -> [String: Any] { (json["cases"] as! [String: Any])[id] as! [String: Any] }

    static func structure(_ id: String) throws -> AtomicStructure {
        let c = Self.case(id)
        let box = c["box"] as! [Double]
        let numbers = c["numbers"] as! [Int]
        let positions = c["positions"] as! [[Double]]
        return try AtomicStructure(box: SIMD3(box[0], box[1], box[2]), atoms: zip(numbers, positions).map {
            AtomicStructure.Atom(atomicNumber: $0.0, position: SIMD3($0.1[0], $0.1[1], $0.1[2]))
        })
    }

    /// abTEM's displacements of every configuration, given to the port.
    static func frozenPhonons(_ id: String) throws -> FrozenPhonons {
        let s = try structure(id)
        let n = max(1, (Self.case(id)["seeds"] as! [Int]).count)
        let d = (0..<n).map { c -> [SIMD3<Double>] in
            let v = doubles("\(id).displacement.\(c)")
            return stride(from: 0, to: v.count, by: 3).map { SIMD3(v[$0], v[$0 + 1], v[$0 + 2]) }
        }
        return try FrozenPhonons(structure: s, displacements: d)
    }

    static func potentials(_ id: String, pi: LobatoParameterization.PiConstant) throws -> [ProjectedPotential] {
        let c = Self.case(id)
        let g = c["gpts"] as! [Int]
        return try frozenPhonons(id).potentials(nx: g[0], ny: g[1], sliceThickness: c["slice_thickness_requested"] as! Double,
                                                parameterization: LobatoParameterization(piConstant: pi))
    }

    static func probe(_ id: String, grid: Grid) throws -> Probe {
        let c = Self.case(id)
        return try Probe(beam: ElectronBeam(energy: c["energy"] as! Double), grid: grid,
                         aperture: Aperture(semiangleCutoff: c["cutoff_mrad"] as! Double))
    }

    static func scan(_ id: String) throws -> GridScan {
        let s = Self.case(id)["scan"] as! [String: Any]
        let a = s["start"] as! [Double], b = s["end"] as! [Double], g = s["gpts"] as! [Int]
        return try GridScan(start: SIMD2(a[0], a[1]), end: SIMD2(b[0], b[1]), nx: g[0], ny: g[1])
    }

    /// The registered detectors, by name, in the pins' order.
    static func detectors() throws -> [(String, any Detector)] {
        try (json["detectors"] as! [[Any]]).map { row in
            let name = row[0] as! String, kind = row[1] as! String, p = row[2] as! [String: Any]
            switch kind {
            case "annular":
                return (name, try AnnularDetector(inner: p["inner"] as! Double, outer: p["outer"] as? Double))
            case "pixelated":
                let m: PixelatedDetector.MaxAngle = ["cutoff": .cutoff, "valid": .valid, "full": .full][p["max_angle"] as! String]!
                return (name, PixelatedDetector(maxAngle: m))
            default:
                return (name, try SegmentedDetector(radialBins: p["nbins_radial"] as! Int,
                                                    azimuthalBins: p["nbins_azimuthal"] as! Int,
                                                    inner: p["inner"] as! Double, outer: p["outer"] as! Double,
                                                    rotation: p["rotation"] as! Double))
            }
        }
    }
}

/// max |a − b| / scale.
/// A size mismatch is reported as +infinity, so the expectation fails (rather than the test process crashing).
func scaledMax(_ a: [Double], _ b: [Double], scale: Double) -> Double {
    guard a.count == b.count else { return .infinity }
    var d = 0.0
    for i in a.indices { d = max(d, abs(a[i] - b[i])) }
    return d / scale
}

func relMax2(_ a: [Double], _ b: [Double]) -> Double {
    guard a.count == b.count else { return .infinity }
    return scaledMax(a, b, scale: max(b.map(abs).max() ?? 0, 1e-300))
}

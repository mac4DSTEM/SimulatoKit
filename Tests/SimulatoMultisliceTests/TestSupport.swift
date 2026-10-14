//
//  TestSupport.swift
//  Role: the pins written by tools/slice-1/pins.py (registration section 4): the JSON and its float64 sidecar (checked
//        against the sha256 the generator recorded), the cases' structures and beams, and the comparison helpers.
//

import CryptoKit
import Foundation
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe

enum Pins {
    static let directory = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures")

    struct Loaded: @unchecked Sendable {  // immutable after load
        let json: [String: Any]
        let blob: Data
    }

    static let loaded: Loaded = {
        let json = try! JSONSerialization.jsonObject(
            with: Data(contentsOf: directory.appendingPathComponent("slice-1-pins.json"))) as! [String: Any]
        let arrays = json["arrays"] as! [String: Any]
        let blob = try! Data(contentsOf: directory.appendingPathComponent(arrays["file"] as! String))
        let digest = SHA256.hash(data: blob).map { String(format: "%02x", $0) }.joined()
        precondition(digest == arrays["sha256"] as! String, "slice-1-arrays.bin does not match its sha256 in the pins")
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
                Double(bitPattern: UInt64(littleEndian: raw.loadUnaligned(fromByteOffset: offset + 8 * i,
                                                                          as: UInt64.self)))
            }
        }
    }

    static func shape(_ name: String) -> [Int] {
        let index = (json["arrays"] as! [String: Any])["index"] as! [String: Any]
        return (index[name] as! [String: Any])["shape"] as! [Int]
    }

    static func complex(_ name: String) -> ComplexField {
        let d = doubles(name), s = shape(name)
        return ComplexField(nx: s[0], ny: s[1], real: stride(from: 0, to: d.count, by: 2).map { d[$0] },
                            imaginary: stride(from: 1, to: d.count, by: 2).map { d[$0] })
    }

    static func hasArray(_ name: String) -> Bool {
        (((json["arrays"] as! [String: Any])["index"] as! [String: Any])[name]) != nil
    }

    static func `case`(_ id: String) -> [String: Any] {
        (json["cases"] as! [String: Any])[id] as! [String: Any]
    }

    static func structure(_ id: String) throws -> AtomicStructure {
        let c = Self.case(id)
        let box = c["box"] as! [Double]
        let numbers = c["numbers"] as! [Int]
        let positions = c["positions"] as! [[Double]]
        return try AtomicStructure(box: SIMD3(box[0], box[1], box[2]),
                                   atoms: zip(numbers, positions).map {
                                       AtomicStructure.Atom(atomicNumber: $0.0, position: SIMD3($0.1[0], $0.1[1], $0.1[2]))
                                   })
    }

    static func potential(_ id: String, pi: LobatoParameterization.PiConstant) throws -> ProjectedPotential {
        let c = Self.case(id)
        let g = c["gpts"] as! [Int]
        return try ProjectedPotential(structure: structure(id), nx: g[0], ny: g[1],
                                      sliceThickness: c["slice_thickness_requested"] as! Double,
                                      parameterization: LobatoParameterization(piConstant: pi))
    }

    static func beam(_ id: String) throws -> ElectronBeam {
        try ElectronBeam(energy: Self.case(id)["energy"] as! Double)
    }

    static func aberrations(_ dict: [String: Any]) -> Aberrations {
        var a = Aberrations()
        for (k, v) in dict {
            let x = v as! Double
            switch k {
            case "C10": a.C10 = x
            case "C12": a.C12 = x
            case "phi12": a.phi12 = x
            case "C30": a.C30 = x
            default: fatalError("aberration \(k) not mapped in the tests")
            }
        }
        return a
    }

    /// The incident wave of a beam of a case (plane wave or probe).
    static func incident(_ id: String, _ label: String, grid: Grid) throws -> Waves {
        let b = (Self.case(id)["beams"] as! [String: Any])[label] as! [String: Any]
        let beam = try beam(id)
        if (b["kind"] as! String) == "plane" { return Waves.planeWave(beam: beam, grid: grid) }
        let p = b["position"] as! [Double]
        let probe = try Probe(beam: beam, grid: grid, aperture: Aperture(semiangleCutoff: b["cutoff_mrad"] as! Double),
                              aberrations: aberrations(b["aberrations"] as! [String: Any]))
        return try probe.build(at: SIMD2(p[0], p[1]))
    }
}

/// max |a − b| / max |b| (0 when both are zero).
func relMax(_ a: [Double], _ b: [Double]) -> Double {
    precondition(a.count == b.count, "relMax: \(a.count) vs \(b.count) values")
    var d = 0.0, m = 0.0
    for i in a.indices {
        d = max(d, abs(a[i] - b[i]))
        m = max(m, abs(b[i]))
    }
    return m > 0 ? d / m : d
}

/// max |a − b| / max |b| for complex fields.
func relMax(_ a: ComplexField, _ b: ComplexField) -> Double {
    precondition(a.nx == b.nx && a.ny == b.ny)
    var d = 0.0, m = 0.0
    for i in a.real.indices {
        let dr = a.real[i] - b.real[i], di = a.imaginary[i] - b.imaginary[i]
        d = max(d, (dr * dr + di * di).squareRoot())
        m = max(m, (b.real[i] * b.real[i] + b.imaginary[i] * b.imaginary[i]).squareRoot())
    }
    return m > 0 ? d / m : d
}

func relDiff(_ a: Double, _ b: Double) -> Double { b == 0 ? abs(a) : abs(a - b) / abs(b) }

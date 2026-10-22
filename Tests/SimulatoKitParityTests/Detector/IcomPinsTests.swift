//
//  IcomPinsTests.swift
//  Role: `RealField.integratingGradient`, the integrated centre of mass of a scan, `integrateDisc` and `tiled` against abTEM
//        1.0.10 (tools/icom/pins.py). abTEM states no tolerance: assert_allclose's float64 default (rtol 1e-7) with atol
//        1e-12 of the largest value. The scan's own centre of mass is held at the 1e-5 the centre-of-mass detector is
//        registered at (Slice4ScanTests), carried through the integration.
//

import CryptoKit
import Foundation
import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoMultislice
import SimulatoProbe

enum IcomPins {
    static let directory = URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent().appendingPathComponent("Fixtures")

    struct Loaded: @unchecked Sendable {
        let json: [String: Any]
        let blob: Data
    }

    static let loaded: Loaded = {
        let json = try! JSONSerialization.jsonObject(
            with: Data(contentsOf: directory.appendingPathComponent("icom-pins.json"))) as! [String: Any]
        let arrays = json["arrays"] as! [String: Any]
        let blob = try! Data(contentsOf: directory.appendingPathComponent(arrays["file"] as! String))
        let digest = SHA256.hash(data: blob).map { String(format: "%02x", $0) }.joined()
        precondition(digest == arrays["sha256"] as! String, "icom-arrays.bin does not match its sha256")
        return Loaded(json: json, blob: blob)
    }()

    static var json: [String: Any] { loaded.json }

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

    /// The complex array's components as two real fields.
    static func gradient(_ name: String, nx: Int, ny: Int, sx: Double, sy: Double) -> (x: RealField, y: RealField) {
        let d = doubles(name)
        let gx = stride(from: 0, to: d.count, by: 2).map { d[$0] }, gy = stride(from: 1, to: d.count, by: 2).map { d[$0] }
        return (RealField(nx: nx, ny: ny, values: gx, samplingX: sx, samplingY: sy),
                RealField(nx: nx, ny: ny, values: gy, samplingX: sx, samplingY: sy))
    }

    /// assert_allclose(actual, desired, rtol: 1e-7, atol: 1e-12 · max|desired|): the worst excess over the bound, ≤ 0 passes.
    static func allcloseExcess(_ a: [Double], _ b: [Double], rtol: Double = 1e-7) -> Double {
        precondition(a.count == b.count)
        let atol = 1e-12 * (b.map(abs).max() ?? 0)
        var worst = -Double.infinity
        for i in a.indices { worst = max(worst, abs(a[i] - b[i]) - (atol + rtol * abs(b[i]))) }
        return worst
    }

    /// max |a − b| / max |b|
    static func relMax(_ a: [Double], _ b: [Double]) -> Double {
        var d = 0.0, m = 0.0
        for i in a.indices { d = max(d, abs(a[i] - b[i])); m = max(m, abs(b[i])) }
        return m > 0 ? d / m : d
    }

    static let scan: [ScanMeasurement] = {
        let s = json["scan"] as! [String: Any]
        let b = vec(s["box"]!)
        let atoms = (s["atoms"] as! [[Any]]).map { row -> AtomicStructure.Atom in
            let p = vec(row[1])
            return .init(atomicNumber: (row[0] as! NSNumber).intValue, position: SIMD3(p[0], p[1], p[2]))
        }
        let structure = try! AtomicStructure(box: SIMD3(b[0], b[1], b[2]), atoms: atoms)
        let n = s["gpts"] as! Int
        let potential = try! Potential(structure: structure, nx: n, ny: n, sliceThickness: (s["sliceThickness"] as! NSNumber).doubleValue)
        let beam = try! ElectronBeam(energy: (s["energy"] as! NSNumber).doubleValue)
        let probe = try! Probe(beam: beam, grid: potential.grid,
                               aperture: Aperture(semiangleCutoff: (s["semiangle"] as! NSNumber).doubleValue))
        let a = vec(s["start"]!), e = vec(s["end"]!), g = s["scanGpts"] as! [Int]
        let grid = try! GridScan(start: SIMD2(a[0], a[1]), end: SIMD2(e[0], e[1]), nx: g[0], ny: g[1])
        return try! STEM.scan(probe: probe, potentials: [potential], scan: grid,
                              detectors: [CenterOfMassDetector(maxAngle: .cutoff, units: .inverseAngstroms),
                                          CenterOfMassDetector(maxAngle: .cutoff, units: .milliradians)])
    }()
}

@Suite("Integrated centre of mass, integrate_gradient, integrate_disc, tile (abTEM 1.0.10)")
struct IcomPinsTests {
    @Test func integratedGradientMatchesAbTEM() throws {
        for case let c as [String: Any] in IcomPins.json["gradients"] as! [Any] {
            let label = c["label"] as! String, shape = c["shape"] as! [Int], s = IcomPins.vec(c["sampling"]!)
            let g = IcomPins.gradient("grad_\(label)_g", nx: shape[0], ny: shape[1], sx: s[0], sy: s[1])
            let out = try RealField.integratingGradient(x: g.x, y: g.y)
            let want = IcomPins.doubles("grad_\(label)_out")
            let excess = IcomPins.allcloseExcess(out.values, want)
            #expect(excess <= 0, "\(label): excess \(excess), relMax \(IcomPins.relMax(out.values, want))")
        }
    }

    @Test func plantedPotentialsAreRecoveredLikeAbTEMDoes() throws {
        for case let c as [String: Any] in IcomPins.json["gradients"] as! [Any] {
            let label = c["label"] as! String, shape = c["shape"] as! [Int], s = IcomPins.vec(c["sampling"]!)
            guard ["even", "odd", "mixed"].contains(label) else { continue }  // the 2 x 3 grid has only Nyquist frequencies
            // the potential less its minimum; relMax divides by its maximum, the potential's range
            let t = IcomPins.doubles("grad_\(label)_T")
            let low = t.min()!
            let want = t.map { $0 - low }
            let g = IcomPins.gradient("grad_\(label)_g", nx: shape[0], ny: shape[1], sx: s[0], sy: s[1])
            let kit = try RealField.integratingGradient(x: g.x, y: g.y).values
            #expect(IcomPins.relMax(IcomPins.doubles("grad_\(label)_out"), want) <= 1e-9, "abTEM \(label)")
            #expect(IcomPins.relMax(kit, want) <= 1e-9, "kit \(label)")
        }
    }

    @Test func integratedCenterOfMassOfAbTEMsCentreOfMassMatches() throws {
        let s = IcomPins.json["scan"] as! [String: Any]
        let sampling = IcomPins.vec(s["sampling"]!), n = s["scanGpts"] as! [Int]
        for (tag, index) in [("inv", 0), ("mrad", 1)] {
            // abTEM's own centre-of-mass image through the kit's integration: the strict bound
            let g = IcomPins.gradient("scan_com_\(tag)", nx: n[0], ny: n[1], sx: sampling[0], sy: sampling[1])
            let icom = try RealField.integratingGradient(x: g.x, y: g.y)
            let want = IcomPins.doubles("scan_icom_\(tag)")
            let excess = IcomPins.allcloseExcess(icom.values, want)
            #expect(excess <= 0, "\(tag): excess \(excess), relMax \(IcomPins.relMax(icom.values, want))")
            // the kit's centre of mass, registered at 1e-5 of abTEM's, carried through
            let kit = try IcomPins.scan[index].integratedCenterOfMass(samplingX: sampling[0], samplingY: sampling[1])
            #expect(kit.nx == n[0] && kit.ny == n[1])
            #expect(IcomPins.relMax(kit.values, want) <= 1e-5, "\(tag): kit COM \(IcomPins.relMax(kit.values, want))")
        }
        // abTEM's one-call chain `integrated_center_of_mass()` is the 1/Å one
        let chain = IcomPins.doubles("scan_icom_chain")
        #expect(IcomPins.allcloseExcess(chain, IcomPins.doubles("scan_icom_inv")) <= 0)
    }

    @Test func integrateDiscMatchesAbTEM() throws {
        let d = IcomPins.json["disc"] as! [String: Any]
        let shape = d["shape"] as! [Int], s = IcomPins.vec(d["sampling"]!)
        let image = RealField(nx: shape[0], ny: shape[1], values: IcomPins.doubles("disc_image"), samplingX: s[0], samplingY: s[1])
        for case let c as [String: Any] in d["cases"] as! [Any] {
            let p = IcomPins.vec(c["position"]!), r = (c["radius"] as! NSNumber).doubleValue
            let sum = try image.integrateDisc(x: p[0], y: p[1], radius: r)
            let mean = try image.integrateDisc(x: p[0], y: p[1], radius: r, mean: true)
            let wantSum = (c["sum"] as! NSNumber).doubleValue, wantMean = (c["mean"] as! NSNumber).doubleValue
            #expect(IcomPins.allcloseExcess([sum], [wantSum]) <= 0, "sum at \(p) r \(r): \(sum) vs \(wantSum)")
            #expect(IcomPins.allcloseExcess([mean], [wantMean]) <= 0, "mean at \(p) r \(r): \(mean) vs \(wantMean)")
        }
    }

    @Test func tileMatchesAbTEMExactly() throws {
        let d = IcomPins.json["disc"] as! [String: Any]
        let shape = d["shape"] as! [Int], s = IcomPins.vec(d["sampling"]!)
        let image = RealField(nx: shape[0], ny: shape[1], values: IcomPins.doubles("disc_image"), samplingX: s[0], samplingY: s[1])
        let t = IcomPins.json["tile"] as! [String: Any]
        let reps = t["repetitions"] as! [Int], tiledShape = t["shape"] as! [Int]
        let out = try image.tiled(repetitionsX: reps[0], repetitionsY: reps[1])
        #expect([out.nx, out.ny] == tiledShape)
        #expect(out.values == IcomPins.doubles("tile_2_3"))
    }
}

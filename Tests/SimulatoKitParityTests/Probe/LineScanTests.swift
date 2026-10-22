//
//  LineScanTests.swift
//  Role: LineScan (points, sampling, endpoint, add_margin, at_position, fractional) and STEM.scan over a LineScan with an
//        annular detector against abTEM 1.0.10 (tools/linescan/pins.py). The positions are NumPy's linspace of the same
//        float64 arithmetic; the scan values hold the other slices' bound (the kit's pi = abTEM's float32 pi: 1e-10).
//

import Foundation
import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoMultislice
import SimulatoProbe

@Suite("LineScan against abTEM 1.0.10")
struct LineScanTests {
    static func make(_ c: [String: Any]) throws -> LineScan {
        let a = c["args"] as! [String: Any]
        func opt(_ k: String) -> Double? { a[k].map { LinescanPins.num($0) } }
        func pt(_ k: String) -> SIMD2<Double>? { (a[k]).map { v in let p = LinescanPins.vec(v); return SIMD2(p[0], p[1]) } }
        let gpts = (a["gpts"] as? NSNumber)?.intValue
        let endpoint = (a["endpoint"] as? NSNumber)?.boolValue ?? true
        if a["atPosition"] != nil {
            return try LineScan.atPosition(center: pt("center")!, extent: opt("extent") ?? 1.0, angle: opt("angle") ?? 0.0,
                                           gpts: gpts, sampling: opt("sampling"), endpoint: endpoint)
        }
        var scan = try LineScan(start: pt("start") ?? .zero, end: pt("end")!, gpts: gpts, sampling: opt("sampling"),
                                endpoint: endpoint, fractional: (a["fractional"] as? NSNumber)?.boolValue ?? false,
                                extent: pt("extent2"))
        if let m = a["margin"] {
            let v = LinescanPins.vec(m)
            scan = try scan.addingMargin(before: v[0], after: v[1])
        }
        return scan
    }

    @Test func constructionMatchesAbTEM() throws {
        var worst = 0.0
        for c in LinescanPins.list("lineScans") {
            let label = c["label"] as! String
            let scan = try Self.make(c)
            let s = LinescanPins.vec(c["start"]!), e = LinescanPins.vec(c["end"]!)
            #expect(scan.gpts == (c["gpts"] as! Int), "\(label) gpts")
            #expect(scan.endpoint == (c["endpoint"] as! Bool), "\(label) endpoint")
            // exact but for the last bit: np.linalg.norm and the kit's sqrt(x² + y²) can differ by one ulp in the margin's direction
            #expect(allClose([scan.start.x, scan.start.y, scan.end.x, scan.end.y], s + e, rtol: 1e-14).passes, "\(label) end points")
            #expect(scan.sampling == LinescanPins.num(c["sampling"]), "\(label) sampling")
            #expect(scan.extent == LinescanPins.num(c["extent"]), "\(label) extent")
            if scan.extent > 0 { #expect(scan.angle == LinescanPins.num(c["angle"]), "\(label) angle") }
            let want = LinescanPins.doubles("line.\(label).positions")
            let got = scan.positions.flatMap { [$0.x, $0.y] }
            #expect(got.count == want.count, "\(label) count")
            let close = allClose(got, want)
            worst = max(worst, close.relativeToMax)
            #expect(close.passes, "\(label) positions")
            #expect(scan.shape == [scan.gpts])
        }
        print("LineScan positions: worst max|Δ|/max|p| = \(worst)")
    }

    @Test func refusals() throws {
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(1, 1)) }
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(1, 1), gpts: 0) }
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(1, 1), sampling: -1) }
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(.nan, 1), gpts: 3) }
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(1, 1), gpts: 3, fractional: true) }
        #expect(throws: SimulatoError.self) { try LineScan(start: SIMD2(1, 1), end: SIMD2(1, 1), sampling: 0.5) }
        let point = try LineScan(start: SIMD2(1, 1), end: SIMD2(1, 1), gpts: 2)
        #expect(throws: SimulatoError.self) { try point.addingMargin(0.1) }
        let line = try LineScan(end: SIMD2(1, 0), gpts: 5)
        #expect(throws: SimulatoError.self) { try line.addingMargin(-0.5) }  // both ends meet: no positions
    }

    /// The scan the pins ran: Al, 64², annular detector 50 to 150 mrad.
    static func run(_ scan: LineScan, pi: LobatoParameterization.PiConstant, backend: Backend = .cpu) throws -> ScanMeasurement {
        let n = LinescanPins.json["gpts"] as! Int
        let potential = try ProjectedPotential(structure: LinescanPins.structure, nx: n, ny: n,
                                               sliceThickness: LinescanPins.num(LinescanPins.json["sliceThickness"]),
                                               parameterization: LobatoParameterization(piConstant: pi))
        let beam = try ElectronBeam(energy: LinescanPins.num(LinescanPins.json["energy"]))
        let probe = try Probe(beam: beam, grid: potential.grid,
                              aperture: Aperture(semiangleCutoff: LinescanPins.num(LinescanPins.json["semiangle"])))
        return try STEM.scan(probe: probe, potentials: [potential], scan: scan, detectors: [AnnularDetector(inner: 50, outer: 150)],
                             backend: backend)[0]
    }

    static func scan(_ spec: [String: Any]) throws -> LineScan {
        let s = LinescanPins.vec(spec["start"]!), e = LinescanPins.vec(spec["end"]!)
        return try LineScan(start: SIMD2(s[0], s[1]), end: SIMD2(e[0], e[1]), gpts: spec["gpts"] as? Int,
                            endpoint: (spec["endpoint"] as? Bool) ?? true)
    }

    @Test func stemScanOverALineEqualsAbTEM() throws {
        for spec in LinescanPins.list("scans") {
            let label = spec["label"] as! String
            let scan = try Self.scan(spec)
            #expect(scan.sampling == LinescanPins.num(spec["sampling"]))
            let want = LinescanPins.doubles("scan.\(label).adf")
            for (pi, bound) in [(LobatoParameterization.PiConstant.abTEMFloat32, 1e-10), (.exact, 1e-6)] {
                let m = try Self.run(scan, pi: pi)
                #expect(m.scanShape == [scan.gpts] && m.valueShape == [])
                var d = 0.0
                for i in want.indices { d = max(d, abs(m.values[i] - want[i])) }
                let rel = d / (want.map(abs).max() ?? 1)
                print("STEM.scan LineScan \(label) pi \(pi): max|Δ|/max|ADF| = \(rel)")
                #expect(rel <= bound, "\(label) \(pi)")
            }
        }
    }

    @Test(.enabled(if: Multislice.metalIsAvailable, "no Metal device: the Metal line scan is not tested here"))
    func metalLineScanEqualsCPU() throws {
        let spec = LinescanPins.list("scans")[0]
        let scan = try Self.scan(spec)
        let cpu = try Self.run(scan, pi: .exact)
        let metal = try Self.run(scan, pi: .exact, backend: .metal)
        #expect(metal.scanShape == cpu.scanShape)
        var d = 0.0
        for i in cpu.values.indices { d = max(d, abs(metal.values[i] - cpu.values[i])) }
        let rel = d / (cpu.values.map(abs).max() ?? 1)
        print("STEM.scan LineScan metal against cpu: max|Δ|/max|ADF| = \(rel)")
        #expect(rel <= 2e-4, "the float32 bound of the other Metal scans")
    }
}

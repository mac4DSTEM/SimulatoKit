//
//  LineInterpolationTests.swift
//  Role: interpolate_line and interpolate_line_at_position on an image and on a centred diffraction pattern against
//        abTEM 1.0.10 (tools/pins/linescan/pins.py), at assert_allclose's float64 defaults (rtol 1e-7, atol 1e-12 of the maximum).
//        The truths that need no reference and the refusals are in SimulatoDetectorTests.
//

import Foundation
import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoProbe

@Suite("interpolate_line against abTEM 1.0.10")
struct LineInterpolationTests {
    static func note(_ what: String, _ c: Closeness) {
        print("\(what): max|Δ|/max|abTEM| = \(c.relativeToMax), assert_allclose ratio \(c.ratio)")
    }

    static var image: RealField {
        let s = LinescanPins.vec((LinescanPins.json["image"] as! [String: Any])["sampling"]!)
        let n = (LinescanPins.json["image"] as! [String: Any])["shape"] as! [Int]
        return RealField(nx: n[0], ny: n[1], values: LinescanPins.doubles("img.array"), samplingX: s[0], samplingY: s[1])
    }

    static var pattern: RealField {
        let s = LinescanPins.vec((LinescanPins.json["pattern"] as! [String: Any])["sampling"]!)
        let n = (LinescanPins.json["pattern"] as! [String: Any])["shape"] as! [Int]
        return RealField(nx: n[0], ny: n[1], values: LinescanPins.doubles("dp.array"), samplingX: s[0], samplingY: s[1])
    }

    @Test func interpolateLineEqualsAbTEM() throws {
        for spec in LinescanPins.list("interpLines") {
            let label = spec["label"] as! String
            let args = spec["args"] as! [String: Any]
            func pt(_ k: String) -> SIMD2<Double>? { args[k].map { let v = LinescanPins.vec($0); return SIMD2(v[0], v[1]) } }
            func opt(_ k: String) -> Double? { args[k].map { LinescanPins.num($0) } }
            let isImage = (spec["target"] as! String) == "image"
            let field = isImage ? Self.image : Self.pattern
            let gpts = args["gpts"] as? Int
            let order = (args["order"] as? Int) ?? 3
            let frac = (args["fractional"] as? Bool) ?? false
            let values: [Double], sampling: Double, start: SIMD2<Double>?, end: SIMD2<Double>?, width: Double?
            if (spec["kind"] as! String) == "line" {
                let endpoint = (args["endpoint"] as? Bool) ?? false
                if isImage {
                    let r = try field.interpolateLine(start: pt("start") ?? .zero, end: pt("end") ?? SIMD2(0, field.samplingX * Double(field.nx)),
                                                      sampling: opt("sampling"), gpts: gpts, width: opt("width") ?? 0,
                                                      margin: opt("margin") ?? 0, order: order, endpoint: endpoint, fractional: frac)
                    (values, sampling, start, end, width) = (r.values, r.sampling, r.start, r.end, r.lineWidth)
                } else {
                    let r = try field.interpolateReciprocalLine(start: pt("start") ?? .zero,
                                                                end: pt("end") ?? SIMD2(0, -Double(field.nx - 1) * field.samplingX),
                                                                sampling: opt("sampling"), gpts: gpts, width: opt("width") ?? 0,
                                                                margin: opt("margin") ?? 0, order: order, endpoint: endpoint,
                                                                fractional: frac)
                    (values, sampling, start, end, width) = (r.values, r.sampling, r.start, r.end, r.lineWidth)
                }
            } else {
                let endpoint = (args["endpoint"] as? Bool) ?? true
                if isImage {
                    let r = try field.interpolateLineAtPosition(center: pt("center")!, angle: opt("angle")!, extent: opt("extent")!,
                                                                gpts: gpts, sampling: opt("sampling"), width: opt("width") ?? 0,
                                                                order: order, endpoint: endpoint)
                    (values, sampling, start, end, width) = (r.values, r.sampling, r.start, r.end, r.lineWidth)
                } else {
                    let r = try field.interpolateReciprocalLineAtPosition(center: pt("center")!, angle: opt("angle")!,
                                                                          extent: opt("extent")!, gpts: gpts, sampling: opt("sampling"),
                                                                          width: opt("width") ?? 0, order: order, endpoint: endpoint)
                    (values, sampling, start, end, width) = (r.values, r.sampling, r.start, r.end, r.lineWidth)
                }
            }
            #expect(values.count == (spec["gpts"] as! Int), "\(label) points")
            #expect(sampling == LinescanPins.num(spec["sampling"]) || abs(sampling - LinescanPins.num(spec["sampling"])) <= 1e-14 * sampling,
                    "\(label) sampling \(sampling) against \(spec["sampling"]!)")
            let s = LinescanPins.vec(spec["start"]!), e = LinescanPins.vec(spec["end"]!)
            #expect(allClose([start!.x, start!.y, end!.x, end!.y], s + e, rtol: 1e-14).passes, "\(label) line")
            #expect(width == (spec["width"] as? NSNumber)?.doubleValue, "\(label) width metadata")
            let c = allClose(values, LinescanPins.doubles("\(label).values"))
            Self.note("interpolate_line \(label)", c)
            #expect(c.passes, "\(label)")
        }
    }
}

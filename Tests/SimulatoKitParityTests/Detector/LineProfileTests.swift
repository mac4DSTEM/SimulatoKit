//
//  LineProfileTests.swift
//  Role: RealSpaceLineProfiles / ReciprocalSpaceLineProfiles (width, interpolate) and interpolate_line /
//        interpolate_line_at_position against abTEM 1.0.10 (tools/pins/linescan/pins.py), at assert_allclose's float64 defaults
//        (rtol 1e-7, atol 1e-12 of the maximum). abTEM stores width and interpolate results as float32, so those
//        agree to float32 rounding (6e-8), inside the same bound.
//

import Foundation
import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoProbe

@Suite("Line profiles against abTEM 1.0.10")
struct LineProfileTests {

    static func note(_ what: String, _ c: Closeness) {
        print("\(what): max|Δ|/max|abTEM| = \(c.relativeToMax), assert_allclose ratio \(c.ratio)")
    }

    @Test func widthOfTheScanProfileAndOfPlantedProfilesEqualsAbTEM() throws {
        for spec in LinescanPins.list("scans") {
            let label = spec["label"] as! String
            let scan = try LineScanTests.scan(spec)
            let profile = RealSpaceLineProfiles(count: scan.gpts, values: LinescanPins.doubles("scan.\(label).adf"),
                                                sampling: scan.sampling, start: scan.start, end: scan.end)
            for h in [0.5, 0.2] {
                let w = try profile.width(height: h)
                let want = LinescanPins.num(spec["width\(h)"])
                let c = allClose(w, [want])
                Self.note("width \(label) \(h)", c)
                #expect(c.passes, "\(label) width at \(h): \(w) against \(want)")
            }
        }
        for spec in LinescanPins.list("planted") {
            let label = spec["label"] as! String
            let shape = (spec["shape"] as! [Int])
            let profile = RealSpaceLineProfiles(ensembleShape: Array(shape.dropLast()), count: shape.last!,
                                                values: LinescanPins.doubles("planted.\(label).array"),
                                                sampling: LinescanPins.num(spec["sampling"]))
            for h in [0.5, 0.25, 0.8] {
                let c = allClose(try profile.width(height: h), LinescanPins.doubles("planted.\(label).width\(h)"))
                Self.note("width \(label) \(h)", c)
                #expect(c.passes, "\(label) width at \(h)")
            }
        }
    }

    @Test func interpolateOfTheScanProfileEqualsAbTEM() throws {
        for spec in LinescanPins.list("scans") {
            let label = spec["label"] as! String
            let scan = try LineScanTests.scan(spec)
            let profile = RealSpaceLineProfiles(count: scan.gpts, values: LinescanPins.doubles("scan.\(label).adf"),
                                                sampling: scan.sampling, start: scan.start, end: scan.end)
            for (tag, entry) in spec["interp"] as! [String: [String: Any]] {
                let args = entry["args"] as! [String: Any]
                let r = try profile.interpolate(sampling: args["sampling"].map { LinescanPins.num($0) }, gpts: args["gpts"] as? Int,
                                                order: (args["order"] as? Int) ?? 3, endpoint: (args["endpoint"] as? Bool) ?? false)
                #expect(r.count == (entry["gpts"] as! Int), "\(label) \(tag) points")
                #expect(r.sampling == LinescanPins.num(entry["sampling"]), "\(label) \(tag) sampling")
                #expect(r.start == scan.start && r.end == scan.end, "the line travels with the profile")
                let c = allClose(r.values, LinescanPins.doubles("scan.\(label).interp.\(tag)"))
                Self.note("interpolate scan \(label) \(tag)", c)
                #expect(c.passes, "\(label) \(tag)")
            }
        }
    }

    @Test func interpolateOfPlantedProfilesEqualsAbTEM() throws {
        for spec in LinescanPins.list("planted") {
            let label = spec["label"] as! String
            let shape = (spec["shape"] as! [Int])
            let profile = RealSpaceLineProfiles(ensembleShape: Array(shape.dropLast()), count: shape.last!,
                                                values: LinescanPins.doubles("planted.\(label).array"),
                                                sampling: LinescanPins.num(spec["sampling"]))
            for (tag, entry) in spec["interp"] as! [String: [String: Any]] {
                let args = entry["args"] as! [String: Any]
                let r = try profile.interpolate(sampling: args["sampling"].map { LinescanPins.num($0) }, gpts: args["gpts"] as? Int,
                                                order: (args["order"] as? Int) ?? 3, endpoint: (args["endpoint"] as? Bool) ?? false)
                #expect(r.ensembleShape + [r.count] == entry["shape"] as! [Int], "\(label) \(tag) shape")
                #expect(r.sampling == LinescanPins.num(entry["sampling"]), "\(label) \(tag) sampling")
                let c = allClose(r.values, LinescanPins.doubles("planted.\(label).interp.\(tag)"))
                Self.note("interpolate \(label) \(tag)", c)
                #expect(c.passes, "\(label) \(tag)")
            }
        }
    }
}

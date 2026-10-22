//
//  LineProfileTruthTests.swift
//  Role: line profiles and line interpolation against truths that need no reference (the spline passes through its
//        knots, a plane is interpolated exactly, a triangle has a known width) and the refusals.
//

import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoProbe

@Suite("Line profiles: truths and refusals")
struct LineProfileTruthTests {
    @Test func theSplineReturnsTheDataAtTheKnots() throws {
        let data = (0..<9).map { Double(($0 * 7) % 5) - 1.5 }
        let profile = RealSpaceLineProfiles(count: 9, values: data, sampling: 0.5)
        for order in 0...5 {
            // the same sampling, any order: the profile back (the padded spline passes through every sample)
            let same = try profile.interpolate(order: order)
            #expect(same.count == 9 && same.sampling == 0.5)
            #expect(zip(same.values, data).allSatisfy { abs($0 - $1) < 1e-12 }, "order \(order)")
        }
    }

    @Test func theWidthOfATriangleIsKnown() throws {
        // [0 1 2 3 4 3 2 1 0] minus half the maximum changes sign between points 1/2 and 2/3, and 5/6 and 6/7... : the first
        // sign change is at 1, the last at 6 (the zero at points 2 and 6 counts through NumPy's sign), width (6 - 1) x 0.25
        let profile = RealSpaceLineProfiles(count: 9, values: [0, 1, 2, 3, 4, 3, 2, 1, 0], sampling: 0.25)
        #expect(try profile.width() == [1.25])
        // two profiles, a different maximum each
        let two = RealSpaceLineProfiles(ensembleShape: [2], count: 4, values: [0, 4, 4, 0, 0, 1, 1, 0], sampling: 2)
        #expect(try two.width() == [4, 4])
    }

    @Test func refusals() throws {
        let profile = RealSpaceLineProfiles(count: 4, values: [1, 2, 3, 4], sampling: 1)
        #expect(throws: SimulatoError.self) { try profile.interpolate(sampling: 0.5, gpts: 8) }
        #expect(throws: SimulatoError.self) { try profile.interpolate(gpts: 8, order: 6) }
        #expect(throws: SimulatoError.self) { try profile.interpolate(gpts: 0) }
        #expect(throws: SimulatoError.self) { try profile.interpolate(sampling: -1) }
        // abTEM raises an IndexError for a profile that never crosses the height (a flat one: the maximum itself is the threshold)
        #expect(throws: SimulatoError.self) { try RealSpaceLineProfiles(count: 3, values: [-1, -1, -1], sampling: 1).width() }
    }

    @Test func aLineScanMeasurementBecomesProfiles() throws {
        let scan = try LineScan(end: SIMD2(2, 0), gpts: 3)
        let one = ScanMeasurement(scanShape: [3], valueShape: [], values: [1, 2, 3])
        let p = try one.lineProfiles(for: scan)
        #expect(p.ensembleShape == [] && p.count == 3 && p.values == [1, 2, 3] && p.sampling == 1 && p.end == SIMD2(2, 0))
        #expect(one[1] == [2])
        // a vector detector (centre of mass): the scan axis moves last, as abTEM's moveaxis does
        let two = ScanMeasurement(scanShape: [3], valueShape: [2], values: [1, 10, 2, 20, 3, 30])
        let q = try two.lineProfiles(for: scan)
        #expect(q.ensembleShape == [2] && q.values == [1, 2, 3, 10, 20, 30])
        let grid = ScanMeasurement(scanShape: [2, 2], valueShape: [], values: [1, 2, 3, 4])
        #expect(throws: SimulatoError.self) { try grid.lineProfiles(for: scan) }
        let pattern = ScanMeasurement(scanShape: [3], valueShape: [2, 2], values: [Double](repeating: 0, count: 12))
        #expect(throws: SimulatoError.self) { try pattern.lineProfiles(for: scan) }
        #expect(throws: SimulatoError.self) { try one.gaussianSourceSize(sigma: 1, samplingX: 1, samplingY: 1) }
    }
}

//
//  PolarMeasurementTests.swift
//  Role: the operations on a SegmentedDetector scan (abTEM's PolarMeasurements.integrate, integrate_radial and
//        differentials, the DPC signal), against numbers abTEM 1.0.10 printed for a planted 2 x 3 scan of 2 radial x 4
//        azimuthal bins (values k^1.5, k the flat index; radial sampling 10 mrad from 5 mrad), and their refusals.
//

import Foundation
import Testing
import SimulatoCore
import SimulatoDetector

@Suite("Polar measurements (segmented detector)")
struct PolarMeasurementTests {
    /// The planted scan and the detector whose bins it fills: inner 5, outer 25 mrad, 2 x 4 bins.
    static func planted() throws -> (ScanMeasurement, SegmentedDetector) {
        let values = (0..<(2 * 3 * 2 * 4)).map { pow(Double($0), 1.5) }
        let m = ScanMeasurement(scanShape: [2, 3], valueShape: [2, 4], values: values)
        return (m, try SegmentedDetector(radialBins: 2, azimuthalBins: 4, inner: 5, outer: 25))
    }

    /// assert_allclose: abTEM 1.0.10 returns float64 sums here; differentials are complex64 (rtol 1e-6).
    static func close(_ a: [Double], _ b: [Double], rtol: Double) -> Bool {
        a.count == b.count && zip(a, b).allSatisfy { abs($0 - $1) <= rtol * abs($1) }
    }

    @Test func detectorRegionsAreSummedAsAbTEM() throws {
        let (m, _) = try Self.planted()
        let r = try m.integrate(detectorRegions: [1, 6])
        #expect(r.scanShape == [2, 3] && r.valueShape == [])
        #expect(Self.close(r.values, [15.696938456699069, 79.38320341483518, 173.28194235161567, 289.31676725154983,
                                      423.81829944857606, 574.5152729585091], rtol: 1e-7))
        #expect(Self.close(try m.integrate(detectorRegions: [3]).values,
                           [5.196152422706632, 36.4828726939094, 82.8190799272728, 140.29611541307906,
                            207.06279240848656, 281.96985654498604], rtol: 1e-7))
        // A region given twice counts twice.
        #expect(Self.close(try m.integrate(detectorRegions: [3, 3]).values,
                           [10.392304845413264, 72.9657453878188, 165.6381598545456, 280.5922308261581,
                            414.1255848169731, 563.9397130899721], rtol: 1e-7))
        // A line scan keeps its one scan axis.
        let line = ScanMeasurement(scanShape: [6], valueShape: [2, 4], values: m.values)
        #expect(try line.integrate(detectorRegions: [1, 6]).scanShape == [6])
        #expect(try line.integrate(detectorRegions: [1, 6]).values == (try m.integrate(detectorRegions: [1, 6])).values)
    }

    @Test func radialLimitsTruncateToBinsAsAbTEM() throws {
        let (m, d) = try Self.planted()
        // (15, 25) is bin 1; (5, 25) both bins (the whole detector); (9, 20) truncates to bin 0.
        #expect(Self.close(try m.integrateRadial(inner: 15, outer: 25, detector: d).values,
                           [52.39753752165015, 198.91933957063137, 399.1700804463722, 641.2493153257982,
                            918.8648676713927, 1227.935103956358], rtol: 1e-7))
        let whole: [Double] = [61.422117069102974, 316.65240586419407, 692.4494883772924, 1156.6954457458824,
                               1694.7699278241307, 2297.6063763826915]
        #expect(Self.close(try m.integrateRadial(inner: 5, outer: 25, detector: d).values, whole, rtol: 1e-7))
        #expect(Self.close(try m.integrate(detector: d).values, whole, rtol: 1e-7))
        // int(−0.5) is 0 in abTEM: an inner angle half a bin below the detector's counts from bin 0.
        #expect(Self.close(try m.integrateRadial(inner: 0, outer: 25, detector: d).values, whole, rtol: 1e-7))
        // Inner above outer, or an inner angle far past the bins: an empty slice, zeros.
        #expect(try m.integrateRadial(inner: 25, outer: 15, detector: d).values == [0, 0, 0, 0, 0, 0])
        #expect(try m.integrateRadial(inner: 1e30, outer: 25, detector: d).values == [0, 0, 0, 0, 0, 0])
        #expect(Self.close(try m.integrateRadial(inner: 9, outer: 20, detector: d).values,
                           [9.024579547452822, 117.73306629356271, 293.2794079309202, 515.446130420084,
                            775.9050601527379, 1069.6712724263334], rtol: 1e-7))
    }

    @Test func differentialsSubtractTheFirstRegionsFromTheSecond() throws {
        let (m, _) = try Self.planted()
        // abTEM differentials(((0, 4), (2, 6)), ((1, 5), (3, 7))): integrate(item 1) − integrate(item 0), complex64.
        let d = try m.differentials(direction1: ([0, 4], [2, 6]), direction2: ([1, 5], [3, 7]))
        #expect(d.direction1.scanShape == [2, 3] && d.direction2.valueShape == [])
        #expect(Self.close(d.direction1.values, [9.525365829467773, 19.809343338012695, 26.11396026611328,
                                                 31.15369415283203, 35.48075866699219, 39.33309555053711], rtol: 1e-6))
        #expect(Self.close(d.direction2.values, [11.53607177734375, 20.705455780029297, 26.79631996154785,
                                                 31.727031707763672, 35.98493194580078, 39.788352966308594], rtol: 1e-6))
    }

    @Test func aSymmetricPatternHasNoDifferential() throws {
        // Four quadrants of equal intensity at every position: both directions are 0.
        let m = ScanMeasurement(scanShape: [3], valueShape: [1, 4], values: [Double](repeating: 2.5, count: 12))
        let d = try m.differentials(direction1: ([0], [2]), direction2: ([1], [3]))
        #expect(d.direction1.values == [0, 0, 0] && d.direction2.values == [0, 0, 0])
        #expect(d.direction1.scanShape == [3])
    }

    @Test func refusals() throws {
        let (m, d) = try Self.planted()
        // Regions outside the bins, or none.
        for bad in [[8], [-1], [Int]()] { #expect(throws: SimulatoError.self) { try m.integrate(detectorRegions: bad) } }
        // Beyond the outer bin (abTEM: "Integration limit exceeded"), below the inner angle, not finite.
        #expect(throws: SimulatoError.self) { try m.integrateRadial(inner: 5, outer: 35, detector: d) }
        #expect(throws: SimulatoError.self) { try m.integrateRadial(inner: -20, outer: 25, detector: d) }
        #expect(throws: SimulatoError.self) { try m.integrateRadial(inner: .nan, outer: 25, detector: d) }
        // A detector whose bins are not the measurement's.
        let other = try SegmentedDetector(radialBins: 1, azimuthalBins: 4, inner: 5, outer: 25)
        #expect(throws: SimulatoError.self) { try m.integrate(detector: other) }
        // A one-value measurement has no bins.
        let annular = ScanMeasurement(scanShape: [2], valueShape: [], values: [1, 2])
        #expect(throws: SimulatoError.self) { try annular.integrate(detectorRegions: [0]) }
        #expect(throws: SimulatoError.self) { try annular.differentials(direction1: ([0], [0]), direction2: ([0], [0])) }
    }
}

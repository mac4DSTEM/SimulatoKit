//
//  PolarMeasurements.swift
//  Role: abTEM's operations on a polar measurement (`PolarMeasurements`, the scan of a `SegmentedDetector`): `integrate`
//        over detector regions, `integrate_radial` between two angles, and `differentials`, the differential phase
//        contrast (DPC) signal of two pairs of regions.
//

import Foundation
import SimulatoCore

extension ScanMeasurement {
    /// The sum of the given detector regions at every scan position (abTEM `PolarMeasurements.integrate(detector_regions=)`).
    /// A region is a flat bin index, radial-major: `azimuthal + radial × azimuthalBins`, as `SegmentedDetector` fills
    /// them. A region given twice counts twice, as in abTEM. Any scan with two value axes is read as radial × azimuthal;
    /// the bins are not checked against a detector (use `integrate(detector:)` or `integrateRadial` for that).
    /// - Returns: one value per scan position (`valueShape` `[]`).
    /// - Throws: `invalidScan` when the measurement does not have two value axes, for no region (abTEM returns zeros) or
    ///   one outside the bins (abTEM would read a negative one from the end).
    public func integrate(detectorRegions: [Int]) throws -> ScanMeasurement {
        let n = try polarBinCount()
        guard !detectorRegions.isEmpty else { throw SimulatoError.invalidScan("no detector region given") }
        for r in detectorRegions where !(0..<n).contains(r) {
            throw SimulatoError.invalidScan("detector region \(r) is outside the \(n) bins")
        }
        return reduced { bins in detectorRegions.reduce(0.0) { $0 + bins[$1] } }
    }

    /// The sum of every bin of `detector` at every scan position (abTEM `PolarMeasurements.integrate()`).
    /// - Throws: `invalidScan` when the measurement's bins are not the detector's.
    public func integrate(detector: SegmentedDetector) throws -> ScanMeasurement {
        try checkBins(detector)
        return reduced { bins in bins.reduce(0, +) }
    }

    /// The sum of the radial bins of `detector` between `inner` and `outer` [mrad], all azimuths (abTEM
    /// `PolarMeasurements.integrate_radial(inner, outer)`): the bins from trunc((inner − detector.inner) / Δ) to
    /// trunc((outer − detector.inner) / Δ), the last excluded, with Δ the radial bin width (truncated toward zero, as
    /// abTEM's `int`, so an inner angle up to one bin below the detector's counts from bin 0). An angle that does not fall
    /// on a bin edge selects the bin it truncates to; inner above outer gives zeros, as abTEM's empty slice.
    /// - Throws: `invalidScan` when the bins are not the detector's, for a non-finite angle, an outer index beyond the last
    ///   bin (abTEM's "Integration limit exceeded") or a negative index (which abTEM would read from the end).
    public func integrateRadial(inner: Double, outer: Double, detector: SegmentedDetector) throws -> ScanMeasurement {
        try checkBins(detector)
        guard inner.isFinite, outer.isFinite else {
            throw SimulatoError.invalidScan("the integration limits must be finite, got (\(inner), \(outer)) mrad")
        }
        let width = (detector.outer - detector.inner) / Double(detector.radialBins)
        let a = ((inner - detector.inner) / width).rounded(.towardZero)
        let b = ((outer - detector.inner) / width).rounded(.towardZero)
        guard a >= 0, b >= 0 else {
            throw SimulatoError.invalidScan("the integration limits (\(inner), \(outer)) mrad fall below the detector's inner angle")
        }
        guard b <= Double(detector.radialBins) else {
            throw SimulatoError.invalidScan("the outer limit \(outer) mrad exceeds the detector's \(detector.outer) mrad")
        }
        // An inner index past the last bin is an empty slice in abTEM (zeros); clamped so it cannot overflow Int.
        let first = Int(min(a, Double(detector.radialBins))), last = Int(b), nAz = detector.azimuthalBins
        return reduced { bins in
            var sum = 0.0
            if first < last { for i in (first * nAz)..<(last * nAz) { sum += bins[i] } }
            return sum
        }
    }

    /// The differential (DPC) signal of two directions (abTEM `PolarMeasurements.differentials(direction_1,
    /// direction_2)`): for each direction the integral of its second regions minus that of its first, which is abTEM's
    /// order (its code subtracts item 0 from item 1). abTEM's complex image is `direction1 + i direction2`; summed here in
    /// float64 where abTEM stores complex64.
    /// - Throws: as `integrate(detectorRegions:)`.
    public func differentials(direction1: ([Int], [Int]), direction2: ([Int], [Int])) throws
        -> (direction1: ScanMeasurement, direction2: ScanMeasurement) {
        func difference(_ d: ([Int], [Int])) throws -> ScanMeasurement {
            let plus = try integrate(detectorRegions: d.1), minus = try integrate(detectorRegions: d.0)
            return ScanMeasurement(scanShape: scanShape, valueShape: [], values: zip(plus.values, minus.values).map { $0 - $1 })
        }
        return (try difference(direction1), try difference(direction2))
    }

    private func polarBinCount() throws -> Int {
        guard valueShape.count == 2 else {
            throw SimulatoError.invalidScan("a polar measurement has radial × azimuthal bins; this one has value shape \(valueShape)")
        }
        return valueCount
    }

    private func checkBins(_ detector: SegmentedDetector) throws {
        _ = try polarBinCount()
        guard valueShape == [detector.radialBins, detector.azimuthalBins] else {
            throw SimulatoError.invalidScan("the measurement has \(valueShape) bins, the detector \([detector.radialBins, detector.azimuthalBins])")
        }
    }

    /// One number per scan position from that position's bins (indexed from 0).
    private func reduced(_ f: (UnsafeBufferPointer<Double>) -> Double) -> ScanMeasurement {
        let v = valueCount, positions = values.count / max(v, 1)
        let out = values.withUnsafeBufferPointer { all in
            (0..<positions).map { p in f(UnsafeBufferPointer(rebasing: all[(p * v)..<((p + 1) * v)])) }
        }
        return ScanMeasurement(scanShape: scanShape, valueShape: [], values: out)
    }
}

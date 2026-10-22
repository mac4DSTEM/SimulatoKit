//
//  LineProfiles.swift
//  Role: abTEM 1.0.10's `RealSpaceLineProfiles` and `ReciprocalSpaceLineProfiles` (measurements.py, `_BaseMeasurement1D`):
//        one or more profiles over a last axis with a sampling, their `width` (the distance between the first and last
//        crossing of a fraction of the maximum) and `interpolate` (a spline through the profile padded periodically by 3
//        points, `map_coordinates`), and what `STEM.scan` returns for a `LineScan`: `ScanMeasurement.lineProfiles(for:)`.
//

import Foundation
import SimulatoCore
import SimulatoProbe

/// One or more profiles along a line: `values` is row-major over `ensembleShape` and then the points of the line (abTEM's
/// last axis).
public protocol LineProfiles: Sendable, Hashable {
    /// The shape of the leading axes (empty for a single profile).
    var ensembleShape: [Int] { get }
    /// The points of each profile.
    var count: Int { get }
    /// `ensembleShape` × `count` values, the points fastest.
    var values: [Double] { get }
    /// The distance between points: Å for real-space profiles, 1/Å for reciprocal-space profiles.
    var sampling: Double { get }
    /// The line the profiles were taken along [Å or 1/Å in the measurement's frame], when known.
    var start: SIMD2<Double>? { get }
    var end: SIMD2<Double>? { get }
    /// The width the profiles were averaged over across the line, when they were.
    var lineWidth: Double? { get }

    /// - Precondition: `values` holds the product of `ensembleShape` × `count` numbers.
    init(ensembleShape: [Int], count: Int, values: [Double], sampling: Double, start: SIMD2<Double>?, end: SIMD2<Double>?,
         lineWidth: Double?)
}

/// Real-space line profiles (abTEM `RealSpaceLineProfiles`): sampling in Å.
public struct RealSpaceLineProfiles: LineProfiles {
    public let ensembleShape: [Int]
    public let count: Int
    public let values: [Double]
    public let sampling: Double
    public let start: SIMD2<Double>?
    public let end: SIMD2<Double>?
    public let lineWidth: Double?

    public init(ensembleShape: [Int] = [], count: Int, values: [Double], sampling: Double, start: SIMD2<Double>? = nil,
                end: SIMD2<Double>? = nil, lineWidth: Double? = nil) {
        precondition(values.count == ensembleShape.reduce(1, *) * count, "RealSpaceLineProfiles: the values do not fill the shape")
        self.ensembleShape = ensembleShape
        self.count = count
        self.values = values
        self.sampling = sampling
        self.start = start
        self.end = end
        self.lineWidth = lineWidth
    }
}

/// Reciprocal-space line profiles (abTEM `ReciprocalSpaceLineProfiles`): sampling in 1/Å.
public struct ReciprocalSpaceLineProfiles: LineProfiles {
    public let ensembleShape: [Int]
    public let count: Int
    public let values: [Double]
    public let sampling: Double
    public let start: SIMD2<Double>?
    public let end: SIMD2<Double>?
    public let lineWidth: Double?

    public init(ensembleShape: [Int] = [], count: Int, values: [Double], sampling: Double, start: SIMD2<Double>? = nil,
                end: SIMD2<Double>? = nil, lineWidth: Double? = nil) {
        precondition(values.count == ensembleShape.reduce(1, *) * count, "ReciprocalSpaceLineProfiles: the values do not fill the shape")
        self.ensembleShape = ensembleShape
        self.count = count
        self.values = values
        self.sampling = sampling
        self.start = start
        self.end = end
        self.lineWidth = lineWidth
    }
}

extension LineProfiles {
    /// The extent of a profile: `sampling` × `count` (abTEM `extent`).
    public var extent: Double { sampling * Double(count) }

    private var profileCount: Int { ensembleShape.reduce(1, *) }

    /// The width of each profile at a fraction `height` of its maximum (abTEM `width`): the profile minus `height` × its
    /// maximum, then the distance between the first and the last sign change, times the sampling. The widths are in
    /// `ensembleShape` order (one number for a single profile). abTEM stores them as float32; these are float64.
    /// - Throws: `invalidScan` for a profile that never crosses the height (abTEM raises an IndexError).
    public func width(height: Double = 0.5) throws -> [Double] {
        var out: [Double] = []
        out.reserveCapacity(profileCount)
        for p in 0..<profileCount {
            let row = values[(p * count)..<((p + 1) * count)]
            let threshold = (row.max() ?? 0) * height
            func sign(_ v: Double) -> Double { v > 0 ? 1 : (v < 0 ? -1 : (v == 0 ? 0 : .nan)) }
            var previous = sign(row[row.startIndex] - threshold)
            var left = -1, right = -1
            for k in 1..<max(count, 1) {
                let current = sign(row[row.startIndex + k] - threshold)
                if current - previous != 0 {  // NumPy: diff != 0, true also for a NaN
                    if left < 0 { left = k - 1 }
                    right = k - 1
                }
                previous = current
            }
            guard left >= 0 else {
                throw SimulatoError.invalidScan("profile \(p) never crosses \(height) of its maximum, so it has no width")
            }
            out.append(Double(right - left) * sampling)
        }
        return out
    }

    /// The profiles on another sampling (abTEM `interpolate`): a spline of `order` through each profile, which is padded by
    /// 3 points of its own periodic continuation, at `gpts` points from the first point (`endpoint` false) or ending on
    /// the last (`endpoint` true) of the original. Give `sampling` or `gpts`, not both; neither keeps the sampling.
    /// - Parameters:
    ///   - sampling: the new sampling; the points are ⌈extent / sampling⌉.
    ///   - gpts: the new number of points; the sampling is extent / gpts.
    ///   - order: the spline order, 0 to 5 (default 3).
    ///   - endpoint: false (default) as abTEM.
    /// - Throws: `invalidScan` for both `sampling` and `gpts`, a non-positive value, or an order outside 0 to 5.
    public func interpolate(sampling newSampling: Double? = nil, gpts: Int? = nil, order: Int = 3,
                            endpoint: Bool = false) throws -> Self {
        if gpts != nil && newSampling != nil { throw SimulatoError.invalidScan("give sampling or gpts to interpolate, not both") }
        guard SplineInterpolation.orders.contains(order) else {
            throw SimulatoError.invalidScan("the spline order must be 0 to 5, got \(order)")
        }
        var s = newSampling, n = gpts
        if s == nil && n == nil { s = sampling }
        if n == nil {
            guard let s, s.isFinite, s > 0 else { throw SimulatoError.invalidScan("sampling must be finite and positive") }
            let c = (extent / s).rounded(.up)
            guard c >= 1, c < 1e9 else { throw SimulatoError.invalidScan("sampling \(s) gives \(c) points") }
            n = Int(c)
        }
        guard let points = n, points >= 1 else { throw SimulatoError.invalidScan("gpts must be at least 1") }
        if s == nil { s = extent / Double(points) }
        // abTEM: wrap-pad by 3, then linspace(3, N + 3, gpts, endpoint) in the padded index space.
        let positions = GridScan.linspace(3.0, Double(count + 6) - 3.0, points, endpoint: endpoint)
        var out = [Double](repeating: 0, count: profileCount * points)
        var padded = [Double](repeating: 0, count: count + 6)
        for p in 0..<profileCount {
            for j in 0..<(count + 6) {
                var m = (j - 3) % count
                if m < 0 { m += count }
                padded[j] = values[p * count + m]
            }
            let c = SplineInterpolation.coefficients(padded, order: order)
            for (k, x) in positions.enumerated() { out[p * points + k] = SplineInterpolation.evaluate(c, at: x, order: order) }
        }
        return Self(ensembleShape: ensembleShape, count: points, values: out, sampling: s!, start: start, end: end,
                    lineWidth: lineWidth)
    }
}

extension ScanMeasurement {
    /// The value at scan position `i` of a one-axis (line) scan.
    public subscript(i: Int) -> [Double] {
        precondition(scanShape.count == 1, "ScanMeasurement: one index needs a line scan; use (ix, iy) for a grid scan")
        return Array(values[(i * valueCount)..<((i + 1) * valueCount)])
    }

    /// The line profiles of a scan over `scan` (what abTEM returns for a `LineScan` with a detector of one number, or of
    /// a vector of numbers whose axis then comes first): the scan axis moves to the last axis, the sampling is the scan's.
    /// - Throws: `invalidScan` when this is not a one-axis scan of `scan`'s number of positions, or a detector value has
    ///   more than one axis (abTEM returns other measurement types for those).
    public func lineProfiles(for scan: LineScan) throws -> RealSpaceLineProfiles {
        guard scanShape == [scan.gpts] else {
            throw SimulatoError.invalidScan("the measurement has scan shape \(scanShape), the line scan \(scan.gpts) positions")
        }
        guard valueShape.count <= 1 else {
            throw SimulatoError.invalidScan("line profiles need detectors of one number or a vector, got value shape \(valueShape)")
        }
        let n = scan.gpts, v = valueCount
        var moved = [Double](repeating: 0, count: values.count)
        for i in 0..<n { for k in 0..<v { moved[k * n + i] = values[i * v + k] } }
        return RealSpaceLineProfiles(ensembleShape: valueShape, count: n, values: moved, sampling: scan.sampling,
                                     start: scan.start, end: scan.end)
    }
}

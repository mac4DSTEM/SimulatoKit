//
//  LineScan.swift
//  Role: STEM scan positions along a straight line (abTEM 1.0.10 `LineScan`, scan.py): NumPy's
//        linspace(start, end, gpts, endpoint) along each axis, the number of points and the sampling kept consistent
//        the way abTEM's setters keep them (`_adjust_gpts`, `_adjust_sampling`), with `at_position` and `add_margin`.
//        `ScanPattern` is what `STEM.scan` takes: the positions and the shape of the scan axes.
//

import Foundation
import SimulatoCore

/// A set of probe positions with the shape of the scan axes in the result: `[nx, ny]` for a `GridScan`, `[n]` for a
/// `LineScan`.
public protocol ScanPattern: Sendable {
    /// The positions [Å], first scan axis slowest.
    var positions: [SIMD2<Double>] { get }
    /// The number of positions along each scan axis; the product is `positions.count`.
    var shape: [Int] { get }
}

extension GridScan: ScanPattern {
    public var shape: [Int] { [nx, ny] }
}

/// Probe positions along a straight line (abTEM `LineScan`).
///
/// The line runs from `start` to `end` [Å]. `gpts` points are placed along it, the last at `end` when `endpoint` is true.
/// As in abTEM, a given `sampling` wins over a given `gpts`: the number of points is ⌈length / sampling⌉ and the realised
/// sampling is the length over that number (over one less with `endpoint`).
public struct LineScan: Sendable, Hashable, ScanPattern {
    public private(set) var start: SIMD2<Double>
    public private(set) var end: SIMD2<Double>
    /// The number of positions.
    public private(set) var gpts: Int
    /// The distance between positions [Å]: length / (gpts − 1) with `endpoint` and more than one point, else length / gpts.
    public private(set) var sampling: Double
    public let endpoint: Bool

    /// - Parameters:
    ///   - start: the first point [Å], or, with `fractional`, as fractions of `extent`.
    ///   - end: the last point [Å] (a fraction of `extent` with `fractional`).
    ///   - gpts: the number of positions; give this or `sampling`.
    ///   - sampling: the largest distance between positions [Å]; give this or `gpts`.
    ///   - endpoint: true (abTEM's default) puts the last position at `end`.
    ///   - fractional: `start` and `end` are fractions of `extent` along each axis (abTEM `fractional=True` with a potential).
    ///   - extent: the extent [Å] of the potential the fractions refer to (its grid's `extentX`, `extentY`).
    /// - Throws: `invalidScan` for non-finite points, neither `gpts` nor `sampling`, a non-positive `gpts` or `sampling`, a
    ///   scan of no positions, or `fractional` without `extent`.
    public init(start: SIMD2<Double> = .zero, end: SIMD2<Double>, gpts: Int? = nil, sampling: Double? = nil,
                endpoint: Bool = true, fractional: Bool = false, extent: SIMD2<Double>? = nil) throws {
        var s = start, e = end
        if fractional {
            guard let extent else { throw SimulatoError.invalidScan("provide the potential's extent for fractional coordinates") }
            s = SIMD2(extent.x * s.x, extent.y * s.y)
            e = SIMD2(extent.x * e.x, extent.y * e.y)
        }
        guard s.x.isFinite, s.y.isFinite, e.x.isFinite, e.y.isFinite else {
            throw SimulatoError.invalidScan("start and end must be finite")
        }
        if let gpts { guard gpts >= 1 else { throw SimulatoError.invalidScan("gpts must be at least 1, got \(gpts)") } }
        if let sampling {
            guard sampling.isFinite, sampling > 0 else {
                throw SimulatoError.invalidScan("sampling must be finite and positive, got \(sampling) Å")
            }
        }
        guard gpts != nil || sampling != nil else {
            throw SimulatoError.invalidScan("provide gpts or sampling for a line scan")
        }
        self.start = s
        self.end = e
        self.endpoint = endpoint
        self.gpts = gpts ?? 0
        self.sampling = sampling ?? 0
        if sampling != nil { try adjustGpts() }
        adjustSampling()
    }

    /// The length of the line [Å] (abTEM `extent`).
    public var extent: Double {
        let d = end - start
        return (d.x * d.x + d.y * d.y).squareRoot()
    }

    /// The unit vector from `start` to `end`.
    public var direction: SIMD2<Double> {
        let d = end - start
        return d / extent
    }

    /// The angle of the line to the first axis [rad] (abTEM's `angle`, which documents degrees and returns `arctan2`).
    public var angle: Double {
        let d = direction
        return atan2(d.y, d.x)
    }

    public var shape: [Int] { [gpts] }

    /// abTEM's `_adjust_gpts` (a sampling is always set when this runs): ⌈length / sampling⌉ points.
    private mutating func adjustGpts() throws {
        let n = (extent / sampling).rounded(.up)
        guard n >= 1, n < 1e9 else {
            throw SimulatoError.invalidScan("the line of \(extent) Å at sampling \(sampling) Å has \(n) positions")
        }
        gpts = Int(n)
        adjustSampling()
    }

    /// abTEM's `_adjust_sampling`.
    private mutating func adjustSampling() {
        sampling = endpoint && gpts > 1 ? extent / Double(gpts - 1) : extent / Double(gpts)
    }

    /// The line with a margin added at both ends (abTEM `add_margin`): `start` moves back along the line by `before` and
    /// `end` forward by `after` [Å]. As abTEM's setters do, the number of points is recomputed from the sampling the line
    /// had before each end moved, so it can change.
    /// - Throws: `invalidScan` for a line of zero length, a non-finite margin or a line left with no positions.
    public func addingMargin(before: Double, after: Double) throws -> LineScan {
        guard before.isFinite, after.isFinite else { throw SimulatoError.invalidScan("the margin must be finite") }
        guard extent > 0 else { throw SimulatoError.invalidScan("a line of zero length has no direction for a margin") }
        var line = self
        let d = direction
        line.start = start - d * before
        try line.adjustGpts()
        line.end = end + d * after
        try line.adjustGpts()
        return line
    }

    /// The same margin at both ends.
    public func addingMargin(_ margin: Double) throws -> LineScan { try addingMargin(before: margin, after: margin) }

    /// The end points of a line of `extent` [Å] centred on `center` at `angle` [deg] to the first axis (abTEM's arithmetic).
    static func endpoints(center: SIMD2<Double>, extent: Double, angle: Double) -> (start: SIMD2<Double>, end: SIMD2<Double>) {
        let a = angle * (Double.pi / 180)
        let d = SIMD2(cos(a), sin(a))
        return (SIMD2(center.x - extent / 2 * d.x, center.y - extent / 2 * d.y),
                SIMD2(center.x + extent / 2 * d.x, center.y + extent / 2 * d.y))
    }

    /// A line centred on `center` [Å] (abTEM `LineScan.at_position`).
    /// - Parameters:
    ///   - extent: the length of the line [Å] (abTEM's default is 1).
    ///   - angle: the angle to the first axis [deg].
    ///   - endpoint: true (default) puts the last position at the end of the line.
    public static func atPosition(center: SIMD2<Double>, extent: Double = 1.0, angle: Double = 0.0, gpts: Int? = nil,
                                  sampling: Double? = nil, endpoint: Bool = true) throws -> LineScan {
        let p = endpoints(center: center, extent: extent, angle: angle)
        return try LineScan(start: p.start, end: p.end, gpts: gpts, sampling: sampling, endpoint: endpoint)
    }

    /// The positions [Å]: linspace(start, end, gpts, endpoint) along each axis.
    public var positions: [SIMD2<Double>] {
        let xs = GridScan.linspace(start.x, end.x, gpts, endpoint: endpoint)
        let ys = GridScan.linspace(start.y, end.y, gpts, endpoint: endpoint)
        return zip(xs, ys).map { SIMD2($0, $1) }
    }
}

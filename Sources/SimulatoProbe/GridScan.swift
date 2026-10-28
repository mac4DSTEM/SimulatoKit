//
//  GridScan.swift
//  Role: STEM scan positions (abTEM `GridScan.get_positions`): along each axis
//        NumPy's linspace(start, end, gpts, endpoint), x the first index.
//

import SimulatoCore

/// A rectangular grid of probe positions.
public struct GridScan: Sendable, Hashable {
    public let start: SIMD2<Double>
    public let end: SIMD2<Double>
    public let nx: Int
    public let ny: Int
    public let endpoint: Bool

    public init(start: SIMD2<Double>, end: SIMD2<Double>, nx: Int, ny: Int, endpoint: Bool = false) throws {
        guard nx >= 1, ny >= 1 else { throw SimulatoError.invalidScan("gpts must be at least 1, got \(nx) x \(ny)") }
        guard start.x.isFinite, start.y.isFinite, end.x.isFinite, end.y.isFinite else {
            throw SimulatoError.invalidScan("start and end must be finite")
        }
        guard end.x - start.x > 0 || end.y - start.y > 0 else {
            throw SimulatoError.invalidScan("scan extent must be positive on at least one axis (abTEM's rule)")
        }
        self.start = start
        self.end = end
        self.nx = nx
        self.ny = ny
        self.endpoint = endpoint
    }

    /// A scan whose pixel size is at most `samplingX` × `samplingY` [Å] (abTEM `GridScan(sampling=)`): ⌈extent / sampling⌉
    /// points per axis, so the realised sampling is the extent over the points (⌈extent / sampling⌉ + 1 with `endpoint`).
    public init(start: SIMD2<Double>, end: SIMD2<Double>, samplingX: Double, samplingY: Double,
                endpoint: Bool = false) throws {
        guard samplingX.isFinite, samplingY.isFinite, samplingX > 0, samplingY > 0 else {
            throw SimulatoError.invalidScan("sampling must be finite and positive, got (\(samplingX), \(samplingY)) Å")
        }
        let ex = end.x - start.x, ey = end.y - start.y
        guard ex.isFinite, ey.isFinite else { throw SimulatoError.invalidScan("start and end must be finite") }
        guard ex >= 0, ey >= 0 else { throw SimulatoError.invalidScan("end must not lie before start") }
        // A zero-length axis has a single point; abTEM adds one point per axis with `endpoint`.
        let px = ex > 0 ? (ex / samplingX).rounded(.up) : 0, py = ey > 0 ? (ey / samplingY).rounded(.up) : 0
        guard px < 1e9, py < 1e9 else { throw SimulatoError.invalidScan("sampling gives too many points") }
        let extra = endpoint ? 1 : 0
        try self.init(start: start, end: end, nx: max(Int(px) + (ex > 0 ? extra : 1), 1),
                      ny: max(Int(py) + (ey > 0 ? extra : 1), 1), endpoint: endpoint)
    }

    /// NumPy's `linspace(start, stop, num, endpoint)`: i · step + start, step = (stop − start) / (num or num − 1); the
    /// last value is stop exactly when endpoint is true.
    public static func linspace(_ a: Double, _ b: Double, _ n: Int, endpoint: Bool) -> [Double] {
        let div = endpoint ? n - 1 : n
        guard div > 0 else { return [a] }
        let step = (b - a) / Double(div)
        var out = (0..<n).map { Double($0) * step + a }
        if endpoint && n > 1 { out[n - 1] = b }
        return out
    }

    /// The positions [Å], index ix * ny + iy.
    public var positions: [SIMD2<Double>] {
        let xs = Self.linspace(start.x, end.x, nx, endpoint: endpoint)
        let ys = Self.linspace(start.y, end.y, ny, endpoint: endpoint)
        return xs.flatMap { x in ys.map { SIMD2(x, $0) } }
    }
}

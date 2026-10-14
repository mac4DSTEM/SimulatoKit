//
//  GridScan.swift
//  Role: STEM scan positions (registration-slice-2.md section 3.3; abTEM `GridScan.get_positions`): along each axis
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
        self.start = start
        self.end = end
        self.nx = nx
        self.ny = ny
        self.endpoint = endpoint
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

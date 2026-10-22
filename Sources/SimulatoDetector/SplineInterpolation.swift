//
//  SplineInterpolation.swift
//  Role: scipy.ndimage's `map_coordinates` as abTEM calls it (orders 0 to 5, mode "constant", prefilter on): B-spline
//        coefficients by the recursive filter with mirror boundaries (scipy's `spline_filter1d`), then the tensor-product
//        B-spline sum at each coordinate; a coordinate outside [0, n − 1] is 0, an index beyond the edge is mirrored.
//        Used by the line profiles' `interpolate` and by `interpolate_line`.
//

import Foundation
import SimulatoCore

enum SplineInterpolation {
    static let orders = 0...5

    /// scipy's filter poles for an order (none below 2).
    static func poles(_ order: Int) -> [Double] {
        switch order {
        case 2: return [sqrt(8.0) - 3]
        case 3: return [sqrt(3.0) - 2]
        case 4: return [sqrt(664 - sqrt(438976.0)) + sqrt(304.0) - 19,
                        sqrt(664 + sqrt(438976.0)) - sqrt(304.0) - 19]
        case 5: return [sqrt(67.5 - sqrt(4436.25)) + sqrt(26.25) - 6.5,
                        sqrt(67.5 + sqrt(4436.25)) - sqrt(26.25) - 6.5]
        default: return []
        }
    }

    /// The coefficients of the line `c[base + k * stride]`, k < n, in place (scipy's `spline_filter1d`, mirror boundaries).
    static func filterLine(_ c: inout [Double], base: Int, step: Int, n: Int, poles: [Double]) {
        guard n > 1, !poles.isEmpty else { return }
        var gain = 1.0
        for z in poles { gain *= (1 - z) * (1 - 1 / z) }
        for k in 0..<n { c[base + k * step] *= gain }
        for z in poles {
            // causal initialisation (mirror)
            var zi = z
            let zn1 = pow(z, Double(n - 1))
            var c0 = c[base] + zn1 * c[base + (n - 1) * step]
            if n > 2 {
                for i in 1..<(n - 1) {
                    c0 += zi * (c[base + i * step] + zn1 * c[base + (n - 1 - i) * step])
                    zi *= z
                }
            }
            c[base] = c0 / (1 - zn1 * zn1)
            for i in 1..<n { c[base + i * step] += z * c[base + (i - 1) * step] }
            // anticausal initialisation (mirror)
            c[base + (n - 1) * step] = (z / (z * z - 1)) * (c[base + (n - 1) * step] + z * c[base + (n - 2) * step])
            for i in Swift.stride(from: n - 2, through: 0, by: -1) {
                c[base + i * step] = z * (c[base + (i + 1) * step] - c[base + i * step])
            }
        }
    }

    /// The coefficients of a 1D array.
    static func coefficients(_ a: [Double], order: Int) -> [Double] {
        var c = a
        filterLine(&c, base: 0, step: 1, n: a.count, poles: poles(order))
        return c
    }

    /// The coefficients of an nx × ny array (row-major), filtered along each axis.
    static func coefficients(_ a: [Double], nx: Int, ny: Int, order: Int) -> [Double] {
        var c = a
        let p = poles(order)
        guard !p.isEmpty else { return c }
        for iy in 0..<ny { filterLine(&c, base: iy, step: ny, n: nx, poles: p) }
        for ix in 0..<nx { filterLine(&c, base: ix * ny, step: 1, n: ny, poles: p) }
        return c
    }

    /// The cardinal B-spline of order n ≥ 1 at t: (1/n!) Σ (−1)^k C(n+1, k) (t + (n+1)/2 − k)₊ⁿ.
    static func bspline(_ n: Int, _ t: Double) -> Double {
        let half = Double(n + 1) / 2
        guard abs(t) < half else { return 0 }
        var sum = 0.0, binomial = 1.0
        for k in 0...(n + 1) {
            let u = t + half - Double(k)
            if u > 0 { sum += (k % 2 == 0 ? 1 : -1) * binomial * pow(u, Double(n)) }
            binomial = binomial * Double(n + 1 - k) / Double(k + 1)
        }
        var factorial = 1.0
        if n >= 2 { for i in 2...n { factorial *= Double(i) } }
        return sum / factorial
    }

    /// The first index and the weights of the order + 1 taps at coordinate `cc`.
    static func taps(_ order: Int, _ cc: Double) -> (first: Int, weights: [Double]) {
        if order == 0 { return (Int((cc + 0.5).rounded(.down)), [1]) }
        let centre = order % 2 == 1 ? cc.rounded(.down) : (cc + 0.5).rounded(.down)
        let first = Int(centre) - order / 2
        return (first, (0...order).map { bspline(order, cc - Double(first + $0)) })
    }

    /// The mirror image of an index outside 0 ..< n.
    @inline(__always) static func mirror(_ i: Int, _ n: Int) -> Int {
        if i >= 0 && i < n { return i }
        if n == 1 { return 0 }
        let period = 2 * (n - 1)
        var m = i % period
        if m < 0 { m += period }
        return m >= n ? period - m : m
    }

    /// The spline through the 1D coefficients at `cc` (0 outside [0, n − 1]).
    static func evaluate(_ c: [Double], at cc: Double, order: Int) -> Double {
        let n = c.count
        guard cc >= 0, cc <= Double(n - 1) else { return 0 }
        let t = taps(order, cc)
        var sum = 0.0
        for (k, w) in t.weights.enumerated() { sum += w * c[mirror(t.first + k, n)] }
        return sum
    }

    /// The spline through the 2D coefficients at (x, y) (0 outside the array).
    static func evaluate(_ c: [Double], nx: Int, ny: Int, x: Double, y: Double, order: Int) -> Double {
        guard x >= 0, x <= Double(nx - 1), y >= 0, y <= Double(ny - 1) else { return 0 }
        let tx = taps(order, x), ty = taps(order, y)
        var sum = 0.0
        for (i, wx) in tx.weights.enumerated() {
            let row = mirror(tx.first + i, nx) * ny
            var rowSum = 0.0
            for (j, wy) in ty.weights.enumerated() { rowSum += wy * c[row + mirror(ty.first + j, ny)] }
            sum += wx * rowSum
        }
        return sum
    }
}

//
//  LineInterpolation.swift
//  Role: abTEM 1.0.10's `Images.interpolate_line` and `interpolate_line_at_position` (measurements.py, `BaseMeasurements2D`
//        and `_interpolate_stack`): a line scan's positions in pixels, the image padded periodically by 2 × order on each
//        side, a spline of the given order at each position (`map_coordinates`, mode "constant"), averaged across a width
//        when one is given. On a `RealField` of an image (Å) or of a diffraction pattern (1/Å, zero frequency at index
//        n / 2, as `Waves.diffractionPattern` and `PixelatedDetector` return it).
//

import Foundation
import SimulatoCore
import SimulatoProbe

extension RealField {
    /// What abTEM calls the measurement's `offset` and `extent` for an image, or for a centred diffraction pattern.
    struct Frame {
        let offset: SIMD2<Double>
        let extent: SIMD2<Double>
    }

    /// An image: the origin at the first pixel, the extent sampling × points.
    var imageFrame: Frame {
        Frame(offset: .zero, extent: SIMD2(samplingX * Double(nx), samplingY * Double(ny)))
    }

    /// A centred pattern: the first pixel at −⌊n / 2⌋ × sampling, and abTEM's `extent` (the lowest limit minus the highest,
    /// which is −(n − 1) × sampling).
    var patternFrame: Frame {
        Frame(offset: SIMD2(-Double(nx / 2) * samplingX, -Double(ny / 2) * samplingY),
              extent: SIMD2(-Double(nx - 1) * samplingX, -Double(ny - 1) * samplingY))
    }

    /// abTEM's `interpolate_line`: the profile values, the line's sampling and the line itself.
    func interpolatedLine(frame: Frame, start: SIMD2<Double>, end: SIMD2<Double>, sampling: Double?, gpts: Int?, width: Double,
                          margin: Double, order: Int, endpoint: Bool, fractional: Bool) throws
        -> (values: [Double], line: LineScan, width: Double?) {
        guard SplineInterpolation.orders.contains(order) else {
            throw SimulatoError.invalidScan("the spline order must be 0 to 5, got \(order)")
        }
        guard width.isFinite, width >= 0 else { throw SimulatoError.invalidScan("the width must be finite and not negative, got \(width)") }
        guard margin.isFinite else { throw SimulatoError.invalidScan("the margin must be finite") }
        let minSampling = Swift.min(samplingX, samplingY)
        let s = (sampling == nil && gpts == nil) ? minSampling : sampling
        var scan = try LineScan(start: start, end: end, gpts: gpts, sampling: s, endpoint: endpoint, fractional: fractional,
                                extent: fractional ? frame.extent : nil)
        if margin != 0 { scan = try scan.addingMargin(margin) }
        let pixels = scan.positions.map { SIMD2(($0.x - frame.offset.x) / samplingX, ($0.y - frame.offset.y) / samplingY) }

        // The image padded periodically by 2 * order, the spline coefficients of the padded array.
        let pad = 2 * order
        let px = nx + 2 * pad, py = ny + 2 * pad
        var padded = [Double](repeating: 0, count: px * py)
        for i in 0..<px {
            var ix = (i - pad) % nx
            if ix < 0 { ix += nx }
            for j in 0..<py {
                var iy = (j - pad) % ny
                if iy < 0 { iy += ny }
                padded[i * py + j] = values[ix * ny + iy]
            }
        }
        let c = SplineInterpolation.coefficients(padded, nx: px, ny: py, order: order)
        let shift = Double(pad)

        guard width > 0 else {
            return (pixels.map { SplineInterpolation.evaluate(c, nx: px, ny: py, x: $0.x + shift, y: $0.y + shift, order: order) },
                    scan, nil)
        }
        let length = scan.extent
        guard length > 0 else { throw SimulatoError.invalidScan("a width needs a line of positive length") }
        let direction = (scan.end - scan.start) / length
        let perpendicular = SIMD2(-direction.y, direction.x)
        let n = (width / minSampling / 2).rounded(.down) * 2 + 1
        let across = GridScan.linspace(-n / 2, n / 2, Int(n), endpoint: true)
        var out = [Double](repeating: 0, count: pixels.count)
        for (k, p) in pixels.enumerated() {
            var sum = 0.0
            for a in across {
                let q = SIMD2(a * perpendicular.x + p.x, a * perpendicular.y + p.y)
                sum += SplineInterpolation.evaluate(c, nx: px, ny: py, x: q.x + shift, y: q.y + shift, order: order)
            }
            out[k] = sum / Double(across.count)
        }
        return (out, scan, width)
    }

    /// An image interpolated along a line (abTEM `Images.interpolate_line`): the real-space profile from `start` to `end`
    /// [Å] at `gpts` points or at `sampling` [Å] (neither: the finer pixel size), averaged across `width` [Å] when given.
    /// - Parameters:
    ///   - start: the first point [Å] (default the origin), or a fraction of the image extent with `fractional`.
    ///   - end: the last point [Å] (abTEM's default, `(0, extent[0])`, is not offered).
    ///   - margin: added to both ends of the line [Å].
    ///   - order: the spline order, 0 to 5 (default 3); positions outside the image wrap around it.
    ///   - endpoint: false (default) leaves `end` out, as abTEM.
    /// - Throws: `invalidScan` for an invalid line (see `LineScan`), an order outside 0 to 5, or a negative width.
    public func interpolateLine(start: SIMD2<Double> = .zero, end: SIMD2<Double>, sampling: Double? = nil, gpts: Int? = nil,
                                width: Double = 0, margin: Double = 0, order: Int = 3, endpoint: Bool = false,
                                fractional: Bool = false) throws -> RealSpaceLineProfiles {
        let r = try interpolatedLine(frame: imageFrame, start: start, end: end, sampling: sampling, gpts: gpts, width: width,
                                     margin: margin, order: order, endpoint: endpoint, fractional: fractional)
        return RealSpaceLineProfiles(count: r.values.count, values: r.values, sampling: r.line.sampling, start: r.line.start,
                                     end: r.line.end, lineWidth: r.width)
    }

    /// An image interpolated along a line centred on `center` (abTEM `interpolate_line_at_position`): a line of `extent`
    /// [Å] at `angle` [deg] to the first axis; `endpoint` defaults to true here, as in abTEM.
    public func interpolateLineAtPosition(center: SIMD2<Double>, angle: Double, extent: Double, gpts: Int? = nil,
                                          sampling: Double? = nil, width: Double = 0, order: Int = 3,
                                          endpoint: Bool = true) throws -> RealSpaceLineProfiles {
        let p = LineScan.endpoints(center: center, extent: extent, angle: angle)
        return try interpolateLine(start: p.start, end: p.end, sampling: sampling, gpts: gpts, width: width, order: order,
                                   endpoint: endpoint)
    }

    /// A centred diffraction pattern (sampling in 1/Å, the zero frequency at index n / 2 of each axis, as
    /// `Waves.diffractionPattern` returns it) interpolated along a line in reciprocal space (abTEM
    /// `DiffractionPatterns.interpolate_line`): `start` and `end` are spatial frequencies [1/Å] and the result is a
    /// reciprocal-space profile. With `fractional`, abTEM's own pattern extent applies (it is negative: −(n − 1) × sampling).
    public func interpolateReciprocalLine(start: SIMD2<Double> = .zero, end: SIMD2<Double>, sampling: Double? = nil,
                                          gpts: Int? = nil, width: Double = 0, margin: Double = 0, order: Int = 3,
                                          endpoint: Bool = false, fractional: Bool = false) throws -> ReciprocalSpaceLineProfiles {
        let r = try interpolatedLine(frame: patternFrame, start: start, end: end, sampling: sampling, gpts: gpts, width: width,
                                     margin: margin, order: order, endpoint: endpoint, fractional: fractional)
        return ReciprocalSpaceLineProfiles(count: r.values.count, values: r.values, sampling: r.line.sampling,
                                           start: r.line.start, end: r.line.end, lineWidth: r.width)
    }

    /// The same for a line centred on `center` [1/Å] (abTEM `interpolate_line_at_position` on a diffraction pattern).
    public func interpolateReciprocalLineAtPosition(center: SIMD2<Double>, angle: Double, extent: Double, gpts: Int? = nil,
                                                    sampling: Double? = nil, width: Double = 0, order: Int = 3,
                                                    endpoint: Bool = true) throws -> ReciprocalSpaceLineProfiles {
        let p = LineScan.endpoints(center: center, extent: extent, angle: angle)
        return try interpolateReciprocalLine(start: p.start, end: p.end, sampling: sampling, gpts: gpts, width: width,
                                             order: order, endpoint: endpoint)
    }
}

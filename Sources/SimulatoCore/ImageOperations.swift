//
//  ImageOperations.swift
//  Role: what abTEM does to a finished image (`Images`): `integrate_gradient` (the Fourier-space integration of a
//        two-component gradient field, `_integrate_gradient_2d`), `integrate_disc` and `tile`. CPU float64 on the kit's own
//        FFT. Rows and columns: x is the first index, arrays are row-major `[ix * ny + iy]`.
//

import Foundation

/// A refusal by one of the image operations (`SimulatoError` covers the grid lengths and samplings).
public enum ImageOperationError: Error, Sendable, Equatable, CustomStringConvertible {
    /// The two components of a gradient differ in shape or sampling.
    case gradientComponentsDiffer
    /// A repetition count below 1, or a tiled image of more than 2^28 values (or one whose size overflows).
    case invalidRepetitions(Int, Int)
    /// A disc position or radius that is not finite, or a radius that is not positive.
    case invalidDisc(String)
    /// `Border.raise` and the disc crosses the image border.
    case discOutsideImage

    public var description: String {
        switch self {
        case .gradientComponentsDiffer: return "the two gradient components must have the same shape and sampling"
        case .invalidRepetitions(let x, let y): return "repetitions must be at least 1 along both axes and give at most 2^28 values, got (\(x), \(y))"
        case .invalidDisc(let why): return "invalid integration disc: \(why)"
        case .discOutsideImage: return "the integration region is outside the image"
        }
    }
}

extension RealField {
    /// The integral of a gradient field (abTEM `Images.integrate_gradient` of an image whose real part is the gradient's x
    /// component and imaginary part its y component): T̂ = (k_x ĝ_x + k_y ĝ_y) / (2πi |k|²) with NumPy's `fftfreq`
    /// frequencies (the Nyquist frequency of an even length is negative), T̂ = 0 at zero frequency, no padding; the real
    /// part of the inverse transform, shifted so its minimum is 0. The result has the sampling of the components.
    /// - Throws: `ImageOperationError.gradientComponentsDiffer` when the two differ in shape or sampling,
    ///   `SimulatoError.invalidSampling` for a non-finite or non-positive sampling, `unsupportedGridLength` for a length
    ///   below 2.
    public static func integratingGradient(x gx: RealField, y gy: RealField) throws -> RealField {
        guard gx.nx == gy.nx, gx.ny == gy.ny, gx.samplingX == gy.samplingX, gx.samplingY == gy.samplingY else {
            throw ImageOperationError.gradientComponentsDiffer
        }
        let values = try integrateGradient(gx: gx.values, gy: gy.values, nx: gx.nx, ny: gx.ny,
                                           samplingX: gx.samplingX, samplingY: gx.samplingY)
        return RealField(nx: gx.nx, ny: gx.ny, values: values, samplingX: gx.samplingX, samplingY: gx.samplingY)
    }

    /// The tiled image (abTEM `Images.tile`): `repetitionsX` × `repetitionsY` copies; the sampling stays.
    /// - Throws: `ImageOperationError.invalidRepetitions` for a count below 1 or a result of more than 2^28 values.
    public func tiled(repetitionsX: Int, repetitionsY: Int) throws -> RealField {
        guard repetitionsX >= 1, repetitionsY >= 1 else { throw ImageOperationError.invalidRepetitions(repetitionsX, repetitionsY) }
        let (mx, o1) = nx.multipliedReportingOverflow(by: repetitionsX)
        let (my, o2) = ny.multipliedReportingOverflow(by: repetitionsY)
        let (total, o3) = mx.multipliedReportingOverflow(by: my)
        guard !o1, !o2, !o3, total <= 1 << 28 else { throw ImageOperationError.invalidRepetitions(repetitionsX, repetitionsY) }
        var out = [Double](repeating: 0, count: mx * my)
        for ix in 0..<mx {
            let row = (ix % nx) * ny
            for iy in 0..<my { out[ix * my + iy] = values[row + iy % ny] }
        }
        return RealField(nx: mx, ny: my, values: out, samplingX: samplingX, samplingY: samplingY)
    }

    /// How `integrateDisc` treats a disc that crosses the image border (abTEM's `border`).
    public enum Border: Sendable, Hashable {
        /// The image is extended by wrapping around to the opposite edge.
        case wrap
        /// A disc crossing the border is refused.
        case raise
    }

    /// The sum (or mean) of the image over a disc of `radius` [Å] centred at (`x`, `y`) [Å] from the image's origin
    /// (abTEM `integrate_disc`, no `interpolate`): the pixels of the bounding box with a soft edge, weight
    /// 1 − clip((r − radius) / mean sampling, 0, 1). The distances are abTEM's, which place the box's pixels on
    /// linspace(0, m d, m) with the endpoint (spacing m d / (m − 1), not d); ported as they are.
    /// With `mean` and a disc whose pixels all have weight 0 the result is NaN, as abTEM's 0 / 0.
    /// - Throws: `ImageOperationError.invalidDisc` for a non-finite position or a radius that is not finite and positive,
    ///   or a disc or position of more than 2^31 pixels across or 2^53 pixels from the origin,
    ///   `discOutsideImage` for `.raise` at the border.
    public func integrateDisc(x: Double, y: Double, radius: Double, mean: Bool = false, border: Border = .wrap) throws -> Double {
        guard x.isFinite, y.isFinite else { throw ImageOperationError.invalidDisc("the position must be finite, got (\(x), \(y)) Å") }
        guard radius.isFinite, radius > 0 else { throw ImageOperationError.invalidDisc("the radius must be finite and positive, got \(radius) Å") }
        guard samplingX.isFinite, samplingY.isFinite, samplingX > 0, samplingY > 0 else {
            throw SimulatoError.invalidSampling(samplingX, samplingY)
        }
        // Int(...) of a huge Double traps: the box must fit 2^31 pixels per axis and the corner 2^53 pixels from the origin.
        guard 2 * radius / samplingX < 0x1p31, 2 * radius / samplingY < 0x1p31,
              abs(x / samplingX) < 0x1p53, abs(y / samplingY) < 0x1p53 else {
            throw ImageOperationError.invalidDisc("the disc of radius \(radius) Å at (\(x), \(y)) Å is too large for the sampling (\(samplingX), \(samplingY)) Å")
        }
        let mx = Int((2 * radius / samplingX).rounded(.up)), my = Int((2 * radius / samplingY).rounded(.up))
        let cornerX = Int((x / samplingX).rounded(.down)) - mx / 2
        let cornerY = Int((y / samplingY).rounded(.down)) - my / 2
        if border == .raise, cornerX < 0 || cornerY < 0 || cornerX + mx > nx || cornerY + my > ny {
            throw ImageOperationError.discOutsideImage
        }
        // abTEM's pixel coordinates inside the box are numpy's linspace(0, m * d, m) with the endpoint (an Images axis has
        // endpoint=True), so the spacing is m d / (m − 1), a little more than d, and the last one is exactly m d
        func coordinate(_ i: Int, _ m: Int, _ d: Double) -> Double {
            if m == 1 { return 0 }
            return i == m - 1 ? Double(m) * d : Double(i) * ((Double(m) * d) / Double(m - 1))
        }
        let px = x - Double(cornerX) * samplingX, py = y - Double(cornerY) * samplingY
        let meanSampling = (samplingX + samplingY) / 2
        var total = 0.0, weights = 0.0
        for i in 0..<mx {
            let ix = ((cornerX + i) % nx + nx) % nx
            let dx = coordinate(i, mx, samplingX) - px
            for j in 0..<my {
                let iy = ((cornerY + j) % ny + ny) % ny
                let dy = coordinate(j, my, samplingY) - py
                let r = (dx * dx + dy * dy).squareRoot()
                let w = 1 - Swift.min(Swift.max((r - radius) / meanSampling, 0), 1)
                total += values[ix * ny + iy] * w
                weights += w
            }
        }
        return mean ? total / weights : total
    }
}

extension ComplexField {
    /// `RealField.integratingGradient` for a field whose real part is the gradient's x component and imaginary part its y
    /// component (abTEM's complex `Images`), with the sampling [Å] of its pixels.
    public func integratedGradient(samplingX: Double, samplingY: Double) throws -> RealField {
        let values = try integrateGradient(gx: real, gy: imaginary, nx: nx, ny: ny, samplingX: samplingX, samplingY: samplingY)
        return RealField(nx: nx, ny: ny, values: values, samplingX: samplingX, samplingY: samplingY)
    }
}

/// abTEM `_integrate_gradient_2d` in float64.
func integrateGradient(gx: [Double], gy: [Double], nx: Int, ny: Int, samplingX: Double, samplingY: Double) throws -> [Double] {
    guard samplingX.isFinite, samplingY.isFinite, samplingX > 0, samplingY > 0 else {
        throw SimulatoError.invalidSampling(samplingX, samplingY)
    }
    let fft = try FFT2D(nx: nx, ny: ny)
    var fx = ComplexField(nx: nx, ny: ny, real: gx, imaginary: [Double](repeating: 0, count: nx * ny))
    var fy = ComplexField(nx: nx, ny: ny, real: gy, imaginary: [Double](repeating: 0, count: nx * ny))
    fft.forward(&fx)
    fft.forward(&fy)
    let kx = fftFrequencies(nx, spacing: samplingX), ky = fftFrequencies(ny, spacing: samplingY)
    var that = ComplexField(nx: nx, ny: ny)
    for ix in 0..<nx {
        for iy in 0..<ny {
            let i = ix * ny + iy
            var k = kx[ix] * kx[ix] + ky[iy] * ky[iy]
            if k == 0 { k = 1e-12 }
            // (a + ib) / (2πi k) = (b − ia) / (2π k)
            let a = fx.real[i] * kx[ix] + fy.real[i] * ky[iy]
            let b = fx.imaginary[i] * kx[ix] + fy.imaginary[i] * ky[iy]
            let c = 2 * Double.pi * k
            that.real[i] = b / c
            that.imaginary[i] = -a / c
        }
    }
    fft.inverse(&that)
    var out = that.real
    let low = out.min() ?? 0
    for i in out.indices { out[i] -= low }
    return out
}

/// NumPy's `fftfreq(n, d)`: k/(n d) for k = 0 ... ⌈n/2⌉ − 1, then −⌊n/2⌋ ... −1.
func fftFrequencies(_ n: Int, spacing d: Double) -> [Double] {
    let positive = (n - 1) / 2 + 1
    let scale = 1.0 / (Double(n) * d)
    return (0..<n).map { i in Double(i < positive ? i : i - n) * scale }
}

//
//  Fields.swift
//  Role: the arrays the kit returns: a complex field (split real and imaginary parts, vDSP's layout) and a real field,
//        both on a grid, x the first index, row-major `[ix * ny + iy]` (abTEM's `[x, y]`).
//

/// A complex array on an nx × ny grid, split into real and imaginary parts.
public struct ComplexField: Sendable, Hashable {
    public let nx: Int
    public let ny: Int
    public var real: [Double]
    public var imaginary: [Double]

    /// A field of zeros.
    public init(nx: Int, ny: Int) {
        self.nx = nx
        self.ny = ny
        real = [Double](repeating: 0, count: nx * ny)
        imaginary = [Double](repeating: 0, count: nx * ny)
    }

    /// - Precondition: both arrays hold nx × ny values.
    public init(nx: Int, ny: Int, real: [Double], imaginary: [Double]) {
        precondition(real.count == nx * ny && imaginary.count == nx * ny, "ComplexField: \(nx) x \(ny) values required")
        self.nx = nx
        self.ny = ny
        self.real = real
        self.imaginary = imaginary
    }

    /// The value at (ix, iy).
    public subscript(ix: Int, iy: Int) -> (re: Double, im: Double) {
        (real[ix * ny + iy], imaginary[ix * ny + iy])
    }

    /// |ψ|² per pixel (abTEM's `abs2`: re² + im²).
    public func squaredMagnitude() -> [Double] {
        var out = [Double](repeating: 0, count: real.count)
        for i in out.indices { out[i] = real[i] * real[i] + imaginary[i] * imaginary[i] }
        return out
    }

    /// The largest |ψ|.
    public var maxMagnitude: Double {
        var m = 0.0
        for i in real.indices { m = Swift.max(m, (real[i] * real[i] + imaginary[i] * imaginary[i]).squareRoot()) }
        return m
    }
}

/// A real array on an nx × ny grid with its sampling (Å per pixel for images, 1/Å per pixel for patterns).
public struct RealField: Sendable, Hashable {
    public let nx: Int
    public let ny: Int
    public var values: [Double]
    public let samplingX: Double
    public let samplingY: Double

    public init(nx: Int, ny: Int, values: [Double], samplingX: Double, samplingY: Double) {
        precondition(values.count == nx * ny, "RealField: \(nx) x \(ny) values required")
        self.nx = nx
        self.ny = ny
        self.values = values
        self.samplingX = samplingX
        self.samplingY = samplingY
    }

    public subscript(ix: Int, iy: Int) -> Double { values[ix * ny + iy] }

    public var sum: Double { values.reduce(0, +) }
    public var max: Double { values.max() ?? 0 }
}

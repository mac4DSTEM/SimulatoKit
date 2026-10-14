//
//  Grid.swift
//  Role: the lateral simulation grid (gpts, extent, sampling), its spatial frequencies (NumPy's fftfreq, bit for bit)
//        and the FFT lengths both backends accept (registration section 3.9).
//

/// A lateral grid of nx × ny points over an extent in Å; x is the first index and arrays are row-major `[ix * ny + iy]`.
public struct Grid: Sendable, Hashable {
    public let nx: Int
    public let ny: Int
    /// Extent along x and y [Å].
    public let extentX: Double
    public let extentY: Double

    /// - Throws: `SimulatoError.unsupportedGridLength` for a length outside `isSupportedLength`, `invalidExtent` for a
    ///   non-finite or non-positive extent.
    public init(nx: Int, ny: Int, extentX: Double, extentY: Double) throws {
        guard Self.isSupportedLength(nx) else { throw SimulatoError.unsupportedGridLength(nx) }
        guard Self.isSupportedLength(ny) else { throw SimulatoError.unsupportedGridLength(ny) }
        guard extentX.isFinite, extentY.isFinite, extentX > 0, extentY > 0 else {
            throw SimulatoError.invalidExtent(extentX, extentY)
        }
        self.nx = nx
        self.ny = ny
        self.extentX = extentX
        self.extentY = extentY
    }

    /// Sampling along x and y [Å per pixel] (extent / gpts, as abTEM's `Grid`).
    public var samplingX: Double { extentX / Double(nx) }
    public var samplingY: Double { extentY / Double(ny) }
    public var count: Int { nx * ny }

    /// Reciprocal-space sampling [1/Å]: 1 / (n · (extent / n)), as abTEM computes it (`core/grid.py`
    /// `reciprocal_space_sampling`); on some grids this differs from 1 / extent by one ulp, enough to move a pixel across a
    /// detector limit placed at a whole number of pixels (slice-2 refuter, finding B1).
    public var reciprocalSamplingX: Double { 1 / (Double(nx) * samplingX) }
    public var reciprocalSamplingY: Double { 1 / (Double(ny) * samplingY) }

    /// The lengths both backends accept: f · 2^m with f ∈ {1, 3, 5, 15} and m ≥ 3 (vDSP's DFT lengths).
    public static func isSupportedLength(_ n: Int) -> Bool {
        guard n >= 8 else { return false }
        var m = n
        var power = 0
        while m % 2 == 0 { m /= 2; power += 1 }
        return [1, 3, 5, 15].contains(m) && power >= 3
    }

    /// The same grid within rounding: equal point counts and extents within 1e−12 relative (exact `==` refuses a grid
    /// whose extent was computed along another path, e.g. 9.600000000000001 against 9.6; refuter finding, 2026-10-15).
    public func matches(_ other: Grid) -> Bool {
        nx == other.nx && ny == other.ny
            && abs(extentX - other.extentX) <= 1e-12 * max(abs(extentX), abs(other.extentX))
            && abs(extentY - other.extentY) <= 1e-12 * max(abs(extentY), abs(other.extentY))
    }

    /// NumPy's `fftfreq(n, d)`: integers [0, 1, ..., ⌈n/2⌉ − 1, −⌊n/2⌋, ..., −1] times 1 / (n d), bit for bit.
    public static func fftfreq(_ n: Int, _ d: Double) -> [Double] {
        let val = 1.0 / (Double(n) * d)
        let positive = (n - 1) / 2 + 1
        return (0..<n).map { i in Double(i < positive ? i : i - n) * val }
    }

    /// Spatial frequencies kx (length nx) and ky (length ny) [1/Å].
    public var frequenciesX: [Double] { Self.fftfreq(nx, samplingX) }
    public var frequenciesY: [Double] { Self.fftfreq(ny, samplingY) }
}

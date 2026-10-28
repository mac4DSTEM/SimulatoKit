//
//  Antialias.swift
//  Role: abTEM's antialias aperture (`antialias.py`, configuration cutoff 0.6666666, taper 0.01): a cosine-tapered disk
//        at two thirds of the Nyquist frequency of the coarser axis, and the size of the diffraction pattern cropped to
//        it.
//

import Foundation

/// The bandlimit every multislice step applies.
public enum Antialias {
    /// The cutoff as a fraction of the sampling frequency 1/Δ (abTEM's `antialias.cutoff`).
    public static let cutoff = 0.6666666
    /// The taper width as a fraction of 1/Δ (abTEM's `antialias.taper`).
    public static let taper = 0.01

    /// The aperture on the grid's spatial frequencies: 1 inside k_c − w, a cosine taper to 0 at k_c, 0 outside.
    public static func aperture(on grid: Grid) -> [Double] {
        let dmax = max(grid.samplingX, grid.samplingY)
        let kc = cutoff / dmax / 2
        let w = taper / dmax
        let kx = grid.frequenciesX, ky = grid.frequenciesY
        var out = [Double](repeating: 0, count: grid.count)
        for ix in 0..<grid.nx {
            for iy in 0..<grid.ny {
                let r = (kx[ix] * kx[ix] + ky[iy] * ky[iy]).squareRoot()
                let value: Double
                if r > kc {
                    value = 0
                } else if r > kc - w {
                    value = 0.5 * (1 + cos(Double.pi * (r - kc + w) / w))
                } else {
                    value = 1
                }
                out[ix * grid.ny + iy] = value
            }
        }
        return out
    }

    /// The pattern size at the antialias cutoff (abTEM's `diffraction_patterns(max_angle="cutoff")`, parity "odd"):
    /// ⌊round((2/3) / max(Δ) × L, 7)⌋, moved by +1 to the grid's parity, then by +1 if even.
    public static func cutoffPatternSize(on grid: Grid) -> (nx: Int, ny: Int) {
        let kcut = 2.0 / 3.0 / max(grid.samplingX, grid.samplingY)
        func size(_ n: Int, _ extent: Double) -> Int {
            var m = Int(floor(roundToDecimals(kcut * extent, 7)))
            if m % 2 != n % 2 { m += 1 }
            if m % 2 == 0 { m += 1 }
            return m
        }
        return (size(grid.nx, Double(grid.nx) * grid.samplingX), size(grid.ny, Double(grid.ny) * grid.samplingY))
    }

    /// NumPy's `np.round(x, decimals)` for a scalar: x · 10^d rounded half to even, divided by 10^d.
    static func roundToDecimals(_ x: Double, _ decimals: Int) -> Double {
        let f = pow(10.0, Double(decimals))
        return (x * f).rounded(.toNearestOrEven) / f
    }
}

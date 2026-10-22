//
//  WavesResample.swift
//  Role: abTEM's `fft_interpolate` for a wave: the Fourier-space crop or zero-pad to a new number of points over the same
//        extent (`core/fft.py`: `fft_interpolation_masks`, `fft_crop`, `fft_interpolate` with normalization "values").
//

import Foundation

extension Waves {
    /// The wave resampled to `rows` × `columns` points over the same extent by Fourier-space cropping (smaller) or
    /// zero-padding (larger), as abTEM's `fft_interpolate(..., normalization="values")` does for `WavesDetector(gpts)`:
    /// FFT, keep the lowest frequencies (the first ⌊m/2⌋ + 1 and last ⌊m/2⌋ of an odd axis, m/2 and m/2 of an even one),
    /// inverse FFT, scaled so the amplitude of the values is kept. The beam, normalisation and tilt are kept.
    /// - Throws: `SimulatoError.unsupportedGridLength` for fewer than 2 points per axis, `SimulatoError.gridTooLarge` when
    ///   the result would not fit in half of the physical memory (or its size overflows).
    public func resampled(rows: Int, columns: Int) throws -> Waves {
        guard Grid.isSupportedLength(rows) else { throw SimulatoError.unsupportedGridLength(rows) }
        guard Grid.isSupportedLength(columns) else { throw SimulatoError.unsupportedGridLength(columns) }
        let limit = Int(ProcessInfo.processInfo.physicalMemory / 2)
        let (points, o1) = rows.multipliedReportingOverflow(by: columns)
        let (bytes, o2) = points.multipliedReportingOverflow(by: 16)
        guard !o1, !o2, bytes <= limit else {
            throw SimulatoError.gridTooLarge(rows: rows, columns: columns, required: o1 || o2 ? .max : bytes, limit: limit)
        }
        let target = try Grid(nx: rows, ny: columns, extentX: grid.extentX, extentY: grid.extentY)
        var f = field
        try FFT2D(nx: grid.nx, ny: grid.ny).forward(&f)
        let mapX = Self.resampleMap(from: grid.nx, to: rows)
        let mapY = Self.resampleMap(from: grid.ny, to: columns)
        var out = ComplexField(nx: rows, ny: columns)
        for (sx, dx) in mapX {
            for (sy, dy) in mapY {
                out.real[dx * columns + dy] = f.real[sx * grid.ny + sy]
                out.imaginary[dx * columns + dy] = f.imaginary[sx * grid.ny + sy]
            }
        }
        try FFT2D(nx: rows, ny: columns).inverse(&out)
        // ifft2 divides by the new point count, abTEM then multiplies by new / old.
        let scale = Double(rows * columns) / Double(grid.count)
        for i in out.real.indices { out.real[i] *= scale; out.imaginary[i] *= scale }
        return Waves(field: out, grid: target, beam: beam, normalization: normalization, tilt: tilt)
    }

    /// (source index, destination index) of every frequency that survives going from n to m points (abTEM's
    /// `_fft_interpolation_masks_1d`): the lowest m (crop) or all n (pad), positive frequencies first.
    package static func resampleMap(from n: Int, to m: Int) -> [(Int, Int)] {
        if m <= n {
            if m == 1 { return [(0, 0)] }
            return cropIndices(n: n, m: m).enumerated().map { ($1, $0) }
        }
        if n == 1 { return [(0, 0)] }
        let positive = n % 2 == 0 ? n / 2 : n / 2 + 1
        return (0..<n).map { i in (i, i < positive ? i : m - n + i) }
    }
}

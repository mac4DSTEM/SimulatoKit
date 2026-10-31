//
//  ImageInterpolation.swift
//  Role: abTEM 1.0.10's `Images.interpolate` with its default method "fft" (measurements.py, and `fft_interpolate` in
//        core/fft.py): the image's spectrum cropped or zero-padded to the new number of points over the same extent, the
//        real part of the inverse transform, scaled by new / old points for normalization "values". The "spline" method
//        is not offered. CPU float64 on the kit's own FFT.
//

import Foundation

extension RealField {
    /// What `interpolated` keeps (abTEM's `normalization`).
    public enum InterpolationNormalization: Sendable, Hashable {
        /// The pixel values (abTEM "values", the default): the result is scaled by new points / old points.
        case values
        /// The sum over the image (abTEM "intensity" or "amplitude"): no scaling.
        case intensity
    }

    /// The image resampled to `rows` × `columns` points over the same extent (abTEM `Images.interpolate(gpts=...)`,
    /// method "fft", periodic boundary): FFT, keep the lowest frequencies (crop) or pad with zeros, inverse FFT, the real
    /// part. The sampling becomes extent / points along each axis.
    /// - Throws: `SimulatoError.unsupportedGridLength` for fewer than 2 points on an axis of the image or of the result
    ///   (a departure: abTEM keeps the zero frequency of a one-point axis; the kit's FFT and grids start at 2 points),
    ///   `SimulatoError.gridTooLarge` when the result would not fit in half of the physical memory.
    public func interpolated(rows: Int, columns: Int, normalization: InterpolationNormalization = .values) throws -> RealField {
        guard Grid.isSupportedLength(nx) else { throw SimulatoError.unsupportedGridLength(nx) }
        guard Grid.isSupportedLength(ny) else { throw SimulatoError.unsupportedGridLength(ny) }
        guard Grid.isSupportedLength(rows) else { throw SimulatoError.unsupportedGridLength(rows) }
        guard Grid.isSupportedLength(columns) else { throw SimulatoError.unsupportedGridLength(columns) }
        let limit = Int(ProcessInfo.processInfo.physicalMemory / 2)
        let (points, o1) = rows.multipliedReportingOverflow(by: columns)
        let (bytes, o2) = points.multipliedReportingOverflow(by: 16)
        guard !o1, !o2, bytes <= limit else {
            throw SimulatoError.gridTooLarge(rows: rows, columns: columns, required: o1 || o2 ? .max : bytes, limit: limit)
        }
        var f = ComplexField(nx: nx, ny: ny, real: values, imaginary: [Double](repeating: 0, count: nx * ny))
        try FFT2D(nx: nx, ny: ny).forward(&f)
        let mapX = Waves.resampleMap(from: nx, to: rows), mapY = Waves.resampleMap(from: ny, to: columns)
        var out = ComplexField(nx: rows, ny: columns)
        for (sx, dx) in mapX {
            for (sy, dy) in mapY {
                out.real[dx * columns + dy] = f.real[sx * ny + sy]
                out.imaginary[dx * columns + dy] = f.imaginary[sx * ny + sy]
            }
        }
        try FFT2D(nx: rows, ny: columns).inverse(&out)
        var result = out.real
        if normalization == .values {
            let scale = Double(points) / Double(nx * ny)
            for i in result.indices { result[i] *= scale }
        }
        return RealField(nx: rows, ny: columns, values: result, samplingX: samplingX * Double(nx) / Double(rows),
                         samplingY: samplingY * Double(ny) / Double(columns))
    }

    /// The image resampled to a sampling of about `samplingX` × `samplingY` [Å] (abTEM `Images.interpolate(sampling=...)`,
    /// method "fft"): ⌈extent / sampling⌉ points along each axis, then as `interpolated(rows:columns:normalization:)`;
    /// the realised sampling is extent / points, at most the one asked for.
    /// - Throws: `SimulatoError.invalidSampling` for a sampling that is not finite and positive, or one giving 2^31 points
    ///   or more; otherwise as `interpolated(rows:columns:normalization:)`.
    public func interpolated(samplingX newX: Double, samplingY newY: Double,
                             normalization: InterpolationNormalization = .values) throws -> RealField {
        guard newX.isFinite, newY.isFinite, newX > 0, newY > 0 else { throw SimulatoError.invalidSampling(newX, newY) }
        // abTEM: int(ceil(extent / sampling)) with extent = sampling × points.
        let rows = (samplingX * Double(nx) / newX).rounded(.up), columns = (samplingY * Double(ny) / newY).rounded(.up)
        guard rows < 0x1p31, columns < 0x1p31 else { throw SimulatoError.invalidSampling(newX, newY) }
        return try interpolated(rows: Int(rows), columns: Int(columns), normalization: normalization)
    }
}

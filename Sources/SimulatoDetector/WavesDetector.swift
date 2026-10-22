//
//  WavesDetector.swift
//  Role: abTEM's `WavesDetector(gpts)`: the complex exit wave at every scan position (ScanWaves), optionally resampled in
//        Fourier space. One potential only (see `WavesDetectorError.severalPotentials`); a refusal before any work when the
//        waves would not fit in half of the physical memory.
//

import Foundation
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe

/// A refusal by the waves detector.
public enum WavesDetectorError: Error, Sendable, Equatable, CustomStringConvertible {
    /// The exit waves of the scan need `required` bytes, more than the `limit` allowed (half of the physical memory by
    /// default).
    case tooLarge(positions: Int, rows: Int, columns: Int, required: Int, limit: Int)
    /// More than one potential (frozen-phonon configuration) was given.
    case severalPotentials(Int)

    public var description: String {
        switch self {
        case .tooLarge(let p, let r, let c, let need, let limit):
            return "\(p) exit waves of \(r) x \(c) points need \(need) bytes, more than the \(limit) bytes allowed (half of the physical memory); scan fewer positions or detect on a smaller `rows` x `columns`"
        case .severalPotentials(let n):
            return "the waves detector takes one potential, got \(n): averaging complex waves over configurations is not meaningful; scan each configuration on its own"
        }
    }
}

/// Detects the complex exit wave (abTEM's `WavesDetector`). With `rows` and `columns` the wave is resampled in Fourier
/// space (abTEM's `gpts`), otherwise it keeps the simulation grid.
public struct WavesDetector: Sendable, Hashable {
    public let rows: Int?
    public let columns: Int?

    /// The exit wave on the simulation grid.
    public init() { rows = nil; columns = nil }

    /// The exit wave resampled to `rows` × `columns` points over the same extent (abTEM's `gpts`).
    /// - Throws: `SimulatoError.unsupportedGridLength` for fewer than 2 points per axis.
    public init(rows: Int, columns: Int) throws {
        guard Grid.isSupportedLength(rows) else { throw SimulatoError.unsupportedGridLength(rows) }
        guard Grid.isSupportedLength(columns) else { throw SimulatoError.unsupportedGridLength(columns) }
        self.rows = rows
        self.columns = columns
    }

    /// The bytes of `positions` complex float64 waves of `rows` × `columns` points (16 bytes a point), saturating.
    public static func requiredBytes(positions: Int, rows: Int, columns: Int) -> Int {
        let (a, o1) = positions.multipliedReportingOverflow(by: rows)
        let (b, o2) = a.multipliedReportingOverflow(by: columns)
        let (c, o3) = b.multipliedReportingOverflow(by: 16)
        return o1 || o2 || o3 ? Int.max : c
    }

    /// Half of the physical memory [bytes], the default limit.
    public static var defaultMemoryLimit: Int { Int(ProcessInfo.processInfo.physicalMemory / 2) }

    /// The detector applied to one exit wave (abTEM's `WavesDetector.detect`).
    public func detect(_ waves: Waves) throws -> Waves {
        guard let r = rows, let c = columns else { return waves }
        return try waves.resampled(rows: r, columns: c)
    }
}

/// The exit waves of a scan, one per position, positions row-major over (scan x, scan y): the values are
/// `[position][row][column]` in two flat arrays.
public struct ScanWaves: Sendable {
    public let scanShape: [Int]
    /// The grid, beam and normalisation of every wave.
    public let template: Waves
    public let real: [Double]
    public let imaginary: [Double]

    public var positionCount: Int { scanShape[0] * scanShape[1] }
    public var rows: Int { template.grid.nx }
    public var columns: Int { template.grid.ny }

    /// - Precondition: the arrays hold positions × rows × columns numbers.
    public init(scanShape: [Int], template: Waves, real: [Double], imaginary: [Double]) {
        var filled = scanShape.count == 2
        if filled {
            let (a, o1) = scanShape[0].multipliedReportingOverflow(by: scanShape[1])
            let (b, o2) = a.multipliedReportingOverflow(by: template.grid.count)
            filled = !o1 && !o2 && real.count == b && imaginary.count == real.count
        }
        precondition(filled, "ScanWaves: the values do not fill the shapes")
        self.scanShape = scanShape
        self.template = template
        self.real = real
        self.imaginary = imaginary
    }

    /// The exit wave at scan position (ix, iy).
    public subscript(ix: Int, iy: Int) -> Waves {
        let n = template.grid.count
        let base = (ix * scanShape[1] + iy) * n
        let field = ComplexField(nx: rows, ny: columns, real: Array(real[base..<(base + n)]),
                                 imaginary: Array(imaginary[base..<(base + n)]))
        return Waves(field: field, grid: template.grid, beam: template.beam, normalization: template.normalization,
                     tilt: template.tilt)
    }
}

extension STEM {
    /// Scans `probe` over `scan` through the potential and returns the complex exit wave at every position
    /// (abTEM's `probe.scan(potential, scan, WavesDetector(gpts))`). `.metal` runs float32 on the GPU and is held to `.cpu`
    /// at `MetalAccuracy.registered.wave`.
    /// - Parameter memoryLimit: refuse above this many bytes of waves (default half of the physical memory).
    /// - Throws: `WavesDetectorError.severalPotentials` for frozen phonons (scan each configuration on its own),
    ///   `WavesDetectorError.tooLarge`.
    public static func scan(probe: Probe, potentials: [ProjectedPotential], scan: GridScan, detector: WavesDetector,
                            backend: Backend = .cpu, options: ExecutionOptions = .standard,
                            memoryLimit: Int = WavesDetector.defaultMemoryLimit) throws -> ScanWaves {
        guard potentials.count == 1 else { throw WavesDetectorError.severalPotentials(potentials.count) }
        guard potentials[0].grid.matches(probe.grid) else { throw SimulatoError.gridMismatch }
        try checkMemory(probe: probe, scan: scan, detector: detector, limit: memoryLimit)
        let potential = potentials[0]
        return try collect(probe: probe, scan: scan, detector: detector, options: options) { incidents in
            try Multislice.runBatch(incidents, through: potential, backend: backend, options: options)
        }
    }

    /// The same through an unbuilt potential: `.cpu` builds it, `.metal` builds each slice on the GPU.
    public static func scan(probe: Probe, potentials: [Potential], scan: GridScan, detector: WavesDetector,
                            backend: Backend = .cpu, options: ExecutionOptions = .standard,
                            memoryLimit: Int = WavesDetector.defaultMemoryLimit) throws -> ScanWaves {
        guard potentials.count == 1 else { throw WavesDetectorError.severalPotentials(potentials.count) }
        guard potentials[0].grid.matches(probe.grid) else { throw SimulatoError.gridMismatch }
        try checkMemory(probe: probe, scan: scan, detector: detector, limit: memoryLimit)
        let potential = potentials[0]
        return try collect(probe: probe, scan: scan, detector: detector, options: options) { incidents in
            try Multislice.runBatch(incidents, through: potential, backend: backend, options: options)
        }
    }

    private static func checkMemory(probe: Probe, scan: GridScan, detector: WavesDetector, limit: Int) throws {
        let r = detector.rows ?? probe.grid.nx, c = detector.columns ?? probe.grid.ny
        let (positions, overflow) = scan.nx.multipliedReportingOverflow(by: scan.ny)
        guard !overflow else {
            throw WavesDetectorError.tooLarge(positions: .max, rows: r, columns: c, required: .max, limit: limit)
        }
        let need = WavesDetector.requiredBytes(positions: positions, rows: r, columns: c)
        guard need <= limit else {
            throw WavesDetectorError.tooLarge(positions: positions, rows: r, columns: c, required: need, limit: limit)
        }
    }

    private static func collect(probe: Probe, scan: GridScan, detector: WavesDetector, options: ExecutionOptions,
                                run: ([Waves]) throws -> [Waves]) throws -> ScanWaves {
        let positions = scan.positions
        let r = detector.rows ?? probe.grid.nx, c = detector.columns ?? probe.grid.ny
        let need = WavesDetector.requiredBytes(positions: positions.count, rows: r, columns: c)
        guard need != .max else {
            throw WavesDetectorError.tooLarge(positions: positions.count, rows: r, columns: c, required: .max, limit: .max)
        }
        var re = [Double](), im = [Double]()
        re.reserveCapacity(need / 16)
        im.reserveCapacity(need / 16)
        var template: Waves?
        // Chunks as the measurement scan: the incident and exit waves of one chunk only, besides the result.
        let chunk = max(1, options.scanChunkBytes / (probe.grid.count * 16 * 2))
        for first in stride(from: 0, to: positions.count, by: chunk) {
            let incidents = try positions[first..<min(first + chunk, positions.count)].map { try probe.build(at: $0) }
            for exit in try run(incidents) {
                let w = try detector.detect(exit)
                re.append(contentsOf: w.field.real)
                im.append(contentsOf: w.field.imaginary)
                if template == nil { template = w }
            }
        }
        guard let t = template else { throw SimulatoError.invalidScan("no scan position") }
        return ScanWaves(scanShape: [scan.nx, scan.ny], template: t, real: re, imaginary: im)
    }
}

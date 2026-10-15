//
//  STEM.swift
//  Role: the STEM driver (registration-slice-2.md sections 3.5, 3.6): every configuration, every scan position, the
//        probe there through that configuration's potential, each detector on the exit wave, the mean over
//        configurations (abTEM's `ensemble_mean`). Positions run as one batch per configuration (.cpu one by one,
//        .metal batched on the GPU).
//

import SimulatoCore
import SimulatoMultislice
import SimulatoProbe

/// One detector's measurement over a scan: values row-major over (scan x, scan y, the detector's shape). Named
/// `ScanMeasurement`, not abTEM's `Measurement`, so it does not collide with Foundation's `Measurement`.
public struct ScanMeasurement: Sendable {
    public let scanShape: [Int]
    public let valueShape: [Int]
    public let values: [Double]

    var valueCount: Int { valueShape.reduce(1, *) }

    /// The values at scan position (ix, iy).
    public subscript(ix: Int, iy: Int) -> [Double] {
        let base = (ix * scanShape[1] + iy) * valueCount
        return Array(values[base..<(base + valueCount)])
    }

    /// The image of a one-number detector (annular), sampled as the scan.
    public func image(samplingX: Double, samplingY: Double) -> RealField? {
        guard valueCount == 1 else { return nil }
        return RealField(nx: scanShape[0], ny: scanShape[1], values: values, samplingX: samplingX, samplingY: samplingY)
    }
}

public enum STEM {
    /// Scans `probe` over `scan` through every potential (one per frozen-phonon configuration; one for a static
    /// structure) and returns each detector's measurement averaged over the potentials. On `.metal` the probes, the
    /// transmissions and the diffraction intensities are computed on the GPU (registration-performance.md 1.6).
    public static func scan(probe: Probe, potentials: [ProjectedPotential], scan: GridScan, detectors: [any Detector],
                            backend: Backend = .cpu) throws -> [ScanMeasurement] {
        guard !potentials.isEmpty else { throw SimulatoError.invalidScan("no potential") }
        for p in potentials where !p.grid.matches(probe.grid) { throw SimulatoError.gridMismatch }
        switch backend {
        case .cpu:
            return try cpuScan(probe: probe, potentials: potentials, scan: scan, detectors: detectors)
        case .metal:
            return try metalScan(probe: probe, count: potentials.count, scan: scan, detectors: detectors) { k, s, n, r, d in
                try Multislice.metalScanIntensities(spectrum: s, positions: n, ramps: r, through: potentials[k],
                                                    beam: probe.beam, deliver: d)
            }
        }
    }

    /// The same through unbuilt potentials: on `.cpu` each is built (`build()`), on `.metal` each slice is built on the GPU.
    public static func scan(probe: Probe, potentials: [Potential], scan: GridScan, detectors: [any Detector],
                            backend: Backend = .cpu) throws -> [ScanMeasurement] {
        guard !potentials.isEmpty else { throw SimulatoError.invalidScan("no potential") }
        for p in potentials where !p.grid.matches(probe.grid) { throw SimulatoError.gridMismatch }
        switch backend {
        case .cpu:
            return try cpuScan(probe: probe, potentials: potentials.map { try $0.build() }, scan: scan, detectors: detectors)
        case .metal:
            return try metalScan(probe: probe, count: potentials.count, scan: scan, detectors: detectors) { k, s, n, r, d in
                try Multislice.metalScanIntensities(spectrum: s, positions: n, ramps: r, through: potentials[k],
                                                    beam: probe.beam, deliver: d)
            }
        }
    }

    static func cpuScan(probe: Probe, potentials: [ProjectedPotential], scan: GridScan,
                        detectors: [any Detector]) throws -> [ScanMeasurement] {
        let positions = scan.positions
        let shapes = try detectors.map { try $0.shape(for: probe.grid, beam: probe.beam) }
        let counts = shapes.map { $0.reduce(1, *) }
        var sums = counts.map { [Double](repeating: 0, count: $0 * positions.count) }
        // Positions in chunks, so the waves of one chunk only are held (slice-2 refuter: 13 MB per position at 480²).
        let chunk = max(1, (256 << 20) / (probe.grid.count * 16 * 2))
        for first in stride(from: 0, to: positions.count, by: chunk) {
            let range = first..<min(first + chunk, positions.count)
            let incidents = try range.map { try probe.build(at: positions[$0]) }
            for potential in potentials {
                let exits = try Multislice.runBatch(incidents, through: potential, backend: .cpu)
                for (offset, exit) in exits.enumerated() {
                    try accumulate(DiffractionIntensity(exit), at: first + offset, detectors, counts, &sums)
                }
            }
        }
        return average(sums, shapes, scan: scan, over: potentials.count)
    }

    typealias Ramps = (Int) throws -> (x: ComplexField, y: ComplexField)

    /// Potentials outer, positions inner: each position still sums its configurations in potential order.
    static func metalScan(probe: Probe, count: Int, scan: GridScan, detectors: [any Detector],
                          run: (Int, ComplexField, Int, Ramps, (Int, [Double]) throws -> Void) throws -> Void)
        throws -> [ScanMeasurement] {
        let positions = scan.positions
        let shapes = try detectors.map { try $0.shape(for: probe.grid, beam: probe.beam) }
        let counts = shapes.map { $0.reduce(1, *) }
        var sums = counts.map { [Double](repeating: 0, count: $0 * positions.count) }
        let spectrum = probe.centredSpectrum()
        for k in 0..<count {
            try run(k, spectrum, positions.count, { try probe.shiftRamps(at: positions[$0]) }) { q, values in
                let intensity = DiffractionIntensity(grid: probe.grid, beam: probe.beam, values: values)
                try accumulate(intensity, at: q, detectors, counts, &sums)
            }
        }
        return average(sums, shapes, scan: scan, over: count)
    }

    static func accumulate(_ intensity: DiffractionIntensity, at q: Int, _ detectors: [any Detector], _ counts: [Int],
                           _ sums: inout [[Double]]) throws {
        for (d, detector) in detectors.enumerated() {
            let v = try detector.detect(intensity)
            precondition(v.count == counts[d], "a detector returned the wrong number of values")
            for k in 0..<v.count { sums[d][q * counts[d] + k] += v[k] }
        }
    }

    static func average(_ sums: [[Double]], _ shapes: [[Int]], scan: GridScan, over n: Int) -> [ScanMeasurement] {
        zip(sums, shapes).map { values, shape in
            ScanMeasurement(scanShape: [scan.nx, scan.ny], valueShape: shape, values: values.map { $0 / Double(n) })
        }
    }
}

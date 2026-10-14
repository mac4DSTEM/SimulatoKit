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

/// One detector's measurement over a scan: values row-major over (scan x, scan y, the detector's shape).
public struct Measurement: Sendable {
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
    /// structure) and returns each detector's measurement averaged over the potentials.
    public static func scan(probe: Probe, potentials: [ProjectedPotential], scan: GridScan, detectors: [any Detector],
                            backend: Backend = .cpu) throws -> [Measurement] {
        guard !potentials.isEmpty else { throw SimulatoError.invalidScan("no potential") }
        for p in potentials where !p.grid.matches(probe.grid) { throw SimulatoError.gridMismatch }
        let positions = scan.positions
        let shapes = try detectors.map { try $0.shape(for: probe.grid, beam: probe.beam) }
        let counts = shapes.map { $0.reduce(1, *) }
        var sums = counts.map { [Double](repeating: 0, count: $0 * positions.count) }
        let incidents = try positions.map { try probe.build(at: $0) }
        for potential in potentials {
            let exits = try Multislice.runBatch(incidents, through: potential, backend: backend)
            for (q, exit) in exits.enumerated() {
                let intensity = try DiffractionIntensity(exit)
                for (d, detector) in detectors.enumerated() {
                    let v = try detector.detect(intensity)
                    precondition(v.count == counts[d], "a detector returned the wrong number of values")
                    for k in 0..<v.count { sums[d][q * counts[d] + k] += v[k] }
                }
            }
        }
        let n = Double(potentials.count)
        return zip(sums, shapes).map { values, shape in
            Measurement(scanShape: [scan.nx, scan.ny], valueShape: shape, values: values.map { $0 / n })
        }
    }
}

//
//  DetectorTruthTests.swift
//  Role: registration-slice-2.md section 5, items 3–8: Parseval, BF + DF, a flat vacuum scan, HAADF rising with Z, DPC
//        symmetry, the pixelated pattern's sum; and the detectors' refusals.
//

import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoMultislice
import SimulatoProbe

@Suite("Detector truth")
struct DetectorTruthTests {
    static func alExit(at position: SIMD2<Double>) throws -> Waves {
        let potentials = try Pins2.potentials("S1", pi: .exact)
        let probe = try Pins2.probe("S1", grid: potentials[0].grid)
        return try Multislice.run(probe.build(at: position), through: potentials[0])[0]
    }

    /// 3, 4 and 8: a detector over all angles collects N Σ|ψ|²; BF + DF = total; segments sum to their annulus.
    @Test(arguments: [SIMD2(0.0, 0.0), SIMD2(1.3, 2.9)])
    func detectorsConserveTheTotal(position: SIMD2<Double>) throws {
        let exit = try Self.alExit(at: position)
        let p = try DiffractionIntensity(exit)
        let n = Double(exit.grid.count)
        let a = p.angularSampling
        let maximum = min(Double(exit.grid.nx / 2) * a.x, Double(exit.grid.ny / 2) * a.y)
        let all = try AnnularDetector(inner: 0, outer: maximum).detect(p)[0]
        #expect(abs(all - n * exit.totalIntensity) <= 1e-12 * all, "Parseval")
        let full = try PixelatedDetector(maxAngle: .full).detect(p).reduce(0, +)
        #expect(abs(full - n * exit.totalIntensity) <= 1e-12 * full, "the pixelated pattern sums to the total")
        let bf = try AnnularDetector(inner: 0, outer: 20).detect(p)[0]
        let df = try AnnularDetector(inner: 20, outer: 120).detect(p)[0]
        let both = try AnnularDetector(inner: 0, outer: 120).detect(p)[0]
        #expect(abs(bf + df - both) <= 1e-12 * both, "BF + DF = total")
        let segments = try SegmentedDetector(radialBins: 3, azimuthalBins: 6, inner: 5, outer: 40, rotation: 0.7).detect(p)
        let ring = try AnnularDetector(inner: 5, outer: 40).detect(p)[0]
        #expect(abs(segments.reduce(0, +) - ring) <= 1e-12 * ring, "segments sum to their annulus")
    }

    /// 5 and 7: in vacuum every image is flat, and the four quadrants of a centred probe are equal.
    @Test func vacuumScanIsFlatAndDPCIsSymmetric() throws {
        let vac = try AtomicStructure(box: SIMD3(8.1, 8.1, 6), atoms: [])
        let p = try ProjectedPotential(structure: vac, nx: 128, ny: 128, sliceThickness: 2)
        let probe = try Probe(beam: ElectronBeam(energy: 200e3), grid: p.grid, aperture: Aperture(semiangleCutoff: 20))
        let scan = try GridScan(start: SIMD2(0, 0), end: SIMD2(4.05, 4.05), nx: 6, ny: 6)
        let detectors: [any Detector] = [try AnnularDetector(inner: 0, outer: 20), try AnnularDetector(inner: 5, outer: 15),
                                         try SegmentedDetector(radialBins: 1, azimuthalBins: 4, inner: 2, outer: 20),
                                         try SegmentedDetector(radialBins: 1, azimuthalBins: 4, inner: 0, outer: 20)]
        let m = try STEM.scan(probe: probe, potentials: [p], scan: scan, detectors: detectors)
        for (k, measurement) in m.enumerated() {
            let count = measurement.valueShape.reduce(1, *)
            for v in 0..<count {
                let column = stride(from: v, to: measurement.values.count, by: count).map { measurement.values[$0] }
                let spread = column.max()! - column.min()!
                #expect(spread <= 1e-12 * max(column.map(abs).max()!, 1e-300), "detector \(k) value \(v) is flat")
            }
        }
        // Amendment 2: the quadrants are equal once the origin pixel is excluded (inner > 0); with inner = 0 the zero
        // frequency (φ = atan2(0, 0) = 0) belongs to the first segment alone, in abTEM's rule as in the port.
        let quadrants = m[2][0, 0]
        #expect(quadrants.max()! - quadrants.min()! <= 1e-12 * quadrants.max()!, "DPC quadrants equal (inner > 0)")
        let withOrigin = m[3][0, 0]
        #expect(withOrigin[0] > withOrigin[1] && abs(withOrigin[1] - withOrigin[2]) <= 1e-12 * withOrigin[1]
                    && abs(withOrigin[2] - withOrigin[3]) <= 1e-12 * withOrigin[2],
                "inner = 0: the first segment holds the origin pixel, the other three are equal")
    }

    /// 6: HAADF on a single atom rises with Z.
    @Test func haadfRisesWithZ() throws {
        let values = try [6, 14, 29, 79].map { try STEMPinsTests.singleAtomHAADF($0, pi: .exact) }
        #expect(zip(values, values.dropFirst()).allSatisfy { $0 < $1 }, "\(values)")
    }

    @Test func refusals() throws {
        #expect(throws: (any Error).self) { try AnnularDetector(inner: -1) }
        #expect(throws: (any Error).self) { try AnnularDetector(inner: 50, outer: 20) }
        #expect(throws: (any Error).self) { try SegmentedDetector(radialBins: 0, azimuthalBins: 4, inner: 0, outer: 20) }
        #expect(throws: (any Error).self) { try SegmentedDetector(radialBins: 1, azimuthalBins: 4, inner: 20, outer: 20) }
        #expect(throws: (any Error).self) { try GridScan(start: SIMD2(0, 0), end: SIMD2(1, 1), nx: 0, ny: 3) }
        let exit = try Self.alExit(at: SIMD2(0, 0))
        #expect(throws: SimulatoError.self) {
            try AnnularDetector(inner: 0, outer: 1000).detect(DiffractionIntensity(exit))
        }
        #expect(!SimulatoDetector.plannedScope.isEmpty)
    }
}

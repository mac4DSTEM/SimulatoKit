//
//  DetectorRefusalsAndVacuumTests.swift
//  Role: the detectors' truths that need no fixture: a flat vacuum scan with DPC symmetry, past review findings each
//        pinned by a number abTEM gave, and the constructors' refusals. The truths on the pinned Al exit wave are in
//        SimulatoKitParityTests/DetectorTruthTests.
//

import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoMultislice
import SimulatoProbe

@Suite("Detectors: vacuum, review findings, refusals")
struct DetectorRefusalsAndVacuumTests {
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
        // The quadrants are equal once the origin pixel is excluded (inner > 0); with inner = 0 the zero
        // frequency (φ = atan2(0, 0) = 0) belongs to the first segment alone, in abTEM's rule as in the port.
        let quadrants = m[2][0, 0]
        #expect(quadrants.max()! - quadrants.min()! <= 1e-12 * quadrants.max()!, "DPC quadrants equal (inner > 0)")
        let withOrigin = m[3][0, 0]
        #expect(withOrigin[0] > withOrigin[1] && abs(withOrigin[1] - withOrigin[2]) <= 1e-12 * withOrigin[1]
                    && abs(withOrigin[2] - withOrigin[3]) <= 1e-12 * withOrigin[2],
                "inner = 0: the first segment holds the origin pixel, the other three are equal")
    }

    /// Past review findings, each pinned by a number abTEM gave.
    @Test func reviewFindingsStayFixed() throws {
        let beam = try ElectronBeam(energy: 200e3)
        let grid = try Grid(nx: 96, ny: 96, extentX: 7.81, extentY: 7.81)
        // B1: abTEM's angular sampling on this grid, to the last bit.
        let vac = try AtomicStructure(box: SIMD3(7.81, 7.81, 2), atoms: [])
        let p0 = try ProjectedPotential(structure: vac, nx: 96, ny: 96, sliceThickness: 2)
        let probe = try Probe(beam: beam, grid: grid, aperture: Aperture(semiangleCutoff: 20))
        let exit = try Multislice.run(probe.build(at: SIMD2(1.1, 2.2)), through: p0)[0]
        let p = try DiffractionIntensity(exit)
        #expect(p.angularSampling.x == 3.211183139222596, "abTEM's Δα on 96 px of 7.81 Å at 200 keV")
        // A: an outer one ulp above a pixel's angle, at a bin count where n_r α / outer rounds up to n_r: the pixel
        // is dropped (abTEM's rule) instead of indexing past the end.
        let angles = Grid.fftfreq(96, 1 / p.angularSampling.x / 96)
        var found: (k: Int, bins: Int)? = nil
        search: for k in 1..<40 {
            for bins in 2...12 where (Double(bins) * angles[k] / angles[k].nextUp).rounded(.down) == Double(bins) {
                found = (k, bins); break search
            }
        }
        let hit = try #require(found, "a pixel whose radial bin rounds up to n_r exists on this grid")
        let segments = try SegmentedDetector(radialBins: hit.bins, azimuthalBins: 4, inner: 0,
                                             outer: angles[hit.k].nextUp).detect(p)
        #expect(segments.count == hit.bins * 4)
        // B2: abTEM's isclose has rtol 1e-5 on outer: max + 5e-4 accepted, max + 3e-3 refused.
        let maximum = Double(48) * p.angularSampling.x
        #expect(throws: Never.self) { try AnnularDetector(inner: 0, outer: maximum + 5e-4).detect(p) }
        #expect(throws: SimulatoError.self) { try AnnularDetector(inner: 0, outer: maximum + 3e-3).detect(p) }
        // inner beyond the default outer is refused (abTEM raises).
        #expect(throws: SimulatoError.self) { try AnnularDetector(inner: 150).detect(p) }
    }

    @Test func refusals() throws {
        #expect(throws: (any Error).self) { try AnnularDetector(inner: -1) }
        #expect(throws: (any Error).self) { try AnnularDetector(inner: 50, outer: 20) }
        #expect(throws: (any Error).self) { try SegmentedDetector(radialBins: 0, azimuthalBins: 4, inner: 0, outer: 20) }
        #expect(throws: (any Error).self) { try SegmentedDetector(radialBins: 1, azimuthalBins: 4, inner: 20, outer: 20) }
        #expect(throws: (any Error).self) { try GridScan(start: SIMD2(0, 0), end: SIMD2(1, 1), nx: 0, ny: 3) }
        #expect(throws: (any Error).self) { try GridScan(start: SIMD2(1, 1), end: SIMD2(0, 1), nx: 3, ny: 3) }
        #expect(throws: Never.self) { try GridScan(start: SIMD2(1, 0), end: SIMD2(0, 1), nx: 3, ny: 3) }
        #expect(!SimulatoDetector.plannedScope.isEmpty)
    }
}

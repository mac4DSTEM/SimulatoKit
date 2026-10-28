//
//  DetectorTruthTests.swift
//  Role: the detectors' truths independent of abTEM: Parseval, BF + DF, a flat vacuum scan, HAADF rising with Z, DPC
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

    /// Section 3.4: the annulus is half-open. With its outer edge exactly on a pixel's angle (the 5th pixel along x),
    /// that pixel belongs to the next annulus, not to both.
    @Test func annulusIsHalfOpenOnAPixel() throws {
        let exit = try Self.alExit(at: SIMD2(1.3, 2.9))
        let p = try DiffractionIntensity(exit)
        let a = p.angularSampling
        let edge = Grid.fftfreq(exit.grid.nx, 1 / a.x / Double(exit.grid.nx))[5]  // the pixel (5, 0)'s angle, exact
        let inner = try AnnularDetector(inner: 0, outer: edge).detect(p)[0]
        let outer = try AnnularDetector(inner: edge, outer: 2 * edge).detect(p)[0]
        let both = try AnnularDetector(inner: 0, outer: 2 * edge).detect(p)[0]
        #expect(abs(inner + outer - both) <= 1e-13 * both, "a pixel on the edge is counted once")
        let pixel = p.values[5 * exit.grid.ny]
        let closedEdgeExcess = try AnnularDetector(inner: 0, outer: edge.nextUp).detect(p)[0] - inner
        #expect(closedEdgeExcess >= pixel, "moving outer past the pixel adds it (it was not in [0, edge))")
    }

    /// Section 3.4: outer = nil is ⌊min(cutoff angles)⌋ (S1's grid: ⌊43 Δα⌋ = ⌊133.14⌋ = 133 mrad).
    @Test func defaultOuterIsTheFlooredCutoffAngle() throws {
        let exit = try Self.alExit(at: SIMD2(0, 0))
        let p = try DiffractionIntensity(exit)
        let c = p.cutoffAngles
        #expect(abs(c.x - 43 * p.angularSampling.x) <= 1e-12 * c.x)
        #expect(try AnnularDetector(inner: 0).outerAngle(p) == 133)
        #expect(try AnnularDetector(inner: 0).outerAngle(p) == (min(c.x, c.y)).rounded(.down))
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
        // The quadrants are equal once the origin pixel is excluded (inner > 0); with inner = 0 the zero
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
        let exit = try Self.alExit(at: SIMD2(0, 0))
        #expect(throws: SimulatoError.self) {
            try AnnularDetector(inner: 0, outer: 1000).detect(DiffractionIntensity(exit))
        }
        #expect(!SimulatoDetector.plannedScope.isEmpty)
    }
}

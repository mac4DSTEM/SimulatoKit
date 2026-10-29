//
//  DetectorTruthTests.swift
//  Role: the detectors' truths independent of abTEM on the pinned Al exit wave: Parseval, BF + DF, the half-open
//        annulus, the default outer angle, HAADF rising with Z, the pixelated pattern's sum, an outer angle past the
//        pattern refused. The fixture-free truths (a flat vacuum scan, DPC symmetry, the refusals) are in SimulatoDetectorTests.
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

    /// 6: HAADF on a single atom rises with Z.
    @Test func haadfRisesWithZ() throws {
        let values = try [6, 14, 29, 79].map { try STEMPinsTests.singleAtomHAADF($0, pi: .exact) }
        #expect(zip(values, values.dropFirst()).allSatisfy { $0 < $1 }, "\(values)")
    }

    @Test func anOuterAngleBeyondThePatternIsRefused() throws {
        let exit = try Self.alExit(at: SIMD2(0, 0))
        #expect(throws: SimulatoError.self) {
            try AnnularDetector(inner: 0, outer: 1000).detect(DiffractionIntensity(exit))
        }
    }
}

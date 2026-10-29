//
//  LineScanTruthTests.swift
//  Role: LineScan and interpolate_line by truths that need no reference: a vacuum line scan is flat, a plane is
//        interpolated exactly, and the refusals. The comparisons with abTEM are in SimulatoKitParityTests
//        (LineScanTests, LineInterpolationTests).
//

import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoMultislice
import SimulatoProbe

@Suite("Line scans and line interpolation: truths and refusals")
struct LineScanTruthTests {
    @Test func aVacuumLineScanIsFlat() throws {
        let vac = try AtomicStructure(box: SIMD3(8, 8, 4), atoms: [])
        let p = try ProjectedPotential(structure: vac, nx: 64, ny: 64, sliceThickness: 2)
        let probe = try Probe(beam: ElectronBeam(energy: 200e3), grid: p.grid, aperture: Aperture(semiangleCutoff: 20))
        let scan = try LineScan(start: SIMD2(1, 1), end: SIMD2(5, 3), gpts: 7)
        let m = try STEM.scan(probe: probe, potentials: [p], scan: scan, detectors: [AnnularDetector(inner: 0, outer: 20)])[0]
        #expect(m.scanShape == [7] && m.valueShape == [])
        let spread = m.values.max()! - m.values.min()!
        #expect(spread <= 1e-12 * m.values.max()!, "a vacuum line scan is flat")
        let profile = try m.lineProfiles(for: scan)
        #expect(profile.count == 7 && profile.end == SIMD2(5, 3))
    }

    @Test func aLinearImageIsInterpolatedExactly() throws {
        // a plane 2 i + 3 j: linear interpolation is exact between pixels (inside, clear of the periodic seam)
        let image = RealField(nx: 20, ny: 16, values: (0..<20).flatMap { i in (0..<16).map { j in 2.0 * Double(i) + 3.0 * Double(j) } },
                              samplingX: 0.5, samplingY: 0.5)
        let r = try image.interpolateLine(start: SIMD2(2, 1), end: SIMD2(5, 3.5), gpts: 7, order: 1, endpoint: true)
        for (k, v) in r.values.enumerated() {
            let t = Double(k) / 6
            let x = (2 + 3 * t) / 0.5, y = (1 + 2.5 * t) / 0.5
            #expect(abs(v - (2 * x + 3 * y)) < 1e-12)
        }
        #expect(r.start == SIMD2(2, 1) && r.end == SIMD2(5, 3.5) && r.lineWidth == nil)
    }

    @Test func lineScanRefusals() throws {
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(1, 1)) }
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(1, 1), gpts: 0) }
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(1, 1), sampling: -1) }
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(.nan, 1), gpts: 3) }
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(1, 1), gpts: 3, fractional: true) }
        #expect(throws: SimulatoError.self) { try LineScan(start: SIMD2(1, 1), end: SIMD2(1, 1), sampling: 0.5) }
        let point = try LineScan(start: SIMD2(1, 1), end: SIMD2(1, 1), gpts: 2)
        #expect(throws: SimulatoError.self) { try point.addingMargin(0.1) }
        let line = try LineScan(end: SIMD2(1, 0), gpts: 5)
        #expect(throws: SimulatoError.self) { try line.addingMargin(-0.5) }  // both ends meet: no positions
    }

    @Test func interpolationRefusals() throws {
        let image = RealField(nx: 4, ny: 4, values: [Double](repeating: 1, count: 16), samplingX: 1, samplingY: 1)
        #expect(throws: SimulatoError.self) { try image.interpolateLine(end: SIMD2(1, 1), gpts: 3, order: 7) }
        #expect(throws: SimulatoError.self) { try image.interpolateLine(end: SIMD2(1, 1), gpts: 3, width: -1) }
        #expect(throws: SimulatoError.self) { try image.interpolateLine(start: SIMD2(1, 1), end: SIMD2(1, 1), gpts: 3, width: 1) }
    }
}

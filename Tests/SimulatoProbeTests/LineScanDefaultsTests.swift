//
//  LineScanDefaultsTests.swift
//  Role: a line scan matched to a probe (abTEM's LineScan(end=None) and match_probe), against numbers abTEM 1.0.10 printed
//        (tools/pins/linescan-defaults/pins.py) for a 20 mrad probe at 200 keV and a 30 mrad probe at 80 keV on a
//        12.3 x 9.7 A grid; and its refusals.
//

import Foundation
import Testing
import SimulatoCore
import SimulatoProbe

@Suite("Line scan defaults matched to a probe")
struct LineScanDefaultsTests {
    static func probe(_ cutoff: Double, _ energy: Double) throws -> Probe {
        try Probe(beam: ElectronBeam(energy: energy), grid: Grid(nx: 64, ny: 48, extentX: 12.3, extentY: 9.7),
                  aperture: Aperture(semiangleCutoff: cutoff))
    }

    /// End points and counts exact; the sampling at the assert_allclose default rtol 1e-7 (abTEM states no tolerance).
    static func check(_ s: LineScan, start: SIMD2<Double>, end: SIMD2<Double>, gpts: Int, sampling: Double) -> Bool {
        s.start == start && s.end == end && s.gpts == gpts && abs(s.sampling - sampling) <= 1e-7 * sampling
    }

    @Test func aMissingEndAndSamplingFollowTheProbeAsAbTEM() throws {
        let p20 = try Self.probe(20, 200e3), p30 = try Self.probe(30, 80e3)
        let end = SIMD2(0.0, 9.7)
        #expect(Self.check(try LineScan(matching: p20), start: .zero, end: end, gpts: 32, sampling: 0.3129032258064516))
        #expect(Self.check(try LineScan(matching: p30), start: .zero, end: end, gpts: 29, sampling: 0.3464285714285714))
        #expect(Self.check(try LineScan(start: SIMD2(1, 2), matching: p20), start: SIMD2(1, 2), end: end, gpts: 26,
                           sampling: 0.3105865418848666))
        #expect(Self.check(try LineScan(start: SIMD2(1, 2), matching: p30), start: SIMD2(1, 2), end: end, gpts: 23,
                           sampling: 0.35293925214189387))
        #expect(Self.check(try LineScan(endpoint: false, matching: p20), start: .zero, end: end, gpts: 32, sampling: 0.303125))
        #expect(Self.check(try LineScan(endpoint: false, matching: p30), start: .zero, end: end, gpts: 29,
                           sampling: 0.33448275862068966))
    }

    @Test func aGivenEndKeepsItsPointAndAGivenSamplingOrCountWins() throws {
        let p20 = try Self.probe(20, 200e3), p30 = try Self.probe(30, 80e3)
        #expect(Self.check(try LineScan(start: SIMD2(1, 1), end: SIMD2(8, 5), matching: p20), start: SIMD2(1, 1),
                           end: SIMD2(8, 5), gpts: 26, sampling: 0.32249030993194194))
        #expect(Self.check(try LineScan(start: SIMD2(1, 1), end: SIMD2(8, 5), matching: p30), start: SIMD2(1, 1),
                           end: SIMD2(8, 5), gpts: 24, sampling: 0.3505329455781978))
        for p in [p20, p30] {
            #expect(Self.check(try LineScan(end: SIMD2(8, 5), gpts: 17, matching: p), start: .zero, end: SIMD2(8, 5),
                               gpts: 17, sampling: 0.5896238207535377))
            #expect(Self.check(try LineScan(sampling: 0.3, matching: p), start: .zero, end: SIMD2(0, 9.7), gpts: 33,
                               sampling: 0.303125))
        }
        // The matched scan is the plain scan with the same end and sampling.
        let s = try LineScan(matching: p20)
        #expect(s == (try LineScan(end: SIMD2(0, 9.7), sampling: 0.99 * LineScan.nyquistSampling(
            semiangleCutoff: 20, wavelength: p20.beam.wavelength))))
    }

    @Test func refusals() throws {
        let p20 = try Self.probe(20, 200e3)
        // abTEM silently replaces this count by the probe-matched 32 positions.
        #expect(throws: SimulatoError.self) { try LineScan(gpts: 10, matching: p20) }
        let open = try Probe(beam: ElectronBeam(energy: 200e3), grid: Grid(nx: 64, ny: 48, extentX: 12.3, extentY: 9.7),
                             aperture: Aperture(semiangleCutoff: .infinity))
        #expect(throws: SimulatoError.self) { try LineScan(matching: open) }
        #expect(throws: SimulatoError.self) { try LineScan(sampling: 0, matching: p20) }
        // A given count or sampling needs no cutoff.
        #expect(try LineScan(end: SIMD2(8, 5), gpts: 4, matching: open).gpts == 4)
        #expect(try LineScan(sampling: 0.5, matching: open).gpts == 20)
    }
}

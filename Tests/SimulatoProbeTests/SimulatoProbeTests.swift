//
//  SimulatoProbeTests.swift
//  Role: SimulatoProbe's refusals and conventions (the probe's numbers are pinned in SimulatoMultisliceTests).
//

import Testing
import SimulatoCore
import SimulatoProbe

@Suite("SimulatoProbe")
struct SimulatoProbeTests {
    @Test func refusalsAndConventions() throws {
        #expect(throws: (any Error).self) { try Aperture(semiangleCutoff: 0) }
        #expect(throws: (any Error).self) { try Aperture(semiangleCutoff: .nan) }
        var a = Aberrations()
        a.defocus = 100
        #expect(a.C10 == -100, "defocus = −C10 (abTEM)")
        a.C30 = .infinity
        let beam = try ElectronBeam(energy: 200e3)
        let grid = try Grid(nx: 64, ny: 64, extentX: 10, extentY: 10)
        #expect(throws: (any Error).self) { try Probe(beam: beam, grid: grid, aperture: Aperture(semiangleCutoff: 20), aberrations: a) }
        let hard = try Aperture(semiangleCutoff: 20, soft: false).values(on: grid, wavelength: beam.wavelength)
        #expect(hard.allSatisfy { $0 == 0 || $0 == 1 })
    }

    @Test func lineScanConventions() throws {
        // endpoint puts the last position at `end`; without it the sampling is length / gpts (NumPy linspace).
        let a = try LineScan(end: SIMD2(3, 4), gpts: 6)
        #expect(a.extent == 5 && a.sampling == 1 && a.positions.first == SIMD2(0, 0) && a.positions.last == SIMD2(3, 4))
        let b = try LineScan(end: SIMD2(3, 4), gpts: 5, endpoint: false)
        #expect(b.sampling == 1 && b.positions.last == SIMD2(2.4, 3.2))
        // a sampling wins over gpts and gives ceil(length / sampling) positions.
        let c = try LineScan(end: SIMD2(3, 4), gpts: 2, sampling: 0.7)
        #expect(c.gpts == 8 && abs(c.sampling - 5.0 / 7.0) < 1e-15)
        // the margin moves both ends along the line and recomputes the count from the sampling.
        let d = try a.addingMargin(0.55)
        #expect(abs(d.start.x + 0.33) < 1e-15 && abs(d.start.y + 0.44) < 1e-15 && abs(d.end.x - 3.33) < 1e-15 && abs(d.end.y - 4.44) < 1e-15)
        // abTEM's setters in turn: the start move gives ceil(5.55 / 1) = 6 positions at 5.55 / 5 = 1.11, the end move ceil(6.1 / 1.11) = 6.
        #expect(d.gpts == 6 && abs(d.sampling - 6.1 / 5.0) < 1e-15)
        let e = try LineScan.atPosition(center: SIMD2(1, 1), extent: 2, angle: 90, gpts: 3)
        #expect(abs(e.start.x - 1) < 1e-15 && abs(e.start.y - 0) < 1e-15 && abs(e.end.y - 2) < 1e-15)
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(1, 0)) }
        // a grid scan is a scan pattern of two axes.
        #expect(try GridScan(start: .zero, end: SIMD2(1, 1), nx: 3, ny: 2).shape == [3, 2])
    }
}

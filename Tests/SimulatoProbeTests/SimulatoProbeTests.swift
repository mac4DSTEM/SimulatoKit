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
}

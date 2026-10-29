//
//  WavesDetectorRefusalTests.swift
//  Role: WavesDetector's typed refusals (too much memory, several or no potentials, a bad resample size) on the
//        aluminium cell of the waves pins, built here without the fixture. The pinned exit waves are in
//        SimulatoKitParityTests/WavesDetectorTests.
//

import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoKitTestSupport
import SimulatoMultislice
import SimulatoProbe

@Suite("WavesDetector refusals")
struct WavesDetectorRefusalTests {
    /// Two Al cells along z on 48 x 48, a/2 slices, 200 keV, 20 mrad: the waves pins' setup.
    static func setup() throws -> (probe: Probe, potential: Potential) {
        let potential = try Potential(structure: aluminiumStructure(cells: (1, 1, 2)), nx: 48, ny: 48, sliceThickness: 2.025,
                                      parameterization: LobatoParameterization(piConstant: .abTEMFloat32))
        let probe = try Probe(beam: ElectronBeam(energy: 200e3), grid: potential.grid, aperture: Aperture(semiangleCutoff: 20))
        return (probe, potential)
    }

    /// The 3 x 2 scan over one cell.
    static func scan() throws -> GridScan { try GridScan(start: SIMD2(0, 0), end: SIMD2(4.05, 4.05), nx: 3, ny: 2) }

    @Test func memoryRefusalIsTyped() throws {
        let (probe, potential) = try Self.setup()
        let scan = try Self.scan()
        // 6 positions x 48 x 48 x 16 bytes = 221 184 bytes.
        let need = WavesDetector.requiredBytes(positions: 6, rows: 48, columns: 48)
        #expect(need == 221_184)
        #expect(throws: WavesDetectorError.tooLarge(positions: 6, rows: 48, columns: 48, required: need, limit: need - 1)) {
            try STEM.scan(probe: probe, potentials: [potential], scan: scan, detector: WavesDetector(), memoryLimit: need - 1)
        }
        // The resampled size counts, not the simulation grid.
        let small = try WavesDetector(rows: 10, columns: 11)
        let ok = try STEM.scan(probe: probe, potentials: [potential], scan: scan, detector: small, memoryLimit: 6 * 10 * 11 * 16)
        #expect(ok.rows == 10 && ok.columns == 11)
        #expect(throws: WavesDetectorError.self) {
            try STEM.scan(probe: probe, potentials: [potential.build()], scan: scan, detector: small, memoryLimit: 6 * 10 * 11 * 16 - 1)
        }
        // Default: half of the physical memory; a huge scan is refused before any work.
        let huge = try GridScan(start: SIMD2(0, 0), end: SIMD2(4, 4), nx: 100_000, ny: 100_000)
        #expect(throws: WavesDetectorError.self) {
            try STEM.scan(probe: probe, potentials: [potential], scan: huge, detector: WavesDetector())
        }
        #expect(WavesDetector.requiredBytes(positions: .max, rows: 2, columns: 2) == .max)
        #expect("\(WavesDetectorError.tooLarge(positions: 1, rows: 2, columns: 3, required: 96, limit: 5))".contains("96 bytes"))
    }

    @Test func frozenPhononsAreRefused() throws {
        let (probe, potential) = try Self.setup()
        #expect(throws: WavesDetectorError.severalPotentials(2)) {
            try STEM.scan(probe: probe, potentials: [potential, potential], scan: Self.scan(), detector: WavesDetector())
        }
        #expect(throws: WavesDetectorError.severalPotentials(0)) {
            try STEM.scan(probe: probe, potentials: [Potential](), scan: Self.scan(), detector: WavesDetector())
        }
        #expect(throws: (any Error).self) { try WavesDetector(rows: 1, columns: 8) }
    }
}

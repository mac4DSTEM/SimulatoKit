//
//  HardeningTests.swift
//  Role: inputs that used to trap late in SimulatoDetector (a line scan's iCOM, a custom scan pattern of the wrong shape,
//        overflowing scans, empty profiles, huge widths) now throw their typed errors before any work.
//

import Foundation
import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoMultislice
import SimulatoProbe

struct BadScan: ScanPattern {
    var positions: [SIMD2<Double>]
    var shape: [Int]
}

@Suite("Detector hardening")
struct DetectorHardeningTests {
    func setup() throws -> (Probe, ProjectedPotential) {
        let grid = try Grid(nx: 16, ny: 16, extentX: 4, extentY: 4)
        let beam = try ElectronBeam(energy: 200e3)
        let probe = try Probe(beam: beam, grid: grid, aperture: Aperture(semiangleCutoff: 20))
        let potential = ProjectedPotential(grid: grid, slices: [PotentialSlice(thickness: 1, values: [Double](repeating: 0, count: 256))])
        return (probe, potential)
    }

    @Test func centreOfMassOfALineScanIsRefused() throws {
        let m = ScanMeasurement(scanShape: [4], valueShape: [2], values: [Double](repeating: 0, count: 8))
        #expect(throws: SimulatoError.self) { try m.gradientImages(samplingX: 1, samplingY: 1) }
        #expect(throws: SimulatoError.self) { try m.integratedCenterOfMass(samplingX: 1, samplingY: 1) }
    }

    @Test func malformedScanPatternIsRefusedBeforeTheSimulation() throws {
        let (probe, potential) = try setup()
        let p = [SIMD2<Double>(1, 1), SIMD2(2, 2), SIMD2(3, 3)]
        let detector = try AnnularDetector(inner: 0, outer: 10)
        for shape in [[2, 2], [], [1, 1, 3]] {
            let scan = BadScan(positions: p, shape: shape)
            #expect(throws: SimulatoError.self) {
                try STEM.scan(probe: probe, potentials: [potential], scan: scan, detectors: [detector])
            }
        }
    }

    @Test func overflowingWaveScanIsTyped() throws {
        let (probe, potential) = try setup()
        let scan = try GridScan(start: SIMD2(0, 0), end: SIMD2(4, 4), nx: 1 << 40, ny: 1 << 40)
        #expect(throws: WavesDetectorError.self) {
            try STEM.scan(probe: probe, potentials: [potential], scan: scan, detector: WavesDetector())
        }
    }

    @Test func emptyProfilesAreRefused() throws {
        let empty = RealSpaceLineProfiles(count: 0, values: [], sampling: 1)
        #expect(throws: SimulatoError.self) { try empty.width() }
        #expect(throws: SimulatoError.self) { try empty.interpolate(gpts: 3) }
        let p = RealSpaceLineProfiles(count: 4, values: [0, 1, 1, 0], sampling: 1)
        #expect(throws: SimulatoError.self) { try p.interpolate(gpts: 1 << 40) }
        #expect(try p.interpolate(gpts: 8).count == 8)
    }

    @Test func hugeWidthIsRefused() throws {
        let f = RealField(nx: 8, ny: 8, values: [Double](repeating: 1, count: 64), samplingX: 1, samplingY: 1)
        #expect(throws: SimulatoError.self) { try f.interpolateLine(end: SIMD2(4, 0), width: 1e300) }
        #expect(throws: SimulatoError.self) { try f.interpolateLine(end: SIMD2(4, 0), width: 1e12) }
        #expect(try f.interpolateLine(end: SIMD2(4, 0), width: 3).count > 0)
    }
}

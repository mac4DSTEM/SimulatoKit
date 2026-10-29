//
//  IntegratedCenterOfMassTests.swift
//  Role: the integrated centre of mass of a scan by its own truth: a planted gradient as a two-value scan is the
//        integral of its two images, and a vacuum scan has no centre of mass and a flat iCOM. The comparison with abTEM
//        is in SimulatoKitParityTests/IcomPinsTests; the line-scan refusal is in HardeningTests.
//

import Foundation
import Testing
import SimulatoCore
import SimulatoDetector
import SimulatoMultislice
import SimulatoProbe

@Suite("Integrated centre of mass: truths")
struct IntegratedCenterOfMassTests {
    @Test func aPlantedGradientIsIntegratedAsItsTwoImages() throws {
        let nx = 12, ny = 10, sx = 0.3, sy = 0.4
        let ex = Double(nx) * sx, ey = Double(ny) * sy
        var t = [Double](repeating: 0, count: nx * ny), gx = t, gy = t
        for (m, n, a, phase) in [(1, 0, 1.0, 0.3), (0, 1, 0.7, 1.1), (1, 2, 0.4, 5.2)] {
            let kx = Double(m) / ex, ky = Double(n) / ey
            for ix in 0..<nx {
                for iy in 0..<ny {
                    let arg = 2 * Double.pi * (kx * Double(ix) * sx + ky * Double(iy) * sy) + phase
                    t[ix * ny + iy] += a * cos(arg)
                    gx[ix * ny + iy] += -a * 2 * Double.pi * kx * sin(arg)
                    gy[ix * ny + iy] += -a * 2 * Double.pi * ky * sin(arg)
                }
            }
        }
        // a centre-of-mass scan holds (x, y) per position, rows by columns
        let scan = ScanMeasurement(scanShape: [nx, ny], valueShape: [2], values: zip(gx, gy).flatMap { [$0, $1] })
        let images = try scan.gradientImages(samplingX: sx, samplingY: sy)
        #expect(images.x.values == gx && images.y.values == gy && images.x.nx == nx && images.x.ny == ny)
        let icom = try scan.integratedCenterOfMass(samplingX: sx, samplingY: sy)
        let direct = try RealField.integratingGradient(x: images.x, y: images.y)
        #expect(icom.values == direct.values && icom.samplingX == sx && icom.samplingY == sy)
        // and it is the planted potential less its minimum
        let low = t.min()!, range = t.max()! - low
        let worst = zip(icom.values, t).map { abs($0 - ($1 - low)) }.max()!
        #expect(worst <= 1e-9 * range, "worst \(worst) of range \(range)")
    }

    @Test func aVacuumScanHasNoCentreOfMass() throws {
        let vac = try AtomicStructure(box: SIMD3(8, 8, 4), atoms: [])
        let p = try ProjectedPotential(structure: vac, nx: 64, ny: 64, sliceThickness: 2)
        let probe = try Probe(beam: ElectronBeam(energy: 200e3), grid: p.grid, aperture: Aperture(semiangleCutoff: 20))
        let scan = try GridScan(start: SIMD2(0, 0), end: SIMD2(4, 4), nx: 4, ny: 4)
        let com = try STEM.scan(probe: probe, potentials: [p], scan: scan,
                                detectors: [CenterOfMassDetector(maxAngle: .cutoff, units: .inverseAngstroms)])[0]
        #expect(com.scanShape == [4, 4] && com.valueShape == [2])
        // the aperture's radius in 1/Å sets the scale: 20 mrad / λ
        let scale = try 20e-3 / ElectronBeam(energy: 200e3).wavelength
        #expect(com.values.allSatisfy { abs($0) <= 1e-12 * scale }, "max |COM| \(com.values.map(abs).max()!)")
        let icom = try com.integratedCenterOfMass(samplingX: 1, samplingY: 1)
        #expect(icom.values.allSatisfy { abs($0) <= 1e-12 * scale })
    }
}

//
//  TruthTests.swift
//  Role: registration section 5, truth independent of abTEM: the far field of one atom against the analytic Lobato
//        projected potential (3a–3c, amendment 1), vacuum (4), the symmetry of Al [001] (5), the probe (6), χ (7) and
//        the thickness series (8). Truth 1 is in BeamAndScatteringTests, truth 2 in PotentialPinsTests.
//

import SimulatoKitTestSupport
import Foundation
import Testing
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe

@Suite("Truth independent of abTEM")
struct TruthTests {
    // MARK: 3, the far field of one Al atom

    @Test func singleAtomConvergesToTheAnalyticProjectedPotential() throws {
        let ff = Pins.json["truth_far_field_Al"] as! [String: Any]
        let grids = ff["grids"] as! [String: [String: Any]]
        var maxAbs: [Int: Double] = [:]
        for n in [512, 1024] {
            let g = grids[String(n)]!
            let offsets = g["pixel_offsets"] as! [Int], r = g["r"] as! [Double], analytic = g["V"] as! [Double]
            let atom = try AtomicStructure(box: SIMD3(25.6, 25.6, 2.0),
                                           atoms: [.init(atomicNumber: 13, position: SIMD3(12.8, 12.8, 1.0))])
            let v = try ProjectedPotential(structure: atom, nx: n, ny: n, sliceThickness: 2.0).slices[0].values
            var relNear = 0.0, absFar = 0.0, worst = 0.0
            for (k, off) in offsets.enumerated() {
                let value = v[(n / 2 + off) * n + n / 2]
                let d = abs(value - analytic[k])
                worst = max(worst, d)
                if r[k] <= 2.0 + 1e-9 { relNear = max(relNear, d / analytic[k]) } else { absFar = max(absFar, d) }
            }
            maxAbs[n] = worst
            if n == 512 {
                #expect(relNear <= 2e-2, "3a: relative error for r <= 2 Å")
                #expect(absFar <= 3e-2, "3b: absolute error for 2 < r <= 4 Å [V·Å]")
            }
        }
        #expect(maxAbs[1024]! / maxAbs[512]! <= 0.5, "3c: halving the pixel halves the error")
    }

    // MARK: 4, vacuum

    @Test func planeWaveThroughVacuumStaysOne() throws {
        let potential = try Pins.potential("C1", pi: .exact)
        let waves = try Multislice.run(.planeWave(beam: try Pins.beam("C1"), grid: potential.grid), through: potential,
                                       exitPlanes: .every(1))
        for w in waves {
            var d = 0.0
            for i in w.field.real.indices { d = max(d, abs(w.field.real[i] - 1), abs(w.field.imaginary[i])) }
            #expect(d <= 1e-14)
        }
    }

    @Test func probeThroughVacuumEqualsTheDefocusedProbe() throws {
        let potential = try Pins.potential("C1", pi: .exact)
        let beam = try Pins.beam("C1")
        var aberrations = Aberrations()
        aberrations.C10 = -50
        aberrations.C30 = 1e4
        aberrations.C12 = 20
        aberrations.phi12 = 0.3
        let position = SIMD2(8.1, 8.1)
        let aperture = try Aperture(semiangleCutoff: 20)
        let probe = try Probe(beam: beam, grid: potential.grid, aperture: aperture, aberrations: aberrations)
        let waves = try Multislice.run(probe.build(at: position), through: potential, exitPlanes: .every(1))
        var thickness = 0.0
        for (q, w) in waves.enumerated().dropFirst() {
            thickness += potential.slices[q - 1].thickness
            var shifted = aberrations
            shifted.C10 += thickness
            let expected = try Probe(beam: beam, grid: potential.grid, aperture: aperture, aberrations: shifted)
                .build(at: position)
            #expect(relMax(w.field, expected.field) <= 1e-12, "after \(thickness) Å of vacuum")
        }
    }

    /// Refuter finding: a grid equal within rounding (an extent 1 ulp off) is accepted; a different grid is refused.
    @Test func gridsEqualWithinRoundingAreAccepted() throws {
        let s = try AtomicStructure(box: SIMD3(9.6, 12.8, 4), atoms: [.init(atomicNumber: 6, position: SIMD3(1, 2, 1))])
        let p = try ProjectedPotential(structure: s, nx: 96, ny: 128, sliceThickness: 2)
        let beam = try ElectronBeam(energy: 200e3)
        let nearly = try Grid(nx: 96, ny: 128, extentX: 9.600000000000001, extentY: 12.8)
        #expect(nearly != p.grid)
        let a = try Multislice.run(.planeWave(beam: beam, grid: nearly), through: p)
        let b = try Multislice.run(.planeWave(beam: beam, grid: p.grid), through: p)
        #expect(a[0].field == b[0].field)
        let other = try Grid(nx: 96, ny: 128, extentX: 9.7, extentY: 12.8)
        #expect(throws: SimulatoError.gridMismatch) { try Multislice.run(.planeWave(beam: beam, grid: other), through: p) }
    }

    // MARK: 5, the symmetry of Al [001]

    /// Truth 5 (amendment 2): slices of a/2 are C-centred, so the pattern lives on h + k even; one slice per cell has
    /// period a/2, so the pattern lives on h and k both even; 4mm symmetry either way. Units of 1/a = 4 pixels (4 cells).
    @Test func aluminiumPatternLivesOnTheSliceLatticeWith4mmSymmetry() throws {
        let beam = try Pins.beam("B1")
        let structure = try Pins.structure("B1")
        for (sliceThickness, bothEven) in [(2.025, false), (4.05, true)] {
            let potential = try ProjectedPotential(structure: structure, nx: 256, ny: 256, sliceThickness: sliceThickness)
            let planes: ExitPlanes = sliceThickness < 3 ? .every(2) : .every(1)
            let waves = try Multislice.run(.planeWave(beam: beam, grid: potential.grid), through: potential,
                                           exitPlanes: planes)
            let n = 256, unit = 4, c = 128
            for w in waves.dropFirst() {
                let p = try w.diffractionPattern(maxAngle: .full)
                let peak = p.max
                var offLattice = 0.0, asymmetry = 0.0
                for ix in 0..<n {
                    for iy in 0..<n {
                        let h = ix - c, k = iy - c
                        let onUnit = h % unit == 0 && k % unit == 0
                        let allowed = onUnit && (bothEven ? (h / unit) % 2 == 0 && (k / unit) % 2 == 0
                                                          : ((h + k) / unit) % 2 == 0)
                        if !allowed { offLattice = max(offLattice, p[ix, iy]) }
                        if h > -c && k > -c {
                            asymmetry = max(asymmetry, abs(p[ix, iy] - p[c + k, c + h]), abs(p[ix, iy] - p[c - h, iy]))
                        }
                    }
                }
                #expect(offLattice <= 1e-20 * peak, "slices of \(sliceThickness) Å: off-lattice intensity")
                #expect(asymmetry <= 1e-12 * peak, "slices of \(sliceThickness) Å: 4mm symmetry")
            }
        }
    }

    // MARK: 6, the probe

    @Test func probeNormalisationSymmetryAndShift() throws {
        let beam = try ElectronBeam(energy: 200e3)
        let grid = try Grid(nx: 128, ny: 128, extentX: 12.8, extentY: 12.8)
        let aperture = try Aperture(semiangleCutoff: 20)
        let plain = try Probe(beam: beam, grid: grid, aperture: aperture)
        let centred = try plain.build(at: SIMD2(6.4, 6.4))
        #expect(abs((try centred.diffractionPattern(maxAngle: .full).sum) - (1)) <= 1e-14, "Σ_k |Ψ|² = 1")
        let intensity = centred.intensity()
        var asym = 0.0
        for dx in -40...40 {
            for dy in -40...40 {
                asym = max(asym, abs(intensity[64 + dx, 64 + dy] - intensity[64 - dx, 64 - dy]))
            }
        }
        #expect(asym <= 1e-12 * intensity.max, "centrosymmetric about a pixel position")

        var aberrations = Aberrations()
        aberrations.defocus = 50
        aberrations.Cs = 1e4
        aberrations.C23 = 300
        aberrations.phi23 = 0.4
        let aberrated = try Probe(beam: beam, grid: grid, aperture: aperture, aberrations: aberrations)
        let base = try aberrated.build(at: SIMD2(3.05, 4.15))
        let (p, q) = (7, -5)
        let moved = try aberrated.build(at: SIMD2(3.05 + Double(p) * grid.samplingX, 4.15 + Double(q) * grid.samplingY))
        var shifted = ComplexField(nx: 128, ny: 128)
        for ix in 0..<128 {
            for iy in 0..<128 {
                let sx = ((ix - p) % 128 + 128) % 128, sy = ((iy - q) % 128 + 128) % 128
                shifted.real[ix * 128 + iy] = base.field.real[sx * 128 + sy]
                shifted.imaginary[ix * 128 + iy] = base.field.imaginary[sx * 128 + sy]
            }
        }
        #expect(relMax(moved.field, shifted) <= 1e-12, "a whole-pixel move is a circular shift")
    }

    // MARK: 7, χ

    @Test func chiEqualsTheClosedForm() {
        var a = Aberrations()
        a.C10 = -73.5
        a.C12 = 12.25
        a.phi12 = 0.7
        a.C23 = 410
        a.phi23 = -0.35
        a.C30 = 2.1e4
        let lambda = 0.025079340317328468
        for (alpha, phi) in [(0.004, 0.3), (0.012, -2.1), (0.027, 1.4), (0.0, 0.0)] {
            let hand = 2 * Double.pi / lambda * (alpha * alpha / 2 * (-73.5 + 12.25 * cos(2 * (phi - 0.7)))
                                                 + alpha * alpha * alpha / 3 * (410 * cos(3 * (phi + 0.35)))
                                                 + alpha * alpha * alpha * alpha / 4 * 2.1e4)
            let chi = a.chi(alpha: alpha, phi: phi, wavelength: lambda)
            #expect(abs(chi - hand) <= 1e-12 * max(abs(hand), 1), "χ at α \(alpha), φ \(phi)")
        }
        #expect(a.defocus == 73.5)
    }

    // MARK: 8, the thickness series

    @Test func thicknessSeriesIsConsistentMonotoneAndConservesIntensity() throws {
        let full = try Pins.structure("B1")
        let potential = try Pins.potential("B1", pi: .exact)
        let beam = try Pins.beam("B1")
        let series = try Multislice.run(.planeWave(beam: beam, grid: potential.grid), through: potential,
                                        exitPlanes: .every(2))
        #expect(series.count == 11)
        for cells in [1, 3, 7] {
            let depth = Double(cells) * 4.05
            let atoms = full.atoms.filter { $0.position.z < depth - 1e-9 }
            let part = try AtomicStructure(box: SIMD3(full.box.x, full.box.y, depth), atoms: atoms)
            let p = try ProjectedPotential(structure: part, nx: 256, ny: 256, sliceThickness: 2.025)
            let w = try Multislice.run(.planeWave(beam: beam, grid: p.grid), through: p)
            #expect(relMax(w[0].field, series[cells].field) <= 1e-13, "\(cells) cells run alone")
        }
        let n = Double(potential.grid.count)
        for q in series.indices {
            if q > 0 {
                #expect(series[q].totalIntensity <= series[q - 1].totalIntensity * (1 + 1e-13), "Σ|ψ|² does not increase (plane \(q))")
            }
            let pattern = try series[q].diffractionPattern(maxAngle: .full)
            #expect(relDiff(pattern.sum, series[q].totalIntensity / n) <= 1e-12, "Parseval, plane \(q)")
            #expect(pattern[128, 128] <= pattern.sum * (1 + 1e-12), "central beam within the total")
        }
    }
}

//
//  CoreChecks.swift
//  Role: SimulatoCore from outside: the beam, constants, errors, grid, antialias, fields, waves and structures. Every
//        number is one a reader can recompute by hand (the comment says how).
//

import Foundation
import SimulatoCore

/// λ at 200 keV [Å] as the kit computes it (CODATA 2014 constants, abTEM's evaluation order).
let lambda200 = 0.025079340317328468

func checkCore(_ r: inout Report) {
    r.section("Core: beam") { r in
        r.equal("SimulatoCore.referenceVersion", SimulatoCore.referenceVersion, "abTEM 1.0.10 (164e644)")
        let beam = try ElectronBeam(energy: 200_000)
        r.equal("ElectronBeam.energy", beam.energy, 200_000)
        r.near("wavelength at 200 keV = 0.025079340317328468 Å", beam.wavelength, lambda200, 1e-16)
        // Textbook: λ = 12.2643 / √(V (1 + 0.97848e-6 V)) Å, V in volts.
        let v = 200_000.0
        r.nearRelative("wavelength against 12.2643 / sqrt(V (1 + 0.97848e-6 V))", beam.wavelength,
                       12.2643 / (v * (1 + 0.97848e-6 * v)).squareRoot(), 2e-5)
        // γ = 1 + E / (m c²), m c² = 510998.95 eV.
        r.nearRelative("relativisticMassFactor = 1 + E / 510998.95 eV", beam.relativisticMassFactor,
                       1 + v / 510_998.95, 1e-8)
        // σ = 2π / (λ E) · (m c² + E) / (2 m c² + E) [1/(V Å)] (Kirkland); 7.288e-4 at 200 keV.
        let sigma = 2 * Double.pi / (beam.wavelength * v) * (510_998.95 + v) / (2 * 510_998.95 + v)
        r.nearRelative("interactionParameter against Kirkland's closed form", beam.interactionParameter, sigma, 1e-7)
        r.check("interactionParameter is 7.29e-4 rad/(V Å) at 200 keV",
                abs(beam.interactionParameter - 7.29e-4) < 1e-6, "got \(beam.interactionParameter)")
        r.check("ElectronBeam is Hashable", try beam == ElectronBeam(energy: 200_000)
                && beam != ElectronBeam(energy: 300_000))
        r.check("higher energy, shorter wavelength", try ElectronBeam(energy: 300_000).wavelength < beam.wavelength)
        r.refuses("ElectronBeam(energy: 0)", .invalidEnergy(0)) { try ElectronBeam(energy: 0) }
        r.refuses("ElectronBeam(energy: -1)", .invalidEnergy(-1)) { try ElectronBeam(energy: -1) }
        r.refusesCase("ElectronBeam(energy: .nan)", "invalidEnergy") { try ElectronBeam(energy: .nan) }
        r.refusesCase("ElectronBeam(energy: .infinity)", "invalidEnergy") { try ElectronBeam(energy: .infinity) }
    }

    r.section("Core: constants") { r in
        let c = PhysicalConstants.self
        r.equal("speedOfLight", c.speedOfLight, 299_792_458.0)
        r.equal("planck (CODATA 2014)", c.planck, 6.62607004e-34)
        r.equal("electronMass (CODATA 2014)", c.electronMass, 9.10938356e-31)
        r.equal("elementaryCharge (CODATA 2014)", c.elementaryCharge, 1.6021766208e-19)
        r.equal("vacuumPermittivity (CODATA 2014)", c.vacuumPermittivity, 8.85418781762039e-12)
        r.equal("bohrRadius [Å]", c.bohrRadius, 0.5291772105638411)
        // κ = 4π ε₀ / (2π a₀ e) = 2 ε₀ / (a₀ e) with ε₀ in F/m, a₀ in Å (1e-10 m), e in C: 0.0208866.
        r.nearRelative("kappa = 2 eps0 / (a0 e)", c.kappa, 2 * c.vacuumPermittivity * 1e-10 / (c.bohrRadius * c.elementaryCharge),
                       1e-7)
    }

    r.section("Core: errors") { r in
        let all: [SimulatoError] = [
            .invalidEnergy(-1), .unsupportedGridLength(7), .invalidSampling(0, 1), .invalidTilt("x"), .invalidExtent(0, 1), .invalidBox(0, 1, 1), .unsupportedElement(0),
            .invalidPosition(index: 3), .invalidSliceThickness(0), .sliceThicknessesDoNotSumToBox(sum: 1, depth: 2),
            .invalidExitPlanes("x"), .gridMismatch, .invalidAperture(0), .invalidAberration("x"), .metalUnavailable("x"),
            .metalFailure("x"), .invalidFrozenPhonons("x"), .invalidScan("x"), .invalidDetector("x"),
        ]
        r.equal("SimulatoError: 19 cases built", all.count, 19)
        r.check("SimulatoError: every description is non-empty and distinct",
                Set(all.map(\.description)).count == 19 && all.allSatisfy { !$0.description.isEmpty })
        r.equal("SimulatoError.invalidEnergy description", SimulatoError.invalidEnergy(-1).description,
                "electron energy must be finite and positive, got -1.0 eV")
        r.check("SimulatoError is Equatable", SimulatoError.invalidPosition(index: 3) == .invalidPosition(index: 3)
                && SimulatoError.invalidPosition(index: 3) != .invalidPosition(index: 4))
        let boxed: any Error = SimulatoError.gridMismatch
        r.check("SimulatoError is an Error and Sendable", boxed is SimulatoError)
    }

    r.section("Core: grid") { r in
        let g = try Grid(nx: 32, ny: 48, extentX: 8, extentY: 12)
        r.check("Grid stores nx, ny, extents", g.nx == 32 && g.ny == 48 && g.extentX == 8 && g.extentY == 12)
        r.check("Grid sampling = extent / gpts", g.samplingX == 0.25 && g.samplingY == 0.25)
        r.equal("Grid.count", g.count, 32 * 48)
        r.near("Grid.reciprocalSamplingX = 1 / extent", g.reciprocalSamplingX, 1.0 / 8, 1e-15)
        r.near("Grid.reciprocalSamplingY = 1 / extent", g.reciprocalSamplingY, 1.0 / 12, 1e-15)
        // fftfreq(8, 0.5): [0, 1, 2, 3, -4, -3, -2, -1] / (8 · 0.5)
        r.equal("Grid.fftfreq(8, 0.5)", Grid.fftfreq(8, 0.5), [0, 0.25, 0.5, 0.75, -1, -0.75, -0.5, -0.25])
        // odd n: [0, 1, 2, -2, -1] / 5
        r.equal("Grid.fftfreq(5, 1) (odd n)", Grid.fftfreq(5, 1), [0, 0.2, 0.4, -0.4, -0.2])
        r.equal("Grid.frequenciesX = fftfreq(nx, samplingX)", g.frequenciesX, Grid.fftfreq(32, 0.25))
        r.equal("Grid.frequenciesY = fftfreq(ny, samplingY)", g.frequenciesY, Grid.fftfreq(48, 0.25))
        r.equal("Grid.frequenciesX[nx / 2] = −Nyquist", g.frequenciesX[16], -2.0)
        // any n ≥ 2 (0.3.0); hasDirectTransform: f · 2^m, f in {1, 3, 5, 15}, m ≥ 3
        let supported = [2, 4, 7, 8, 9, 12, 16, 24, 45, 64, 100, 120, 360, 480, 512, 960]
        let unsupported = [-8, 0, 1]
        r.check("Grid.hasDirectTransform: vDSP's lengths only", [8, 24, 120, 960].allSatisfy(Grid.hasDirectTransform)
                && [4, 7, 12, 100, 360].allSatisfy { !Grid.hasDirectTransform($0) })
        r.check("Grid.isSupportedLength accepts \(supported)", supported.allSatisfy(Grid.isSupportedLength))
        r.check("Grid.isSupportedLength refuses \(unsupported)", unsupported.allSatisfy { !Grid.isSupportedLength($0) })
        let nearly = try Grid(nx: 32, ny: 48, extentX: 8 * (1 + 1e-13), extentY: 12)
        let wider = try Grid(nx: 32, ny: 48, extentX: 8.001, extentY: 12)
        let coarser = try Grid(nx: 40, ny: 48, extentX: 8, extentY: 12)
        r.check("Grid.matches: itself and within 1e-12 relative", g.matches(g) && g.matches(nearly))
        r.check("Grid.matches refuses another extent or size", !g.matches(wider) && !g.matches(coarser))
        r.check("Grid is Hashable", g == (try Grid(nx: 32, ny: 48, extentX: 8, extentY: 12)) && Set([g, g]).count == 1)
        r.refuses("Grid nx = 1", .unsupportedGridLength(1)) { try Grid(nx: 1, ny: 8, extentX: 1, extentY: 1) }
        r.refuses("Grid ny = 0", .unsupportedGridLength(0)) { try Grid(nx: 8, ny: 0, extentX: 1, extentY: 1) }
        r.check("Grid nx = 7, ny = 12 (any length)", (try? Grid(nx: 7, ny: 12, extentX: 1, extentY: 1)) != nil)
        let fromSampling = try Grid(extentX: 8, extentY: 12, samplingX: 0.3, samplingY: 0.25)
        r.check("Grid(extent:sampling:) = ceil(extent / sampling)", fromSampling.nx == 27 && fromSampling.ny == 48)
        r.refusesCase("Grid sampling 0", "invalidSampling") { try Grid(extentX: 8, extentY: 12, samplingX: 0, samplingY: 1) }
        r.refuses("Grid extent 0", .invalidExtent(0, 1)) { try Grid(nx: 8, ny: 8, extentX: 0, extentY: 1) }
        r.refusesCase("Grid extent NaN", "invalidExtent") { try Grid(nx: 8, ny: 8, extentX: 1, extentY: .nan) }
        r.refusesCase("Grid extent negative", "invalidExtent") { try Grid(nx: 8, ny: 8, extentX: -1, extentY: 1) }
    }

    r.section("Core: antialias") { r in
        r.equal("Antialias.cutoff = 0.6666666", Antialias.cutoff, 0.6666666)
        r.equal("Antialias.taper = 0.01", Antialias.taper, 0.01)
        // 32 × 32 over 32 Å: Δ = 1, k = i / 32; k_c = 0.6666666 / 2; w = 0.01. 1 where |(i, j)| / 32 ≤ k_c − w.
        let g = try Grid(nx: 32, ny: 32, extentX: 32, extentY: 32)
        let a = Antialias.aperture(on: g)
        r.equal("aperture: one value per pixel", a.count, 1024)
        r.equal("aperture at k = 0", a[0], 1)
        r.equal("aperture at the Nyquist corner (16, 16)", a[16 * 32 + 16], 0)
        var inside = 0, ones = 0, between = 0
        let kc = 0.6666666 / 2
        for i in -16..<16 {
            for j in -16..<16 {
                let rr = Double(i * i + j * j).squareRoot() / 32
                if rr <= kc - 0.01 { inside += 1 } else if rr <= kc { between += 1 }
                let v = a[((i + 32) % 32) * 32 + (j + 32) % 32]
                if v == 1 { ones += 1 }
            }
        }
        r.equal("aperture: pixels at 1 = pixels inside k_c − w (counted by hand)", ones, inside)
        r.check("aperture: a taper exists and takes values in (0, 1)", between > 0
                && a.contains { $0 > 0 && $0 < 1 } && a.allSatisfy { $0 >= 0 && $0 <= 1 })
        // ⌊round(2/3 · 32, 7)⌋ = 21; parity of 32 → 22; even → 23.
        let n = Antialias.cutoffPatternSize(on: g)
        r.check("cutoffPatternSize on 32² (21 → 22 → 23)", n.nx == 23 && n.ny == 23, "got \(n)")
        // 2/3 · 64 = 42.67 → 42; parity of 64 is even; even → 43
        r.check("cutoffPatternSize on 64²",
                Antialias.cutoffPatternSize(on: try Grid(nx: 64, ny: 64, extentX: 64, extentY: 64)) == (43, 43))
        // the coarser axis decides: 32 pts over 32 Å (Δ 1) and 32 over 16 Å (Δ 0.5): k_cut = 2/3 for both.
        let m = Antialias.cutoffPatternSize(on: try Grid(nx: 32, ny: 32, extentX: 32, extentY: 16))
        r.check("cutoffPatternSize uses the coarser sampling for both axes", m.nx == 23 && m.ny == 11, "got \(m)")
    }

    r.section("Core: fields") { r in
        var f = ComplexField(nx: 2, ny: 3)
        r.check("ComplexField(nx:ny:) is zeros", f.real == [Double](repeating: 0, count: 6) && f.imaginary == f.real)
        f = ComplexField(nx: 2, ny: 3, real: [3, 0, 0, 0, 0, 1], imaginary: [4, 0, 0, 0, 0, 0])
        r.check("ComplexField subscript [ix * ny + iy]", f[0, 0] == (3, 4) && f[1, 2] == (1, 0))
        r.equal("ComplexField.squaredMagnitude", f.squaredMagnitude(), [25, 0, 0, 0, 0, 1])
        r.equal("ComplexField.maxMagnitude (3-4-5)", f.maxMagnitude, 5)
        f.real[1] = 2
        r.equal("ComplexField.real is mutable", f.real[1], 2)
        r.check("ComplexField is Hashable and Equatable", f == f && f != ComplexField(nx: 2, ny: 3))
        let rf = RealField(nx: 2, ny: 2, values: [1, 5, -3, 2], samplingX: 0.5, samplingY: 0.25)
        r.check("RealField stores values and sampling", rf.nx == 2 && rf.ny == 2 && rf.samplingX == 0.5 && rf.samplingY == 0.25)
        r.equal("RealField subscript [ix, iy]", rf[1, 0], -3)
        r.equal("RealField.sum", rf.sum, 5)
        r.equal("RealField.max", rf.max, 5)
        var rm = rf
        rm.values[0] = 9
        r.check("RealField.values is mutable, RealField is Hashable", rm.max == 9 && rm != rf)
    }

    r.section("Core: waves") { r in
        let beam = try ElectronBeam(energy: 200_000)
        let g = try Grid(nx: 32, ny: 32, extentX: 32, extentY: 32)
        let w = Waves.planeWave(beam: beam, grid: g)
        r.check("planeWave: ψ = 1 + 0i at every pixel", w.field.real.allSatisfy { $0 == 1 } && w.field.imaginary.allSatisfy { $0 == 0 })
        r.check("planeWave: .values normalisation, its grid and beam", w.normalization == .values && w.grid == g && w.beam == beam)
        r.equal("Waves.totalIntensity of a plane wave = nx ny", w.totalIntensity, 1024)
        let i = w.intensity()
        r.check("Waves.intensity: ones, sampled as the grid", i.values.allSatisfy { $0 == 1 } && i.samplingX == 1 && i.nx == 32)
        // ψ / N, FFT: all intensity at k = 0, value 1, at (nx/2, ny/2) of the shifted pattern.
        let full = try w.diffractionPattern(maxAngle: .full)
        r.check("diffractionPattern(.full): 32 × 32, sum 1, peak 1 at (16, 16)",
                full.nx == 32 && abs(full.sum - 1) < 1e-12 && abs(full[16, 16] - 1) < 1e-12)
        r.near("diffractionPattern sampling = 1 / extent", full.samplingX, 1.0 / 32, 1e-15)
        let cut = try w.diffractionPattern()
        r.check("diffractionPattern(): cropped to the cutoff, 23 × 23, peak 1 at (11, 11)",
                cut.nx == 23 && cut.ny == 23 && abs(cut[11, 11] - 1) < 1e-12 && abs(cut.sum - 1) < 1e-12)
        let probeLike = Waves(field: ComplexField(nx: 32, ny: 32, real: [Double](repeating: 1, count: 1024),
                                                  imaginary: [Double](repeating: 0, count: 1024)),
                              grid: g, beam: beam, normalization: .reciprocalSpace)
        // .reciprocalSpace: no division by N, so the pattern of ψ = 1 sums to N² = 1024² (k = 0 carries N²).
        let pattern = try probeLike.diffractionPattern(maxAngle: .full)
        r.near("diffractionPattern with .reciprocalSpace is not divided by N (peak N²)", pattern[16, 16], 1024.0 * 1024.0, 1e-6)
        r.check("Waves.Normalization and Waves.MaxAngle cases", Waves.Normalization.values != .reciprocalSpace
                && Waves.MaxAngle.cutoff != .full)
        r.check("Waves is Hashable", w == Waves.planeWave(beam: beam, grid: g) && w != probeLike)
    }

    r.section("Core: structures") { r in
        let atom = AtomicStructure.Atom(atomicNumber: 13, position: SIMD3(1, 2, 3))
        r.check("AtomicStructure.Atom stores Z and position", atom.atomicNumber == 13 && atom.position == SIMD3(1, 2, 3))
        let s = try AtomicStructure(box: SIMD3(4, 5, 6), atoms: [atom])
        r.check("AtomicStructure stores box and atoms", s.box == SIMD3(4, 5, 6) && s.atoms == [atom])
        // ASE order: x outermost, then y, then z; (1, 2, 3) shifted by multiples of the box (4, 5, 6).
        let rep = s.repeated(2, 1, 2)
        r.check("repeated(2, 1, 2): box doubles along x and z", rep.box == SIMD3(8, 5, 12))
        r.equal("repeated(2, 1, 2): positions in ASE order", rep.atoms.map(\.position),
                [SIMD3(1, 2, 3), SIMD3(1, 2, 9), SIMD3(5, 2, 3), SIMD3(5, 2, 9)])
        r.check("AtomicStructure and Atom are Hashable", s == (try AtomicStructure(box: SIMD3(4, 5, 6), atoms: [atom]))
                && Set([atom, atom]).count == 1)
        r.refusesCase("AtomicStructure box with a zero side", "invalidBox") { try AtomicStructure(box: SIMD3(0, 1, 1), atoms: []) }
        r.refusesCase("AtomicStructure box with NaN", "invalidBox") { try AtomicStructure(box: SIMD3(1, .nan, 1), atoms: []) }
        r.refuses("AtomicStructure atom 1 with a non-finite position", .invalidPosition(index: 1)) {
            try AtomicStructure(box: SIMD3(4, 5, 6), atoms: [atom, .init(atomicNumber: 13, position: SIMD3(.infinity, 0, 0))])
        }
    }
}

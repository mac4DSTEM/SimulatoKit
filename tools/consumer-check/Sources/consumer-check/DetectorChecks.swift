//
//  DetectorChecks.swift
//  Role: SimulatoDetector from outside: DiffractionIntensity, the three detectors, the public `Detector` protocol
//        adopted by a detector defined here, and the STEM driver with frozen-phonon averaging.
//

import Foundation
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe
import SimulatoDetector

/// A detector defined outside the kit: the sum of everything (shape []), conforming to the public protocol.
struct TotalDetector: Detector {
    func shape(for grid: Grid, beam: ElectronBeam) throws -> [Int] { [] }
    func detect(_ intensity: DiffractionIntensity) throws -> [Double] { [intensity.values.reduce(0, +)] }
}

/// A custom detector with a two-number value: the pattern's centre-pixel and the pixel count, in a [2] shape.
struct CentreAndCountDetector: Detector {
    func shape(for grid: Grid, beam: ElectronBeam) throws -> [Int] { [2] }
    func detect(_ intensity: DiffractionIntensity) throws -> [Double] { [intensity.values[0], Double(intensity.grid.count)] }
}

/// Compiles only while the type conforms to the public `Detector` protocol from outside the package.
func isDetector(_ d: any Detector) -> Bool { true }

func checkDetector(_ r: inout Report) {
    r.section("Detector: intensity") { r in
        let beam = try ElectronBeam(energy: 200_000)
        let g = try Grid(nx: 32, ny: 32, extentX: 8, extentY: 8)
        let probe = try Probe(beam: beam, grid: g, aperture: Aperture(semiangleCutoff: 20)).build(at: SIMD2(4, 4))
        let p = try DiffractionIntensity(probe)
        r.check("DiffractionIntensity keeps grid, beam; 1024 values", p.grid == g && p.beam == beam && p.values.count == 1024)
        r.near("a probe's intensity sums to 1 (Σ_k |Ψ|² = 1)", p.values.reduce(0, +), 1, 1e-12)
        let plane = try DiffractionIntensity(Waves.planeWave(beam: beam, grid: g))
        r.check("a plane wave's intensity is 1 at k = 0 only (ψ / N)", abs(plane.values[0] - 1) < 1e-12
                && plane.values.dropFirst().allSatisfy { $0 < 1e-20 })
        // Δα = (1/8 Å⁻¹) λ · 10³ mrad
        r.near("angularSampling = λ / L · 10³ = 3.135 mrad", p.angularSampling.x, lambda200 / 8 * 1e3, 1e-12)
        r.equal("angularSampling x = y on a square grid", p.angularSampling.x, p.angularSampling.y)
        // cutoff gpts: ⌊2/3 · 8 / 0.25⌋ = ⌊21.33⌋ = 21 → grid parity (32) → 22; angle = ⌊22/2⌋ Δα = 11 Δα
        r.near("cutoffAngles = 11 Δα (the cutoff gpts 22 at the grid's parity)", p.cutoffAngles.x, 11 * p.angularSampling.x, 1e-12)
    }

    r.section("Detector: annular") { r in
        let beam = try ElectronBeam(energy: 200_000)
        let g = try Grid(nx: 32, ny: 32, extentX: 8, extentY: 8)
        let p = try DiffractionIntensity(Probe(beam: beam, grid: g, aperture: Aperture(semiangleCutoff: 20)).build(at: SIMD2(4, 4)))
        let da = p.angularSampling.x
        let maximum = 16 * da
        let all = try AnnularDetector(inner: 0, outer: maximum)
        let defaults = try AnnularDetector()
        r.check("AnnularDetector stores inner and outer; defaults 0 and nil", all.inner == 0 && all.outer == maximum
                && defaults.inner == 0 && defaults.outer == nil)
        r.equal("shape is [] (one number)", try all.shape(for: g, beam: beam), [])
        r.near("over all angles the annular detector equals the total (1)", try all.detect(p)[0], 1, 1e-12)
        r.equal("outerAngle(p) is the given outer", try all.outerAngle(p), maximum)
        // outer = nil → ⌊cutoff angle⌋ = ⌊11 · 3.1349⌋ = 34 mrad
        let auto = try AnnularDetector()
        r.equal("outerAngle(p) with outer nil = ⌊min cutoff angle⌋ = 34 mrad", try auto.outerAngle(p), 34)
        r.near("annular with outer nil also gets everything of a 20 mrad probe", try auto.detect(p)[0], 1, 1e-12)
        let bf = try AnnularDetector(inner: 0, outer: 10), adf = try AnnularDetector(inner: 10, outer: maximum)
        r.near("BF (0–10) + ADF (10–max) = total", try bf.detect(p)[0] + adf.detect(p)[0], 1, 1e-12)
        // a flat 20 mrad disk in k: the share inside 10 mrad is the area ratio (10/20)² = 0.25, up to pixelation
        let share = try bf.detect(p)[0]
        r.check("the share of a flat 20 mrad disk inside 10 mrad is the area ratio 0.25 (0.15 ... 0.35)", share > 0.15 && share < 0.35, "got \(share)")
        r.check("AnnularDetector is a Detector and Hashable", isDetector(all) && all == all && all != adf)
        r.refusesCase("outer beyond the maximum simulated angle", "invalidDetector") { try AnnularDetector(inner: 0, outer: maximum * 1.5).detect(p) }
        r.refusesCase("inner negative", "invalidDetector") { try AnnularDetector(inner: -1) }
        r.refusesCase("inner NaN", "invalidDetector") { try AnnularDetector(inner: .nan) }
        r.refusesCase("outer below inner", "invalidDetector") { try AnnularDetector(inner: 10, outer: 5) }
        r.refusesCase("inner above the automatic outer (34 mrad)", "invalidDetector") { try AnnularDetector(inner: 50).detect(p) }
    }

    r.section("Detector: pixelated and segmented") { r in
        let beam = try ElectronBeam(energy: 200_000)
        let g = try Grid(nx: 32, ny: 32, extentX: 8, extentY: 8)
        let p = try DiffractionIntensity(Probe(beam: beam, grid: g, aperture: Aperture(semiangleCutoff: 20)).build(at: SIMD2(4, 4)))
        // .full 32 × 32; .cutoff 22 × 22 (grid parity); .valid ⌊22 / √2⌋ = 15 → parity of 32 → 16
        let full = PixelatedDetector(maxAngle: .full), cutoff = PixelatedDetector(maxAngle: .cutoff), valid = PixelatedDetector()
        r.equal("PixelatedDetector default is .valid", valid.maxAngle, .valid)
        r.equal("shape(.full)", try full.shape(for: g, beam: beam), [32, 32])
        r.equal("shape(.cutoff)", try cutoff.shape(for: g, beam: beam), [22, 22])
        r.equal("shape(.valid)", try valid.shape(for: g, beam: beam), [16, 16])
        r.near("the full pattern sums to 1", try full.detect(p).reduce(0, +), 1, 1e-12)
        r.equal("detect(.valid) has 16 × 16 values", try valid.detect(p).count, 256)
        r.near("the valid pattern still holds the whole 20 mrad probe", try valid.detect(p).reduce(0, +), 1, 1e-12)
        r.near("the cutoff pattern holds it too", try cutoff.detect(p).reduce(0, +), 1, 1e-12)
        // fftshift: the zero frequency is at index (n/2, n/2) = (8, 8) of the 16 × 16 pattern
        let v = try valid.detect(p)
        r.check("the pattern is fftshifted: the brightest pixel is not at the corner", v.indices.max { v[$0] < v[$1] }! != 0)
        r.check("PixelatedDetector is a Detector and Hashable; MaxAngle has three cases", full != cutoff && isDetector(valid)
                && PixelatedDetector.MaxAngle.valid != .cutoff)

        let seg = try SegmentedDetector(radialBins: 2, azimuthalBins: 4, inner: 0, outer: 25)
        r.check("SegmentedDetector stores its parameters; rotation defaults to 0", seg.radialBins == 2 && seg.azimuthalBins == 4
                && seg.inner == 0 && seg.outer == 25 && seg.rotation == 0)
        r.equal("shape is [radialBins, azimuthalBins]", try seg.shape(for: g, beam: beam), [2, 4])
        let bins = try seg.detect(p)
        r.near("the 8 bins between 0 and 25 mrad hold the whole 20 mrad probe", bins.reduce(0, +), 1, 1e-12)
        // flat 20 mrad disk: the ring 0 ... 12.5 mrad holds the area ratio (12.5/20)² = 0.39, up to pixelation
        let inner = bins[0..<4].reduce(0, +)
        r.check("bins are non-negative; the inner ring (0 ... 12.5 mrad) holds the area ratio 0.39 (0.25 ... 0.55)",
                bins.allSatisfy { $0 >= 0 } && inner > 0.25 && inner < 0.55, "got \(inner)")
        let rotated = try SegmentedDetector(radialBins: 2, azimuthalBins: 4, inner: 0, outer: 25, rotation: Double.pi / 4)
        r.near("a rotated segmented detector holds the same total", try rotated.detect(p).reduce(0, +), 1, 1e-12)
        r.check("SegmentedDetector is a Detector and Hashable", isDetector(seg) && seg == seg && seg != rotated)
        r.refusesCase("0 radial bins", "invalidDetector") { try SegmentedDetector(radialBins: 0, azimuthalBins: 4, inner: 0, outer: 25) }
        r.refusesCase("outer = inner", "invalidDetector") { try SegmentedDetector(radialBins: 1, azimuthalBins: 4, inner: 5, outer: 5) }
        r.refusesCase("NaN rotation", "invalidDetector") {
            try SegmentedDetector(radialBins: 1, azimuthalBins: 4, inner: 0, outer: 5, rotation: .nan)
        }
    }

    r.section("Detector: STEM driver") { r in
        let beam = try ElectronBeam(energy: 200_000)
        let empty = try AtomicStructure(box: SIMD3(8, 8, 4), atoms: [])
        let vacuum = try ProjectedPotential(structure: empty, nx: 32, ny: 32, sliceThickness: 2)
        let probe = try Probe(beam: beam, grid: vacuum.grid, aperture: Aperture(semiangleCutoff: 20))
        let scan = try GridScan(start: SIMD2(1, 1), end: SIMD2(5, 5), nx: 2, ny: 2)
        let custom = CentreAndCountDetector()
        let ms = try STEM.scan(probe: probe, potentials: [vacuum], scan: scan,
                               detectors: [AnnularDetector(), TotalDetector(), custom, PixelatedDetector(maxAngle: .cutoff)])
        r.equal("one measurement per detector", ms.count, 4)
        let typed: [ScanMeasurement] = ms  // the bare name, beside Foundation's Measurement (renamed for 0.1.0)
        r.equal("the result type is ScanMeasurement", String(reflecting: type(of: typed[0])), "SimulatoDetector.ScanMeasurement")
        // A consumer can build an intensity and test its own detector on it.
        let given = DiffractionIntensity(grid: vacuum.grid, beam: beam, values: [Double](repeating: 1, count: vacuum.grid.count))
        r.equal("DiffractionIntensity(grid:beam:values:): a custom detector on a known intensity",
                try TotalDetector().detect(given), [Double(vacuum.grid.count)])
        let lazy = try Potential(structure: empty, nx: 32, ny: 32, sliceThickness: 2)
        let viaLazy = try STEM.scan(probe: probe, potentials: [lazy], scan: scan, detectors: [AnnularDetector()])
        r.check("STEM.scan through a Potential equals the scan through its ProjectedPotential",
                maxDifference(viaLazy[0].values, ms[0].values) == 0)
        if Multislice.metalIsAvailable {
            let gpu = try STEM.scan(probe: probe, potentials: [lazy], scan: scan, detectors: [AnnularDetector()], backend: .metal)
            r.check("STEM.scan through a Potential on .metal: vacuum reads 1 within 1e-5", gpu[0].values.allSatisfy { abs($0 - 1) < 1e-5 })
        } else {
            r.skip("STEM.scan through a Potential on .metal", "no Metal device")
        }
        r.check("scanShape [2, 2] and the value shapes", ms.allSatisfy { $0.scanShape == [2, 2] }
                && ms[0].valueShape == [] && ms[1].valueShape == [] && ms[2].valueShape == [2] && ms[3].valueShape == [22, 22])
        r.equal("values: scan positions × the detector's numbers", ms.map(\.values.count), [4, 4, 8, 4 * 22 * 22])
        // vacuum: the 20 mrad probe stays inside the antialias cutoff, so every detector sees all of it
        r.check("vacuum: the annular detector reads 1 at every position", ms[0].values.allSatisfy { abs($0 - 1) < 1e-9 })
        r.check("vacuum: the custom total detector (a user type) agrees with the annular one", maxDifference(ms[0].values, ms[1].values) < 1e-12)
        r.check("custom two-number detector: [centre pixel, 1024] per position", ms[2][1, 0].count == 2 && ms[2][1, 0][1] == 1024
                && ms[2][0, 0][0] > 0 && ms[2][1, 1][0] == ms[2].values[6])
        r.equal("subscript(ix, iy) of a one-number detector", ms[0][1, 1].count, 1)
        let image = try unwrap(ms[0].image(samplingX: 4, samplingY: 4))
        r.check("image(samplingX:samplingY:) of a one-number detector: 2 × 2, sampled as given, sum 4", image.nx == 2 && image.ny == 2
                && image.samplingX == 4 && abs(image.sum - 4) < 1e-8)
        r.check("image of a pixelated detector is nil", ms[3].image(samplingX: 1, samplingY: 1) == nil)

        // frozen phonons with σ = 0: the average over configurations equals the static scan
        let al = try AtomicStructure(box: SIMD3(8, 8, 4), atoms: [.init(atomicNumber: 13, position: SIMD3(3, 3, 2))])
        let staticPotential = try ProjectedPotential(structure: al, nx: 32, ny: 32, sliceThickness: 2)
        let phonons = try FrozenPhonons(structure: al, sigmas: [13: 0], configurations: 2, seed: 3)
        let single = try STEM.scan(probe: probe, potentials: [staticPotential], scan: scan, detectors: [AnnularDetector()])
        let averaged = try STEM.scan(probe: probe, potentials: phonons.potentials(nx: 32, ny: 32, sliceThickness: 2), scan: scan,
                                     detectors: [AnnularDetector()])
        r.check("σ = 0 frozen phonons: the mean over 2 configurations equals the static scan", maxDifference(single[0].values, averaged[0].values) < 1e-12)
        r.check("an Al atom: the annular total is ≤ 1 and above 0.9 at every position",
                single[0].values.allSatisfy { $0 <= 1 + 1e-9 && $0 > 0.9 }, "got \(single[0].values)")
        let adf = try STEM.scan(probe: probe, potentials: [staticPotential], scan: scan, detectors: [AnnularDetector(inner: 20, outer: nil)])
        r.check("an Al atom scatters some intensity beyond 20 mrad; far from the vacuum it is zero up to the aperture",
                adf[0].values.allSatisfy { $0 > 0 && $0 < 0.5 })
        let shifted = try STEM.scan(probe: probe, potentials: [staticPotential], scan: scan, detectors: [])
        r.equal("no detectors → no measurements", shifted.count, 0)

        if Multislice.metalIsAvailable {
            let gpu = try STEM.scan(probe: probe, potentials: [staticPotential], scan: scan, detectors: [AnnularDetector()], backend: .metal)
            r.check("STEM.scan(.metal) agrees with .cpu (1e-4)", maxDifference(gpu[0].values, single[0].values) < 1e-4)
        } else {
            r.skip("STEM.scan(.metal)", "Multislice.metalIsAvailable is false (no Metal device on this machine)")
        }

        // An empty literal needs its type since 0.1.0 (two overloads: [ProjectedPotential] and [Potential]).
        r.refusesCase("STEM.scan with no potential", "invalidScan") {
            try STEM.scan(probe: probe, potentials: [] as [ProjectedPotential], scan: scan, detectors: [AnnularDetector()])
        }
        r.refusesCase("STEM.scan with no Potential", "invalidScan") {
            try STEM.scan(probe: probe, potentials: [] as [Potential], scan: scan, detectors: [AnnularDetector()])
        }
        let other = try ProjectedPotential(structure: AtomicStructure(box: SIMD3(6, 6, 4), atoms: []), nx: 32, ny: 32, sliceThickness: 2)
        r.refuses("STEM.scan with a potential on another grid", .gridMismatch) {
            try STEM.scan(probe: probe, potentials: [other], scan: scan, detectors: [AnnularDetector()])
        }
        r.refusesCase("STEM.scan whose detector refuses the pattern", "invalidDetector") {
            try STEM.scan(probe: probe, potentials: [vacuum], scan: scan, detectors: [AnnularDetector(inner: 0, outer: 1000)])
        }
        r.check("SimulatoDetector.plannedScope names the detectors", SimulatoDetector.plannedScope.contains("annular"))
    }
}

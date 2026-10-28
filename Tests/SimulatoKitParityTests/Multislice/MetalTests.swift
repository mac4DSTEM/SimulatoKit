//
//  MetalTests.swift
//  Role: the Metal path (float32, MPSGraph) against the CPU path (float64) on the
//        same potential, for every pinned case and beam and one 512 × 512, 100-slice crystal. Disabled with a named
//        reason when no Metal device exists.
//

import SimulatoKitTestSupport
import Testing
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe

@Suite("Metal against the CPU path",
       .enabled(if: Multislice.metalIsAvailable,
                "no Metal device (MTLCreateSystemDefaultDevice() returned nil): the Metal path is not tested here"))
struct MetalTests {
    /// The total-intensity bound per slice (MPSGraph's float32 FFT round trip loses up to 1.43e−6
    /// of Σ|ψ|² per round trip on the non-zero frequencies; tools/metal-fft-gain).
    static let totalIntensityPerSlice = MetalAccuracy.registered.totalIntensityPerSlice

    struct Worst {
        var wave = 0.0, intensity = 0.0, pattern = 0.0, total = 0.0, totalPerSlice = 0.0
    }

    func compare(_ incident: Waves, _ potential: ProjectedPotential, _ planes: ExitPlanes, label: String) throws -> Worst {
        var worst = Worst()
        let planeIndices = try planes.planes(sliceCount: potential.slices.count)
        let cpu = try Multislice.run(incident, through: potential, exitPlanes: planes, backend: .cpu)
        let gpu = try Multislice.run(incident, through: potential, exitPlanes: planes, backend: .metal)
        #expect(cpu.count == gpu.count)
        for ((a, b), plane) in zip(zip(gpu, cpu), planeIndices) {
            let w = relMax(a.field, b.field)
            let i = relMax(a.intensity().values, b.intensity().values)
            let p = relMax(try a.diffractionPattern().values, try b.diffractionPattern().values)
            let t = relDiff(a.totalIntensity, b.totalIntensity)
            worst.wave = max(worst.wave, w)
            worst.intensity = max(worst.intensity, i)
            worst.pattern = max(worst.pattern, p)
            worst.total = max(worst.total, t)
            if plane >= 0 { worst.totalPerSlice = max(worst.totalPerSlice, t / Double(plane + 1)) }
            #expect(w <= MetalAccuracy.registered.wave, "\(label): exit wave")
            #expect(i <= MetalAccuracy.registered.intensity, "\(label): intensity")
            #expect(p <= MetalAccuracy.registered.pattern, "\(label): diffraction pattern")
            #expect(t <= Self.totalIntensityPerSlice * Double(plane + 1), "\(label): total intensity")
        }
        return worst
    }

    @Test(arguments: CaseBeam.all)
    func metalMatchesTheCPUPathOnThePinnedCases(r: CaseBeam) throws {
        let c = Pins.case(r.id)
        let potential = try Pins.potential(r.potentialCase, pi: .exact)
        let incident = try Pins.incident(r.id, r.label, grid: potential.grid)
        let worst = try compare(incident, potential, .afterSlices(c["exit_planes"] as! [Int]), label: "\(r.id)/\(r.label)")
        print("Metal vs CPU, \(r.id)/\(r.label): wave \(worst.wave), intensity \(worst.intensity), pattern \(worst.pattern), total \(worst.total), per slice \(worst.totalPerSlice)")
    }

    @Test func metalMatchesTheCPUPathOn512Squared100Slices() throws {
        let cell = try AtomicStructure(box: SIMD3(4.05, 4.05, 4.05), atoms: [
            .init(atomicNumber: 13, position: SIMD3(0, 0, 0)), .init(atomicNumber: 13, position: SIMD3(0, 2.025, 2.025)),
            .init(atomicNumber: 13, position: SIMD3(2.025, 0, 2.025)), .init(atomicNumber: 13, position: SIMD3(2.025, 2.025, 0)),
        ])
        let potential = try ProjectedPotential(structure: cell.repeated(8, 8, 50), nx: 512, ny: 512, sliceThickness: 2.025)
        #expect(potential.slices.count == 100)
        let beam = try ElectronBeam(energy: 300e3)
        let plane = try compare(.planeWave(beam: beam, grid: potential.grid), potential, .every(25), label: "Al 512² plane")
        var ab = Aberrations()
        ab.defocus = 30
        let probe = try Probe(beam: beam, grid: potential.grid, aperture: Aperture(semiangleCutoff: 25), aberrations: ab)
        let focused = try compare(probe.build(at: SIMD2(15.3, 16.9)), potential, .every(25), label: "Al 512² probe")
        for (name, w) in [("plane", plane), ("probe", focused)] {
            print("Metal vs CPU, 512² x 100 slices, \(name): wave \(w.wave), intensity \(w.intensity), pattern \(w.pattern), total \(w.total), per slice \(w.totalPerSlice)")
        }
    }
}

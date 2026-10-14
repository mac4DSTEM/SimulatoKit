//
//  MetalTests.swift
//  Role: registration section 4, second table: the Metal path (float32, MPSGraph) against the CPU path (float64) on the
//        same potential, for every pinned case and beam and one 512 × 512, 100-slice crystal. Skipped with a named
//        reason when no Metal device exists.
//

import XCTest
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe

final class MetalTests: XCTestCase {
    override func setUpWithError() throws {
        guard Multislice.metalIsAvailable else {
            throw XCTSkip("no Metal device (MTLCreateSystemDefaultDevice() returned nil): the Metal path is not tested here")
        }
    }

    struct Worst {
        var wave = 0.0, intensity = 0.0, pattern = 0.0, total = 0.0, totalPerSlice = 0.0
    }

    /// The total-intensity bound per slice (amendment 2, Gate D: MPSGraph's float32 FFT round trip loses up to 1.43e−6
    /// of Σ|ψ|² per round trip on the non-zero frequencies; tools/metal-fft-gain).
    static let totalIntensityPerSlice = 1.5e-6

    func compare(_ incident: Waves, _ potential: ProjectedPotential, _ planes: ExitPlanes, label: String,
                 worst: inout Worst) throws {
        let planeIndices = try planes.planes(sliceCount: potential.slices.count)
        let cpu = try Multislice.run(incident, through: potential, exitPlanes: planes, backend: .cpu)
        let gpu = try Multislice.run(incident, through: potential, exitPlanes: planes, backend: .metal)
        XCTAssertEqual(cpu.count, gpu.count)
        for ((a, b), plane) in zip(zip(gpu, cpu), planeIndices) {
            let w = relMax(a.field, b.field)
            let i = relMax(a.intensity().values, b.intensity().values)
            let p = relMax(try a.diffractionPattern().values, try b.diffractionPattern().values)
            let t = relDiff(a.totalIntensity, b.totalIntensity)
            worst.wave = max(worst.wave, w)
            worst.intensity = max(worst.intensity, i)
            worst.pattern = max(worst.pattern, p)
            worst.total = max(worst.total, t)
            XCTAssertLessThanOrEqual(w, 1e-4, "\(label): exit wave")
            XCTAssertLessThanOrEqual(i, 2e-4, "\(label): intensity")
            XCTAssertLessThanOrEqual(p, 2e-4, "\(label): diffraction pattern")
            XCTAssertLessThanOrEqual(t, Self.totalIntensityPerSlice * Double(plane + 1), "\(label): total intensity")
            if plane >= 0 { worst.totalPerSlice = max(worst.totalPerSlice, t / Double(plane + 1)) }
        }
    }

    func testMetalMatchesTheCPUPathOnThePinnedCases() throws {
        var worst = Worst()
        for (id, label) in MultislicePinsTests.runs {
            let c = Pins.case(id)
            let potential = try Pins.potential(MultislicePinsTests.potentialCase(id), pi: .exact)
            let incident = try Pins.incident(id, label, grid: potential.grid)
            try compare(incident, potential, .afterSlices(c["exit_planes"] as! [Int]), label: "\(id)/\(label)",
                        worst: &worst)
        }
        print("Metal vs CPU, pinned cases: wave \(worst.wave), intensity \(worst.intensity), pattern \(worst.pattern), total \(worst.total), per slice \(worst.totalPerSlice)")
    }

    func testMetalMatchesTheCPUPathOn512Squared100Slices() throws {
        let cell = try AtomicStructure(box: SIMD3(4.05, 4.05, 4.05), atoms: [
            .init(atomicNumber: 13, position: SIMD3(0, 0, 0)), .init(atomicNumber: 13, position: SIMD3(0, 2.025, 2.025)),
            .init(atomicNumber: 13, position: SIMD3(2.025, 0, 2.025)), .init(atomicNumber: 13, position: SIMD3(2.025, 2.025, 0)),
        ])
        let crystal = cell.repeated(8, 8, 50)
        let potential = try ProjectedPotential(structure: crystal, nx: 512, ny: 512, sliceThickness: 2.025)
        XCTAssertEqual(potential.slices.count, 100)
        let beam = try ElectronBeam(energy: 300e3)
        var worst = Worst()
        try compare(.planeWave(beam: beam, grid: potential.grid), potential, .every(25), label: "Al 512² plane",
                    worst: &worst)
        var ab = Aberrations()
        ab.defocus = 30
        let probe = try Probe(beam: beam, grid: potential.grid, aperture: Aperture(semiangleCutoff: 25), aberrations: ab)
        try compare(probe.build(at: SIMD2(15.3, 16.9)), potential, .every(25), label: "Al 512² probe", worst: &worst)
        print("Metal vs CPU, 512² x 100 slices: wave \(worst.wave), intensity \(worst.intensity), pattern \(worst.pattern), total \(worst.total), per slice \(worst.totalPerSlice)")
    }
}

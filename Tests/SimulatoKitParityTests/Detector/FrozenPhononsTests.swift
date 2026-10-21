//
//  FrozenPhononsTests.swift
//  Role: registration-slice-2.md: frozen phonons against abTEM with its displacements given (F1: SrTiO₃, four
//        configurations, the averaged plane-wave pattern), and truth 1 (the kit's generator by its statistics) and
//        truth 2 (σ = 0 is slice 1).
//

import SimulatoKitTestSupport
import Foundation
import Testing
import SimulatoCore
import SimulatoMultislice

@Suite("Frozen phonons")
struct FrozenPhononsTests {
    @Test(arguments: [LobatoParameterization.PiConstant.abTEMFloat32, .exact])
    func averagedPlaneWavePatternEqualsAbTEM(pi: LobatoParameterization.PiConstant) throws {
        let potentials = try Pins2.potentials("F1", pi: pi)
        #expect(potentials.count == 4)
        let bound = pi == .abTEMFloat32 ? 1e-11 : 1e-7
        #expect(relMax2(potentials[0].slices[0].values, Pins2.doubles("F1.V.config0.slice0")) <= bound, "V, config 0")
        let beam = try ElectronBeam(energy: Pins2.case("F1")["energy"] as! Double)
        var mean: [Double] = []
        let totals = Pins2.case("F1")["per_config_total"] as! [Double]
        for (c, p) in potentials.enumerated() {
            let exit = try Multislice.run(.planeWave(beam: beam, grid: p.grid), through: p)[0]
            #expect(abs(exit.totalIntensity - totals[c]) / totals[c] <= (pi == .abTEMFloat32 ? 1e-10 : 1e-6), "Σ|ψ|², config \(c)")
            let pattern = try exit.diffractionPattern().values
            mean = mean.isEmpty ? pattern : zip(mean, pattern).map(+)
        }
        mean = mean.map { $0 / 4 }
        #expect(relMax2(mean, Pins2.doubles("F1.pattern_mean")) <= (pi == .abTEMFloat32 ? 1e-10 : 1e-6), "averaged pattern")
    }

    /// Truth 1: the kit's generator, by its statistics.
    @Test func generatorStatistics() throws {
        var g = GaussianGenerator(seed: 1)
        let n = 1_000_000
        var s1 = 0.0, s2 = 0.0, s4 = 0.0, lag = 0.0, previous = 0.0
        for i in 0..<n {
            let z = g.nextNormal()
            s1 += z; s2 += z * z; s4 += z * z * z * z
            if i > 0 { lag += z * previous }
            previous = z
        }
        let mean = s1 / Double(n), variance = s2 / Double(n) - mean * mean
        #expect(abs(mean) <= 5e-3, "mean")
        #expect(abs(variance - 1) <= 7e-3, "variance")
        #expect(abs(s4 / Double(n) - 3) <= 0.05, "fourth moment")
        #expect(abs(lag / Double(n - 1)) <= 5e-3, "lag-1 correlation")
        var a = GaussianGenerator(seed: 1), b = GaussianGenerator(seed: 1), c = GaussianGenerator(seed: 2)
        let sa = (0..<16).map { _ in a.nextNormal() }, sb = (0..<16).map { _ in b.nextNormal() }
        let sc = (0..<16).map { _ in c.nextNormal() }
        #expect(sa == sb, "the same seed, the same stream")
        #expect(sa != sc, "another seed, another stream")
        // Consecutive seeds share no configuration (the refuter found seed + c did).
        let one = try AtomicStructure(box: SIMD3(5, 5, 5), atoms: [.init(atomicNumber: 6, position: SIMD3(1, 1, 1))])
        let e1 = try FrozenPhonons(structure: one, sigmas: [6: 0.1], configurations: 3, seed: 10)
        let e2 = try FrozenPhonons(structure: one, sigmas: [6: 0.1], configurations: 3, seed: 11)
        #expect(Set(e1.displacements.map { $0[0].x }).isDisjoint(with: Set(e2.displacements.map { $0[0].x })))
    }

    /// Truth 1: σ per axis and the directions.
    @Test func displacementsHaveTheirSigmaAndRespectTheDirections() throws {
        let atoms = (0..<10_000).map { AtomicStructure.Atom(atomicNumber: 13, position: SIMD3(Double($0 % 100) * 0.1, 1, 1)) }
        let s = try AtomicStructure(box: SIMD3(10, 10, 10), atoms: atoms)
        let fp = try FrozenPhonons(structure: s, sigmas: [13: 0.1], configurations: 2, seed: 42, directions: [.x, .y])
        for d in fp.displacements {
            for axis in 0..<2 {
                let v = d.map { $0[axis] }
                let m = v.reduce(0, +) / Double(v.count)
                let sd = (v.map { ($0 - m) * ($0 - m) }.reduce(0, +) / Double(v.count)).squareRoot()
                #expect(abs(sd - 0.1) <= 0.003, "σ along axis \(axis)")
            }
            #expect(d.allSatisfy { $0.z == 0 }, "z is not displaced")
        }
        #expect(fp.displacements[0] != fp.displacements[1], "configurations differ")
        let again = try FrozenPhonons(structure: s, sigmas: [13: 0.1], configurations: 2, seed: 42, directions: [.x, .y])
        #expect(again.displacements == fp.displacements, "seeded: reproducible")
        let zOnly = try FrozenPhonons(structure: s, sigmas: [13: 0.1], configurations: 1, seed: 42, directions: [.z])
        #expect(zOnly.displacements[0].allSatisfy { $0.x == 0 && $0.y == 0 }, "only z is displaced")
        #expect(zOnly.displacements[0].map(\.z) == (try FrozenPhonons(structure: s, sigmas: [13: 0.1], configurations: 1,
                                                                       seed: 42)).displacements[0].map(\.z),
                "the streams do not depend on the directions")
    }

    /// Truth 2: σ = 0 is slice 1, to the last bit; the mean of identical configurations is the single result.
    @Test func zeroSigmaIsSliceOne() throws {
        let s = try Pins2.structure("S2")
        let fp = try FrozenPhonons(structure: s, sigmas: [13: 0], configurations: 3, seed: 9)
        let a = try fp.potentials(nx: 128, ny: 128, sliceThickness: 2.025)
        let b = try ProjectedPotential(structure: s, nx: 128, ny: 128, sliceThickness: 2.025)
        for p in a { #expect(p.slices == b.slices, "every configuration equals the static potential") }
    }

    @Test func refusals() throws {
        let s = try AtomicStructure(box: SIMD3(5, 5, 5), atoms: [.init(atomicNumber: 6, position: SIMD3(1, 1, 1))])
        #expect(throws: (any Error).self) { try FrozenPhonons(structure: s, displacements: []) }
        #expect(throws: (any Error).self) { try FrozenPhonons(structure: s, displacements: [[]]) }
        #expect(throws: (any Error).self) { try FrozenPhonons(structure: s, displacements: [[SIMD3(.nan, 0, 0)]]) }
        #expect(throws: (any Error).self) { try FrozenPhonons(structure: s, sigmas: [6: -0.1], configurations: 1, seed: 0) }
        #expect(throws: (any Error).self) { try FrozenPhonons(structure: s, sigmas: [6: 0.1], configurations: 0, seed: 0) }
        #expect(throws: (any Error).self) { try FrozenPhonons(structure: s, sigmas: [13: 0.1], configurations: 1, seed: 0) }
    }
}

//
//  BeamAndScatteringTests.swift
//  Role: registration sections 3.1–3.3 against the pins (λ, σ, κ, the float32 π, f(k²) under both π choices) and
//        truth 1 (Kirkland's tabulated λ and σ).
//

import SimulatoKitTestSupport
import Testing
import SimulatoCore
import SimulatoMultislice

@Suite("Beam and scattering factor")
struct BeamAndScatteringTests {
    @Test func kappaAndFloat32PiEqualThePins() {
        #expect(PhysicalConstants.kappa == Pins.json["kappa"] as! Double)
        #expect(LobatoParameterization.float32Pi == Pins.json["pi_float32"] as! Double)
        #expect(LobatoParameterization.float32PiSquared == Pins.json["pi2_float32"] as! Double)
    }

    @Test func wavelengthAndSigmaEqualAbTEM() throws {
        for row in Pins.json["beam"] as! [[String: Any]] {
            let beam = try ElectronBeam(energy: row["energy"] as! Double)
            #expect(relDiff(beam.wavelength, row["wavelength"] as! Double) <= 1e-15, "λ at \(beam.energy)")
            #expect(relDiff(beam.interactionParameter, row["sigma"] as! Double) <= 1e-15, "σ at \(beam.energy)")
        }
    }

    /// Truth 1: Kirkland (2010) Table 2.2, four significant figures.
    @Test(arguments: [(100e3, 0.03701, 9.244e-4), (200e3, 0.02508, 7.288e-4), (300e3, 0.01969, 6.526e-4)])
    func wavelengthAndSigmaMatchKirklandsTable(energy: Double, lambda: Double, sigma: Double) throws {
        let beam = try ElectronBeam(energy: energy)
        #expect(abs(beam.wavelength - lambda) <= 0.5e-5, "λ at \(energy)")
        #expect(abs(beam.interactionParameter - sigma) <= 0.5e-7, "σ at \(energy)")
    }

    @Test(arguments: [("H", 1), ("C", 6), ("O", 8), ("Al", 13), ("Si", 14), ("Ti", 22), ("Sr", 38), ("Au", 79), ("U", 92)])
    func projectedScatteringFactorEqualsAbTEMUnderBothPiChoices(symbol: String, z: Int) throws {
        let k2 = Pins.json["f_k2"] as! [Double]
        let f = (Pins.json["f"] as! [String: [String: [Double]]])[symbol]!
        let abtem = try LobatoParameterization(piConstant: .abTEMFloat32).projectedScatteringFactor(atomicNumber: z, k2: k2)
        let exact = try LobatoParameterization(piConstant: .exact).projectedScatteringFactor(atomicNumber: z, k2: k2)
        for i in k2.indices {
            #expect(relDiff(abtem[i], f["abtem"]![i]) <= 1e-14, "\(symbol) k2 \(k2[i])")
            #expect(relDiff(exact[i], f["exact"]![i]) <= 1e-14, "\(symbol) exact k2 \(k2[i])")
            #expect(relDiff(exact[i], f["abtem"]![i]) <= 1e-7, "\(symbol) exact vs abTEM")
        }
    }

    @Test func refusals() {
        #expect(throws: (any Error).self) { try ElectronBeam(energy: 0) }
        #expect(throws: (any Error).self) { try ElectronBeam(energy: .nan) }
        #expect(throws: (any Error).self) { try LobatoParameterization.parameters(atomicNumber: 104) }
        #expect(throws: (any Error).self) { try LobatoParameterization.parameters(atomicNumber: 0) }
        #expect(throws: (any Error).self) { try Grid(nx: 1, ny: 128, extentX: 10, extentY: 10) }
        #expect(throws: Never.self) { try Grid(nx: 100, ny: 127, extentX: 10, extentY: 10) }
        #expect(throws: (any Error).self) { try Grid(nx: 128, ny: 128, extentX: 0, extentY: 10) }
        #expect(throws: Never.self) { try Grid(nx: 120, ny: 160, extentX: 10, extentY: 10) }
    }
}

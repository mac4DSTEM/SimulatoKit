//
//  BeamAndScatteringTests.swift
//  Role: registration sections 3.1–3.3 against the pins (λ, σ, κ, the float32 π, f(k²) under both π choices) and
//        truth 1 (Kirkland's tabulated λ and σ).
//

import XCTest
import SimulatoCore
import SimulatoMultislice

final class BeamAndScatteringTests: XCTestCase {
    func testKappaAndFloat32PiEqualThePins() {
        XCTAssertEqual(PhysicalConstants.kappa, Pins.json["kappa"] as! Double)
        XCTAssertEqual(LobatoParameterization.float32Pi, Pins.json["pi_float32"] as! Double)
        XCTAssertEqual(LobatoParameterization.float32PiSquared, Pins.json["pi2_float32"] as! Double)
    }

    func testWavelengthAndSigmaEqualAbTEM() throws {
        for row in Pins.json["beam"] as! [[String: Any]] {
            let beam = try ElectronBeam(energy: row["energy"] as! Double)
            XCTAssertLessThanOrEqual(relDiff(beam.wavelength, row["wavelength"] as! Double), 1e-15, "λ at \(beam.energy)")
            XCTAssertLessThanOrEqual(relDiff(beam.interactionParameter, row["sigma"] as! Double), 1e-15,
                                     "σ at \(beam.energy)")
        }
    }

    /// Truth 1: Kirkland (2010) Table 2.2, four significant figures.
    func testWavelengthAndSigmaMatchKirklandsTable() throws {
        let table: [(Double, Double, Double)] = [(100e3, 0.03701, 9.244e-4), (200e3, 0.02508, 7.288e-4),
                                                 (300e3, 0.01969, 6.526e-4)]
        for (e, lambda, sigma) in table {
            let beam = try ElectronBeam(energy: e)
            XCTAssertEqual(beam.wavelength, lambda, accuracy: 0.5e-5, "λ at \(e)")
            XCTAssertEqual(beam.interactionParameter, sigma, accuracy: 0.5e-7, "σ at \(e)")
        }
    }

    func testProjectedScatteringFactorEqualsAbTEMUnderBothPiChoices() throws {
        let k2 = Pins.json["f_k2"] as! [Double]
        let symbols = ["H": 1, "C": 6, "O": 8, "Al": 13, "Si": 14, "Ti": 22, "Sr": 38, "Au": 79, "U": 92]
        let f = Pins.json["f"] as! [String: [String: [Double]]]
        for (symbol, z) in symbols {
            let abtem = try LobatoParameterization(piConstant: .abTEMFloat32).projectedScatteringFactor(atomicNumber: z, k2: k2)
            let exact = try LobatoParameterization(piConstant: .exact).projectedScatteringFactor(atomicNumber: z, k2: k2)
            for i in k2.indices {
                XCTAssertLessThanOrEqual(relDiff(abtem[i], f[symbol]!["abtem"]![i]), 1e-14, "\(symbol) k2 \(k2[i])")
                XCTAssertLessThanOrEqual(relDiff(exact[i], f[symbol]!["exact"]![i]), 1e-14, "\(symbol) exact k2 \(k2[i])")
                XCTAssertLessThanOrEqual(relDiff(exact[i], f[symbol]!["abtem"]![i]), 1e-7, "\(symbol) exact vs abTEM")
            }
        }
    }

    func testRefusals() {
        XCTAssertThrowsError(try ElectronBeam(energy: 0))
        XCTAssertThrowsError(try ElectronBeam(energy: .nan))
        XCTAssertThrowsError(try LobatoParameterization.parameters(atomicNumber: 104))
        XCTAssertThrowsError(try LobatoParameterization.parameters(atomicNumber: 0))
        XCTAssertThrowsError(try Grid(nx: 100, ny: 128, extentX: 10, extentY: 10))
        XCTAssertThrowsError(try Grid(nx: 128, ny: 128, extentX: 0, extentY: 10))
        XCTAssertNoThrow(try Grid(nx: 120, ny: 160, extentX: 10, extentY: 10))
    }
}

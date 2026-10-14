//
//  MinimalCIFTests.swift
//  Role: the example's CIF reader (tier 2): the shipped Al CIF expands to the four fcc atoms, operations with
//        fractions and signs parse, a non-orthogonal cell and a truncated loop are refused.
//

import Foundation
import Testing
import SimulatoCore
@testable import SimulatoExample

@Suite("The example CIF reader")
struct MinimalCIFTests {
    @Test func aluminiumExpandsToFourAtoms() throws {
        let url = URL(fileURLWithPath: #filePath).deletingLastPathComponent().deletingLastPathComponent()
            .deletingLastPathComponent().appendingPathComponent("Examples/SimulatoExample/Al.cif")
        let (name, cell) = try MinimalCIF.read(String(contentsOf: url, encoding: .utf8))
        #expect(name == "Al")
        #expect(cell.box == SIMD3(4.0495, 4.0495, 4.0495))
        #expect(cell.atoms.count == 4)
        #expect(cell.atoms.allSatisfy { $0.atomicNumber == 13 })
        let expected: [SIMD3<Double>] = [SIMD3(0, 0, 0), SIMD3(0, 0.5, 0.5), SIMD3(0.5, 0, 0.5), SIMD3(0.5, 0.5, 0)]
        for (atom, p) in zip(cell.atoms, expected) {
            #expect(abs((atom.position.x) - (p.x * 4.0495)) <= 1e-12)
            #expect(abs((atom.position.y) - (p.y * 4.0495)) <= 1e-12)
            #expect(abs((atom.position.z) - (p.z * 4.0495)) <= 1e-12)
        }
    }

    @Test func operationsParse() throws {
        let p = SIMD3(0.1, 0.2, 0.3)
        #expect(try MinimalCIF.apply("-y,x-y,z+1/3", p) == SIMD3(-0.2, 0.1 - 0.2, 0.3 + 1.0 / 3.0))
        #expect(try MinimalCIF.apply("1/2-x, 1/2+y, -z", p) == SIMD3(0.5 - 0.1, 0.5 + 0.2, -0.3))
        #expect(throws: (any Error).self) { try MinimalCIF.apply("x,y", p) }
    }

    @Test func refusals() {
        let monoclinic = """
        data_m
        _cell_length_a 5
        _cell_length_b 6
        _cell_length_c 7
        _cell_angle_alpha 90
        _cell_angle_beta 105.3
        _cell_angle_gamma 90
        loop_
        _atom_site_label
        _atom_site_fract_x
        _atom_site_fract_y
        _atom_site_fract_z
        Mg1 0 0 0
        """
        #expect(throws: (any Error).self) { try MinimalCIF.read(monoclinic) }
        let truncated = monoclinic.replacingOccurrences(of: "105.3", with: "90").replacingOccurrences(of: "Mg1 0 0 0", with: "Mg1 0 0")
        #expect(throws: (any Error).self) { try MinimalCIF.read(truncated) }
    }
}

//
//  CrystalSliceThicknessTests.swift
//  Role: `Slicing.crystalSliceThicknesses` against abTEM 1.0.10 `crystal_slice_thicknesses`, the numbers printed by
//        tools/pins/crystal-slices/pins.py held inline at rtol 1e-7 (assert_allclose default; abTEM states none).
//

import Testing
import SimulatoCore
import SimulatoMultislice
import SimulatoKitTestSupport

@Suite("Crystal slice thicknesses")
struct CrystalSliceThicknessTests {
    static func planted(_ z: [Double], depth: Double) throws -> AtomicStructure {
        try AtomicStructure(box: SIMD3(4, 4, depth),
                            atoms: z.map { AtomicStructure.Atom(atomicNumber: 6, position: SIMD3(1, 2, $0)) })
    }

    static func expectClose(_ a: [Double], _ b: [Double]) {
        #expect(a.count == b.count, "slice count")
        guard a.count == b.count else { return }
        for (x, y) in zip(a, b) { #expect(abs(x - y) <= 1e-7 * abs(y), "\(x) against \(y)") }
    }

    @Test func plantedStructuresMatchAbTEM() throws {
        let four = try Self.planted([0.0, 0.6, 1.0, 2.05, 3.99], depth: 4)
        Self.expectClose(try Slicing.crystalSliceThicknesses(four),
                         [0.4, 0.6000000000000001, 1.0, 1.8000000000000003, 0.20000000000000018])
        Self.expectClose(try Slicing.crystalSliceThicknesses(four, tolerance: 0.5), [0.5, 0.5, 1.0, 1.5, 0.5])
        // 0.3 / 0.2 and 0.6 / 0.2 fall just below 1.5 and 3 in float64: bins 1 and 2, as NumPy's floor gives.
        let six = try Self.planted([0.3, 0.6, 1.8, 3.0, 4.2, 5.4], depth: 6)
        Self.expectClose(try Slicing.crystalSliceThicknesses(six, tolerance: 0.2),
                         [0.20000000000000004, 0.19999999999999996, 1.4000000000000001, 1.2, 1.1999999999999997,
                          1.2000000000000002, 0.6000000000000005])
    }

    @Test func aluminiumAtACommensurateToleranceMatchesAbTEM() throws {
        let al = try aluminiumStructure(cells: (1, 1, 1))
        Self.expectClose(try Slicing.crystalSliceThicknesses(al, tolerance: 0.45), [1.7999999999999998, 2.2500000000000004])
        let potential = try Potential(structure: al, nx: 16, ny: 16,
                                      sliceThicknesses: try Slicing.crystalSliceThicknesses(al, tolerance: 0.45))
        #expect(potential.sliceCount == 2)
    }

    @Test func refusesWhatAbTEMRefuses() throws {
        // abTEM asserts that the thicknesses sum to the depth: an incommensurate tolerance, an atom below 0, a box
        // thinner than the tolerance.
        let six = try Self.planted([0.3, 0.6, 1.8, 3.0, 4.2, 5.4], depth: 6)
        #expect(throws: SimulatoError.self) { try Slicing.crystalSliceThicknesses(six, tolerance: 0.7) }
        #expect(throws: SimulatoError.self) { try Slicing.crystalSliceThicknesses(try aluminiumStructure(cells: (1, 1, 1))) }
        #expect(throws: SimulatoError.self) {
            try Slicing.crystalSliceThicknesses(try Self.planted([-0.3, 1.0], depth: 4))
        }
        #expect(throws: SimulatoError.self) { try Slicing.crystalSliceThicknesses(try Self.planted([0.05], depth: 0.1)) }
    }

    @Test func refusesABadTolerance() throws {
        let four = try Self.planted([0.0, 1.0], depth: 4)
        for t in [0, -0.2, .nan, .infinity, 1e-300] {
            #expect(throws: SimulatoError.self) { try Slicing.crystalSliceThicknesses(four, tolerance: t) }
        }
    }
}

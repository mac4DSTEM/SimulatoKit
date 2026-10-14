//
//  PotentialPinsTests.swift
//  Role: registration sections 3.4–3.6 against the pins: the slicing (thicknesses, the slice of every atom), the
//        projected potential of every slice (stored arrays and per-slice sums and maxima), the bandlimited transmission
//        function, and truth 2 (the k = 0 identity) on every slice of every case.
//

import XCTest
import SimulatoCore
import SimulatoMultislice

final class PotentialPinsTests: XCTestCase {
    static let cases = ["A1", "A2", "B1", "C1", "D1"]

    func testSlicingEqualsAbTEM() throws {
        for id in Self.cases {
            let c = Pins.case(id)
            let potential = try Pins.potential(id, pi: .abTEMFloat32)
            XCTAssertEqual(potential.slices.map(\.thickness), c["slice_thickness"] as! [Double], "\(id) thicknesses")
            XCTAssertEqual(potential.sliceOfAtom, c["slice_of_atom"] as! [Int], "\(id) slice of each atom")
        }
    }

    func testPotentialEqualsAbTEMWithItsPi() throws {
        for id in Self.cases {
            let c = Pins.case(id)
            let potential = try Pins.potential(id, pi: .abTEMFloat32)
            let maxima = c["V_max_abs"] as! [Double], sums = c["V_sum"] as! [Double]
            let globalMax = maxima.max()!
            for (j, slice) in potential.slices.enumerated() {
                let m = slice.values.map(abs).max()!
                XCTAssertLessThanOrEqual(abs(m - maxima[j]), 1e-11 * max(globalMax, 1e-300), "\(id) slice \(j) max |V|")
                let s = slice.values.reduce(0, +)
                XCTAssertLessThanOrEqual(abs(s - sums[j]), 1e-11 * max(abs(sums[j]), 1) , "\(id) slice \(j) ΣV")
                if Pins.hasArray("\(id).V.\(j)") {
                    XCTAssertLessThanOrEqual(relMax(slice.values, Pins.doubles("\(id).V.\(j)")), 1e-11, "\(id) V slice \(j)")
                }
            }
        }
    }

    func testPotentialWithTheExactPiIsWithinTheRegisteredDifference() throws {
        for id in Self.cases {
            let potential = try Pins.potential(id, pi: .exact)
            for j in 0..<2 where Pins.hasArray("\(id).V.\(j)") {
                let d = relMax(potential.slices[j].values, Pins.doubles("\(id).V.\(j)"))
                XCTAssertLessThanOrEqual(d, 1e-7, "\(id) exact-π V slice \(j)")
                XCTAssertGreaterThan(d, 0, "\(id): the exact π must differ from abTEM's float32 π somewhere")
            }
        }
    }

    func testBandlimitedTransmissionEqualsAbTEM() throws {
        for id in Self.cases where Pins.hasArray("\(id).t.0") {
            let potential = try Pins.potential(id, pi: .abTEMFloat32)
            let t = try Multislice.transmissionFunctions(potential, beam: Pins.beam(id))
            XCTAssertLessThanOrEqual(relMax(t[0], Pins.complex("\(id).t.0")), 1e-11, "\(id) t slice 0")
        }
    }

    /// Truth 2: Σ V ΔxΔy = Σ_atoms f(0) for every slice (S(0) = ΔxΔy and the bilinear weights sum to one).
    func testIntegratedPotentialEqualsTheSumOfFAtZero() throws {
        for id in Self.cases {
            for pi in [LobatoParameterization.PiConstant.exact, .abTEMFloat32] {
                let potential = try Pins.potential(id, pi: pi)
                let structure = try Pins.structure(id)
                let area = potential.grid.samplingX * potential.grid.samplingY
                for (j, slice) in potential.slices.enumerated() {
                    let zs = zip(structure.atoms, potential.sliceOfAtom).filter { $0.1 == j }.map(\.0.atomicNumber)
                    let expected = try zs.reduce(0.0) {
                        $0 + (try LobatoParameterization(piConstant: pi).projectedScatteringFactor(atomicNumber: $1, k2: 0))
                    }
                    let integral = slice.values.reduce(0, +) * area
                    if zs.isEmpty {
                        XCTAssertEqual(slice.values.map(abs).max()!, 0, "\(id) empty slice \(j) is zero exactly")
                    } else {
                        XCTAssertLessThanOrEqual(relDiff(integral, expected), 1e-12, "\(id) slice \(j) k = 0 identity")
                    }
                }
            }
        }
    }

    /// Section 3.4, the top-face rule: an atom within 1e−10 Å of the top face is put at z = 0 (slice 0), so its potential
    /// equals that of the same atom at z = 0.
    func testTopFaceAtomGoesToTheFirstSlice() throws {
        let top = try AtomicStructure(box: SIMD3(8, 8, 4), atoms: [.init(atomicNumber: 29, position: SIMD3(2.3, 5.1, 4 - 1e-11))])
        let bottom = try AtomicStructure(box: SIMD3(8, 8, 4), atoms: [.init(atomicNumber: 29, position: SIMD3(2.3, 5.1, 0))])
        let a = try ProjectedPotential(structure: top, nx: 64, ny: 64, sliceThickness: 1)
        let b = try ProjectedPotential(structure: bottom, nx: 64, ny: 64, sliceThickness: 1)
        XCTAssertEqual(a.sliceOfAtom, [0])
        XCTAssertEqual(a.slices[0].values, b.slices[0].values)
        XCTAssertEqual(a.slices[3].values.map(abs).max(), 0)
    }

    /// Section 3.4, the nudged edges: 30 slices of 0.1 Å accumulate to 0.30000000000000004 after three, so without the
    /// −1e−12 an atom at z = 0.3 would fall into slice 2; abTEM puts it into slice 3 (the next slice).
    func testAtomOnADriftedSliceBoundaryGoesToTheNextSlice() throws {
        XCTAssertGreaterThan(0.1 + 0.1 + 0.1, 0.3, "the premise: the cumulative sum overshoots")
        let s = try AtomicStructure(box: SIMD3(8, 8, 3), atoms: [.init(atomicNumber: 6, position: SIMD3(1, 1, 0.3))])
        let p = try ProjectedPotential(structure: s, nx: 64, ny: 64, sliceThickness: 0.1)
        XCTAssertEqual(p.slices.count, 30)
        XCTAssertEqual(p.sliceOfAtom, [3])
    }

    func testRefusals() throws {
        let s = try AtomicStructure(box: SIMD3(10, 10, 2), atoms: [.init(atomicNumber: 104, position: SIMD3(1, 1, 1))])
        XCTAssertThrowsError(try ProjectedPotential(structure: s, nx: 64, ny: 64, sliceThickness: 1))
        let ok = try AtomicStructure(box: SIMD3(10, 10, 2), atoms: [.init(atomicNumber: 13, position: SIMD3(1, 1, 1))])
        XCTAssertThrowsError(try ProjectedPotential(structure: ok, nx: 64, ny: 64, sliceThickness: 0))
        XCTAssertThrowsError(try ProjectedPotential(structure: ok, nx: 64, ny: 64, sliceThicknesses: [1, 0.5]))
        XCTAssertThrowsError(try ProjectedPotential(structure: ok, nx: 100, ny: 64, sliceThickness: 1))
        XCTAssertThrowsError(try AtomicStructure(box: SIMD3(0, 10, 2), atoms: []))
        XCTAssertThrowsError(try AtomicStructure(box: SIMD3(10, 10, 2), atoms: [.init(atomicNumber: 13, position: SIMD3(.nan, 1, 1))]))
    }
}

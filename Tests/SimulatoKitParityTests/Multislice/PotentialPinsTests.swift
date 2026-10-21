//
//  PotentialPinsTests.swift
//  Role: registration sections 3.4–3.6 against the pins: the slicing (thicknesses, the slice of every atom), the
//        projected potential of every slice (stored arrays and per-slice sums and maxima), the bandlimited transmission
//        function, truth 2 (the k = 0 identity) on every slice of every case, and the planted slicing rules.
//

import SimulatoKitTestSupport
import Testing
import SimulatoCore
import SimulatoMultislice

@Suite("Projected potential")
struct PotentialPinsTests {
    static let cases = ["A1", "A2", "B1", "C1", "D1"]

    @Test(arguments: cases)
    func slicingEqualsAbTEM(id: String) throws {
        let c = Pins.case(id)
        let potential = try Pins.potential(id, pi: .abTEMFloat32)
        #expect(potential.slices.map(\.thickness) == c["slice_thickness"] as! [Double], "\(id) thicknesses")
        #expect(potential.sliceOfAtom == c["slice_of_atom"] as! [Int], "\(id) slice of each atom")
    }

    @Test(arguments: cases)
    func potentialEqualsAbTEMWithItsPi(id: String) throws {
        let c = Pins.case(id)
        let potential = try Pins.potential(id, pi: .abTEMFloat32)
        let maxima = c["V_max_abs"] as! [Double], sums = c["V_sum"] as! [Double]
        let globalMax = maxima.max()!
        for (j, slice) in potential.slices.enumerated() {
            let m = slice.values.map(abs).max()!
            #expect(abs(m - maxima[j]) <= 1e-11 * max(globalMax, 1e-300), "\(id) slice \(j) max |V|")
            let s = slice.values.reduce(0, +)
            #expect(abs(s - sums[j]) <= 1e-11 * max(abs(sums[j]), 1), "\(id) slice \(j) ΣV")
            if Pins.hasArray("\(id).V.\(j)") {
                #expect(relMax(slice.values, Pins.doubles("\(id).V.\(j)")) <= 1e-11, "\(id) V slice \(j)")
            }
        }
    }

    @Test(arguments: ["A1", "A2", "B1", "D1"])
    func potentialWithTheExactPiIsWithinTheRegisteredDifference(id: String) throws {
        let potential = try Pins.potential(id, pi: .exact)
        for j in 0..<2 where Pins.hasArray("\(id).V.\(j)") {
            let d = relMax(potential.slices[j].values, Pins.doubles("\(id).V.\(j)"))
            #expect(d <= 1e-7, "\(id) exact-π V slice \(j)")
            #expect(d > 1e-9, "\(id): the exact π must differ from abTEM's float32 π (predicted ~5e−8)")
        }
    }

    @Test(arguments: ["A1", "A2", "B1", "D1"])
    func bandlimitedTransmissionEqualsAbTEM(id: String) throws {
        let potential = try Pins.potential(id, pi: .abTEMFloat32)
        let t = try Multislice.transmissionFunctions(potential, beam: Pins.beam(id))
        #expect(relMax(t[0], Pins.complex("\(id).t.0")) <= 1e-11, "\(id) t slice 0")
    }

    /// Truth 2: Σ V ΔxΔy = Σ_atoms f(0) for every slice (S(0) = ΔxΔy and the bilinear weights sum to one).
    @Test(arguments: cases, [LobatoParameterization.PiConstant.exact, .abTEMFloat32])
    func integratedPotentialEqualsTheSumOfFAtZero(id: String, pi: LobatoParameterization.PiConstant) throws {
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
                #expect(slice.values.map(abs).max()! == 0, "\(id) empty slice \(j) is zero exactly")
            } else {
                #expect(relDiff(integral, expected) <= 1e-12, "\(id) slice \(j) k = 0 identity")
            }
        }
    }

    /// Section 3.4, the top-face rule: an atom within 1e−10 Å of the top face is put at z = 0 (slice 0), so its potential
    /// equals that of the same atom at z = 0.
    @Test func topFaceAtomGoesToTheFirstSlice() throws {
        let top = try AtomicStructure(box: SIMD3(8, 8, 4), atoms: [.init(atomicNumber: 29, position: SIMD3(2.3, 5.1, 4 - 1e-11))])
        let bottom = try AtomicStructure(box: SIMD3(8, 8, 4), atoms: [.init(atomicNumber: 29, position: SIMD3(2.3, 5.1, 0))])
        let a = try ProjectedPotential(structure: top, nx: 64, ny: 64, sliceThickness: 1)
        let b = try ProjectedPotential(structure: bottom, nx: 64, ny: 64, sliceThickness: 1)
        #expect(a.sliceOfAtom == [0])
        #expect(a.slices[0].values == b.slices[0].values)
        #expect(a.slices[3].values.map(abs).max() == 0)
    }

    /// Section 3.4, the nudged edges: 30 slices of 0.1 Å accumulate to 0.30000000000000004 after three, so without the
    /// −1e−12 an atom at z = 0.3 would fall into slice 2; abTEM puts it into slice 3 (the next slice).
    @Test func atomOnADriftedSliceBoundaryGoesToTheNextSlice() throws {
        #expect(0.1 + 0.1 + 0.1 > 0.3, "the premise: the cumulative sum overshoots")
        let s = try AtomicStructure(box: SIMD3(8, 8, 3), atoms: [.init(atomicNumber: 6, position: SIMD3(1, 1, 0.3))])
        let p = try ProjectedPotential(structure: s, nx: 64, ny: 64, sliceThickness: 0.1)
        #expect(p.slices.count == 30)
        #expect(p.sliceOfAtom == [3])
    }

    /// Section 3.4: a requested thickness that does not divide the depth gives n = ⌈Lz / s⌉ equal slices of Lz / n.
    @Test func requestedThicknessGivesCeilingManyEqualSlices() throws {
        #expect(try Slicing.thicknesses(depth: 5, requested: 2) == [5.0 / 3, 5.0 / 3, 5.0 / 3])
        #expect(try Slicing.thicknesses(depth: 40.495, requested: 2).count == 21)
        #expect(try Slicing.thicknesses(depth: 4, requested: 4) == [4])
        let s = try AtomicStructure(box: SIMD3(8, 8, 5), atoms: [.init(atomicNumber: 6, position: SIMD3(1, 1, 3.4))])
        let p = try ProjectedPotential(structure: s, nx: 64, ny: 64, sliceThickness: 2)
        #expect(p.slices.count == 3)
        #expect(p.sliceOfAtom == [2], "z = 3.4 Å is in the third slice of 5/3 Å")
    }

    @Test func refusals() throws {
        let s = try AtomicStructure(box: SIMD3(10, 10, 2), atoms: [.init(atomicNumber: 104, position: SIMD3(1, 1, 1))])
        #expect(throws: (any Error).self) { try ProjectedPotential(structure: s, nx: 64, ny: 64, sliceThickness: 1) }
        let ok = try AtomicStructure(box: SIMD3(10, 10, 2), atoms: [.init(atomicNumber: 13, position: SIMD3(1, 1, 1))])
        #expect(throws: (any Error).self) { try ProjectedPotential(structure: ok, nx: 64, ny: 64, sliceThickness: 0) }
        #expect(throws: (any Error).self) { try ProjectedPotential(structure: ok, nx: 64, ny: 64, sliceThicknesses: [1, 0.5]) }
        #expect(throws: (any Error).self) { try ProjectedPotential(structure: ok, nx: 1, ny: 64, sliceThickness: 1) }
        #expect(throws: (any Error).self) { try AtomicStructure(box: SIMD3(0, 10, 2), atoms: []) }
        // Refuter finding: a thickness so small that ⌈Lz / s⌉ overflows is a named refusal, not a trap.
        #expect(throws: SimulatoError.invalidSliceThickness(1e-320)) {
            try ProjectedPotential(structure: ok, nx: 64, ny: 64, sliceThickness: 1e-320)
        }
        #expect(throws: SimulatoError.invalidSliceThickness(1e-7)) { try Slicing.thicknesses(depth: 2, requested: 1e-7) }
        #expect(throws: (any Error).self) {
            try AtomicStructure(box: SIMD3(10, 10, 2), atoms: [.init(atomicNumber: 13, position: SIMD3(.nan, 1, 1))])
        }
    }
}

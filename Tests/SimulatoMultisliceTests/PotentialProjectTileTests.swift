//
//  PotentialProjectTileTests.swift
//  Role: `ProjectedPotential.project`, `.tiled` and `Potential.project` against abTEM 1.0.10 `PotentialArray.project`,
//        `.tile` and `Potential.project`: the numbers printed by tools/pins/potential-project-tile/pins.py held inline,
//        sums at rtol 1e-7 (assert_allclose default; abTEM states none), shapes, extents and thicknesses exact.
//        The planted potential is built from outside through the public `PotentialSlice` initialiser.
//

import Foundation
import Testing
import SimulatoCore
import SimulatoMultislice
import SimulatoKitTestSupport

@Suite("Potential project and tile")
struct PotentialProjectTileTests {
    /// abTEM's planted PotentialArray: 3 slices on 4 x 5, values k^1.3 (k the flat index), thicknesses 1, 1.5, 0.5 Å,
    /// extent 2 x 2.5 Å.
    static func planted() throws -> ProjectedPotential {
        let grid = try Grid(nx: 4, ny: 5, extentX: 2.0, extentY: 2.5)
        let slices = [1.0, 1.5, 0.5].enumerated().map { s, dz in
            PotentialSlice(thickness: dz, values: (0..<20).map { pow(Double(s * 20 + $0), 1.3) })
        }
        return ProjectedPotential(grid: grid, slices: slices)
    }

    /// abTEM's checksum: the sum of value · (flat index + 1), slice by slice, row-major.
    static func weighted(_ values: [Double]) -> Double {
        values.enumerated().reduce(0) { $0 + $1.element * Double($1.offset + 1) }
    }

    static func close(_ a: Double, _ b: Double) -> Bool { abs(a - b) <= 1e-7 * abs(b) }

    @Test func projectSumsTheSlices() throws {
        let p = try Self.planted().project()
        let abTEM = [170.0992068579204, 178.26246245022838, 186.96329271595343, 195.98402939465967, 205.25906320258483,
                     214.75256861125797, 224.4411280205199, 234.30787148190774, 244.33987262942208, 254.52679454798368,
                     264.8601082341936, 275.3326063892693, 285.9380832382715, 296.6711140312202, 307.5268975939282,
                     318.50114048503565, 329.5899695928622, 340.7898647593841, 352.0976058713826, 363.510230635932]
        #expect(p.nx == 4 && p.ny == 5 && p.samplingX == 0.5 && p.samplingY == 0.5)
        #expect(p.values.count == abTEM.count)
        for (a, b) in zip(p.values, abTEM) { #expect(Self.close(a, b), "\(a) against \(b)") }
    }

    @Test(arguments: [
        // (x, y, z) repetitions, shape (slices, nx, ny), extent, weighted, sum, projected weighted
        ((2, 3, 1), (3, 8, 15), (4.0, 7.5), 7810729.876804397, 31462.523464463502, 2021200.3158708385),
        ((2, 1, 2), (6, 8, 5), (4.0, 2.5), 3002294.8578713303, 20975.015642975668, 457231.794640888),
        ((1, 1, 3), (9, 4, 5), (2.0, 2.5), 1611947.7293367172, 15731.261732231753, 185611.2286583485),
    ])
    func tiledMatchesAbTEM(reps: (Int, Int, Int), shape: (Int, Int, Int), extent: (Double, Double),
                           weighted: Double, sum: Double, projectedWeighted: Double) throws {
        let original = try Self.planted()
        let t = try original.tiled(repetitionsX: reps.0, repetitionsY: reps.1, repetitionsZ: reps.2)
        #expect(t.slices.count == shape.0 && t.grid.nx == shape.1 && t.grid.ny == shape.2)
        #expect(t.grid.extentX == extent.0 && t.grid.extentY == extent.1)
        #expect(t.grid.samplingX == 0.5 && t.grid.samplingY == 0.5)
        let thicknesses = [Double](repeating: 0, count: reps.2).flatMap { _ in [1.0, 1.5, 0.5] }
        #expect(t.slices.map(\.thickness) == thicknesses && t.thickness == 3.0 * Double(reps.2))
        #expect(t.sliceOfAtom.isEmpty)
        let flat = t.slices.flatMap(\.values)
        #expect(Self.close(Self.weighted(flat), weighted))
        #expect(Self.close(flat.reduce(0, +), sum))
        #expect(Self.close(Self.weighted(t.project().values), projectedWeighted))
    }

    @Test func tilingOnceIsTheIdentity() throws {
        let original = try Self.planted()
        #expect(try original.tiled(repetitionsX: 1, repetitionsY: 1) == ProjectedPotential(grid: original.grid,
                                                                                             slices: original.slices))
    }

    @Test func aluminiumProjectionMatchesAbTEM() throws {
        let potential = try Potential(structure: aluminiumStructure(cells: (1, 1, 2)), nx: 16, ny: 16, sliceThickness: 2.025,
                                      parameterization: LobatoParameterization(piConstant: .abTEMFloat32))
        let p = try potential.project()
        #expect(p == (try potential.build().project()), "the lazy projection is the built one")
        #expect(potential.sliceCount == 4)
        #expect(Self.close(p.sum, 35209.17454763371), "\(p.sum)")
        #expect(Self.close(p.max, 907.726975700707), "\(p.max)")
        #expect(Self.close(Self.weighted(p.values), 4215022.430595551))
    }

    @Test func tilingRefusesBadRepetitions() throws {
        let p = try Self.planted()
        // Below 1, a slice past 2^28 values, past the slice ceiling (3 · 400 000), past 2^31 values in all, overflow.
        for (x, y, z) in [(0, 1, 1), (1, 0, 1), (1, 1, 0), (-1, 2, 2), (1 << 20, 1 << 20, 1), (1, 1, 400_000),
                          (1_000, 1_000, 200), (1, 1, Int.max), (Int.max / 2, 1, 1)] {
            #expect(throws: ImageOperationError.self) { try p.tiled(repetitionsX: x, repetitionsY: y, repetitionsZ: z) }
        }
        // A potential without slices is refused the same way before any grid is built.
        let empty = ProjectedPotential(grid: p.grid, slices: [])
        #expect(throws: ImageOperationError.self) { try empty.tiled(repetitionsX: Int.max / 2, repetitionsY: 1) }
        #expect(try empty.tiled(repetitionsX: 2, repetitionsY: 2, repetitionsZ: 3).slices.isEmpty)
    }
}

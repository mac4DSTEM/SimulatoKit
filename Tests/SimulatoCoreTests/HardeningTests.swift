//
//  HardeningTests.swift
//  Role: inputs that used to trap late in SimulatoCore (huge discs, tilings and resamples) now throw their typed errors.
//

import Foundation
import Testing
import SimulatoCore

@Suite("Core hardening")
struct CoreHardeningTests {
    func image() -> RealField {
        RealField(nx: 4, ny: 4, values: [Double](repeating: 1, count: 16), samplingX: 1, samplingY: 1)
    }

    @Test func hugeDiscIsRefused() throws {
        let f = image()
        #expect(throws: ImageOperationError.self) { try f.integrateDisc(x: 1, y: 1, radius: 1e300) }
        #expect(throws: ImageOperationError.self) { try f.integrateDisc(x: 1e300, y: 1, radius: 1) }
        #expect(throws: ImageOperationError.self) { try f.integrateDisc(x: 1, y: -1e300, radius: 1) }
        // a valid disc is unchanged: the 3 x 3 box of ones sums to at most 9
        #expect(try f.integrateDisc(x: 2, y: 2, radius: 1) > 0)
    }

    @Test func hugeTilingIsRefused() throws {
        let f = image()
        #expect(throws: ImageOperationError.invalidRepetitions(.max, 2)) { try f.tiled(repetitionsX: .max, repetitionsY: 2) }
        #expect(throws: ImageOperationError.self) { try f.tiled(repetitionsX: 1 << 20, repetitionsY: 1 << 20) }
        #expect(try f.tiled(repetitionsX: 2, repetitionsY: 3).values.count == 8 * 12)
    }

    @Test func hugeResampleIsRefused() throws {
        let grid = try Grid(nx: 4, ny: 4, extentX: 1, extentY: 1)
        let w = Waves(field: ComplexField(nx: 4, ny: 4), grid: grid, beam: try ElectronBeam(energy: 200e3), normalization: .values)
        #expect(throws: SimulatoError.self) { try w.resampled(rows: 1 << 40, columns: 1 << 30) }
        #expect(throws: SimulatoError.self) { try w.resampled(rows: 1 << 20, columns: 1 << 20) }
        #expect(try w.resampled(rows: 8, columns: 8).grid.nx == 8)
    }
}

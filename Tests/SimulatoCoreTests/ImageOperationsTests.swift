//
//  ImageOperationsTests.swift
//  Role: the image operations by their own truth: a known smooth potential recovered from its gradient by
//        `integratingGradient` up to its minimum (the planted-truth test), tiling, the soft-edged disc, and the refusals.
//        The comparison with abTEM's numbers is in SimulatoKitParityTests (icom-pins.json).
//

import Foundation
import Testing
@testable import SimulatoCore

@Suite("Image operations")
struct ImageOperationsTests {
    /// A potential of periodic cosines inside the grid's band and its analytic gradient.
    static func planted(nx: Int, ny: Int, sx: Double, sy: Double) -> (potential: RealField, gx: RealField, gy: RealField) {
        let ex = Double(nx) * sx, ey = Double(ny) * sy
        var t = [Double](repeating: 0, count: nx * ny), gx = t, gy = t
        for (m, n, a, phase) in [(1, 0, 1.0, 0.3), (0, 1, 0.7, 1.1), (1, 1, 0.5, 2.0), (2, -1, 0.9, 4.0), (1, 2, 0.4, 5.2)]
        where 2 * abs(m) < nx && 2 * abs(n) < ny {
            let kx = Double(m) / ex, ky = Double(n) / ey
            for ix in 0..<nx {
                for iy in 0..<ny {
                    let arg = 2 * Double.pi * (kx * Double(ix) * sx + ky * Double(iy) * sy) + phase
                    t[ix * ny + iy] += a * cos(arg)
                    gx[ix * ny + iy] += -a * 2 * Double.pi * kx * sin(arg)
                    gy[ix * ny + iy] += -a * 2 * Double.pi * ky * sin(arg)
                }
            }
        }
        func field(_ v: [Double]) -> RealField { RealField(nx: nx, ny: ny, values: v, samplingX: sx, samplingY: sy) }
        return (field(t), field(gx), field(gy))
    }

    @Test(arguments: [(16, 12, 0.3, 0.5), (15, 11, 0.4, 0.4), (20, 9, 0.25, 0.6), (64, 64, 0.1, 0.1), (45, 49, 0.2, 0.3)])
    func integratingAGradientRecoversThePotentialUpToItsMinimum(nx: Int, ny: Int, sx: Double, sy: Double) throws {
        let p = Self.planted(nx: nx, ny: ny, sx: sx, sy: sy)
        let t = try RealField.integratingGradient(x: p.gx, y: p.gy)
        let low = p.potential.values.min()!
        let range = p.potential.values.max()! - low
        var worst = 0.0
        for i in t.values.indices { worst = Swift.max(worst, abs(t.values[i] - (p.potential.values[i] - low))) }
        #expect(worst <= 1e-9 * range, "worst \(worst) of range \(range)")
        #expect(t.values.min() == 0)
        #expect(t.samplingX == sx && t.samplingY == sy)
    }

    @Test func aConstantGradientCarriesNoPeriodicPart() throws {
        // a uniform gradient has only the zero frequency, which abTEM sets to 0: the integral is flat
        let g = RealField(nx: 8, ny: 6, values: [Double](repeating: 2.5, count: 48), samplingX: 1, samplingY: 1)
        let t = try RealField.integratingGradient(x: g, y: g)
        #expect(t.values.allSatisfy { abs($0) < 1e-12 })
    }

    @Test func complexGradientIsTheSameAsTwoComponents() throws {
        let p = Self.planted(nx: 15, ny: 11, sx: 0.4, sy: 0.4)
        let c = ComplexField(nx: 15, ny: 11, real: p.gx.values, imaginary: p.gy.values)
        #expect(try c.integratedGradient(samplingX: 0.4, samplingY: 0.4) == RealField.integratingGradient(x: p.gx, y: p.gy))
    }

    @Test func gradientRefusals() throws {
        let a = RealField(nx: 4, ny: 4, values: [Double](repeating: 0, count: 16), samplingX: 1, samplingY: 1)
        let b = RealField(nx: 4, ny: 3, values: [Double](repeating: 0, count: 12), samplingX: 1, samplingY: 1)
        let c = RealField(nx: 4, ny: 4, values: [Double](repeating: 0, count: 16), samplingX: 2, samplingY: 1)
        let d = RealField(nx: 1, ny: 4, values: [Double](repeating: 0, count: 4), samplingX: 1, samplingY: 1)
        let e = RealField(nx: 4, ny: 4, values: [Double](repeating: 0, count: 16), samplingX: 0, samplingY: 0)
        #expect(throws: ImageOperationError.gradientComponentsDiffer) { _ = try RealField.integratingGradient(x: a, y: b) }
        #expect(throws: ImageOperationError.gradientComponentsDiffer) { _ = try RealField.integratingGradient(x: a, y: c) }
        #expect(throws: SimulatoError.unsupportedGridLength(1)) { _ = try RealField.integratingGradient(x: d, y: d) }
        #expect(throws: SimulatoError.invalidSampling(0, 0)) { _ = try RealField.integratingGradient(x: e, y: e) }
    }

    @Test func frequenciesAreNumPys() {
        #expect(fftFrequencies(4, spacing: 0.5) == [0, 0.5, -1, -0.5])
        #expect(fftFrequencies(5, spacing: 1) == [0, 0.2, 0.4, -0.4, -0.2])
    }

    @Test func tilingRepeatsTheImage() throws {
        let f = RealField(nx: 2, ny: 3, values: [1, 2, 3, 4, 5, 6], samplingX: 0.5, samplingY: 0.25)
        let t = try f.tiled(repetitionsX: 2, repetitionsY: 2)
        #expect(t.nx == 4 && t.ny == 6 && t.samplingX == 0.5 && t.samplingY == 0.25)
        #expect(t.values == [1, 2, 3, 1, 2, 3, 4, 5, 6, 4, 5, 6, 1, 2, 3, 1, 2, 3, 4, 5, 6, 4, 5, 6])
        #expect(throws: ImageOperationError.invalidRepetitions(0, 2)) { _ = try f.tiled(repetitionsX: 0, repetitionsY: 2) }
    }

    @Test func discOfAConstantImage() throws {
        let f = RealField(nx: 20, ny: 16, values: [Double](repeating: 3, count: 320), samplingX: 0.2, samplingY: 0.25)
        #expect(abs(try f.integrateDisc(x: 2, y: 2, radius: 0.7, mean: true) - 3) < 1e-12)
        // periodic: the same disc moved by one extent along either axis gives the same sum
        let r = RealField(nx: 20, ny: 16, values: (0..<320).map { sin(Double($0)) }, samplingX: 0.2, samplingY: 0.25)
        let a = try r.integrateDisc(x: 0.1, y: 0.2, radius: 0.5)
        let b = try r.integrateDisc(x: 0.1 + 4.0, y: 0.2 - 4.0, radius: 0.5)
        #expect(abs(a - b) < 1e-12)
        #expect(throws: ImageOperationError.discOutsideImage) { _ = try r.integrateDisc(x: 0.1, y: 0.2, radius: 0.5, border: .raise) }
        #expect(try r.integrateDisc(x: 2, y: 2, radius: 0.5, border: .raise) == r.integrateDisc(x: 2, y: 2, radius: 0.5))
        #expect(throws: ImageOperationError.self) { _ = try r.integrateDisc(x: .nan, y: 0, radius: 1) }
        #expect(throws: ImageOperationError.self) { _ = try r.integrateDisc(x: 0, y: 0, radius: 0) }
    }
}

//
//  SimulatoCoreTests.swift
//  Role: SimulatoCore's own rules: the grid lengths (vDSP's), NumPy's fftfreq, the FFT against a direct DFT and its
//        round trip, the antialias aperture, the pattern crop sizes, the structure's repeat order.
//

import Foundation
import Testing
@testable import SimulatoCore

@Suite("SimulatoCore")
struct SimulatoCoreTests {
    @Test func everyLengthFromTwoIsSupported() {
        #expect((2...1100).allSatisfy(Grid.isSupportedLength))
        #expect(!Grid.isSupportedLength(1) && !Grid.isSupportedLength(0) && !Grid.isSupportedLength(-8))
        let direct = (1...1100).filter(Grid.hasDirectTransform)
        #expect(direct == [8, 16, 24, 32, 40, 48, 64, 80, 96, 120, 128, 160, 192, 240, 256, 320, 384, 480, 512, 640, 768, 960, 1024])
    }

    /// Any length against a direct DFT (vDSP lengths, odd, prime, and sizes 2 and 7 below vDSP's floor), both ways.
    @Test(arguments: [(2, 3), (7, 11), (13, 100), (97, 6), (45, 49), (127, 20), (8, 8), (100, 90), (1000, 3)])
    func fFTOfAnyLengthMatchesADirectDFT(nx: Int, ny: Int) throws {
        var f = ComplexField(nx: nx, ny: ny)
        for i in 0..<(nx * ny) {
            f.real[i] = sin(Double(i) * 0.37) + 0.1 * Double(i % 7)
            f.imaginary[i] = cos(Double(i) * 0.11) - 0.05 * Double(i % 5)
        }
        let original = f
        let fft = try FFT2D(nx: nx, ny: ny)
        fft.forward(&f)
        var worst = 0.0
        var scale = 0.0
        for (k1, k2) in [(0, 0), (1, 1 % ny), (nx / 2, ny / 3), (nx - 1, ny - 1), (nx / 3, ny - 1)] {
            var re = 0.0, im = 0.0
            for j1 in 0..<nx {
                for j2 in 0..<ny {
                    let ph = -2 * Double.pi * (Double((k1 * j1) % nx) / Double(nx) + Double((k2 * j2) % ny) / Double(ny))
                    let a = original.real[j1 * ny + j2], b = original.imaginary[j1 * ny + j2]
                    re += a * cos(ph) - b * sin(ph)
                    im += a * sin(ph) + b * cos(ph)
                }
            }
            worst = max(worst, abs(re - f.real[k1 * ny + k2]), abs(im - f.imaginary[k1 * ny + k2]))
            scale = max(scale, abs(re), abs(im))
        }
        // numpy.testing.assert_allclose's default rtol is 1e-7; the transform is good to ~1e-13 of the largest term.
        #expect(worst <= 1e-7 * max(scale, 1))
        fft.inverse(&f)
        var back = 0.0
        for i in f.real.indices { back = max(back, abs(f.real[i] - original.real[i]), abs(f.imaginary[i] - original.imaginary[i])) }
        #expect(back <= 1e-7)
    }

    @Test func fFTFreqIsNumPys() {
        #expect(Grid.fftfreq(8, 0.5) == [0, 0.25, 0.5, 0.75, -1, -0.75, -0.5, -0.25])
        #expect(Grid.fftfreq(15, 1).prefix(9).map { $0 * 15 } == [0, 1, 2, 3, 4, 5, 6, 7, -7])
    }

    @Test func fFTMatchesADirectDFTAndRoundTrips() throws {
        let nx = 24, ny = 40
        var f = ComplexField(nx: nx, ny: ny)
        for i in 0..<(nx * ny) {
            f.real[i] = sin(Double(i) * 0.37) + 0.1 * Double(i % 7)
            f.imaginary[i] = cos(Double(i) * 0.11) - 0.05 * Double(i % 5)
        }
        let original = f
        let fft = try FFT2D(nx: nx, ny: ny)
        fft.forward(&f)
        var worst = 0.0
        for (k1, k2) in [(0, 0), (1, 3), (5, 17), (23, 39), (12, 20)] {
            var re = 0.0, im = 0.0
            for j1 in 0..<nx {
                for j2 in 0..<ny {
                    let ph = -2 * Double.pi * (Double(k1 * j1) / Double(nx) + Double(k2 * j2) / Double(ny))
                    let a = original.real[j1 * ny + j2], b = original.imaginary[j1 * ny + j2]
                    re += a * cos(ph) - b * sin(ph)
                    im += a * sin(ph) + b * cos(ph)
                }
            }
            worst = max(worst, abs(re - f.real[k1 * ny + k2]), abs(im - f.imaginary[k1 * ny + k2]))
        }
        #expect(worst <= 1e-11)
        fft.inverse(&f)
        var back = 0.0
        for i in f.real.indices { back = max(back, abs(f.real[i] - original.real[i]), abs(f.imaginary[i] - original.imaginary[i])) }
        #expect(back <= 1e-14)
    }

    @Test func antialiasApertureAndCutoffSizes() throws {
        let grid = try Grid(nx: 128, ny: 128, extentX: 10, extentY: 10)
        let a = Antialias.aperture(on: grid)
        #expect(a[0] == 1)
        #expect(a[64 * 128 + 64] == 0, "the corner frequency is cut")
        #expect(a.allSatisfy { $0 >= 0 && $0 <= 1 })
        #expect(Antialias.cutoffPatternSize(on: grid).nx == 87, "128 px of 0.078125 Å: 85.33 → 85 → 86 → 87")
        let odd = try Grid(nx: 120, ny: 160, extentX: 11.715, extentY: 15.62)
        let s = Antialias.cutoffPatternSize(on: odd)
        #expect(s.nx == 81 && s.ny == 107)
    }

    @Test func repeatFollowsASEOrder() throws {
        let s = try AtomicStructure(box: SIMD3(1, 2, 3), atoms: [.init(atomicNumber: 1, position: SIMD3(0.1, 0.2, 0.3)),
                                                                   .init(atomicNumber: 2, position: SIMD3(0.5, 0.5, 0.5))])
        let r = s.repeated(2, 1, 2)
        #expect(r.box == SIMD3(2, 2, 6))
        #expect(r.atoms.map(\.atomicNumber) == [1, 2, 1, 2, 1, 2, 1, 2])
        #expect(r.atoms[2].position == SIMD3(0.1, 0.2, 3.3))
        #expect(r.atoms[4].position == SIMD3(1.1, 0.2, 0.3))
    }
}

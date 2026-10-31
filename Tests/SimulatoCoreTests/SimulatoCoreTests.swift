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

    /// Lengths off vDSP take the mixed-radix FFT, not Bluestein (issue 4: 250 = 2 · 5³ ran three padded 512-point DFTs
    /// per line, 506 = 2 · 11 · 23 the same at 1024), unless a large prime factor makes Bluestein the cheaper
    /// (`MixedRadixPlan.prefersBluestein`): 97, 127, 509 and 1009 are prime, 2 · 53 and 2 · 257 have a factor above √n.
    @Test func lengthsPickTheirFFTPath() throws {
        func paths(_ n: Int) throws -> FFTPath {
            let p = try FFT2D(nx: n, ny: 2).paths
            #expect(try FFT2D(nx: 2, ny: n).paths.y == p.x)
            return p.x
        }
        #expect(try paths(256) == .vDSP && paths(480) == .vDSP)
        for n in [2, 3, 6, 7, 11, 13, 45, 47, 49, 90, 100, 250, 496, 500, 506, 1000, 1001] {
            #expect(try paths(n) == .mixedRadix, "\(n)")
        }
        for n in [97, 106, 127, 509, 514, 1009] { #expect(try paths(n) == .bluestein, "\(n)") }
    }

    /// pocketfft's `cost_guess` and the Bluestein test of `make_cfft_plan`, re-weighted for vDSP, at worked values.
    @Test func mixedRadixCostRuleIsPocketffts() {
        #expect(MixedRadixPlan.factorization(506) == [2, 11, 23])
        #expect(MixedRadixPlan.factorization(1000) == [2, 4, 5, 5, 5])
        #expect(MixedRadixPlan.factorization(1) == nil)
        #expect(abs(MixedRadixPlan.costGuess(506) - (2 + 1.1 * 11 + 1.1 * 23) * 506) < 1e-9)
        #expect(MixedRadixPlan.costGuess(1024) == 20 * 1024)
        #expect(!MixedRadixPlan.prefersBluestein(47, bluesteinLength: 128))     // below 50
        #expect(!MixedRadixPlan.prefersBluestein(506, bluesteinLength: 1024))   // cost 19 936, Bluestein 25 600
        #expect(MixedRadixPlan.prefersBluestein(97, bluesteinLength: 256))      // cost 10 350, Bluestein 5 120
        #expect(!MixedRadixPlan.prefersBluestein(1001, bluesteinLength: 2048))  // 13 ≤ √1001
    }

    /// Every element of the transform at the lengths of issue 4 (and odd, generic-radix and Bluestein neighbours) against a
    /// separable direct DFT, at numpy.testing.assert_allclose's float64 rtol 1e-7 of the largest term, both ways.
    @Test(arguments: [(250, 250), (250, 506), (500, 21), (49, 250), (75, 3), (506, 90), (256, 250), (250, 480),
                      (97, 11), (13, 1001), (496, 127), (47, 2 * 3 * 3 * 11)])
    func fFTAtTheMixedRadixLengthsMatchesADirectDFTEverywhere(nx: Int, ny: Int) throws {
        var f = ComplexField(nx: nx, ny: ny)
        for i in 0..<(nx * ny) {
            f.real[i] = sin(Double(i) * 0.37) + 0.1 * Double(i % 7)
            f.imaginary[i] = cos(Double(i) * 0.11) - 0.05 * Double(i % 5)
        }
        let original = f
        // Reference: a direct DFT along y for every row, then along x for every column (exact integer phases).
        func directDFT(_ re: [Double], _ im: [Double]) -> ([Double], [Double]) {
            let n = re.count
            var outRe = [Double](repeating: 0, count: n), outIm = outRe
            for k in 0..<n {
                var sr = 0.0, si = 0.0
                for j in 0..<n {
                    let ph = -2 * Double.pi * Double((k * j) % n) / Double(n)
                    let c = cos(ph), s = sin(ph)
                    sr += re[j] * c - im[j] * s
                    si += re[j] * s + im[j] * c
                }
                outRe[k] = sr
                outIm[k] = si
            }
            return (outRe, outIm)
        }
        var refRe = original.real, refIm = original.imaginary
        for ix in 0..<nx {
            let (a, b) = directDFT(Array(refRe[(ix * ny)..<(ix * ny + ny)]), Array(refIm[(ix * ny)..<(ix * ny + ny)]))
            for iy in 0..<ny { refRe[ix * ny + iy] = a[iy]; refIm[ix * ny + iy] = b[iy] }
        }
        for iy in 0..<ny {
            let (a, b) = directDFT((0..<nx).map { refRe[$0 * ny + iy] }, (0..<nx).map { refIm[$0 * ny + iy] })
            for ix in 0..<nx { refRe[ix * ny + iy] = a[ix]; refIm[ix * ny + iy] = b[ix] }
        }
        let fft = try FFT2D(nx: nx, ny: ny)
        fft.forward(&f)
        var worst = 0.0, scale = 0.0
        for i in f.real.indices {
            worst = max(worst, abs(f.real[i] - refRe[i]), abs(f.imaginary[i] - refIm[i]))
            scale = max(scale, abs(refRe[i]), abs(refIm[i]))
        }
        #expect(worst <= 1e-7 * scale)
        fft.inverse(&f)
        var back = 0.0
        for i in f.real.indices { back = max(back, abs(f.real[i] - original.real[i]), abs(f.imaginary[i] - original.imaginary[i])) }
        #expect(back <= 1e-7)
    }

    /// The inverse against a direct inverse DFT (e^{+2πi jk/n} / n) at lengths with two generic passes in a row
    /// (121 = 11², 169 = 13²) and a generic pass after radix 2 (198 = 2 · 3² · 11), so the u_{p−j} twiddles are checked
    /// in both directions, not only by the round trip.
    @Test(arguments: [(121, 169), (169, 198)])
    func inverseFFTAtGenericLengthsMatchesADirectInverseDFT(nx: Int, ny: Int) throws {
        var f = ComplexField(nx: nx, ny: ny)
        for i in 0..<(nx * ny) {
            f.real[i] = sin(Double(i) * 0.29) - 0.1 * Double(i % 3)
            f.imaginary[i] = cos(Double(i) * 0.13) + 0.07 * Double(i % 4)
        }
        let original = f
        func directInverse(_ re: [Double], _ im: [Double]) -> ([Double], [Double]) {
            let n = re.count
            var outRe = [Double](repeating: 0, count: n), outIm = outRe
            for k in 0..<n {
                var sr = 0.0, si = 0.0
                for j in 0..<n {
                    let ph = 2 * Double.pi * Double((k * j) % n) / Double(n)
                    let c = cos(ph), s = sin(ph)
                    sr += re[j] * c - im[j] * s
                    si += re[j] * s + im[j] * c
                }
                outRe[k] = sr / Double(n)
                outIm[k] = si / Double(n)
            }
            return (outRe, outIm)
        }
        var refRe = original.real, refIm = original.imaginary
        for ix in 0..<nx {
            let (a, b) = directInverse(Array(refRe[(ix * ny)..<(ix * ny + ny)]), Array(refIm[(ix * ny)..<(ix * ny + ny)]))
            for iy in 0..<ny { refRe[ix * ny + iy] = a[iy]; refIm[ix * ny + iy] = b[iy] }
        }
        for iy in 0..<ny {
            let (a, b) = directInverse((0..<nx).map { refRe[$0 * ny + iy] }, (0..<nx).map { refIm[$0 * ny + iy] })
            for ix in 0..<nx { refRe[ix * ny + iy] = a[ix]; refIm[ix * ny + iy] = b[ix] }
        }
        let fft = try FFT2D(nx: nx, ny: ny)
        #expect(fft.paths.x == .mixedRadix && fft.paths.y == .mixedRadix)
        fft.inverse(&f)
        var worst = 0.0, scale = 0.0
        for i in f.real.indices {
            worst = max(worst, abs(f.real[i] - refRe[i]), abs(f.imaginary[i] - refIm[i]))
            scale = max(scale, abs(refRe[i]), abs(refIm[i]))
        }
        #expect(worst <= 1e-7 * scale)
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

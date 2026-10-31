//
//  FFTBenchmarkTests.swift
//  Role: a local timing of the float64 FFT2D at the lengths of issue 4 (250², 506²) beside their vDSP neighbours and
//        lengths near the generic-radix / Bluestein boundary (97, 106, 127, 251, 509, 1001), the mixed-radix passes
//        forced at Bluestein lengths, and the vDSP DFT at the Bluestein padded length, so the share of the time spent in
//        vDSP against the rest is visible and `MixedRadixPlan.prefersBluestein` can be re-checked.
//        Skipped unless SIMULATO_FFT_BENCH=1; it asserts nothing about speed (CI machines differ).
//

import Accelerate
import Foundation
import Testing
@testable import SimulatoCore

@Suite("FFT benchmark", .serialized, .enabled(if: ProcessInfo.processInfo.environment["SIMULATO_FFT_BENCH"] == "1"))
struct FFTBenchmarkTests {
    /// Median wall time in ms of `repeats` calls of `body`, after one warm-up call; `setup` runs untimed before each call.
    private func medianMs(repeats: Int, setup: () -> Void = {}, _ body: () -> Void) -> Double {
        setup()
        body()
        var times: [Double] = []
        for _ in 0..<repeats {
            setup()
            let t0 = DispatchTime.now().uptimeNanoseconds
            body()
            times.append(Double(DispatchTime.now().uptimeNanoseconds - t0) / 1e6)
        }
        return times.sorted()[times.count / 2]
    }

    @Test(arguments: [97, 106, 127, 250, 251, 256, 500, 506, 509, 512, 1001, 1009])
    func forwardFFT2D(n: Int) throws {
        var f = ComplexField(nx: n, ny: n)
        for i in 0..<(n * n) {
            f.real[i] = sin(Double(i) * 0.37)
            f.imaginary[i] = cos(Double(i) * 0.11)
        }
        let fft = try FFT2D(nx: n, ny: n)
        let seed = f
        // The field is reset to the seed before each call (untimed), so the values stay bounded.
        let ms = medianMs(repeats: 21, setup: { f = seed; f.real[0] += 0 }) { fft.forward(&f) }
        print("FFTBENCH fft2 n=\(n) path=\(fft.paths.x) \(String(format: "%.3f", ms)) ms")
    }

    /// The mixed-radix plan forced at lengths the cost rule sends to Bluestein (2 n lines, the work of an n² fft2), to
    /// set beside `forwardFFT2D` at the same n when the rule is re-checked.
    @Test(arguments: [97, 106, 127, 251, 509])
    func mixedRadixAtBluesteinLengths(n: Int) throws {
        let plan = try #require(MixedRadixPlan(length: n))
        var re = (0..<(n * n)).map { sin(Double($0) * 0.37) }, im = (0..<(n * n)).map { cos(Double($0) * 0.11) }
        let seedRe = re, seedIm = im
        let work = UnsafeMutablePointer<Lane2>.allocate(capacity: plan.workCount)
        defer { work.deallocate() }
        let ms = medianMs(repeats: 21, setup: { re = seedRe; im = seedIm }) {
            re.withUnsafeMutableBufferPointer { r in im.withUnsafeMutableBufferPointer { i in
                plan.transformLines(count: n, forward: true, stride: n, r.baseAddress!, i.baseAddress!, work: work)
                plan.transformColumns(count: n, forward: true, stride: n, r.baseAddress!, i.baseAddress!, work: work)
            } }
        }
        print("FFTBENCH forced-mixed n=\(n) \(String(format: "%.3f", ms)) ms")
    }

    /// The vDSP DFT alone at a padded length, called 4 n times (two per line, 2 n lines): the floor of the Bluestein path.
    @Test(arguments: [(250, 512), (506, 1024)])
    func vDSPAtPaddedLength(n: Int, m: Int) throws {
        let setup = try #require(vDSP_DFT_zop_CreateSetupD(nil, vDSP_Length(m), .FORWARD))
        defer { vDSP_DFT_DestroySetupD(setup) }
        var ar = [Double](repeating: 0.5, count: m), ai = [Double](repeating: 0.25, count: m)
        var br = [Double](repeating: 0, count: m), bi = [Double](repeating: 0, count: m)
        let ms = medianMs(repeats: 21) {
            for _ in 0..<(4 * n) { vDSP_DFT_ExecuteD(setup, &ar, &ai, &br, &bi) }
        }
        print("FFTBENCH vdsp m=\(m) x\(4 * n) \(String(format: "%.3f", ms)) ms")
    }
}

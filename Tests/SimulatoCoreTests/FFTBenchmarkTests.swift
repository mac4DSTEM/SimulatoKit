//
//  FFTBenchmarkTests.swift
//  Role: a local timing of the float64 FFT2D at the lengths of issue 4 (250², 506²) beside their vDSP neighbours, and of
//        the vDSP DFT at the Bluestein padded length, so the share of the time spent in vDSP against the rest is visible.
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

    @Test(arguments: [250, 256, 500, 506, 512])
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

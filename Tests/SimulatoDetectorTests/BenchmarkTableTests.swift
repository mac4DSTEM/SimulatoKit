//
//  BenchmarkTableTests.swift
//  Role: the per-function table of docs/benchmarks.md (tools/benchmarks/run.sh). One warm-up, then five timed runs per
//        entry; prints `BENCH|function|size|backend|median|min|max` in seconds. Reported, never gated. Enabled only by
//        SIMULATO_BENCH_TABLE=<comma list of grid sizes> (for example 128,256); a size above 512 is meant for the nightly
//        review. Run under ../macEM-suite/tools/memcap.sh by run.sh.
//

import Foundation
import Testing
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe
import SimulatoDetector

private let benchSizes: [Int] = (ProcessInfo.processInfo.environment["SIMULATO_BENCH_TABLE"] ?? "")
    .split(separator: ",").compactMap { Int($0) }

@Suite("Benchmark table (SIMULATO_BENCH_TABLE=128,256)", .enabled(if: !benchSizes.isEmpty, "set SIMULATO_BENCH_TABLE"))
struct BenchmarkTableTests {
    static func seconds(_ d: Duration) -> Double {
        Double(d.components.seconds) + Double(d.components.attoseconds) * 1e-18
    }

    /// One warm-up, then five runs.
    static func report(_ function: String, _ size: Int, _ backend: String, _ body: () throws -> Void) rethrows {
        try body()
        let clock = ContinuousClock()
        var t: [Double] = []
        for _ in 0..<5 { t.append(seconds(try clock.measure(body))) }
        t.sort()
        print(String(format: "BENCH|%@|%d|%@|%.6f|%.6f|%.6f", function, size, backend, t[2], t[0], t[4]))
    }

    /// Al [001] at 200 keV, 5 slices of a/2 x 2 = 10 slices; `cells` unit cells across so the sampling is about 0.1 Å.
    static func structure(size n: Int) throws -> AtomicStructure {
        let cell = try AtomicStructure(box: SIMD3(4.05, 4.05, 4.05), atoms: [
            .init(atomicNumber: 13, position: SIMD3(0, 0, 0)), .init(atomicNumber: 13, position: SIMD3(0, 2.025, 2.025)),
            .init(atomicNumber: 13, position: SIMD3(2.025, 0, 2.025)), .init(atomicNumber: 13, position: SIMD3(2.025, 2.025, 0)),
        ])
        let cells = max(1, Int((Double(n) * 0.1 / 4.05).rounded()))
        return cell.repeated(cells, cells, 5)
    }

    @Test(arguments: benchSizes)
    func table(n: Int) throws {
        let beam = try ElectronBeam(energy: 200e3)
        let structure = try Self.structure(size: n)
        let potential = try Potential(structure: structure, nx: n, ny: n, sliceThickness: 2.025)
        let built = try potential.build()
        let wave = Waves.planeWave(beam: beam, grid: potential.grid)

        // FFT: forward + inverse of one complex field (vDSP lengths for n; Bluestein for n - 6).
        for m in [n, n - 6] {
            let fft = try FFT2D(nx: m, ny: m)
            var f = ComplexField(nx: m, ny: m, real: (0..<m * m).map { sin(Double($0)) }, imaginary: [Double](repeating: 0, count: m * m))
            try Self.report("fft2+ifft2", m, "cpu") {
                fft.forward(&f)
                fft.inverse(&f)
            }
        }
        try Self.report("potential build (infinite)", n, "cpu") { _ = try potential.build() }
        try Self.report("multislice 10 slices (plane wave, built potential)", n, "cpu") {
            _ = try Multislice.run(wave, through: built)
        }
        try Self.report("transmission functions", n, "cpu") { _ = try Multislice.transmissionFunctions(built, beam: beam) }
        let probe = try Probe(beam: beam, grid: potential.grid, aperture: Aperture(semiangleCutoff: 20))
        try Self.report("probe build", n, "cpu") { _ = try probe.build(at: SIMD2(1, 1)) }
        let scan = try GridScan(start: SIMD2(0, 0), end: SIMD2(4.05, 4.05), nx: 8, ny: 8)
        let detector = try AnnularDetector(inner: 40, outer: 90)
        try Self.report("STEM scan 8x8, annular", n, "cpu") {
            _ = try STEM.scan(probe: probe, potentials: [built], scan: scan, detectors: [detector])
        }
        let phonons = try FrozenPhonons(structure: structure, sigmas: [13: 0.08], configurations: 4, seed: 1)
        let ensemble = try phonons.lazyPotentials(nx: n, ny: n, sliceThickness: 2.025)
        try Self.report("frozen phonons 4 x multislice, built per configuration", n, "cpu") {
            for p in ensemble { _ = try Multislice.run(wave, through: p.build()) }
        }
        if Multislice.metalIsAvailable {
            try Self.report("multislice 10 slices (plane wave, built potential)", n, "metal") {
                _ = try Multislice.run(wave, through: built, backend: .metal)
            }
            try Self.report("multislice 10 slices (plane wave, potential built on the GPU)", n, "metal") {
                _ = try Multislice.run(wave, through: potential, backend: .metal)
            }
            try Self.report("STEM scan 8x8, annular", n, "metal") {
                _ = try STEM.scan(probe: probe, potentials: [built], scan: scan, detectors: [detector], backend: .metal)
            }
        }
    }
}

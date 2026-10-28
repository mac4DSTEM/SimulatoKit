//
//  BenchmarkTests.swift
//  Role: the speed table (reported, not asserted): Al [001] at 300 keV, 100 slices of a/2, about 0.047 Å per
//        pixel, a plane wave and a 64 × 64 HAADF scan over one cell, `.metal` against `.cpu`, from the structure to the
//        numbers. Enabled only by SIMULATO_BENCHMARK=1; the CPU scan at 512² and 1024² is timed on a subset of
//        positions and extrapolated per position (stated in the output).
//

import Foundation
import Testing
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe
import SimulatoDetector

@Suite("Benchmark (SIMULATO_BENCHMARK=1)",
       .enabled(if: ["1", "metal"].contains(ProcessInfo.processInfo.environment["SIMULATO_BENCHMARK"] ?? ""),
                "set SIMULATO_BENCHMARK=1 to run (SIMULATO_BENCHMARK=metal: the Metal side only)"),
       .enabled(if: Multislice.metalIsAvailable, "no Metal device"))
struct BenchmarkTests {
    static func seconds(_ d: Duration) -> Double {
        Double(d.components.seconds) + Double(d.components.attoseconds) * 1e-18
    }

    @Test(arguments: [(256, 3, 256), (512, 6, 64), (1024, 12, 16)])
    func speed(n: Int, cells: Int, cpuPositions: Int) throws {
        let cell = try AtomicStructure(box: SIMD3(4.05, 4.05, 4.05), atoms: [
            .init(atomicNumber: 13, position: SIMD3(0, 0, 0)), .init(atomicNumber: 13, position: SIMD3(0, 2.025, 2.025)),
            .init(atomicNumber: 13, position: SIMD3(2.025, 0, 2.025)), .init(atomicNumber: 13, position: SIMD3(2.025, 2.025, 0)),
        ])
        let structure = cell.repeated(cells, cells, 50)
        let beam = try ElectronBeam(energy: 300e3)
        let clock = ContinuousClock()
        // The plane wave, structure to exit wave.
        var times: [String: Double] = [:]
        let metalOnly = ProcessInfo.processInfo.environment["SIMULATO_BENCHMARK"] == "metal"
        for backend in metalOnly ? [Backend.metal] : [Backend.metal, .cpu] {
            let t = try clock.measure {
                let p = try Potential(structure: structure, nx: n, ny: n, sliceThickness: 2.025)
                _ = try Multislice.run(.planeWave(beam: beam, grid: p.grid), through: p, backend: backend)
            }
            times["plane \(backend)"] = Self.seconds(t)
        }
        // The scan: 64 × 64 over one cell, HAADF 70–200 mrad, probe 25 mrad.
        let full = try GridScan(start: SIMD2(0, 0), end: SIMD2(4.05, 4.05), nx: 64, ny: 64)
        let subset = try GridScan(start: SIMD2(0, 0), end: SIMD2(4.05, 4.05 * Double(cpuPositions) / 4096),
                                  nx: 64, ny: cpuPositions / 64 == 0 ? 1 : cpuPositions / 64)
        let cpuCount = subset.positions.count
        let detector = try AnnularDetector(inner: 70, outer: 200)
        let base = Multislice.metalPeakAllocatedBytes(reset: true)
        var peak = 0
        let tm = try clock.measure {
            let p = try Potential(structure: structure, nx: n, ny: n, sliceThickness: 2.025)
            let probe = try Probe(beam: beam, grid: p.grid, aperture: Aperture(semiangleCutoff: 25))
            _ = try STEM.scan(probe: probe, potentials: [p], scan: full, detectors: [detector], backend: .metal)
            peak = Multislice.metalPeakAllocatedBytes() - base
        }
        if metalOnly {
            print(String(format: "BENCH %d²: plane wave metal %.3f s; scan 64x64 metal %.2f s; peak device memory of the Metal scan %d MB above the start",
                         n, times["plane metal"]!, Self.seconds(tm), peak >> 20))
            return
        }
        let tc = try clock.measure {
            let p = try Potential(structure: structure, nx: n, ny: n, sliceThickness: 2.025)
            let probe = try Probe(beam: beam, grid: p.grid, aperture: Aperture(semiangleCutoff: 25))
            _ = try STEM.scan(probe: probe, potentials: [p], scan: subset, detectors: [detector], backend: .cpu)
        }
        let cpuScan = Self.seconds(tc) * 4096 / Double(cpuCount)
        print(String(format: "BENCH %d²: plane wave metal %.3f s, cpu %.3f s (x%.1f); scan 64x64 metal %.2f s, cpu %.2f s %@ (x%.1f); peak device memory of the Metal scan %d MB above the start",
                     n, times["plane metal"]!, times["plane cpu"]!, times["plane cpu"]! / times["plane metal"]!,
                     Self.seconds(tm), cpuScan, cpuCount == 4096 ? "measured" : "extrapolated from \(cpuCount) positions",
                     cpuScan / Self.seconds(tm), peak >> 20))
    }
}

//
//  WavesResampleTests.swift
//  Role: the Fourier-space resample behind WavesDetector(gpts) (abTEM's fft_interpolate); the pinned comparison with abTEM
//        is in SimulatoKitParityTests/WavesDetectorTests.
//

import Foundation
import Testing
import SimulatoCore

@Suite("Waves resampling")
struct WavesResampleTests {
    /// A smooth complex wave: a few low harmonics, so it survives a crop to 9 points and a pad.
    func wave(nx: Int, ny: Int) throws -> Waves {
        let grid = try Grid(nx: nx, ny: ny, extentX: 6, extentY: 8)
        var f = ComplexField(nx: nx, ny: ny)
        for ix in 0..<nx {
            for iy in 0..<ny {
                let u = 2 * Double.pi * Double(ix) / Double(nx), v = 2 * Double.pi * Double(iy) / Double(ny)
                f.real[ix * ny + iy] = 1 + 0.3 * cos(u + 2 * v) + 0.1 * sin(3 * u)
                f.imaginary[ix * ny + iy] = 0.2 * sin(2 * u - v)
            }
        }
        return Waves(field: f, grid: grid, beam: try ElectronBeam(energy: 200e3), normalization: .reciprocalSpace)
    }

    func maxDiff(_ a: ComplexField, _ b: ComplexField) -> Double {
        zip(zip(a.real, b.real), zip(a.imaginary, b.imaginary)).map { hypot($0.0 - $0.1, $1.0 - $1.1) }.max() ?? 0
    }

    @Test func padThenCropIsTheIdentity() throws {
        for (nx, ny) in [(16, 20), (15, 21), (16, 21)] {
            let w = try wave(nx: nx, ny: ny)
            let padded = try w.resampled(rows: nx + 7, columns: ny + 10)
            #expect(padded.grid.nx == nx + 7 && padded.grid.extentY == 8 && padded.beam == w.beam)
            let back = try padded.resampled(rows: nx, columns: ny)
            #expect(maxDiff(back.field, w.field) < 1e-12, "\(nx) x \(ny)")
        }
    }

    @Test func aBandlimitedWaveKeepsItsValuesWhenCropped() throws {
        let big = try wave(nx: 24, ny: 30)
        let small = try big.resampled(rows: 9, columns: 11)
        let direct = try wave(nx: 9, ny: 11)
        #expect(maxDiff(small.field, direct.field) < 1e-12, "amplitude is kept (abTEM's 'values' normalisation)")
    }

    @Test func aConstantStaysConstantAndMixedAxesWork() throws {
        let plane = Waves.planeWave(beam: try ElectronBeam(energy: 200e3), grid: try Grid(nx: 12, ny: 7, extentX: 3, extentY: 3))
        let r = try plane.resampled(rows: 5, columns: 18)
        #expect(r.field.real.allSatisfy { abs($0 - 1) < 1e-13 } && r.field.imaginary.allSatisfy { abs($0) < 1e-13 })
        #expect(r.normalization == .values)
        #expect(throws: SimulatoError.unsupportedGridLength(1)) { try plane.resampled(rows: 1, columns: 4) }
    }

    @Test func frequencyMapsFollowAbTEMsMasks() {
        // crop 8 -> 5: first 3 (0, 1, 2), last 2 (6, 7); pad 5 -> 8: 0, 1, 2 stay, 3, 4 go to 6, 7; even 4 -> 6: 0, 1 and 4, 5.
        #expect(Waves.resampleMap(from: 8, to: 5).map { [$0.0, $0.1] } == [[0, 0], [1, 1], [2, 2], [6, 3], [7, 4]])
        #expect(Waves.resampleMap(from: 5, to: 8).map { [$0.0, $0.1] } == [[0, 0], [1, 1], [2, 2], [3, 6], [4, 7]])
        #expect(Waves.resampleMap(from: 4, to: 6).map { [$0.0, $0.1] } == [[0, 0], [1, 1], [2, 4], [3, 5]])
        #expect(Waves.resampleMap(from: 6, to: 4).map { [$0.0, $0.1] } == [[0, 0], [1, 1], [4, 2], [5, 3]])
    }
}

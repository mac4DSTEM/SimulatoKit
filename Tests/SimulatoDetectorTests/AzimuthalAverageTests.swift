//
//  AzimuthalAverageTests.swift
//  Role: the azimuthal average of a centred diffraction pattern (abTEM's DiffractionPatterns.azimuthal_average), against
//        numbers abTEM 1.0.10 printed (tools/pins/azimuthal/pins.py) for planted 9 x 8 and 8 x 6 patterns, values
//        (k + 1)^1.3 with k the flat index, at 200 keV; and its refusals.
//

import Foundation
import Testing
import SimulatoCore
import SimulatoDetector

@Suite("Azimuthal average")
struct AzimuthalAverageTests {
    static func pattern(_ nx: Int, _ ny: Int, _ sx: Double, _ sy: Double) -> RealField {
        RealField(nx: nx, ny: ny, values: (0..<(nx * ny)).map { pow(Double($0) + 1, 1.3) }, samplingX: sx, samplingY: sy)
    }

    static let beam = try! ElectronBeam(energy: 200e3)

    /// assert_allclose defaults (rtol 1e-7, atol 0): abTEM states no tolerance for the azimuthal average.
    static func close(_ a: [Double], _ b: [Double], rtol: Double = 1e-7) -> Bool {
        a.count == b.count && zip(a, b).allSatisfy { abs($0 - $1) <= rtol * abs($1) }
    }

    @Test func stepRingsMatchAbTEMOnAnOddEvenPattern() throws {
        let p = Self.pattern(9, 8, 0.05, 0.06)
        let a = try p.azimuthalAverage(beam: Self.beam)
        #expect(a.count == 5 && a.ensembleShape == [])
        #expect(abs(a.sampling - 0.05000000000000001) <= 1e-7 * 0.05)
        #expect(Self.close(a.values, [109.31059072723754, 110.0733829078507, 111.05800473544934, 113.83482973524923,
                                      116.75029096542004]))
        let b = try p.azimuthalAverage(beam: Self.beam, maxAngle: 9, radialSampling: 0.7, width: 1.5)
        #expect(b.count == 11 && abs(b.sampling - 0.035) <= 1e-7 * 0.035)
        #expect(Self.close(b.values, [109.71729296135845, 110.60442943790507, 110.97479454458214, 112.91164715276787,
                                      115.2149506722666, 115.6100936043095, 116.37171872576145, 116.20544576168864,
                                      116.3659299808196, 111.20251111350319, 109.37593155122985]))
    }

    @Test func gaussianRingsMatchAbTEM() throws {
        let p = Self.pattern(9, 8, 0.05, 0.06)
        let a = try p.azimuthalAverage(beam: Self.beam, weighting: .gaussian)
        #expect(Self.close(a.values, [109.52747242355026, 110.19514603212777, 111.79352436215716, 114.5200943034101,
                                      117.77733984988313]))
        let b = try p.azimuthalAverage(beam: Self.beam, maxAngle: 6, radialSampling: 1.3, weighting: .gaussian, width: 0.8)
        #expect(b.count == 4 && abs(b.sampling - 0.06500000000000002) <= 1e-7 * 0.065)
        #expect(Self.close(b.values, [109.39182869735909, 110.29310407876528, 112.87926982867883, 117.66284648473913]))
    }

    @Test func evenPatternMatchesAbTEM() throws {
        let p = Self.pattern(8, 6, 0.04, 0.04)
        let a = try p.azimuthalAverage(beam: Self.beam)
        #expect(a.count == 4 && abs(a.sampling - 0.039999999999999994) <= 1e-7 * 0.04)
        #expect(Self.close(a.values, [76.36746409305131, 76.84996447773611, 77.6798227805934, 75.22480102040664]))
        #expect(Self.close(try p.azimuthalAverage(beam: Self.beam, width: 0.5).values,
                           [76.08612050106568, 76.6139088690763, 77.74902662196723, 77.65019303476393]))
    }

    @Test func aRingWithNoPixelIsZero() throws {
        // A 1 x 1 pattern: one pixel at the centre; rings beyond a step width of it hold nothing.
        let p = RealField(nx: 1, ny: 1, values: [3], samplingX: 0.05, samplingY: 0.05)
        let a = try p.azimuthalAverage(beam: Self.beam, maxAngle: 5, width: 0.5)
        #expect(a.count == 4 && a.values == [3, 0, 0, 0])
    }

    @Test func invalidArgumentsAreRefused() {
        let p = Self.pattern(9, 8, 0.05, 0.06)
        #expect(throws: SimulatoError.self) { try p.azimuthalAverage(beam: Self.beam, radialSampling: 0) }
        #expect(throws: SimulatoError.self) { try p.azimuthalAverage(beam: Self.beam, radialSampling: .nan) }
        #expect(throws: SimulatoError.self) { try p.azimuthalAverage(beam: Self.beam, width: 0) }
        #expect(throws: SimulatoError.self) { try p.azimuthalAverage(beam: Self.beam, width: -1) }
        #expect(throws: SimulatoError.self) { try p.azimuthalAverage(beam: Self.beam, maxAngle: 0) }
        #expect(throws: SimulatoError.self) { try p.azimuthalAverage(beam: Self.beam, maxAngle: .infinity) }
        #expect(throws: SimulatoError.self) { try p.azimuthalAverage(beam: Self.beam, maxAngle: 1e9) }
    }
}

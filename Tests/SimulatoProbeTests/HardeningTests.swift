//
//  HardeningTests.swift
//  Role: a line scan with an explicit huge `gpts` is refused, as the sampling path refuses it.
//

import Testing
import SimulatoCore
import SimulatoProbe

@Suite("Probe hardening")
struct ProbeHardeningTests {
    @Test func hugeGptsIsRefused() throws {
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(1, 0), gpts: 2_000_000_000) }
        #expect(throws: SimulatoError.self) { try LineScan(end: SIMD2(1, 0), gpts: .max) }
        #expect(try LineScan(end: SIMD2(1, 0), gpts: 5).gpts == 5)
        // a line of zero length has NaN direction and angle
        let z = try LineScan(start: SIMD2(1, 1), end: SIMD2(1, 1), gpts: 3)
        #expect(z.direction.x.isNaN && z.angle.isNaN)
    }
}

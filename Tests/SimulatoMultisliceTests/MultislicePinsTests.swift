//
//  MultislicePinsTests.swift
//  Role: registration sections 3.7, 3.8 and 3.10 against the pins, CPU path: the incident probe, the exit waves, the
//        diffraction patterns and the per-plane scalars of every case and beam, with abTEM's π (bounds 1e−10, the probe
//        1e−12) and with the exact π (1e−6).
//

import XCTest
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe

final class MultislicePinsTests: XCTestCase {
    static let runs: [(String, String)] = [("A1", "plane"), ("A2", "plane"), ("B1", "plane"), ("B2", "probe"),
                                           ("C1", "plane"), ("C1", "probe"), ("D1", "plane"), ("D1", "probe")]

    /// B2 uses B1's structure and potential.
    static func potentialCase(_ id: String) -> String { id == "B2" ? "B1" : id }

    func run(_ id: String, _ label: String, pi: LobatoParameterization.PiConstant,
             backend: Backend = .cpu) throws -> (waves: [Waves], beam: [String: Any]) {
        let c = Pins.case(id)
        let b = (c["beams"] as! [String: Any])[label] as! [String: Any]
        let potential = try Pins.potential(Self.potentialCase(id), pi: pi)
        let incident = try Pins.incident(id, label, grid: potential.grid)
        let waves = try Multislice.run(incident, through: potential,
                                       exitPlanes: .afterSlices(c["exit_planes"] as! [Int]), backend: backend)
        XCTAssertEqual(waves.count, (c["exit_planes"] as! [Int]).count, "\(id)/\(label) plane count")
        return (waves, b)
    }

    func testIncidentProbeEqualsAbTEM() throws {
        for (id, label) in Self.runs where label == "probe" {
            let potential = try Pins.potential(Self.potentialCase(id), pi: .exact)
            let probe = try Pins.incident(id, label, grid: potential.grid)
            XCTAssertLessThanOrEqual(relMax(probe.field, Pins.complex("\(id)/\(label).incident")), 1e-12, "\(id) probe")
        }
    }

    func testExitWavesPatternsAndScalarsEqualAbTEM() throws {
        for (id, label) in Self.runs {
            let (waves, b) = try run(id, label, pi: .abTEMFloat32)
            try compare(id, label, waves, b, bound: 1e-10)
        }
    }

    func testExactPiStaysWithinTheRegisteredDifference() throws {
        for (id, label) in Self.runs {
            let (waves, b) = try run(id, label, pi: .exact)
            try compare(id, label, waves, b, bound: 1e-6)
        }
    }

    func compare(_ id: String, _ label: String, _ waves: [Waves], _ b: [String: Any], bound: Double) throws {
        let totals = b["total_intensity"] as! [Double], sums = b["pattern_sum"] as! [Double]
        let maxima = b["exit_max_abs"] as! [Double], centres = b["pattern_center"] as! [Double]
        for (q, w) in waves.enumerated() {
            let pattern = try w.diffractionPattern(maxAngle: .cutoff)
            XCTAssertEqual([pattern.nx, pattern.ny], b["pattern_gpts"] as! [Int], "\(id)/\(label) pattern size")
            XCTAssertLessThanOrEqual(relDiff(w.totalIntensity, totals[q]), bound, "\(id)/\(label) plane \(q) Σ|ψ|²")
            XCTAssertLessThanOrEqual(relDiff(pattern.sum, sums[q]), bound, "\(id)/\(label) plane \(q) Σ pattern")
            XCTAssertLessThanOrEqual(relDiff(w.field.maxMagnitude, maxima[q]), bound, "\(id)/\(label) plane \(q) max|ψ|")
            XCTAssertLessThanOrEqual(abs(pattern[pattern.nx / 2, pattern.ny / 2] - centres[q]), bound * pattern.max,
                                     "\(id)/\(label) plane \(q) pattern centre")
            if (b["stored_planes"] as! [Int]).contains(q) {
                XCTAssertLessThanOrEqual(relMax(w.field, Pins.complex("\(id)/\(label).exit.\(q)")), bound,
                                         "\(id)/\(label) exit wave plane \(q)")
            }
            if (b["stored_pattern_planes"] as! [Int]).contains(q) {
                XCTAssertLessThanOrEqual(relMax(pattern.values, Pins.doubles("\(id)/\(label).pattern.\(q)")), bound,
                                         "\(id)/\(label) pattern plane \(q)")
            }
        }
    }
}

//
//  MultislicePinsTests.swift
//  Role: registration sections 3.7, 3.8 and 3.10 against the pins, CPU path: the incident probe, the exit waves, the
//        diffraction patterns and the per-plane scalars of every case and beam, with abTEM's π (bounds 1e−10, the probe
//        1e−12) and with the exact π (1e−6).
//

import Testing
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe

/// One pinned case and beam, e.g. "B2/probe".
struct CaseBeam: Sendable, CustomTestStringConvertible {
    let id: String
    let label: String
    var testDescription: String { "\(id)/\(label)" }
    /// B2 uses B1's structure and potential.
    var potentialCase: String { id == "B2" ? "B1" : id }

    static let all = [CaseBeam(id: "A1", label: "plane"), CaseBeam(id: "A2", label: "plane"),
                      CaseBeam(id: "B1", label: "plane"), CaseBeam(id: "B2", label: "probe"),
                      CaseBeam(id: "C1", label: "plane"), CaseBeam(id: "C1", label: "probe"),
                      CaseBeam(id: "D1", label: "plane"), CaseBeam(id: "D1", label: "probe")]
}

@Suite("Multislice against abTEM")
struct MultislicePinsTests {
    func run(_ r: CaseBeam, pi: LobatoParameterization.PiConstant) throws -> (waves: [Waves], beam: [String: Any]) {
        let c = Pins.case(r.id)
        let b = (c["beams"] as! [String: Any])[r.label] as! [String: Any]
        let potential = try Pins.potential(r.potentialCase, pi: pi)
        let incident = try Pins.incident(r.id, r.label, grid: potential.grid)
        let waves = try Multislice.run(incident, through: potential, exitPlanes: .afterSlices(c["exit_planes"] as! [Int]))
        #expect(waves.count == (c["exit_planes"] as! [Int]).count, "\(r.id)/\(r.label) plane count")
        return (waves, b)
    }

    @Test(arguments: CaseBeam.all.filter { $0.label == "probe" })
    func incidentProbeEqualsAbTEM(r: CaseBeam) throws {
        let potential = try Pins.potential(r.potentialCase, pi: .exact)
        let probe = try Pins.incident(r.id, r.label, grid: potential.grid)
        #expect(relMax(probe.field, Pins.complex("\(r.id)/\(r.label).incident")) <= 1e-12, "\(r.id) probe")
    }

    @Test(arguments: CaseBeam.all)
    func exitWavesPatternsAndScalarsEqualAbTEM(r: CaseBeam) throws {
        let (waves, b) = try run(r, pi: .abTEMFloat32)
        try compare(r, waves, b, bound: 1e-10)
    }

    @Test(arguments: CaseBeam.all)
    func exactPiStaysWithinTheRegisteredDifference(r: CaseBeam) throws {
        let (waves, b) = try run(r, pi: .exact)
        try compare(r, waves, b, bound: 1e-6)
    }

    func compare(_ r: CaseBeam, _ waves: [Waves], _ b: [String: Any], bound: Double) throws {
        let totals = b["total_intensity"] as! [Double], sums = b["pattern_sum"] as! [Double]
        let maxima = b["exit_max_abs"] as! [Double], centres = b["pattern_center"] as! [Double]
        let tag = "\(r.id)/\(r.label)"
        for (q, w) in waves.enumerated() {
            let pattern = try w.diffractionPattern(maxAngle: .cutoff)
            #expect([pattern.nx, pattern.ny] == b["pattern_gpts"] as! [Int], "\(tag) pattern size")
            #expect(relDiff(w.totalIntensity, totals[q]) <= bound, "\(tag) plane \(q) Σ|ψ|²")
            #expect(relDiff(pattern.sum, sums[q]) <= bound, "\(tag) plane \(q) Σ pattern")
            #expect(relDiff(w.field.maxMagnitude, maxima[q]) <= bound, "\(tag) plane \(q) max|ψ|")
            #expect(abs(pattern[pattern.nx / 2, pattern.ny / 2] - centres[q]) <= bound * pattern.max,
                    "\(tag) plane \(q) pattern centre")
            if (b["stored_planes"] as! [Int]).contains(q) {
                #expect(relMax(w.field, Pins.complex("\(tag).exit.\(q)")) <= bound, "\(tag) exit wave plane \(q)")
            }
            if (b["stored_pattern_planes"] as! [Int]).contains(q) {
                #expect(relMax(pattern.values, Pins.doubles("\(tag).pattern.\(q)")) <= bound, "\(tag) pattern plane \(q)")
            }
        }
    }
}

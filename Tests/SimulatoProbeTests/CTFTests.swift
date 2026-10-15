//
//  CTFTests.swift
//  Role: registration-slice-3.md: the CTF kernels, the HRTEM images of a SrTiO₃ exit wave and the closed forms against
//        abTEM (section 4), and the truth of section 5 (the identity CTF, the envelopes' closed forms, a flat vacuum image,
//        weak-phase contrast at Scherzer, the textbook values).
//

import CryptoKit
import Foundation
import Testing
import SimulatoCore
import SimulatoMultislice
import SimulatoProbe

enum Pins3 {
    static let directory = URL(fileURLWithPath: #filePath).deletingLastPathComponent().appendingPathComponent("Fixtures")

    struct Loaded: @unchecked Sendable {  // immutable after load
        let json: [String: Any]
        let blob: Data
    }

    static let loaded: Loaded = {
        let json = try! JSONSerialization.jsonObject(
            with: Data(contentsOf: directory.appendingPathComponent("slice-3-pins.json"))) as! [String: Any]
        let arrays = json["arrays"] as! [String: Any]
        let blob = try! Data(contentsOf: directory.appendingPathComponent(arrays["file"] as! String))
        let digest = SHA256.hash(data: blob).map { String(format: "%02x", $0) }.joined()
        precondition(digest == arrays["sha256"] as! String, "slice-3-arrays.bin does not match its sha256 in the pins")
        return Loaded(json: json, blob: blob)
    }()

    static var json: [String: Any] { loaded.json }

    static func doubles(_ name: String) -> [Double] {
        let index = (json["arrays"] as! [String: Any])["index"] as! [String: Any]
        let entry = index[name] as! [String: Any]
        let offset = entry["offset"] as! Int
        let count = (entry["count"] as! Int) * ((entry["kind"] as! String) == "c16" ? 2 : 1)
        return loaded.blob.withUnsafeBytes { raw in
            (0..<count).map { i in
                Double(bitPattern: UInt64(littleEndian: raw.loadUnaligned(fromByteOffset: offset + 8 * i, as: UInt64.self)))
            }
        }
    }

    static func `case`(_ id: String) -> [String: Any] { (json["cases"] as! [String: Any])[id] as! [String: Any] }

    static func ctf(_ id: String) throws -> ContrastTransferFunction {
        let c = Self.case(id)
        var a = Aberrations()
        for (k, v) in c["aberrations"] as! [String: Double] {
            let path: WritableKeyPath<Aberrations, Double> = [
                "C10": \.C10, "C12": \.C12, "phi12": \.phi12, "C21": \.C21, "phi21": \.phi21, "C23": \.C23, "phi23": \.phi23,
                "C30": \.C30, "C32": \.C32, "phi32": \.phi32, "C34": \.C34, "phi34": \.phi34, "C41": \.C41, "phi41": \.phi41,
                "C43": \.C43, "phi43": \.phi43, "C45": \.C45, "phi45": \.phi45, "C50": \.C50, "C52": \.C52, "phi52": \.phi52,
                "C54": \.C54, "phi54": \.phi54, "C56": \.C56, "phi56": \.phi56,
            ][k]!
            a[keyPath: path] = v
        }
        return try ContrastTransferFunction(aberrations: a, semiangleCutoff: (c["cutoff_mrad"] as? Double) ?? .infinity,
                                            softAperture: c["soft"] as! Bool, focalSpread: c["focal_spread"] as! Double,
                                            angularSpread: c["angular_spread"] as! Double)
    }

    static func grid(_ id: String) throws -> Grid {
        let c = Self.case(id)
        let g = c["gpts"] as! [Int], e = c["extent"] as! [Double]
        return try Grid(nx: g[0], ny: g[1], extentX: e[0], extentY: e[1])
    }

    /// The SrTiO₃ [001] exit wave of H1–H3 (3 × 3 × 5 cells, 120², 300 keV), by the port.
    static func strontiumTitanateExit(pi: LobatoParameterization.PiConstant) throws -> Waves {
        let a0 = 3.905
        let cell = try AtomicStructure(box: SIMD3(a0, a0, a0), atoms: [
            .init(atomicNumber: 38, position: SIMD3(0, 0, 0)), .init(atomicNumber: 22, position: SIMD3(0.5, 0.5, 0.5) * a0),
            .init(atomicNumber: 8, position: SIMD3(0.5, 0.5, 0) * a0), .init(atomicNumber: 8, position: SIMD3(0.5, 0, 0.5) * a0),
            .init(atomicNumber: 8, position: SIMD3(0, 0.5, 0.5) * a0),
        ])
        let p = try ProjectedPotential(structure: cell.repeated(3, 3, 5), nx: 120, ny: 120, sliceThickness: 1.9525,
                                       parameterization: LobatoParameterization(piConstant: pi))
        return try Multislice.run(.planeWave(beam: ElectronBeam(energy: 300e3), grid: p.grid), through: p)[0]
    }
}

func relMax3(_ a: [Double], _ b: [Double]) -> Double {
    guard a.count == b.count else { return .infinity }
    var d = 0.0, m = 0.0
    for i in a.indices { d = max(d, abs(a[i] - b[i])); m = max(m, abs(b[i])) }
    return m > 0 ? d / m : d
}

func interleaved(_ f: ComplexField) -> [Double] { zip(f.real, f.imaginary).flatMap { [$0.0, $0.1] } }

@Suite("CTF against abTEM")
struct CTFTests {
    @Test(arguments: ["H1", "H2", "H3", "K1", "K2"])
    func kernelsEqualAbTEM(id: String) throws {
        let c = Pins3.case(id)
        let beam = try ElectronBeam(energy: c["energy"] as! Double)
        let h = try Pins3.ctf(id).kernel(on: Pins3.grid(id), wavelength: beam.wavelength)
        #expect(relMax3(interleaved(h), Pins3.doubles("\(id).kernel")) <= 1e-12, "\(id) kernel")
    }

    @Test(arguments: [LobatoParameterization.PiConstant.abTEMFloat32, .exact])
    func hrtemImagesEqualAbTEM(pi: LobatoParameterization.PiConstant) throws {
        let exit = try Pins3.strontiumTitanateExit(pi: pi)
        for id in ["H1", "H2", "H3"] {
            let image = try exit.applying(Pins3.ctf(id)).intensity().values
            #expect(relMax3(image, Pins3.doubles("\(id).image")) <= (pi == .abTEMFloat32 ? 1e-10 : 1e-6), "\(id) image")
        }
    }

    @Test func closedFormsEqualAbTEMAndTheTextbook() throws {
        for row in Pins3.json["closed_forms"] as! [[String: Double]] {
            let s = try ContrastTransferFunction.scherzerDefocus(Cs: row["Cs"]!, energy: row["energy"]!)
            let r = try ContrastTransferFunction.pointResolution(Cs: row["Cs"]!, energy: row["energy"]!)
            #expect(abs(s - row["scherzer"]!) <= 1e-15 * abs(row["scherzer"]!), "Scherzer \(row)")
            #expect(abs(r - row["point_resolution"]!) <= 1e-15 * row["point_resolution"]!, "point resolution \(row)")
        }
        for row in Pins3.json["nyquist"] as! [[String: Double]] {
            let n = try ContrastTransferFunction.nyquistSampling(semiangleCutoff: row["semiangle"]!, energy: row["energy"]!)
            #expect(abs(n - row["sampling"]!) <= 1e-15 * row["sampling"]!, "Nyquist \(row)")
        }
        // Truth 5: the textbook values, four figures.
        #expect(abs(try ContrastTransferFunction.scherzerDefocus(Cs: 1.2e7, energy: 200e3) - 671.9) < 0.05)
        #expect(abs(try ContrastTransferFunction.pointResolution(Cs: 1.2e7, energy: 200e3) - 2.370) < 0.0005)
        #expect(abs(try ContrastTransferFunction.nyquistSampling(semiangleCutoff: 20, energy: 200e3) - 0.3135) < 0.00005)
    }
}

@Suite("CTF truth")
struct CTFTruthTests {
    /// 1: no aberrations, no envelopes, no aperture: the image is |ψ|².
    @Test func identityCTFLeavesTheImage() throws {
        let exit = try Pins3.strontiumTitanateExit(pi: .exact)
        let image = try exit.applying(ContrastTransferFunction()).intensity().values
        #expect(relMax3(image, exit.intensity().values) <= 1e-14)
    }

    /// 2: the envelopes' closed forms for round aberrations.
    @Test(arguments: [0.005, 0.010, 0.020, 0.040])
    func envelopesEqualTheirClosedForms(alpha: Double) throws {
        let lambda = try ElectronBeam(energy: 200e3).wavelength
        var a = Aberrations()
        a.C10 = -400; a.C30 = 1.2e7; a.C50 = 3e9
        let beta = 0.6, delta = 35.0
        let hand = exp(-pow(beta / 1e3 / 2, 2) * pow(2 * Double.pi / lambda * (a.C10 * alpha + a.C30 * pow(alpha, 3)
                                                                            + a.C50 * pow(alpha, 5)), 2))
        let es = ContrastTransferFunction.spatialEnvelope(alpha: alpha, phi: 0.7, wavelength: lambda, aberrations: a,
                                                          angularSpread: beta)
        #expect(abs(es - hand) <= 1e-14, "E_s at \(alpha)")
        #expect(a.dChiDPhi(alpha: alpha, phi: 0.7, wavelength: lambda) == 0, "round aberrations: no azimuthal term")
        let et = ContrastTransferFunction.temporalEnvelope(alpha: alpha, wavelength: lambda, focalSpread: delta)
        #expect(abs(et - exp(-pow(Double.pi * delta * alpha * alpha / (2 * lambda), 2))) <= 1e-14, "E_t at \(alpha)")
        #expect(et >= 0 && et <= 1 && es >= 0 && es <= 1, "in [0, 1] (amendment 2: the envelope underflows to 0 at large α)")
        let edge = (2 * lambda / (Double.pi * delta)).squareRoot()  // π Δ α² / (2λ) = 1
        #expect(abs(ContrastTransferFunction.temporalEnvelope(alpha: edge, wavelength: lambda, focalSpread: delta) - exp(-1))
                <= 1e-14)
        #expect(ContrastTransferFunction.temporalEnvelope(alpha: 0, wavelength: lambda, focalSpread: delta) == 1)
    }

    /// 3: vacuum is flat under any CTF.
    @Test func vacuumImageIsFlat() throws {
        let vac = try AtomicStructure(box: SIMD3(9.6, 14.4, 4), atoms: [])
        let p = try ProjectedPotential(structure: vac, nx: 96, ny: 160, sliceThickness: 2)
        let exit = try Multislice.run(.planeWave(beam: ElectronBeam(energy: 80e3), grid: p.grid), through: p)[0]
        var a = Aberrations()
        a.C10 = 35; a.C12 = 12; a.phi12 = 0.3; a.C23 = 90; a.phi23 = 1.9; a.C30 = -4e5
        let ctf = try ContrastTransferFunction(aberrations: a, semiangleCutoff: 35, focalSpread: 25, angularSpread: 0.8)
        let image = try exit.applying(ctf).intensity().values
        #expect(image.map { abs($0 - 1) }.max()! <= 1e-14)
    }

    /// 4: weak-phase contrast at Scherzer: the atom is dark and the image follows 1 + 2 IFFT[φ̂ sin χ A].
    @Test func weakPhaseContrastAtScherzer() throws {
        let beam = try ElectronBeam(energy: 300e3)
        let atom = try AtomicStructure(box: SIMD3(12.8, 12.8, 2), atoms: [.init(atomicNumber: 6, position: SIMD3(6.4, 6.4, 1))])
        let p = try ProjectedPotential(structure: atom, nx: 128, ny: 128, sliceThickness: 2)
        let exit = try Multislice.run(.planeWave(beam: beam, grid: p.grid), through: p)[0]
        var a = Aberrations()
        a.Cs = 1.2e7
        a.defocus = try ContrastTransferFunction.scherzerDefocus(Cs: 1.2e7, energy: 300e3)
        let cutoff = beam.wavelength / (try ContrastTransferFunction.pointResolution(Cs: 1.2e7, energy: 300e3)) * 1e3
        let ctf = try ContrastTransferFunction(aberrations: a, semiangleCutoff: cutoff)
        let image = try exit.applying(ctf).intensity()
        // (a) the atom is dark: the minimum is at the atom's pixel.
        let minIndex = image.values.indices.min { image.values[$0] < image.values[$1] }!
        #expect(minIndex == 64 * 128 + 64, "the image's minimum is at the atom")
        // (b) the weak-phase prediction: φ = σ V (bandlimited as the transmission), I − 1 ≈ 2 IFFT[φ̂ sin χ A].
        let sigma = beam.interactionParameter
        let v = p.slices[0].values
        let phase = ComplexField(nx: 128, ny: 128, real: v.map { sigma * $0 }, imaginary: [Double](repeating: 0, count: v.count))
        var kernel = try ctf.kernel(on: p.grid, wavelength: beam.wavelength)  // exp(−iχ) A; sin χ A = −Im(kernel)
        for i in kernel.real.indices { kernel.real[i] = -kernel.imaginary[i]; kernel.imaginary[i] = 0 }
        let predicted = try Self.filter(phase, kernel, grid: p.grid).real.map { 2 * $0 }
        let measured = image.values.map { $0 - 1 }
        let scale = measured.map(abs).max()!
        let deviation = zip(predicted, measured).map { abs($0 - $1) }.max()! / scale
        #expect(deviation <= 0.1, "weak-phase prediction within 0.1 of max |I − 1| (got \(deviation))")
        print("weak-phase deviation \(deviation)")
    }

    /// IFFT(FFT(f) · h), with the test-only DFT below (the kit's FFT is package-internal).
    static func filter(_ f: ComplexField, _ h: ComplexField, grid: Grid) throws -> ComplexField {
        let spectrum = try Self.fft(f, grid: grid, inverse: false)
        var product = spectrum
        for i in product.real.indices {
            product.real[i] = spectrum.real[i] * h.real[i] - spectrum.imaginary[i] * h.imaginary[i]
            product.imaginary[i] = spectrum.real[i] * h.imaginary[i] + spectrum.imaginary[i] * h.real[i]
        }
        return try Self.fft(product, grid: grid, inverse: true)
    }

    /// A direct 2-D DFT by separable 1-D DFTs (test-only; 128² is small).
    static func fft(_ f: ComplexField, grid: Grid, inverse: Bool) throws -> ComplexField {
        let nx = grid.nx, ny = grid.ny, s = inverse ? 1.0 : -1.0
        func dft(_ re: [Double], _ im: [Double]) -> ([Double], [Double]) {
            let n = re.count
            var r = [Double](repeating: 0, count: n), m = r
            for k in 0..<n {
                for j in 0..<n {
                    let a = s * 2 * Double.pi * Double((k * j) % n) / Double(n)
                    r[k] += re[j] * cos(a) - im[j] * sin(a)
                    m[k] += re[j] * sin(a) + im[j] * cos(a)
                }
            }
            return (r, m)
        }
        var out = f
        for ix in 0..<nx {
            let (r, m) = dft(Array(out.real[(ix * ny)..<((ix + 1) * ny)]), Array(out.imaginary[(ix * ny)..<((ix + 1) * ny)]))
            for iy in 0..<ny { out.real[ix * ny + iy] = r[iy]; out.imaginary[ix * ny + iy] = m[iy] }
        }
        for iy in 0..<ny {
            let (r, m) = dft((0..<nx).map { out.real[$0 * ny + iy] }, (0..<nx).map { out.imaginary[$0 * ny + iy] })
            for ix in 0..<nx { out.real[ix * ny + iy] = r[ix]; out.imaginary[ix * ny + iy] = m[ix] }
        }
        if inverse {
            let k = 1 / Double(nx * ny)
            out.real = out.real.map { $0 * k }
            out.imaginary = out.imaginary.map { $0 * k }
        }
        return out
    }

    @Test func refusals() {
        #expect(throws: (any Error).self) { try ContrastTransferFunction(focalSpread: -1) }
        #expect(throws: (any Error).self) { try ContrastTransferFunction(angularSpread: -0.4) }
        #expect(throws: (any Error).self) { try ContrastTransferFunction(semiangleCutoff: 0) }
        #expect(throws: (any Error).self) { try ContrastTransferFunction(focalSpread: .nan) }
        #expect(throws: (any Error).self) { try ContrastTransferFunction.nyquistSampling(semiangleCutoff: 0, energy: 200e3) }
        #expect(throws: (any Error).self) { try ContrastTransferFunction.nyquistSampling(semiangleCutoff: -5, energy: 200e3) }
    }
}

@Suite("CTF conventions as parameters (registration-generality.md)")
struct CTFConventionTests {
    @Test func theClosedFormsTakeTheirConstants() throws {
        let lambda = try ElectronBeam(energy: 300e3).wavelength
        let cs = 1.2e7
        #expect(try ContrastTransferFunction.scherzerDefocus(Cs: cs, energy: 300e3, factor: 4.0 / 3)
                == (4.0 / 3 * cs * lambda).squareRoot())
        #expect(try ContrastTransferFunction.scherzerDefocus(Cs: -cs, energy: 300e3, factor: 1)
                == -(cs * lambda).squareRoot())
        #expect(try ContrastTransferFunction.pointResolution(Cs: cs, energy: 300e3, divisor: 1)
                == pow(pow(lambda, 3) * cs, 0.25))
        #expect(try ContrastTransferFunction.nyquistSampling(semiangleCutoff: 20, energy: 300e3, oversampling: 2)
                == 1 / (2 * 20 / lambda * 1e-3))
        #expect(throws: (any Error).self) { try ContrastTransferFunction.scherzerDefocus(Cs: cs, energy: 300e3, factor: 0) }
        #expect(throws: (any Error).self) { try ContrastTransferFunction.pointResolution(Cs: cs, energy: 300e3, divisor: -1) }
        #expect(throws: (any Error).self) {
            try ContrastTransferFunction.nyquistSampling(semiangleCutoff: 20, energy: 300e3, oversampling: .nan)
        }
    }
}

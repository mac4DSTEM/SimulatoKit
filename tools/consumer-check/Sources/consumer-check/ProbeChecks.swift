//
//  ProbeChecks.swift
//  Role: SimulatoProbe from outside: the aperture, the aberrations and their derivatives, the probe, the scan grid, the
//        contrast transfer function with its closed forms, and the CTF applied to waves.
//

import Foundation
import SimulatoCore
import SimulatoProbe

func checkProbe(_ r: inout Report) {
    r.section("Probe: aperture") { r in
        let g = try Grid(nx: 32, ny: 32, extentX: 8, extentY: 8)
        r.refuses("Aperture(0)", .invalidAperture(0)) { try Aperture(semiangleCutoff: 0) }
        r.refusesCase("Aperture(NaN)", "invalidAperture") { try Aperture(semiangleCutoff: .nan) }
        let none = try Aperture(semiangleCutoff: .infinity)
        r.check("Aperture: .infinity is allowed and is 1 everywhere", none.values(on: g, wavelength: lambda200).allSatisfy { $0 == 1 })
        let hard = try Aperture(semiangleCutoff: 20, soft: false)
        let softDefault = try Aperture(semiangleCutoff: 20)
        r.check("Aperture stores semiangleCutoff and soft; default soft is true", hard.semiangleCutoff == 20 && !hard.soft
                && softDefault.soft)
        // hard edge: 1 where α = |k| λ ≤ 20 mrad; counted by hand on the fftfreq grid
        let values = hard.values(on: g, wavelength: lambda200)
        var count = 0
        let k = Grid.fftfreq(32, 0.25)
        for i in 0..<32 { for j in 0..<32 where (k[i] * k[i] + k[j] * k[j]).squareRoot() * lambda200 <= 0.020 { count += 1 } }
        r.equal("hard aperture: pixels at 1 = pixels with α ≤ 20 mrad", Int(values.reduce(0, +)), count)
        r.check("hard aperture is 0 or 1", values.allSatisfy { $0 == 0 || $0 == 1 })
        let soft = try Aperture(semiangleCutoff: 20).values(on: g, wavelength: lambda200)
        r.check("soft aperture: 1 at k = 0, values in [0, 1], at least one edge pixel strictly between",
                soft[0] == 1 && soft.allSatisfy { $0 >= 0 && $0 <= 1 } && soft.contains { $0 > 0 && $0 < 1 })
        r.check("Aperture is Hashable", hard == hard && hard != softDefault)
    }

    r.section("Probe: aberrations") { r in
        var ab = Aberrations()
        r.check("Aberrations(): every coefficient is zero", ab.C10 == 0 && ab.C30 == 0 && ab.C56 == 0 && ab.phi12 == 0
                && ab.defocus == 0 && ab.Cs == 0 && ab == Aberrations())
        ab.defocus = 50
        r.equal("defocus = −C10", ab.C10, -50)
        ab.Cs = 1e7
        r.equal("Cs = C30", ab.C30, 1e7)
        let lam = lambda200, alpha = 0.01
        // χ = (2π/λ) Σ α^{n+1}/(n+1) C_n cos(...): defocus Δf: −π Δf α² / λ ... (α²/2) C10 · 2π/λ
        var def = Aberrations(); def.C10 = -100
        r.near("chi: pure C10 = (2π/λ)(α²/2) C10", def.chi(alpha: alpha, phi: 0.3, wavelength: lam),
               2 * Double.pi / lam * alpha * alpha / 2 * -100, 1e-9)
        var cs = Aberrations(); cs.C30 = 1e7
        r.near("chi: pure C30 = (2π/λ)(α⁴/4) C30", cs.chi(alpha: alpha, phi: 1.1, wavelength: lam),
               2 * Double.pi / lam * pow(alpha, 4) / 4 * 1e7, 1e-9)
        var ast = Aberrations(); ast.C12 = 30; ast.phi12 = 0.4
        r.near("chi: C12 = (2π/λ)(α²/2) C12 cos(2(φ − φ12))", ast.chi(alpha: alpha, phi: 0.9, wavelength: lam),
               2 * Double.pi / lam * alpha * alpha / 2 * 30 * cos(2 * (0.9 - 0.4)), 1e-9)
        var high = Aberrations(); high.C56 = 1e9; high.phi56 = 0.2; high.C41 = 1e8; high.C23 = 1e5; high.phi23 = 0.1; high.C45 = 1e8
        high.C34 = 1e6; high.phi34 = 0.3; high.C50 = 1e10; high.C32 = 1e6; high.C21 = 1e4; high.C43 = 1e7; high.C52 = 1e8; high.C54 = 1e9
        r.check("chi: all orders finite and nonzero", high.chi(alpha: alpha, phi: 0.5, wavelength: lam).isFinite
                && high.chi(alpha: alpha, phi: 0.5, wavelength: lam) != 0)
        r.near("dChiDAlpha: pure C10 = (2π/λ) C10 α", def.dChiDAlpha(alpha: alpha, phi: 0.3, wavelength: lam),
               2 * Double.pi / lam * -100 * alpha, 1e-9)
        r.near("dChiDAlpha: pure C30 = (2π/λ) C30 α³", cs.dChiDAlpha(alpha: alpha, phi: 0.3, wavelength: lam),
               2 * Double.pi / lam * 1e7 * pow(alpha, 3), 1e-9)
        r.equal("dChiDPhi: pure C10 has no azimuthal dependence", def.dChiDPhi(alpha: alpha, phi: 0.3, wavelength: lam), 0)
        // against central differences of chi (h = 1e-6): ∂χ/∂α, and ∂χ/∂φ = α · dChiDPhi
        var mix = high; mix.C10 = -100; mix.C30 = 1e7; mix.C12 = 30; mix.phi12 = 0.4; mix.C32 = 1e6; mix.phi32 = 0.2
        let h = 1e-6, a0 = 0.012, p0 = 0.7
        let da = (mix.chi(alpha: a0 + h, phi: p0, wavelength: lam) - mix.chi(alpha: a0 - h, phi: p0, wavelength: lam)) / (2 * h)
        let dp = (mix.chi(alpha: a0, phi: p0 + h, wavelength: lam) - mix.chi(alpha: a0, phi: p0 - h, wavelength: lam)) / (2 * h)
        let analyticA = mix.dChiDAlpha(alpha: a0, phi: p0, wavelength: lam)
        let analyticP = mix.dChiDPhi(alpha: a0, phi: p0, wavelength: lam) * a0
        r.check("dChiDAlpha equals the central difference of chi (rel 1e-5)", abs(da - analyticA) <= 1e-5 * abs(da), "\(da) vs \(analyticA)")
        r.check("dChiDPhi · α equals the central difference of chi in φ (rel 1e-5)", abs(dp - analyticP) <= 1e-5 * abs(dp), "\(dp) vs \(analyticP)")
        r.check("Aberrations is Hashable", ab == ab && ab != Aberrations())
    }

    r.section("Probe: probe") { r in
        let beam = try ElectronBeam(energy: 200_000)
        let g = try Grid(nx: 32, ny: 32, extentX: 8, extentY: 8)
        let ap = try Aperture(semiangleCutoff: 20)
        let probe = try Probe(beam: beam, grid: g, aperture: ap)
        r.check("Probe stores beam, grid, aperture, aberrations", probe.beam == beam && probe.grid == g && probe.aperture == ap
                && probe.aberrations == Aberrations())
        let w = try probe.build(at: SIMD2(4, 4))
        r.check("build: .reciprocalSpace normalisation on the same grid", w.normalization == .reciprocalSpace && w.grid == g && w.beam == beam)
        // Σ_k |Ψ(k)|² = 1: the pattern of a .reciprocalSpace wave is |FFT ψ|² undivided, summing to 1; Σ_r |ψ|² = 1 / (nx ny).
        r.near("Σ_k |Ψ(k)|² = 1 (diffraction pattern, whole grid)", try w.diffractionPattern(maxAngle: .full).sum, 1, 1e-12)
        r.near("Σ_r |ψ|² = 1 / (nx ny)", w.totalIntensity, 1.0 / 1024, 1e-15)
        // centred on pixel (16, 16): 4 Å / 0.25 Å per pixel
        let intensity = w.intensity()
        let peak = intensity.values.indices.max { intensity.values[$0] < intensity.values[$1] }!
        r.equal("the probe peaks at the pixel of its position (16, 16)", peak, 16 * 32 + 16)
        let shifted = try probe.build(at: SIMD2(2, 6))
        let sh = shifted.intensity()
        r.equal("a probe at (2, 6) Å peaks at pixel (8, 24)", sh.values.indices.max { sh.values[$0] < sh.values[$1] }!, 8 * 32 + 24)
        var ab = Aberrations(); ab.defocus = 100
        let defocused = try Probe(beam: beam, grid: g, aperture: ap, aberrations: ab).build(at: SIMD2(4, 4))
        r.near("a defocused probe is normalised the same way", try defocused.diffractionPattern(maxAngle: .full).sum, 1, 1e-12)
        r.check("defocus lowers the peak intensity", defocused.intensity().max < intensity.max)
        var bad = Aberrations(); bad.C30 = .nan
        r.refusesCase("Probe with a NaN aberration", "invalidAberration") { try Probe(beam: beam, grid: g, aperture: ap, aberrations: bad) }
        r.refuses("build(at:) with a NaN position", .invalidPosition(index: 0)) { try probe.build(at: SIMD2(.nan, 1)) }
        let narrower = try Probe(beam: beam, grid: g, aperture: Aperture(semiangleCutoff: 10))
        r.check("Probe is Hashable", probe == probe && probe != narrower)
    }

    r.section("Probe: scan") { r in
        // linspace(0, 2, 2, endpoint: false) = [0, 1]; (0, 4, 2) = [0, 2]; index ix * ny + iy
        let s = try GridScan(start: SIMD2(0, 0), end: SIMD2(2, 4), nx: 2, ny: 2)
        r.check("GridScan stores start, end, gpts, endpoint (default false)", s.start == SIMD2(0, 0) && s.end == SIMD2(2, 4)
                && s.nx == 2 && s.ny == 2 && !s.endpoint)
        r.equal("positions, endpoint false", s.positions, [SIMD2(0, 0), SIMD2(0, 2), SIMD2(1, 0), SIMD2(1, 2)])
        let e = try GridScan(start: SIMD2(0, 0), end: SIMD2(2, 4), nx: 2, ny: 3, endpoint: true)
        r.equal("positions, endpoint true (3 along y: 0, 2, 4)", e.positions,
                [SIMD2(0, 0), SIMD2(0, 2), SIMD2(0, 4), SIMD2(2, 0), SIMD2(2, 2), SIMD2(2, 4)])
        r.equal("linspace(0, 1, 5, endpoint: true)", GridScan.linspace(0, 1, 5, endpoint: true), [0, 0.25, 0.5, 0.75, 1])
        r.equal("linspace(0, 1, 4, endpoint: false)", GridScan.linspace(0, 1, 4, endpoint: false), [0, 0.25, 0.5, 0.75])
        r.equal("linspace with one point and an endpoint = [start]", GridScan.linspace(3, 9, 1, endpoint: true), [3])
        let line = try GridScan(start: SIMD2(1, 0), end: SIMD2(1, 3), nx: 1, ny: 3)
        r.equal("a scan with extent on one axis only is allowed", line.positions, [SIMD2(1, 0), SIMD2(1, 1), SIMD2(1, 2)])
        r.check("GridScan is Hashable", s == s && s != e)
        r.refusesCase("GridScan with 0 points", "invalidScan") { try GridScan(start: .zero, end: SIMD2(1, 1), nx: 0, ny: 2) }
        r.refusesCase("GridScan with a NaN start", "invalidScan") { try GridScan(start: SIMD2(.nan, 0), end: SIMD2(1, 1), nx: 2, ny: 2) }
        r.refusesCase("GridScan without extent", "invalidScan") { try GridScan(start: SIMD2(1, 1), end: SIMD2(1, 1), nx: 2, ny: 2) }
    }

    r.section("Probe: contrast transfer function") { r in
        let cs = 1.2e7  // 1.2 mm in Å
        // Scherzer: sign(Cs) √(1.5 |Cs| λ)
        r.near("scherzerDefocus(Cs: 1.2 mm, 200 keV) = √(1.5 Cs λ)", try ContrastTransferFunction.scherzerDefocus(Cs: cs, energy: 200_000),
               (1.5 * cs * lambda200).squareRoot(), 1e-9)
        r.check("scherzerDefocus is 671.9 Å for Cs = 1.2 mm at 200 keV",
                abs(try ContrastTransferFunction.scherzerDefocus(Cs: cs, energy: 200_000) - 671.9) < 0.1)
        r.check("scherzerDefocus: negative Cs gives a negative defocus, Cs = 0 gives 0",
                try ContrastTransferFunction.scherzerDefocus(Cs: -cs, energy: 200_000) < 0
                && ContrastTransferFunction.scherzerDefocus(Cs: 0, energy: 200_000) == 0)
        r.near("pointResolution = (λ³ |Cs| / 6)^¼", try ContrastTransferFunction.pointResolution(Cs: cs, energy: 200_000),
               pow(pow(lambda200, 3) * cs / 6, 0.25), 1e-12)
        r.check("pointResolution(Cs: 1.2 mm, 200 keV) is 2.37 Å",
                abs(try ContrastTransferFunction.pointResolution(Cs: cs, energy: 200_000) - 2.37) < 0.01)
        r.near("nyquistSampling(30 mrad) = λ / (4 α)", try ContrastTransferFunction.nyquistSampling(semiangleCutoff: 30, energy: 200_000),
               lambda200 / (4 * 0.030), 1e-12)
        r.refusesCase("nyquistSampling(0 mrad)", "invalidAperture") { try ContrastTransferFunction.nyquistSampling(semiangleCutoff: 0, energy: 200_000) }
        r.refuses("scherzerDefocus at an invalid energy", .invalidEnergy(-5)) { try ContrastTransferFunction.scherzerDefocus(Cs: cs, energy: -5) }
        r.refusesCase("pointResolution at energy 0", "invalidEnergy") { try ContrastTransferFunction.pointResolution(Cs: cs, energy: 0) }

        // envelopes: Δ = 50 Å at α = 20 mrad: exp(−((π/(2λ)) Δ α²)²)
        let want = exp(-pow(0.5 * Double.pi / lambda200 * 50 * 0.02 * 0.02, 2))
        r.near("temporalEnvelope = exp(−((π/2λ) Δ α²)²)", ContrastTransferFunction.temporalEnvelope(alpha: 0.02, wavelength: lambda200, focalSpread: 50),
               want, 1e-14)
        r.equal("temporalEnvelope(α = 0) = 1", ContrastTransferFunction.temporalEnvelope(alpha: 0, wavelength: lambda200, focalSpread: 50), 1)
        var ab = Aberrations(); ab.C10 = -100
        let dchi = ab.dChiDAlpha(alpha: 0.02, phi: 0, wavelength: lambda200)
        r.near("spatialEnvelope = exp(−(β/2)² (∂χ/∂α)²) for a rotationally symmetric lens",
               ContrastTransferFunction.spatialEnvelope(alpha: 0.02, phi: 0, wavelength: lambda200, aberrations: ab, angularSpread: 0.5),
               exp(-pow(0.0005 / 2, 2) * dchi * dchi), 1e-12)
        r.equal("spatialEnvelope with β = 0 is 1", ContrastTransferFunction.spatialEnvelope(alpha: 0.02, phi: 0, wavelength: lambda200,
                                                                                              aberrations: ab, angularSpread: 0), 1)

        let plain = try ContrastTransferFunction()
        r.check("CTF defaults: no aberrations, no aperture (∞), soft, no spreads", plain.aberrations == Aberrations()
                && plain.semiangleCutoff == .infinity && plain.softAperture && plain.focalSpread == 0 && plain.angularSpread == 0)
        var ctf = try ContrastTransferFunction(aberrations: ab, semiangleCutoff: 30, softAperture: false, focalSpread: 40, angularSpread: 0.2)
        r.check("CTF stores its parameters; they are mutable; Hashable", ctf.semiangleCutoff == 30 && !ctf.softAperture
                && ctf.focalSpread == 40 && ctf.angularSpread == 0.2 && ctf != plain)
        ctf.focalSpread = 0
        let g = try Grid(nx: 32, ny: 32, extentX: 8, extentY: 8)
        let hk = try ctf.kernel(on: g, wavelength: lambda200)
        r.check("kernel: 1 + 0i at k = 0; zero outside the 30 mrad hard aperture; |H| ≤ 1",
                hk.real[0] == 1 && hk.imaginary[0] == 0 && hk[16, 16] == (0, 0)
                && hk.squaredMagnitude().allSatisfy { $0 <= 1 + 1e-12 })
        // a pure phase (no aperture, no envelope): |H| = 1 and H = exp(−iχ)
        let pure = try ContrastTransferFunction(aberrations: ab)
        let hp = try pure.kernel(on: g, wavelength: lambda200)
        r.check("kernel of a pure phase has |H| = 1", hp.squaredMagnitude().allSatisfy { abs($0 - 1) < 1e-12 })
        let k = Grid.fftfreq(32, 0.25)
        let alpha1 = (k[1] * k[1]).squareRoot() * lambda200
        let chi1 = ab.chi(alpha: alpha1, phi: 0, wavelength: lambda200)
        r.check("kernel at (1, 0) = exp(−iχ)", abs(hp.real[32] - cos(chi1)) < 1e-12 && abs(hp.imaginary[32] + sin(chi1)) < 1e-12)

        let beam = try ElectronBeam(energy: 200_000)
        let plane = Waves.planeWave(beam: beam, grid: g)
        let same = try plane.applying(plain)
        r.check("applying(default CTF) leaves a plane wave at 1", maxDifference(same.field.real, plane.field.real) < 1e-12
                && maxDifference(same.field.imaginary, plane.field.imaginary) < 1e-12)
        let probe = try Probe(beam: beam, grid: g, aperture: Aperture(semiangleCutoff: 20)).build(at: SIMD2(3, 5))
        let blurred = try probe.applying(pure)
        r.near("a pure-phase CTF conserves the total intensity (Parseval)", blurred.totalIntensity, probe.totalIntensity, 1e-15)
        r.check("applying keeps grid, beam, normalisation; the pure-phase CTF changes the wave",
                blurred.grid == g && blurred.beam == beam && blurred.normalization == .reciprocalSpace && blurred != probe)
        r.refusesCase("CTF with an infinite aberration", "invalidAberration") {
            var bad = Aberrations(); bad.C10 = .infinity
            return try ContrastTransferFunction(aberrations: bad)
        }
        r.refuses("CTF with cutoff 0", .invalidAperture(0)) { try ContrastTransferFunction(semiangleCutoff: 0) }
        r.refusesCase("CTF with a negative focal spread", "invalidAberration") { try ContrastTransferFunction(focalSpread: -1) }
        r.refusesCase("CTF with a NaN angular spread", "invalidAberration") { try ContrastTransferFunction(angularSpread: .nan) }
    }

    r.section("Probe: namespace") { r in
        r.check("SimulatoProbe namespace type exists", String(describing: SimulatoProbe.self) == "SimulatoProbe")
    }
}

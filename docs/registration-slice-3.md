# SimulatoKit slice 3: the contrast transfer function, partial coherence and HRTEM images (registration)

Registered 2026-10-15 (about 01:35 CEST), before any slice-3 Swift, by the SimulatoKit session under the owner's standing
instruction for the night (a kit that closes a slice before 05:00 with budget left takes its next slice). Tier 1: the
numbers are HRTEM image intensities, the transfer function a microscopist reads off (Scherzer defocus, point resolution) and
the envelopes that set the information limit. Slices 1 and 2 hold unchanged; this file adds to them.

## 1. What slice 3 makes

In SimulatoProbe:
- `ContrastTransferFunction` (abTEM's `CTF`): aberrations (slice 1's `Aberrations`), an objective aperture (semiangle cutoff,
  soft or hard; none by default), a focal spread Δ (Å) and an angular spread β (mrad); `kernel(on:wavelength:)` gives the
  complex transfer function on a grid.
- `SpatialEnvelope` and `TemporalEnvelope` values (also usable alone), and the closed forms `scherzerDefocus(Cs:energy:)`,
  `pointResolution(Cs:energy:)`, `nyquistSampling(semiangleCutoff:energy:)`.
- `Waves.applying(_ ctf:)` (abTEM's `apply_ctf`): the wave in Fourier space times the kernel, back to real space; its
  `intensity()` is the HRTEM image.

Not in slice 3: the CTF's Wiener filter and phase flipping, point-spread functions, profiles and plots, aperture shapes other
than the disk, distributions of parameters (focal series as ensembles; a caller loops), partial coherence beyond the
quasi-coherent envelopes (no incoherent focal integration), detector MTF and noise.

## 2. The reference

abTEM 1.0.10 at `164e644`, the same venv and configuration (float64, numpy FFT, eager). `transfer.py`: `CTF`
(`_evaluate_from_angular_grid`: the aberration array, then × the spatial envelope if β ≠ 0, × the temporal envelope if Δ ≠ 0,
× the aperture if the cutoff is finite), `SpatialEnvelope._evaluate_from_angular_grid`, `TemporalEnvelope.
_evaluate_from_angular_grid`, `Aberrations`, `Aperture`, `scherzer_defocus`, `point_resolution`, `nyquist_sampling`;
`waves.py` `Waves.apply_ctf`; `transform.py` `ReciprocalSpaceMultiplication` (FFT, × kernel, IFFT).

## 3. The definitions, pinned

On the grid's spatial frequencies, α = |k| λ (rad), φ = atan2(ky, kx), as slice 1's probe (3.8).

1. **Aberrations:** exp(−i χ(α, φ)), χ as slice 1's 3.8 (all orders to C56).
2. **Spatial envelope** (β in mrad, used as β / 10³ rad): E_s = exp(−sign(β) (β/2)² ((∂χ/∂α)² + (∂χ/∂φ)²)) with, as abTEM
   writes them,
   ∂χ/∂α = (2π/λ) [ (C12 cos 2(φ−φ12) + C10) α + (C23 cos 3(φ−φ23) + C21 cos(φ−φ21)) α² + (C34 cos 4(φ−φ34) + C32 cos 2(φ−φ32)
   + C30) α³ + (C45 cos 5(φ−φ45) + C43 cos 3(φ−φ43) + C41 cos(φ−φ41)) α⁴ + (C56 cos 6(φ−φ56) + C54 cos 4(φ−φ54) + C52
   cos 2(φ−φ52) + C50) α⁵ ],
   ∂χ/∂φ = −(2π/λ) [ ½ (2 C12 sin 2(φ−φ12)) α + ⅓ (3 C23 sin 3(φ−φ23) + C21 sin(φ−φ21)) α² + ¼ (4 C34 sin 4(φ−φ34) + 2 C32
   sin 2(φ−φ32)) α³ + ⅕ (5 C45 sin 5(φ−φ45) + 3 C43 sin 3(φ−φ43) + C41 sin(φ−φ41)) α⁴ + ⅙ (6 C56 sin 6(φ−φ56) + 4 C54
   sin 4(φ−φ54) + 2 C52 sin 2(φ−φ52)) α⁵ ].
   (abTEM's ∂χ/∂φ is the derivative divided by α; it is ported as written.)
3. **Temporal envelope** (Δ in Å, the 1/e half-width of the focal spread): E_t = exp(−((½ π / λ) Δ α²)²).
4. **Aperture:** slice 1's (soft by default; none when the cutoff is infinite).
5. **The kernel:** H = exp(−iχ) × E_s (if β ≠ 0) × E_t (if Δ ≠ 0) × aperture (if the cutoff is finite), in that order.
6. **Applying it** (`apply_ctf`): ψ' = IFFT₂(FFT₂(ψ) · H); the image is |ψ'|².
7. **Closed forms:** Scherzer defocus = sign(Cs) √(1.5 |Cs| λ) (abTEM's `defocus` convention, defocus = −C10);
   point resolution = (λ³ |Cs| / 6)^{1/4}; Nyquist sampling = 1 / (4 α_c / λ · 10⁻³) (α_c in mrad), all in Å.

## 4. Pins (agreement with abTEM)

`tools/slice-3/pins.py` writes `Tests/SimulatoProbeTests/Fixtures/slice-3-pins.json` and a float64 sidecar before any
slice-3 Swift, re-deriving every quantity from this text in NumPy and asserting it equals abTEM's.

| Case | Exit wave | CTF |
|---|---|---|
| H1 | SrTiO₃ [001], 3 × 3 × 5 cells, 120 × 120, 300 keV plane wave (slice 1's pipeline, abTEM's π) | Cs 1.2 × 10⁷ Å (1.2 mm), defocus Scherzer, cutoff 20 mrad soft |
| H2 | as H1 | as H1 plus Δ 30 Å and β 0.5 mrad |
| H3 | as H1 | C10 −80 Å, C12 50 Å at 0.4 rad, C21 300 Å at 1.1, C23 200 Å at 0.1, C30 1.2 × 10⁷ Å, Δ 40 Å, β 0.3 mrad, no aperture |
| K1 | none (kernel only), 96 × 160 over 9.6 × 14.4 Å, 80 keV | every coefficient C10 ... C56 with angles, Δ 25 Å, β 0.8 mrad, cutoff 35 mrad hard |
| K2 | none, 128 × 128 over 12.8 Å, 200 keV | β < 0 (abTEM's sign(β) branch), Δ 0, cutoff infinite |

Recorded: the kernels (complex), the images |ψ'|² of H1–H3, and the closed forms at Cs ∈ {−5 × 10⁵, 10⁶, 1.2 × 10⁷} Å and
80, 200, 300 keV, and the Nyquist sampling at 10, 20, 30 mrad.

Agreement required: kernels within 1e−12 × max |H| (CPU, float64); images within 1e−10 × max with the port's own exit wave
(`.abTEMFloat32`) and 1e−6 with `.exact`; the closed forms within 1e−15 relative.

## 5. Truth (independent of abTEM)

1. **The identity CTF.** No aberrations, no envelopes, no aperture: the image equals |ψ|² within 1e−14 × max.
2. **The envelopes' closed forms.** For round aberrations only (C10, C30, C50), ∂χ/∂φ = 0 and ∂χ/∂α = (2π/λ)(C10 α + C30 α³
   + C50 α⁵): E_s and E_t equal the hand formulas within 1e−14 at α = 5, 10, 20, 40 mrad; E_t(0) = E_s(0) = 1; both lie in
   (0, 1] everywhere for β, Δ > 0; E_t(α) = 1/e where (π Δ α² / (2λ)) = 1.
3. **Vacuum is flat.** A plane wave through vacuum, with any CTF of section 3: the image is 1 at every pixel within 1e−14
   (H(0) = 1: χ(0) = 0, the envelopes and the soft aperture are 1 at α = 0).
4. **Weak-phase contrast at Scherzer.** A single carbon atom in a 2 Å slice at 300 keV, Cs 1.2 mm at Scherzer defocus,
   cutoff at the point-resolution angle: (a) the image has its minimum at the atom (atoms dark: sin χ < 0 across the passband
   for positive Cs at Scherzer); (b) the image minus 1 equals the weak-phase prediction 2 IFFT[φ̂ · sin χ · A] (φ = σ V, the
   exit phase, A the aperture) within 0.1 × max |I − 1| (the bound holds the second-order term φ², predicted ~0.04).
5. **The closed forms.** Scherzer defocus at Cs = 1.2 mm, 200 keV = √(1.5 × 1.2 × 10⁷ × λ) = 671.9 Å (textbook ~67 nm);
   point resolution (λ³ Cs / 6)^{1/4} = 2.37 Å; the Nyquist sampling at 20 mrad, 200 keV = λ / (4 α_c) = 0.3135 Å; each by
   its formula within 1e−15, and to four figures as stated.

## 6. Gates and predicted outcomes

1. Pins (section 4). 2. Truth (section 5). 3. Every new test broken once (planned: E_t without the square; E_t with λ
multiplied instead of divided; β not converted from mrad; the ½ of ∂χ/∂φ's first term dropped (its sign is squared
away, so the sign is not a mutant); ∂χ/∂α missing the C30 term; sign(β) dropped; the soft aperture made hard; the kernel
conjugated (exp(+iχ)); the envelopes' order of application exchanged with the aperture (equivalent: real factors commute,
so instead the temporal envelope applied twice); apply_ctf without the inverse FFT's 1/N; Scherzer's 3/2 replaced by
4/3; point resolution with a cube root). 4. One Gate B refuter on the frozen commit (other aberrations, envelopes, grids, signs,
focal series). 5. `tools/test.sh` green with the stamp; a cold release build with zero warnings.

A Metal path for the CTF is not registered: one FFT pair per image does not pay. The Metal performance lane (after this
slice) is registered separately if it moves a number.

## 7. What is not claimed

The envelopes are the quasi-coherent approximation (abTEM's): valid for a weak, slowly varying object and small spreads,
not a substitute for an incoherent sum over defocus and illumination angles. No detector, no noise, no MTF: a simulated
HRTEM image is still not a micrograph.

## 8. Licence

As slices 1 and 2: abTEM's algorithms (GPL-3.0-or-later), ported.

## 9. Amendments

### 2026-10-15 ~01:40 CEST, after the pin generator ran, before any slice-3 Swift

`tools/slice-3/pins.py`: 35 self-checks pass (kernels bit-equal to the text, images within 5.7e−16), deterministic,
1.54 MB; the closed forms gave the predicted 671.9 Å, 2.37 Å and 0.3135 Å. One case is changed:

1. **K2's negative angular spread is dropped, and a negative spread is refused by the port.** With β < 0, abTEM's
   −sign(β) turns the envelope into exp(+(β/2)²(…)), which overflows to infinity at large angles (the generator's comparison
   met NaN): there is no meaningful number to pin. K2 now has β = 0.4 mrad (no aperture, no focal spread); the port refuses
   β < 0 and Δ < 0 with a named error. Section 4's table is read with that change.

# SimulatoKit slice 2: frozen phonons, STEM scans and detectors (registration)

Registered 2026-10-15 (about 01:20 CEST), before any slice-2 Swift, by the SimulatoKit session on the owner's instruction of
00:55 ("we have time and usage to build slice 2 everywhere"). Tier 1: the numbers are images and patterns a scientist reads
(HAADF, ADF and BF intensities, a 4D-STEM pattern, DPC segments) and the thermal averaging that makes them comparable to
experiment. Slice 1's definitions (`docs/registration-slice-1.md`, sections 3.1–3.10 and amendments 1–3) hold unchanged;
this file adds to them. A change after code starts is a dated amendment (section 9).

## 1. What slice 2 makes

- **SimulatoMultislice, frozen phonons.** `FrozenPhonons`: an ensemble of configurations of one structure, each the
  structure with every atom displaced by a Gaussian of standard deviation σ_Z (Å, per element) along the chosen axes
  (default x, y and z). Two sources of displacements: **given** (an array per configuration, e.g. read from abTEM's
  pins, so the port can be checked against abTEM exactly) and **generated** by the kit's own seeded generator (section
  3.2, registered by its statistics, not by abTEM's bits). `ProjectedPotential` gains a constructor from a configuration.
- **SimulatoProbe, scans.** `GridScan` (start, end, gpts, endpoint; positions as abTEM's).
- **SimulatoDetector** (now depending on SimulatoMultislice and SimulatoProbe). The STEM driver
  `STEM.scan(probe:potentials:scan:detectors:backend:)`: for every configuration and every position, the probe at that
  position through that configuration's potential, each detector applied to the exit wave, the measurement averaged over
  configurations (abTEM's `ensemble_mean=True`). `AnnularDetector` (inner, outer mrad: BF, ADF, HAADF), `PixelatedDetector` (the pattern per
  position, max angle "cutoff", "valid" or "full"), `SegmentedDetector` (radial × azimuthal bins between inner and outer,
  with a rotation: DPC). Each detects one exit wave; the driver assembles images (gpts of the scan) or a 4D cube.
- **Metal.** The STEM driver's `.metal` backend propagates a batch of probe positions at once (one MPSGraph over a
  [batch, nx, ny] tensor, the transmission broadcast); detection is on the CPU from the read-back exit waves. The CPU
  path (float64) stays the reference.

Not in slice 2: partial coherence (spatial and temporal envelopes), the CTF and HRTEM images (slice 3), `FlexibleAnnular`
detectors, detector offsets, resampling of pixelated patterns, line scans and custom scans beyond a list of positions,
anisotropic σ (a 3-vector per element), PRISM, the Metal performance lane (after slice 2).

## 2. The reference

abTEM 1.0.10 at `164e644`, the same venv and configuration as slice 1 (float64, numpy FFT, eager). What is used:
`inelastic/phonons.py` (`FrozenPhonons`, `validate_seeds`, `randomize`), `potentials/iam.py` (`_prepare_atoms`: displace,
then wrap, then the top-face rule), `scan.py` (`GridScan.get_positions`), `detectors.py` (`AnnularDetector`,
`PixelatedDetector`, `SegmentedDetector`), `measurements.py` (`_annular_detector_mask`, `integrate_radial`,
`_polar_detector_bins`, `_radial_binning`, `polar_binning`), `waves.py` (`Probe.scan`, `diffraction_patterns` with
`parity`, `_gpts_within_angle`, `rectangle_cutoff_angles`, `cutoff_angles`).

## 3. The definitions, pinned

### 3.1 Frozen phonons (abTEM's)

Configuration c of a structure of N atoms: seed s_c; r = `default_rng(s_c).normal(size=(N, 3))`; for each axis a in the
directions, position[:, a] += σ_{Z(atom)} · r[:, a]. Seeds: an integer seed with n configurations makes n seeds by
`default_rng(seed).integers(2³¹ − 1)` drawn until n distinct (`validate_seeds`); a tuple gives them directly. The displaced
structure then goes through slice 1's section 3.4 (wrap, top-face rule, slicing): displacement first, wrap after
(`_prepare_atoms`, periodic). σ is per element (a dictionary) or one value for all.

### 3.2 The kit's own generator (not abTEM's bits)

Bit-reproducing NumPy's PCG64 and its ziggurat normal sampler is not in scope. The kit's generator is SplitMix64 seeding a
xoshiro256** stream; normals by the Box–Muller transform, both outputs used. Same seed, same displacements, on every
platform (integer arithmetic and one `log`, `sqrt`, `cos`, `sin` per pair). Registered by its statistics (section 5).

### 3.3 Scan positions

`GridScan(start, end, gpts, endpoint)`: positions along each axis `linspace(start, end, gpts, endpoint)` (NumPy's, x first,
meshgrid "ij"); the image is gpts.x × gpts.y, x the first index.

### 3.4 Detectors (on one exit wave ψ of a probe, normalisation "reciprocal space", so no 1/N)

Angular sampling Δα = λ / L per axis (mrad: × 10³). The full unshifted pattern P(k) = |FFT₂(ψ)|².

- **Annular** (`integrate_radial`): Σ P over the pixels with inner² ≤ α² < outer² (α the pixel's angle in mrad on the
  `fftfreq` grid; half-open). Outer must not exceed the pattern's maximum angle min(⌊n/2⌋ Δα) beyond `isclose` (atol 1e−5);
  `outer = nil` means ⌊min(cutoff angles)⌋ (the cutoff angle of an axis is ⌊n'/2⌋ Δα, n' the cutoff gpts of slice 1's
  3.10 with parity "same").
- **Pixelated** (`diffraction_patterns(max_angle=..., parity="same")`): "cutoff" crops to the cutoff gpts with the grid's
  parity (slice 1's rule, parity "same" instead of "odd"); "valid" to ⌊n'/√2⌋ (abTEM's `safe_floor_int`) moved to the
  grid's parity; "full" keeps all; then `fftshift`. abTEM's default is "valid".
- **Segmented** (`polar_binning`): the pattern cropped to max_angle = outer (gpts 2⌈outer/Δα⌉ + 1, moved to the grid's
  parity), `fftshift`ed; each pixel's α and φ = atan2(ky, kx) on the cropped grid's `fftfreq` (then shifted); φ ← (φ −
  rotation) mod 2π; a pixel with inner ≤ α < outer is in radial bin ⌊n_r (α − inner)/(outer − inner)⌋ and azimuthal bin
  clip(⌊n_φ φ / 2π⌋, 0, n_φ − 1); each bin sums its pixels. Output n_r × n_φ per position.

### 3.5 Averaging

The measurement of a scan is the mean over configurations of each detector's value at each position (abTEM's
`ensemble_mean`), in configuration order.

### 3.6 Metal batching

A batch of B positions: ψ of shape [B, nx, ny] complex64; per slice ψ ← IFFT(FFT(ψ · t) · P) over axes 1 and 2, t and P
broadcast. The batch size is chosen so the batch's waves fit in 512 MB; results in position order.

## 4. Pins (agreement with abTEM)

`tools/slice-2/pins.py` writes `Tests/SimulatoDetectorTests/Fixtures/slice-2-pins.json` and a float64 sidecar before any
slice-2 Swift, re-deriving every quantity from this text in NumPy and asserting it equals abTEM's (a failure stops it).

| Case | Structure | Grid | Beam and scan | Configurations |
|---|---|---|---|---|
| F1 | SrTiO₃ [001], 2 × 2 × 4 cells | 120 × 120 | plane wave 200 keV, pattern averaged | 4, σ Sr 0.08, Ti 0.06, O 0.07 Å, seed 7 |
| S1 | Al [001] fcc, 2 × 2 × 3 cells (8.1 × 8.1 × 12.15 Å) | 128 × 128 | probe 200 keV, 20 mrad, defocus 0; GridScan (0,0)–(4.05,4.05), 6 × 6, endpoint False | 3, σ Al 0.085 Å, seed 3 |
| S2 | as S1, σ = 0 | 128 × 128 | as S1 | 1 (the frozen-phonon machinery equals slice 1) |
| S3 | vacuum 8.1 × 8.1 × 6 Å | 128 × 128 | as S1 | 1 |
| S4 | single atoms C, Si, Cu, Au, each at the centre of its own 8 × 8 × 2 Å box | 128 × 128 | probe 200 keV 25 mrad at the atom | 1 |

Detectors on S1–S4: Annular BF (0–20), ADF (30–60), HAADF (70–130 mrad), Annular with outer = nil; Pixelated "cutoff",
"valid" and "full" (patterns stored for two positions only); Segmented n_r = 2, n_φ = 4, inner 0, outer 20, rotation 0
and rotation 0.3 rad.

Recorded: the seeds abTEM derives, every configuration's displacements, V of slice 0 for configuration 0 (F1, S1), the
scan positions, each detector's averaged image, the pixelated patterns of two positions, the segmented values, the
per-configuration HAADF of one position (to localise a disagreement).

Agreement required, CPU path (float64) with `.abTEMFloat32` and abTEM's displacements given (relative to the array's max):
displacements and positions bit-equal (they are read, not computed); V 1e−11; annular images and segmented values 1e−10;
pixelated patterns 1e−10; the averaged F1 pattern 1e−10. With `.exact`: 1e−6. A difference beyond a bound is a failure of
the port.

Metal (float32, batched) against the CPU path, same potentials: annular images 2e−4 × max, segmented 2e−4 × max, pixelated
2e−4 × max, up to 100 slices (slice 1's bounds; the per-slice Σ|ψ|² bound of amendment 2, accepted by the owner, applies to
each position's exit wave).

## 5. Truth (independent of abTEM)

1. **The generator's statistics.** 10⁶ normals from seed 1: mean within 5 × 10⁻³ (5 σ/√n), variance within 1 ± 7 × 10⁻³,
   fourth moment 3 ± 0.05, the lag-1 correlation within ± 5 × 10⁻³; two seeds give different streams; the same seed the
   same stream. The displacements of 10⁴ atoms with σ = 0.1 Å have a per-axis standard deviation within 0.1 ± 0.003 Å and
   the axes not displaced (directions "xy") are untouched exactly.
2. **σ = 0 is slice 1.** A frozen-phonon ensemble with σ = 0 gives exactly the slice-1 potential and exit waves (S2,
   equality to the last bit), and averaging n identical configurations gives the single result within 1e−15 relative.
3. **A detector over all angles collects the total.** Annular 0 to ≥ the cutoff angle equals N Σ_r |ψ|² (Parseval, the
   probe's normalisation) within 1e−12 relative, at every position.
4. **BF + DF = total.** Annular (0, a) + annular (a, b) = annular (0, b) exactly (half-open masks, disjoint), and the
   segments of a segmented detector sum to the annular detector with the same inner and outer within 1e−12.
5. **An atom-free scan is flat.** In vacuum (S3) every detector's image is constant over the positions within 1e−12
   relative (the probe's shift is a phase ramp; |Ψ(k)| does not depend on the position).
6. **HAADF rises with Z.** S4: the HAADF signal of a probe on a single atom increases strictly from C to Si to Cu to Au.
7. **DPC symmetry.** In vacuum the four segments of a centred, aberration-free probe (rotation 0) are equal within 1e−12
   (the soft aperture is symmetric under the quarter turns of a square grid).
8. **The pixelated pattern sums to the total**: Σ pattern ("full") = N Σ|ψ|² within 1e−12.

## 6. Gates and predicted outcomes

1. Pins (section 4) on every case. 2. Truth (section 5). 3. Metal against CPU (section 4), the tests enabled only with a
device and the reason named. 4. Every new test broken once by a mutant (planned: σ applied twice; displacement after the
wrap instead of before; the directions ignored; the seed rule's distinctness dropped; positions with endpoint inverted;
the annular mask closed (≤ outer); inner applied as > instead of ≥; the probe pattern divided by N; outer = nil without
the floor; the pixelated "valid" without √2; the segment's φ from atan2(kx, ky); the rotation added instead of
subtracted; the radial bin rounded instead of floored; the average a sum; on Metal, the batch's FFT over axes 0 and 1).
5. One Gate B refuter on the frozen commit (other structures, σ, seeds, detectors at the cutoff, scans off the cell, many
configurations). 6. `tools/test.sh` green with the stamp; a cold release build with zero warnings.

## 7. What is not claimed

Frozen phonons with Einstein-model σ are an approximation to thermal diffuse scattering; σ is the caller's (no Debye–Waller
table in slice 2). Agreement with abTEM uses abTEM's displacements; the kit's own generator agrees with abTEM only in
distribution. No partial coherence: a simulated HAADF image is not yet a micrograph. Detector efficiency and noise are not
modelled.

## 8. Licence

As slice 1: the algorithms are abTEM's (GPL-3.0-or-later), ported; xoshiro256** and SplitMix64 are public domain
(Blackman and Vigna).

## 9. Amendments

### 2026-10-15 01:45, after the pin generator ran, before any slice-2 Swift

`tools/slice-2/pins.py`: 98 self-checks pass (largest text-against-abTEM difference 1.1e−14 outside the segmented
detector), deterministic (sha256 equal on two runs), 1.48 MB. Two readings are settled here.

1. **abTEM's segmented detector accumulates in float32.** `_radial_binning` sums the segments into a `float32` array with
   numba `fastmath` (`measurements.py` 3056–3061), even in float64 mode; summing the same pattern over abTEM's own bins in
   float64 equals the text exactly. The segmented values therefore agree with abTEM only to float32 rounding (measured
   2e−7 × max). The port sums in float64. Pins: the port against the text's float64 sum (stored as `<case>.<det>.float64`)
   within 1e−10, and against abTEM's output within 1e−6 × the scale below.
2. **The scale of a detector bound.** An image that is pure rounding noise (vacuum ADF and HAADF outside a 20 mrad probe,
   ~1e−32) cannot be judged relative to its own maximum. Every detector bound of section 4 is relative to the larger of the
   image's maximum and the scan's total intensity (the maximum of the annular image with outer = nil), stored per detector
   as `scale`.

### 2026-10-15 02:05, at the first slice-2 test run: truth 7 corrected, and the kit generator's seeds

1. **Truth 7 was false as written, and is corrected.** "In vacuum the four segments of a centred, aberration-free probe are
   equal" fails with inner = 0: the zero frequency has φ = atan2(0, 0) = 0, so it belongs to the first segment alone, and it
   carries the probe's largest value (abTEM's rule and the port's alike; the pins against abTEM passed). The truth is now:
   with inner > 0 (2 mrad, which excludes only the origin), the four segments are equal within 1e−12; with inner = 0 the
   first is larger and the other three are equal within 1e−12. My first reading, that −π/2 mod 2π rounds into the third
   quadrant, was checked and is wrong (it lands in the fourth, as it should).
2. **The kit generator's configurations** (section 3.2, made concrete): configuration c uses the generator seeded with
   seed + c (wrapping); three normals per atom in atom order (x, y, z); an axis not in the directions is set to zero after
   drawing (so the streams do not depend on the directions); an element absent from σ is not displaced.

All other slice-2 tests passed on the first run: frozen phonons (F1, four configurations, the averaged pattern), every
detector on S1–S3 under both π choices, the scan positions bit-equal, S4's HAADF on C, Si, Cu and Au, and the batched Metal
scan against the CPU scan.

# SimulatoKit slice 1: projected potential and multislice for a plane wave and a convergent probe (registration)

Registered 2026-10-15 (night of 2026-10-14), before any Swift code, by the SimulatoKit session (the coordinator's side
session; the owner's direction in the suite `ROADMAP.md`, "Kits", phase 5: SimulatoKit slice 1 is a projected potential
and a multislice image of a thin crystal against abTEM at a pinned commit, Metal first with a CPU reference and a
registered float32-against-float64 tolerance). Tier 1: every number here is one a scientist reads or one the other kits
will use as truth (a wavelength, a projected potential, an exit wave, an image or diffraction intensity). What is written
here is fixed before the code; a change after code starts is an amendment, dated, with its reason (section 9).

## 1. What slice 1 makes

Four products, one per module; SimulatoCore and SimulatoMultislice and SimulatoProbe carry slice 1, SimulatoDetector is an
empty module with its planned scope (slice 2). Public API only; nothing an outside user needs is `package`.

- **SimulatoCore.** `PhysicalConstants` (CODATA 2014, the values abTEM reads from ASE; section 3.1); `ElectronBeam` (energy
  in eV, the relativistic wavelength λ in Å, the interaction parameter σ in 1/(V·Å)); `AtomicStructure` (an orthogonal box
  Lx, Ly, Lz in Å and atoms (Z, Cartesian x, y, z in Å); `repeated(nx, ny, nz)`; a minimal local type, **marked for
  replacement** by DiffractoKit's `CrystalStructure` when DiffractoKit has its first tag; `git -C ../DiffractoKit tag` is
  empty on 2026-10-15 00:00); `Grid` (gpts nx × ny, extent, sampling; the FFT sizes of section 3.9); `ComplexField`
  (a complex array on the grid, x the first index, row-major `[ix * ny + iy]`, abTEM's `[x, y]` layout); `Waves` (a
  complex field with its energy: `intensity()`, `diffractionPattern(maxAngle:)`); the antialias aperture (3.6); the CPU FFT
  (Accelerate `vDSP.DiscreteFourierTransform<Double>`, package).
- **SimulatoMultislice.** `LobatoParameterization` (the table of 3.3, the scaled parameters, the projected scattering
  factor, the π choice of 3.3); the slicing of 3.4; `ProjectedPotential` (the slices of 3.5); the transmission function and
  the Fresnel propagator (3.6); `Multislice` (3.7) with exit planes; one API with two backends, `.cpu` (float64,
  Accelerate) and `.metal` (float32, MPSGraph), the backend a parameter, the CPU path the reference.
- **SimulatoProbe.** `Aberrations` (abTEM's polar coefficients C10 ... C56 and angles; `defocus` = −C10); `Aperture` (the
  soft edge of 3.8); `Probe` (energy, grid, semiangle cutoff, aberrations; built at a position, normalised as abTEM does).
- **SimulatoDetector.** Empty: bright field, annular (ADF, HAADF), pixelated (4D) and segmented detectors are slice 2,
  with frozen phonons.

Not in slice 1: frozen phonons, partial coherence, detectors that integrate over a scan, STEM scans, the finite
projection integrals (`projection="finite"`), the Kirkland, Peng and other parameterisations, PRISM, tilt, non-orthogonal
cells, a CIF importer in the library (the example executable carries a minimal reader, marked for replacement by
DiffractoKit's `CIFImport`).

## 2. The reference

**abTEM 1.0.10** (GPL-3.0-or-later; `github.com/abtem/abtem`, tag `v1.0.10` = commit
`164e644ff411167db23ec4ef26aeb4fe8f384b26`, 2026-07-06), installed from PyPI into the gitignored `References/abtem-venv`
(Python 3.14.7) and the same tag cloned into `References/abtem-164e644` so the definitions below cite file:line. The venv
holds abTEM 1.0.10, ASE 3.29.0, NumPy 2.5.4, SciPy 1.18.1, numba 0.68.0, pyfftw 0.15.1 (`References/abtem-venv-freeze.txt`).
The oracle runs as `References/abtem-venv/bin/python -I -B` with abTEM's configuration set to
`precision: float64`, `fft: numpy`, `dask.lazy: false` (abTEM's defaults are float32 and FFTW; float64 makes the oracle a
float64 reference, and NumPy's pocketfft makes it deterministic, with no FFTW plan chosen by timing). An unpinned copy is
never a reference. Agreement with abTEM is agreement, not truth; truth is section 5. Kirkland, *Advanced Computing in
Electron Microscopy* (2nd ed., 2010), is the text for the method; Dr. Probe is a second oracle for a later slice, not
this one.

What of abTEM is used, and for what (paths under `abtem/`):

- `core/energy.py`: `energy2wavelength`, `energy2sigma` (3.2).
- `parametrizations/__init__.py` `LobatoParametrization` (its `scaled_parameters`) and `parametrizations/functions/
  lobato.py` `projected_scattering_factor` (3.3); `parametrizations/data/lobato.json` (the table).
- `integrals.py` `ScatteringFactorProjectionIntegrals.integrate_on_grid`, `superpose_deltas`, `sinc` (3.5): the default
  `Potential(..., parametrization="lobato", projection="infinite")`.
- `slicing.py` `_validate_slice_thickness`, `SliceIndexedAtoms`; `potentials/iam.py` `_FieldBuilderFromAtoms._prepare_atoms`
  (3.4).
- `potentials/iam.py` `PotentialArray._transmission_function`; `multislice.py` `_fresnel_propagator_array`,
  `conventional_multislice_step`, `multislice_and_detect` (3.6, 3.7); `antialias.py` `antialias_aperture` and the
  configuration `antialias.cutoff` 0.6666666, `antialias.taper` 0.01 (`core/abtem.yaml`).
- `waves.py` `PlaneWave`, `Probe._calculate_array`, `_WavesNormalization`, `Waves.diffraction_patterns`,
  `_antialias_cutoff_gpts`; `transfer.py` `Aperture`, `soft_aperture`, `Aberrations._evaluate_from_angular_grid`;
  `core/fft.py` `fft_shift_kernel`, `fft_crop` (3.8, 3.10).

## 3. The definitions, pinned

Units: lengths in Å, energies in eV, angles in mrad where abTEM uses mrad, potentials in V (projected: V·Å), spatial
frequencies k in 1/Å without 2π (abTEM's `fftfreq`). Arrays are `[x, y]`.

### 3.1 Constants

ASE 3.29.0's CODATA 2014 values, which abTEM reads: h = 6.62607004e−34 J s, mₑ = 9.10938356e−31 kg, c = 299792458 m/s,
e = 1.6021766208e−19 C, ε₀ = 8.85418781762039e−12 F/m, Bohr a₀ = 0.5291772105638411 Å. abTEM's `kappa`
(`core/constants.py`) = 4π ε₀ / (2π a₀ e) in ASE units = **0.0208865737082965** (the pin generator records the value and
the Swift constant is that double).

### 3.2 Wavelength and interaction parameter

λ(E) = h c / √(E (2 mₑ c² / e + E)) / e × 10¹⁰ Å (E in eV); σ(E) = 2π m(E) e λ / h² in abTEM's units, m(E) = mₑ (1 + e E /
(mₑ c²)), which is σ = 2π mₑ (1 + eE/mₑc²) e λ / h² with λ in m and the result converted to 1/(V·Å). The pin generator
records λ and σ at 60, 80, 100, 120, 200, 300 keV.

### 3.3 The parameterisation: Lobato and Van Dyck (2014), abTEM's default

I. Lobato, D. Van Dyck, Acta Cryst. A 70, 636–649 (2014). Five pairs (aᵢ, bᵢ) per element, Z = 1 ... 103, from abTEM's
`lobato.json` (the paper's table, written into a generated Swift source with its provenance; section 4). Scaled for the
projected scattering factor (`scaled_parameters`): Aᵢ = π² aᵢ / bᵢ^{3/2} / κ, Bᵢ = 2π / √bᵢ. The projected scattering
factor of an infinite projection, a function of k² (1/Å²):

f(k²) = 8π Σᵢ [ Aᵢ / Bᵢ / (q + Bᵢ²) + Aᵢ Bᵢ / (q + Bᵢ²)² ],  q = 4π² k²,

in V·Å² (its integral over the plane is the projected potential's). **The π choice.** abTEM evaluates the π and π² of
this formula as float32 constants (`lobato.py`, `pi = np.array(np.pi, dtype=np.float32)`), even in float64 mode; measured
on Al at 2026-10-15 00:05, its f differs from the exact-π f by 2.8e−8 to 5.3e−8 relative and equals the float32-π
evaluation to 4e−16. The port's default is the exact π (`.exact`); `.abTEMFloat32` reproduces abTEM's numbers (the
float32-rounded π = 3.1415927410125732 and π² = 9.869604110717773, used exactly where abTEM uses them). The scaling of A
and B uses the exact π in both, as abTEM does. No Debye–Waller damping in slice 1 (abTEM's `sigmas` is empty by default).

### 3.4 The structure and its slices

The box must be orthogonal (abTEM refuses otherwise). Before slicing, every atom is wrapped into the box (fractional
coordinate mod 1, abTEM's `atoms.wrap(eps=0)`), and an atom with z > Lz − 1e−10 is set to z = 0. Slice thicknesses: a
single value s gives n = ⌈Lz / s⌉ slices of Lz / n each (`_validate_slice_thickness`); a list must sum to Lz within
`np.isclose` (rtol 1e−5, atol 1e−8). An atom is in slice j when e_{j−1} ≤ z < e_j with e_j = (Σ_{i≤j} sᵢ) − 1e−12
(`np.digitize` on the nudged cumulative edges, `SliceIndexedAtoms`); an atom past the last edge goes into the last slice.
The infinite projection takes no margin and no padding.

### 3.5 The projected potential of a slice

For each element present in the slice, with its atoms' (x, y) in pixel units u = x / Δx, v = y / Δy:

1. Bilinear deltas (`superpose_deltas`, `round_positions=False`): i = ⌊u⌋, j = ⌊v⌋, fx = u − i, fy = v − j; add
   (1 − fx)(1 − fy) at (i, j), fx (1 − fy) at (i + 1, j), (1 − fx) fy at (i, j + 1), fx fy at (i + 1, j + 1), indices mod
   (nx, ny). Note abTEM's weight expression `1 + xy − y − x` is the same number.
2. D = FFT₂(deltas) (unnormalised forward transform, NumPy's sign: e^{−2πi k·r}).
3. Multiply by f(k²) / S(k), with k² = kx² + ky² on the `fftfreq` grid and S abTEM's compensation (`integrals.py`
   `sinc`): S(k) = sin(κ) / κ × Δx Δy with κ = √((kx Δx)² + (ky Δy)²) (a dimensionless frequency in cycles per pixel,
   **not** multiplied by π: abTEM's definition, ported as it is), S(0) = Δx Δy.
4. V_element = Re(IFFT₂(·)) (NumPy's normalised inverse, 1/(nx ny)); the slice's V is the sum over elements.

V is in V·Å. A slice with no atoms has V = 0 exactly.

### 3.6 Transmission, propagation, antialiasing

- Antialias aperture A(k) (`antialias_aperture`): k_c = 0.6666666 / max(Δx, Δy) / 2, taper w = 0.01 / max(Δx, Δy);
  A = 1 for |k| ≤ k_c − w; A = ½ (1 + cos(π (|k| − k_c + w) / w)) for k_c − w < |k| ≤ k_c; A = 0 for |k| > k_c.
- Transmission (`_transmission_function`, then `AntialiasAperture.bandlimit`): t = exp(i σ V), then
  t ← IFFT₂(FFT₂(t) · A).
- Fresnel propagator over Δz (order 1, `_fresnel_propagator_array`): P(k) = exp(−i π λ Δz k²) · A(k).

### 3.7 Multislice

ψ ← IFFT₂(FFT₂(ψ · t_j) · P_{Δz_j}) for j = 1 ... n: transmit, then propagate by that slice's thickness, the last slice
included (`conventional_multislice_step`, `transpose=False`). Exit planes: the wave after slice j for each requested j
(abTEM's `exit_planes`, given as every m-th slice, the last always included). The incident plane wave is ψ = 1 at every
pixel (`PlaneWave(normalize=False)`); the incident probe is section 3.8's.

### 3.8 The probe

On the reciprocal grid, α = |k| λ (rad), φ = atan2(ky, kx):

- Aperture (`Aperture`, `soft=True`, abTEM's default): a(k) = clip((α_c − α) / d + ½, 0, 1) with d = √((cos φ Δα_x)² +
  (sin φ Δα_y)²), Δα = λ / (n Δ) the angular sampling in rad, α_c the semiangle cutoff in rad; a(0) = 1.
- Aberrations (`Aberrations._evaluate_from_angular_grid`): χ = (2π/λ) [ ½ α² (C10 + C12 cos 2(φ − φ12)) + ⅓ α³ (C21
  cos(φ − φ21) + C23 cos 3(φ − φ23)) + ¼ α⁴ (C30 + C32 cos 2(φ − φ32) + C34 cos 4(φ − φ34)) + ⅕ α⁵ (C41 cos(φ − φ41) + C43
  cos 3(φ − φ43) + C45 cos 5(φ − φ45)) + ⅙ α⁶ (C50 + C52 cos 2(φ − φ52) + C54 cos 4(φ − φ54) + C56 cos 6(φ − φ56)) ], the
  factor exp(−i χ). Coefficients in Å, angles in rad; defocus = −C10.
- Position (`fft_shift_kernel`): the factor exp(−2πi (kx' u + ky' v)) with (u, v) the position in pixels and k' the
  `fftfreq` of the grid at unit spacing.
- Ψ(k) = shift · a · exp(−iχ); normalised so Σ_k |Ψ(k)|² = 1 (`normalize()`, reciprocal space); ψ = IFFT₂(Ψ). Hence
  Σ_r |ψ|² = 1 / (nx ny) (NumPy's inverse), abTEM's convention, kept.

### 3.9 Grids

nx, ny each of the form f · 2^m with f ∈ {1, 3, 5, 15} and n ≥ 8 (the lengths `vDSP.DiscreteFourierTransform` accepts;
the Metal path is held to the same set so both backends accept the same grids). Any other size is a refusal with a named
error. abTEM accepts any size; the pins use only these.

### 3.10 Outputs

- Intensity: |ψ|² per pixel.
- Diffraction pattern (`Waves.diffraction_patterns`, `max_angle="cutoff"`, `fftshift=True`): for a plane wave (metadata
  `normalization: "values"`) ψ / (nx ny) first, then Ψ = FFT₂, then cropped to the antialias cutoff gpts (n' = ⌊(2/3) /
  max(Δ) × L⌋ per axis, with abTEM's rounding `floor(round(·, 7))`, adjusted to n's parity by +1), the crop keeping the
  first ⌈n'/2⌉ and the last ⌊n'/2⌋ frequencies, then |Ψ|², then `fftshift`. For a probe no division (its metadata carries
  no `normalization`). `maxAngle: .full` keeps the whole grid.

## 4. Pins (agreement with abTEM)

A generator, `tools/slice-1/pins.py` (run in `References/abtem-venv`), writes `Tests/SimulatoMultisliceTests/Fixtures/
slice-1-pins.json` (and, if the arrays are large, binary float64 little-endian sidecars named in it, each with its sha256)
before any Swift exists. It also writes `Sources/SimulatoMultislice/LobatoTable.swift` from `lobato.json`, with the file's
sha256 in the header. It records the versions (abTEM, ASE, NumPy, SciPy, Python) and the configuration. The cases
(monoclinic is DiffractoKit's; here the order is the task's):

| Case | Structure | Grid | Slices | Beam |
|---|---|---|---|---|
| A1 single atom | Au at (5.0, 5.0, 1.0) in a 10 × 10 × 2 Å box | 128 × 128 | 1 | plane wave 200 keV |
| A2 single atoms | C at (4.3, 6.1, 0.5) and Si at (12.7, 3.2, 1.5), off-grid, 16 × 16 × 2 Å | 160 × 160 | 1 | plane wave 300 keV |
| B1 Al [001] | fcc Al a = 4.05 Å, cubic cell 4 × 4 × 10 (16.2 × 16.2 × 40.5 Å), 40 atoms per layer pair | 256 × 256 | 20 of 2.025 Å, exit planes every 2 (one per cell) | plane wave 200 keV |
| B2 Al [001] probe | as B1 | 256 × 256 | 20, exit planes every 2 | probe 200 keV (below) |
| C1 vacuum | no atoms, 16.2 × 16.2 × 20 Å | 256 × 256 | 10 of 2 Å | plane wave and probe (as B2) |
| D1 SrTiO₃ [001] | perovskite a = 3.905 Å, 3 × 4 × 5 cells, two elements per slice, non-square grid | 120 × 160 | 10 of 1.9525 Å | plane wave 300 keV and probe 300 keV α_c 25 mrad at (5.3, 7.9) Å |

B2's probe: α_c = 20 mrad, C10 = −50 Å (defocus 50 Å), C30 = 1.0e4 Å (1 µm), C12 = 20 Å at φ12 = 0.3 rad, at the
off-grid position (8.37, 7.91) Å. C1's probe: the same aberrations at the box centre.

Recorded per case: λ, σ; per slice the projected potential V (both π choices: abTEM's own output, and the generator's
re-evaluation with the exact π, section 3.3); the transmission function after bandlimiting (first slice); the incident
wave; the exit wave at every exit plane; the intensity and the diffraction pattern (cutoff) at every exit plane; for
probes the incident probe and its Σ|Ψ(k)|². The generator also re-derives each quantity from this text in NumPy (written
from the registration, not from abTEM's code paths) and asserts it equals abTEM's to the tolerances below; a failure stops
the generator. That proves the text of section 3 is complete before the Swift is written.

Agreement required, CPU path (float64) with `.abTEMFloat32` against abTEM float64 (relative to the array's max |·|):

| Quantity | Bound | Predicted |
|---|---|---|
| λ, σ, κ | 1e−15 relative | equal or 1 ulp |
| f(k²) at the pinned k² | 1e−14 relative | ~4e−16 |
| V per slice | 1e−11 × max V | ~1e−14 |
| bandlimited t | 1e−11 × max | ~1e−14 |
| incident probe | 1e−12 × max | ~1e−15 |
| exit wave per exit plane | 1e−10 × max | ~1e−13 |
| intensity, diffraction pattern | 1e−10 × max | ~1e−13 |

With `.exact` (the default) against abTEM: f within 1e−7 relative; V within 1e−7 × max V; exit waves, intensities and
patterns within 1e−6 × max (predicted ~5e−8: the phase σ V changes by σ · ΔV). A difference beyond a bound is a failure of
the port, never a tolerance to widen.

Agreement required, Metal path (float32, MPSGraph) against the CPU path (float64), on the same potential (built by the
CPU path in float64, cast to float32 transmission functions once), for every case above and the refuter's cases, up to
100 slices on grids up to 512 × 512:

| Quantity | Bound | Predicted |
|---|---|---|
| exit wave | 1e−4 × max |ψ| | ≤ 3e−5 |
| intensity | 2e−4 × max | ≤ 5e−5 |
| diffraction pattern (cutoff) | 2e−4 × max | ≤ 5e−5 |
| total intensity Σ|ψ|² | 1e−5 relative; amended to 1.5e−6 × slices (amendment 2) | ≤ 3e−6 (measured 8e−7 per slice) |

These are float32 rounding through 2n FFTs: about 1e−7 × log₂(N) per transform, √(2n) accumulation. A Metal result beyond
the bound is a failure of the Metal path.

## 5. Truth (independent of abTEM)

1. **Wavelength and σ** against Kirkland's tabulated values (2010, Table 2.2, four significant figures): λ = 0.03701 Å at
   100 keV, 0.02508 at 200, 0.01969 at 300; σ = 9.244e−4 at 100 keV, 7.288e−4 at 200, 6.526e−4 at 300 (1/(V·Å)). Equal to
   the stated digits.
2. **The k = 0 identity of 3.5.** Σ_pixels V Δx Δy = Σ_atoms f(0) for every slice, to 1e−12 relative: S(0) = ΔxΔy and the
   bilinear weights sum to one, so the method conserves each atom's integrated projected potential exactly.
3. **The far field of one atom.** A single Al atom on a grid point at the centre of a 25.6 × 25.6 Å box, 512 × 512: the
   FFT potential at 1 ≤ r ≤ 4 Å equals the analytic Lobato projected potential V(r) = 2 Σᵢ [2 Aᵢ/Bᵢ K₀(Bᵢ r) + Aᵢ r K₁(Bᵢ r)]
   (abTEM's `projected_potential`, its real-space closed form, evaluated by SciPy `kn` in the generator) within 1e−2
   relative. The bound holds the aliasing of f's 1/k² tail and the compensation S; it is a truth about the method, not
   about abTEM's code.
4. **Vacuum.** With no atoms: the plane wave's exit wave equals 1 at every pixel within 1e−14; the probe after a vacuum slab
   of thickness T equals the probe built with C10 + T within 1e−12 × max (P_Δz = exp(−iπλΔz k²) = exp(−iχ) with C10 = Δz),
   for T = 20 Å and every exit plane of C1.
5. **The symmetry of Al [001]** (corrected by amendment 2). With the extent a whole number of cells and the sampling
   dividing a/2 (B1: a/2 = 32 pixels): with slices of a/2 every pattern has its weight only on (h, k) with h + k even (in
   units of 1/a), and with one slice per cell only on h and k both even; off-lattice intensity ≤ 1e−20 × max. The pattern
   has the 4mm symmetry: I(h, k) = I(k, h) = I(−h, k) within 1e−12 × max.
6. **The probe.** In vacuum Σ_k |Ψ(k)|² = 1 within 1e−14; an aberration-free probe (C10 = 0, α_c 20 mrad) is
   centrosymmetric about its position, |ψ(r₀ + d)|² = |ψ(r₀ − d)|² within 1e−12 × max, when r₀ is on a pixel; a probe
   placed at r₀ + (p Δx, q Δy) equals the probe at r₀ circularly shifted by (p, q) pixels within 1e−12 × max.
7. **χ.** For given α, φ and coefficients, the phase equals the closed form of 3.8 written by hand in the test (C10, C12 with
   φ12, C30, C23 with φ23) within 1e−12 relative.
8. **The thickness series.** In B1, the exit wave at exit plane k of the 10-cell run equals the exit wave of a separate
   k-cell run sliced the same way within 1e−13 × max; the total intensity Σ|ψ|² never increases with thickness (the
   propagator has modulus ≤ 1, the transmission is unitary before bandlimiting, the bandlimit removes) — checked plane by
   plane; and, for the plane wave, the central beam's intensity of the pattern differs from 1 by less than the summed
   intensity of all other beams plus the loss, within 1e−12 (Σ over the pattern equals the total intensity, Parseval).

## 6. Gates and predicted outcomes (registered now)

1. Pins (section 4): every case, every quantity, within its bound, CPU `.abTEMFloat32` and CPU `.exact`.
2. Truth (section 5): every item.
3. Metal (section 4, second table): every case, when a Metal device exists; otherwise the tests skip with the named
   reason "no Metal device (MTLCreateSystemDefaultDevice() returned nil)".
4. Every new test broken once by a mutant aimed at its rule (planned: λ without the relativistic term; σ with the rest mass;
   A and B swapped in the scaling; the compensation S dropped; S with π inside (sin(πκ)/(πκ)); the bilinear weights' fx
   and fy swapped; the propagator's sign flipped; the propagator over the whole thickness; the antialias cutoff ½ instead
   of ⅔; the transmission not bandlimited; propagate before transmit; the top-boundary wrap (z > Lz − 1e−10 → 0) dropped;
   the slice-edge nudge dropped; the probe aperture hard; the probe normalised in real space; χ's sign flipped; the shift
   kernel's sign flipped; the crop's parity rule dropped; the plane wave's 1/N dropped in the pattern; on Metal, the
   propagator not applied, the inverse FFT unscaled). The proof is recorded in `docs/slice-1-mutants-<date>.md`.
5. One Gate B refuter on the frozen slice-1 commit (a second only on a named doubt), asked to compare the port with abTEM on
   structures and thicknesses the fixtures do not hold (other elements, 80 and 300 keV, thicker crystals, other grids,
   atoms on slice boundaries and on the box's top face, probes near the box edge, large aberrations) and to try to break
   the Metal tolerance; Gate D only if a discrepancy's cause is unproven.
6. `tools/test.sh` green with the stamp; a cold build with zero warnings; the example executable runs on a CIF and writes a
   PNG and the numbers.

## 7. What is not claimed

Agreement with abTEM and the truths of section 5 do not make a simulated image match a micrograph: no frozen phonons, no
partial coherence, no detector, no inelastic scattering, no thermal diffuse scattering. The independent-atom model with
an infinite projection per slice is abTEM's default and an approximation (atoms near a slice boundary are projected whole
into one slice). abTEM's compensation function S (3.5) is ported as abTEM defines it; whether it is the right correction
for bilinear deltas is not judged here. The float32 Metal path is registered against the float64 CPU path only.

## 8. Licence

abTEM is GPL-3.0-or-later; the definitions of section 3 are abTEM's algorithms, ported, and the Lobato table is taken from
abTEM's data file (the paper's numbers). SimulatoKit is therefore GPL-3.0 (LICENSE, NOTICE crediting abTEM and Lobato &
Van Dyck); an app built on it must be GPL-compatible.

## 9. Amendments

### 2026-10-15 00:10, after the pin generator ran, before any Swift: what the generator found

`tools/slice-1/pins.py` re-derives every quantity of section 3 from this text in NumPy and asserts it equals abTEM's
float64 output. 310 self-checks passed; the largest differences between the text and abTEM: V 1.5e−14 × max, the
bandlimited t 2.2e−16, the incident probe 3.0e−16, exit waves 1.0e−14, patterns 6.6e−16, f and λ 0 (bit-equal), the k = 0
identity of truth 2 3.7e−16. The exact-π V differs from abTEM's by at most 5.7e−8 × max, inside the registered 1e−7.
The generator is deterministic (the arrays' sha256 equal on two runs). Three readings are settled here.

1. **Section 3.10, the crop, read in full.** `diffraction_patterns` has `parity="odd"` by default, and its cutoff size goes
   through two parity steps: n' = ⌊round((2/3) / max(Δ) × L, 7)⌋ is first moved by +1 to the grid's own parity
   (`_antialias_cutoff_gpts`, parity "same"), then by +1 again if it is even (parity "odd"). So the pattern is always odd,
   with the zero frequency at its centre after `fftshift`; the crop keeps the first ⌊n'/2⌋ + 1 and the last ⌊n'/2⌋
   frequencies (`_fft_interpolation_masks_1d`). Example: 128 px of 0.078125 Å gives 85.33 → 85 → 86 → 87. Section 3.10's
   "adjusted to n's parity by +1" was half of the rule. Also section 3.7: abTEM's integer `exit_planes = m` includes the
   entrance plane, so the planes are (incident, after slice m, 2m, ..., n) and B1's ten cells give eleven planes.
2. **Truth 3 is replaced by a three-part form.** Measured, the FFT method's error against the analytic V(r) is an absolute
   aliasing error of order 1e−2 V·Å (the f ∝ 1/k² tail beyond the Nyquist frequency is missing), which grows relative to
   V(r) as V(r) decays (27 % at 4 Å on the 512² grid), so "1e−2 relative for 1 ≤ r ≤ 4 Å" was the wrong form for a true
   statement. The truth is now: (a) on 512² at 0.05 Å, within 2e−2 relative for 0.25 ≤ r ≤ 2 Å (measured 9.5e−3); (b) within
   3e−2 V·Å absolute for 2 < r ≤ 4 Å (measured 1.8e−2; this bound was set after the measurement, so it is weak on its own);
   (c) halving the pixel (1024², 0.025 Å) at least halves the largest absolute error over 0.25 ≤ r ≤ 4 Å (measured ratio
   0.28). (c) carries the weight: it is the method converging to the closed form, which abTEM's code cannot supply.
3. **What the fixtures hold.** To keep the repository small (13.3 MB instead of 30), the stored arrays are: V of slices 0
   and 1 and the bandlimited t of slice 0 for A1, A2, B1 and D1 (B2 and C1 reuse B1's and vacuum's, recomputed); exit waves at
   the last plane (A1, A2, D1 their only plane); B1's patterns at planes 1, 5 and 10 (cells 1, 5, 10); every other plane is
   pinned by its scalars (Σ|ψ|², Σ pattern, max |ψ|, the pattern's centre). Intensities are not stored: they are |ψ|² of
   the stored exit waves. Grids and cases are unchanged.

### 2026-10-15 00:20, at the first full test run: truth 5 corrected, and the Metal path's total intensity (Gate D)

**Truth 5 was wrong physics, and is corrected.** The registration said the projected potential of Al [001] has period
a/2, so every exit wave lives on the 2/a lattice. The projection has period a/2, but a 2.025 Å slice holds only half the
atoms: (0, 0) and (½, ½) at z = 0, (0, ½) and (½, 0) at z = a/2. Each slice is C-centred with translation (a/2, a/2), and
the propagator between the two kinds of slice makes the (h + k even, h and k odd) beams real: the dynamical coupling to
the first-order Laue zone ((1 1 ±1) is allowed). Measured: up to 4.9e−5 of the pattern's maximum at ten cells, growing with
thickness. The truth is now: (a) with slices of a/2, all intensity on (h, k) (in units of 1/a) with h + k even, off it
≤ 1e−20 × max; (b) with one slice per cell (thickness a, each slice the whole cell's projection, period a/2), all intensity
on h and k both even, off it ≤ 1e−20 × max; (c) the 4mm symmetry as registered (each slice is 4mm about the origin). The
port did not change; the test did.

**Metal: the total intensity drifts, cause proven (Gate D).** The registered bound for Σ|ψ|² (Metal against CPU,
1e−5 relative, predicted ≤ 3e−6) failed: 1.8e−5 on B2 (20 slices) and 8.4e−5 on the 512², 100-slice probe; the exit wave
(worst 5.0e−5 against 1e−4), the intensity (9.7e−5 against 2e−4) and the pattern (1.24e−4 against 2e−4) held. The
drift is linear in the number of slices, about 8e−7 per slice at 512². Gate D, `tools/metal-fft-gain/`: MPSGraph's
float32 FFT round trip, with no multiplication in between, loses 1.4e−7 to 1.4e−6 of Σ|ψ|² per round trip (by grid and
field), linearly over 100 round trips, for both of its scaling modes; the zero frequency passes exactly; after a fitted
scalar the remainder is 1e−13 of the energy. So the loss is a near-uniform gain of about 1 − 5e−7 on the non-zero
frequencies of each MPSGraph FFT pair, a property of the framework on this OS and device, not of the port (one step is one
round trip). The bound is amended to **1.5e−6 × (number of slices)** relative (the measured worst per round trip, a delta
at 256², is 1.43e−6), with the cause recorded. Not done tonight, and why: a scalar compensation would over-correct
DC-heavy fields (the plane wave's exit), and a frequency gain map or a compute-shader FFT with twiddles rounded from
float64 is an acceleration lane of its own. **The owner accepted this amended bound on 2026-10-15** (the board's waiting
list); the compute-shader FFT stays an optional acceleration lane, not a correction owed.

### 2026-10-15 01:40, Gate B: the refuter on the frozen commit `11bcdb4`, and the fixes it required

One Gate B refuter (a Sonnet subagent) compared the port with abTEM 1.0.10 (float64, as the pins) on 247 inputs the
fixtures do not hold: H to Lr (Z = 1 ... 103 sampled, mixtures of seven elements), 20 to 1000 keV, up to 306 slices,
grids from 8 to 960 points including 24 × 960 and anisotropic sampling, atoms on slice edges (±1e−12), at the top-face
threshold (±1e−13), negative, far outside the box, at 1e9 Å, probes at the corner and outside the box with C10 to ±1e6 Å
and every order to C56, semiangles 0.5 mrad to infinity, uneven slicings and explicit lists, rotated SrTiO₃ and Pnma
CaTiO₃ (720 atoms). **It could not break the registered agreement:** worst V 6.7e−14 × max (bound 1e−11), incident probe
1.65e−15 (1e−12), exit wave 1.3e−13 (1e−10), intensity 1.65e−13, pattern 7.3e−14 (1e−10); the exact π within its 1e−7 /
1e−6. Metal against CPU passed every case up to 102 slices (worst wave 4.4e−5, intensity 8.8e−5, pattern 1.59e−4 = 80 %
of its bound, Σ|ψ|² 9.5e−7 per slice = 63 %).

Fixed (tier 2, a test each, each test shown to fail without its fix):
1. A requested thickness so small that ⌈Lz / s⌉ overflows (1e−320) trapped in `Int(n)`; now `invalidSliceThickness`, with a
   ceiling of `Slicing.maximumSliceCount` = 10⁶ slices.
2. `Multislice.run` compared grids with exact `==`, so a grid whose extent was computed along another path (9.600000000000001
   against 9.6) was refused; now `Grid.matches` (equal counts, extents within 1e−12 relative).

Recorded, not changed (no number of the registration moves):
3. **Atoms past the last edge of a short explicit list.** A list may sum short of Lz within `isclose` (up to 1e−5 × Lz); an
   atom in that gap is dropped by abTEM (`label_to_index(..., max_label=n−1)`, `core/utils.py`) and put into the last slice
   by the port, as section 3.4 registered. Kept: the port never drops an atom of the box silently.
4. **abTEM drops an atom at a tiny negative x or y in `Potential.build()` but not in `multislice(potential)`**: `wrap(eps=0)`
   rounds a fractional −1e−20 to exactly 1.0, and `pad_atoms` → `atoms_in_cell` (`atoms.py` 954–956) requires < 1. The port
   matches abTEM's multislice path. An abTEM inconsistency, reported here; ASE builders produce such −1e−16 coordinates.
5. **`pbc`.** The port is periodic on all three axes (section 3.4); abTEM wraps only the axes ASE marks periodic and crops
   the rest. Agreement holds for `pbc=True`; an ASE `Atoms` built with `pbc=False` (ASE's default) differs.
6. **`ExitPlanes.afterSlices` always includes the last slice** (section 3.7 as ported); abTEM's tuple `(1, 3)` returns only
   those planes.
7. **Metal beyond the registered 100 slices.** At 304 slices (240², probe) the pattern error reaches 3.9e−4 against the 2e−4
   bound, growing about 1.3e−6 per slice (float32 rounding through 2n FFTs; the Σ|ψ|² gain stays inside its bound). The
   Metal bounds are registered for up to 100 slices; a longer run needs its own bound or a better FFT (the acceleration
   lane of amendment 2).
8. **Metal memory.** All transmission functions are held in float64 and again in float32: 960² × 100 slices peaked at
   3.8 GB. A performance lane (transmissions built per slice on the GPU) is planned.

The refuter's scripts, results and logs are in the session's scratchpad (not committed): 247 cases, `results.jsonl`.

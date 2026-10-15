# SimulatoKit status

## State (2026-10-15, afternoon: the Metal performance lane and the 0.1.0 candidate)

- **The Metal performance lane, part 2, is done** (`docs/registration-performance.md`, registered at `1b34eb6`, code
  `5300102`, memory fix and mutants `88b30da`): `Potential` (abTEM's lazy potential, `build()` bit-equal to the eager
  one); the potential built on the GPU from host-merged deltas; the transmissions on the GPU for both forms; one compiled
  MPSGraph step encoded per slice into one command buffer with ping-pong wave buffers; a per-step renormalisation of
  Σ|·|²; STEM probes and intensities on the GPU. Speed, structure to numbers, 100 slices: plane wave 6.6× / 33× / 190×,
  64 × 64 HAADF scan 31× / 28× / 37× at 256² / 512² / 1024²; device memory 1.0 to 1.7 GB. The gain-map attempt was
  rejected with its measurements (amendment 1); the renormalisation brought the measured Σ|ψ|² drift from 9.3e−7 to
  1.4e−7 per slice on Al and at most 5e−7 on heavy elements (slice 1's accepted 1.5e−6 stays the registered bound).
  Break-once 18 of 18 (`docs/performance-mutants-2026-10-15.md`). Gate B (amendment 3): about 450 cases, the bounds not
  broken on physical inputs (worst pattern 79 %); three defects fixed (no slice trapped, amplitude 1e8 gave NaN, the
  executable cache grew without bound). 0.0.1 candidate `bcde1b5`, gated on a clean clone.
- **0.0.1 candidate** at the lane's gated commit; **0.1.0 candidate** after it: `docs/parity.md` filled against abTEM's
  documented API (24 done, 34 partial, 64 absent), a README without the suite, the example (a CIF in, an HRTEM image and a
  HAADF scan out as PNGs and a CSV), `tools/consumer-check` (the whole public API from a separate package, negative
  controls), `.spi.yml`, `.github/workflows/ci.yml`, the hook refusing `pkill -f` and `killall`, `docs/consumers-0.0.1.md`
  and `docs/consumers-0.1.0.md` (`Measurement` renamed `ScanMeasurement`; `DiffractionIntensity.init(grid:beam:values:)`
  public). The owner tags; the commands are on the board row.
- **Platform floor**: declared macOS 27; the code uses nothing newer than macOS 15 (`Synchronization.Mutex` is the only
  macOS 15 API; MPSGraph's FFT needs 14).

## State (2026-10-15, end of the first session)

- **Repo created 2026-10-14/15 night** by the SimulatoKit session (the coordinator's side session). The owner created
  `github.com/mac4DSTEM/SimulatoKit` and pushed through `79dc601`, and made it public on 2026-10-15 morning; the later
  commits wait for the owner's push (board, waiting list). No tag.
- **Slice 3 is ported, refuted and green** (2026-10-15, owner's standing instruction): `ContrastTransferFunction`
  (abTEM's CTF: aberrations, spatial and temporal envelopes, aperture), `Waves.applying` for HRTEM images, Scherzer
  defocus, point resolution, Nyquist sampling. Registration: `docs/registration-slice-3.md` (three amendments). Pins: 35
  self-checks. Break-once 14 of 15 (one equivalent, `docs/slice-3-mutants-2026-10-15.md`). Gate B: about 6700 kernels, 738
  images and 2184 closed forms, worst 2.2e−15; one refusal fixed.
- **Performance lane, part 1** (`235d661`, tier 2): the CPU batch runs positions on at most 8 threads, bit for bit equal to
  sequential runs; 36 positions 0.068 s against 0.228 s.
- **Slice 2 is ported, refuted and green** (2026-10-15, owner's instruction of 00:55): frozen phonons (abTEM's
  displacements given, or the kit's seeded `GaussianGenerator`), `GridScan`, the STEM driver averaged over configurations,
  `AnnularDetector` (BF, ADF, HAADF), `PixelatedDetector` (4D-STEM), `SegmentedDetector` (DPC), the scan batched on Metal.
  Registration: `docs/registration-slice-2.md` (three amendments). Pins: 98 self-checks against abTEM. Break-once 22 of 23
  (one equivalent, `docs/slice-2-mutants-2026-10-15.md`). Gate B: 1843 comparisons passed; its trap, two convention
  differences and four refusals fixed (amendment 3).
- **Slice 1 is ported, refuted and green**: the projected potential (Lobato, abTEM's default infinite projection) and
  multislice for a plane wave and a convergent probe, on the CPU (float64, Accelerate) and on Metal (float32, MPSGraph),
  against abTEM 1.0.10 (`164e644`) in float64. Registration: `docs/registration-slice-1.md` (sections 1–8 and three
  amendments in section 9).
- **Gate** at `eb18296`: `tools/test.sh` green and stamped, 63 test functions in 14 suites, no warnings. (Slice 2's gate
  was at `e2afa7c`, 52 functions, cold build 0 warnings. Slice 1's gate was at `a61b736`: 37 test functions in 9 suites (Swift Testing, the kits'
  standard since 2026-10-15; parameterised per case and beam), the Metal suite run on the M5 Pro; a cold release build with
  0 warnings.) Slice 1 break-once: 40 of 40 mutants killed (`docs/slice-1-mutants-2026-10-15.md`). Gate B: the refuter could not
  break the agreement on 247 new inputs (amendment 3); its two fixes are in.
- **The name** passed the clash check on 2026-10-14 night: no Apple framework, no Swift Package Index package and no
  GitHub Swift repository is called SimulatoKit or Simulato (one JavaScript model-based testing tool, GannettDigital/simulato,
  is unrelated); no electron-microscopy software of that name.

## What agrees with what (measured)

| Gate | Bound | Measured |
|---|---|---|
| λ, σ, κ, the float32 π against abTEM | 1e−15 | equal on every pinned energy |
| f(k²) against abTEM, `.abTEMFloat32` | 1e−14 | passes (the generator's text form is bit-equal) |
| V per slice, `.abTEMFloat32` | 1e−11 × max | passes on A1, A2, B1, C1, D1 |
| exit waves, patterns, per-plane scalars, `.abTEMFloat32` | 1e−10 × max | passes on all eight case-beams |
| the same with `.exact` (the default) | 1e−6 × max (V 1e−7) | passes; V differs from abTEM's by ~5e−8 |
| truth 1–8 (Kirkland's table, the k = 0 identity, the far field and its convergence, vacuum, the slice lattice and 4mm, the probe, χ, the thickness series) | section 5 | all pass |
| Metal against CPU, wave / intensity / pattern | 1e−4 / 2e−4 / 2e−4 × max | worst 5.0e−5 / 9.7e−5 / 1.24e−4 (512², 100 slices); 8.5e−6 / 1.7e−5 / 4.0e−5 on the pinned cases |
| Metal against CPU, Σ\|ψ\|² | 1.5e−6 × slices (amendment 2) | 8.4e−7 to 9.2e−7 per slice |
| refuter, 247 new inputs, CPU `.abTEMFloat32` | as above | V 6.7e−14, probe 1.65e−15, exit 1.3e−13, pattern 7.3e−14 |
| refuter, Metal to 102 slices | as above | wave 4.4e−5, pattern 1.59e−4 (80 % of the bound), Σ\|ψ\|² 9.5e−7 per slice |
| slice 2: frozen phonons, every detector, `.abTEMFloat32` | V 1e−11, annular and pixelated 1e−10 × scale | ~1e−14 (pins and refuter) |
| slice 2: segmented against abTEM | 1e−6 (abTEM sums in float32) | 3.8e−7; against the float64 sum 2e−15 |
| slice 2: Metal batched scan against CPU | 2e−4 × scale | 8.6e−6 (refuter, 480², 225 positions) |
| slice 3: CTF kernels, HRTEM images, closed forms | 1e−12, 1e−10, 1e−15 | 3.7e−16, 2.2e−15, bit-equal (refuter) |
| performance: GPU-built V against the CPU's | 1e−5 × max V | 3.1e−7 (A1, A2, B1, D1) |
| performance: GPU probe against `Probe.build` | 1e−5 × max | 5.0e−7; Σ\|Ψ(k)\|² − 1 = −5e−7 |
| performance: Metal (renormalised) against CPU, 512² × 100 slices | slice 1's | wave 1.6e−5, pattern 6.5e−5, Σ\|ψ\|² 2.2e−8 to 1.4e−7 per slice |
| performance: refuter, about 450 cases, both forms | slices 1 and 2's | wave 4.6e−5, pattern 1.58e−4, scans 3.4e−5, Σ\|ψ\|² ≤ 4.9e−7 per slice |
| performance: the same at 1024² × 100 slices | slice 1's | wave 1.3e−5, pattern 3.9e−5, 8.4e−8 per slice |
| performance: GPU-built scans (S1, S3) against CPU | 2e−4 × scale | 4.9e−6; bit-equal on a second run |

## Known limits

- **MPSGraph's float32 FFT is not exactly unitary** (Gate D, amendment 2, `tools/metal-fft-gain/`): each round trip loses
  1.4e−7 to 1.4e−6 of Σ|ψ|², depending on the data (a fixed gain map does not correct it: performance amendment 1). The
  per-step renormalisation keeps the Metal path within 5e−7 per slice of the CPU path's (measured; 1.4e−7 on Al); the
  registered bound stays 1.5e−6 per slice, accepted by the owner. Re-measure after OS updates.
- **Weak high-angle signals on Metal** sit at float32's floor (about 3e−10 of the probe's intensity, absolute): a single H
  atom's 85–100 mrad image is off by 2.3e−4 of its own maximum (bound 2e−4), as it was before the performance lane
  (Gate D). A floor for the bound is proposed to the owner (performance amendment 3, item 6); use `.cpu` for such signals.
- **Unphysical phases** (σV above about 200 rad in one pixel column) break the Metal bounds.
- **The Metal bounds are registered to 100 slices** (grids to 1024² since the performance lane). Before the lane, at 304
  slices the pattern error reached 3.9e−4 (bound 2e−4); not re-measured with the renormalisation. A 64 × 64 scan at 1024²
  is GPU-bound at about 100 ms per position.
- **Differences from abTEM at the edges**, recorded in amendment 3 and `docs/parity.md`: always periodic (abTEM honours
  ASE's `pbc`), an atom in the gap of a short slice list kept (abTEM drops it), `.afterSlices` always adds the last slice,
  and abTEM's own `build()` / `multislice()` disagreement on atoms at a tiny negative x or y.
- **Orthogonal cells only; `AtomicStructure` is a local type** marked for replacement by DiffractoKit's cell types at
  DiffractoKit's first tag (no tag on 2026-10-15). The example's CIF reader likewise.
- **Partial coherence is the quasi-coherent envelope approximation** (abTEM's), not an incoherent focal or angular sum;
  no detector MTF or noise: a simulated image is still not a micrograph. Frozen phonons use Einstein-model σ given by the
  caller (no Debye–Waller table).
- **Frozen phonons agree with abTEM only when abTEM's displacements are given**; the kit's generator agrees in distribution.
- **Segmented detectors** sum in float64 where abTEM sums in float32: they agree with abTEM to ~2e−7.
- **A `ProjectedPotential` on `.metal`** uploads every slice's V (float32, 4 bytes per pixel per slice); a `Potential`
  uploads only the atoms' merged deltas, so prefer `Potential` for large Metal runs. CPU scans use up to 8 threads.

## Hand-off

Next (the board's `lanes/simulatokit` row holds the exact prompt): the owner's tags 0.0.1 and 0.1.0; then the kit loop's
asks from the apps and the consumer kits, if any, as 0.1.x; then slice 4, exit-wave and focal-series reconstruction
(abTEM has none: the reference is the literature's closed forms and a simulated focal series from this kit), or what
the apps ask first; the generality pass is 0.2.0.

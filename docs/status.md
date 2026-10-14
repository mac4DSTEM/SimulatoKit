# SimulatoKit status

## State (2026-10-15, end of the first session)

- **Repo created 2026-10-14/15 night** by the SimulatoKit session (the coordinator's side session). The owner created
  `github.com/mac4DSTEM/SimulatoKit` (private until 0.1.0) and pushed through `79dc601`; the later commits wait for the
  owner's push (board, waiting list). No tag.
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

## Known limits

- **MPSGraph's float32 FFT is not exactly unitary** (Gate D, amendment 2, `tools/metal-fft-gain/`): each round trip loses
  1.4e−7 to 1.4e−6 of Σ|ψ|² on the non-zero frequencies, so the Metal path's total intensity drifts ~8e−7 per slice below
  the CPU path's; the amended bound of 1.5e−6 per slice is accepted by the owner (2026-10-15). Re-measure after OS updates.
- **Metal is not yet faster.** At 512² and 102 slices the example took 2.81 s on the CPU and 2.35 s on Metal, both
  dominated by the potential, which the CPU path builds for both. One MPSGraph run per slice, with host round trips for
  the transmission functions, is the slice-1 design (correctness first); a performance lane (one executable for all
  slices, the potential on the GPU) is planned.
- **The Metal bounds hold to 100 slices.** At 304 slices the pattern error reaches 3.9e−4 (bound 2e−4), growing about
  1.3e−6 per slice; Metal memory holds every transmission twice (960² × 100 slices peaked at 3.8 GB). Amendment 3.
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
- **The potential is built on the CPU** for both backends (the remaining performance lane); CPU scans use up to 8 threads.

## Hand-off

Next: the consumers and the 0.1.0 preparation (docs/consumers-0.1.0.md for mac4DSTEM, macTEM and the other kits' tests; DocC;
the public API pass; a consumer-check package like DiffractoKit's), and the rest of the Metal performance lane (the
potential on the GPU, one MPSGraph executable over all slices, measured against the CPU path, and an attempt to tighten
slice 1's amendment-2 intensity bound). The board's `lanes/simulatokit` row holds the exact prompt.

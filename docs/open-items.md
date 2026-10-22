# Open items

A plain list, most important first. Benchmarks: `docs/benchmarks.md`; coverage: `docs/parity.md`.

## Slower than Python at a typical size

1. **Any-length FFT off the vDSP lengths (Bluestein), 250² and 506²**: 3.8 ms against NumPy's 0.57 ms at 250², 20 ms against 9 ms
   at 506² (kit CPU float64 against abTEM float64, 2026-10-21, machine loaded). Cause: each line costs three power-of-two
   FFTs of 2n to 4n points; NumPy's FFT runs mixed radix on 250 = 2·5³ directly. Fix candidate: a mixed-radix plan for lengths
   that factor into 2, 3, 5 and 7, with Bluestein for the rest, or the smallest f · 2^m ≥ 2n − 1 (f in 1, 3, 5, 15) as the
   padded length. Numbers are right (rtol 1e−7); only speed is open.
2. **vDSP lengths between 186² and 256², `fft2+ifft2`**: the kit is 0.3 to 0.5 of NumPy at 192², 256² and 512² in the loaded
   run (1.25 ms against 0.35 ms at 192²) but equal in a quiet run (1.06 ms against 1.04 ms at 256²). Re-measure in the nightly
   review before changing anything (diagnose first, `../CLAUDE.md`); the column gather of `FFT2D.transform` is the suspect.
3. **Multislice and frozen phonons at 192², kit CPU float64 against abTEM float64**: 0.79x and 0.90x; neighbours (128², 256²) are
   faster than abTEM, so this looks like the load of the run. Re-measure in the nightly review.
4. **Probe build at 192² to 512², best kit against abTEM float32**: 0.67x to 0.83x. abTEM float32 computes the aberration
   function in float32; the kit is float64 only. A Metal or float32 CPU path for the probe is the work item if it matters;
   a scan builds its probes on the GPU already.

## Work items by priority

- **Study, next:** the study rows still missing in `docs/parity.md` (the 0.4.0 rows `LineScan`, integrated centre of mass and `WavesDetector` are done).
- **Study, with DiffractoKit:** non-orthogonal and rotated cells, `periodic=False`, `orthogonalize_cell` and
  `standardize_cell` through `CrystalStructure` and `CIFImport`, replacing the local `AtomicStructure` stand-in. Finite
  projection, the Kirkland, Peng and Ewald parametrisations, Debye–Waller tables.
- **Study, last:** the public `PotentialSlice` initialiser, `crystal_slice_thicknesses`, `PotentialArray.project`, tiled potentials.
- **Metal:** the Metal bounds hold to 100 slices and grids to 1024²; re-measure past 300 slices with the renormalisation
  (ROADMAP); the benchmark table has no 1024² row yet (nightly).
- **Poisson noise:** `samples > 1` and the NumPy random stream are not offered (the draws are the kit's own).
- **Tilt ensembles:** a distribution of tilts (abTEM's tilt axis) is not offered; loop over tilts. A `.metal` batch needs one tilt.

0.4.0 notes (line scans, iCOM, waves):

- **Fast test target for SimulatoDetector:** `Package.swift` has none; the scan, refusal and Metal tests of `LineScan`, `WavesDetector` and iCOM live in `SimulatoKitParityTests` (skipped by the per-session gate). Add `SimulatoDetectorTests` (dependencies SimulatoDetector, SimulatoMultislice, SimulatoProbe, SimulatoCore, SimulatoKitTestSupport) and move `memoryRefusalIsTyped` and `frozenPhononsAreRefused` there.
- **WavesDetector with frozen phonons:** one potential per call; abTEM's `ensemble_mean=False` for waves is a loop over `STEM.scan` calls, no batch API. abTEM 1.0.10's own `probe.scan(..., WavesDetector(gpts=...))` raises a broadcast error (the detector does not declare its output shape); the pins use `probe.multislice` then `WavesDetector(gpts).detect`, abTEM's per-wave code.
- **Line profiles:** `ScanMeasurement.lineProfiles(for:)` takes detectors of one number or a vector; abTEM returns other measurement types for segmented and pixelated line scans, which the kit leaves as the raw `ScanMeasurement`. `ReciprocalSpaceLineProfiles` has no producer in a scan yet (the azimuthal average is not ported); `interpolateReciprocalLine` produces them. `LineScan` has no `end=None` or probe-matched sampling default: pass `end` and `gpts` or `sampling`.
- **iCOM:** no calibration beyond abTEM's (the image is in units of the centre of mass times Å, minimum 0); a phase or potential scaling is not in abTEM and not built; `integrate_disc(interpolate:)` and the DPC measurement of a `SegmentedDetector` are not built.
- **Lobato π:** the pins of 0.4.0 use `LobatoParameterization(piConstant: .abTEMFloat32)` like the other slices; with the default `.exact` they differ by 2.6e-8 of the maximum.

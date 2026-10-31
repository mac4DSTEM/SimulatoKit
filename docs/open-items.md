# Open items

A plain list, most important first. Each item has an issue; the issue holds the detail. Benchmarks: `docs/benchmarks.md`;
coverage: `docs/parity.md`.

## Slower than Python at a typical size

- **Any-length FFT with a prime factor above 7 (Bluestein)**: lengths off vDSP that factor into 2, 3, 5 and 7 run a
  mixed-radix FFT since #42 (250² 0.36 ms against NumPy 0.24 ms, 500² 1.58 ms against 1.52 ms); 506² (2 × 11 × 23) stays on
  Bluestein, 1.6x slower than NumPy. Numbers are right (rtol 1e−7), only speed is open. Fix candidate: generic odd-radix
  passes (#4).
- **Rows measured on a loaded machine**: `fft2+ifft2` at 192² to 512², multislice and frozen phonons at 192², the probe build
  at 192² to 512² against abTEM float32. Re-measure in a quiet run before changing anything (#5).

## Work items

- **Detector and noise follow-ups from 0.4.0**: line profiles for segmented and
  pixelated line scans, `LineScan` defaults, the "spline" method of `Images.interpolate`, tilt distributions (#9).
  `integrate_disc(interpolate:)` and `Images.interpolate` (method "fft") are in; the DPC measurement of a `SegmentedDetector`
  is under Polar measurements. The azimuthal average of a single pattern is in (`RealField.azimuthalAverage`), over a
  pixelated scan it is not.
- **Potential API**: `CrystalPotential` (a unit potential repeated, frozen phonons per cell), its API undecided (#7).
- **Cells**: non-orthogonal and rotated cells, `periodic=False`, `orthogonalize_cell` in `AtomicStructure` (the kit keeps
  its own structure type: `docs/adr/0001-standalone-structure-types.md`). Finite projection, the Kirkland, Peng and Ewald
  parametrisations, Debye–Waller tables (`docs/parity.md`).
- **Metal**: the bounds hold to 100 slices and grids to 1024²; re-measure past 300 slices with the renormalisation, and add
  a 1024² row to the benchmark table (#5).
- **Polar measurements**: the DPC measurement of a `SegmentedDetector` is in (`integrate`, `integrateRadial`,
  `differentials`); `PolarMeasurements.to_diffraction_patterns` and `azimuthal_limits` remain (#9).

## Notes

- **Lobato π**: the pins use `LobatoParameterization(piConstant: .abTEMFloat32)`; with the default `.exact` results differ
  by about 2.6e−8 of the maximum.
- **`WavesDetector` pins**: abTEM 1.0.10's own `probe.scan(..., WavesDetector(gpts=...))` raises a broadcast error (the
  detector does not declare its output shape); the pins use `probe.multislice` then `WavesDetector(gpts).detect`, abTEM's
  per-wave code.

# Open items

A plain list, most important first. Each item has an issue; the issue holds the detail. Benchmarks: `docs/benchmarks.md`;
coverage: `docs/parity.md`.

## Slower than Python at a typical size

- **Any-length FFT off the vDSP lengths (Bluestein)**: 2x to 7x slower than NumPy at 250² and 506²; numbers are right
  (rtol 1e−7), only speed is open. Fix candidate: a mixed-radix plan for lengths that factor into 2, 3, 5 and 7, Bluestein
  for the rest (#4).
- **Rows measured on a loaded machine**: `fft2+ifft2` at 192² to 512², multislice and frozen phonons at 192², the probe build
  at 192² to 512² against abTEM float32. Re-measure in a quiet run before changing anything (#5).

## Work items

- **Detector and noise follow-ups from 0.4.0**: `WavesDetector` with frozen phonons, line profiles for segmented and
  pixelated line scans, the azimuthal average, `LineScan` defaults, `integrate_disc(interpolate:)`, the DPC measurement of a
  `SegmentedDetector`, Poisson noise with `samples > 1`, tilt distributions (#9).
- **Potential API**: the public `PotentialSlice` initialiser, `crystal_slice_thicknesses`, `PotentialArray.project`, tiled
  potentials (#7).
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

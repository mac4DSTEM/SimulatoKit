# Consumers of SimulatoKit 0.1.0

SimulatoKit is the truth source for other packages' tests: simulated exit waves, diffraction patterns, 4D-STEM cubes and
images with known inputs, computed by a port of abTEM 1.0.10 that is pinned against abTEM in float64. This note says what
a consumer imports, which numbers it can rely on, and how to pin them.

## Pinning

```swift
.package(url: "https://github.com/mac4DSTEM/SimulatoKit", exact: "0.1.0")
```

Products: `SimulatoCore` (beam, grid, structure, fields, waves, the CPU FFT), `SimulatoMultislice` (`Potential`,
`ProjectedPotential`, `Multislice`, `FrozenPhonons`), `SimulatoProbe` (`Probe`, `Aperture`, `Aberrations`, `GridScan`,
`ContrastTransferFunction`), `SimulatoDetector` (`STEM`, `AnnularDetector`, `PixelatedDetector`, `SegmentedDetector`,
`DiffractionIntensity`, `ScanMeasurement`). GPL-3.0-or-later: a package that links it is GPL-compatible; a package that only
uses it in its test target and ships none of it should still state the licence of its test dependencies.

**Use the CPU path for test truth.** `.cpu` is float64 and agrees with abTEM to 1e−10 of the array's maximum or better
(exit waves, patterns, detectors; `docs/status.md`). `.metal` is float32 and agrees with `.cpu` within 1e−4 (waves) and
2e−4 (intensities, patterns, detector images) of the maximum, up to 100 slices; use it for speed, not as a reference.
Fix the parameterisation: `LobatoParameterization(piConstant: .abTEMFloat32)` reproduces abTEM's numbers bit-closely;
the default `.exact` differs from abTEM by about 5e−8 relative.

**Store what you compare.** A consumer test should store the simulated numbers it uses (a fixture with the SimulatoKit
tag and the inputs), or regenerate them at a pinned tag, so a SimulatoKit bump never moves a consumer's expected value
silently. A bump is a reviewed change: re-run, diff the fixture, record the difference.

## For a diffraction package (cells, scattering, pattern geometry, indexing)

- Zone-axis patterns: `Waves.planeWave` through a `Potential` of the oriented crystal (orthogonal cells; rotate the
  structure before slicing), `Multislice.run`, then `Waves.diffractionPattern(maxAngle:)` (abTEM's
  `diffraction_patterns`, odd-sized, zero frequency at the centre, sampled at 1/extent). Spot positions sit on the
  reciprocal lattice of the projected cell exactly (truth 5 of slice 1); intensities are dynamical, not kinematic, so
  compare positions and symmetry, and compare intensities only to a kinematic model in the thin limit (one thin slice,
  small σV).
- Convergent-beam patterns: a `Probe` (semiangle, aberrations) through the crystal and `PixelatedDetector` or the
  pattern of the exit wave.
- What it cannot give yet: non-orthogonal cells (until the cell types are shared with a crystallography package), a
  detector MTF or noise, inelastic scattering.

## For a 4D-STEM package (virtual images, centre of mass and DPC, strain, ptychography)

- A cube: `STEM.scan(probe:potentials:scan:detectors:backend:)` with `PixelatedDetector(maxAngle: .cutoff)` (or
  `.valid`, `.full`) gives one diffraction pattern per scan position, row-major over (scan x, scan y, kx, ky) in
  `ScanMeasurement.values`, with `ScanMeasurement.scanShape` and `valueShape`; averaged over frozen-phonon configurations when
  several potentials are given (`FrozenPhonons.lazyPotentials`).
- Virtual-image truth: `AnnularDetector(inner:outer:)` on the same scan is abTEM's annular integration of the same
  patterns (half-open masks: BF + DF = the total exactly, slice 2 truth 4); a virtual image computed from the cube with
  the same masks must equal it to rounding.
- DPC truth: `SegmentedDetector(radialBins:azimuthalBins:inner:outer:rotation:)`; in vacuum the four segments of a
  centred aberration-free probe are equal (truth 7), and the centre of mass of a vacuum pattern is the centre.
- Strain truth: build the strained structure explicitly (scale `AtomicStructure.box` and the positions) and scan it;
  the Bragg disks move by the known reciprocal strain.
- Scale: 64 × 64 positions on 256² to 1024² run on `.metal` in seconds to minutes (`docs/status.md`, the speed table);
  `.cpu` runs positions on up to 8 threads.

## For an app

Use the products directly; the example executable (`Examples/SimulatoExample`) is the minimal end-to-end use: a CIF in,
an HRTEM image and a HAADF scan out. `AtomicStructure` is a local, orthogonal-cell type that will be replaced by a shared
crystal type in a later minor version (a breaking change, announced in that version's consumer note).

## API changes from 0.0.1

- **Breaking:** `Measurement` is renamed `ScanMeasurement` (it collided with Foundation's `Measurement` in any file
  importing both; the consumer check found it). Replace the name; nothing else changes.
- `STEM.scan`, `Multislice.run` and `Multislice.runBatch` take either `[ProjectedPotential]` / `ProjectedPotential` or
  `[Potential]` / `Potential`; an empty array literal now needs its type (`[] as [ProjectedPotential]`).
- `DiffractionIntensity.init(grid:beam:values:)` is public (a consumer can test its own `Detector` on a known intensity).
- Everything else is documentation, CI, the example and the consumer check.

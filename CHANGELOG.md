# Changelog

All notable changes to SimulatoKit. Each release's Migration paragraph is the record for consumers.

No release dates are given: the repository has no git tags yet, so a date for each version could not be verified.
Anything not stated here, in the README or in the git history is left out.

## 0.4.0

The study release: line scans and their profiles, integrated centre of mass, and exit waves per scan position (suite task 23). Pinned against abTEM 1.0.10 float64 (`tools/linescan/pins.py`, `tools/icom/pins.py`, `tools/waves/pins.py`); `docs/parity.md` goes to 37 present, 36 partial, 55 missing of 128 rows (from 28, 36, 59 of 123).

### Added
- `LineScan(start:end:gpts:sampling:endpoint:fractional:extent:)` with `addingMargin` and `atPosition` (abTEM `LineScan`, `.add_margin`, `.at_position`); `STEM.scan` takes any `ScanPattern` (`GridScan` or `LineScan`) and returns `scanShape` `[n]` for a line.
- `RealSpaceLineProfiles` and `ReciprocalSpaceLineProfiles` with `width(height:)` and `interpolate(sampling:gpts:order:endpoint:)`; `ScanMeasurement.lineProfiles(for:)` turns a line-scan measurement into profiles.
- `RealField.interpolateLine`, `interpolateLineAtPosition` and, for centred diffraction patterns, `interpolateReciprocalLine`, `interpolateReciprocalLineAtPosition` (abTEM `interpolate_line`, `interpolate_line_at_position`: spline orders 0 to 5, width, margin, fractional coordinates).
- `RealField.integratingGradient(x:y:)` and `ComplexField.integratedGradient(samplingX:samplingY:)`: abTEM's `Images.integrate_gradient`, the Fourier-space integral of a two-component gradient field (NumPy `fftfreq`, zero frequency set to 0, no padding, minimum shifted to 0), CPU float64 on the kit's FFT, any length from 2.
- `ScanMeasurement.integratedCenterOfMass(samplingX:samplingY:)` and `gradientImages(...)`: the integrated centre of mass (iCOM / DPC) of a `CenterOfMassDetector` scan, as abTEM's `integrated_center_of_mass`.
- `RealField.integrateDisc(x:y:radius:mean:border:)` (abTEM `integrate_disc`) and `RealField.tiled(repetitionsX:repetitionsY:)` (`Images.tile`).
- `ImageOperationError`.
- `WavesDetector` (abTEM's `WavesDetector(gpts)`): `STEM.scan(probe:potentials:scan:detector:backend:options:memoryLimit:)`
  returns a `ScanWaves`, the complex exit wave at every scan position (`[position][row][column]`, subscript `scanWaves[ix, iy]`
  gives a `Waves`). `WavesDetector(rows:columns:)` resamples each wave in Fourier space (crop or zero-pad over the same
  extent), as abTEM's `fft_interpolate`. Pinned against abTEM 1.0.10 float64 at 1.8e-15 of the wave maximum (assert_allclose
  rtol 1e-7, atol 1e-12 of the maximum); `.metal` (float32) holds `MetalAccuracy.registered.wave` (measured 3.6e-7).
- `Waves.resampled(rows:columns:)` (SimulatoCore): abTEM's `fft_interpolate` for a wave, odd and even sizes, crop and pad.
- `WavesDetectorError`: `tooLarge` (positions x rows x columns x 16 bytes above half of the physical memory, or the
  `memoryLimit` given) and `severalPotentials` (frozen phonons).

### Migration
`STEM.scan(scan:)` is now generic over `ScanPattern`; a caller passing a `GridScan` needs no change. `ScanMeasurement.scanShape` can have one entry for a line scan: code that indexes `scanShape[1]` must check `scanShape.count` first, and the `(ix, iy)` subscript and `image(...)` need two scan axes. `STEM.scan` with a `WavesDetector` returns `ScanWaves`, a different type from the measurement of the other detectors (one potential only: a frozen-phonon list is refused with `WavesDetectorError.severalPotentials`; scan each configuration alone). A review of the new surface turned inputs that trapped into throws before the tag: `SimulatoError.gridTooLarge(rows:columns:required:limit:)` is a new case (a consumer that switches exhaustively over `SimulatoError` adds it; `Waves.resampled` throws it above half of the physical memory), `LineScan` caps explicit `gpts` below 1e9, `LineProfiles.interpolate` caps its points and the profile count times points at 2^31, `RealField.tiled` caps the result at 2^28 elements, `integrateDisc` refuses a disc or position beyond what an `Int` holds, and `ScanMeasurement.integratedCenterOfMass` and `gradientImages` refuse a one-axis scan. Everything else is additive.

## 0.3.1

Documentation only; no source change and no number moves.

- The 0.3.0 entry no longer claims Bluestein runs at NumPy's speed; it is 2x to 7x slower (`docs/open-items.md` 1), as measured.
- README: install line `from: "0.3.1"` (was `"0.1.0"`); the "no noise" sentence now names `ScanMeasurement.poissonNoise`.
- The 0.3.0 Migration paragraph now lists the defaulted `tilt:` argument on `Probe.init`, `Waves.planeWave` and `Multislice.propagator`, not only `Waves.init`.
- Added `docs/status.md`.

**Migration.** None.

## 0.3.0

- Any grid size: `Grid`, `Potential` and every FFT accept any length of at least 2, on `.cpu` and `.metal` (before:
  f · 2^m with f in 1, 3, 5, 15 and m ≥ 3). vDSP's DFT runs where it has the length, Bluestein's algorithm (own
  implementation; ImagoKit's FFT is not a dependency, the kits stay standalone) elsewhere. Bluestein is slower than NumPy's
  mixed-radix FFT: 2x to 7x on the sizes measured (3.8 ms against 0.57 ms at 250², `docs/open-items.md` 1); the numbers are right, only speed is open. Added `Grid.hasDirectTransform(_:)`. A 90 x 75 grid is pinned against abTEM (rtol 1e−7, plane wave and
  probe).
- Added `Potential(structure:sampling:sliceThickness:)` (a sampling or an (x, y) pair), `Grid(extentX:extentY:samplingX:samplingY:)`
  and `GridScan(start:end:samplingX:samplingY:endpoint:)`: ⌈extent / sampling⌉ points (one more with `endpoint`), as abTEM.
- Added beam tilt [mrad]: `Waves.tilt`, `Waves.planeWave(…tilt:)`, `Probe(…tilt:)` and `Multislice.propagator(…tilt:)`. As in
  abTEM it enters the Fresnel propagator only; exit waves keep it. Pinned against abTEM (plane wave and probe); `.metal`
  matches `.cpu` at the registered bound. A `.metal` batch must share one tilt (`SimulatoError.invalidTilt`).
- Added `CenterOfMassDetector`, `ScanMeasurement.gaussianSourceSize(…)`, `ScanMeasurement.poissonNoise(…)` and a public
  `ScanMeasurement.init`. The Poisson draw uses the kit's generator, not NumPy's stream.
- `docs/parity.md` rewritten from abTEM's source (28 present, 36 partial, 59 missing of 123 rows, each marked study, useful or
  not planned); `docs/benchmarks.md` and `tools/benchmarks/` added; `docs/open-items.md` added.
- Added `SimulatoError.invalidSampling` and `.invalidTilt`; `.unsupportedGridLength` is now thrown only for a length below 2.
- No number moves for a grid or call that worked in 0.2.1.

**Migration.** Nothing is required. `SimulatoError` gained two cases, so an exhaustive `switch` over it needs a `default` or the
two new cases. `Grid.isSupportedLength(_:)` is now true for every n ≥ 2: code that used it to pick vDSP lengths should call
`Grid.hasDirectTransform(_:)`. The `Waves` initialiser has a new defaulted `tilt:` argument, so a type that mirrors its
signature must add it. The same defaulted `tilt:` (zero) is new on `Probe.init`, `Waves.planeWave` and `Multislice.propagator`: a call
that worked in 0.2.1 compiles and gives the same numbers, and a type that mirrors any of those signatures must add it.

## 0.2.1

- Added `MetalAccuracy.detectorFloor` (1e-5) and `MetalAccuracy.detectorBound(imageMaximum:incidentIntensity:)`: the
  bound for Metal detector images is 2e-4 x max(the image's maximum, 1e-5 x the incident intensity).
- No number moves.

## 0.2.0

- Every choice the kit makes is now a parameter, and every default is the previous value. On the defaults no number
  moves: 14 bit digests of both backends and both potential forms are identical on 0.1.0 and 0.2.0, and every
  fixture is unchanged.
- Added `ExecutionOptions` (`options:` on `Multislice.run`, `Multislice.runBatch` and `STEM.scan`),
  `MetalAccuracy.registered`, `Grid.matches(_:relativeTolerance:)`, `maximumSliceCount:` on the slicing functions and
  both `Potential` initialisers, and the constants of `ContrastTransferFunction.scherzerDefocus`, `pointResolution`
  and `nyquistSampling` as parameters.
- Added a DocC catalogue with one article per module.
- Fixed an invalid `import a, b` line in the README quick start.
- Lowered the platform floor to macOS 15 and iOS 18 (was macOS 27). Numbers are validated on macOS only; no test runs
  on iOS yet.

**Migration.** Nothing, unless you want another choice: call sites written against 0.1.0 compile unchanged. One
source break: an unapplied reference to a changed function no longer has its old type
(`grids.allSatisfy(grid.matches)` fails to type-check; write `grids.allSatisfy { grid.matches($0) }`), and a protocol
requirement that mirrors an old signature must add the new parameter.

## 0.1.0

- The first API meant for use by other packages; added a README that stands on its own, `docs/parity.md` (abTEM's
  public API against SimulatoKit), a consumer-check package that uses the whole public API from outside, CI and
  the Swift Package Index configuration.
- `DiffractionIntensity.init(grid:beam:values:)` is now public.
- `STEM.scan`, `Multislice.run` and `Multislice.runBatch` accept either `[ProjectedPotential]` / `ProjectedPotential`
  or `[Potential]` / `Potential`.
- Rewrote the example executable as a first-ten-minutes walkthrough (a CIF in; an HRTEM image and a STEM HAADF scan
  out).

**Migration.** `Measurement` is renamed `ScanMeasurement` (it collided with Foundation's `Measurement`): replace the
name. An empty array literal passed to the functions above now needs its type (`[] as [ProjectedPotential]`).

## 0.0.1

- The first version: potentials (Lobato parameterisation, infinite projection), multislice on `.cpu` (float64) and
  `.metal` (float32), probes and aberrations, frozen phonons, scan grids, annular, pixelated and segmented detectors,
  the STEM driver, and the contrast transfer function, checked against abTEM 1.0.10 in float64.
- Added `Potential` (the sliced structure, built lazily; on `.metal` it is built on the GPU) and the Metal performance
  work. `.metal` now renormalises every step, so its numbers moved within their documented bounds (the total-intensity
  error per slice from about 9e-7 to at most 1.4e-7); `.cpu` numbers did not move.

Unverified: the notes for 0.0.1 and 0.1.0 describe changes relative to the earlier untagged `main`; the exact
commit of each version could not be confirmed without tags.

Consumer notes for 0.0.1 to 0.2.1 are archived.

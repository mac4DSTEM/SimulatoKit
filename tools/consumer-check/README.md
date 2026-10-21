# SimulatoKit consumer check

A separate SwiftPM package (`tools 6.2`, macOS 27, Swift 6) that depends on SimulatoKit by path (`../..`) and uses, from
outside the package, the public API of all four modules: SimulatoCore, SimulatoMultislice, SimulatoProbe and
SimulatoDetector. `package` access does not cross a package boundary, so the executable compiles only while every symbol it
names is `public`: a `public` reverted to `package` by mistake turns the build red here, which the kit's own tests (inside
the package, where `package` is visible) never notice.

## What it covers

326 checks (the count printed by `run.sh`), about 0.1 s to run, under 20 s to build in release.

- **SimulatoCore**: `ElectronBeam` (λ at 200 keV = 0.025079340317328468 Å as the kit computes it, against 12.2643/√(V(1+0.97848e-6 V));
  σ against Kirkland's closed form; the mass factor), `PhysicalConstants`, all 17 `SimulatoError` cases, `Grid` (`fftfreq` by hand,
  `isSupportedLength`, `matches`, the refusals), `Antialias` (the aperture counted pixel by pixel, `cutoffPatternSize` 21 → 22 → 23),
  `ComplexField`, `RealField`, `Waves` (plane wave stays 1, both `Normalization`s and `MaxAngle`s), `AtomicStructure` (ASE order of
  `repeated`).
- **SimulatoMultislice**: `LobatoParameterization` (f(0) = 16π Σ A/B³, both `PiConstant`s), `Slicing`, `ProjectedPotential` (Σ V = f(0)/ΔxΔy,
  periodic wrap, every init), `ExitPlanes` (`every(m)` plane indices), `Multislice` (`run`, `runBatch`, `transmissionFunctions`,
  `propagator`, `maximumCPUThreads`, a vacuum plane wave stays 1, batch = single runs, exact today; D12 allows a tolerance), `Backend.metal` only if
  `Multislice.metalIsAvailable` (otherwise a SKIP line with the reason, and the `metalUnavailable` refusal is checked),
  `FrozenPhonons` (σ = 0 equals the static potential, directions, seeds, statistics, every refusal), `GaussianGenerator`.
- **SimulatoProbe**: `Aperture`, `Aberrations` (`chi` closed forms, `dChiDAlpha` and `dChiDPhi` against central differences),
  `Probe` (Σ_k|Ψ(k)|² = 1, the peak at its position), `GridScan` (`positions`, `linspace`), `ContrastTransferFunction` (Scherzer
  defocus, point resolution, Nyquist sampling, both envelopes, the kernel, `Waves.applying`).
- **SimulatoDetector**: `DiffractionIntensity`, `AnnularDetector` (over all angles equals the total; BF + ADF = total),
  `PixelatedDetector` (all three `MaxAngle`s), `SegmentedDetector`, the public `Detector` protocol adopted by two detectors defined
  in the consumer, `STEM.scan` with frozen-phonon averaging (σ = 0 equals the static scan), `ScanMeasurement` (`subscript`, `image`).

Each check prints `PASS` or `FAIL`; a refusal is checked for the case it throws. The process exits non-zero on any failure.

The section "Multislice: the lazy Potential" in `MultisliceChecks.swift` checks the lazy `Potential`, the `Potential` overloads of
`Multislice.run`, `runBatch` and `STEM.scan` (CPU, and Metal when a device exists), `FrozenPhonons.lazyPotentials` and
`PotentialSlice(thickness:values:)`; the detector section builds a `DiffractionIntensity(grid:beam:values:)` directly.

## Running it

    tools/consumer-check/run.sh

builds in **release** (`swift build -c release -j 4`; the kit is built as a consumer builds it, with no unsafe flags), requires a
warning-free build whose log names no test or example target, runs the checks, and prints `consumer-check: all passed — N checks ...`.

## The second consumer: the README's Quick start

    tools/consumer-check/readme.sh      # also the last step of run.sh

A second separate package, `readme/` (same identity trick, same four products), whose `main.swift` is **generated**: `readme.sh`
copies the first ```` ```swift ```` block after the `## Quick start` heading of `README.md` verbatim, then appends the fixed
`readme/checks.swift.txt`, which asserts on the README's own `hrtem`, `haadf` and `scan` (HRTEM image finite with mean in
(0.5, 1.5); HAADF values finite, non-negative, not all equal, `count == scan.nx * scan.ny`). It builds in release (`-j 4`, no
warnings allowed), runs, and prints `readme-consumer: the README's quick start builds and runs (N checks)`. If the quick start
stops defining `hrtem`, `haadf` or `scan` the script fails with a message saying so. The generated `main.swift` is gitignored
(`readme/.gitignore`); the README code is never edited by the script, so a quick start that stops compiling turns this red.

## Public declarations a consumer cannot use fully (by design or worth a look)

- `PotentialSlice` has no public initialiser: a consumer can read slices but build a `ProjectedPotential` from slices only by reusing
  another potential's `PotentialSlice` values.
- `DiffractionIntensity.init(grid:beam:values:)` is `package`: a custom `Detector` can be run through `STEM.scan` but not unit-tested on
  an intensity the consumer builds itself (only `DiffractionIntensity(_ waves:)` is public).
- `FFT2D`, `ASEUnits`, `Waves.cropIndices` and the Metal batch are `package`: intended.
- `Measurement` (SimulatoDetector) was ambiguous with `Foundation.Measurement` in a file that imports both: renamed `ScanMeasurement` for 0.1.0.

# Notes: WavesDetector (suite task 23, part 3)

For the session that merges CHANGELOG, parity.md, open-items and status. Nothing of those four was edited here.

## CHANGELOG bullets (Added)

- `WavesDetector` (abTEM's `WavesDetector(gpts)`): `STEM.scan(probe:potentials:scan:detector:backend:options:memoryLimit:)`
  returns a `ScanWaves`, the complex exit wave at every scan position (`[position][row][column]`, subscript `scanWaves[ix, iy]`
  gives a `Waves`). `WavesDetector(rows:columns:)` resamples each wave in Fourier space (crop or zero-pad over the same
  extent), as abTEM's `fft_interpolate`. Pinned against abTEM 1.0.10 float64 at 1.8e-15 of the wave maximum (assert_allclose
  rtol 1e-7, atol 1e-12 of the maximum); `.metal` (float32) holds `MetalAccuracy.registered.wave` (measured 3.6e-7).
- `Waves.resampled(rows:columns:)` (SimulatoCore): abTEM's `fft_interpolate` for a wave, odd and even sizes, crop and pad.
- `WavesDetectorError`: `tooLarge` (positions x rows x columns x 16 bytes above half of the physical memory, or the
  `memoryLimit` given) and `severalPotentials` (frozen phonons).

## parity.md rows

Row 160, same columns (abTEM | kit | state | priority):

| `WavesDetector(gpts)` | `WavesDetector(rows:columns:)`, `STEM.scan(probe:potentials:scan:detector:)` returning `ScanWaves` | **present**: the exit wave per position, `gpts` as a Fourier crop or pad; pinned at 1.8e-15 of the maximum on a 3 x 2 scan, six `gpts` (crop, pad, mixed, odd, even), `.metal` at 3.6e-7 against `.cpu`; one potential only (a frozen-phonon list is refused: averaging complex waves is not abTEM's `ensemble_mean` of intensities; scan each configuration alone); refuses above half of the physical memory; `to_cpu`, `url` do not apply | study |

New row (core, near the `Waves` rows):

| `fft_interpolate(array, new_shape)` for waves | `Waves.resampled(rows:columns:)` | **present**: crop or pad in Fourier space, `normalization="values"`; the 2D case only; pinned through the WavesDetector rows | study |

Row 94 (`Waves.multislice(potential, detectors)`): the "no `detectors=`" remark still holds; add "a `WavesDetector` over a scan is `STEM.scan(... detector:)`".

## open-items

Remove `WavesDetector` from "Study, next" (leaving `LineScan`, integrated centre of mass and DPC to the other sessions). Add: "**WavesDetector with frozen phonons:** one potential per call; abTEM's `ensemble_mean=False` for waves (one wave set per configuration) is a loop over `STEM.scan` calls, no batch API."
Also: "**Package.swift has no fast test target for SimulatoDetector:** the WavesDetector scan and refusal tests live in SimulatoKitParityTests (skipped by the per-session gate); add `SimulatoDetectorTests` (dependencies SimulatoDetector, SimulatoMultislice, SimulatoProbe, SimulatoCore, SimulatoKitTestSupport) and move `memoryRefusalIsTyped` and `frozenPhononsAreRefused` (they need no abTEM fixture beyond the cell) there."

## Deviations from abTEM

1. abTEM 1.0.10's `probe.scan(potential, scan, WavesDetector(gpts=...))` raises `ValueError: could not broadcast input array from shape (3,2,33,40) into shape (3,2,48,48)` (the detector declares no output shape, so the measurement is preallocated on the simulation grid). `tools/waves/pins.py` takes the exit waves of the scan (`probe.multislice`) and applies `WavesDetector(gpts).detect(waves)`, abTEM's per-wave code. Without `gpts` the scan works and is pinned as is. The kit supports `gpts` in the scan.
2. Frozen phonons: refused (`WavesDetectorError.severalPotentials`) instead of averaged; the other detectors average intensities, a mean of complex waves has no physical meaning.
3. The pins use `LobatoParameterization(piConstant: .abTEMFloat32)`: abTEM uses a float32 pi in its Lobato parameters even in float64. With the kit's default `.exact` the exit waves differ by 2.6e-8 of the maximum (cause diagnosed: the same difference as every other pin; with `.abTEMFloat32` it drops to 1e-15).
4. `to_cpu`, `url` are not offered.
5. `rows` and `columns` of at least 2 only (abTEM allows 1).

## Measured agreement (2026-10-22, CPU float64 against abTEM 1.0.10 float64, 3 x 2 scan, 48 x 48 grid, Al 1x1x2)

max |a - b| / max |b| per case: full 1.3e-15, (33,40) 1.2e-15, (81,96) 1.8e-15, (30,71) 1.5e-15, (64,25) 1.5e-15, (10,11) 9e-16.
Worst assert_allclose ratio |a - b| / (atol + rtol |b|): 4.4e-4 (limit 1).
Metal against CPU, max |difference| / max |wave|: full 3.6e-7, (30,71) 3.3e-7, (81,96) 3.5e-7 (registered 1e-4).

## New public names

`WavesDetector` (`init()`, `init(rows:columns:)`, `rows`, `columns`, `detect(_:)`, `requiredBytes(positions:rows:columns:)`, `defaultMemoryLimit`), `ScanWaves` (`scanShape`, `template`, `real`, `imaginary`, `positionCount`, `rows`, `columns`, `subscript(ix, iy)`), `WavesDetectorError`, `STEM.scan(probe:potentials:scan:detector:...)` (both `[ProjectedPotential]` and `[Potential]`), `Waves.resampled(rows:columns:)`; `package`: `Waves.resampleMap`.

## Files

`Sources/SimulatoCore/WavesResample.swift`, `Sources/SimulatoDetector/WavesDetector.swift`, `tools/waves/pins.py`, `Tests/SimulatoKitParityTests/Fixtures/waves-{pins.json,arrays.bin}`, `Tests/SimulatoKitParityTests/Detector/WavesDetectorTests.swift`, `Tests/SimulatoCoreTests/WavesResampleTests.swift`. No shared file edited.

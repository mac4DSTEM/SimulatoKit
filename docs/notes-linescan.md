# Notes for the merge session: LineScan, line profiles, interpolate_line (suite task 23, part 1)

Branch `worktree-agent-a81aa0277a1843f5b` (worktree `agent-a334cb8dfa7ab1467`), three commits on top of `c924afb`-era `main`. Delete this file after the merge.
CHANGELOG.md, docs/parity.md, docs/status.md, docs/open-items.md, README.md and Package.swift are untouched; their text is below.

## New public names

SimulatoProbe (`Sources/SimulatoProbe/LineScan.swift`):
- `protocol ScanPattern { positions: [SIMD2<Double>]; shape: [Int] }`, conformed by `GridScan` (`shape` = `[nx, ny]`) and `LineScan`.
- `struct LineScan` (`start`, `end`, `gpts`, `sampling`, `endpoint`, `extent`, `direction`, `angle` in rad, `shape`, `positions`),
  `init(start:end:gpts:sampling:endpoint:fractional:extent:) throws`, `addingMargin(_:)`, `addingMargin(before:after:)`,
  `LineScan.atPosition(center:extent:angle:gpts:sampling:endpoint:)`, `LineScan.endpoints(center:extent:angle:)`.

SimulatoDetector:
- `protocol LineProfiles` (`ensembleShape`, `count`, `values`, `sampling`, `start`, `end`, `lineWidth`, `extent`,
  `width(height:)`, `interpolate(sampling:gpts:order:endpoint:)`), `struct RealSpaceLineProfiles`, `struct ReciprocalSpaceLineProfiles`.
- `ScanMeasurement.lineProfiles(for: LineScan)`, `ScanMeasurement[i]` (one-axis scan).
- `RealField.interpolateLine(start:end:sampling:gpts:width:margin:order:endpoint:fractional:)`,
  `RealField.interpolateLineAtPosition(center:angle:extent:gpts:sampling:width:order:endpoint:)` (an image, Å),
  `RealField.interpolateReciprocalLine(...)`, `RealField.interpolateReciprocalLineAtPosition(...)` (a centred diffraction pattern, 1/Å,
  zero frequency at index n / 2, as `Waves.diffractionPattern` returns it).
- `STEM.scan(... scan: some ScanPattern ...)` in both overloads (was `GridScan`; source compatible).
- `ScanMeasurement.scanShape` may now have one axis; the `(ix, iy)` subscript and `image(...)` need two; `gaussianSourceSize` on a
  one-axis scan throws `invalidScan`.

Internal: `SplineInterpolation` (scipy's `map_coordinates` for orders 0 to 5, mode constant, prefilter mirror).

## CHANGELOG bullets (new top version, minor bump; suggest 0.4.0 together with the other agents' work)

- `LineScan(start:end:gpts:sampling:endpoint:fractional:extent:)` with `addingMargin` and `atPosition` (abTEM `LineScan`, `.add_margin`, `.at_position`); `STEM.scan` takes any `ScanPattern` (`GridScan` or `LineScan`) and returns `scanShape` `[n]` for a line.
- `RealSpaceLineProfiles` and `ReciprocalSpaceLineProfiles` with `width(height:)` and `interpolate(sampling:gpts:order:endpoint:)`; `ScanMeasurement.lineProfiles(for:)` turns a line-scan measurement into profiles.
- `RealField.interpolateLine`, `interpolateLineAtPosition` and, for centred diffraction patterns, `interpolateReciprocalLine`, `interpolateReciprocalLineAtPosition` (abTEM `interpolate_line`, `interpolate_line_at_position`: spline orders 0 to 5, width, margin, fractional coordinates).
- Migration: `STEM.scan(scan:)` is now generic over `ScanPattern`; callers passing a `GridScan` need no change. `ScanMeasurement.scanShape` can have one entry for a line scan; code that indexes `scanShape[1]` must check `scanShape.count` first.

## docs/parity.md rows (columns: abTEM | SimulatoKit | State | Priority)

Row 147 (Scans) becomes:
`| `LineScan(start, end, gpts, sampling, endpoint=True, fractional, potential)`, `.at_position`, `.add_margin` | `LineScan(start:end:gpts:sampling:endpoint:fractional:extent:)`, `.addingMargin`, `LineScan.atPosition`, `.positions`, `.shape` | **present**: pinned on 13 lines against abTEM (positions bit-equal, the end points to the last bit after a margin); a sampling wins over gpts and a margin recomputes the count end by end, as abTEM's setters do; `end` and one of `gpts`/`sampling` are required (abTEM's `end=None` and `match_probe` defaults not offered); `potential=` is the extent `extent:`; `STEM.scan` over a line: 5e-15 of the ADF with abTEM's float32 π, 1.1e-7 with the exact π, `.metal` 1.5e-6 of `.cpu` | study |`

Row 149 (`BaseScan.get_positions` ...) becomes: `| `BaseScan.get_positions`, `.limits`, `.num_positions`, `.shape` | `GridScan.positions`, `.nx`, `.ny`, `LineScan.positions`, `.gpts`, `ScanPattern.shape` | **partial**: positions and the shape; `limits` not offered | useful |`

Row 171 (Measurements) becomes (and `.metal` is covered by the line-scan row):
`| `RealSpaceLineProfiles`, `ReciprocalSpaceLineProfiles`, `_BaseMeasurement1D.width`, `.interpolate` | `RealSpaceLineProfiles`, `ReciprocalSpaceLineProfiles`, `LineProfiles.width(height:)`, `.interpolate(sampling:gpts:order:endpoint:)`, `ScanMeasurement.lineProfiles(for:)` | **present**: width to 4.5e-8 and interpolate to 5.5e-8 of the maximum against abTEM (abTEM stores both as float32; `assert_allclose` rtol 1e-7 holds at 0.6 of its bound), on the scan's profile and on planted profiles with ensemble axes, orders 0 to 5, both endpoints; `lineProfiles` only for detectors of one number or a vector (abTEM returns other measurement types for segmented and pixelated); `RealSpaceLineProfiles.tile`, `ReciprocalSpaceLineProfiles.angular_extent`, `show` and the lazy path not offered | useful |`

Row 190 (Measurements) becomes two rows:
`| `interpolate_line`, `interpolate_line_at_position` | `RealField.interpolateLine`, `.interpolateLineAtPosition` (images), `.interpolateReciprocalLine`, `.interpolateReciprocalLineAtPosition` (centred patterns) | **present**: 38 lines on a 24 × 18 image and a 21 × 20 pattern against abTEM, orders 0 to 5, width, margin, fractional, both endpoints: 3.5e-14 of the maximum, `assert_allclose` ratio below 1.2e-5 of its bound; `end` is required (abTEM's `(0, extent[0])` default not offered); a pattern is assumed centred (zero frequency at n / 2); complex and lazy arrays not offered; in SimulatoDetector, not SimulatoCore, because `LineScan` lives in SimulatoProbe | study |`
`| `integrate_disc`, `Images.integrate_gradient`, `Images.tile` | | **missing**: planned: SimulatoCore | study |`

Summary line of parity.md and status.md: three present rows more (LineScan, line profiles, interpolate_line) and one partial unchanged; adjust the counts after the other agents' rows.

## docs/open-items.md

- In "Work items by priority", remove `LineScan` from "Study, next" (leaving integrated centre of mass and DPC, `WavesDetector`).
- Add: **Line profiles from segmented and pixelated line scans**: `ScanMeasurement.lineProfiles(for:)` takes detectors of one number or a vector; abTEM returns `PolarMeasurements` and `DiffractionPatterns` ensembles for the others, which the kit leaves as the raw `ScanMeasurement`. `ReciprocalSpaceLineProfiles` has no producer in a scan yet (the azimuthal average is not ported); `interpolateReciprocalLine` produces them.
- Add: **`LineScan` defaults**: `end=None` and the probe-matched sampling (`match_probe`, 0.99 × Nyquist) are not offered; pass `end` and `gpts` or `sampling`.

## Measured agreement (abTEM 1.0.10, numpy 2.5.4, scipy 1.18.1, float64; `tools/linescan/pins.py`, deterministic: a second run is byte-identical)

| What | Measured | Bound used |
|---|---|---|
| `LineScan.positions`, 13 constructions (gpts, sampling, endpoint, both, one point, margin, asymmetric margin, `at_position`, fractional) | 0 (bit-equal); end points after a margin differ by one ulp (np.linalg.norm against sqrt(x² + y²)), tested at rtol 1e-14 | `assert_allclose` rtol 1e-7 |
| `STEM.scan` over two lines, ADF 50 to 150 mrad, Al 64² | 4.8e-15 of the maximum with `PiConstant.abTEMFloat32`; 1.1e-7 with `.exact` | 1e-10 and 1e-6 as the other slices' scan pins |
| `.metal` against `.cpu`, the same scan | 1.5e-6 of the maximum | 2e-4 (the other Metal scan tests) |
| `width`, 4 profiles × 2 to 3 heights | 0 to 4.5e-8 of the maximum (float32 storage in abTEM), ratio to the bound 0.46 | rtol 1e-7, atol 1e-12 max |
| `interpolate`, 2 scan profiles × 9 argument sets + 3 planted ensembles × 8 | up to 5.5e-8 of the maximum (the same float32 storage), ratio 0.59 | rtol 1e-7, atol 1e-12 max |
| `interpolate_line` and `_at_position`, 19 lines × (image, pattern) | 3.5e-14 of the maximum (spline arithmetic), ratio below 1.2e-5 | rtol 1e-7, atol 1e-12 max |

The 1D `width` and `interpolate` results are float32 in abTEM (`_calculate_widths` and `_interpolate` allocate float32); the kit returns float64, so the
difference is abTEM's rounding (half a float32 ulp, 6e-8), not the spline. Nothing was loosened; every bound is assert_allclose's float64 default.

## Deviations from abTEM

- `LineScan`: `end` is required and one of `gpts`/`sampling` must be given (abTEM leaves them for `match_probe`); `potential=` is `extent:`; `angle` is in radians (abTEM's docstring says degrees, `arctan2` returns radians); the margin is `addingMargin` returning a new value (a struct), not an in-place mutator.
- `interpolateLine` family: `end` required; `RealField` only (no complex, no lazy); a pattern is assumed centred. Placed in SimulatoDetector (parity.md said SimulatoCore) because it uses `LineScan`.
- `lineProfiles(for:)` for value shapes of at most one axis; the scan axis moves last as in abTEM's `moveaxis`.
- Errors: every refusal is `SimulatoError.invalidScan` (no new case, to keep `SimulatoError.swift` free of conflicts): a profile that never crosses the width height (abTEM: IndexError), both `gpts` and `sampling` to `interpolate` (abTEM: RuntimeError), an order outside 0 to 5, a negative width.
- `STEM.swift` was edited in five signature tokens (`GridScan` to `some ScanPattern`), the `ScanMeasurement` initialiser precondition and subscript, and `average`'s `scanShape`; `ScanOperations.swift` got one guard. The WavesDetector agent's edits to the same file may need a hand merge on those lines.
- Tests live in `SimulatoKitParityTests` (14 tests) and `SimulatoProbeTests` (1 test, `lineScanConventions`); there is no fast test target on the Detector side (Package.swift is not edited here).
- `Tests/Fixtures/README.md`: two rows added for `linescan-pins.json` and `linescan-arrays.bin`.

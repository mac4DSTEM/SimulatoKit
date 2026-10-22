# Notes for the merge: integrated centre of mass (suite task 23, part 2)

Code: `Sources/SimulatoCore/ImageOperations.swift` (new), `Sources/SimulatoDetector/IntegratedCenterOfMass.swift` (new).
Pins: `tools/icom/pins.py` -> `Tests/SimulatoKitParityTests/Fixtures/icom-pins.json`, `icom-arrays.bin`.
Tests: `Tests/SimulatoCoreTests/ImageOperationsTests.swift` (fast, planted truth, refusals),
`Tests/SimulatoKitParityTests/Detector/IcomPinsTests.swift` (abTEM). No shared file was edited.

## New public names
- `RealField.integratingGradient(x:y:) throws -> RealField` (abTEM `Images.integrate_gradient`)
- `ComplexField.integratedGradient(samplingX:samplingY:) throws -> RealField` (real part = x, imaginary part = y, as abTEM's complex `Images`)
- `RealField.tiled(repetitionsX:repetitionsY:) throws -> RealField` (`Images.tile`)
- `RealField.integrateDisc(x:y:radius:mean:border:) throws -> Double`, `RealField.Border { wrap, raise }` (`integrate_disc`)
- `ImageOperationError { gradientComponentsDiffer, invalidRepetitions(Int, Int), invalidDisc(String), discOutsideImage }`
  (its own enum, so `SimulatoError.swift` stays untouched; the grid-length and sampling refusals are `SimulatoError`)
- `ScanMeasurement.gradientImages(samplingX:samplingY:) throws -> (x: RealField, y: RealField)`
- `ScanMeasurement.integratedCenterOfMass(samplingX:samplingY:) throws -> RealField` (`integrated_center_of_mass`)

## CHANGELOG bullets (Added)
- `RealField.integratingGradient(x:y:)` and `ComplexField.integratedGradient(samplingX:samplingY:)`: abTEM's `Images.integrate_gradient`, the Fourier-space integral of a two-component gradient field (NumPy `fftfreq`, zero frequency set to 0, no padding, minimum shifted to 0), CPU float64 on the kit's FFT, any length from 2.
- `ScanMeasurement.integratedCenterOfMass(samplingX:samplingY:)` and `gradientImages(...)`: the integrated centre of mass (iCOM / DPC) of a `CenterOfMassDetector` scan, as abTEM's `integrated_center_of_mass`.
- `RealField.integrateDisc(x:y:radius:mean:border:)` (abTEM `integrate_disc`) and `RealField.tiled(repetitionsX:repetitionsY:)` (`Images.tile`).
- `ImageOperationError`.

## docs/parity.md rows to change
Row 183 (`DiffractionPatterns.center_of_mass`): replace the tail "`integrated_center_of_mass` (DPC) missing" by "`integrated_center_of_mass` is `ScanMeasurement.integratedCenterOfMass` (see below)"; state stays **partial**.

New / changed rows (Python | Swift | state | use):

| `DiffractionPatterns.integrated_center_of_mass()` | `ScanMeasurement.integratedCenterOfMass(samplingX:samplingY:)` of a `CenterOfMassDetector` scan | **pinned**: abTEM's centre of mass through the kit 2e−16 (rtol 1e−7 default); the kit's own centre of mass 5e−8 (its COM is registered at 1e−5); no scaling beyond abTEM's (units of the COM times Å, minimum 0) | study |
| `Images.integrate_gradient` | `RealField.integratingGradient(x:y:)`, `ComplexField.integratedGradient(samplingX:samplingY:)` | **pinned** within assert_allclose float64 defaults (rtol 1e−7, atol 1e−12 of the max); measured 8e−16 of the maximum on even, odd, mixed, 2×3 and random (Nyquist) fields; planted potentials recovered to 1e−9 of their range | study |
| `integrate_disc` (`Images.integrate_disc`) | `RealField.integrateDisc(x:y:radius:mean:border:)` | **pinned** same tolerance; measured 2e−15 relative; `interpolate` not offered | study |
| `Images.tile` | `RealField.tiled(repetitionsX:repetitionsY:)` | **pinned** exact | useful |

Row 190: remove `integrate_disc`, `Images.integrate_gradient`, `Images.tile` from it; `interpolate_line`, `interpolate_line_at_position` stay missing (the LineScan session owns them).

## open-items.md
In the "Study, next" line drop "integrated centre of mass and DPC (`Images.integrate_gradient`)". Add: "iCOM: no calibration beyond abTEM's (the image is in units of the centre of mass times Å, minimum 0); a phase/potential scaling (electrons per mrad to radians) is not in abTEM and not built; `integrate_disc(interpolate:)` and the DPC measurement of a `SegmentedDetector` (differences to a two-value COM-like gradient) are not built."

## Deviations from abTEM
- `integrate_disc`: abTEM places the box's pixels on `linspace(0, m d, m)` with the endpoint (an `Images` axis is endpoint=True), so the spacing is m d/(m−1), not d. Ported as it is (the pins were wrong until this was copied); `interpolate=` and a non-zero image offset are not offered.
- `integrateDisc` rejects a non-finite position and a non-positive radius (abTEM would return garbage or crash).
- Single image only (abTEM maps over ensemble axes; a `ScanMeasurement` here has one scan and one value set).
- abTEM's `integrate_gradient` also relabels the metadata ("iCOM", "arb. unit"); the kit has no metadata on `RealField`.
- iCOM is not scaled to a phase/potential: abTEM does not either, so there is nothing to port; the units are those of the `CenterOfMassDetector` (1/Å or mrad) times Å.

## Measured agreement (float64, abTEM 1.0.10, CPU)
- integrate_gradient: max |Δ|/max|abTEM| 3e−16 (16×12), 5e−16 (15×11), 8e−16 (20×9), 6e−17 (2×3), 4e−16 / 6e−16 (random 14×10 / 13×9, Nyquist included).
- iCOM of a 12×10 Al scan, abTEM's COM as input: 2e−16 (1/Å), 3e−16 (mrad). With the kit's COM (registered at 1e−5): 4.8e−8 (1/Å), 3.8e−8 (mrad).
- integrate_disc: 7 discs (wrapping, corner, r from 0.08 to 2 Å, sum and mean): worst relative 1.6e−15. tile: bit exact.
- Planted truth (both abTEM's pin and the kit): recovers the potential less its minimum to 1e−9 of the range (the unit test uses 64×64, 45×49 and others).

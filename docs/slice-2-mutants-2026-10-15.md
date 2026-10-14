# Slice 2 break-once proof (2026-10-15)

Registration-slice-2.md section 6, gate 4. `tools/slice-2/mutants.py` (slice 1's runner, slice 2's list): one textual change
at a time on a clone, `swift test -c release --filter`; KILLED when the named suites fail, a test-process crash after the
run started included (the runner was taught that on this day; before, such a run was mislabelled BUILD-FAILED).

## First pass, on `3eea476`: 14 killed, 4 survived, 2 mislabelled

- **N02** (the directions ignored on x) survived: the test only left z out. Fixed: a z-only ensemble must leave x and y at
  zero, and its z stream must equal the all-axes stream's (the streams do not depend on the directions).
- **N07** (the annular mask closed at outer) survived: no tested pixel lay exactly on an outer angle. Fixed: an annulus
  whose edge is the 5th pixel's angle, bit for bit; that pixel is counted once.
- **N10** (outer = nil without the floor) survived: no pixel lies between 133 and 133.14 mrad on S1's grid. Fixed: a test of
  `outerAngle` itself (133 exactly).
- **N11** (pixelated "valid" without √2) and **N18** (cutoff gpts without the grid's parity) crashed the test process on a
  size precondition in the comparison helper; the helper now reports +∞, and both are killed by `detectorsEqualAbTEM`.
- **N19** (the segmented pattern not cropped at outer) survived and is **equivalent**: the crop (2⌈outer/Δα⌉ + 1 points,
  grid parity) keeps every pixel with α < outer, so binning the uncropped pattern sums the same pixels. Recorded, not chased.

## Second pass, on `345171f` (N02, N07, N10, N11, N18, N19): 5 killed, N19 equivalent

| Mutant | Rule broken | Killed by |
|---|---|---|
| N01 | σ applied twice (displacement doubled) | averagedPlaneWavePatternEqualsAbTEM, detectorsEqualAbTEM |
| N02 | the directions ignored on x | displacementsHaveTheirSigmaAndRespectTheDirections |
| N03 | every configuration from the same seed | displacementsHaveTheirSigmaAndRespectTheDirections |
| N04 | Box–Muller's angle over half the circle (wrong statistics) | displacementsHaveTheirSigmaAndRespectTheDirections, generatorStatistics |
| N05 | the spare normal equal to the first (correlated pairs) | generatorStatistics |
| N06 | linspace with endpoint inverted | detectorsEqualAbTEM, scanPositionsEqualAbTEM |
| N07 | the annular mask closed at outer | annulusIsHalfOpenOnAPixel |
| N08 | inner as > instead of ≥ | detectorsConserveTheTotal, detectorsEqualAbTEM |
| N09 | the probe's pattern divided by N | detectorsConserveTheTotal, detectorsEqualAbTEM, singleAtomHAADFEqualsAbTEM |
| N10 | outer = nil without the floor | defaultOuterIsTheFlooredCutoffAngle |
| N11 | pixelated 'valid' without √2 | detectorsEqualAbTEM |
| N12 | the segment's φ from atan2(kx, ky) | detectorsEqualAbTEM |
| N13 | the rotation added instead of subtracted | detectorsEqualAbTEM |
| N14 | the radial bin rounded instead of floored | detectorsEqualAbTEM |
| N15 | the average over configurations a sum | detectorsEqualAbTEM |
| N16 | Metal: the batch's forward FFT over axes 0 and 1 | metalScanMatchesTheCPUScan |
| N17 | an outer beyond the maximum angle accepted | refusals |
| N18 | the cutoff gpts without the grid's parity | detectorsEqualAbTEM |
| N19 | the segmented pattern not cropped at outer | equivalent (see above) |
| N20 | the module's scope emptied | refusals |
| N21 | every element scattered as carbon | haadfRisesWithZ |
| N22 | the scan weighted by position | vacuumScanIsFlatAndDPCIsSymmetric |
| N23 | every configuration offset by 1e−9 Å | averagedPlaneWavePatternEqualsAbTEM, zeroSigmaIsSliceOne |

N21–N23 were added after the second pass because three test functions had not yet been reached by any mutant; all three
killed on `345171f`.

Result: 22 of 23 killed, one equivalent (N19). Every slice-2 test function appears in the last column.

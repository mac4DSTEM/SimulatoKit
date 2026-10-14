# Slice 1 break-once proof (2026-10-15)

Registration section 6, gate 4: every new test is broken once by a mutant aimed at its rule. `tools/slice-1/mutants.py`
applies one textual change at a time to a clone of the commit and runs the named suites with `swift test -c release
--filter`; a mutant is KILLED when those tests fail (a build failure does not count).

## First run, on the frozen XCTest commit `11bcdb4`: 38 of 40 killed

Two survived, and each exposed a weak test, not a wrong rule:

- **M24** (the exact π replaced by abTEM's float32 π): the test only asserted that the exact-π potential differs from
  abTEM's at all, which FFT rounding (~1e−14) satisfies. Fixed: the difference must exceed 1e−9 (predicted ~5e−8).
- **M27** (⌊Lz / s⌋ slices instead of ⌈Lz / s⌉): every pinned depth is a whole number of slices. Fixed: a planted test with
  depth 5 Å and s = 2 Å (three slices of 5/3 Å), 40.495 / 2 (21 slices), and an atom placed by the uneven slicing.

## Second run, on `b168c5c` (the tests in Swift Testing, the two fixes in): 40 of 40 killed

| Mutant | Rule broken | Killed by |
|---|---|---|
| M01 | λ without the relativistic term | wavelengthAndSigmaEqualAbTEM, wavelengthAndSigmaMatchKirklandsTable |
| M02 | σ with the rest mass | wavelengthAndSigmaEqualAbTEM, wavelengthAndSigmaMatchKirklandsTable |
| M03 | A and B swapped | bandlimitedTransmissionEqualsAbTEM, potentialEqualsAbTEMWithItsPi, potentialWithTheExactPiIsWithinTheRegisteredDifference, projectedScatteringFactorEqualsAbTEMUnderBothPiChoices |
| M04 | the compensation S dropped | bandlimitedTransmissionEqualsAbTEM, potentialEqualsAbTEMWithItsPi, potentialWithTheExactPiIsWithinTheRegisteredDifference |
| M05 | S with π inside | bandlimitedTransmissionEqualsAbTEM, potentialEqualsAbTEMWithItsPi, potentialWithTheExactPiIsWithinTheRegisteredDifference |
| M06 | the bilinear weights' fx and fy swapped | bandlimitedTransmissionEqualsAbTEM, potentialEqualsAbTEMWithItsPi, potentialWithTheExactPiIsWithinTheRegisteredDifference |
| M07 | the propagator's sign flipped | exactPiStaysWithinTheRegisteredDifference, exitWavesPatternsAndScalarsEqualAbTEM, probeThroughVacuumEqualsTheDefocusedProbe |
| M08 | the propagator over the whole thickness | exactPiStaysWithinTheRegisteredDifference, exitWavesPatternsAndScalarsEqualAbTEM, probeThroughVacuumEqualsTheDefocusedProbe, thicknessSeriesIsConsistentMonotoneAndConservesIntensity |
| M09 | the antialias cutoff 1/2 | bandlimitedTransmissionEqualsAbTEM, exactPiStaysWithinTheRegisteredDifference, exitWavesPatternsAndScalarsEqualAbTEM |
| M10 | the transmission not bandlimited | bandlimitedTransmissionEqualsAbTEM, exactPiStaysWithinTheRegisteredDifference, exitWavesPatternsAndScalarsEqualAbTEM |
| M11 | propagate before transmit (part 1: transmit removed here, re-added after the inverse by M11b's second edit) | exactPiStaysWithinTheRegisteredDifference, exitWavesPatternsAndScalarsEqualAbTEM |
| M12 | the top-face rule dropped | topFaceAtomGoesToTheFirstSlice |
| M13 | the slice-edge nudge dropped | atomOnADriftedSliceBoundaryGoesToTheNextSlice, potentialEqualsAbTEMWithItsPi, slicingEqualsAbTEM |
| M14 | the probe aperture hard | exactPiStaysWithinTheRegisteredDifference, exitWavesPatternsAndScalarsEqualAbTEM, incidentProbeEqualsAbTEM |
| M15 | the probe normalised in real space | exactPiStaysWithinTheRegisteredDifference, exitWavesPatternsAndScalarsEqualAbTEM, incidentProbeEqualsAbTEM, probeNormalisationSymmetryAndShift |
| M16 | χ's sign flipped | exactPiStaysWithinTheRegisteredDifference, exitWavesPatternsAndScalarsEqualAbTEM, incidentProbeEqualsAbTEM |
| M17 | the shift kernel's sign flipped | exactPiStaysWithinTheRegisteredDifference, exitWavesPatternsAndScalarsEqualAbTEM, incidentProbeEqualsAbTEM, probeNormalisationSymmetryAndShift |
| M18 | the crop's first parity step dropped | antialiasApertureAndCutoffSizes, exitWavesPatternsAndScalarsEqualAbTEM |
| M19 | the plane wave's 1/N dropped in the pattern | exactPiStaysWithinTheRegisteredDifference, exitWavesPatternsAndScalarsEqualAbTEM, probeNormalisationSymmetryAndShift, thicknessSeriesIsConsistentMonotoneAndConservesIntensity |
| M20 | Metal: the propagator not applied | metalMatchesTheCPUPathOn512Squared100Slices, metalMatchesTheCPUPathOnThePinnedCases |
| M21 | Metal: the inverse FFT scaled by 1/sqrt(N) | metalMatchesTheCPUPathOn512Squared100Slices, metalMatchesTheCPUPathOnThePinnedCases |
| M22 | abTEM's float32 π replaced by the exact π | bandlimitedTransmissionEqualsAbTEM, exitWavesPatternsAndScalarsEqualAbTEM, potentialEqualsAbTEMWithItsPi |
| M23 | Swift's Float.pi (rounded toward zero) for the float32 π | kappaAndFloat32PiEqualThePins, projectedScatteringFactorEqualsAbTEMUnderBothPiChoices |
| M24 | the exact π replaced by abTEM's | potentialWithTheExactPiIsWithinTheRegisteredDifference, projectedScatteringFactorEqualsAbTEMUnderBothPiChoices |
| M25 | a zero energy accepted | refusals |
| M26 | the bilinear weights no longer sum to one | bandlimitedTransmissionEqualsAbTEM, integratedPotentialEqualsTheSumOfFAtZero, potentialEqualsAbTEMWithItsPi, potentialWithTheExactPiIsWithinTheRegisteredDifference |
| M27 | ⌊Lz / s⌋ slices instead of ⌈Lz / s⌉ | requestedThicknessGivesCeilingManyEqualSlices |
| M28 | slice thicknesses that do not sum to the depth accepted | refusals |
| M29 | the aperture below one inside the cutoff | planeWaveThroughVacuumStaysOne, probeThroughVacuumEqualsTheDefocusedProbe |
| M30 | an anisotropic propagator (breaks 4mm) | aluminiumPatternLivesOnTheSliceLatticeWith4mmSymmetry, probeThroughVacuumEqualsTheDefocusedProbe |
| M31 | atoms shifted along x by depth (breaks the slice lattice) | aluminiumPatternLivesOnTheSliceLatticeWith4mmSymmetry, singleAtomConvergesToTheAnalyticProjectedPotential |
| M32 | a propagator with gain (Σ|ψ|² grows) | planeWaveThroughVacuumStaysOne, probeThroughVacuumEqualsTheDefocusedProbe, thicknessSeriesIsConsistentMonotoneAndConservesIntensity |
| M33 | the probe off its position along y | probeNormalisationSymmetryAndShift |
| M34 | χ's second-order term with α² | chiEqualsTheClosedForm |
| M35 | lengths vDSP refuses accepted | supportedLengthsAreVDSPs |
| M36 | fftfreq's Nyquist bin positive | fFTFreqIsNumPys |
| M37 | the forward FFT inverse along y | fFTMatchesADirectDFTAndRoundTrips |
| M38 | the repeat's axes swapped | repeatFollowsASEOrder |
| M39 | a zero semiangle accepted | refusalsAndConventions |
| M40 | the placeholder's scope emptied | plannedScopeIsStated |

Every test function of the library's suites is killed by at least one mutant (checked by name against the table). The
example's CIF reader (tier 2, `MinimalCIFTests`: four atoms, operations, refusals) was checked instead by three sampled
mutants on 2026-10-15, all killed: the operations ignored, the 90° check removed, fractions inverted.

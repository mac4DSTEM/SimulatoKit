# Performance lane, part 2: break-once proof (2026-10-15)

registration-performance.md section 6, gate 3. `tools/performance/mutants.py` (slice 1's runner) on a clone of the frozen
`88b30da`, one mutant at a time, `swift test -c release -j 4 --filter 'MetalGPUTests|MetalTests|STEMPinsTests|PotentialTests'`
(P15: `PotentialTests`), on the M5 Pro with a Metal device. **15 of 15 killed**; every new test function of the lane fails
under at least one mutant (the benchmark is not a gate and is not mutated).

| Mutant | Rule broken | Killed by |
|---|---|---|
| P01 | σ replaced in the GPU phase | metalMatchesTheCPUPathOn512Squared100Slices, metalMatchesTheCPUPathOnThePinnedCases, metalScanMatchesTheCPUScan, metalScanThroughUnbuiltPotentialsMatchesTheCPUScan, potentialOnMetalMatchesTheCPUPath, potentialOnMetalMatchesTheCPUPathAt1024Squared, renormalizationTightensTheTotalIntensity, runsAreDeterministic, unevenSlicesMatchTheCPUPath |
| P02 | the kernel of the next element | gpuPotentialMatchesTheCPUPotential, potentialOnMetalMatchesTheCPUPath, unevenSlicesMatchTheCPUPath |
| P03 | the delta weights' x and y swapped | gpuPotentialMatchesTheCPUPotential, metalScanThroughUnbuiltPotentialsMatchesTheCPUScan, potentialOnMetalMatchesTheCPUPath, potentialOnMetalMatchesTheCPUPathAt1024Squared, renormalizationTightensTheTotalIntensity, runsAreDeterministic, unevenSlicesMatchTheCPUPath |
| P04 | the first thickness's propagator for every slice | unevenSlicesMatchTheCPUPath |
| P05 | the GPU probe's x ramp sign flipped | metalScanMatchesTheCPUScan, metalScanThroughUnbuiltPotentialsMatchesTheCPUScan, vacuumAndProbesOnTheGPU |
| P06 | the renormalisation squared (s² for s) | renormalizationTightensTheTotalIntensity |
| P07 | the scan intensity as Re² only | metalScanMatchesTheCPUScan, metalScanThroughUnbuiltPotentialsMatchesTheCPUScan |
| P08 | the wave read one slice early | metalMatchesTheCPUPathOn512Squared100Slices, metalMatchesTheCPUPathOnThePinnedCases, metalScanMatchesTheCPUScan, metalScanThroughUnbuiltPotentialsMatchesTheCPUScan, potentialOnMetalMatchesTheCPUPath, potentialOnMetalMatchesTheCPUPathAt1024Squared, renormalizationTightensTheTotalIntensity, runsAreDeterministic, unevenSlicesMatchTheCPUPath |
| P09 | the host merge overwrites instead of summing | runsAreDeterministic, unevenSlicesMatchTheCPUPath |
| P10 | the GPU probe's spectrum normalised 1e−4 off | metalScanMatchesTheCPUScan, metalScanThroughUnbuiltPotentialsMatchesTheCPUScan, vacuumAndProbesOnTheGPU |
| P11 | only the first element summed | gpuPotentialMatchesTheCPUPotential, potentialOnMetalMatchesTheCPUPath, unevenSlicesMatchTheCPUPath |
| P12 | the GPU transmission not bandlimited | metalMatchesTheCPUPathOn512Squared100Slices, metalMatchesTheCPUPathOnThePinnedCases, potentialOnMetalMatchesTheCPUPath, potentialOnMetalMatchesTheCPUPathAt1024Squared, renormalizationTightensTheTotalIntensity, runsAreDeterministic, unevenSlicesMatchTheCPUPath |
| P13 | Potential.build puts slice 0's atoms in every slice | detectorsEqualAbTEM, gpuPotentialMatchesTheCPUPotential, metalScanThroughUnbuiltPotentialsMatchesTheCPUScan, potentialOnMetalMatchesTheCPUPath, potentialOnMetalMatchesTheCPUPathAt1024Squared, renormalizationTightensTheTotalIntensity, unevenSlicesMatchTheCPUPath |
| P14 | the padding carries weight | gpuPotentialMatchesTheCPUPotential, metalScanThroughUnbuiltPotentialsMatchesTheCPUScan, potentialOnMetalMatchesTheCPUPath, potentialOnMetalMatchesTheCPUPathAt1024Squared, runsAreDeterministic, unevenSlicesMatchTheCPUPath |
| P15 | the eager potential drops the parameterisation it was given | buildEqualsTheEagerPotentialBitForBit |

P15 was added after the first fourteen: `buildEqualsTheEagerPotentialBitForBit` compares two paths that share `build()`,
so P13 moved both and the test passed; P15 breaks the eager path's delegation instead.

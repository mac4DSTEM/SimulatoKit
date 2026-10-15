# Generality pass: break-once proof (2026-10-15)

registration-generality.md section 3. `tools/generality/mutants.py` (slice 1's runner) on clones of the frozen commits
(`69568e1`; G14 again after its test was strengthened, G15–G17 on the Gate B fixes), one mutant at a time,
`swift test -c release -j 4 --filter 'GeneralityTests|CTFConventionTests'`. **17 of 17 killed** after two first-pass
survivors were dealt with:

- **G14** (a padded Metal batch read one wave short) survived: the batching test never checked the default run's own
  length. The test now checks the count and the agreement with the CPU path; G14 is killed.
- **G16** (first form: a NaN or negative grid tolerance taken as given) survived as an equivalent mutant: the clamp
  could not change a result. The clamp was removed; G16 now breaks the equal-extents short-circuit and is killed.

Not run, equivalent by design: ignoring `maximumCPUThreads` or `scanChunkBytes` changes no number (both are tested bit
for bit); G12 and G13 instead make those options change a number, and the tests catch it.

| Mutant | Rule broken | Killed by |
|---|---|---|
| G01 | the Metal batch bytes ignored | theMetalBatchPlanFollowsItsOptions |
| G02 | the exact batch shapes ignored | theExecutableCacheFollowsItsLimit, theMetalBatchPlanFollowsItsOptions |
| G03 | the executable cache limit ignored | theExecutableCacheFollowsItsLimit |
| G04 | the renormalisation forced on | renormalizationOffDriftsMore |
| G05 | the grid tolerance ignored | theGridToleranceIsAParameter |
| G06 | the slice ceiling ignored | theSliceCeilingIsAParameter |
| G07 | the Scherzer factor ignored | theClosedFormsTakeTheirConstants |
| G08 | the point-resolution divisor ignored | theClosedFormsTakeTheirConstants |
| G09 | the Nyquist oversampling ignored | theClosedFormsTakeTheirConstants |
| G10 | a registered Metal bound moved | theDefaultsAreTheRegisteredValues |
| G11 | the thread count no longer clamped | theDefaultsAreTheRegisteredValues |
| G12 | one thread returns the last incident unpropagated | oneCPUThreadEqualsTheDefaultBitForBit |
| G13 | a one-position chunk skips the first position | theScanChunkIsAParameter |
| G14 | a padded batch read one wave short | metalBatchingAndShapesChangeRoundingOnly |
| G15 | the cache limit not applied on a hit | edgeValuesAreSafe |
| G16 | equal extents no longer match whatever the tolerance | edgeValuesAreSafe |
| G17 | a cache limit set through the var not clamped | edgeValuesAreSafe |

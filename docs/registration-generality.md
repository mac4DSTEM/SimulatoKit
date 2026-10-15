# SimulatoKit 0.2.0: the generality pass (registration)

Registered 2026-10-15 ~14:05 CEST on `f3bacd6` (0.1.0), before any code of the pass. Tier 1: the defaults must reproduce
every number bit for bit, and a policy that becomes a parameter can move numbers when it is changed.

## 1. What the pass makes

Every choice the kit made for its first users becomes a parameter whose default is that choice, documented where it is
declared. No default changes.

| Policy | Today | Parameter (default = today) |
|---|---|---|
| CPU threads for a batch or scan | `Multislice.maximumCPUThreads` = 8 | `ExecutionOptions.maximumCPUThreads` |
| Metal batch memory | 512 MiB of waves | `ExecutionOptions.metalBatchBytes` |
| CPU scan chunk | 256 MiB of exit waves | `ExecutionOptions.scanChunkBytes` |
| Metal per-step renormalisation | on | `ExecutionOptions.metalRenormalization` |
| Metal executable cache | 12, least recently used | `ExecutionOptions.metalExecutableCacheLimit` |
| Metal batch shapes | powers of two, zero-padded | `ExecutionOptions.metalBatchShapes` (`.powersOfTwo`, `.exact`) |
| The registered Metal bounds (in tests only) | wave 1e−4, intensity, pattern, detectors 2e−4, Σ\|ψ\|² 1.5e−6 per slice, to 100 slices | `MetalAccuracy.registered`, public, read by the tests |
| Grid match tolerance | 1e−12 relative | `Grid.matches(_:relativeTolerance:)` |
| Slice-count ceiling | 10⁶ | `maximumSliceCount:` on `Slicing.thicknesses`, `Slicing.validated` and the `Potential` initialisers |
| Scherzer defocus | sign(Cs) √(1.5 \|Cs\| λ) (abTEM) | `scherzerDefocus(Cs:energy:factor:)`, factor 1.5 (4/3 is the other convention in use) |
| Point resolution | (λ³ \|Cs\| / 6)^¼ (abTEM) | `pointResolution(Cs:energy:divisor:)`, divisor 6 |
| Nyquist sampling | 1 / (4 α_c / λ) (abTEM) | `nyquistSampling(semiangleCutoff:energy:oversampling:)`, 4 |

`ExecutionOptions` is passed as `options:` (default `.standard`) to `Multislice.run`, `Multislice.runBatch` and
`STEM.scan`, for both potential forms. Not parameters, and why: abTEM's own conventions that define agreement with it
(the antialias cutoff 2/3 and taper, the slice-edge nudge 1e−12 and the top-face rule 1e−10 Å, NumPy's `isclose` in the
slicing and the detector's outer angle), the π choice (already `LobatoParameterization.piConstant`), the frozen-phonon
seed rule (displacements can be given), and the Metal detector floor (amendment 3, item 6, of the performance
registration: the owner has not decided).

Also in the pass: no application or suite name in the sources, test comments, README, DocC, `Package.swift` or the
example (a `git grep -i` sweep as the proof); a DocC catalogue per module with one article each; a second consumer that
builds the README's quick start verbatim; the macOS 15 floor measured again (the declaration stays macOS 27).

## 2. The proof that the defaults move nothing

1. Every existing fixture and test passes unchanged (`tools/test.sh`), the abTEM pins bit-close as before.
2. **Bit digests before and after**: a scratch test (`DefaultsDigest`, in the session's scratchpad, not committed)
   hashes the float64 bits of V, plane-wave runs at every third plane, a 37-probe batch and a 7 × 5 scan with annular,
   pixelated and segmented detectors, through both potential forms on `.cpu` and `.metal`, and the three CTF closed
   forms. Recorded on `f3bacd6`; the same test on the pass's code must give the same 14 digests.

## 3. Tests of the parameter paths (tier 2)

`maximumCPUThreads` 1 equals the default bit for bit; `scanChunkBytes` of one position equals the default bit for bit;
`metalBatchBytes` small enough for batches of one or two, and `.exact` shapes, within 1e−6 of the default (the GPU's FFT
may pick another kernel for another batch); `metalExecutableCacheLimit` 2 holds at most 2; `metalRenormalization` off
drifts more than on; `MetalAccuracy.registered` equals the registered numbers; `Grid.matches` refuses a 1e−10 offset by
default and accepts it at 1e−9; `maximumSliceCount` 10 refuses 11 slices; the three closed forms at their alternative
constants equal their textbook formulas. Each test broken once (`tools/generality/mutants.py`, planned: the threads
option ignored, the batch bytes ignored, the cache limit ignored, renormalisation forced on, a tolerance parameter
ignored, the ceiling ignored, each closed form's parameter ignored). One Gate B refuter on the frozen commit: a default
that moved, an option that changes a number it should not, a refusal that differs.

## 4. Amendments

### 2026-10-15 ~14:45 CEST, at the first full run

**The defaults moved nothing.** The 14 bit digests of section 2 (`DefaultsDigest`, scratch) are identical on `f3bacd6` and
on the pass's code; every fixture passes unchanged (`tools/test.sh`, 84 test functions). One addition to section 1: the
Metal batching of `runBatch` and of scans now comes from one function, `batchPlan` (a package hook,
`Multislice.metalBatchPlan`, lets a test read it), so `metalBatchBytes` and `metalBatchShapes`, which change rounding only,
have an observable effect to test. Two planned mutants are equivalent by design and are not run: ignoring
`maximumCPUThreads` or `scanChunkBytes` changes no number (both are tested bit for bit against the default). The macOS
floor, measured again on this code: it builds at macOS 15; at macOS 14 only `Synchronization.Mutex` fails (four uses in
`MetalEngine.swift`); the declaration stays macOS 27.

### 2026-10-15 ~15:15 CEST, Gate B: the refuter on the frozen `69568e1`, and the fixes it required

One refuter (a Sonnet subagent; `GRProbes.swift` and `GRCrash.swift` in its scratch clone, not committed): the 14
digests identical to 0.1.0's; every changed public declaration only gained trailing defaulted arguments, and every old
call site compiles and gives bit-identical results; `.cpu` bit-equal across 13 thread and 15 scan-chunk values; the
renormalisation option reaching every entry point; no cache cross-contamination between options; four threads with
different options, no mismatch; refusals identical across backends; the CTF constants at 0, NaN, ±∞ and negative all
refused; worst Metal numbers 8.4e−7 (batch options) and 4.0e−6 (scans), all inside `MetalAccuracy.registered`.

Fixed (tier 2, a test each in `edgeValuesAreSafe`, broken once by G15–G17):
1. **A cache limit of −1 or `Int.min` set through the public `var` trapped** in the cache's trim: every count of
   `ExecutionOptions` is now clamped to at least 1 on every set, and the engine clamps again.
2. **The cache limit held only on inserts**: seven entries stayed under limit 1 after a run of hits; it now holds on hits.
3. **`Slicing.thicknesses(..., maximumCount: Int.max)` trapped** in `Int(n)` at n = 2^63: the count is capped at 2^53.
4. **`Grid.matches` with a negative or NaN tolerance refused a grid compared with itself**: equal extents now match whatever
   the tolerance (a negative or NaN one then accepts nothing else, as 0 would).
5. **The `ProjectedPotential` initialisers had no `maximumSliceCount:`**; they have it now, as `Potential`'s.

Recorded, not changed: a slice count over the ceiling is refused as `invalidSliceThickness` (whose text names a
thickness), and a bad `oversampling` as `invalidAperture` (a new error case would break consumers' exhaustive switches;
0.3.0 may add one); unapplied references to the changed functions (`grids.allSatisfy(grid.matches)`) no longer
type-check, a source break named in `docs/consumers-0.2.0.md`; a zero-slice `runBatch` on `.metal` returns the incident
rounded to float32. The tests that count cached executables now hold the Metal run lock around the run and the count
(`Multislice.withMetalExclusive`), so a run of another suite in parallel cannot change the count between the two.

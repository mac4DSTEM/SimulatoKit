# Consumers of SimulatoKit 0.2.0

The generality pass (`docs/registration-generality.md`): every choice the kit made for its first users is now a
parameter, and every default is the old value. **No number moves** on the defaults: 14 bit digests of both backends and
both potential forms are identical on 0.1.0 and 0.2.0, and every fixture is unchanged.

## What a consumer changes

Nothing, unless it wants another choice. All new parameters have defaults, so 0.1.0 call sites compile unchanged:

- `ExecutionOptions` (new, `options:` on `Multislice.run`, `Multislice.runBatch` and `STEM.scan`): `maximumCPUThreads`
  (8), `metalBatchBytes` (512 MiB), `scanChunkBytes` (256 MiB), `metalRenormalization` (on), `metalExecutableCacheLimit`
  (12), `metalBatchShapes` (`.powersOfTwo`, or `.exact`). The CPU options never change a number; the Metal ones change
  rounding only, within `MetalAccuracy.registered`, except the renormalisation, which changes the agreement
  (`docs/registration-performance.md`, amendment 1).
- `MetalAccuracy.registered` (new): the Metal-against-CPU bounds the kit tests (wave 1e−4, intensity, pattern and
  detectors 2e−4 of the maximum, Σ|ψ|² 1.5e−6 per slice, to 100 slices and 1024²). A consumer comparing its own Metal
  results with a CPU fixture can use these instead of copying numbers.
- `Grid.matches(_:relativeTolerance:)` (1e−12), `maximumSliceCount:` on `Slicing.thicknesses`, `Slicing.validated` and
  both `Potential` initialisers (10⁶).
- `ContrastTransferFunction.scherzerDefocus(Cs:energy:factor:)` (1.5, abTEM's; 4/3 the other convention),
  `pointResolution(Cs:energy:divisor:)` (6), `nyquistSampling(semiangleCutoff:energy:oversampling:)` (4); each refuses a
  non-positive or non-finite constant.

A consumer that only calls these functions changes nothing. One source break: an unapplied reference to a changed
function no longer has its old type (`grids.allSatisfy(grid.matches)` fails to type-check; write
`grids.allSatisfy { grid.matches($0) }`), and a protocol requirement that mirrors an old signature must add the parameter.
`ExecutionOptions` raises any count below 1 to 1, in the initialiser and on every later set.

## Also new

- A DocC catalogue per module, each with one article, including "Metal on Apple GPUs: pitfalls and design" in
  SimulatoMultislice (autorelease pools in long loops, bounded executable caches, one encoded step against a graph
  loop, MPSGraph's data-dependent FFT loss and the renormalisation that corrects it).
- `tools/consumer-check/readme.sh`: a second consumer that builds and runs the README's code blocks verbatim.
- The README states the choices and their defaults; the README's quick start had an invalid `import a, b` line (found by
  the second consumer) and is fixed.

## Not changed

The Metal detector bound for very weak signals (performance amendment 3, item 6) waits for a decision.

## The platform floor (owner's decision, 2026-10-15 13:58)

`Package.swift` now declares **macOS 15 and iOS 18** (was macOS 27), so a consumer on macOS 15 or later, or on
iPadOS/iOS 18 or later, can depend on the kit. Every module builds for both, the Metal path included (MPSGraph's FFT
exists from macOS 14 / iOS 17); no availability branch was added. The proof at the new floor: `tools/test.sh` green (85
test functions), the consumer check 326 of 326 and the README consumer 7 of 7, the 14 bit digests of both backends
identical to 0.1.0's, the Metal timing unchanged (512², 100 slices: 30 to 43 ms, against 36 to 48 ms before), and
`xcodebuild` builds of all four libraries for an iOS device and the iOS simulator without a warning. The numbers are
validated on macOS only; on iOS the same code runs but no test runs there yet.

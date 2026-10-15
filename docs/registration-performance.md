# SimulatoKit performance lane, part 2: the Metal path on the GPU end to end (registration)

Registered 2026-10-15 11:15 CEST (`1b34eb6`), before any code of this lane, on `aa0cae3`. Tier 1: the Metal path's numbers move
(the potential and the transmissions are computed in float32 on the GPU instead of in float64 on the CPU). The CPU path
does not change and stays the reference; every bound below is slice 1's or slice 2's registered bound unless it is named
as new. Part 1 (`235d661`, tier 2) was the CPU batch on up to 8 threads.

## 1. What this lane makes

1. **`Potential`**, abTEM's lazy `Potential` (`potentials/iam.py`): the grid, the slice thicknesses, the wrapped atoms with
   their slice and the parameterisation, without the arrays. `build()` gives the `ProjectedPotential` (abTEM's
   `PotentialArray`) on the CPU; `ProjectedPotential(structure:...)` becomes `Potential(...).build()`, the same
   operations in the same order (bit for bit; tier 3 on its own). `FrozenPhonons.lazyPotentials(...)` gives one per
   configuration.
2. **One compiled MPSGraph executable per step**, encoded once per slice into one command buffer and committed once, in
   place of one synchronous `graph.run` per slice with host round trips. The executable's shapes do not depend on the
   number of slices, so one compilation (about 10 ms, measured below) serves every potential on that grid and batch.
   Measured before this registration (scratch prototype, M5 Pro, macOS 27, 100 slices, the potential built in the step):
   a `for` loop over slices inside one graph ran 16 ms at 512², but compiled in 0.65 to 1.9 s per shape, and its shapes
   carry the slice count; the encoded step ran 20 ms at 512² and 87 ms at 1024² (103 ms for the loop) and compiled in
   10 ms. The encoded step is the design; the loop is rejected for its compile time.
3. **The projected potential on the GPU** (for a `Potential` on `.metal`), per slice, in the step: the bilinear deltas of
   section 3.5 of slice 1 computed on the host in float64 (the same weights and wrapped indices as the CPU path), merged
   per element and pixel in float64 (contributions to one pixel summed on the host, so the GPU scatter sees each index
   once and is deterministic), cast to float32, padded to the largest slice with weight 0 at index 0, scattered into
   [E, nx, ny] (E the structure's elements in increasing Z), FFT over the last two axes, times the element's kernel
   f(k²)/S(k) (computed on the host in float64 by the CPU path's code, cast to float32 once per grid and element), summed
   over elements, inverse FFT, real part. Then t = IFFT(FFT(exp(iσV)) · A) and ψ ← IFFT(FFT(ψ · t) · P) as before.
4. **A `ProjectedPotential` on `.metal`** (built on the CPU in float64): V of every slice uploaded once as float32
   [nx, ny] per slice; the step computes t from V on the GPU. No float64 or float32 copy of the transmission functions is
   held on the host or on the device (slice-1 amendment 3, item 8: 960² × 100 slices peaked at 3.8 GB).
5. **The propagator** for each distinct thickness computed on the host in float64 by the CPU path's code and cast once;
   the aperture A likewise.
6. **STEM on `.metal`**, for both potential forms: the probe built on the GPU as IFFT(Ψ₀ · r_x ⊗ r_y), with
   Ψ₀ = aperture · exp(−iχ) / ‖·‖ computed once on the host in float64 (the CPU path's `Probe.build` without the shift; its
   norm does not depend on the position) and the shift's two one-dimensional ramps r_x[ix] = exp(−2πi p_x[ix] u_x),
   r_y[iy] = exp(−2πi p_y[iy] u_y) computed per position on the host in float64 exactly as `Probe.build` computes them
   (so no float32 phase of thousands of radians), cast to float32; the diffraction intensity |FFT(ψ)|² (÷ N² for a
   plane wave's normalisation) computed on the GPU and read back widened to float64; the detectors read it on the host
   as now. The positions run in batches whose waves fit in 512 MB (slice 2, section 3.6).
7. **The FFT gain compensation** (an attempt; amendment 2 of slice 1, `tools/metal-fft-gain/`): per grid shape and FFT
   layout, at first use, the round-trip gain G(k) = Σ_m Re(conj(Y_m) Z_m) / Σ_m |Y_m|² with Y_m = FFT(ψ_m),
   Z_m = FFT(IFFT(FFT(ψ_m))) on the GPU, for M = 16 fields ψ_m of seeded uniform random phase (the kit's `GaussianGenerator`,
   seed fixed), the sums in float64 on the host; C(k) = 1 / G(k) multiplies P, A and the kernels (one multiplication per
   round trip, folded into arrays the step already multiplies by). The ratio cancels the forward transform's own gain,
   so C corrects the round trip's. It is kept only if it tightens the measured per-slice drift of Σ|ψ|²; either way the
   registered bound (1.5e−6 × slices, accepted by the owner) does not change, and the new measured value is reported.
8. **Memory measured**: the device's `currentAllocatedSize` before and after, the peak over the run, for 512² and 960²
   × 100 slices and a 64-position batch.

## 2. The reference

The CPU path of slices 1 and 2 (float64, Accelerate), itself pinned against abTEM 1.0.10. Agreement with it is the gate.

## 3. Bounds (Metal against the CPU path, relative to the array's max |·|)

| Quantity | Bound | Predicted | Source |
|---|---|---|---|
| exit wave, intensity, pattern (cutoff), per exit plane | 1e−4, 2e−4, 2e−4 | ≤ 6e−5, 1.2e−4, 1.5e−4 | slice 1, up to 100 slices, grids to 512² |
| Σ\|ψ\|² | 1.5e−6 × slices | ≤ 1e−6 per slice; ≤ 2e−7 with C | slice 1, amendment 2 |
| annular, segmented, pixelated over a scan | 2e−4 × max | ≤ 3e−5 | slice 2 |
| the same at 1024² (new) | the slice-1 and slice-2 bounds | wave ≤ 8e−5 | new: slice 1 registered grids to 512² |
| V per slice built on the GPU against the CPU's V (new) | 1e−5 × max V | ≤ 2e−6 | new |
| the probe built on the GPU against `Probe.build` (new) | 1e−5 × max \|ψ\| | ≤ 1e−6 | new |
| a `Potential` on `.metal` against its `build()` on `.cpu` | the slice-1 bounds | as the first row | the same check as slice 1 |

A Metal result beyond a bound is a failure of the Metal path, never a tolerance to widen.

## 4. Truth (independent of the CPU path)

1. **The k = 0 identity on the GPU.** Σ_pixels V Δx Δy of each GPU slice = Σ_atoms f(0) within 1e−5 relative (slice 1,
   truth 2, at float32).
2. **Vacuum.** A plane wave through a GPU-built potential with no atoms stays 1 at every pixel within 1e−6; a GPU-built
   probe's Σ_k |Ψ(k)|² = 1 within 1e−6.
3. **Determinism.** Two runs of the same input give bit-equal exit waves and scan values (the merged scatter).
4. **The gain map**, if kept: a round trip with C applied preserves Σ|ψ|² of a random-phase field not used to measure it
   within 2e−7 per round trip, against about 5e−7 without C.

## 5. Speed (reported, not gated)

End to end, `.metal` against `.cpu` (8 threads for a batch), Al [001] fcc at 300 keV, 100 slices of a/2, 256², 512²
and 1024² (cells repeated to keep about 0.05 Å per pixel), the incident plane wave run and a 64 × 64 HAADF scan, from the
structure to the numbers (potential included). The CPU scans at 512² and 1024² are timed on a subset of positions and
extrapolated per position (the CPU's cost is linear in positions; the subset is stated). A benchmark test, enabled only
by `SIMULATO_BENCHMARK=1`, prints the table; the numbers go in `docs/status.md`. Predicted: the plane-wave run 20× to 50×
faster at 512², the scan 5× to 10×.

## 6. Gates and predicted outcomes

1. The bounds of section 3 on every pinned case and beam of slice 1 (both potential forms), the 512² × 100-slice crystal,
   a 1024² × 100-slice crystal, the slice-2 scans S1–S4 (both potential forms), and the truth of section 4. The Metal
   tests run only with a device, the reason named.
2. `Potential.build()` equals `ProjectedPotential(structure:...)` bit for bit on every pinned case (and the existing pins
   pass unchanged).
3. Every new test broken once by a mutant (`tools/performance/mutants.py`, one at a time on a clone). Planned: σ dropped
   from the GPU phase; the kernel of the next element; x and y of the delta weights swapped; the propagator of the first
   thickness for every slice; the shift ramp's sign flipped; the gain map applied twice; the intensity as Re² only; the
   exit plane one slice late; the host merge dropped (a duplicate index overwritten instead of summed); the plane wave's
   ÷ N² dropped.
4. One Gate B refuter on the frozen commit: other structures and elements, uneven slicings, many elements per slice,
   atoms stacked in a column (the merge), 1024², non-square and odd-factor grids (120 × 160, 240 × 480), exit planes,
   batches across the 512 MB split, frozen-phonon ensembles, scans at the corner.
5. `tools/test.sh` green with the stamp; no new warning.

## 7. What is not claimed

No speed bound is gated. No change to the CPU path or to any number of a CPU result. The Metal bounds hold to 100 slices
(slice-1 amendment 3, item 7); past that a run needs its own bound. The CTF stays on the CPU (slice 3, section 6). A
compute-shader FFT with float64-rounded twiddles is not part of this lane.

## 8. Amendments

### 2026-10-15 ~11:30 CEST, at the first test run (written by 11:51, `88b30da`): the gain map rejected, the per-step renormalisation in its place

**The gain map (1.7) is rejected; the cause is measured.** Prototype (scratch, M5 Pro, macOS 27): the map G(k) measured
from 16 random-phase fields, then 100 round trips of fields not used to measure it, Σ|ψ|² change per round trip without
and with C = 1/G:

| Field | 256² without / with C | 512² | 1024² |
|---|---|---|---|
| random phase | −9.4e−7 / −5.8e−9 | −8.2e−7 / −2.7e−9 | −1.3e−6 / −2.9e−9 |
| weak phase (0.3 rad) | −3.2e−8 / −2.3e−7 | −1.7e−8 / −1.1e−7 | −1.6e−8 / +1.0e−9 |
| a delta | 0 / +9.5e−7 | 0 / +9.5e−7 | 0 / +1.4e−6 |

The loss is not a fixed gain per frequency: it depends on the data (a delta passes every butterfly exactly; a
weak-phase field, mostly the zero frequency, loses little), so a map measured on one kind of field over-corrects
another, as amendment 2 of slice 1 predicted for a scalar. Not kept.

**In its place, the per-step renormalisation** (in `MetalGraphs.renormalized`, on by default): after each FFT pair the
device measures E_x = Σ|x|² before the forward transform, E_X = Σ|X|² after it, E_Y = Σ|X · P|² (or |X · A|² for the
transmission) and E_y = Σ|y|² after the inverse, and scales y by √(E_Y E_x / (E_X E_y)); so the round trip's Σ|·|² follows
Parseval, Σ|y|² = Σ|P · FFT(x)|² / N, with the forward transform's measured loss removed, for each wave of a batch and
for the transmission. 1 where a product is 0. No array of the CPU path changes. Measured against the CPU path (512²,
100 slices, Al, plane wave and probe, both potential forms), per slice, without / with: plane 2.8e−7 / 1.8e−8 (atoms),
2.9e−7 / 1.4e−7 (values); probe 9.2e−7 / 2.2e−8, 9.3e−7 / 1.05e−7; and the errors fall with it (probe wave 5.4e−5 /
1.6e−5, pattern 1.3e−4 / 6.5e−5). At 1024² × 100 slices: 8.4e−8 per slice, wave 1.3e−5. **The new measured bound is
1.4e−7 per slice** (worst over these runs), against 9.3e−7 before; slice 1's accepted 1.5e−6 × slices stays the
registered bound, and this lane's test guards the renormalised drift at 3e−7 per slice (`renormalizationTightensTheTotalIntensity`).

**Mutants (gate 3), as run:** P01 σ replaced in the GPU phase; P02 the next element's kernel; P03 the delta weights'
x and y swapped; P04 the first thickness's propagator for every slice (a test with uneven slices added for it,
`unevenSlicesMatchTheCPUPath`); P05 the GPU probe's x ramp sign; P06 the renormalisation squared (for "the gain map
applied twice"); P07 the scan intensity as Re² only; P08 the wave read one slice early (for "one slice late"); P09 the
host merge overwriting (a column of eight gold atoms in one slice added to `runsAreDeterministic`); P10 the GPU probe's
spectrum normalised 1e−4 off (for "the plane wave's ÷ N² dropped": the scans run probes only); P11 only the first
element summed; P12 the GPU transmission not bandlimited; P13 `build()` putting slice 0's atoms everywhere; P14 the
padding carrying weight. The slice-1 and slice-2 mutants that named the old per-slice graph (M20, M21, N16) now name the
step in `MetalEngine.swift`; M15, M17 and M33 name `build(at:)` by a longer text (the probe's new helpers repeat a line);
N03 and N17 had gone stale at slice 2's refuter fixes and name the current lines.

### 2026-10-15 ~11:50 CEST (`88b30da`), the speed and memory of section 5, and a memory fix

Measured by `BenchmarkTests` (`SIMULATO_BENCHMARK=1`; log `bench-2026-10-15.log` in the session's scratchpad), M5 Pro,
macOS 27, Al [001] at 300 keV, 100 slices of a/2, about 0.047 Å per pixel, structure to numbers (the potential
included). Other sessions ran heavy tests on the machine at the same time, so the CPU times are, if anything, slow.

| Grid | Plane wave: Metal / CPU | Speed-up | 64 × 64 HAADF scan: Metal / CPU | Speed-up | Metal peak memory |
|---|---|---|---|---|---|
| 256² | 0.089 s / 0.585 s | 6.6× | 16.6 s / 516 s (from 256 positions) | 31× | 0.98 GB |
| 512² | 0.133 s / 4.45 s | 33× | 90.2 s / 2521 s (from 64 positions) | 28× | 1.3 GB |
| 1024² | 0.281 s / 53.5 s | 190× | 408 s / 15157 s (from 64 positions) | 37× | 1.7 GB |

The CPU scan runs 8 threads; its time is the subset's, scaled by 4096 / subset (linear in positions). The predictions
(20× to 50× for the plane wave at 512², 5× to 10× for the scan) were low for the scan. The renormalisation costs 20 % to
33 % of the GPU time (512², 100 slices, one wave: 36 ms without, 43 to 48 ms with). At 1024² the scan is GPU-bound at
about 100 ms per position (21 positions per 512 MB batch); the host's share (ramps, widening, detectors) is a few ms.

**Memory, fixed before this table.** The first run peaked at 5.6, 19.8 and 73 GB of device memory: each batch's buffers
(autoreleased Objective-C objects) lived until the scan returned. Every slice step and every batch now runs in its own
autorelease pool; the peaks above are after the fix (the earlier per-slice design held 960² × 100 slices at 3.8 GB).

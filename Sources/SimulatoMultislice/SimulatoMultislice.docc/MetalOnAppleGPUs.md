# Metal on Apple GPUs: pitfalls and design

Five problems met while running a long multislice computation on MPSGraph, how each was fixed, and how SimulatoMultislice's `.metal` path is organised.

## Overview

The `.metal` backend of ``Multislice`` runs the whole propagation on the GPU in float32 with MPSGraph. The first working version was correct and slow, and then used far more memory than the data justified. The five findings below are not specific to electron microscopy. They apply to any Swift program that compiles MPSGraph executables and runs them in a loop, or that relies on float32 FFTs on Apple GPUs.

All measurements were made on an M5 Pro, macOS 27, 2026-10-15.

### 1. Autoreleased objects in a long loop are not freed until the pool drains

MPSGraph and Metal return many autoreleased objects: `MPSGraphTensorData`, `MTLBuffer`, `MPSCommandBuffer`. In a Swift loop that runs for a long time without returning to the run loop, the enclosing autorelease pool never drains, and every iteration's buffers stay alive.

A 64 x 64 probe scan, propagated in batches, held every batch's buffers until the scan ended. The device's allocated memory peaked at 73 GB for 1024 x 1024 grids and 19.8 GB at 512 x 512. Wrapping each batch, and within it each per-slice encode, in `autoreleasepool` brought the 1024 x 1024 peak to 1.7 GB.

```swift
for start in stride(from: 0, to: items.count, by: batchSize) {
    try autoreleasepool {
        // allocate tensor data, encode, commit, wait, read back
    }
}
```

Check any loop that creates these objects per iteration.

### 2. Compiled executables hold device memory

A compiled `MPSGraphExecutable` owns device memory. One compiled for a batch of 21 waves at 1024 x 1024 held about 240 MB. A cache keyed by batch size, left unbounded, grew to 2.9 GB over twelve different batch sizes.

Two measures keep it small:

- Bound the cache. SimulatoMultislice keeps the 12 most recently used executables and drops the oldest. The limit is ``ExecutionOptions/metalExecutableCacheLimit``.
- Round shapes. A remainder batch (the last, smaller batch of a run) is rounded up to a power of two and zero-padded, so a run compiles a few shapes instead of one per remainder. The padded waves are never read back. This is ``MetalBatchShapes/powersOfTwo``; ``MetalBatchShapes/exact`` compiles the exact size.

### 3. Encode a compiled step per iteration instead of compiling a graph loop

A propagation over many slices can be written as a graph `for` loop (a single large graph) or as one compiled step per slice, encoded repeatedly into one command buffer. The second form compiled in 10 ms; the graph loop compiled in 0.65 to 1.9 s. Run time was the same.

The step is compiled once for a shape, then encoded once per slice into one `MPSCommandBuffer`, with the wave kept in two device buffers that alternate (ping-pong). The results are passed to `encode` explicitly. Without preallocated results every iteration's output stays allocated until the command buffer completes: 6.4 GB for 100 slices at 1024 x 1024.

```swift
var current = waveData          // MPSGraphTensorData
var other = makeTensorData()    // same shape, allocated once
let commandBuffer = MPSCommandBuffer(from: queue)
for slice in 0..<sliceCount {
    try autoreleasepool {
        _ = executable.encode(to: commandBuffer,
                              inputs: inputs(for: slice, wave: current),
                              results: [other],
                              executionDescriptor: nil)
        swap(&current, &other)
    }
}
commandBuffer.commit()
commandBuffer.waitUntilCompleted()
```

Where the wave is needed at intermediate planes, the command buffer is committed and read at those slices and a new one is started.

### 4. MPSGraph's float32 FFT loses energy, and the loss depends on the data

A forward and inverse FFT pair on the device does not preserve Parseval's relation exactly. A delta function loses nothing. A random-phase field loses about 1e-6 of its energy per round trip, and a multislice step does one round trip per slice and per transmission.

A fixed gain map, one correction per frequency measured once on random fields, does not fix it: because the loss depends on the data it over-corrects other fields, and it made fields dominated by the zero frequency worse. It was dropped.

What works is a renormalisation inside each step from energies measured on the device itself. Let E_x be Σ|x|² before the forward transform, E_X after it, E_Y after multiplication by the kernel (the propagator), and E_y after the inverse transform. The output y is scaled by √((E_Y / E_X) · (E_x / E_y)), with the scale set to 1 where E_X or E_y is zero. The result follows Parseval with the forward transform's own loss removed. Measured against the float64 reference, the drift of Σ|ψ|² per slice fell from 9.3e-7 to at most 5e-7 (1.4e-7 on the aluminium case that was tracked most closely). The renormalisation is controlled by ``ExecutionOptions/metalRenormalization``.

Compute the factor as a product of ratios, not as a product of energies. E_Y · E_x overflows float32 for a wave of amplitude 1e8, whereas each ratio stays near 1.

### 5. Float32 has an absolute noise floor on detector sums

Sums of diffraction intensity on the Metal path carry an absolute error of about 3e-10 of the probe's total intensity. A bound stated relative to the image's own maximum fails once the signal is weak enough: below about 1e-5 of the total intensity the floor exceeds the relative bound. A single hydrogen atom seen on a high-angle annular detector is such a case.

A bound for such sums needs an absolute floor next to the relative term, for example 2e-4 times the larger of the image's maximum and 1e-5 of the incident intensity. The bounds in ``MetalAccuracy`` are relative to the array's maximum, so a signal below that floor can exceed them without any defect in the algorithm.

## How the .metal path is built

**The potential is built on the GPU.** For a ``Potential``, the host first computes each atom's bilinear contributions to the grid and merges them per pixel in float64 into a list of deltas per slice and element. The device step then evaluates V = Re IFFT(Σₑ FFT(δₑ) · Kₑ), where Kₑ is the element's projected scattering factor. Merging on the host makes the scatter deterministic, so two runs give bit-equal results. A ``ProjectedPotential`` supplied as values is uploaded as float32 instead.

**Transmission is computed on the device.** t = IFFT(FFT(exp(iσV)) · A) with A the antialias aperture, evaluated once per slice and fed to the step.

**Constants are float64 on the host.** Every array that the step multiplies by (kernels, aperture, propagators, the probe's spectrum and the shift ramps of a scan) is computed on the host in float64 by the same code as the CPU path and cast to float32 once. The device only does transforms, products and reductions.

**The wave stays on the device.** The wave sits in two alternating device buffers. Scans build the probes there from a centred spectrum and per-position ramps, propagate them in batches whose waves fit a memory budget (512 MB by default, ``ExecutionOptions/metalBatchBytes``), and return the diffraction intensity of each position.

**One Metal run at a time.** The executables and the command queue are shared, so a lock serialises Metal runs within a process.

**The CPU path is the reference.** Float64 Accelerate code computes the same quantities, and every Metal result is checked against it. The registered bounds, relative to the array's maximum, are:

| Quantity | Bound |
|---|---|
| exit wave, intensity, diffraction pattern | 1e-4, 2e-4, 2e-4 |
| Σ\|ψ\|² | 1.5e-6 per slice |
| annular, segmented and pixelated detectors over a scan | 2e-4 |
| a potential slice built on the GPU against the CPU's | 1e-5 of the maximum |
| a probe built on the GPU against the CPU's | 1e-5 of the maximum |

The bounds hold up to 100 slices and grids up to 1024 x 1024. Beyond that a run needs its own bound. A result outside a bound counts as a failure of the Metal path and the bound is not widened. These bounds are what ``MetalAccuracy/registered`` holds.

## Using it

```swift
import SimulatoCore
import SimulatoMultislice

let beam = try ElectronBeam(energy: 300e3)
let potential = try Potential(structure: structure, nx: 512, ny: 512, sliceThickness: 2.0)
let incident = Waves.planeWave(beam: beam, grid: potential.grid)

if Multislice.metalIsAvailable {
    let exit = try Multislice.run(incident, through: potential, backend: .metal)[0]
    _ = exit
}

// Tuning for a memory-constrained machine:
var options = ExecutionOptions.standard
options.metalExecutableCacheLimit = 4
options.metalBatchShapes = .powersOfTwo
let exits = try Multislice.runBatch([incident], through: potential, backend: .metal, options: options)
```

Without a Metal device, `.metal` throws `SimulatoError.metalUnavailable` (SimulatoCore); check ``Multislice/metalIsAvailable`` first.

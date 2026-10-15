# Grids, waves and conventions

Array layout, frequency order, the two wave normalisations and the units used throughout SimulatoKit.

## Overview

Every array the library returns or accepts follows a handful of conventions inherited from abTEM and NumPy. They are fixed, so a result can be compared element by element with abTEM's.

### Units

| Quantity | Unit |
|---|---|
| lengths, extents, sampling, structure coordinates, aberration coefficients | Å |
| beam energy | eV |
| angles (aperture semiangle, detector limits, angular spread) | mrad |
| aberration orientation angles (`phi12`, ...) and detector `rotation` | rad |
| reciprocal-space coordinates | 1/Å |
| projected potential | V·Å |
| interaction parameter σ (``ElectronBeam/interactionParameter``) | 1/(V·Å) |

The phase a slice imprints on the wave is σ·V, which is dimensionless.

### Layout: x first, row-major

A grid has `nx` points along x and `ny` along y. Every 2-D array is flat and row-major with x as the first index:

```swift
let grid = try Grid(nx: 128, ny: 64, extentX: 12.0, extentY: 6.0)
var image = [Double](repeating: 0, count: grid.count)
let ix = 5, iy = 7
image[ix * grid.ny + iy] = 1          // abTEM's array[x, y]
```

``ComplexField`` and ``RealField`` offer `subscript(ix, iy)` that does this arithmetic. A ``ComplexField`` stores real and imaginary parts as two separate arrays.

Lengths `nx` and `ny` must be f · 2^m with f in {1, 3, 5, 15} and m of at least 3 (see ``Grid/isSupportedLength(_:)``); both backends accept exactly these. Sampling is extent divided by point count: ``Grid/samplingX``.

### Frequency order

Reciprocal-space arrays use NumPy's `fftfreq` order, unshifted. For n = 8 the integer indices are `0, 1, 2, 3, -4, -3, -2, -1`, multiplied by 1/(n·d). ``Grid/fftfreq(_:_:)`` reproduces NumPy bit for bit, and ``Grid/frequenciesX`` and ``Grid/frequenciesY`` give kx and ky in 1/Å.

The forward FFT is unnormalised with e^(−2πi k·r); the inverse carries 1/(nx·ny), as in NumPy's `fft2` and `ifft2`.

Arrays meant for display are the exception. ``Waves/diffractionPattern(maxAngle:)`` returns a shifted pattern with the zero frequency at index `(m/2, m/2)` of an `m`-sized axis (odd for the default `.cutoff`, `nx / 2` for `.full`). Detectors that return pixel arrays (`PixelatedDetector` in SimulatoDetector) do the same.

### Two normalisations

``Waves/Normalization`` decides how the diffraction pattern is scaled:

- ``Waves/Normalization/values``: the wave is defined by its real-space values. The plane wave from ``Waves/planeWave(beam:grid:)`` has ψ = 1 at every pixel. The pattern is computed from ψ / (nx·ny), so it sums to about 1.
- ``Waves/Normalization/reciprocalSpace``: the wave is defined by its spectrum, Σₖ|Ψ(k)|² = 1. This is a probe. The pattern is |FFT(ψ)|² with no division, and a real-space sum satisfies Σᵣ|ψ|² = 1 / (nx·ny).

The normalisation travels with the wave through multislice and the contrast transfer function, so `intensity()` and `diffractionPattern()` stay consistent.

```swift
let beam = try ElectronBeam(energy: 300e3)       // eV
let plane = Waves.planeWave(beam: beam, grid: grid)
let pattern = try plane.diffractionPattern()      // odd-sized, cropped to the antialias cutoff
print(pattern.nx, pattern.samplingX)              // samplingX in 1/Å
```

### Bandlimit

Each multislice step is bandlimited by ``Antialias`` at two thirds of the Nyquist frequency of the coarser axis, with a cosine taper. The default diffraction pattern is cropped to that cutoff (``Antialias/cutoffPatternSize(on:)``); detector angles beyond it contain no signal.

### Constants

``PhysicalConstants`` holds the CODATA 2014 values that ASE 3.29.0 supplies to abTEM 1.0.10, so ``ElectronBeam/wavelength`` and ``ElectronBeam/interactionParameter`` equal abTEM's to the last bit.

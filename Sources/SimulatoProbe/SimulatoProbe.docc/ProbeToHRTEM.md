# From a probe to an HRTEM image

Build a probe, choose aberrations, and turn an exit wave into a high-resolution TEM image with a contrast transfer function.

## Overview

SimulatoProbe covers both illumination modes. A convergent probe for STEM is built in reciprocal space and shifted to a scan position. HRTEM uses plane-wave illumination, and the objective lens is applied to the exit wave afterwards as a contrast transfer function (CTF). Units follow the conventions of SimulatoCore (Grids, waves and conventions): lengths in Å, angles in mrad, aberration coefficients in Å.

### Aberrations

``Aberrations`` holds the polar coefficients of orders 1 to 5 in Krivanek's notation with abTEM's symbols: `C10`, `C12` with `phi12`, `C21`/`C23`, `C30`/`C32`/`C34`, up to `C56`. Each coefficient with an orientation `phi..` is a length in Å and an angle in rad. The phase is

χ(α, φ) = (2π / λ) Σₙ αⁿ⁺¹ / (n + 1) Σₘ Cₙₘ cos(m (φ − φₙₘ))

and ``Aberrations/chi(alpha:phi:wavelength:)`` evaluates it. Two accessors cover the usual cases: ``Aberrations/defocus`` with abTEM's sign, defocus = −C10, and ``Aberrations/Cs`` which is C30.

```swift
var aberrations = Aberrations()
aberrations.Cs = 1e7           // 1 mm in Å
aberrations.defocus = 300      // Å
aberrations.C12 = 20           // two-fold astigmatism, Å
aberrations.phi12 = 0.3        // rad
```

### A probe

A ``Probe`` combines a beam, a grid, an ``Aperture`` (semiangle in mrad, with abTEM's one-pixel soft edge by default) and aberrations. `build(at:)` returns a `Waves` at a position in Å, normalised so that Σₖ|Ψ(k)|² = 1 (the `reciprocalSpace` normalisation).

```swift
import SimulatoCore
import SimulatoProbe

let beam = try ElectronBeam(energy: 300e3)
let grid = try Grid(nx: 256, ny: 256, extentX: 20, extentY: 20)
let probe = try Probe(beam: beam, grid: grid,
                      aperture: Aperture(semiangleCutoff: 20),
                      aberrations: aberrations)
let waves = try probe.build(at: SIMD2(10, 10))
let intensity = waves.intensity()          // RealField, sampled in Å
```

A ``GridScan`` lists positions in Å. Along each axis it uses NumPy's `linspace(start, end, n, endpoint)`; `endpoint` defaults to `false`. `positions` is flat, index `ix * ny + iy`.

```swift
let scan = try GridScan(start: SIMD2(0, 0), end: SIMD2(20, 20), nx: 32, ny: 32)
let positions = scan.positions             // 1024 entries
```

### The contrast transfer function

``ContrastTransferFunction`` multiplies the exit wave's spectrum by a kernel

H(k) = exp(−iχ) · Eₛ · Eₜ · aperture

with each factor present only when its parameter is set.

| Parameter | Effect |
|---|---|
| `aberrations` | the phase exp(−iχ) |
| `semiangleCutoff` [mrad], `softAperture` | objective aperture; `.infinity` for none |
| `focalSpread` Δ [Å] | temporal envelope exp(−((π/2λ) Δ α²)²) |
| `angularSpread` β [mrad] | spatial envelope from the gradient of χ (quasi-coherent approximation) |

The initialiser throws `SimulatoError` for a non-finite coefficient, a non-positive aperture or negative spreads.

Apply the function to an exit wave with `Waves.applying(_:)`. The intensity of the result is the HRTEM image.

```swift
import SimulatoMultislice

let potential = try Potential(structure: structure, nx: 256, ny: 256, sliceThickness: 2.0)
let exit = try Multislice.run(.planeWave(beam: beam, grid: potential.grid),
                              through: potential)[0]

var lens = Aberrations()
lens.Cs = 1e7
lens.defocus = try ContrastTransferFunction.scherzerDefocus(Cs: lens.Cs, energy: 300e3)
let ctf = try ContrastTransferFunction(aberrations: lens, focalSpread: 30)
let image = try exit.applying(ctf).intensity()     // the HRTEM image
```

### Closed forms

For a given Cs and beam energy:

- ``ContrastTransferFunction/scherzerDefocus(Cs:energy:factor:)``: sign(Cs) √(factor |Cs| λ) in Å, factor 1.5 by default (abTEM's; 4/3 is the other convention in use) (assign it to `Aberrations.defocus`).
- ``ContrastTransferFunction/pointResolution(Cs:energy:divisor:)``: (λ³ |Cs| / divisor)^¼ in Å, divisor 6 by default (abTEM's).
- ``ContrastTransferFunction/nyquistSampling(semiangleCutoff:energy:oversampling:)``: the pixel size in Å that samples a semiangle given in mrad at the Nyquist rate.

The envelope functions are available on their own as ``ContrastTransferFunction/temporalEnvelope(alpha:wavelength:focalSpread:)`` and ``ContrastTransferFunction/spatialEnvelope(alpha:phi:wavelength:aberrations:angularSpread:)``, and ``ContrastTransferFunction/kernel(on:wavelength:)`` returns H on a grid in FFT order.

The CTF is evaluated on the CPU in float64 on every backend.

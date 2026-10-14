# SimulatoKit

The macEM suite's electron-scattering simulation library for Apple platforms: a Swift package of numbers, never UI, that
mac4DSTEM, macTEM and the other kits use as their source of simulated truth, and that other people can build their own
apps on. It is the Mac-native port of [abTEM](https://github.com/abtem/abtem): the same numbers on the same inputs, checked
against abTEM 1.0.10 in float64. Every path runs on the CPU in float64 (Accelerate) as the reference; the multislice loop
also runs on the GPU in float32 (Metal, MPSGraph) behind the same API, within a registered tolerance.

| Module | Scope | State |
|---|---|---|
| `SimulatoCore` | the beam (λ, σ), grids, fields, waves, intensity and diffraction patterns, the antialias aperture, structures | slice 1 |
| `SimulatoMultislice` | the Lobato parameterisation, slicing, the projected potential, transmission, Fresnel propagation, multislice on `.cpu` and `.metal`, frozen phonons, batched propagation | slices 1, 2 |
| `SimulatoProbe` | the soft aperture, polar aberrations C10 ... C56, the probe at a position, scan grids | slices 1, 2 |
| `SimulatoDetector` | annular (BF, ADF, HAADF), pixelated (4D-STEM) and segmented (DPC) detectors, the STEM driver averaged over frozen phonons | slice 2 |

```swift
import SimulatoCore
import SimulatoMultislice

let cell = try AtomicStructure(box: SIMD3(4.05, 4.05, 4.05), atoms: [
    .init(atomicNumber: 13, position: SIMD3(0, 0, 0)), .init(atomicNumber: 13, position: SIMD3(0, 2.025, 2.025)),
    .init(atomicNumber: 13, position: SIMD3(2.025, 0, 2.025)), .init(atomicNumber: 13, position: SIMD3(2.025, 2.025, 0)),
])
let potential = try ProjectedPotential(structure: cell.repeated(4, 4, 10), nx: 256, ny: 256, sliceThickness: 2.025)
let beam = try ElectronBeam(energy: 200e3)
let exit = try Multislice.run(.planeWave(beam: beam, grid: potential.grid), through: potential, backend: .metal)[0]
let pattern = try exit.diffractionPattern()      // |Ψ(k)|², cropped at the antialias cutoff, zero frequency centred
```

`swift run -c release SimulatoExample <file.cif>` simulates an image and a pattern from a CIF and writes PNGs and the
numbers (`Examples/SimulatoExample`).

## Platform

macOS 27 or later (Swift 6.2 tools, Swift 6 language mode, strict concurrency). The Metal path needs a Metal device; without
one `.metal` throws `SimulatoError.metalUnavailable` and the tests skip it with that reason. Grid lengths are f · 2^m with
f ∈ {1, 3, 5, 15} and m ≥ 3 (8 ... 1024 and beyond), the lengths Accelerate's DFT accepts.

## Licence

GPL-3.0-or-later (`LICENSE`, `NOTICE`). The algorithms are ported from abTEM (GPL-3.0-or-later), and the Lobato and Van Dyck
(2014) parameters come from abTEM's data file. An app built on SimulatoKit must be GPL-compatible.

## Build and test

```sh
tools/test.sh                      # the gate: swift test -c release, stamps .build/last-green on a clean tree
tools/last-green.sh                # the last green commit and whether HEAD equals it
```

The oracle (abTEM 1.0.10 in a Python venv) is not in the repository; `docs/registration-slice-1.md` section 2 says how it is
pinned, and `tools/slice-1/pins.py` regenerates the fixtures from it. What agrees with abTEM, and to what tolerance, is in
`docs/registration-slice-1.md`; what of abTEM's API is ported is in `docs/parity.md`; the state is in `docs/status.md`.

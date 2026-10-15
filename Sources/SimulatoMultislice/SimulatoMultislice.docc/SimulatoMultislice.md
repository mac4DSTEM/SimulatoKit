# ``SimulatoMultislice``

Projected potentials, slicing, frozen phonons and the multislice propagation of a wave through them, on a float64 CPU path or a float32 Metal path.

## Overview

A structure is cut into slices along the beam (``Slicing``, ``Potential``). Each slice becomes a projected potential from the Lobato parameterisation of the electron scattering factors (``LobatoParameterization``, ``ProjectedPotential``). ``Multislice`` then alternates transmission and Fresnel propagation, bandlimited at two thirds of Nyquist.

One API serves both backends, selected by ``Backend``. `.cpu` is float64 on Accelerate and is the reference. `.metal` is float32 on MPSGraph and agrees with the reference within registered bounds. Execution parameters of the Metal path and the CPU batch (thread count, batch memory, executable cache) are gathered in ``ExecutionOptions``.

``Potential`` keeps the atoms and slice assignment without the arrays: `.cpu` builds it with ``Potential/build()``, `.metal` builds each slice on the GPU.

See <doc:MetalOnAppleGPUs> for how the Metal path is built and for memory and accuracy pitfalls that apply to any MPSGraph program.

## Topics

### Essentials

- <doc:MetalOnAppleGPUs>
- ``SimulatoMultislice/SimulatoMultislice``

### Propagating waves

- ``Multislice``
- ``Backend``
- ``ExitPlanes``

### Execution options

- ``ExecutionOptions``
- ``MetalBatchShapes``
- ``MetalAccuracy``

### Potentials and slicing

- ``Potential``
- ``ProjectedPotential``
- ``PotentialSlice``
- ``Slicing``
- ``LobatoParameterization``
- ``LobatoParameterization/PiConstant``

### Frozen phonons

- ``FrozenPhonons``
- ``FrozenPhonons/Directions``
- ``GaussianGenerator``

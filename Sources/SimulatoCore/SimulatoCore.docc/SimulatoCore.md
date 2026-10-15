# ``SimulatoCore``

The shared foundation of SimulatoKit: electron beams, lateral grids, structures, complex and real fields, waves and the antialias aperture.

## Overview

SimulatoCore holds the types every other SimulatoKit module builds on. All conventions (array layout, FFT ordering, normalisation, constants) follow abTEM 1.0.10. The module also contains the float64 CPU FFT that serves as the numerical reference for the Metal path.

Read <doc:GridsWavesConventions> before indexing into any array the library returns.

## Topics

### Essentials

- <doc:GridsWavesConventions>
- ``SimulatoCore/SimulatoCore``

### Beam and constants

- ``ElectronBeam``
- ``PhysicalConstants``

### Grids and fields

- ``Grid``
- ``ComplexField``
- ``RealField``

### Waves

- ``Waves``
- ``Waves/Normalization``
- ``Waves/MaxAngle``
- ``Antialias``

### Structures

- ``AtomicStructure``
- ``AtomicStructure/Atom``

### Errors

- ``SimulatoError``

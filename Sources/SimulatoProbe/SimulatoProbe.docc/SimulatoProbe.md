# ``SimulatoProbe``

Convergent probes, the polar aberration series, scan positions and the contrast transfer function for HRTEM images.

## Overview

SimulatoProbe describes the optics in front of and behind the specimen. A ``Probe`` is built from a beam, a grid, an ``Aperture`` and ``Aberrations`` at any position, normalised in reciprocal space. A ``GridScan`` lists probe positions. A ``ContrastTransferFunction`` with its partial-coherence envelopes turns an exit wave into an HRTEM image.

The module depends only on SimulatoCore; the specimen enters as a `Waves` exit wave computed with SimulatoMultislice.

## Topics

### Essentials

- <doc:ProbeToHRTEM>
- ``SimulatoProbe/SimulatoProbe``

### Probes and apertures

- ``Probe``
- ``Aperture``

### Aberrations

- ``Aberrations``

### Scans

- ``GridScan``

### Contrast transfer

- ``ContrastTransferFunction``

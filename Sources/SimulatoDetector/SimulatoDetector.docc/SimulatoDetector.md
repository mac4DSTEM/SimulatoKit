# ``SimulatoDetector``

Detectors for scanned-probe simulations and the STEM driver that applies them over a scan and over frozen-phonon configurations.

## Overview

``STEM`` runs a probe over a scan through one or more potentials, passes each exit wave's diffraction intensity (``DiffractionIntensity``) to every detector and averages the result over the potentials. Annular detectors give bright-field, annular dark-field and high-angle annular dark-field images, the pixelated detector gives 4D-STEM data, and the segmented detector gives DPC-style signals. Any type conforming to ``Detector`` can be used in the same call.

## Topics

### Essentials

- <doc:SimulatingSTEMScan>
- ``SimulatoDetector/SimulatoDetector``

### Running a scan

- ``STEM``
- ``ScanMeasurement``

### Detectors

- ``Detector``
- ``AnnularDetector``
- ``PixelatedDetector``
- ``SegmentedDetector``

### Detector input

- ``DiffractionIntensity``

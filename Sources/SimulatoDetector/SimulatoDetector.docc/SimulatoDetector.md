# ``SimulatoDetector``

Detectors for scanned-probe simulations and the STEM driver that applies them over a scan and over frozen-phonon configurations.

## Overview

``STEM`` runs a probe over a scan through one or more potentials, passes each exit wave's diffraction intensity (``DiffractionIntensity``) to every detector and averages the result over the potentials. Annular detectors give bright-field, annular dark-field and high-angle annular dark-field images, the pixelated detector gives 4D-STEM data, the segmented detector gives DPC-style signals, the centre-of-mass detector gives the data for integrated centre of mass (iCOM), and ``WavesDetector`` keeps the complex exit wave at every position. A scan runs over a grid (`GridScan`) or a line (`LineScan`); a line scan's measurement becomes line profiles. Any type conforming to ``Detector`` can be used in the same call.

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
- ``CenterOfMassDetector``
- ``WavesDetector``
- ``ScanWaves``
- ``WavesDetectorError``

### Line profiles

- ``LineProfiles``
- ``RealSpaceLineProfiles``
- ``ReciprocalSpaceLineProfiles``

### Detector input

- ``DiffractionIntensity``

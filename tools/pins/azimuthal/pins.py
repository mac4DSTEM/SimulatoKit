"""SimulatoKit azimuthal average: DiffractionPatterns.azimuthal_average (abTEM 1.0.10, numpy).

Prints the numbers that Tests/SimulatoDetectorTests/AzimuthalAverageTests.swift holds inline (no fixture file): planted
centred patterns (fftshift=True, as PixelatedDetector and Waves.diffractionPattern return them) of 9 x 8 and 8 x 6
pixels, values (k + 1)^1.3 with k the flat index, float64, at 200 keV. abTEM builds the angular coordinates in float32
and sums in float64; the test holds the profiles at rtol 1e-7 (assert_allclose defaults; abTEM states no tolerance).

Run from the repo root:  .venv/bin/python -I -B tools/pins/azimuthal/pins.py
"""

import numpy as np
from abtem.measurements import DiffractionPatterns


def pattern(nx, ny, sampling):
    arr = (np.arange(nx * ny, dtype=np.float64).reshape(nx, ny) + 1.0) ** 1.3
    return DiffractionPatterns(arr, sampling=sampling, fftshift=True, metadata={"energy": 200e3})


def show(name, p):
    print(name, "count:", p.array.shape[-1], "sampling:", repr(p.sampling))
    print(name, "values:", p.array.ravel().tolist())


a = pattern(9, 8, (0.05, 0.06))
b = pattern(8, 6, (0.04, 0.04))
show("odd default", a.azimuthal_average())
show("odd max 9 rs 0.7 w 1.5", a.azimuthal_average(max_angle=9.0, radial_sampling=0.7, width=1.5))
show("odd gaussian", a.azimuthal_average(weighting_function="gaussian"))
show("odd gaussian max 6 rs 1.3 w 0.8", a.azimuthal_average(max_angle=6.0, radial_sampling=1.3,
                                                           weighting_function="gaussian", width=0.8))
show("even default", b.azimuthal_average())
show("even width 0.5", b.azimuthal_average(width=0.5))
one = DiffractionPatterns(np.array([[3.0]]), sampling=0.05, fftshift=True, metadata={"energy": 200e3})
show("one pixel max 5 w 0.5", one.azimuthal_average(max_angle=5.0, width=0.5))
print("max angles default:", -min(min(a.angular_limits)), -min(min(b.angular_limits)))

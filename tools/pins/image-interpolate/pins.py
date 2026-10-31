"""SimulatoKit image interpolation: Images.interpolate (method "fft") and integrate_disc(interpolate=) (abTEM 1.0.10).

Prints the numbers that Tests/SimulatoCoreTests/ImageInterpolationTests.swift holds inline (no fixture file): a planted
5 x 6 image at sampling (0.3, 0.25) A, values sin(0.7 i) + cos(0.4 j)^2 + 0.05 i j, resampled up, down and mixed by
Fourier-space crop or zero-pad, with both normalizations, and disc integrals of the image interpolated to 0.1 A. abTEM
runs in float64 here (precision "float64"); the test holds the values at rtol 1e-7 (assert_allclose defaults; abTEM
states no tolerance).

Run from the repo root:  .venv/bin/python -I -B tools/pins/image-interpolate/pins.py
"""

import numpy as np
import abtem
from abtem.measurements import Images, integrate_disc

abtem.config.set({"device": "cpu", "precision": "float64"})

i, j = np.meshgrid(np.arange(5), np.arange(6), indexing="ij")
arr = np.sin(0.7 * i) + np.cos(0.4 * j) ** 2 + 0.05 * i * j
img = Images(arr, sampling=(0.3, 0.25))


def show(name, im):
    print(name, "shape", im.shape, "sampling", repr(im.sampling))
    print(name, "values", im.array.ravel().tolist())


show("sampling 0.2", img.interpolate(sampling=0.2))
show("gpts 3 x 4", img.interpolate(gpts=(3, 4)))
show("gpts 9 x 7 intensity", img.interpolate(gpts=(9, 7), normalization="intensity"))
show("sampling 0.17 x 0.4", img.interpolate(sampling=(0.17, 0.4)))
print("disc (0.7, 0.8) r 0.45 sum", repr(float(integrate_disc(img, (0.7, 0.8), 0.45, interpolate=0.1))))
print("disc (0.7, 0.8) r 0.45 mean", repr(float(integrate_disc(img, (0.7, 0.8), 0.45, return_mean=True, interpolate=0.1))))
print("disc (0.1, 1.4) r 0.5 sum wrap", repr(float(integrate_disc(img, (0.1, 1.4), 0.5, interpolate=0.1))))
print("disc (0.7, 0.8) r 0.45 sum, no interpolation", repr(float(integrate_disc(img, (0.7, 0.8), 0.45))))
print("interpolated to 0.1:", img.interpolate(0.1).shape, repr(img.interpolate(0.1).sampling))

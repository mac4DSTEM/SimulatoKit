"""SimulatoKit polar measurements: PolarMeasurements.integrate, integrate_radial and differentials (abTEM 1.0.10, numpy).

Prints the numbers that Tests/SimulatoDetectorTests/PolarMeasurementTests.swift holds inline (no fixture file): a planted
2 x 3 scan of 2 radial x 4 azimuthal bins, values k^1.5 with k the flat index, radial sampling 10 mrad from 5 mrad (the
bins of SegmentedDetector(2, 4, inner=5, outer=25)). Sums are float64; differentials are complex64 in abTEM, so the test
holds them at rtol 1e-6 and the float64 sums at rtol 1e-7 (assert_allclose defaults; abTEM states no tolerance).

Run from the repo root:  .venv/bin/python -I -B tools/pins/polar/pins.py
"""

import numpy as np
from abtem.core.axes import ScanAxis
from abtem.measurements import PolarMeasurements

arr = np.arange(2 * 3 * 2 * 4, dtype=np.float64).reshape(2, 3, 2, 4) ** 1.5
pm = PolarMeasurements(
    arr,
    radial_sampling=10.0,
    azimuthal_sampling=np.pi / 2,
    radial_offset=5.0,
    azimuthal_offset=0.0,
    ensemble_axes_metadata=[
        ScanAxis(label="x", sampling=0.5, units="A"),
        ScanAxis(label="y", sampling=0.5, units="A"),
    ],
)
print("regions [1, 6]:", pm.integrate(detector_regions=[1, 6]).array.ravel().tolist())
print("region 3:", pm.integrate(detector_regions=3).array.ravel().tolist())
print("radial 15-25:", pm.integrate(radial_limits=(15.0, 25.0)).array.ravel().tolist())
print("radial 5-25:", pm.integrate(radial_limits=(5.0, 25.0)).array.ravel().tolist())
print("radial 9-20:", pm.integrate(radial_limits=(9.0, 20.0)).array.ravel().tolist())
print("all:", pm.integrate().array.ravel().tolist())
d = pm.differentials(((0, 4), (2, 6)), ((1, 5), (3, 7)), return_complex=True)
print("differentials real:", d.array.real.ravel().tolist())
print("differentials imag:", d.array.imag.ravel().tolist())
try:
    pm.integrate(radial_limits=(5.0, 35.0))
except RuntimeError as e:
    print("outer beyond the bins:", e)

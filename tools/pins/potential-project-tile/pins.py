"""SimulatoKit potential projection and tiling: PotentialArray.project and PotentialArray.tile (abTEM 1.0.10, numpy).

Prints the numbers that Tests/SimulatoMultisliceTests/PotentialProjectTileTests.swift holds inline (no fixture file):
a planted PotentialArray of 3 slices on a 4 x 5 grid (values k^1.3 with k the flat index, thicknesses 1.0, 1.5, 0.5 A,
extent 2.0 x 2.5 A), its projection and its tilings, and the projection of a built aluminium potential. Sums are float64;
the test holds them at rtol 1e-7 (assert_allclose default; abTEM states no tolerance) and the shapes, extents and
thicknesses exactly.

Run from the repo root:  .venv/bin/python -I -B tools/pins/potential-project-tile/pins.py
"""

import numpy as np
from ase.build import bulk
import abtem
from abtem.potentials.iam import PotentialArray

abtem.config.set({"precision": "float64", "fft": "numpy", "dask.lazy": False, "diagnostics.progress_bar": False})

arr = np.arange(3 * 4 * 5, dtype=np.float64).reshape(3, 4, 5) ** 1.3
pa = PotentialArray(arr, slice_thickness=(1.0, 1.5, 0.5), extent=(2.0, 2.5))


def weighted(a):
    """An order-sensitive checksum: sum of value * (flat index + 1)."""
    flat = np.asarray(a, dtype=np.float64).ravel()
    return float(np.sum(flat * (np.arange(flat.size) + 1)))


p = pa.project()
print("project values:", [repr(float(v)) for v in p.array.ravel()])
print("project sampling:", p.sampling)

for reps in [(2, 3), (2, 1, 2), (1, 1, 3)]:
    t = pa.tile(reps)
    print(f"tile {reps}: shape={t.array.shape} extent={t.extent} sampling={t.sampling}")
    print(f"  slice_thickness={t.slice_thickness} thickness={t.thickness!r}")
    print(f"  weighted={weighted(t.array)!r} sum={float(t.array.sum())!r}")
    print(f"  project weighted={weighted(t.project().array)!r}")

atoms = bulk("Al", cubic=True) * (1, 1, 2)
potential = abtem.Potential(atoms, gpts=(16, 16), slice_thickness=2.025)
proj = potential.project()
print("Al project: shape", proj.array.shape, "sum", repr(float(proj.array.sum())), "max", repr(float(proj.array.max())),
      "weighted", repr(weighted(proj.array)))
built = potential.build()
print("Al built slices:", built.num_slices, "built project equals lazy project:",
      bool(np.array_equal(built.project().array, proj.array)))

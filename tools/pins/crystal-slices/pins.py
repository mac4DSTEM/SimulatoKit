"""SimulatoKit crystal slice thicknesses: abtem.slicing.crystal_slice_thicknesses (abTEM 1.0.10).

Prints the numbers that Tests/SimulatoMultisliceTests/CrystalSliceThicknessTests.swift holds inline (no fixture file):
the slice thicknesses of planted structures at several tolerances, and the cases abTEM refuses (its assertion that the
thicknesses sum to the cell depth). The thicknesses are float64 differences of bin centres; the test holds them at
rtol 1e-7 (assert_allclose default for float64; abTEM states no tolerance).

Run from the repo root:  .venv/bin/python -I -B tools/pins/crystal-slices/pins.py
"""

import numpy as np
from ase import Atoms
from ase.build import bulk
from abtem.slicing import crystal_slice_thicknesses


def planted(z, depth):
    n = len(z)
    positions = np.column_stack([np.full(n, 1.0), np.full(n, 2.0), np.array(z, dtype=np.float64)])
    return Atoms("C" * n, positions=positions, cell=[4.0, 4.0, depth], pbc=True)


cases = [
    ("planted 4 A, tol 0.2", planted([0.0, 0.6, 1.0, 2.05, 3.99], 4.0), 0.2),
    ("planted 4 A, tol 0.5", planted([0.0, 0.6, 1.0, 2.05, 3.99], 4.0), 0.5),
    ("planted 6 A, tol 0.2", planted([0.3, 0.6, 1.8, 3.0, 4.2, 5.4], 6.0), 0.2),
    ("planted 6 A, tol 0.7", planted([0.3, 0.6, 1.8, 3.0, 4.2, 5.4], 6.0), 0.7),
    ("Al cubic, tol 0.45", bulk("Al", cubic=True), 0.45),
    ("Al cubic, tol 0.2", bulk("Al", cubic=True), 0.2),
    ("planted 4 A, z below 0, tol 0.2", planted([-0.3, 1.0], 4.0), 0.2),
    ("planted 0.1 A, tol 0.2", planted([0.05], 0.1), 0.2),
]

for name, atoms, tol in cases:
    try:
        t = crystal_slice_thicknesses(atoms, tolerance=tol)
        print(f"{name}: z={atoms.positions[:, 2].tolist()} depth={atoms.cell[2, 2]!r}")
        print("  thicknesses:", [repr(float(v)) for v in t], "sum:", repr(float(np.sum(t))))
    except AssertionError:
        print(f"{name}: z={atoms.positions[:, 2].tolist()} depth={atoms.cell[2, 2]!r}")
        print("  refused: AssertionError (the thicknesses do not sum to the depth)")

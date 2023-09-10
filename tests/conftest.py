"""
Global pytest fixtures
"""

# Use this file if you need to share any fixtures
# across multiple modules
# More information at
# https://docs.pytest.org/en/stable/how-to/fixtures.html#scope-sharing-fixtures-across-classes-modules-packages-or-session

import MDAnalysis as mda
import numpy as np
import pytest

from mda_tui.__main__ import MDA


@pytest.fixture()
def app():
    return MDA()


@pytest.fixture(scope="session")
def universe_filenames(tmp_path_factory):
    """
    Write a dummy universe to file for testing later.

    Note, we use `tmp_path_factory` because this fixture requies `module` scope
    so we can read the file after the fixture has been created, and
    `tmp_path` has function-level scope.
    """
    n_atoms = 1
    n_frames = 1
    u = mda.Universe.empty(n_atoms, trajectory=True)
    coordinates = np.empty((n_frames, u.atoms.n_atoms, 3))
    coordinates[0] = [5, 5, 5]
    u.load_new(coordinates, order="fac")
    dims = np.asarray(
        [
            [10, 10, 10, 90, 90, 90],
        ],
    )
    u.dimensions = dims

    tmp_pdb = tmp_path_factory.getbasetemp() / "dummy_universe.pdb"
    tmp_xtc = tmp_path_factory.getbasetemp() / "dummy_universe.xtc"
    u.atoms.write(tmp_pdb.as_posix())
    with mda.Writer(tmp_xtc.as_posix()) as f:
        for _ts in u.trajectory:
            u.dimensions = dims
            f.write(u.atoms)

    return tmp_pdb, tmp_xtc

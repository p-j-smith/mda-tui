import pathlib

import MDAnalysis as mda
import numpy as np
from numpy.testing import assert_allclose
import pytest
from textual.widgets import (
    Input,
    Select,
)

from mda_tui.widgets import (
    TopologyReaderSelector,
    TrajectoryReaderSelector,
    TrajectoryWriterSelector,
    TransformationSelector,
    transformations,
)


@pytest.mark.asyncio()
async def test_translate(app, universe_filenames: tuple[pathlib.Path, pathlib.Path]):
    pdb, xtc = universe_filenames
    xtc_output = xtc.parent / "translated_trajectory.xtc"
    translate_by = 10

    async with app.run_test() as pilot:
        topology_reader_input = pilot.app.query_one(TopologyReaderSelector).query_one(Input)
        trajectory_reader_input = pilot.app.query_one(TrajectoryReaderSelector).query_one(Input)
        trajectory_writer_output = pilot.app.query_one(TrajectoryWriterSelector).query_one(Input)
        transformation_selector_input = pilot.app.query_one(TransformationSelector).query_one(
            "#transformation",
            Select,
        )
        translate_transformation = pilot.app.query_one(transformations.Translate)

        topology_reader_input.value = pdb.as_posix()
        trajectory_reader_input.value = xtc.as_posix()
        trajectory_writer_output.value = xtc_output.as_posix()
        transformation_selector_input.value = translate_transformation
        translate_transformation.query_one("#translate_x", Input).value = str(translate_by)
        translate_transformation.query_one("#translate_y", Input).value = str(translate_by)
        translate_transformation.query_one("#translate_z", Input).value = str(translate_by)

        pilot.app.run_transformation()

    u = mda.Universe(pdb.as_posix(), xtc.as_posix())
    u_translated = mda.Universe(pdb.as_posix(), xtc_output.as_posix())
    diff = (u.trajectory.timeseries() - u_translated.trajectory.timeseries()).flatten()
    assert_allclose(np.abs(diff), translate_by, atol=1e-5)


@pytest.mark.asyncio()
async def test_center_in_box(app, universe_filenames: tuple[pathlib.Path, pathlib.Path]):
    pdb, xtc = universe_filenames
    xtc_output = xtc.parent / "centered_trajectory.xtc"
    ag = "all"

    async with app.run_test() as pilot:
        topology_reader_input = pilot.app.query_one(TopologyReaderSelector).query_one(Input)
        trajectory_reader_input = pilot.app.query_one(TrajectoryReaderSelector).query_one(Input)
        trajectory_writer_output = pilot.app.query_one(TrajectoryWriterSelector).query_one(Input)
        transformation_selector_input = pilot.app.query_one(TransformationSelector).query_one(
            "#transformation",
            Select,
        )
        center_transformation = pilot.app.query_one(transformations.CenterInBox)

        topology_reader_input.value = pdb.as_posix()
        trajectory_reader_input.value = xtc.as_posix()
        trajectory_writer_output.value = xtc_output.as_posix()
        transformation_selector_input.value = center_transformation
        center_transformation.query_one("#ag", Input).value = ag

        pilot.app.run_transformation()

    u = mda.Universe(pdb.as_posix(), xtc.as_posix())
    u_centered = mda.Universe(pdb.as_posix(), xtc_output.as_posix())
    box = u.dimensions[:3]
    with pytest.raises(AssertionError):
        assert_allclose(box / 2, np.mean(u.atoms.positions, axis=0))
    assert_allclose(box / 2, np.mean(u_centered.atoms.positions, axis=0))


@pytest.mark.asyncio()
async def test_wrap(app, universe_filenames: tuple[pathlib.Path, pathlib.Path]):
    pdb, xtc = universe_filenames
    xtc_output = xtc.parent / "wrapped_trajectory.xtc"
    ag = "all"

    async with app.run_test() as pilot:
        topology_reader_input = pilot.app.query_one(TopologyReaderSelector).query_one(Input)
        trajectory_reader_input = pilot.app.query_one(TrajectoryReaderSelector).query_one(Input)
        trajectory_writer_output = pilot.app.query_one(TrajectoryWriterSelector).query_one(Input)
        transformation_selector_input = pilot.app.query_one(TransformationSelector).query_one(
            "#transformation",
            Select,
        )
        wrap_transformation = pilot.app.query_one(transformations.Wrap)

        topology_reader_input.value = pdb.as_posix()
        trajectory_reader_input.value = xtc.as_posix()
        trajectory_writer_output.value = xtc_output.as_posix()
        transformation_selector_input.value = wrap_transformation
        wrap_transformation.query_one("#ag", Input).value = ag

        pilot.app.run_transformation()

    u = mda.Universe(pdb.as_posix(), xtc.as_posix())
    u_wrapped = mda.Universe(pdb.as_posix(), xtc_output.as_posix())
    with pytest.raises(AssertionError):
        assert_allclose(u.atoms.wrap(inplace=False), u.atoms.positions)
    assert_allclose(u.atoms.wrap(inplace=False), u_wrapped.atoms.positions)

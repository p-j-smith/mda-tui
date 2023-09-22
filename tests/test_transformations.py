import pathlib

import MDAnalysis as mda
import numpy as np
from numpy.testing import assert_allclose
from packaging.specifiers import SpecifierSet
import pytest
from textual.widgets import (
    Input,
    Select,
    Switch,
)

from mda_tui.widgets import (
    TopologyReaderSelector,
    TrajectoryReaderSelector,
    TrajectoryWriterSelector,
    TransformationSelector,
    transformations,
)

mda_supported_versions = SpecifierSet(">2.6.1", prereleases=True)


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

        await pilot.app.run_transformation()
        await pilot.app.workers.wait_for_complete()

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

        await pilot.app.run_transformation()
        await pilot.app.workers.wait_for_complete()

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

        await pilot.app.run_transformation()
        await pilot.app.workers.wait_for_complete()

    u = mda.Universe(pdb.as_posix(), xtc.as_posix())
    u_wrapped = mda.Universe(pdb.as_posix(), xtc_output.as_posix())
    with pytest.raises(AssertionError):
        assert_allclose(u.atoms.wrap(inplace=False), u.atoms.positions)
    assert_allclose(u.atoms.wrap(inplace=False), u_wrapped.atoms.positions)


@pytest.mark.asyncio()
async def test_unwrap(app, universe_filenames: tuple[pathlib.Path, pathlib.Path]):
    pdb, xtc = universe_filenames
    xtc_output = xtc.parent / "unwrapped_trajectory.xtc"
    ag = "all"

    async with app.run_test() as pilot:
        topology_reader_input = pilot.app.query_one(TopologyReaderSelector).query_one(Input)
        trajectory_reader_input = pilot.app.query_one(TrajectoryReaderSelector).query_one(Input)
        trajectory_writer_output = pilot.app.query_one(TrajectoryWriterSelector).query_one(Input)
        transformation_selector_input = pilot.app.query_one(TransformationSelector).query_one(
            "#transformation",
            Select,
        )
        unwrap_transformation = pilot.app.query_one(transformations.Unwrap)

        topology_reader_input.value = pdb.as_posix()
        trajectory_reader_input.value = xtc.as_posix()
        trajectory_writer_output.value = xtc_output.as_posix()
        transformation_selector_input.value = unwrap_transformation
        unwrap_transformation.query_one("#ag", Input).value = ag

        await pilot.app.run_transformation()
        await pilot.app.workers.wait_for_complete()

    u = mda.Universe(pdb.as_posix(), xtc.as_posix())
    atom_1_position, atom_2_position = u.bonds[0].atoms.positions
    with pytest.raises(AssertionError):
        assert all(np.abs(atom_1_position - atom_2_position) < (u.dimensions[:3] / 2))

    u_unwrapped = mda.Universe(pdb.as_posix(), xtc_output.as_posix())
    atom_1_position, atom_2_position = u_unwrapped.bonds[0].atoms.positions
    assert all(np.abs(atom_1_position - atom_2_position) < (u_unwrapped.dimensions[:3] / 2))


@pytest.mark.xfail(
    mda.__version__ not in mda_supported_versions,
    reason="bug in NoJump transformation is fixed in MDAnalysis >2.6.1",
    strict=True,
)
@pytest.mark.asyncio()
async def test_nojump(app, universe_filenames: tuple[pathlib.Path, pathlib.Path]):
    pdb, xtc = universe_filenames
    xtc_output = xtc.parent / "nojump_trajectory.xtc"

    async with app.run_test() as pilot:
        topology_reader_input = pilot.app.query_one(TopologyReaderSelector).query_one(Input)
        trajectory_reader_input = pilot.app.query_one(TrajectoryReaderSelector).query_one(Input)
        trajectory_writer_output = pilot.app.query_one(TrajectoryWriterSelector).query_one(Input)
        transformation_selector_input = pilot.app.query_one(TransformationSelector).query_one(
            "#transformation",
            Select,
        )
        nojump_transformation = pilot.app.query_one(transformations.NoJump)
        nojump_transformation.query_one(Switch).value = True

        topology_reader_input.value = pdb.as_posix()
        trajectory_reader_input.value = xtc.as_posix()
        trajectory_writer_output.value = xtc_output.as_posix()
        transformation_selector_input.value = nojump_transformation

        await pilot.app.run_transformation()
        await pilot.app.workers.wait_for_complete()

    u = mda.Universe(pdb.as_posix(), xtc.as_posix())
    distances_moved = np.diff(u.trajectory.timeseries(), axis=0)
    with pytest.raises(AssertionError):
        assert np.all(np.abs(distances_moved) < (u.dimensions[:3] / 2))

    u_nojump = mda.Universe(pdb.as_posix(), xtc_output.as_posix())
    distances_moved = np.diff(u_nojump.trajectory.timeseries(), axis=0)
    assert np.all(np.abs(distances_moved) < (u_nojump.dimensions[:3] / 2))

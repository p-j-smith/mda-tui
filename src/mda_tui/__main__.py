import pathlib

import MDAnalysis as mda
import numpy as np
from textual import on, work
from textual.app import App, ComposeResult
from textual.widgets import (
    Button,
    Header,
    ProgressBar,
    Select,
    TabbedContent,
    TabPane,
)

from mda_tui.widgets import (
    MDARun,
    TopologyReaderSelector,
    TrajectoryReaderSelector,
    TrajectoryWriterSelector,
    TransformationSelector,
)


class MDA(App):
    """A TUI for MDAnalysis."""

    CSS_PATH = "mda-tui.tcss"

    def __init__(self):
        super().__init__()
        self.title = "MDA TUI"
        self.sub_title = "Trajectory transformations your terminal!"

    def compose(self) -> ComposeResult:
        """Set up the layout of the app"""
        yield Header()
        with TabbedContent():
            with TabPane("Universe"):
                yield TopologyReaderSelector()
                yield TrajectoryReaderSelector()
            with TabPane("Transformations"):
                yield TransformationSelector()
            with TabPane("Output"):
                yield TrajectoryWriterSelector()
        yield MDARun(id="run")

    @on(Button.Pressed, "TopologyReaderSelector #openFile")
    def open_topology_input(self) -> None:
        """Show the `FileOpen` dialog when the button is pushed."""
        topology_widget = self.query_one(TopologyReaderSelector)
        self.push_screen(
            topology_widget.launch_dialogue(),
            callback=topology_widget.show_selected_file,
        )

    @on(Button.Pressed, "TrajectoryReaderSelector #openFile")
    def open_trajectory_input(self) -> None:
        """Show the `FileOpen` dialog when the button is pushed."""
        trajectory_widget = self.query_one(TrajectoryReaderSelector)
        self.push_screen(
            trajectory_widget.launch_dialogue(),
            callback=trajectory_widget.show_selected_file,
        )

    @on(Button.Pressed, "TrajectoryWriterSelector #openFile")
    def open_trajectory_output(self) -> None:
        """Show the `FileOpen` dialog when the button is pushed."""
        trajectory_widget = self.query_one(TrajectoryWriterSelector)
        self.push_screen(
            trajectory_widget.launch_dialogue(),
            callback=trajectory_widget.show_selected_file,
        )

    @on(Button.Pressed, "MDARun #run-button")
    async def run_transformation(self) -> None:
        """Perform the transformation"""

        # Validate input
        topology_selector = self.query_one(TopologyReaderSelector)
        valid_topology_result = topology_selector.validate()
        if not valid_topology_result.is_valid:
            for failure in valid_topology_result.failures:
                self.notify(failure.description, severity="error", timeout=10)
            return

        trajectory_selector = self.query_one(TrajectoryReaderSelector)
        valid_trajectory_result = trajectory_selector.validate()
        if not valid_trajectory_result.is_valid:
            for failure in valid_trajectory_result.failures:
                self.notify(failure.description, severity="error", timeout=10)
            return

        output_trajectory_selector = self.query_one(TrajectoryWriterSelector)
        valid_output_trajectory_result = output_trajectory_selector.validate()
        if not valid_output_trajectory_result.is_valid:
            for failure in valid_output_trajectory_result.failures:
                self.notify(failure.description, severity="error", timeout=10)
            return

        # Load universe
        topology = topology_selector.file
        trajectory = trajectory_selector.file
        u = mda.Universe(topology.as_posix(), trajectory.as_posix())

        # Apply transformation
        transformation_wrapper = (
            self.query_one(TransformationSelector).query_one("#transformation", Select).value
        )
        if transformation_wrapper is None:
            self.notify(
                "No transformation selected - please select a transformation",
                severity="error",
                timeout=20,
            )
            return
        valid_transformation_wrapper_results = transformation_wrapper.validate()
        all_valid = all(result.is_valid for result in valid_transformation_wrapper_results)
        if not all_valid:
            for result in valid_transformation_wrapper_results:
                for failure in result.failures:
                    self.notify(failure.description, severity="error", timeout=10)
            return
        transformation = transformation_wrapper.setup_transformation(universe=u)
        u.trajectory.add_transformations(transformation)

        # Write transformed trajectory
        (start, stop, step) = u.trajectory.check_slice_indices(*self.query_one(MDARun).slice)
        output_trajectory: pathlib.Path = self.query_one(TrajectoryWriterSelector).file
        output_trajectory.parent.mkdir(parents=True, exist_ok=True)

        await self.mount(ProgressBar(np.arange(start, stop, step).size, id="pbar"))
        self._run_transformation(u, start, stop, step, output_trajectory)

    @work(exclusive=True, thread=True)
    def _run_transformation(self, u, start, stop, step, output_trajectory):
        with mda.Writer(output_trajectory.as_posix()) as f:
            for _ in u.trajectory[start:stop:step]:
                f.write(u.atoms)
                self.call_from_thread(self._update_progressbar)

    def _update_progressbar(self):
        self.query_one(ProgressBar).advance(1)

    async def on_worker_state_changed(self, event) -> None:
        """Delete the progress bar once the worker has finished."""
        if not event.worker.is_finished:
            return
        output_trajectory = event.worker._work.args[-1]
        self.notify(f"Finished writing transfomed trajectory to {output_trajectory}", timeout=20)
        await self.query_one(ProgressBar).remove()


def main():
    app = MDA()
    app.run()


if __name__ == "__main__":
    main()

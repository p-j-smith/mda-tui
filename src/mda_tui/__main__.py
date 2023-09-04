import MDAnalysis as mda
from textual import on
from textual.app import App, ComposeResult
from textual.widgets import (
    Button,
    Header,
    Input,
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

# TODO: Add validation for Input widgets


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

    @on(Button.Pressed, "MDARun #run")
    def run_transformation(self) -> None:
        """Perform the transformation"""

        # Load universe
        trajectory = self.query_one(TopologyReaderSelector).query_one(Input).value
        topology = self.query_one(TrajectoryReaderSelector).query_one(Input).value
        u = mda.Universe(topology, trajectory)

        # Apply transformation
        transformation = self.query_one(TransformationSelector).query_one(Select).value
        u.trajectory.add_transformations(transformation)

        # Write transformed trajectory
        start = self.query_one(MDARun).query_one("#start", Input).value
        stop = self.query_one(MDARun).query_one("#stop", Input).value
        step = self.query_one(MDARun).query_one("#step", Input).value
        output_trajectory = self.query_one(TrajectoryWriterSelector).query_one(Input).value
        with mda.Writer(output_trajectory) as f:
            for _ts in u.trajectory[start:stop:step]:
                f.write(u.atoms)


def main():
    app = MDA()
    app.run()


if __name__ == "__main__":
    main()

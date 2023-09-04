from textual import on
from textual.app import App, ComposeResult
from textual.widgets import (
    Button,
    Header,
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


def main():
    app = MDA()
    app.run()


if __name__ == "__main__":
    main()

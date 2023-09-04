import pathlib

import MDAnalysis as mda
from textual import (
    app,
    on,
    widgets,
)
from textual.containers import (
    Horizontal,
)
from textual.widgets import (
    Button,
    Header,
    Input,
    TabbedContent,
    TabPane,
)
import textual_fspicker


class MDA(app.App[str]):
    """A TUI for MDAnalysis."""

    CSS_PATH = "mda-tui.tcss"

    def __init__(self):
        super().__init__()
        self.title = "MDA TUI"
        self.sub_title = "Trajectory transformations your terminal!"

    def compose(self) -> app.ComposeResult:
        """Set up the layout of the app"""
        yield Header()
        with TabbedContent():
            with TabPane("Universe"):
                yield MDATopologyPicker()
                yield MDATrajectoryPicker()
            yield TabPane("Transformations")
            yield TabPane("Output")
        yield MDARun(id="run")

    @on(widgets.Button.Pressed, "MDATopologyPicker #openFile")
    def open_topology(self) -> None:
        """Show the `FileOpen` dialog when the button is pushed."""
        filters = create_file_filters(extensions=sorted(mda._PARSERS.keys()))
        self.push_screen(
            textual_fspicker.FileOpen(
                ".",
                filters=textual_fspicker.Filters(*filters),
            ),
            callback=lambda path: show_selected_file(
                widget=self.query_one(MDATopologyPicker).query_one(Input),
                path=path,
            ),
        )

    @on(widgets.Button.Pressed, "MDATrajectoryPicker #openFile")
    def open_trajectory(self) -> None:
        """Show the `FileOpen` dialog when the button is pushed."""
        filters = create_file_filters(extensions=sorted(mda._READERS.keys()))
        self.push_screen(
            textual_fspicker.FileOpen(
                ".",
                filters=textual_fspicker.Filters(*filters),
            ),
            callback=lambda path: show_selected_file(
                widget=self.query_one(MDATopologyPicker).query_one(Input),
                path=path,
            ),
        )


class MDATopologyPicker(Horizontal):
    """Widget for selecting a topology to load"""

    def compose(self) -> app.ComposeResult:
        """Create the layout for the topology picker"""
        yield Button("browse", id="openFile")
        yield Input(placeholder="select topology file")


class MDATrajectoryPicker(Horizontal):
    """Widget for selecting a trajectory to load"""

    def compose(self) -> app.ComposeResult:
        """Create the layout for the trajectory picker"""
        yield Button("browse", id="openFile")
        yield Input(placeholder="select trajectory file")


class MDARun(Horizontal):
    """Widget for selecting run parameters"""

    def compose(self) -> app.ComposeResult:
        """Create the layout for setting run parameters"""
        yield Input(placeholder="start", id="start")
        yield Input(placeholder="stop", id="stop")
        yield Input(placeholder="step", id="step")
        yield Button("run", id="run")


def create_file_filters(extensions: list[str]):
    """Create a list of file filters from on a list of file extensions.

    This is a naive approach to creating filters - each filter will show a single file
    type (with the same extension as the name of the filter)

    e.g. the 'PSF' filter will show all files with a '.psf' extension.

    The exception to this is the 'ANY' filter, which shows all files.
    """
    filters = [
        (extension, lambda file, e=extension: file.suffix.lower() == f".{e.lower()}")
        for extension in extensions
    ]
    return [("Any", lambda _: True), *filters]


def show_selected_file(widget: Input, path: pathlib.Path | None) -> None:
    """Show the file that was selected by the user.

    Args:
        widget: Input widget to update the value of
        path: The file to show.
    """
    if path is None:  # action was cancelled
        return
    widget.value = str(path)


def main():
    app = MDA()
    app.run()


if __name__ == "__main__":
    main()

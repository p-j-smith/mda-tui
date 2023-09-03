import MDAnalysis as mda
from textual import (
    app,
    containers,
    on,
    widgets,
)
import textual_fspicker


class MDA(app.App[str]):
    """A TUI for MDAnalysis."""

    # TODO: set default width for buttons and Input
    # TODO: fix spacing / layout
    # TODO: use selected file as Input text
    DEFAULT_CSS = """
    widgets.Input {
        width: 30%;
    }
    """

    def __init__(self):
        super().__init__()
        self.title = "MDA TUI"
        self.sub_title = "On-the-fly transformations and trajectory analysis in your terminal!"

    def compose(self) -> app.ComposeResult:
        """Set up the layout of the app"""
        yield widgets.Header()
        with widgets.TabbedContent():
            with widgets.TabPane("Universe", id="universe"):
                yield MDATopologyPicker()
                yield MDATrajectoryPicker()
            yield widgets.TabPane("Transformations", id="transformations")
            yield widgets.TabPane("Output", id="output")
        with containers.Horizontal():
            yield widgets.Input(placeholder="Start")
            yield widgets.Input(placeholder="Stop")
            yield widgets.Input(placeholder="Step")
        yield widgets.Button("Run")
        yield widgets.Footer()

    @on(widgets.Button.Pressed, "MDATopologyPicker #browse")
    def open_topology(self) -> None:
        """Show the `FileOpen` dialog when the button is pushed."""
        filters = [
            (parser, lambda file, p=parser: file.suffix.lower() == f".{p.lower()}")
            for parser in sorted(mda._PARSERS.keys())
        ]
        self.push_screen(
            textual_fspicker.FileOpen(
                ".",
                filters=textual_fspicker.Filters(
                    ("Any", lambda _: True),
                    *filters,
                ),
            ),
            callback=lambda _: "",
        )

    @on(widgets.Button.Pressed, "MDATrajectoryPicker #browse")
    def open_trajectory(self) -> None:
        """Show the `FileOpen` dialog when the button is pushed."""
        filters = [
            (reader, lambda file, r=reader: file.suffix.lower() == f".{r.lower()}")
            for reader in sorted(mda._READERS.keys())
        ]
        self.push_screen(
            textual_fspicker.FileOpen(
                ".",
                filters=textual_fspicker.Filters(
                    ("Any", lambda _: True),
                    *filters,
                ),
            ),
            callback=lambda _: "",
        )


class MDATopologyPicker(containers.Horizontal):
    """Widget for selecting a topology to load"""

    DEFAULT_CSS = """
    MDATopologyPicker widgets.Input {
        width: 85%;
    }
    MDATopologyPicker widgets.Button {
        width: 15%;
    }
    """

    def compose(self) -> app.ComposeResult:
        """Create the layout for the topology picker"""
        yield widgets.Button("browse", id="browse")
        yield widgets.Input(placeholder="select topology file", id="input")


class MDATrajectoryPicker(containers.Horizontal):
    """Widget for selecting a trajectory to load"""

    DEFAULT_CSS = """
    MDATrajectoryPicker widgets.Input {
        width: 85%;
    }
    MDATrajectoryPicker widgets.Button {
        width: 15%;
    }
    """

    def compose(self) -> app.ComposeResult:
        """Create the layout for the trajectory picker"""
        yield widgets.Button("browse", id="browse")
        yield widgets.Input(placeholder="select trajectory file", id="input")


def main():
    app = MDA()
    app.run()


if __name__ == "__main__":
    main()

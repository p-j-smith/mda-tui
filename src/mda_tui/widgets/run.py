from textual.app import ComposeResult
from textual.containers import Horizontal
from textual.widgets import (
    Button,
    Input,
)


class MDARun(Horizontal):
    """Widget for selecting run parameters"""

    def compose(self) -> ComposeResult:
        """Create the layout for setting run parameters"""
        yield Input(placeholder="start", id="start")
        yield Input(placeholder="stop", id="stop")
        yield Input(placeholder="step", id="step")
        yield Button("run", id="run")

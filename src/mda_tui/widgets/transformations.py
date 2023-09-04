from textual.app import ComposeResult
from textual.containers import Horizontal
from textual.widgets import OptionList


class TransformationSelector(Horizontal):
    """Widget for selecting a transformation to apply"""

    def compose(self) -> ComposeResult:
        """Create the layout for the transformation selector"""
        yield OptionList()

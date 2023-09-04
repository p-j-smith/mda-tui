from MDAnalysis.transformations import (
    TransformationBase,
    center_in_box,
    translate,
)
from textual.app import ComposeResult
from textual.containers import Horizontal
from textual.widgets import Select


class TransformationSelector(Horizontal):
    """Widget for selecting a transformation to apply"""

    def compose(self) -> ComposeResult:
        """Create the layout for the transformation selector"""
        options = [
            ("Translate coordinates by a given vector", translate),
            ("Center atoms / molecules", center_in_box),
        ]
        select: Select[TransformationBase] = Select(
            options=options,
            prompt="select transformation",
        )
        yield select

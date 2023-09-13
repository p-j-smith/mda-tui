from textual import on
from textual.app import ComposeResult
from textual.containers import Horizontal
from textual.widgets import (
    Button,
    Input,
)

from mda_tui.validators import IntegerOrNoneValidator


class MDARun(Horizontal):
    """Widget for selecting run parameters"""

    def compose(self) -> ComposeResult:
        """Create the layout for setting run parameters"""
        start = Input(
            placeholder="start",
            validators=IntegerOrNoneValidator(
                failure_description="'start' must be an integer or empty",
            ),
            id="start",
        )
        stop = Input(
            placeholder="stop",
            validators=IntegerOrNoneValidator(
                failure_description="'stop' must be an integer or empty",
            ),
            id="stop",
        )
        step = Input(
            placeholder="step",
            validators=IntegerOrNoneValidator(
                failure_description="'step' must be an integer or empty",
            ),
            id="step",
        )
        run = Button("run", id="run-button")

        # Define tooltips
        start.tooltip = "first frame to transform (0-based index). If empty, defaults to the first frame of the trajectory."
        stop.tooltip = "final frame to transform (non-inclusive). If empty, defaults to the final frame of the trajectory."
        step.tooltip = "number of frames to skip between each transformed frame. If empty, defaults to no frame skipped."
        run.tooltip = "run the transformation. This will read the input trajectory, apply the selected transformation to the selected frames and write the transformed trajecotry to the output file."

        yield start
        yield stop
        yield step
        yield run

    @on(Input.Changed)
    def show_invalid_reasons(self, event: Input.Changed) -> None:
        if not event.validation_result.is_valid:
            for failure in event.validation_result.failures:
                self.notify(failure.description, severity="error", timeout=3)
            return

    @property
    def slice(self) -> tuple[int | None, int | None, int | None]:  # noqa:  A003
        """Get the start, stop, and step values as integers or None"""

        start = self.query_one("#start", Input).value
        stop = self.query_one("#stop", Input).value
        step = self.query_one("#step", Input).value

        start = None if not start else int(start)
        stop = None if not stop else int(stop)
        step = None if not step else int(step)

        return start, stop, step

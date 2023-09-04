from textual import on
from textual.app import ComposeResult
from textual.containers import Horizontal
from textual.validation import Integer, ValidationResult
from textual.widgets import (
    Button,
    Input,
)


class IntegerOrNone(Integer):
    """Check that input is empty or can be cast to an integer"""

    def __init__(self, failure_description):
        super().__init__(failure_description=failure_description)

    def validate(self, value: str) -> ValidationResult:
        value = value if value else "0"
        return super().validate(value) if "." not in value else self.failure()


class MDARun(Horizontal):
    """Widget for selecting run parameters"""

    def compose(self) -> ComposeResult:
        """Create the layout for setting run parameters"""
        yield Input(
            placeholder="start",
            validators=IntegerOrNone(failure_description="'start' must be an integer or empty"),
            id="start",
        )
        yield Input(
            placeholder="stop",
            validators=IntegerOrNone(failure_description="'stop' must be an integer or empty"),
            id="stop",
        )
        yield Input(
            placeholder="step",
            validators=IntegerOrNone(failure_description="'step' must be an integer or empty"),
            id="step",
        )
        yield Button("run", id="run")

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

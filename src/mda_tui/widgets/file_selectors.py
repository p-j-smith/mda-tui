from pathlib import Path
from typing import Callable, ClassVar

import MDAnalysis as mda
from textual.app import ComposeResult
from textual.containers import (
    Horizontal,
)
from textual.validation import Validator
from textual.widgets import (
    Button,
    Input,
)
import textual_fspicker
from typing_extensions import TypeAlias

from mda_tui.validators import (
    FileExtensionValidator,
    FileValidator,
)
from mda_tui.widgets.file_dialogues import FileOpen, FileSave

FilterFunction: TypeAlias = Callable[[Path], bool]
FileFilters: TypeAlias = list[tuple[str, FilterFunction]]


class FileSelector(Horizontal):
    """Widget for selecting a file using a dialogue box or Input widget

    Attributes:
        filters: List of filters to use when selecting files
        input_placeholder: String to show in empty Input widget
        input_id: Textual ID to use for the Input widget
        button_id: Textual ID to use for the Button widget
    """

    filters: FileFilters | None = None
    placeholder: str = "select file"
    input_id: str | None = None
    button_id: str | None = "openFile"
    validators: Validator | list[Validator] | None = None

    def compose(self) -> ComposeResult:
        """Create the layout for the file selector"""
        yield Button("browse", id=self.button_id)
        yield Input(placeholder=self.placeholder, validators=self.validators, id=self.input_id)

    @classmethod
    def create_file_filters(cls, extensions: list[str]) -> FileFilters:
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

    def launch_dialogue(self):
        """Open a file dialogue box"""
        return FileOpen(
            ".",
            filters=textual_fspicker.Filters(*self.filters),
        )

    def show_selected_file(self, path):
        """Display the selected file in the Input widget"""
        self.query_one(Input).value = "" if path is None else str(path)

    @property
    def file(self) -> Path | None:
        return Path(self.query_one(Input).value)

    def validate(self):
        input = self.query_one(Input)  # noqa: A001
        return input.validate(input.value)


class TopologyReaderSelector(FileSelector):
    """Widget for selecting a topology to load"""

    filters: FileFilters = FileSelector.create_file_filters(
        extensions=sorted(mda._PARSERS.keys()),
    )
    placeholder: str = "select topology file"
    validators: ClassVar = [
        FileValidator(failure_description="Input topology file does not exist"),
        FileExtensionValidator(
            valid_extensions=sorted(mda._PARSERS.keys()),
            failure_description="Unknown input topology format",
        ),
    ]
    tooltip = "select a valid topology file to load."


class TrajectoryReaderSelector(FileSelector):
    """Widget for selecting a trajectory to load"""

    filters: FileFilters = FileSelector.create_file_filters(
        extensions=sorted(mda._READERS.keys()),
    )
    placeholder: str = "select trajectory file"
    validators: ClassVar = [
        FileValidator(failure_description="Input trajectory file does not exist"),
        FileExtensionValidator(
            valid_extensions=sorted(mda._READERS.keys()),
            failure_description="Unknown input trajectory format",
        ),
    ]
    tooltip = "select a valid trajectory file to load."


class TrajectoryWriterSelector(FileSelector):
    """Widget for selecting a trajectory to write"""

    filters: FileFilters = FileSelector.create_file_filters(
        extensions=sorted(mda._MULTIFRAME_WRITERS.keys()),
    )
    placeholder: str = "select transformed trajectory file"
    validators: ClassVar = [
        FileExtensionValidator(
            valid_extensions=sorted(mda._MULTIFRAME_WRITERS.keys()),
            failure_description="Unknown output trajectory format",
        ),
    ]
    tooltip = "select file to write the transformed trajetory to. Any existing file of the same name will be overwritten."

    def launch_dialogue(self):
        """Open a file dialogue box"""
        return FileSave(
            location=".",
            filters=textual_fspicker.Filters(*self.filters),
        )

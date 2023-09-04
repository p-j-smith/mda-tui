from pathlib import Path
from typing import Callable, ClassVar

import MDAnalysis as mda
from textual.app import ComposeResult
from textual.containers import (
    Horizontal,
)
from textual.validation import Failure, ValidationResult, Validator
from textual.widgets import (
    Button,
    Input,
)
import textual_fspicker
from typing_extensions import TypeAlias

from mda_tui.widgets.file_dialogues import FileOpen, FileSave

FilterFunction: TypeAlias = Callable[[Path], bool]
FileFilters: TypeAlias = list[tuple[str, FilterFunction]]


class File(Validator):
    """Check that input is a valid file."""

    def __init__(
        self,
        failure_description: str | None = None,
    ) -> None:
        super().__init__(failure_description=failure_description)

    class InvalidFile(Failure):
        """Indicates that the file does not exist."""

    def validate(self, value: str) -> ValidationResult:
        """Validates that `value` is a valid file (i.e. it exists).

        Args:
            value: The value to validate.

        Returns:
            The result of the validation.
        """
        invalid_file = ValidationResult.failure([File.InvalidFile(self, value)])
        file = Path(value).resolve()
        if not file.exists() or not file.is_file():
            return invalid_file
        return self.success()


class FileExtension(Validator):
    """Check that input has a valid file extension."""

    def __init__(
        self,
        failure_description: str | None = None,
        valid_extensions: list[str | Path] | None = None,
    ) -> None:
        super().__init__(failure_description=failure_description)
        self.valid_extensions = valid_extensions

    class InvalidFileExtension(Failure):
        """Indicates that the file extension is invalid."""

    def validate(self, value: str) -> ValidationResult:
        """Validates that `value` is a valid file extension.

        Args:
            value: The value to validate.

        Returns:
            The result of the validation.
        """
        invalid_extension = ValidationResult.failure(
            [FileExtension.InvalidFileExtension(self, value)],
        )
        extension = Path(value).suffix.lstrip(".")
        if self.valid_extensions is not None and extension.upper() not in self.valid_extensions:
            return invalid_extension
        return self.success()


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
        File(failure_description="Input topology file does not exist"),
        FileExtension(
            valid_extensions=sorted(mda._PARSERS.keys()),
            failure_description="Unknown input topology format",
        ),
    ]


class TrajectoryReaderSelector(FileSelector):
    """Widget for selecting a trajectory to load"""

    filters: FileFilters = FileSelector.create_file_filters(
        extensions=sorted(mda._READERS.keys()),
    )
    placeholder: str = "select trajectory file"
    validators: ClassVar = [
        File(failure_description="Input trajectory file does not exist"),
        FileExtension(
            valid_extensions=sorted(mda._READERS.keys()),
            failure_description="Unknown input trajectory format",
        ),
    ]


class TrajectoryWriterSelector(FileSelector):
    """Widget for selecting a trajectory to write"""

    filters: FileFilters = FileSelector.create_file_filters(
        extensions=sorted(mda._MULTIFRAME_WRITERS.keys()),
    )
    placeholder: str = "select transformed trajectory file"
    validators: ClassVar = [
        FileExtension(
            valid_extensions=sorted(mda._MULTIFRAME_WRITERS.keys()),
            failure_description="Unknown output trajectory format",
        ),
    ]

    def launch_dialogue(self):
        """Open a file dialogue box"""
        return FileSave(
            location=".",
            filters=textual_fspicker.Filters(*self.filters),
        )

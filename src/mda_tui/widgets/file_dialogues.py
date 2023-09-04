import pathlib

from textual_fspicker.file_dialog import BaseFileDialog
from textual_fspicker.path_filters import Filters


class FileOpen(BaseFileDialog):
    """A file opening dialog.

    This is the same as textual_fspicker.FileOpen except the 'select_button' text
    is changed from 'Open' to 'Select'.
    """

    def __init__(
        self,
        location: str | pathlib.Path = ".",
        title: str = "Open",
        *,
        select_button="Select",
        filters: Filters | None = None,
        must_exist: bool = True,
    ) -> None:
        """Initialise the `FileOpen` dialog.

        Args:
            location: Optional starting location.
            title: Optional title.
            filters: Optional filters to show in the dialog.
            must_exist: Flag to say if the file must exist.
        """
        super().__init__(location, title, select_button=select_button, filters=filters)
        self._must_exist = must_exist
        """Must the file exist?"""

    def _should_return(self, candidate: pathlib.Path) -> bool:
        """Perform the final checks on the chosen file.

        Args:
            candidate: The file to check.
        """
        if self._must_exist and not candidate.exists():
            self._set_error("The file must exist")
            return False
        return True


class FileSave(BaseFileDialog):
    """A file save dialog.

    This is the same as textual_fspicker.FileSave except the 'select_button' text
    is changed from 'Save' to 'Select'.
    """

    def __init__(
        self,
        location: str | pathlib.Path = ".",
        title: str = "Save as",
        *,
        filters: Filters | None = None,
        can_overwrite: bool = True,
        select_button="Select",
    ) -> None:
        """Initialise the `FileOpen` dialog.

        Args:
            location: Optional starting location.
            title: Optional title.
            filters: Optional filters to show in the dialog.
            can_overwrite: Flag to say if an existing file can be overwritten.
        """
        super().__init__(location, title, select_button=select_button, filters=filters)
        self._can_overwrite = can_overwrite
        """Can an existing file be overwritten?"""

    def _should_return(self, candidate: pathlib.Path) -> bool:
        """Perform the final checks on the chosen file.

        Args:
            candidate: The file to check.
        """
        if candidate.exists() and not self._can_overwrite:
            self._set_error("Overwrite is not allowed")
            return False
        return True

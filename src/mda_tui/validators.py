from pathlib import Path

import MDAnalysis as mda
from textual.validation import (
    Failure,
    Integer,
    Number,
    ValidationResult,
    Validator,
)


class NumberOrNoneValidator(Number):
    """Check that input is empty or is a number"""

    def __init__(self, failure_description):
        super().__init__(failure_description=failure_description)

    def validate(self, value: str) -> ValidationResult:
        value = value if value else "0"
        return super().validate(value)


class IntegerOrNoneValidator(Integer):
    """Check that input is empty or can be cast to an integer"""

    def __init__(self, failure_description):
        super().__init__(failure_description=failure_description)

    def validate(self, value: str) -> ValidationResult:
        value = value if value else "0"
        return super().validate(value) if "." not in value else self.failure()


class FileValidator(Validator):
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
        invalid_file = ValidationResult.failure([FileValidator.InvalidFile(self, value)])
        file = Path(value).resolve()
        if not file.exists() or not file.is_file():
            return invalid_file
        return self.success()


class FileExtensionValidator(Validator):
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
            [FileExtensionValidator.InvalidFileExtension(self, value)],
        )
        extension = Path(value).suffix.lstrip(".")
        if self.valid_extensions is not None and extension.upper() not in self.valid_extensions:
            return invalid_extension
        return self.success()


class AtomSelectionValidator(Validator):
    """Check that input is a valid selection string"""

    def __init__(
        self,
        failure_description: str | None = None,
    ) -> None:
        super().__init__(failure_description=failure_description)

    class InvalidAtomSelection(Failure):
        """Indicates that the selection string is invalid."""

    def validate(self, value: str) -> ValidationResult:
        """Validates that `value` is a valid selection strin (i.e. it returns an AtomGroup).

        Args:
            value: The value to validate.

        Returns:
            The result of the validation.
        """
        failure = AtomSelectionValidator.InvalidAtomSelection(self, value)
        invalid_selection = ValidationResult.failure([failure])
        value = value if value else "all"
        try:
            mda.core.selection.Parser.parse(selectstr=value, selgroups={})
        except mda.exceptions.SelectionError as e:
            failure.description = f"Invalid atom selection - {e}"
            return invalid_selection
        except mda.exceptions.NoDataError as e:
            failure.description = f"Invalid atom selection - {e}"
            return invalid_selection
        except AttributeError as e:
            failure.description = f"Invalid atom selection - {e}"
            return invalid_selection
        else:
            return self.success()

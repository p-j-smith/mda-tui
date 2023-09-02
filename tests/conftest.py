"""
Global pytest fixtures
"""

# Use this file if you need to share any fixtures
# across multiple modules
# More information at
# https://docs.pytest.org/en/stable/how-to/fixtures.html#scope-sharing-fixtures-across-classes-modules-packages-or-session

import pytest

from mda_tui.data.files import MDANALYSIS_LOGO


@pytest.fixture()
def mdanalysis_logo_text() -> str:
    """Example fixture demonstrating how data files can be accessed"""
    with open(MDANALYSIS_LOGO, encoding="utf8") as f:  # noqa: PTH123
        logo_text = f.read()
    return logo_text

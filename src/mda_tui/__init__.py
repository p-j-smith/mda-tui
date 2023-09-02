"""
MDAnalysis TUI
A textual user interface (TUI) for MDAnalysis
"""

# Add imports here

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

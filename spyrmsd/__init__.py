"""
Python RMSD tool with symmetry correction.
"""

from ._version import get_versions
from .due import Doi, due

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

# This will print latest Zenodo version
due.cite(
    Doi("10.5281/zenodo.3631876"),
    path="spyrmsd",
    description="spyrmsd",
)

due.cite(
    Doi("10.1186/s13321-020-00455-2"),
    path="spyrmsd",
    description="spyrmsd",
)

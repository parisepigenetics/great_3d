"""The GREAT 3D package for multi-3D (epi)genomics."""

# Add imports here
from .transcriptome_3d import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

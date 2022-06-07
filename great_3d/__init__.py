"""The GREAT 3D suite for integrative multi-omics 3D epigenomics.
"""

# Add imports here
from .transcriptome_3d import *
from .great_pdb import *

# Handle versioneer
from ._version import get_versions
versions = get_versions()
__version__ = versions['version']
__git_revision__ = versions['full-revisionid']
del get_versions, versions

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions

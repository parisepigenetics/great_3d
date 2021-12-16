"""
Unit and regression test for the great_3d package.
"""

# Import package, test suite, and other packages as needed
import sys

import pytest

import great_3d


def test_great_3d_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "great_3d" in sys.modules

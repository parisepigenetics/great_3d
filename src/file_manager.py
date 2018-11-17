#!/usr/bin/env python3
"""This a module to manage files and open or read it properly"""

from pathlib import Path
import sys


def check_file(path):
    """Function check if a path exist, and if it point on a file.

    ARGUMENTS:

    path = path to the file, str.

    RETURN:
    the path if it's ok.
    """
    file = Path(path)
    if file.exists():
        if file.is_file():
            return file
        sys.exit("Your path: '{}' doesn't point on a file!".format(path))
    sys.exit("Your Path: '{}' is wrong, file doesn't exist!".format(path))

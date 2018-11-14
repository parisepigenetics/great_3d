#!/usr/bin/env python3
"""This a module to manage files and open or read it properly"""

from pathlib import Path
import sys

def checkFile(path):
	file = Path(path)
	if(file.exists()):
		if(file.is_file()):
			return file
		else:
			sys.exit("Your path: '{}' doesn't point on a file!".format(path))
	else:
		sys.exit("Your Path: '{}' is wrong, file doesn't exist!".format(path))


"""Module with the main classes/functions to interface PDB formt with GREAT."""

import numpy as np
import pandas as pd
from Bio import PDB
from scipy import stats
from multiprocessing import cpu_count, Pool, Process
from functools import partial
# Bottleneck can speed up SOME functions on numpy arrays https://bottleneck.readthedocs.io/en/latest/intro.html
import bottleneck

import time
import pprint  # for testing only!

def timing(f):
    """Wrapper to time functions.py

    Works as a decorator. Taken from https://stackoverflow.com/questions/5478351/python-time-measure-function
    """
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print("{:s} function took {:.3f} ms".format(f.__name__, (time2 - time1) * 1000.0))
        return ret
    return wrap


# TODO under developemnt for parsing (IFF REQUIRED) .pdb genome structure files.
# parser = PDB.PDBParser()
# io = PDB.PDBIO()
# struct = parser.get_structure("chr3", "chr3.pdb")
#
# for model in struct:
#     for chain in model:
#         for residue in chain:
#             for atom in residue:
#                 na =  atom.get_full_id()
#                 x,y,z = atom.get_coord()
#                 print(na[0], na[3][1], x, y, z)

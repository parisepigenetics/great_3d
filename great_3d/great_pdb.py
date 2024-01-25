"""Module to interface classes/functions of genome 3D PDB with GREAT."""


import numpy as np
import pandas as pd
from Bio import PDB
from scipy import stats
from multiprocessing import cpu_count, Pool, Process
from functools import partial
# Bottleneck can speed up SOME functions on numpy arrays https://bottleneck.readthedocs.io/en/latest/intro.html
import bottleneck

from great_3d import timing
import pprint  # for testing only!


@timing
# TODO under developemnt for parsing (IFF REQUIRED) .pdb genome structure files.
# TODO create a PDB genome class perhaps
parser = PDB.PDBParser()
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

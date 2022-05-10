"""Module with the main classes/functions to interface PDB formt with GREAT.
"""

import numpy as np
import pandas as pd
import scipy
from scipy import stats
from scipy.spatial import distance_matrix
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


def aaa():
    """
    """
    pass

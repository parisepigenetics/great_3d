# -*- coding: utf-8  -*-

"""Main multi epigenetics viewer moddule.
"""

import math
import time
import numpy as np
import scipy
from multiprocessing import cpu_count, Pool
import pandas as pd
from scipy.spatial import distance_matrix


def matrix_distance(data):
    # TODO parallelise the calculation of the distance matrix.
    """Gets a position filename, return the distance matrix as a Pandas DataFrame.
    """
    df = pd.read_csv(data, sep='\t')
    geneNames = df.index
    ndarray = scipy.spatial.distance.pdist(df)
    matrix_uni = scipy.spatial.distance.squareform(ndarray)
    matrix_dist = pd.DataFrame(matrix_uni)
    matrix_dist.columns = geneNames
    matrix_dist.index = geneNames
    return matrix_dist

def dico_matrix(matrix):
    '''Get a distance matrix (Pandas DataFrame) and return a dictionary of sorted genes with each gene as a key.
    '''
    sortingDict = {}
    for gene_name in matrix.index:
        sortingDict[gene_name] = matrix[gene_name].sort_values()
    return sortingDict

def parallelise_sorting_matrix(df) :
    '''Running the dico_matrix function using Multiprocessing
    '''
    #cpus = cpu_count() - 1
    cpus = 1
    with Pool(cpus) as pool:
        # Split the data frame
        df_split = np.array_split(df, cpus)
        sorted_list_df = pool.map(dico_matrix, df_split)
        pool.close()
        pool.join()
    return sorted_list_df

## Nice wrapper to time functions. Works as a decorator.
# Taken from https://stackoverflow.com/questions/5478351/python-time-measure-function
def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print('{:s} function took {:.3f} ms'.format(f.__name__, (time2-time1)*1000.0))
        return ret
    return wrap

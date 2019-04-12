# -*- coding: utf-8  -*-

"""Main multi epigenetics viewer moddule.
"""

import math
import numpy as np
import scipy
from multiprocessing import Pool
import pandas as pd
from scipy.spatial import distance_matrix

def matrix_distance(data):
    """Gets a position filename, return the distance matrix as a Pandas DataFrame.
    """
    df = pd.read_csv(data, sep='\t')
    geneNames = df.index
    ndarray = scipy.spatial.distance.pdist(df)
    matrix_uni = scipy.spatial.distance.squareform(ndarray)
    matrix_dist = pd.DataFrame(matrix_uni)
    matrix_dist.columns = geneNames
    matrix_dist.index = geneNames
    #for i in range(len(matrix_dist)):  # No need to do this copmplicated loop.
    #    matrix_dist = matrix_dist.rename(columns={i:df.index[i]})
    #    matrix_dist = matrix_dist.rename(index={i:df.index[i]})
    return matrix_dist

def dico_matrix(matrix):
    '''Get a distance matrix (Pandas DataFrame) and return a dictionary of sorted genes with each gene as a key.
    '''
    sortingDict = {}
    # TODO put multiprocessing and parallelise this for loop ONLY (for the moment).
    for i in range(len(matrix)):  # Loop over the data frame row names (i.e. gene names).
        sortingDict[matrix.index[i]] = matrix[matrix.index[i]].sort_values()
    return sortingDict

# If you want put it here.
#print('le nombre de processeur utilisé est de : ', multiprocessing.cpu_count())

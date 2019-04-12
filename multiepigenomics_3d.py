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
    return matrix_dist

def dico_matrix(matrix):
    '''Get a distance matrix (Pandas DataFrame) and return a dictionary of sorted genes with each gene as a key.
    '''
    sortingDict = {}
    for gene_name in matrix.index:
        sortingDict[gene_name] = matrix[gene_name].sort_values()
    return sortingDict

def parallelization_dico_matrix(matrix) :
    '''Running the dico_matrix function using Multiprocessing
    '''
    with Pool(5) as p :
    	sortingDict_multi = p.map(dico_matrix,[matrix])
    return sortingDict_multi

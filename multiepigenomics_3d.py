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


def calculate_distance(position_file):
    # TODO parallelise the calculation of the distance matrix.
    """Get a position file, return the distance matrix as a Pandas DataFrame.
    """
    df = pd.read_csv(position_file, sep='\t')
    geneNames = df.index
    ndarray = scipy.spatial.distance.pdist(df)
    matrix_uni = scipy.spatial.distance.squareform(ndarray)
    matrix_dist = pd.DataFrame(matrix_uni)
    matrix_dist.columns = geneNames
    matrix_dist.index = geneNames
    return matrix_dist


def sorting_distances(dist_df):
    '''Get a distance matrix (Pandas DataFrame) and return a dictionary of sorted genes with each gene as a key.
    '''
    sortedDict = {}
    for gene_name in dist_df.index:
        sortedDict[gene_name] = dist_df.loc[gene_name,].sort_values()
    return sortedDict


def sum_Correlation(distance_dico , correlation_matrix):
    """Gets the Dictionnary of the closest_gene and the correlation Matrix, return a Dictionnary of the sum of correlation for each gene
    """
    dico_sum_Correlation = {}
    for gene , closest_gene in distance_dico.items():
        sum_correlation = abs(correlation_matrix.loc[gene,closest_gene].sum())
        dico_sum_Correlation[gene] = sum_correlation
    return dico_sum_Correlation



#================ OBSOLETE functions ======================

def sorting_dists_parallel(dist_df) :
    '''Get a distance matrix (pandas data frame), return a dictionary of genes and their sorted neighbours.
    '''
    cpus = cpu_count() - 1
    sorted_dict = {}
    with Pool(cpus) as pool:
        # Split the data frame
        df_split = np.array_split(dist_df, cpus)
        sorted_dfs = pool.map(sorting_distances, df_split)
        pool.close()
        pool.join()
    # Put together all the results.
    for d in sorted_dfs: sorted_dict.update(d)
    return sorted_dict

def closest_gene_name(dico_closest_gene):
    '''Get the dictionnary of sorted genes and return only the names of the N closest gene with each gene as a key
    '''
    dico_name_closest_gene = {}
    for gene_name in dico_closest_gene :
        closest_genes = list(pd.DataFrame(dico_closest_gene[gene_name][1:]).index)
        dico_name_closest_gene[gene_name] = closest_genes
    return dico_name_closest_gene

def correlation_matrix(gene_expression_file):
    """Gets the gene Expression filename, return the Correlation matrix as a Pandas DataFrame.
    """
    gene_expression_df = pd.read_csv(gene_expression_file , sep ='\t').transpose()
    return gene_expression_df.corr(method = 'spearman') #or Pearson

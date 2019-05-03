# -*- coding: utf-8  -*-

"""Module with all the multi-omics epigenetic analysis functions.
"""

import math
import time
import numpy as np
import scipy
from multiprocessing import cpu_count, Pool , Process
from functools import partial
import pandas as pd
from scipy.spatial import distance_matrix
from scipy.stats import pearsonr , kendalltau , spearmanr


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


def sum_correlation(sorted_dists, ge_file, no_genes , type_correlation):
    """Gets the dictionnary of the closest genes, the gene expression file and the number of genes we need to compute correlation.
    Return a dictionnary of the sum of correlations for each gene
    """
    # Read the GE file
    geDF = pd.read_csv(ge_file, sep='\t')
    # Convert the GE data frame to a dictionary.
    geD = geDF.T.to_dict('list')
    correlation_sums = {}
    # Here is the actual calculation.
    for gene_ref, closest_genes in sorted_dists.items():
        # TODO check if we gain time when we paralelise this for loop!
        selected_genes = list(closest_genes[1:no_genes + 1].index)
        ref_GE = geD[gene_ref]
        if(type_correlation == 'pearson'):
            correlations = [pearsonr(ref_GE, geD[s])[0] for s in selected_genes]
        elif(type_correlation == 'kendall'):
            correlations = [kendalltau(ref_GE, geD[s])[0] for s in selected_genes]
        elif(type_correlation == 'spearman'):
            correlations = [spearmanr(ref_GE, geD[s])[0] for s in selected_genes]

        # NOTE incase we need to use the absolute values we can flip to that.  correlations = [abs(pearsonr(ref_GE, geD[s])[0]) for s in selected_genes]
        # NOTE The pearsonr method returns the p-value for a 2 tauiled-correlation test, so perhpas we can use it in a later analysis.
        # TODO the function need to be modular in the sense that we can pass different correlation functions to it. Perhaps use decorators.
        sum_correlation = sum(correlations)
        correlation_sums[gene_ref] = sum_correlation
    return correlation_sums

#================ OBSOLETE functions ======================

def sumCor_mp(sorted_dists, ge_file, no_genes , type_correlation) :
    cpus = cpu_count() - 1
    dict = {}
    with Pool(cpus) as pool:
        dict = pool.map(sum_correlation, (sorted_dists , ge_file , no_genes , type_correlation))
        pool.close()
        pool.join()
    return dict


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
        closest_genes = list(pd.DataFrame(dico_closest_gene[gene_name][1:]).index) # This computes correlations between the gene of reference and ALL the other genes.
        dico_name_closest_gene[gene_name] = closest_genes
    return dico_name_closest_gene

def correlation_matrix(gene_expression_file):
    """Gets the gene Expression filename, return the Correlation matrix as a Pandas DataFrame.
    """
    gene_expression_df = pd.read_csv(gene_expression_file , sep ='\t').transpose()
    return gene_expression_df.corr(method = 'spearman') #or Pearson

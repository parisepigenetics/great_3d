#! /usr/bin/python3
"""
program that creates distance and correlation matrix
"""

import scipy.spatial


def distance_matrix (position_data_frame):
    """
    function that creates a distance matrix
    """
    return scipy.spatial.distance.cdist(position_data_frame, position_data_frame, metric='euclidean')


def correlation_matrix (expression_data_frame):
    """
    function that creates a correlation matrix
    """
    print("processing the creation of the correlation matrix this may take several minutes")
    return expression_data_frame.corr(method = 'spearman')


def close_genes_correlation (dict_matrix, corr_matrix, nbr_gen, overr):
    """
    function that gives the correlation sum of the nbr_gen closest genes for each gene
    """
    SUM_CORR = {}
    for gene in overr :
        closeGenes = dict_matrix.loc[gene,].sort_values()[1:nbr_gen+1].index.tolist()
        gene3Dcorr = abs(corr_matrix.loc[gene,closeGenes].sum())
        SUM_CORR[gene]=gene3Dcorr
    return(SUM_CORR)


if __name__ == "__main__":

	print("this programm doesn't run independently sorry !")

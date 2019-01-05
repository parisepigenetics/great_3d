#! /usr/bin/python3
"""
Genes position and expression files Parsing
"""
import argparse
import pandas as pd
import scipy.spatial
import multiprocessing as mp
from functools import partial


def gene_position (position_file):
    """
    
    """
    return pd.read_csv(position_file, delim_whitespace = True)


def gene_expression (expression_file):
    """

    """
    return pd.read_csv(expression_file, delim_whitespace = True)


def overlaping_genes (position_data_frame, expression_data_frame):
    """
    
    """
    return expression_data_frame.intersection(position_data_frame)


def distance_matrix (position_data_frame):
    """
    
    """
    return scipy.spatial.distance.cdist(position_data_frame, position_data_frame, metric='euclidean')


def correlation_matrix (expression_data_frame):
    """
    
    """
    return expression_data_frame.corr(method = 'spearman')


def close_genes_correlation (list_data_nbr, gene_name):
    """

    """
    closeGenes = list_data_nbr[0].loc[gene_name,].sort_values()[1:list_data_nbr[2]+1].index.tolist()
    gene3Dcorr = abs(list_data_nbr[1].loc[gene_name,closeGenes].sum())
    #dataframe.loc[row,col]
    #list[0]=dis  list[1]=corr  list[2]=n
    list_data_nbr[-1][gene_name]=gene3Dcorr
    return(list_data_nbr[-1])



if __name__ == "__main__":


    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("pos_file", help="the file containing the x y z coordinates of the genes", type=str)
    PARSER.add_argument("exp_file", help="the file containing the genes expression data", type=str)
    PARSER.add_argument("nbr_genes", help="the number of close genes", type=int)
    ARGS = PARSER.parse_args()
    POS_FILE = ARGS.pos_file
    EXP_FILE = ARGS.exp_file
    NBR_GEN = ARGS.nbr_genes
    
    NB_CPU = mp.cpu_count()
    POOL = mp.Pool(NB_CPU)


    GENE_POS = gene_position(POS_FILE)
    #print(GENE_POS)

    GENE_EXP = gene_expression(EXP_FILE)
    #print(GENE_EXP)

    overr=overlaping_genes(GENE_POS.index, GENE_EXP.index)
    overr=overr.tolist()
    #print(len(overr))
    #for i in overr :
    #    print(GENE_POS.loc[i,].tolist())

    GEN_POS = pd.DataFrame(GENE_POS.loc[overr])
    GENE_EXP = pd.DataFrame(GENE_EXP.loc[overr])

    DIST_MATRIX = pd.DataFrame(distance_matrix(GEN_POS[['X','Y','Z']]))
    DIST_MATRIX.index = GEN_POS.index
    DIST_MATRIX.columns = GEN_POS.index
    print(DIST_MATRIX)

    CORR_MATRIX = correlation_matrix(GENE_EXP.transpose())


    dic = {}


    #print(len(overr))
    func = partial(close_genes_correlation, [DIST_MATRIX, CORR_MATRIX, NBR_GEN, dic])
    dic = POOL.map(func, overr)

    #print(len(dic))
    dic = dic [0]
    #print(type(dic))
    dic2 = pd.DataFrame([dic], columns=dic.keys())
    print(dic2)
    #dic2 = pd.DataFrame.from_dict(dic, orient='index')
    #print(dic2)


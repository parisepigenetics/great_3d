#! /usr/bin/python3
"""
Genes position and expression files Parsing
"""
import argparse
import pandas as pd
import scipy.spatial

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


if __name__ == "__main__":

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("pos_file", help="the file containing the x y z coordinates of the genes", type=str)
    PARSER.add_argument("exp_file", help="the file containing the genes expression data", type=str)
    ARGS = PARSER.parse_args()
    POS_FILE = ARGS.pos_file
    EXP_FILE = ARGS.exp_file
    
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

    print(correlation_matrix(GENE_EXP))

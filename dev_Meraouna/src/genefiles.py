#! /usr/bin/python3
"""
gene expression and position files management program
"""

import argparse
import pandas as pd


def gene_data (gene_file):
    """
    function that reads a gene position or gene expression file and
    returns a dataframe
    """
    return pd.read_csv(gene_file, delim_whitespace = True)


def overlaping_genes (position_data_frame, expression_data_frame):
    """
    function that returns the overlaping genes between the two data frames
    """
    return expression_data_frame.intersection(position_data_frame)


if __name__ == "__main__":

    PARSER = argparse.ArgumentParser()

    PARSER.add_argument("pos_file", help="the file containing the x y z coordinates of the genes", type=str)
    PARSER.add_argument("exp_file", help="the file containing the genes expression data", type=str)
    PARSER.add_argument("nbr_genes", help="the number of close genes to select", type=int)

    ARGS = PARSER.parse_args()
    
    POS_FILE = ARGS.pos_file
    EXP_FILE = ARGS.exp_file
    NBR_GEN = ARGS.nbr_genes

    GENE_POS = gf.gene_data(POS_FILE)
    print(GENE_POS)
    GENE_EXP = gf.gene_data(EXP_FILE)
    print(GENE_EXP)

    overr = gf.overlaping_genes(GENE_POS.index, GENE_EXP.index)
    overr = overr.tolist()
    print(overr)

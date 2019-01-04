#! /usr/bin/python3
"""
Genes position and expression files Parsing
"""
import argparse
import pandas as pd

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



if __name__ == "__main__":

    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("pos_file", help="the file containing the x y z coordinates of the genes", type=str)
    PARSER.add_argument("exp_file", help="the file containing the genes expression data", type=str)
    ARGS = PARSER.parse_args()
    POS_FILE = ARGS.pos_file
    EXP_FILE = ARGS.exp_file
    
    GENE_POS = gene_position(POS_FILE)
    print(GENE_POS)

    GENE_EXP = gene_expression(EXP_FILE)
    print(GENE_EXP)

    overr=overlaping_genes(GENE_POS.index, GENE_EXP.index)
    overr=overr.tolist()
    print(type(overr))
    print(len(overr))
    for i in overr :
        print(GENE_POS.loc[i,].tolist())

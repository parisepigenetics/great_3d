#! /usr/bin/python3
"""
the main program of the project
"""

import argparse
import genefiles as gf
import distcorr as dc
import visual as vsl
import pandas as pd
import numpy as np
import time
from os import makedirs, path


def complete_time(start, end):
    """
    function that gives the time duration in hours minutes and seconds
    """
    total = end - start
    hrs = total // 3600
    minu = (total % 3600) // 60
    sec = (total % 3600) % 60
    if total < 60:
        return "{:2.0f} sec".format(sec)
    elif total < 3600:
        return "{:2.0f} min {:2.0f} sec".format(minu, sec)
    else:
        return "{:3.0f} h {:2.0f} min {:2.0f} sec".format(hrs, minu, sec)


if __name__ == "__main__":
    
    START = time.time()

    PARSER = argparse.ArgumentParser()

    PARSER.add_argument("pos_file", help="the file containing the x y z coordinates of the genes", type=str)
    PARSER.add_argument("exp_file", help="the file containing the genes expression data", type=str)
    PARSER.add_argument("nbr_genes", help="the number of close genes to select", type=int)

    ARGS = PARSER.parse_args()

    POS_FILE = ARGS.pos_file
    EXP_FILE = ARGS.exp_file
    NBR_GEN = ARGS.nbr_genes

    # get the data from the files
    GENE_POS = gf.gene_data(POS_FILE)
    GENE_EXP = gf.gene_data(EXP_FILE)

    # get the overlaping genes
    OVERR = gf.overlaping_genes(GENE_POS.index, GENE_EXP.index)
    OVERR = OVERR.tolist()

    # keeps only the overlaping genes
    GEN_POS = pd.DataFrame(GENE_POS.loc[OVERR])
    GENE_EXP = pd.DataFrame(GENE_EXP.loc[OVERR])

    # calculating the distance matrix
    DIST_MATRIX = pd.DataFrame(dc.distance_matrix(GEN_POS[['X','Y','Z']]))
    DIST_MATRIX.index = GEN_POS.index
    DIST_MATRIX.columns = GEN_POS.index

    # calculating the correlation matrix
    CORR_MATRIX = dc.correlation_matrix(GENE_EXP.transpose())

    # getting the correlations sum of the closest genes for each gene
    SUM_CORR = dc.close_genes_correlation (DIST_MATRIX, CORR_MATRIX, NBR_GEN, OVERR)

    SUM_CORR = pd.DataFrame([SUM_CORR], columns=SUM_CORR.keys())
    SUM_CORR.rename(index = {0: "sum_corr"}, inplace = True)

    # joining to get the position and correlations sum data frame
    TRANSMAP3D = SUM_CORR.transpose().join(GEN_POS[['X','Y','Z']], how='outer')

    if not path.exists("result"):
        makedirs("result")

    # ploting the result
    END = vsl.visualisation_3d (TRANSMAP3D, EXP_FILE.split("/")[-1].split(".")[0])

    print("duration of the programm : ", complete_time(START, END))

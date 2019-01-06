#!/usr/bin/env python3
"""Program to print a 3D map."""

import argparse
import subprocess
import sys
import time
import pandas
from scipy.spatial.distance import squareform, pdist
import file_manager
import map3d as mp


if __name__ == '__main__':
    # Arguments
    PARSER = argparse.ArgumentParser(
        description="3D transcription map, main file")
    PARSER.add_argument('gene3Dpos', type=str,
                        help='path to the file contains 3D positions of genes')
    PARSER.add_argument('geneTranscriptionTable', type=str,
                        help='path to the file contains expression of genes')
    PARSER.add_argument('NCloserGenes', type=int,
                        help='number of closer genes to take')
    PARSER.add_argument('colorMap', type=str, help="color map for plot.\
        Here the name:\
        https://matplotlib.org/examples/color/colormaps_reference.html\
        advice: use 'magma'.")
    PARSER.add_argument('maxDist', type=float,
                        help="maximal distance autorized for two gene")
    ARGS = PARSER.parse_args()

    # Start of the program
    START = time.time()
    GENE_3D_POS = file_manager.check_file(ARGS.gene3Dpos)
    GENE_EXPRESSION = file_manager.check_file(ARGS.geneTranscriptionTable)
    DF_3DPOS = pandas.read_table(GENE_3D_POS, header=0)
    DF_GENE_EXPR = pandas.read_table(GENE_EXPRESSION, header=0)
    GENE_IN_3DPOS = list(DF_3DPOS.index)
    GENE_IN_EXPR = list(DF_GENE_EXPR.index)
    OVERLAP_GENE = set(GENE_IN_EXPR).intersection(GENE_IN_3DPOS)
    if not OVERLAP_GENE:
        sys.exit("There are no genes in common between two files.")
    # to have the same genes in the two matrix
    DF_3DPOS = DF_3DPOS.loc[list(OVERLAP_GENE)]
    DF_GENE_EXPR = DF_GENE_EXPR.loc[list(OVERLAP_GENE)]
    DF_GENE_EXPR.to_csv("../data/for_corr.csv", "\t")
    # save matrix to use R cor function, rpy2 not stable
    subprocess.check_call(['Rscript', 'correlation.R'], shell=False)
    # called cor R function because it is faster than pandas.corr
    DF_CORR = pandas.read_table("../data/correlation.csv")
    # DF_CORR = pandas.read_pickle("../data/corr.pckl")
    # read the R code result
    DF_DIST = pandas.DataFrame(data=squareform(pdist(DF_3DPOS.drop('chr', 1),
                                                     metric="euclidean")),
                               index=DF_3DPOS.index,
                               columns=DF_3DPOS.index)
    print("{} secondes to create matrix".format(time.time() - START))

    #build and display the 3D map
    START = time.time()
    TRANSMAP = mp.TranscripMap3D(DF_3DPOS, DF_DIST, DF_CORR)
    TRANSMAP.create_plot_dic(ARGS.NCloserGenes, ARGS.maxDist)
    print("{} secondes to create PLOT_DIC".format(time.time() - START))
    TRANSMAP.display(ARGS.colorMap)

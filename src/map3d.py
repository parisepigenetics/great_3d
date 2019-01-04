#!/usr/bin/env python3
"""Program to print a 3D map"""

import argparse
import subprocess
import time
import pandas
# import numpy as np # oour les xs, ys, zs
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import squareform, pdist
import file_manager


if __name__ == '__main__':
    START = time.time()
    PARSER = argparse.ArgumentParser(
        description="3D transcription map, main file")
    PARSER.add_argument('gene3Dpos', type=str,
                        help='path to the file contains 3D positions of genes')
    PARSER.add_argument('geneTranscriptionTable', type=str,
                        help='path to the file contains expression of genes')
    PARSER.add_argument('NCloserGenes', type=int,
                        help='number of closer genes to take')
    PARSER.add_argument('colorMap', type=str, help="color map for plot")
    ARGS = PARSER.parse_args()
    GENE_3D_POS = file_manager.check_file(ARGS.gene3Dpos)
    GENE_EXPRESSION = file_manager.check_file(ARGS.geneTranscriptionTable)
    DF_3DPOS = pandas.read_table(GENE_3D_POS, header=0)
    DF_GENE_EXPR = pandas.read_table(GENE_EXPRESSION, header=0)
    GENE_IN_3DPOS = list(DF_3DPOS.index)
    GENE_IN_EXPR = list(DF_GENE_EXPR.index)
    OVERLAP_GENE = set(GENE_IN_EXPR).intersection(GENE_IN_3DPOS)
    # mettre une erreur si rien ne s'overlap
    # to have the same genes in the two matrix
    DF_3DPOS = DF_3DPOS.loc[list(OVERLAP_GENE)]
    DF_GENE_EXPR = DF_GENE_EXPR.loc[list(OVERLAP_GENE)]
    DF_GENE_EXPR.to_csv("../data/for_corr.csv", "\t")
    # called cor R function because it is faster than pandas.corr
    subprocess.check_call(['Rscript', 'correlation.R'], shell=False)
    # read the R code result
    DF_CORR = pandas.read_table("../data/correlation.csv")
    # DF_CORR = pandas.read_pickle("../data/corr.pckl")
    DF_DIST = pandas.DataFrame(data=squareform(pdist(DF_3DPOS.drop('chr', 1),
                                                     metric="euclidean")),
                               index=DF_3DPOS.index,
                               columns=DF_3DPOS.index)
    print("{} secondes to create matrix".format(time.time() - START))
    # caculate the distance matrix
    DIST_LIST = [pandas.DataFrame]*len(DF_DIST)
    # dictionary with gene name for key and list of N closer genes for values
    DIST_DIC = {}
    # dictionary with gene name for key and list[x,y,z,correlation]
    PLOT_DIC = {}
    for i, val in enumerate(DF_DIST):
        DIST_LIST[i] = DF_DIST.iloc[[i], :]
        tmp_array = DIST_LIST[i].values.argsort()  # sorted array
        tmp_index = [-1]*ARGS.NCloserGenes
        count = ARGS.NCloserGenes-1
        j = 0
        while count >= 0:
            if tmp_array[0][j] not in list(range(i-2, i+3)):  # we don't take gene linearly closer
                tmp_index[count] = tmp_array[0][j]
                count = count - 1
            j += 1
        DIST_DIC[DIST_LIST[i].index[0]] = list(DIST_LIST[i].iloc[:, tmp_index].columns)
        sum_corr = 0
        for closer_gene in DIST_DIC[DIST_LIST[i].index[0]]:
            sum_corr = sum_corr + abs(DF_CORR[DIST_LIST[i].index[0]][closer_gene])
        coord = DF_3DPOS.loc[DIST_LIST[i].index[0]].values[1:]
        PLOT_DIC[DIST_LIST[i].index[0]] = [coord[0], coord[1], coord[2], sum_corr]
    print("{} secondes to create PLOT_DIC".format(time.time() - START))
    FIG = plt.figure()
    AX = FIG.add_subplot(111, projection='3d')
    XS = [-1]*len(PLOT_DIC)
    YS = [-1]*len(PLOT_DIC)
    ZS = [-1]*len(PLOT_DIC)
    COLOR = [-1]*len(PLOT_DIC)
    for i, gene in enumerate(PLOT_DIC):
        XS[i] = PLOT_DIC[gene][0]
        YS[i] = PLOT_DIC[gene][1]
        ZS[i] = PLOT_DIC[gene][2]
        COLOR[i] = PLOT_DIC[gene][3]
    P = AX.scatter(XS, YS, ZS, c=COLOR, cmap=ARGS.colorMap, marker="o")
    AX.set_xlabel('X Label')
    AX.set_ylabel('Y Label')
    AX.set_zlabel('Z Label')
    FIG.colorbar(P)
    print("{} secondes to Plot".format(time.time() - START))
    plt.show()


# Tu as un tableau de coordonnées x, y, z et un tableau d'expression des gènes
# 1/ Il faut calculer la matrice de corrélation (tableau exp des gènes)
# 2/ calcul des distances avec le tableau des coord
# 3/ pour un gène donné, sélection des N gènes les plus proches puis tu sommes
# la correlation des N gènes les plus proches
# (utiliser la matrice de corrélation) cela te donne un chiffre X
# 5/ Plot des coordonnées avec en couleur le chiffre X
# tester le NaN
# errors catcher
# ajouer label des points, titre graphe, meilleures différences dans les couleurs.

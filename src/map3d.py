#!/usr/bin/env python3
"""Program to print a 3D map.
Comments in this code explains the above line"""

import argparse
import subprocess
import sys
import time
import pandas
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.distance import squareform, pdist
import file_manager


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
        https://matplotlib.org/examples/color/colormaps_reference.html")
    PARSER.add_argument('maxDist', type=float, help="maximal distance autorized for two gene")
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

    # caculate the distance matrix
    DIST_LIST = [pandas.DataFrame]*len(DF_DIST)
    # List which contains for one index one row of the DF_DIST
    DIST_DIC = {}
    # dictionary which contains {"gene_id": [N closer gene_id ]}
    PLOT_DIC = {}
    # dictionary with gene name for key and list[x,y,z,correlation]
    for i, val in enumerate(DF_DIST):
        selec_gene = DF_DIST.iloc[[i], :]
        # array contains one selected gene and his distances with other genes.
        sorted_array_index_by_dist = selec_gene.values.argsort()
        # index in DF_DIST of closer genes of selec_gene sorted by distance
        index_of_closer_genes = [-1]*ARGS.NCloserGenes
        count = ARGS.NCloserGenes-1
        j = 0

        while count >= 0:
            if j > len(sorted_array_index_by_dist[0]) - 1:
                sys.exit("There are too much NaN or maxDist is too big.")
            closer_gene = sorted_array_index_by_dist[0][j]
            if closer_gene not in list(range(i-2, i+3)):
                # we don't take gene linearly closer,
                # so if index is -2 or 2 we don't take it
                # because close index means close gene linearly.
                if DF_DIST.iloc[[i],[closer_gene]].values[0][0] < ARGS.maxDist\
                or\
                pandas.isna(DF_CORR.loc[selec_gene.index[0], selec_gene.columns[closer_gene]])\
                ==False:
                    # we verify if value is not NaN and if the gene is not to far
                    index_of_closer_genes[count] = closer_gene
                    # save the index of closer_gene which satisfate conditions
                    count = count - 1
            j += 1

        DIST_DIC[selec_gene.index[0]] = list(selec_gene.iloc[:, index_of_closer_genes].columns)
        # dictionary which contains {"gene_id": [N closer gene_id ]}
        # calculate the transcription score
        sum_corr = 0
        for closer_genes in DIST_DIC[selec_gene.index[0]]:
            sum_corr = sum_corr + abs(DF_CORR[selec_gene.index[0]][closer_genes])
        coord = DF_3DPOS.loc[selec_gene.index[0]].values[1:]
        PLOT_DIC[selec_gene.index[0]] = [coord[0], coord[1], coord[2], sum_corr]

    print("{} secondes to create PLOT_DIC".format(time.time() - START))
    FIG = plt.figure()
    AX = FIG.add_subplot(111, projection='3d')
    XS = [-1]*len(PLOT_DIC)
    YS = [-1]*len(PLOT_DIC)
    ZS = [-1]*len(PLOT_DIC)
    COLOR = [-1]*len(PLOT_DIC)
    LABELS = [""]*len(PLOT_DIC)
    for i, gene in enumerate(PLOT_DIC):
        XS[i] = PLOT_DIC[gene][0]
        YS[i] = PLOT_DIC[gene][1]
        ZS[i] = PLOT_DIC[gene][2]
        COLOR[i] = PLOT_DIC[gene][3]
        LABELS[i] = gene
    P = AX.scatter(XS, YS, ZS, c=COLOR, cmap=ARGS.colorMap, marker="o")
    AX.set_xlabel('X Label')
    AX.set_ylabel('Y Label')
    AX.set_zlabel('Z Label')
    FIG.colorbar(P)
    print("{} secondes to Plot".format(time.time() - START))
    plt.title("3D TRANSCRIPTION MAP", loc="left")
    plt.show()


# Tu as un tableau de coordonnées x, y, z et un tableau d'expression des gènes
# 1/ Il faut calculer la matrice de corrélation (tableau exp des gènes)
# 2/ calcul des distances avec le tableau des coord
# 3/ pour un gène donné, sélection des N gènes les plus proches puis tu sommes
# la correlation des N gènes les plus proches
# (utiliser la matrice de corrélation) cela te donne un chiffre X
# 5/ Plot des coordonnées avec en couleur le chiffre X
# ajouer label des points. diviser le code?? revérifier l'overlap et save image.
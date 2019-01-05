#! /usr/bin/python3
"""
Genes position and expression files Parsing
"""
import argparse
import pandas as pd
import scipy.spatial
import multiprocessing as mp
from functools import partial
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from os import makedirs, path
import ploting as pltg
import numpy as np


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


#def close_genes_correlation (list_data_nbr, gene_name):
    #"""

    #"""
    #closeGenes = list_data_nbr[0].loc[gene_name,].sort_values()[1:list_data_nbr[2]+1].index.tolist()
    #gene3Dcorr = abs(list_data_nbr[1].loc[gene_name,closeGenes].sum())
    #dataframe.loc[row,col]
    #list[0]=dis  list[1]=corr  list[2]=n
    #list_data_nbr[-1][gene_name]=gene3Dcorr
    #return(list_data_nbr[-1])


def close_genes_correlation (dict_matrix, corr_matrix, nbr_gen, overr):
    """
    """
    SUM_CORR = {}
    for gene in overr :
        closeGenes = dict_matrix.loc[gene,].sort_values()[1:nbr_gen+1].index.tolist()
        gene3Dcorr = abs(corr_matrix.loc[gene,closeGenes].sum())
        SUM_CORR[gene]=gene3Dcorr
    return(SUM_CORR)

def visualisation_3d (data_frame, file_name):
    """
    """
    visual = plt.figure().gca(projection='3d')

    visual.scatter(data_frame['X'], data_frame['Y'], data_frame['Z'], c = data_frame['sum_corr'], cmap="jet") #jet blue -> red
    visual.set_xlabel('X')
    visual.set_ylabel('Y')
    visual.set_zlabel('Z')

    plt.savefig("result/"+file_name+"_fig.pdf")
    plt.show()

######################################################  2D

    names = np.array(data_frame.index.tolist())
    x = np.array(data_frame['X'].values)
    #print(" x values : ", x[1])
    y = np.array(data_frame['Y'].values)
    c = np.array(data_frame['sum_corr'].values)

    z = np.array(data_frame['Z'].values)

    norm = plt.Normalize(1,4)
    cmap = plt.cm.jet

    #fig,ax = plt.subplots()
    #sc = plt.scatter(x,y,c=c, s=100, cmap="jet")

    fig = plt.figure(figsize = (16,10))
    ax = fig.add_subplot(111, projection = '3d')
    sc = ax.scatter(x, y, z, c=c, cmap="jet", depthshade = False, picker = True)


    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="-|>"))
    annot.set_visible(False)

    def update_annot(ind):

        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = ""
        if (len(ind["ind"]) > 1):
            for n in ind["ind"] :
                text = text+" "+names[n]
        else :
            text = names[ind["ind"]]

        annot.set_text(text)
        annot.get_bbox_patch().set_facecolor(cmap(c[ind["ind"][0]]))


    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)

    plt.show()

    #pltg.visualize3DData(data_frame.reset_index().values)



if __name__ == "__main__":


    PARSER = argparse.ArgumentParser()
    PARSER.add_argument("pos_file", help="the file containing the x y z coordinates of the genes", type=str)
    PARSER.add_argument("exp_file", help="the file containing the genes expression data", type=str)
    PARSER.add_argument("nbr_genes", help="the number of close genes to select", type=int)
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
    print("dist matrix :")
    print(DIST_MATRIX)

    CORR_MATRIX = correlation_matrix(GENE_EXP.transpose())
    print("corr matrix :")
    print(CORR_MATRIX)


    #dic = {}

    #func = partial(close_genes_correlation, [DIST_MATRIX, CORR_MATRIX, NBR_GEN, dic])
    #dic = POOL.map(func, overr)

    #dic = dic [0]
    #dic2 = pd.DataFrame([dic], columns=dic.keys())

    SUM_CORR = close_genes_correlation (DIST_MATRIX, CORR_MATRIX, NBR_GEN, overr)

    SUM_CORR = pd.DataFrame([SUM_CORR], columns=SUM_CORR.keys())
    SUM_CORR.rename(index = {0: "sum_corr"}, inplace = True)
    print("sum corr matrix :")
    print(SUM_CORR)

    TRANSMAP3D = SUM_CORR.transpose().join(GEN_POS[['X','Y','Z']], how='outer')
    print("transmap 3d matrix :")
    print(TRANSMAP3D)

    if not path.exists("result"):
        makedirs("result")

    visualisation_3d (TRANSMAP3D, EXP_FILE.split("/")[-1].split(".")[0])

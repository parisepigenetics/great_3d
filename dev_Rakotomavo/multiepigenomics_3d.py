#!/usr/bin/python
# -*- coding: utf-8  -*-

import pandas as pd
import numpy as np
import scipy
from scipy.spatial import distance_matrix
from multiprocessing import Pool
import math

def distance_matrice(data):
    '''
    Function that reads a gene position or a gene expression file
    and return the distance Matrice into a DataFrame
    '''
    dt = pd.read_csv(data, sep = '\t')
    ndarray = scipy.spatial.distance.pdist(dt)
    matrix_uni = scipy.spatial.distance.squareform(ndarray)
    matrice_dist = pd.DataFrame(matrix_uni)
    for i in range(len(matrice_dist)) :
        matrice_dist= matrice_dist.rename(columns= {i:dt.index[i]});
        matrice_dist= matrice_dist.rename(index= {i:dt.index[i]})
    return matrice_dist

def dico_matrice(matrice):
    '''
    Function that return All distance of genes
    '''
    dico = {}
    for i in range(len(matrice)):
        #dico[dt.index[i]] =  matrice[dt.index[i]].sort_values()[1:10].index
        dico[matrice.index[i]] =   matrice[matrice.index[i]].sort_values()[1:]
    return dico

def dico_N_matrice(dico ,  nb_genes):
    '''
    Function that return the N closest Genes, Because
    The multiprocessing Function doesn t work on a
    function with more than One argument
    '''
    dico_N_genes = {}
    for i in range(len(dico)):
        a = list(dico.keys())[i]
        b = list(dico.values())[i][1:nb_genes+1]
        dico_N_genes[a] = b
    return dico_N_genes

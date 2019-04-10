#!/usr/bin/python
# -*- coding: utf-8  -*-

import pandas as pd
import numpy as np
import scipy
from scipy.spatial import distance_matrix
import math

def distance_matrice(data):
    import pandas as pd
    dt = pd.read_csv(data, sep = '\t')
    ndarray = scipy.spatial.distance.pdist(dt)
    matrix_uni = scipy.spatial.distance.squareform(ndarray)
    matrice_dist = pd.DataFrame(matrix_uni)
    for i in range(len(matrice_dist)) :
        matrice_dist= matrice_dist.rename(columns= {i:dt.index[i]});
        matrice_dist= matrice_dist.rename(index= {i:dt.index[i]})
    return matrice_dist

def dico_matrice(matrice):
    dico = {}
    for i in range(len(matrice)):
        #dico[dt.index[i]] =  matrice[dt.index[i]].sort_values()[1:10].index
        dico[matrice.index[i]] =   matrice[matrice.index[i]].sort_values()[1:11]
    return dico

'''if __name__ == '__main__':
    #print(dico('1000genes'))
    matrice = gd.dico('test_GE_Dist_Miara.tab')
    from multiprocessing import Pool
    with Pool(5) as p :
        x= p.map(gd.dico_matrice,[matrice])
    print(x)'''

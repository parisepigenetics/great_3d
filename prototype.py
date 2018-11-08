#! /usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np


def calculate_distances(coord):

    for i in range(coord.shape[0]):
        new_dist =  np.sqrt((coord.X.shift()-coord['X'][i+1:])**2 +\
                            (coord.Y.shift()-coord['Y'][i+1:])**2 +\
                            (coord.Z.shift()-coord['Z'][i+1:])**2)
        new_dist.name = coord.index[i]
            # new_dist.sort_values(ascending=False, na_position='first')
        
        try:
            dist
            dist = pd.concat([dist,new_dist], axis=1)

        except:
            dist = new_dist
    return dist

if __name__ == "__main__":

    GENE_COORDINATES_FILE = 'data/toy_ex_coord.txt'
    GENE_EXPRESSION_FILE = 'data/toy_ex_profiles.txt'

    COORD = pd.read_csv(GENE_COORDINATES_FILE, sep='\t')

    ### Main calculations
    #####################
    DIST = calculate_distances(COORD)

    df.sort_values(by='PF3D7_0100100', ascending=False, na_position='first')


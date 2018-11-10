"""
.. module:: parsing
   :synopsis: This module implements all the functions to parse either input
                or additional necessary files.
"""

import numpy as np

def calculate_transcription_map(coordinates, correlations, n, gene):
    coordinates_copy = coordinates.copy()
    coordinates_copy.loc[gene] = np.NaN
    distances =  np.sqrt((coordinates_copy['X'] - coordinates['X'][gene])**2 +\
                         (coordinates_copy['Y'] - coordinates['Y'][gene])**2 +\
                         (coordinates_copy['Z'] - coordinates['Z'][gene])**2)
    distances = distances.sort_values(ascending=True, na_position='last')
    close_genes = distances[:n].index
    transcription_map = sum(abs(correlations[gene][close_genes]))
    return transcription_map

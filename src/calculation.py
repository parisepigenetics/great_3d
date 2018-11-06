"""
.. module:: calculation
   :synopsis: 
"""

# Third-party modules
import numpy as np

def calculate_distances_matrix(genes):
    nb_genes = len(genes)
    # This sets a numpy matrix of shape query * query which will contain all
    # distances between all pairs of genes coordinates
    distances = np.empty((nb_genes, nb_genes), dtype=object)
    # Filling the matrix afterwards with "NaN" is faster
    distances.fill(np.nan)

    for i in range(nb_genes):
        for j in range(i+2, nb_genes-1):
            # Calculate distance between two genes
            distances[i, j] = genes[i].calculate_distance(genes[j])
    return distances

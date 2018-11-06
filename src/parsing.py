"""
.. module:: parsing
   :synopsis: This module implements all the functions to parse either input
                or additional necessary files.
"""

# Third-party modules
import numpy as np
from src.gene import Gene

def parse_gene_coordinates(gene_coordinates_file):
    genes = []
    with open(gene_coordinates_file, "r") as file:
        # The header is skipped
        file.readline()
        for line in file:
            # Elements of the current line are stored into a list
            ele = line[:-1].split("\t")
            # A gene object is added in the genes list
            genes.append(Gene(ele[0], ele[1],
                np.array([float(ele[2]), float(ele[3]), float(ele[4])])))
    return genes

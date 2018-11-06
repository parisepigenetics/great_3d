"""
.. module:: Gene
   :synopsis: This module implements the Gene class.
"""

# Third-party modules
import numpy as np

class Gene:
    """
    .. class:: Gene

      This class groups informations about a gene.

    Attributes:
        name: Name of the gene
        chr: Chromosome corresponding
        coords: X, Y, Z coordinates
    """
    __slots__ = ("name", "chr", "coords")

    def __init__(self, name, chromosome, coords):
        self.name = name
        self.chr = chromosome
        self.coords = coords

    def __str__(self):
        return self.name + " " + self.chr + " " + str(self.coords)

    def calculate_distance(self, gene):
        """
            Calculate Euclidian distance between two residues with THE most efficient method.
            Formula: distance = sqrt((xa-xb)**2 + (ya-yb)**2 + (za-zb)**2)

            Args:
                residue (object): An object of the Residue class.

            Returns:
                dist (float): The calculated distance.
        """
        a_min_b = self.coords - gene.coords
        distance = np.sqrt(np.einsum('i,i->', a_min_b, a_min_b))
        return distance

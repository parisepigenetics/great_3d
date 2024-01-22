"""Module with the main classes/functions for the base calculations of the 3D
analysis of transcriptoms."""

import numpy as np
import pandas as pd
from scipy import stats, spatial
# from multiprocessing import cpu_count, Pool, Process
# from functools import partial
# from scipy.spatial import distance_matrix

import pprint  # for testing only!
from great_3d import timing


@timing
def calculate_distance(df, distance_metric="seuclidean"):
    """Take a genes position dataFrame, return the distance matrix as a pandas data frame.
    - Args:
    - `df`: The genes position dataFrame.
    - `distance_metric`: The numpy distance metric to use (default: "seuclidean").

    - Returns:
    - The distance matrix as a pandas data frame.
    """
    geneNames = df.index
    da = df[["X", "Y", "Z"]].to_numpy()
    ndarray = spatial.distance.pdist(da, distance_metric, V=None)  # Look at other available implementations from numpy
    matrix_uni = spatial.distance.squareform(ndarray)
    matrix_dist = pd.DataFrame(matrix_uni)
    matrix_dist.columns = geneNames
    matrix_dist.index = geneNames
    return matrix_dist


@timing
def sorting_distances(dist_df):
    """Take a genes distance matrix (Pandas DataFrame), return a dictionary of
    gene names sorted by distance with each gene as a key.
    - Args:
    - `dist_df`: The genes distance matrix (Pandas DataFrame).

    - Returns:
    - A dictionary of gene names sorted by distance.
    """
    return {
        gene_name: dist_df.loc[gene_name, :].sort_values()
        for gene_name in dist_df.index
    }


@timing
def sum_correlation(dists_sorted, ge_file, no_genes, correlation_type):
    """Take the dictionnary of the closest genes, the gene expression file, the
    number of genes we need to compute correlation and the correlation type.

    - Args:
    - `dists_sorted`: The dictionary of the closest genes.
    - `ge_file`: The gene expression file.
    - `no_genes`: The number of genes to compute correlation.
    - `correlation_type`: The correlation type.

    - Returns:
    - A complex dictionary containing the sum of correlations, the sum of
    absolute correlations, and the neighboring genes for each gene.
    """
    # Read the GE file
    geDF = pd.read_table(ge_file)
    # Convert the GE data frame to a dictionary
    # geD = geDF.T.to_dict("list")
    # FIXME If a gene has zero expression we give it zero correlation immediately
    correlation_sums = {}
    # Here is the actual calculation
    # if correlation_type == "pearson":
    #     corr_funct = stats.pearsonr
    # elif correlation_type == "spearman":
    #     corr_funct = stats.spearmanr
    # elif correlation_type == "kendall":
    #     corr_funct = stats.kendalltau
    # else:
    #     raise NotImplementedError("Desired correlation function not yet impelmented.")
    # Compute correlation matrix once
    geNumpy = geDF.to_numpy()
    corrMnp = np.corrcoef(geNumpy)
    corrMat = pd.DataFrame(corrMnp, index=geDF.index, columns=geDF.index)
    for gene_ref, sorted_genes in dists_sorted.items():
        # TODO check if we gain time when we paralelise this for loop!
        selected_genes = list(sorted_genes[1: no_genes+1].index)
        # Select all proximal correlations from the coorelation matrix.
        correlation = corrMat.loc[gene_ref, selected_genes]
        # correlation = np.nan_to_num(correlation)  # Correction in case of NAN values
        sum_correlationA = sum(abs(correlation))
        sum_correlation = sum(correlation)
        correlation_sums[gene_ref] = [sum_correlation, sum_correlationA, selected_genes]
    return correlation_sums

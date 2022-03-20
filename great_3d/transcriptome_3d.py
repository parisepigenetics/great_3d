"""Module with the main classes/functions for the 3D analysis of transcriptome.
"""

import numpy as np
import pandas as pd
import scipy
from multiprocessing import cpu_count, Pool, Process
from functools import partial
from scipy.spatial import distance_matrix
from scipy.stats import pearsonr, kendalltau, spearmanr
import plotly
import plotly.graph_objs as go

import time
import pprint  # for testing only!


def timing(f):
    """Wrapper to time functions.py

    Works as a decorator. Taken from https://stackoverflow.com/questions/5478351/python-time-measure-function
    """
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print("{:s} function took {:.3f} ms".format(f.__name__, (time2 - time1) * 1000.0))
        return ret
    return wrap


## Actual functions
@timing
def calculate_distance(position_file):
    # TODO parallelise the calculation of the distance matrix.
    """Take a genes position file, return the distance matrix as a pandas data frame.
    """
    df = pd.read_csv(position_file, sep="\t")
    position_file.seek(0)
    del df["chr"]
    geneNames = df.index
    da = df.to_numpy()
    ndarray = scipy.spatial.distance.pdist(da, 'seuclidean', V=None)
    matrix_uni = scipy.spatial.distance.squareform(ndarray)
    matrix_dist = pd.DataFrame(matrix_uni)
    matrix_dist.columns = geneNames
    matrix_dist.index = geneNames
    return matrix_dist


@timing
def sorting_distances(dist_df):
    """Take a genes distance matrix (Pandas DataFrame), return a dictionary of gene names sorted by distance with each gene as a key.
    """
    sortedDict = {}
    for gene_name in dist_df.index:
        sortedDict[gene_name] = dist_df[gene_name].sort_values()
    return sortedDict


@timing
def sum_correlation(dists_sorted, ge_file, no_genes, correlation_type):
    """Take the dictionnary of the closest genes, the gene expression file, the number of genes we need to compute correlation and the correlation type.

    Return two dictionnaries: One wit the sum of correlations and one with the sum of absolute correlations for each gene
    """
    # Read the GE file
    geDF = pd.read_csv(ge_file, sep="\t")
    # Convert the GE data frame to a dictionary.
    geD = geDF.T.to_dict("list")
    correlation_sums = {}
    # FIXME If a gene has zero expression we give it zero correlation immediately.
    # Here is the actual calculation.
    for gene_ref, sorted_genes in dists_sorted.items():
        # TODO check if we gain time when we paralelise this for loop!
        selected_genes = list(sorted_genes[1 : no_genes].index)
        ref_GE = geD[gene_ref]
        # Select the desired correlation.
        if correlation_type == "pearson":
            #correlationA = [abs(pearsonr(ref_GE, geD[s])[0]) for s in selected_genes]
            correlation = [pearsonr(ref_GE, geD[s])[0] for s in selected_genes]
        elif correlation_type == "kendall":
            #correlationA = [abs(kendalltau(ref_GE, geD[s])[0]) for s in selected_genes]
            correlation = [kendalltau(ref_GE, geD[s])[0] for s in selected_genes]
        elif correlation_type == "spearman":
            #correlationA = [abs(spearmanr(ref_GE, geD[s])[0]) for s in selected_genes]
            correlation = [spearmanr(ref_GE, geD[s])[0] for s in selected_genes]
        # NOTE Each method returns the p-value for a 2 taill-correlation test, so perhpas we can use it in a later analysis
        # TODO the function need to be modular in the sense that we can pass different correlation functions to it. Perhaps use decorators
        correlation = np.nan_to_num(correlation)
        #correlationA = np.nan_to_num(correlationA)  # trick???
        sum_correlationA = sum(abs(correlation))
        sum_correlation = sum(correlation)
        correlation_sums[gene_ref] = [sum_correlation, sum_correlationA, selected_genes]
    return correlation_sums


def generate_genome_3D(genome_coords_file, position_file, correlation_dict):
    """Generate the trace for the 3D genome.
    Genome coordinates file should contain a column named "chr" with the different chromosome names

    Return: A list of all chromosome traces
    """
    traces = []
    colors = plotly.colors.qualitative.Vivid + plotly.colors.qualitative.Safe
    # Build the genome contur trace
    #TODO look at the line 3D plots for a possible alternative https://plotly.com/python/3d-line-plots/
    pos_dt = pd.read_csv(genome_coords_file, sep="\t")
    chroms = pos_dt.loc[:, "chr"].tolist()
    ## Generate the chromosome traces
    chromosomes = list(set(chroms))
    for i in range(len(chromosomes)):
        ch = chromosomes[i]
        # Take the dataframe slice that corresponds to each chromosome
        dfChrom = pos_dt[pos_dt["chr"] == ch]
        traceChr = go.Scatter3d(
            x=dfChrom["X"],
            y=dfChrom["Y"],
            z=dfChrom["Z"],
            name=ch,
            hoverinfo="text",
            mode="lines",
            opacity=0.5,
            line=dict(width=4, color=colors[i]),
            showlegend=True)
        traces.append(traceChr)
    ## Generate the genes traces
    nearGenes = []
    pos_df = pd.read_csv(position_file, sep="\t")
    for gene_ref in correlation_dict:
        pos_df.loc[gene_ref, "Corr"] = correlation_dict[gene_ref][0]
        pos_df.loc[gene_ref, "CorrA"] = correlation_dict[gene_ref][1]
        # Append the list of selectes genes with a string
        nearGenes.append('<br>'.join(correlation_dict[gene_ref][2]))
    corr = pos_df.loc[:, "Corr"].tolist()
    corrA = pos_df.loc[:, "CorrA"].tolist()
    chroms = pos_df.loc[:, "chr"].tolist()
    # Construct the hover text list
    htext = [f"{n}<br>{m}<br>{c:3f}<br>{s}" for n, m, c, s in zip(pos_df.index, chroms, corr, nearGenes)]
    htextA = [f"{n}<br>{m}<br>{c:3f}<br>{s}" for n, m, c, s in zip(pos_df.index, chroms, corrA, nearGenes)]
    ## Construct the orrelation traces
    X = pos_df["X"]
    Y = pos_df["Y"]
    Z = pos_df["Z"]
    ## Correlation trace
    traceCorr = go.Scatter3d(
        x=X,
        y=Y,
        z=Z,
        name="Corr",
        hoverinfo="text",
        hovertext=htext,
        mode="markers",
        opacity=0.75,
        marker=dict(size=3, color=corr, colorscale="RdYlBu_r", showscale=True),
        showlegend=True)
    ## Absolute correlation trace
    traceCorrA = go.Scatter3d(
        x=X,
        y=Y,
        z=Z,
        name="CorrABS",
        visible='legendonly',
        hoverinfo="text",
        hovertext=htextA,
        mode="markers",
        opacity=0.75,
        marker=dict(size=3, color=corrA, colorscale="Hot_r", showscale=True),
        showlegend=True)
    traces.append(traceCorr)
    traces.append(traceCorrA)
    ## Extra traces
    # Significant correlation trace
    #signifGenes = get_significant_corr_genes(correlation_dict)
    # User specified trace
    #userGenes = userGenes
    return traces




def visualise_3D_plotly(traces):
    """Render figure laout for Plotly.
    """
    # Set layout elements, size, margins, legend(s)
    layout = go.Layout(
        plot_bgcolor="#FFF",
        autosize=False,
        width=1600,
        height=1200,
        margin=dict(l=1, r=1, b=1, t=50),
        legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.1),
        modebar={"orientation": "h", "bgcolor": "salmon", "color": "white", "activecolor": "#9ED3CD"},
    )
    # Actual construction of the graph
    fig = go.Figure(data=traces, layout=layout)
    # Remove the axis projections (i.e. spikes), do other things with axes and set the grid colour and lines.
    fig.update_scenes(
        xaxis_spikethickness=1,
        yaxis_spikethickness=1,
        zaxis_spikethickness=1,
        xaxis_spikesides=False,
        yaxis_spikesides=False,
        zaxis_spikesides=False,
        xaxis_spikecolor="#666",
        yaxis_spikecolor="#666",
        zaxis_spikecolor="#666",
        xaxis_title="X",
        yaxis_title="Y",
        zaxis_title="Z",
        xaxis_backgroundcolor="rgb(255, 255, 255)",
        yaxis_backgroundcolor="rgb(255, 255, 255)",
        zaxis_backgroundcolor="rgb(255, 255, 255)",
        xaxis_gridcolor="#ccc",
        yaxis_gridcolor="#ccc",
        zaxis_gridcolor="#ccc",
    )
    return fig



def canvas(with_attribution=True):
    """Placeholder function to show example docstring (NumPy format).

    Replace this function and doc string for your own project.

    Parameters
    ----------
    with_attribution : bool, Optional, default: True
        Set whether or not to display who the quote is from.

    Returns
    -------
    quote : str
        Compiled string including quote and optional attribution.
    """
    quote = "The code is but a canvas to our imagination."
    if with_attribution:
        quote += "\n\t- Adapted from Henry David Thoreau"
    return quote


if __name__ == "__main__":
    # Do something if this file is invoked on its own
    print(canvas())

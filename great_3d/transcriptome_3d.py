"""Module with the main classes/functions for the 3D analysis of transcriptome.
"""

import numpy as np
import pandas as pd
import scipy
from multiprocessing import cpu_count, Pool, Process
from functools import partial
from scipy.spatial import distance_matrix
from scipy import stats
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
        sortedDict[gene_name] = dist_df.loc[gene_name,:].sort_values()
    return sortedDict


@timing
def sum_correlation(dists_sorted, ge_file, no_genes, correlation_type):
    """Take the dictionnary of the closest genes, the gene expression file, the number of genes we need to compute correlation and the correlation type.

    Return a complex dictionnaries containing: The sum of correlations, the sum of absolute correlations and the neighbouring genes for each gene.
    """
    # Read the GE file
    geDF = pd.read_csv(ge_file, sep="\t")
    # Convert the GE data frame to a dictionary
    geD = geDF.T.to_dict("list")
    correlation_sums = {}
    # FIXME If a gene has zero expression we give it zero correlation immediately
    # Here is the actual calculation
    for gene_ref, sorted_genes in dists_sorted.items():
        # TODO check if we gain time when we paralelise this for loop!
        selected_genes = list(sorted_genes[1 : no_genes].index)
        ref_GE = geD[gene_ref]
        # Select the desired correlation
        if correlation_type == "pearson":
            correlation = [stats.pearsonr(ref_GE, geD[s])[0] for s in selected_genes]
        elif correlation_type == "kendall":
            correlation = [stats.kendalltau(ref_GE, geD[s])[0] for s in selected_genes]
        elif correlation_type == "spearman":
            correlation = [stats.spearmanr(ref_GE, geD[s])[0] for s in selected_genes]
        # NOTE Each method returns the p-value for a 2 taill-correlation test, so perhpas we can use it in a later analysis
        # TODO the function need to be modular in the sense that we can pass different correlation functions to it. Perhaps use decorators
        correlation = np.nan_to_num(correlation)
        sum_correlationA = sum(abs(correlation))
        sum_correlation = sum(correlation)
        correlation_sums[gene_ref] = [sum_correlation, sum_correlationA, selected_genes]
    return correlation_sums


def generate_genome_3D(genome_coords_file, position_file, correlation_dict, user_genes):
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
            x=dfChrom.loc[:,"X"],
            y=dfChrom.loc[:,"Y"],
            z=dfChrom.loc[:,"Z"],
            name=ch,
            hoverinfo="text",
            mode="lines",
            opacity=0.5,
            line=dict(width=3, color=colors[i]),
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
    ## Construct the correlation traces
    X = pos_df.loc[:,"X"]
    Y = pos_df.loc[:,"Y"]
    Z = pos_df.loc[:,"Z"]
    ## Correlation trace
    traceCorr = go.Scatter3d(
        x=X,
        y=Y,
        z=Z,
        ids=pos_df.index.values,
        name="Corr",
        hoverinfo="text",
        hovertext=htext,
        mode="markers",
        opacity=0.7,
        marker=dict(size=4, color=corr, colorscale="RdYlBu_r", showscale=True),
        showlegend=True)
    traces.append(traceCorr)
    ## Absolute correlation trace
    traceCorrA = go.Scatter3d(
        x=X,
        y=Y,
        z=Z,
        ids=pos_df.index.values,
        name="CorrABS",
        visible='legendonly',
        hoverinfo="text",
        hovertext=htextA,
        mode="markers",
        opacity=0.7,
        marker=dict(size=4, color=corrA, colorscale="Hot_r", showscale=True),
        showlegend=True)
    traces.append(traceCorrA)
    ## Extra traces
    # Significant correlation trace
    signifGenes = get_significant_corr_genes(correlation_dict)
    posDF_sign = pos_df.loc[signifGenes,:]
    nearGenes = []
    for gene_ref in signifGenes:
        posDF_sign.loc[gene_ref, "Corr"] = correlation_dict[gene_ref][0]
        # Append the list of selectes genes with a string
        nearGenes.append('<br>'.join(correlation_dict[gene_ref][2]))
    corr = posDF_sign.loc[:,"Corr"].tolist()
    chroms = posDF_sign.loc[:,"chr"].tolist()
    # Construct the hover text list
    htext = [f"{n}<br>{m}<br>{c:3f}<br>{s}" for n, m, c, s in zip(posDF_sign.index, chroms, corr, nearGenes)]
    print(posDF_sign.head())
    print(posDF_sign.shape)
    traceSign = go.Scatter3d(
        x=posDF_sign.loc[:,"X"],
        y=posDF_sign.loc[:,"Y"],
        z=posDF_sign.loc[:,"Z"],
        ids=posDF_sign.index.values,
        name="Sign_Corr",
        hoverinfo="text",
        hovertext=htext,
        mode="markers",
        opacity=0.7,
        marker=dict(size=4, color=corr, colorscale="RdYlBu_r", showscale=True),
        showlegend=True)
    traces.append(traceSign)
    # User specified trace
    if user_genes is not None:
        with user_genes as fh:
            userGenes = [l.rstrip() for l in fh]
        posDF_user = pos_df.loc[userGenes,:]
        nearGenes = []
        for gene_ref in userGenes:
            posDF_sign.loc[gene_ref, "Corr"] = correlation_dict[gene_ref][0]
            # Append the list of selectes genes with a string
            nearGenes.append('<br>'.join(correlation_dict[gene_ref][2]))
        corr = posDF_user.loc[:,"Corr"].tolist()
        chroms = posDF_user.loc[:,"chr"].tolist()
        htext = [f"{n}<br>{m}<br>{c:3f}<br>{s}" for n, m, c, s in zip(posDF_user.index, chroms, corr, nearGenes)]
        print(posDF_user.shape)
        traceUser = go.Scatter3d(
            x=posDF_user.loc[:,"X"],
            y=posDF_user.loc[:,"Y"],
            z=posDF_user.loc[:,"Z"],
            name="User_Genes",
            hoverinfo="text",
            hovertext=htext,
            mode="markers",
            opacity=0.7,
            marker=dict(size=4, color=corr, colorscale="RdYlBu_r", showscale=True),
            showlegend=True)
        traces.append(traceUser)
    return traces


def get_significant_corr_genes(corSumsDict, coef=2):
    """Calculate statistical tests and return a list of significantly correlated genes from a sums_of_correlations dictionary
    """
    # Compute the MAD for correlation sums
    corrs = [i[0] for i in list(corSumsDict.values())]
    mad = stats.median_abs_deviation(corrs)
    med = np.median(corrs)
    thresP = med + coef*mad
    thresN = med - coef*mad
    print(med, mad, thresP, thresN)
    signGenes = []
    for k, v in corSumsDict.items():
        if v[0] > thresP or (v[0] <= thresN and v[0] <= 0):  # We need both conditions that's why we cannot use lambda
            signGenes.append(k)
    #print([(k, corSumsDict[k][0]) for k in signGenes])
    return signGenes


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
        hoverlabel=dict(bgcolor='rgba(255,255,255,0.7)', font=dict(color='black')),
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

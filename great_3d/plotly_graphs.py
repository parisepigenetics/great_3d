"""Module for generating and visualizing 3D genome plots using Plotly."""

import plotly
import plotly.graph_objs as go
import pandas as pd
import numpy as np

# import pprint  # for testing only!
from great_3d import timing
from scipy import stats


def get_significant_corr_genes(corSumsDict, coef=2):
    """Calculate statistical tests and return a list of significantly
    correlated genes from a sums_of_correlations dictionary.
    - Args:
    - `corSumsDict`: The dictionary of sums_of_correlations.
    - `coef`: Optional. The coefficient for determining significance (default: 2).

    - Returns:
    - A list of significantly correlated genes.
    """
    # Compute the MAD for correlation sums
    corrs = [i[0] for i in list(corSumsDict.values())]
    mad = stats.median_abs_deviation(corrs)
    med = np.median(corrs)
    thresP = med + coef*mad
    thresN = med - coef*mad
    return [
        k
        for k, v in corSumsDict.items()
        if v[0] > thresP or (v[0] <= thresN and v[0] <= 0)
    ]


@timing
def generate_genome_3D(genome_coords_file, position_df, correlation_dict, user_genes):
    """Generate the trace for the 3D genome.

    Genome coordinates file should contain a column named "chr" with the different chromosome names

    - Args:
    - `genome_coords_file`: The file path to the genome coordinates file.
    - `position_df`: The gene positions DataFrame.
    - `correlation_dict`: The dictionary of gene/genes correlation values.
    - `user_genes`: The user-specified genes.

    - Returns:
    - A list of all plotted traces.
    """
    traces = []
    colours = plotly.colors.qualitative.Pastel + plotly.colors.qualitative.Safe + plotly.colors.qualitative.Vivid
    # Build the genome contur trace
    # TODO look at the line 3D plots for a possible alternative https://plotly.com/python/3d-line-plots/
    pos_dt = pd.read_table(genome_coords_file)
    chroms = pos_dt.loc[:, "chr"].tolist()
    # Generate the chromosome traces
    chromosomes = list(set(chroms))
    for i in range(len(chromosomes)):
        ch = chromosomes[i]
        # Take the dataframe slice that corresponds to each chromosome
        dfChrom = pos_dt[pos_dt["chr"] == ch]
        traceChr = go.Scatter3d(
            x=dfChrom.loc[:, "X"],
            y=dfChrom.loc[:, "Y"],
            z=dfChrom.loc[:, "Z"],
            name=ch,
            hoverinfo="text",
            mode="lines",
            opacity=0.5,
            line=dict(width=2, color=colours[-(i+1)]),  # index from the end of colours
            showlegend=True)
        traces.append(traceChr)
    # Generate the genes traces
    # FIXME find a way to deal with gene ovelarps (perhaps introduce a jitter as a quick fix.)
    nearGenes = []
    pos_df = position_df
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
    # Construct the correlation traces
    X = pos_df.loc[:, "X"]
    Y = pos_df.loc[:, "Y"]
    Z = pos_df.loc[:, "Z"]
    # Correlation trace
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
    # Absolute correlation trace
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
        opacity=0.5,
        marker=dict(size=4, color=corrA, colorscale="Hot_r", showscale=True),
        showlegend=True)
    traces.append(traceCorrA)
    # Extra traces
    # Significant correlation trace
    signifGenes = get_significant_corr_genes(correlation_dict)
    # print("Significant Genes:\n")
    # pp = pprint.PrettyPrinter(indent=2)
    # pp.pprint(signifGenes)  #FIXME remove these two comment printing lines
    posDF_sign = pos_df.loc[signifGenes, :]
    nearGenes = []
    for gene_ref in signifGenes:
        posDF_sign.loc[gene_ref, "Corr"] = correlation_dict[gene_ref][0]
        # Append the list of selectes genes with a string
        nearGenes.append('<br>'.join(correlation_dict[gene_ref][2]))
    corr = posDF_sign.loc[:, "Corr"].tolist()
    chroms = posDF_sign.loc[:, "chr"].tolist()
    # Construct the hover text list
    htext = [f"{n}<br>{m}<br>{c:3f}<br>{s}" for n, m, c, s in zip(posDF_sign.index, chroms, corr, nearGenes)]
    traceSign = go.Scatter3d(
        x=posDF_sign.loc[:, "X"],
        y=posDF_sign.loc[:, "Y"],
        z=posDF_sign.loc[:, "Z"],
        ids=posDF_sign.index.values,
        name="Sign_Corr",
        hoverinfo="text",
        hovertext=htext,
        mode="markers",
        opacity=0.5,
        marker=dict(size=7, color=corr, colorscale="RdYlBu_r", showscale=True),
        showlegend=True)
    traces.append(traceSign)
    # User specified trace
    if user_genes is not None:
        with user_genes as fh:
            userGenes = [l.rstrip() for l in fh]
        # print(f"Overlap Significant AND user:\n")
        # interSectGenes = list(set(signifGenes) & set(userGenes))  #FIXME use this to the dash app
        # pp.pprint(interSectGenes)
        posDF_user = pos_df.loc[userGenes, :]
        nearGenes = []
        for gene_ref in userGenes:
            posDF_sign.loc[gene_ref, "Corr"] = correlation_dict[gene_ref][0]
            # Append the list of selectes genes with a string
            nearGenes.append('<br>'.join(correlation_dict[gene_ref][2]))
        corr = posDF_user.loc[:, "Corr"].tolist()
        chroms = posDF_user.loc[:, "chr"].tolist()
        htext = [f"{n}<br>{m}<br>{c:3f}<br>{s}" for n, m, c, s in zip(posDF_user.index, chroms, corr, nearGenes)]
        traceUser = go.Scatter3d(
            x=posDF_user.loc[:, "X"],
            y=posDF_user.loc[:, "Y"],
            z=posDF_user.loc[:, "Z"],
            name="User_Genes",
            hoverinfo="text",
            hovertext=htext,
            mode="markers",
            opacity=0.8,
            marker=dict(size=4, color=corr, colorscale="RdYlBu_r", showscale=True),
            showlegend=True)
        traces.append(traceUser)
    return traces


def visualise_3D_plotly(traces):
    """Render figure laout for Plotly.

    - Args:
    - `traces`: The list of traces.

    - Returns:
    - The Plotly figure.
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

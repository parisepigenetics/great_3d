"""Module with the multi-omics epigenetic analysis classes/functions.
"""

import math
import time
import numpy as np
import scipy
from multiprocessing import cpu_count, Pool , Process
from functools import partial
import pandas as pd
from scipy.spatial import distance_matrix
from scipy.stats import pearsonr , kendalltau , spearmanr

import plotly
import plotly.graph_objs as go

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import pprint  # for testing only!

## Nice wrapper to time functions. Works as a decorator.
# Taken from https://stackoverflow.com/questions/5478351/python-time-measure-function
def timing(f):
    def wrap(*args):
        time1 = time.time()
        ret = f(*args)
        time2 = time.time()
        print('{:s} function took {:.3f} ms'.format(f.__name__, (time2-time1)*1000.0))
        return ret
    return wrap


def calculate_distance(position_file):
    # TODO parallelise the calculation of the distance matrix.
    """Get a position file, return the distance matrix as a Pandas DataFrame.
    """
    df = pd.read_csv(position_file, sep='\t')
    position_file.seek(0)

    del df[' chr']
    geneNames = df.index
    ndarray = scipy.spatial.distance.pdist(df)
    matrix_uni = scipy.spatial.distance.squareform(ndarray)
    matrix_dist = pd.DataFrame(matrix_uni)
    matrix_dist.columns = geneNames
    matrix_dist.index = geneNames
    return matrix_dist


def sorting_distances(dist_df):
    '''Take a distance matrix (Pandas DataFrame) and return a dictionary of sorted neighboring genes with each gene as a key.
    '''
    sortedDict = {}
    for gene_name in dist_df.index:
        sortedDict[gene_name] = dist_df.loc[gene_name,].sort_values()
    return sortedDict


def sum_correlation(dists_sorted, ge_file, no_genes , correlation_type):
    """Take the dictionnary of the closest genes, the gene expression file and the number of genes we need to compute correlation.
    Return a dictionnary of the sum of correlations for each gene
    """
    # Read the GE file
    geDF = pd.read_csv(ge_file, sep='\t')
    # Convert the GE data frame to a dictionary.
    geD = geDF.T.to_dict('list')
    correlation_sums = {}
    #FIXME If a gene has zero expression we give it zero correlation immediately.
    # Here is the actual calculation.
    for gene_ref, closest_genes in dists_sorted.items():
        # TODO check if we gain time when we paralelise this for loop!
        selected_genes = list(closest_genes[1:no_genes + 1].index)
        ref_GE = geD[gene_ref]
        if(correlation_type == 'pearson'):
            correlations = [abs(pearsonr(ref_GE, geD[s])[0]) for s in selected_genes]
        elif(correlation_type == 'kendall'):
            correlations = [abs(kendalltau(ref_GE, geD[s])[0]) for s in selected_genes]
        elif(correlation_type == 'spearman'):
            correlations = [abs(spearmanr(ref_GE, geD[s])[0]) for s in selected_genes]
        # NOTE in case we need to use the absolute values we can flip to that.  correlations = [abs(pearsonr(ref_GE, geD[s])[0]) for s in selected_genes]
        # NOTE Each method returns the p-value for a 2 taill-correlation test, so perhpas we can use it in a later analysis.
        # TODO the function need to be modular in the sense that we can pass different correlation functions to it. Perhaps use decorators.
        correlations = np.nan_to_num(correlations)
        sum_correlation = sum(correlations)
        correlation_sums[gene_ref] = sum_correlation
    return correlation_sums


def visualization_3D_plotly(position_file, correlation_dict, outfile_name, autoopen = True):
    """Take correlation dictionnary, gene postion file.
    Return HTML outfile with an interactive plot of the 3D gene correlation
    """
    pos_dt = pd.read_csv(position_file, sep='\t')
    pos_dt['Corr'] = ""
    # Adding the correlation column into the DataFrame
    for gene_ref in correlation_dict:
        pos_dt.loc[gene_ref, 'Corr'] = correlation_dict[gene_ref]
    # Get the coordinates and the ploting vaiable(s)
    x = pos_dt.loc[:,"X"].tolist()
    y = pos_dt.loc[:,"Y"].tolist()
    z = pos_dt.loc[:,"Z"].tolist()
    corr = pos_dt.loc[:,"Corr"].tolist()
    corr2 = [x**2 for x in corr]
    # Construct the hover text list
    htext = [f"{n}<br>{c:3f}" for n, c in zip(pos_dt.index,  corr)]
    htext2 = [f"{n}<br>{c:3f}" for n, c in zip(pos_dt.index,  corr2)]
    # Construct the trace(s)
    trace1 = go.Scatter3d(x=x, y=y, z=z, name = "Correlation",
    hoverinfo="text", hovertext = htext, mode='markers',
    marker=dict(size=3, color=corr, colorscale='Portland', cmin =1, cmax = max(corr), opacity=0.67, showscale=True),
    showlegend=True)
    trace2 = go.Scatter3d(x=x, y=y, z=z, name = "Correlation2",
    hoverinfo="text", hovertext = htext2, mode='markers',
    marker=dict(size=3, color=corr2, colorscale='Portland', cmin =2, cmax = max(corr2), opacity=0.67, showscale=True),
    showlegend=True)
    # Combine multiple tracks
    data = [trace1, trace2]
    # Set layout elements, size, margins, legend(s)
    layout = go.Layout(plot_bgcolor="#FFF", autosize=False, width=1600, height=1200,
    margin=dict(l=1, r=1, b=1, t=50),  # leave a small margin to the right for the modebar
    legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.1),
    title_text="GREAT 3D transcriptome", title_x=0.5,
    modebar={'orientation':'h', 'bgcolor':'salmon', 'color':'white', 'activecolor': '#9ED3CD'})
    # Actual construction of the graph
    fig = go.Figure(data=data, layout=layout)
    # Remove the axis projections (i.e. spikes) and do other things with axes
    fig.update_scenes(xaxis_spikethickness=1, yaxis_spikethickness=1, zaxis_spikethickness=1,
    xaxis_spikesides=False, yaxis_spikesides=False, zaxis_spikesides=False,
    xaxis_spikecolor="#666", yaxis_spikecolor="#666", zaxis_spikecolor="#666",
    xaxis_title="X", yaxis_title="Y", zaxis_title="Z",
    xaxis_backgroundcolor="rgb(255, 255, 255)", yaxis_backgroundcolor="rgb(255, 255, 255)", zaxis_backgroundcolor="rgb(255, 255, 255)",
    xaxis_gridcolor="#ccc", yaxis_gridcolor="#ccc", zaxis_gridcolor="#ccc")
    # Set the grid colour and lines
    #fig.update_xaxes(showline=True, linewidth=2, linecolor='black', gridcolor='white')
    return plotly.offline.plot(fig, filename = outfile_name, auto_open=autoopen, config={'displaylogo':False})




def visualization_3D_plotly_line(position_file,correlation_dict, outfile_name):
    """Gets the dictionnary of the closest genes, the gene postion file and
    return an HTML outfile with a plot of the 3D gene correlation by using Plotly
    """
    # FIXME the programm isn't running when we use the file with the chromosome,
    # The problem appears on the pdist ????
    pos_dt = pd.read_csv(position_file , sep ='\t')
    pos_dt['corr'] = ""
    # Adding the correlation columns into the genePos data Frame
    for gene_ref in correlation_dict:
        pos_dt['corr'][gene_ref] = correlation_dict[gene_ref]
    # 3D VISUALIZATION_3D
    x = pos_dt.iloc[:,1].tolist()
    y = pos_dt.iloc[:,2].tolist()
    z = pos_dt.iloc[:,3].tolist()
    corr = pos_dt.iloc[:,4].tolist()
    trace1 = go.Scatter3d(x=x, y=y, z=z, text = pos_dt.index,
    hoverinfo = 'text', mode='markers',
    marker=dict(size=4, color=corr, colorscale='Reds',  # choose a colorscale
    opacity=0.5,  # Transparency
    showscale = False),
    showlegend = True)
    #LINE PLOT
    #Putting a none inside Y columns when the chromosoe changes
    #TODO find the number of chromosomes, create the data list before and append a trace every time you have a new chromosome (you create the trace for each chromosoem in a loop here.) Do something similar for the colour.
    #FIXME fix the line BUG for multi-chromosomes NOW!
    gap = []
    for i in range(len(pos_dt.index)-1) :
        if((pos_dt[' chr'][i]) != (pos_dt[' chr'][i+1])):
            #print('{}lol'.format(i+3))
            chr_change = i+1
            gap.append(chr_change)
    for gene_index in gap :
        pos_dt['Z'][pos_dt.index[gene_index]] = None
        pos_dt['Y'][pos_dt.index[gene_index]] = None
        pos_dt['X'][pos_dt.index[gene_index]] = None
    x2 = pos_dt.iloc[:,1].tolist()
    y2 = pos_dt.iloc[:,2].tolist()
    z2 = pos_dt.iloc[:,3].tolist()
    corr = pos_dt.iloc[:,4].tolist()
    trace2 = go.Scatter3d(x=x2, y=y2, z=z2, text = pos_dt.index,
    hoverinfo = 'text', mode='lines',
    line = dict(color = 'steelblue', width=5),
    showlegend = False , connectgaps = True, opacity=0.5)
    #TODO (Check the line.shape of plottly and it seems it does not allow. PEHAPS it is possible to interpolate and create a curve instead of segments. check interpolate.interp2d(x, y, z, kind='cubic') from scipy.interpolate.)
    data = [trace1, trace2]
    layout = go.Layout(margin=dict(l=0, r=0, b=0, t=0), legend=dict(yanchor="top", y=0.99, xanchor="left",
    x=0.01))
    fig = go.Figure(data=data, layout=layout)
    return plotly.offline.plot(fig, filename = outfile_name, auto_open=False)


def visualization_3D_mtp(position_file,correlation_dict):
    """Gets the dictionnary of the closest genes, the gene postion file and
    return a plot of the 3D gene correlation by using Mathplotlib
    """
    pos_dt = pd.read_csv(position_file, sep='\t')
    fig = plt.figure()
    ax = fig.add_subplot(111 , projection='3d')
    pos_dt['corr'] = ""
    # Adding the correlation columns into the genePos data Frame
    for gene_ref in correlation_dict:
        pos_dt['corr'][gene_ref] = correlation_dict[gene_ref]
    # Collect the values.
    x = pos_dt.iloc[:,1].tolist()
    y = pos_dt.iloc[:,2].tolist()
    z = pos_dt.iloc[:,3].tolist()
    corr = pos_dt.iloc[:,4].tolist()
    # 3D VISUALIZATION
    ax.scatter(x, y, z, c='r', cmap = corr, marker = 'o')
    ax.set_xlabel('X'), ax.set_ylabel('Y'), ax.set_zlabel('Z')
    plt.show()



###############################################################################
#=================================== OBSOLETE functions ======================#
def sumCor_mp(sorted_dists, ge_file, no_genes , type_correlation) :
    cpus = cpu_count() - 1
    dict = {}
    with Pool(cpus) as pool:
        dict = pool.map(sum_correlation, (sorted_dists , ge_file , no_genes , type_correlation))
        pool.close()
        pool.join()
    return dict

def sorting_dists_parallel(dist_df) :
    '''Get a distance matrix (pandas data frame), return a dictionary of genes and their sorted neighbours.
    '''
    cpus = cpu_count() - 1
    sorted_dict = {}
    with Pool(cpus) as pool:
        # Split the data frame
        df_split = np.array_split(dist_df, cpus)
        sorted_dfs = pool.map(sorting_distances, df_split)
        pool.close()
        pool.join()
    # Put together all the results.
    for d in sorted_dfs: sorted_dict.update(d)
    return sorted_dict

def closest_gene_name(dico_closest_gene):
    '''Get the dictionnary of sorted genes and return only the names of the N closest gene with each gene as a key
    '''
    dico_name_closest_gene = {}
    for gene_name in dico_closest_gene :
        closest_genes = list(pd.DataFrame(dico_closest_gene[gene_name][1:]).index) # This computes correlations between the gene of reference and ALL the other genes.
        dico_name_closest_gene[gene_name] = closest_genes
    return dico_name_closest_gene

def correlation_matrix(gene_expression_file):
    """Gets the gene Expression filename, return the Correlation matrix as a Pandas DataFrame.
    """
    gene_expression_df = pd.read_csv(gene_expression_file , sep ='\t').transpose()
    return gene_expression_df.corr(method = 'spearman') #or Pearson

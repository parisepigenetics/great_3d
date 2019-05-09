# -*- coding: utf-8  -*-

"""Module with all the multi-omics epigenetic analysis functions.
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
import plotly.plotly as py
import plotly.graph_objs as go
#plotly.tools.set_credentials_file(username='miara1502', api_key='LM3BdwFIOFpJmq3DP6CQ')
# NOTE: we have to create an account to run the programm correctyl
# NOTE: we can visualize the Plot on our own account

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import pprint

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
    geneNames = df.index
    ndarray = scipy.spatial.distance.pdist(df)
    matrix_uni = scipy.spatial.distance.squareform(ndarray)
    matrix_dist = pd.DataFrame(matrix_uni)
    matrix_dist.columns = geneNames
    matrix_dist.index = geneNames
    return matrix_dist


def sorting_distances(dist_df):
    '''Get a distance matrix (Pandas DataFrame) and return a dictionary of sorted genes with each gene as a key.
    '''
    sortedDict = {}
    for gene_name in dist_df.index:
        sortedDict[gene_name] = dist_df.loc[gene_name,].sort_values()
    return sortedDict


def sum_correlation(sorted_dists, ge_file, no_genes , type_correlation):
    """Gets the dictionnary of the closest genes, the gene expression file and the number of genes we need to compute correlation.
    Return a dictionnary of the sum of correlations for each gene
    """
    # Read the GE file
    geDF = pd.read_csv(ge_file, sep='\t')
    # Convert the GE data frame to a dictionary.
    geD = geDF.T.to_dict('list')
    correlation_sums = {}
    #FIXME If a gene has zero expression we give it zero correlation immediately.
    # Here is the actual calculation.
    for gene_ref, closest_genes in sorted_dists.items():
        # TODO check if we gain time when we paralelise this for loop!
        selected_genes = list(closest_genes[1:no_genes + 1].index)
        ref_GE = geD[gene_ref]
        if(type_correlation == 'pearson'):
            correlations = [abs(pearsonr(ref_GE, geD[s])[0]) for s in selected_genes]
        elif(type_correlation == 'kendall'):
            correlations = [abs(kendalltau(ref_GE, geD[s])[0]) for s in selected_genes]
        elif(type_correlation == 'spearman'):
            correlations = [abs(spearmanr(ref_GE, geD[s])[0]) for s in selected_genes]

        # NOTE incase we need to use the absolute values we can flip to that.  correlations = [abs(pearsonr(ref_GE, geD[s])[0]) for s in selected_genes]
        # NOTE The pearsonr method returns the p-value for a 2 tauiled-correlation test, so perhpas we can use it in a later analysis.
        # TODO the function need to be modular in the sense that we can pass different correlation functions to it. Perhaps use decorators.
        sum_correlation = sum(correlations)
        correlation_sums[gene_ref] = sum_correlation
    return correlation_sums


def visualization_3D_plotly(position_file, correlation_dict):
    """Gets the dictionnary of the closest genes, the gene postion file and
    return a plot of the 3D gene correlation by using Plotly
    """
    pos_dt = pd.read_csv(position_file, sep='\t')
    pos_dt['corr'] = ""
    # Adding the correlation columns into the genePos data Frame
    for gene_ref in correlation_dict:
        pos_dt['corr'][gene_ref] = correlation_dict[gene_ref]
    # Taken from https://plot.ly/python/3d-scatter-plots/ and adjusted with our values
    # 3D VISUALIZATION_3D
    x = pos_dt.iloc[:,0].tolist()
    y = pos_dt.iloc[:,1].tolist()
    z = pos_dt.iloc[:,2].tolist()
    corr = pos_dt.iloc[:,3].tolist()
    trace1 = go.Scatter3d(
        x=x,
        y=y,
        z=z,
        text = pos_dt.index ,
        hoverinfo = 'text' ,

        mode='markers',
        marker=dict(
            size=6,
            color=corr,
            colorscale='Reds', # choose a colorscale
            opacity=0.8
        ),
        showlegend = False
    )

    data = [trace1]
    layout = go.Layout(
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=0
        )
    )
    fig = go.Figure(data=data, layout=layout)
    py.iplot(fig, filename='3d-scatter-colorscale_3D_TRANSMAP')


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
    x = pos_dt.iloc[:,0].tolist()
    y = pos_dt.iloc[:,1].tolist()
    z = pos_dt.iloc[:,2].tolist()
    corr = pos_dt.iloc[:,3].tolist()
    # 3D VISUALIZATION
    ax.scatter(x , y , z , c='r' , cmap = corr , marker = 'o')
    ax.set_xlabel('X label') , ax.set_ylabel('Y label') , ax.set_zlabel('Z label')
    plt.show()




#================ OBSOLETE functions ======================

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

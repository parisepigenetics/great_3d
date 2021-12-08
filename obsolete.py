###############################################################################
#============ OBSOLETE functions from multiEpigenomics3D =====================#
###############################################################################

# OBSOLETE ALL the job is done by the visualization_3D_plotly now
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
    # find the number of chromosomes, create the data list before and append a trace every time you have a new chromosome (you create the trace for each chromosoem in a loop here.) Do something similar for the colour.
    # fix the line BUG for multi-chromosomes NOW!
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

# OBSOLETE All visualisation with plotly now.
# Kept as REFERENCE for matplotlib functionalites
def visualization_3D_mtp(position_file,correlation_dict):
    """Gets the dictionnary of the closest genes, the gene postion file and
    return a plot of the 3D gene correlation by using Mathplotlib
    """
    import matplotlib.pyplot as plt
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

# OLD FUNCTIONS for calculations
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

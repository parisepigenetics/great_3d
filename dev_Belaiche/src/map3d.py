#!/usr/bin/env python3
"""Module to manage 3D Transcription map. """

import sys
import pandas
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class TranscripMap3D:
    """Class to manage 3D Transcription map. """

    def __init__(self, pos3DMatrix, DistMatrix, CorrMatrix):
        """To init the class with a matrix:
        - of 3D position of each gene,
        - of Distances,
        - of correlation of transcription."""

        self.df_pos = pos3DMatrix
        self.df_dist = DistMatrix
        self.df_corr = CorrMatrix
        self.plot_dic = {}
        self.number_of_save = 0

    def create_plot_dic(self, n_closer_genes, max_dist):
        """Function to create a dictionary which is used to plot the genome.
        The dictionary is built like this:
        dictionary with gene name for key and list[x,y,z,correlation]

        Arguments are given by the main.py and the command line."""

        # List which contains for one index one row of the df_dist
        dist_dic = {}
        # dictionary which contains {"gene_id": [N closer gene_id ]}
        self.plot_dic = {}
        # dictionary with gene name for key and list[x,y,z,correlation]
        for i, val in enumerate(self.df_dist):
            selec_gene = self.df_dist.iloc[[i], :]
            # array contains one selected gene and
            # his distances with other genes.
            sorted_array_index_by_dist = selec_gene.values.argsort()
            # index in df_dist of closer genes of selec_gene sorted by distance
            index_of_closer_genes = [-1]*n_closer_genes
            count = n_closer_genes-1
            j = 0

            while count >= 0:
                if j > len(sorted_array_index_by_dist[0]) - 1:
                    sys.exit("There are too much NaN\
                    or max_dist is too small.")
                closer_gene = sorted_array_index_by_dist[0][j]
                if closer_gene not in list(range(i-2, i+3)):
                    # we don't take gene linearly closer,
                    # so if index is -2 or 2 we don't take it
                    # because close index means close gene linearly.
                    if self.df_dist.iloc[[i], [closer_gene]].values[0][0]\
                    < max_dist\
                        or\
                        pandas.isna(self.df_corr.loc[selec_gene.index[0],
                                    selec_gene.columns[closer_gene]])\
                            is False:
                        # we verify if value is not NaN and
                        # if the gene is not to far
                        index_of_closer_genes[count] = closer_gene
                        # save the index of closer_gene
                        # which satisfate conditions
                        count = count - 1
                j += 1

            dist_dic[selec_gene.index[0]] =\
            list(selec_gene.iloc[:, index_of_closer_genes].columns)
            # dictionary which contains {"gene_id": [N closer gene_id ]}
            # calculate the transcription score
            sum_corr = 0
            for closer_genes in dist_dic[selec_gene.index[0]]:
                sum_corr = sum_corr +\
                abs(self.df_corr[selec_gene.index[0]][closer_genes])
            coord = self.df_pos.loc[selec_gene.index[0]].values[1:]
            self.plot_dic[selec_gene.index[0]] =\
            [coord[0], coord[1], coord[2], sum_corr]

    def display(self, color_map):
        """Function to display the 3D map with matplotlib."""

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        xs = [-1]*len(self.plot_dic)
        ys = [-1]*len(self.plot_dic)
        zs = [-1]*len(self.plot_dic)
        color = [-1]*len(self.plot_dic)
        size = [-1]*len(self.plot_dic)
        names = [""]*len(self.plot_dic)
        for i, gene in enumerate(self.plot_dic):
            xs[i] = self.plot_dic[gene][0]
            ys[i] = self.plot_dic[gene][1]
            zs[i] = self.plot_dic[gene][2]
            color[i] = self.plot_dic[gene][3]
            names[i] = gene
        mini=min(color)
        maxi=max(color)
        colors = [(float(i)-mini)/(maxi-mini) for i in color]
        plot = ax.scatter(xs, ys, zs, c=colors, cmap=color_map, marker="o")
        ax.set_xlabel('X axe')
        ax.set_ylabel('Y axe')
        ax.set_zlabel('Z axe')
        fig.colorbar(plot)
        genome_name = names[0].split("_")[0]
        plt.title("3D TRANSCRIPTION MAP OF {}".format(genome_name), loc="left")
        annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="->"))
        annot.set_visible(False)
        # https://stackoverflow.com/questions/7908636/possible-to-make-labels-\
        #appear-when-hovering-over-a-point-in-matplotlib?fbclid=IwAR0XLYQ8h4wJa\
        #SBPzA5OGMmJLLcbZjdtQrC8KtVUX6FSTHsIxVW-HQF2zbA
        
        def update_annot(ind):
            pos = plot.get_offsets()[ind["ind"][0]]
            annot.xy = pos
            text = "{}".format(" ".join([names[n] for n in ind["ind"]]))
            annot.set_text(text)
            annot.get_bbox_patch().set_alpha(0.4)
        
        def onclick(event):
            vis = annot.get_visible()
            if event.inaxes == ax:
                cont, ind = plot.contains(event)
                if cont:
                    update_annot(ind)
                    annot.set_visible(True)
                    fig.canvas.draw_idle()
                else:
                    if vis:
                        annot.set_visible(False)
                        fig.canvas.draw_idle()
        
        def save(event):
            global number_of_save
            if event.inaxes == ax:
                if event.key == "ctrl+alt+s":
                    self.number_of_save+=1
                    fig.savefig("../results/3DMAP_{}_{}.pdf".format(genome_name, self.number_of_save), bbox_inches='tight')
                    print("Figure saved.")

        fig.canvas.mpl_connect("button_press_event", onclick)
        fig.canvas.mpl_connect("key_press_event", save)
        plt.show()


if __name__ == '__main__':
    import map3d
    print(help(map3d))

#! /usr/bin/python3
"""
script that plots the genes 3D position coloured by their expression score
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import time


def visualisation_3d (data_frame, file_name):
    """
    the function that plots
    """
    names = np.array(data_frame.index.tolist())
    x = np.array(data_frame['X'].values)
    y = np.array(data_frame['Y'].values)
    z = np.array(data_frame['Z'].values)
    c = np.array(data_frame['sum_corr'].values)

    # the figure window
    fig = plt.figure(figsize = (16,10))
    # adding axes to the figure
    ax = fig.add_subplot(111, projection = '3d')
    # ploting the genes
    sc = ax.scatter(x, y, z, c=c, cmap="jet", depthshade = False, picker = True)

    # title and labels
    ax.text2D(0.05, 0.95, "3D Transmap of "+file_name, transform=ax.transAxes)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    # the annotation (appearing box and arrow on hover)
    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="-|>"))
    # setting it invisible (untill the hover)
    annot.set_visible(False)

    # the cmap color bar
    fig.colorbar(sc, shrink=0.5, aspect=5)

    # text of the annotation
    def update_annot(ind):

        # setting the position (the indice)
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        # choosing the text of the annotation
        text = ""
        if (len(ind["ind"]) > 1):
            for n in ind["ind"] :
                text = text+" "+"['"+names[n]+"']"
        else :
            text = names[ind["ind"]]

        annot.set_text(text)

    # hovering event management
    def hover(event):
        vis = annot.get_visible()
        # if the mouse on the figure
        if event.inaxes == ax:
            # cont (true/false), ind the indice (the gene indice)
            cont, ind = sc.contains(event)
            # cont true -> mouse hovering a point of sc (a gene)
            if cont:
                # call update_annot with the gene indice
                update_annot(ind)
                # set the annotation visible
                annot.set_visible(True)
                # redraw the figure
                fig.canvas.draw_idle()
            # cont false -> mouse not hovering a gene
            else:
                # if it came out of hovering a gene
                if vis:
                    # set the annotation invisible
                    annot.set_visible(False)
                    # redraw the figure
                    fig.canvas.draw_idle()

    # connect the callback function the the event manager in order to receive events
    fig.canvas.mpl_connect("motion_notify_event", hover)
    
    end = time.time()
    
    # save the plot as pdf
    plt.savefig("result/"+file_name+"_fig.pdf")
    # show the plot
    plt.show()

    return(end)


if __name__ == "__main__":

    print("this programm doesn't run independently sorry !")

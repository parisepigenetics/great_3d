#! /usr/bin/python3
"""
programm that creates distance and correlation matrix
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import time


def visualisation_3d (data_frame, file_name):
    """
    """
    names = np.array(data_frame.index.tolist())
    x = np.array(data_frame['X'].values)
    y = np.array(data_frame['Y'].values)
    z = np.array(data_frame['Z'].values)
    c = np.array(data_frame['sum_corr'].values)

    norm = plt.Normalize(1,4)
    cmap = plt.cm.jet

    fig = plt.figure(figsize = (16,10))
    ax = fig.add_subplot(111, projection = '3d')
    sc = ax.scatter(x, y, z, c=c, cmap="jet", depthshade = False, picker = True)

    ax.text2D(0.05, 0.95, "3D Transmap of "+file_name, transform=ax.transAxes)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')

    annot = ax.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                        bbox=dict(boxstyle="round", fc="w"),
                        arrowprops=dict(arrowstyle="-|>"))
    annot.set_visible(False)

    fig.colorbar(sc, shrink=0.5, aspect=5)

    def update_annot(ind):

        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = ""
        if (len(ind["ind"]) > 1):
            for n in ind["ind"] :
                text = text+" "+"['"+names[n]+"']"
        else :
            text = names[ind["ind"]]

        annot.set_text(text)

    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", hover)
    
    end = time.time()
    

    plt.savefig("result/"+file_name+"_fig.pdf")
    plt.show()

    return(end)


if __name__ == "__main__":

    print("this programm doesn't run independently sorry !")

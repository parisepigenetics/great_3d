"""
.. module:: visualisation
   :synopsis: This module implements all the functions to parse either input
                or additional necessary files.
"""

# Third-party modules
import matplotlib.cm as cmx
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.colors as col

def visualize_4d_genome(transcription_map, coordinates, colorsMap='YlOrRd'):
    cs = transcription_map
    cm = plt.get_cmap(colorsMap)
    cNorm = col.Normalize(vmin=min(cs), vmax=max(cs))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=cm)
    fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    ax = Axes3D(fig)
    ax.scatter(coordinates['X'], coordinates['Y'], coordinates['Z'],
               c=scalarMap.to_rgba(cs), marker='o')
    scalarMap.set_array(cs)
    fig.colorbar(scalarMap)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()

"""
.. module:: visualisation
   :synopsis: This module implements the 3D visualization part
"""

# Third-party modules
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cmx
import matplotlib.pyplot as plt
import matplotlib.colors as col

def visualize_4d_genome(transcription_map, coordinates, colors_map='YlOrRd'):
    """
    Plot a 3D visualization of the genes colored by their transcription map.

    Args:
        transcription_map (list): Transcription map for the studied genes
        coordinates (Pandas Dataframe): X, Y, Z Coordinates of the genes
        colors_map (str): Chosen colors of the color bar. Default: YlOrRd (Yellow or Red)
    """
    col_map = plt.get_cmap(colors_map)
    col_norm = col.Normalize(vmin=min(transcription_map), vmax=max(transcription_map))
    scalar_map = cmx.ScalarMappable(norm=col_norm, cmap=col_map)
    fig = plt.figure()
    axes = Axes3D(fig)
    axes.scatter(coordinates['X'], coordinates['Y'], coordinates['Z'],\
                 c=scalar_map.to_rgba(transcription_map), marker='o')
    scalar_map.set_array(transcription_map)
    fig.colorbar(scalar_map)
    axes.set_xlabel('X Label')
    axes.set_ylabel('Y Label')
    axes.set_zlabel('Z Label')
    plt.show()

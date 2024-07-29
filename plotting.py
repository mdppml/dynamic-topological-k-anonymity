import os
os.environ["MKL_CBWR"] = "AUTO"
import numpy as np
import matplotlib.pyplot as plt
import time
import cechmate as cm
from persim import plot_diagrams
import tadasets
from collections import defaultdict
from itertools import combinations


def plot_cech(points, maxdim=1):
    """
    Plots the given points and their Cech complex up to a specified dimension.

    Args:
        points (np.ndarray): An array of point coordinates.
        maxdim (int): The maximum dimension for the Cech complex. Defaults to 1.

    Returns:
        list: List of persistence diagrams for each dimension.
    """
    plt.scatter(points[:, 0], points[:, 1])
    plt.title("Generated 15x2 Points")
    plt.show()
    tt=time.time()
    max_radius=25
    
    cech = cm.Cech(maxdim=1)  # Go up to 1D homology
    simplices=cech.build(points)
    dgmscech = cech.diagrams()
    
    plot_diagrams(dgmscech, labels=['$H_0$ Cech', '$H_1$ Cech'])
    plt.show()
    return dgmscech 
  
def plot_weighted_barcodes(shifted_results_2, sorted_vertices, visual limit):
    """
    Plots weighted barcodes for given results.

    Args:
        shifted_results_2 (dict): The processed history of vertex connected components.
        sorted_vertices (list): List of vertices sorted by a specific criteria.
        visual_limit (float): Maximum limit for visualization on x-axis.
    """
    fig, ax = plt.subplots()
    
    for idx, vertex in enumerate(sorted_vertices):
        history = shifted_results_2[vertex]
        last_death_time = 0
        for component, death_time in history:
            thickness = len(component) * 0.5
    
            if death_time == float('inf'):
                ax.plot([last_death_time, visual_limit], [idx, idx], linewidth=thickness * 2, solid_capstyle='butt', color='red')
                ax.scatter(visual_limit, idx, marker='x', color='white', s=thickness * 10, zorder=5)
            else:
                ax.plot([last_death_time, death_time], [idx, idx], linewidth=thickness * 2, solid_capstyle='butt', color='red')
    
            last_death_time = death_time
    
    ax.set_xlabel("Filtration Parameter")
    ax.set_ylabel("Sorted Vertex")
    ax.set_yticks(range(len(sorted_vertices)))
    ax.set_yticklabels(sorted_vertices)
    ax.set_xlim([0, visual_limit])
    x_limits = plt.xlim()
    
    plt.show()



def plot_first_graph(ax, shifted_results_2):
    """
    Plots the first graph with adjusted filtration parameters and sorted vertices - components.

    Args:
        ax (matplotlib.axes._subplots.AxesSubplot): The axes on which to plot.
        shifted_results_2 (dict): The processed history of vertex connected components.
    """
    max_death_time = 0.0
    sorted_vertices = sorted(shifted_results_2.keys(), key=lambda k: shifted_results_2[k][-1][1], reverse=True)

    for vertex in sorted_vertices:
        for _, death_time in shifted_results_2[vertex]:
            if death_time != float('inf'):
                max_death_time = max(max_death_time, death_time)

    visual_limit = max_death_time * 1.1

    for idx, vertex in enumerate(sorted_vertices):
        history = shifted_results_2[vertex]
        last_death_time = 0
        for component, death_time in history:
            thickness = len(component) * 0.5

            if death_time == float('inf'):
                ax.plot([last_death_time, x_limits[1]], [idx, idx], linewidth=thickness * 2, solid_capstyle='butt', color='#DAA520')
                ax.scatter(x_limits[1], idx, marker='x', color='white', s=thickness * 10, zorder=5)
            else:
                ax.plot([last_death_time, death_time], [idx, idx], linewidth=thickness * 2, solid_capstyle='butt', color='#DAA520')

            last_death_time = death_time

    ax.set_xlabel("Filtration Parameter")
    ax.set_xlim([0, x_limits[1]])

def plot_second_graph(ax, holes_data, n):
    """
    Plots the second graph showing the life of holes over the filtration process.

    Args:
        ax (matplotlib.axes._subplots.AxesSubplot): The axes on which to plot.
        holes_data (list): List of dictionaries containing detailed information about each hole.
        n (int): Number of points or vertices in the graph.
    """
    y = n + 1
    spacing = 0.1
    color_map = {1: 'green', 2: 'orange', 3: 'blue'}

    for hole_data in sorted(holes_data, key=lambda x: x['death'][0], reverse=True):
        dim = hole_data['dim']
        birth, birth_vertices = hole_data['birth']
        death, *death_vertices = hole_data['death']
        vertices_changes = hole_data['vertices_changes']

        current_y = y
        color = color_map.get(dim, 'black')

        last_change_point = birth  # Initialize to birth value
        current_vertices = birth_vertices  # Initialize to vertices at birth

        for vertices, changes in vertices_changes:
            for new_vertices, change_point in changes:
                # Update thickness before plotting
                thickness = len(current_vertices) * 0.1
                bar_width = change_point - last_change_point
                plt.barh(current_y, bar_width, left=last_change_point, height=thickness, color=color)

                # Add label on top of the bar
                label_x = last_change_point + bar_width / 2
                label_y = current_y + thickness / 2

                last_change_point = change_point  # Update for the next iteration
                current_vertices = new_vertices  # Update the current vertices for the next iteration

        # Draw the final bar
        thickness = len(current_vertices) * 0.1  # Update thickness one final time
        bar_width = death - last_change_point
        plt.barh(current_y, bar_width, left=last_change_point, height=thickness, color=color)

        # Add label on top of the final bar
        label_x = last_change_point + bar_width / 2
        label_y = current_y + thickness / 2

        y += 1 + spacing  

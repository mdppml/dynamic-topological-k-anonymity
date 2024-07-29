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

from data_generation import generate_cluster, generate_loop
from plotting import plot_cech, plot_weighted_barcodes, plot_first_graph, plot_second_graph
from utilities import get_topological_info,cech_complex_filtration, get_skeleton, get_weighted information, get_holes_data

n=15
np.random.seed(43)

# Generating synthetic data points
cluster1 = generate_loop(center=[0, 0], num_points=5, radius=2)
cluster2 = generate_cluster(center=[4, 5], num_points=5, radius=3)
cluster3 = generate_cluster(center=[-5, 5], num_points=5, radius=2)

# Combining clusters and loops
points = np.vstack([cluster1, cluster2,cluster3])

# Plot Cech complex and compute persistence diagrams
dgmscech=plot_cech(points)

# Extract topological features from the persistence diagrams
unique_points, zero_dim_count, one_dim_count, births, deaths, simplices = get_topological_info(points, dgmscech)
    
lists = []
skeleton1=get_skeleton(simplices,1)
for simplex in skeleton1:  # Up to 1-dimensional skeleton
    vertices, filtration_value = simplex
    if len(vertices) == 2:  # Check if the simplex is an edge
        lists.append(simplex)
del skeleton1

# Sort the list of edges by their filtration value
lists = sorted(lists, key=lambda x: x[1])

result = cech_complex_filtration(lists)

# Get weighted information for plotting
shifted_results_2,sorted_vertices,visual_limit = get_weighted_information()
plot_weighted_barcodes(shifted_results_2, sorted_vertices, visual limit)

diag = []

# Collecting and sorting birth-death pairs from persistence diagrams
for dim, arr in enumerate(dgmscech):
    # iterate through each birth, death pair in that dimension
    for birth, death in arr:
        # create a tuple (dim, (birth, death)) and append it to diag
        diag.append((dim, (birth, death)))

diag.sort(key=lambda x: (x[0], x[1][0]))

# Get filtration data about the holes and where they form
holes_data=get_holes_data(diag)

if x_limits[1]<3:
    x_limits=(0,3)

fig, ax = plt.subplots()
plot_first_graph(ax, shifted_results_2)
plot_second_graph(ax, holes_data, n)
ax.set_xlim(0, 3)
ax.set_aspect('auto')
plt.show()

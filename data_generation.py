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

# Minimum distance parameter to ensure generated points are not too close
min_distance = 0.2

def generate_cluster(center, num_points, radius=1.0, min_distance=min_distance):
        """
    Generates a cluster of points around a center within a specified radius.

    Ensures that no two points are closer than a specified minimum distance.

    Parameters:
    - center: A tuple or list representing the x and y coordinates of the cluster center.
    - num_points: The number of points to generate in the cluster.
    - radius: The maximum distance of any point from the center. Default is 1.0.
    - min_distance: The minimum allowable distance between any two points.

    Returns:
    - A numpy array of shape (num_points, 2), each row representing x and y coordinates of a point.
    """
    points = []
    while len(points) < num_points:
        angle = np.random.rand() * 2 * np.pi
        radii = radius * np.sqrt(np.random.rand())
        x = center[0] + radii * np.cos(angle)
        y = center[1] + radii * np.sin(angle)

        # Check if the point is too close to existing points
        too_close = False
        for px, py in points:
            if np.sqrt((px - x)**2 + (py - y)**2) < min_distance:
                too_close = True
                break

        if not too_close:
            points.append([x, y])

    return np.array(points)

def generate_loop(center, num_points, radius=1.0):
     """
    Generates points in a loop (circle) shape.

    Parameters:
    - center: A tuple or list representing the x and y coordinates of the loop's center.
    - num_points: The number of points to generate along the loop.
    - radius: The radius of the loop.

    Returns:
    - A numpy array of shape (num_points, 2), each row representing x and y coordinates of a point.
    """
    angles = np.linspace(0, 2*np.pi, num_points, endpoint=False)
    x = center[0] + radius * np.cos(angles)
    y = center[1] + radius * np.sin(angles)
    return np.column_stack((x, y))


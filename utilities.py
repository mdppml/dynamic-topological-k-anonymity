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

# Number of unique points
def get_topological_info(points,dgmscech):
  """
    Extracts unique points and topological information from a set of points and 
    their corresponding persistence diagrams (dgmscech).

    Args:
        points (np.ndarray): Array of point coordinates.
        dgmscech (list): List of persistence diagrams for each dimension.

    Returns:
        tuple: Contains unique points, counts of 0D and 1D homology classes, 
        birth times, death times, and sorted simplices.
    """
  unique_points = np.unique(points, axis=0)
  
  # Count 0D and 1D homology classes
  zero_dim_count = sum(1 for dim, _ in dgmscech[0])
  one_dim_count = sum(1 for dim, _ in dgmscech[1])
  
  simplices = sorted(simplices, key=lambda x: x[1])
  births = [sub_array[0] for sub_array in dgmscech[1]]
  deaths = [sub_array[1] for sub_array in dgmscech[1]]
  return unique_points, zero_dim_count, one_dim_count, births, deaths, simplices



def update_history(hist, new_comp, time):
  """
    Updates the history of a vertex's connected components over time.

    Args:
        hist (dict): History dictionary to update.
        new_comp (list): New connected component.
        time (float): Time at which the update happens.
    """
    for v in new_comp:
        hist[v].append((new_comp[:], time))

def cech_complex_filtration(arr):
     """
    Constructs a Cech complex filtration from an array of simplices.

    Args:
        arr (list): List of tuples representing simplices and their filtration values.

    Returns:
        defaultdict: A history of the connected components over time for each vertex.
    """
    parent = {}
    rank = {}
    history = defaultdict(list)

    def find(v):
        if parent[v] != v:
            parent[v] = find(parent[v])
        return parent[v]

    def union(x, y, time):
        root_x = find(x)
        root_y = find(y)

        if root_x != root_y:
            if rank[root_x] < rank[root_y]:
                root_x, root_y = root_y, root_x
            parent[root_y] = root_x
            rank[root_x] = max(rank[root_x], rank[root_y] + 1)

            new_component = [v for v in parent.keys() if find(v) == root_x]
            update_history(history, new_component, time)

    # Initialize disjoint sets and history
    for a, (vertices, time) in enumerate(arr):
        for v in vertices:
            if v not in parent:
                parent[v] = v
                rank[v] = 0
                history[v].append(([v], time))

    # Create the filtration
    for vertices, time in arr:
        x, y = vertices
        union(x, y, time)

    return history


def get_skeleton(simplices, dim):
    """Retrieve the skeleton up to the given dimension from a list of simplices.

    Args:
        simplices (list): List of simplices as produced by Cech.build()
        dim (int): The maximum dimension of the skeleton

    Returns:
        list: A list of simplices up to the given dimension
    """
    return [simplex for simplex in simplices if len(simplex[0]) <= dim + 1]
    
def get_filtration(simplices):
    return [simplex for simplex in simplices]



shifted_results = {}

def shift_floats(vertex_history):
    if len(vertex_history) == 0:
        return []

    shifted_floats = [vertex_history[i][1] for i in range(1, len(vertex_history))] + [float('inf')]

    shifted_history = [(lst, val) for (lst, _), val in zip(vertex_history, shifted_floats)]

    return shifted_history

# Process each vertex's history
for k, v in result.items():
    shifted_results[k] = shift_floats(v)


def prune_vertex_history(vertex_dict):
    # Collect all unique paths across all vertices
    all_paths = {}
    for vertex, history in vertex_dict.items():
        for path, value in history:
            all_paths[tuple(sorted(path))] = value

    new_vertex_dict = {}
    for vertex, history in vertex_dict.items():
        new_history = []
        for path, value in history:
            if tuple(sorted(path)) in all_paths:
                new_history.append((path, value))
                del all_paths[tuple(sorted(path))]  # Remove this path to prevent it from appearing in another vertex's history
        new_vertex_dict[vertex] = new_history

    return new_vertex_dict



def get_weighted_information(shifted_results, dgmscech):
    """
    Calculates weighted information for visualization from the shifted results and 
    persistence diagrams.

    Args:
        shifted_results (dict): Dictionary of shifted results per vertex.
        dgmscech (list): List of persistence diagrams.

    Returns:
        tuple: Contains updated shifted results, sorted vertices, and visual limit for plotting.
    """
  shifted_results_2=prune_vertex_history(shifted_results.copy())
  
  # Initialize max_death_time
  max_death_time = 0.0
  
  # Sort vertices by the death time of their last component
  sorted_vertices = sorted(shifted_results_2.keys(), key=lambda k: shifted_results_2[k][-1][1], reverse=True)
  
  # Find the maximum finite death time
  for vertex in sorted_vertices:
      for _, death_time in shifted_results_2[vertex]:
          if death_time != float('inf'):
              max_death_time = max(max_death_time, death_time)
  
  # Add 10% to the maximum finite death time for x-axis limit
  visual_limit = max(max_death_time * 1.1,dgmscech[1][len(dgmscech[1])-1][1])
  return shifted_results_2,sorted_vertices,visual_limit

def get_simplices_in_filtration_range(simplices, birth, death, dim=None):
    result=[]
    for simplex, value in simplices:
        if (dim is None or len(simplex) - 1 == dim) and birth <= value <= death:
            result.append((simplex, value))
    return result

def has_cycle(graph):
    g = defaultdict(list)
    for u, v in graph:
        g[u].append(v)
        g[v].append(u)

    visited = set()
    def dfs(u, parent):
        visited.add(u)
        for v in g[u]:
            if v not in visited:
                if dfs(v, u):
                    return True
            elif v != parent:
                return True
        return False

    for u in g.keys():
        if u not in visited:
            if dfs(u, None):
                return True
    return False

from itertools import combinations
from heapq import heappop, heappush

def is_subset(list1, list2):
    return set(list1).issubset(set(list2))

def find_shortest_loop(edges,component,faces,death_vertex):
    shortest_loop = None
    shortest_distance = math.inf

    # Generate all combinations of 4 edges to check for loops
    for comb in combinations(edges, 4):
        points = {}

        # Count the occurrences of each vertex
        for edge in comb:
            a, b = edge[0]
            if a not in points:
                points[a] = 0
            if b not in points:
                points[b] = 0
            points[a] += 1
            points[b] += 1
        comb_set=[]

        for edge in comb:
           comb_set.extend(edge[0])

        # Remove duplicates by converting to a set and back to a list
        comb_set = list(set(comb_set))



        # Check if the combination forms a loop (each vertex appears exactly twice)
        if all(val == 2 for val in points.values()):
            skip = False
            for lst in faces:
                if is_subset(lst, comb_set):
                    skip=True

            if skip == False:
                # Calculate the total distance of the loop
                total_distance = sum(edge[1] for edge in comb)
                first=False
                second=False
                for k in range(len(comb)):
                    if component[0] in comb[k][0]:
                        first = True
                        break
                for k in range(len(comb)):
                    if component[1] in comb[k][0]:
                        second = True
                        break
                if first == True and second == True:
                    if total_distance < shortest_distance:
                        shortest_distance = total_distance
                        shortest_loop = comb

    if shortest_loop is not None:
        # Reformat the loop to the format [a, b, c, d]
        loop_vertices = []
        for edge in shortest_loop:
            loop_vertices.extend(edge[0])
        loop_vertices = list(set(loop_vertices))
        return loop_vertices
    else:
        return death_vertex


def quality_check(simplices, death, death_vertex, component, dim):
    # Step 1: Get all simplices
    all_simplices = get_simplices_in_filtration_range(simplices, 0, death-0.0001)
    # Step 2: Filter simplices
    filtered_simplices = [s for s in all_simplices if all(v in death_vertex for v in s[0]) and len(s[0]) == 2]
    filtered_faces = [s[0] for s in all_simplices if all(v in death_vertex for v in s[0]) and len(s[0]) > dim]
    death_vertex2 = find_shortest_loop(filtered_simplices,component, filtered_faces, death_vertex)
    return death_vertex2

def get_holes_data(diag): 
    """
    Extracts detailed information about the holes (homology classes) from the persistence diagram.

    Args:
        diag (list): Sorted list of tuples containing dimension and birth-death pairs.

    Returns:
        list: List of dictionaries containing detailed information about each hole.
    """
    holes_data = []
    for dim, (birth, death) in diag:
        if dim >0:
            # Capture death simplices
            birth_simplices = get_simplices_in_filtration_range(simplices, birth, birth)
            death_simplices = get_simplices_in_filtration_range(simplices, death, death)
            # Capture all simplices between birth and death
            all_simplices = get_simplices_in_filtration_range(simplices, birth+0.00001, death-0.000001)
    
    
            union_set = set()
            # Loop through the array and perform the union operation
            for inner_list, _ in death_simplices:
                union_set |= set(inner_list)
            # Convert the set back to a list and sort it
            death_vertices = sorted(list(union_set))
            death_vertices= quality_check(simplices, death, death_vertices, death_simplices[0][0], dim)
            birth_vertices = death_vertices.copy()
            important_vertices = death_vertices.copy()
            intermediate = []
            intermediate.append((death_vertices.copy(),[death_simplices[0]]))
            # Important vertices are a copy of death vertices
    
            for simplex, _ in reversed(all_simplices):
    
                if set(simplex).issubset(set(important_vertices)):
                    # Find other simplices with the same filtration value
                    matching_simplices = [(s, f) for s, f in all_simplices if s != simplex and f == _]
                    # Identify new vertices not in important_vertices
                    new_vertices = [v for v in matching_simplices[0][0] if v not in important_vertices]
    
                    # Add new vertices to important_vertices and birth_vertices
                    important_vertices += new_vertices
                    birth_vertices += new_vertices
    
                    # Record the new state of important_vertices in intermediate
                    intermediate.append((important_vertices.copy(),matching_simplices))
            # Print the intermediate results
    
            birth_vertices = (birth_vertices,[birth_simplices[0]])
            death_vertices = (death_vertices.copy(),[death_simplices[0]])
            # Print the final state of important_vertices and birth_vertices
    
            hole_data = {
                'dim' : dim,
                'birth': birth,
                'death': death,
                'vertices_changes': intermediate,
            }
            holes_data.append(hole_data)
  return holes_data


def find_shared_edges(s1, s2):
    """
    Finds shared edges between two sets of simplices.

    Args:
        s1 (set): A set of simplices.
        s2 (set): Another set of simplices.

    Returns:
        set: A set of edges that are common to both s1 and s2.
    """
    edges1 = {frozenset(e) for e in combinations(s1, 2)}
    edges2 = {frozenset(e) for e in combinations(s2, 2)}
    return edges1 & edges2


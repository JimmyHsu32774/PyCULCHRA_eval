import math
import numpy as np
from scipy.spatial import KDTree

#Author: Caroline Rhoads (C1)
def find_close_pairs(coords: list[list[float]], cutoff: float, excluded_pairs: set[tuple[int, int]]) -> np.ndarray:
    """C1: Implement the contract above; document your reasoning.
    Implemented by: Caroline Rhoads
    
    Inputs:
        coords: list[list[float]] = a list of 3D coordinates for atoms, where each inner list contains three floats representing the x, y, and z coordinates of an atom.
        cutoff: float = the distance threshold below which atoms are considered to be in close proximity.
        excluded_pairs: set[tuple[int, int]] = a set of atom index pairs that should be excluded from consideration.
    Outputs:
        clash_arr: ndarray = an array of shape (m, 2) containing the indices of two atoms that are closer than the cutoff distance and are not in the excluded_pairs set.
    Reasoning:
        1. Input validation: The function first checks if the cutoff is a finite positive number and if the coords input is a 2D array with shape (n, 3) and contains finite values. If any of these checks fail, a ValueError is raised.
        2. Distance calculation: The function iterates through all unique pairs of atoms (i, j) where i < j. For each pair, it checks if the pair is in the excluded_pairs set.
        3. If the pair is not excluded, it calculates the squared distance between the two atoms using the Euclidean distance formula. If the squared distance is less than the square of the cutoff, the pair is added to the clash_list.
        4. Finally, the function returns the clash_list containing all pairs of atoms that are closer than the cutoff distance and not excluded.
    """
    
    # Input validation:
       
    try:
        cutoff = float(cutoff)
    except (ValueError, TypeError) as error:
        raise ValueError("cutoff passed into find_close_pairs couldn't be converted into float type") from error
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff passed into find_close_pairs is not a finite positive number")
    
    arr = np.asarray(coords, dtype = float)    #np.asarray() doesn't copy the data over if it's already in an array which we won't need
    
    if arr.ndim != 2:
        raise ValueError("the dimension of the atom list passed into find_close_pairs is not 2")
    if arr.shape[1] != 3:
        raise ValueError("the shape of the atom list passed into find_close_pairs is not n x 3")
    if not np.isfinite(arr).all():
        raise ValueError("the entries in the atom list passed into find_close_pairs must be finite")

    # build sciPy KDTree for efficient search using query_pairs()
    tree = KDTree(arr)
    clash_arr = tree.query_pairs(r = cutoff, output_type='ndarray') # all other parameters are default, this stores the output as an ndarray of shape (m, 2) where m is the number of pairs found
    if len(clash_arr) == 0:
        return np.empty((0, 2), dtype=int) # if no pairs are found, return an empty array of shape (0, 2)

    # query_pairs() uses <= cutoff, but we want < cutoff, so we filter out the pairs that are exactly at the cutoff distance:
    angstroms = np.linalg.norm(arr[clash_arr[:, 0]] - arr[clash_arr[:, 1]], axis=1) # finds euclidean distance between each pair of atoms in clash_arr
    clash_arr = clash_arr[angstroms < cutoff]   # returns a list of bools aligned with clash_arr.
    
    # Remove excluded pairs
    clash_arr = np.array([
        pair for pair in clash_arr
        if (pair[0], pair[1]) not in excluded_pairs
        and (pair[1], pair[0]) not in excluded_pairs
    ]).reshape(-1, 2)

    return clash_arr
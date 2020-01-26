"""
measure.py

This module is for functions that take measurements.
"""
import numpy as np
#import atom_data
from .atom_data import atomic_weights

def calculate_molecular_mass(symbols):
    """Calculate the mass of a molecule.
    
    Parameters
    ----------
    symbols : list
        A list of elements.
    
    Returns
    -------
    mass : float
        The mass of the molecule
    """

    mass = 0.0

    for atom_name in symbols:
        if atom_name in atomic_weights.keys():
            mass = mass + atomic_weights[atom_name]
        else:
            raise KeyError("Element " + atom_name + " is not stored in the atom_data dict!!!" ) 
    return mass

def calculate_center_of_mass(symbols, coordinates):
    """
    Calculate the center of mass of a molecule.
    
    The center of mass is weighted by each atom's weight.
    
    Parameters
    ----------
    symbols : list
        A list of elements for the molecule
    coordinates : np.ndarray
        The coordinates of the molecule.
    
    Returns
    -------
    center_of_mass: np.ndarray
        The center of mass of the molecule.
    Notes
    -----
    The center of mass is calculated with the formula
    
    .. math:: \\vec{R}=\\frac{1}{M} \\sum_{i=1}^{n} m_{i}\\vec{r_{}i}
    
    """

    sum_of_mass_times_coordinates = np.array([0.,0.,0.])
    sum_of_mass = calculate_molecular_mass(symbols)

    for atom_name, atom_coordinate in zip(symbols,coordinates):
        sum_of_mass_times_coordinates += atomic_weights[atom_name] * atom_coordinate
    center_of_mass = sum_of_mass_times_coordinates / sum_of_mass
        
    return center_of_mass


def calculate_distance(rA, rB):
    """
    Calculate the distance between two points.

    Parameters
    ----------
    rA, rB : np.ndarray
        The coordinates of each point.

    Returns
    -------
    distance : float
        The distance between the two points.

    Examples
    --------
    >>> r1 = np.array([0,0,0])
    >>> r2 = np.array([0,0.1,0])
    >>> calculate_distance(r1,r2)
    0.1
    """
    # This function calculates the distance between two points given as numpy arrays.
    
    if not isinstance(rA, np.ndarray) or not isinstance(rB, np.ndarray):
        raise TypeError("Input must be type np.ndarray for calculate_distance!")

    distance_vector = (rA - rB)
    distance_scalar = np.linalg.norm(distance_vector)
    return distance_scalar

def calculate_angle(rA, rB, rC, degrees=False):
    # Calculate the angle between three points. Answer is given in radians by default, but can be given in degrees
    # by setting degrees=True
    AB = rB - rA
    BC = rB - rC
    theta=np.arccos(np.dot(AB, BC)/(np.linalg.norm(AB)*np.linalg.norm(BC)))

    if degrees:
        return np.degrees(theta)
    else:
        return theta


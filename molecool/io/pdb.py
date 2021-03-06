"""
pdb.py

This module handles pdb type of data.
"""
import numpy as py

def open_pdb(file_location):
    """
    Open and read coordinates and symbols from a pbd file.

    The pbd file must specify the atom elements in the last column, and follow the conventions outlined in the PDB format specification.

    Parameters
    ----------
    file_location : str
        The directory of the file.

    Returns
    -------
    symbols : str
        The name of the element.

    coordinates : np.ndarray
        3D list containing the xyz coordinates of the atom.

    Examples
    --------
    >>> file_location = ""
    >>> open_pdb(file_location)
    >>> 
    """
    # This function reads in a pdb file and returns the atom names and coordinates.
    with open(file_location) as f:
        data = f.readlines()

    coordinates = []
    symbols = []
    for line in data:
        if 'ATOM' in line[0:6] or 'HETATM' in line[0:6]:
            symbols.append(line[76:79].strip())
            atom_coordinates = [float(x) for x in line[30:55].split()]
            coordinates.append(atom_coordiantes)

    # Convert list to numpy array
    coordinates = np.array(coordinates)

    return symbols, coordinates


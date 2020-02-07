# Some utility classes to represent a PDB structure
import numpy as np

class Atom:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type):
        self.type = type
        self.coords = (0.0, 0.0, 0.0)

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.type

class Residue:
    """
    A simple class for an amino acid residue
    """

    def __init__(self, type, number):
        self.type = type
        self.number = number
        self.atoms = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return "{0} {1}".format(self.type, self.number)
    # find and return the alpha carbon for the residue
    def get_ca(self):
        for atom in self.atoms:
            if atom.type == "CA":
                return atom

class ActiveSite:
    """
    A simple class for an active site
    """

    def __init__(self, name):
        self.name = name
        self.residues = []

    # Overload the __repr__ operator to make printing simpler.
    def __repr__(self):
        return self.name
    # Find the length,width,height of the active site
    def get_shape(self):
        x_coords = list()
        y_coords = list()
        z_coords = list()
        for residue in self.residues:
            ca = residue.get_ca()
            x_coords.append(ca.coords[0])
            y_coords.append(ca.coords[1])
            z_coords.append(ca.coords[2])
        length = max(x_coords) - min(x_coords)
        width = max(y_coords) - min(y_coords)
        height = max(z_coords) - min(z_coords)
        return (length,width,height)


        

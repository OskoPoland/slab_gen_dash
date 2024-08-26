import numpy as np
from ase import Atoms
from ase.neighborlist import natural_cutoffs
from scipy.spatial.distance import pdist, squareform


class Graph:
    """
    Class with methods for generating simple graph represtations of an atomic
    configuration stored as an adjacency matrix in a np.ndarray.
    """

    def __init__(self, atoms: Atoms, scaling_factor: float = 1.05):
        self._atoms = atoms
        self._mult = scaling_factor
        self._adj_mat = None
        self._nat_cut = None

        n = atoms.get_pbc().sum()
        if n == 3:
            self.atoms_type = 'crystal'
        elif n == 2:
            self.atoms_type = 'surface'
        elif n == 1:
            self.atoms_type = 'chain'
        elif n == 0:
            self.atoms_type = 'molecule'

        self.gen_adj_matrix()

    def gen_adj_matrix(self) -> None:
        """Generates an adjacency matrix representation of bonds within a unit cell.

        Bond distances are based on sets of natural cutoffs provided by ASE. Bonds
        across periodic boundaries are determined using the minimum image convetion.

        Parameters
        ----------
        scaling_factor (float, optional): Scalar for natural cutoff distances

        Returns
        -------
        graph (np.ndarray, shape=NxN): Graph representation of bonds within a unit
        cell stored as a symmetric adjacency matrix, where N = # of atoms in the
        system.

        This method uses numpy broadcasting to accelerate the following loop:

            for i in range(nat):
                tmp = []
                ri = fpos[i, :]
                for j in range(nat):
                    dx = ri - fpos[j, :]
                    dx = dx - np.rint(dx)
                    dr = np.matmul(cell, dx)
                    dr_ = np.linalg.norm(dr)
                    tmp.append(np.linalg.norm(dr_))
                d.append(tuple(tmp))
            d = np.array(d)
        """

        # -- Construct a symmetric NxN matrix of natural cutoffs for atom pairs
        mask = None
        if self.atoms_type in ['molecule', 'chain']:
            self._nat_cut = natural_cutoffs(self._atoms, mult=self._mult)
            mask = np.array([[a+b for a in self._nat_cut]
                            for b in self._nat_cut])

        # -- Construct a symmetric NxN matrix of atom pair distances
        # using the minimum image convention
        fpos = self._atoms.get_scaled_positions()
        cell = np.array(self._atoms.get_cell())

        # -- Vectorized version of the above double for loop
        # Compute r[i, :] - r[j, :]
        x = fpos[:, np.newaxis, :] - fpos[np.newaxis, :, :]
        # Minimum image convention in fractional coordinates
        if self.atoms_type in ['surface', 'crystal']:
            x = x - np.rint(x)
        # Convert back to cartesian coordinates
        x = x @ cell
        # Compute the distances
        d = np.linalg.norm(x, axis=2)

        if mask is None:
            d_min = np.amin(d[np.where(np.abs(d) > 1.E-3)]) * self._mult
            x = range(len(self._atoms))
            mask = np.array([[d_min for a in x] for b in x])

        # -- Build the adjaceny matrix graph representation
        N = len(self._atoms)
        self._adj_mat = (d < mask).astype(np.int8) - np.eye(N, dtype=np.int8)

    @property
    def adj_mat(self):
        """Returns the adjaceny matrix."""
        return self._adj_mat

    @property
    def coord_num(self):
        """Returns the total number of bonds to each atom."""
        return self._adj_mat.sum(axis=0)

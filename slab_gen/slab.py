from ase import Atoms
import numpy as np
from itertools import product as iter_prod
from sys import float_info
from ase.io import read, write
from ase.visualize import view
from ase.build import bulk, surface, fcc100, fcc110, fcc111, fcc211, make_supercell, add_adsorbate
from ase.constraints import FixAtoms
from ase.spacegroup import get_spacegroup
from .graph import Graph
from .npencoder import NpEncoder
import os, errno, json

class Slab:
    def __init__(self, 
                 fname : str, 
                 hkl : list[int],
                 nlay : int = 1,
                 vac_height : float = 10., 
                 symprec=1E-5):
        """Description: Build a slab of the material with the hkl orientation 

        Args:
            fname (str): Location of the desired ?unit cell?
            hkl (list[int]): Chosen miller indicies
            nlay (int, optional): # of unit cells per layer normal to the interface . Defaults to 1.
            vac_height (float, optional): Height of the vaccuum for ASE GUI visualization (Angstrom) -- ASK STEVE. Defaults to 10..
            symprec (_type_, optional): --ASK STEVE--. Defaults to 1E-5.
        """

        self.fname = self.validate_fname(fname)
        self.hkl = hkl
        self.bulk = read(self.fname)

        sg = get_spacegroup(self.bulk, symprec=symprec)
        self.spacegroup = f'{sg.symbol} ({sg.no})'
        self.h, self.k, self.l = hkl
        # -- Determine whether i needed for specific geometry?
        # -- TODO: Clarify comment with steve
        if sg.no >= 143 and sg.no <= 194:
            self.i = -(self.h + self.k)

        # -- deferred attributes
        self.slab = None
        self.slab_center = None
        self.nlay = nlay
        self.vac_height = vac_height
        
        self.gen_slab()
     
    @classmethod
    def from_json(self, path : str) -> Atoms:
        if not os.path.exists(path):
            raise FileNotFoundError(
                errno.ENOENT,
                os.strerror(errno.ENOENT),
                path
            )
        
        with open(path, 'r') as json_data:
            data = json.load(json_data)
            
            params = {
                'fname' : data['fname'],
                'hkl' : data['hkl'],
                'nlay' : data['nlay'],
                'vac_height' : data['vac_height']
            }
            
            return Slab(**params)
    
    def to_json(self, name : str = None, out : str = None) -> None:
        """Serializes the slab object as json

        Args:
            name (str, optional): Name of output json file. Defaults to None.
            out (str, optional): Location of output file. Defaults to None. If None then
            file is output to current working directory
        """
        data = self.serialize()

        if out is not None:
            if os.path.exists(out):
                os.chdir(out)

        default = 'slab'
        with open(f'{name if name is not None else default}.json', 'w') as out:
            json.dump(data, out, indent=4, cls=NpEncoder)
            
    def serialize(self) -> dict[str, any]:
        """Serializes the class as a dictionary
        """
        return {
            'fname' : self.fname,
            'hkl' : self.hkl,
            'nlay' : self.nlay,
            'vac_height' : self.vac_height
        }

    @classmethod
    def deserialize(self, data : dict[str, any]) -> Atoms:
        if data is None: return None
        
        return Slab(**data)

    def validate_fname(self, fname : str) -> str:
        if not os.path.exists(fname):
            raise FileNotFoundError(
                errno.ENOENT,
                os.strerror(errno.ENOENT),
                fname
            )
        else:
            return fname
    
    def export(self, ftype : str, out : str | None = None, name : str | None = None) -> bool:
        """
        Uses ase.write() to return the interface as a chosen filetype
        """
        default = 'img'
        if out is not None:
            if os.path.exists(out): # Default to cwd if no path
                os.chdir(out)
        
        if self.slab is None:
            print('Interface has not yet been generated')
            return False

        fname = f'{name if name is not None else default}.{ftype}'
        write(fname, self.slab)
        return True

    def fix_layers(self, n_layers: int = 2, scaling_factor: float = 1.05) -> None:
        """
        Fixes all atoms of the slab except the top n_layers. A "layer" in this 
        context consists of the set of atoms that are undercoordinated at the 
        surface. Layers are determined iteratively by successively masking the
        top-most layer. This method allows for atoms along a corrugated surface
        to be treated as one single layer.

        Parameters
        ----------
        n_layers: int, 
        """
        # -- Copy the z coordinates of all atoms
        z = self.slab.positions[:, 2]

        # -- Generate a graph of the structure
        graph = Graph(self.slab)
        bonds_graph = graph.gen_adj_matrix(scaling_factor=scaling_factor)
        max_coord = np.max(np.sum(bonds_graph, axis=0))

        exclude = []
        for n in range(n_layers):

            # -- Select atom indices to exclude
            for i, row in enumerate(bonds_graph):
                s = np.sum(row)
                if s < max_coord and s > 0 and z[i] > self.slab_center:
                    exclude.append(i)

            # -- Update the graph by zeroing out the previous layer
            for i in exclude:
                bonds_graph[:, i] = 0
                bonds_graph[i, :] = 0

        # -- Set the atom indices to be fixed
        fix_bottom_layers = FixAtoms(
            indices=[atom.index for atom in self.slab if atom.index not in exclude])
        self.slab.set_constraint(fix_bottom_layers)

    def gen_slab(self):
        """
        Generates a basic slab, adds vaccuum, and applies centering
        """
        self.slab = surface(self.bulk, self.hkl, self.nlay)
        self.slab.center(vacuum=self.vac_height, axis=2)

        z = self.slab.positions[:, 2]
        self.slab_center = 0.5 * (np.max(z) + np.min(z))

        self.sort_atomic_pos()

    def get_cell_param(self, cell: np.ndarray) -> dict:
        """
        Assumes cell vectors are stored columnwise.
        """
        a1 = cell[:, 0]
        a2 = cell[:, 1]
        a3 = cell[:, 2]
        a1_, a2_, a3_ = np.linalg.norm(cell, axis=0)

        c_dict = {}
        c_dict['a1'] = cell[:, 0]
        c_dict['a1_'] = np.linalg.norm(a1)
        c_dict['a2'] = cell[:, 1]
        c_dict['a2_'] = np.linalg.norm(a2)
        c_dict['a3'] = cell[:, 2]
        c_dict['a3_'] = np.linalg.norm(a3)

        conv = 180 / np.pi
        c_dict['alpha'] = np.arccos(np.dot(a2, a3) / (a2_ * a3_)) * conv
        c_dict['beta'] = np.arccos(np.dot(a1, a3) / (a1_ * a3_)) * conv
        c_dict['gamma'] = np.arccos(np.dot(a1, a2) / (a1_ * a2_)) * conv

        return c_dict

    def get_super(self, s):
        """Returns a new Slab object containing the super cell.
        s : array-like
            Can be a tuple 3 integer scalars or an integer transformation matrix.
        """
        from copy import deepcopy
        super_surf = deepcopy(self)

        if len(s) == 3 and isinstance(s[0], int):
            s = np.array(s)
            if np.any(s <= 0):
                raise ValueError("'s' 3-tuple elements must be >0.")
            print('Generating a surface super cell with lattice vector scalars:')
            print(s)
            super_surf.slab *= s

        elif len(s) == 3 and isinstance(s[0], (list, tuple, np.ndarray)):
            diag = np.array([s[i][i] for i in range(3)])
            if np.any(diag == 0):
                raise ValueError("'s' diagonal elements cannot = 0.")
            print('Generating a surface super cell with the transformation matrix:')
            print(np.array(s).astype(int))
            super_surf.slab = make_supercell(super_surf.slab, s)

        else:
            raise ValueError(
                "'s' must be a 3-tuple of integers or a 3x3 integer transformation matrix.")

        super_surf.sort_atomic_pos()

        return super_surf

    def print_cell(self, cell: np.ndarray):

        # cell = np.array(self.slab.cell).T

        c_dict = self.get_cell_param(cell)
        a1 = c_dict['a1']
        a2 = c_dict['a2']
        a3 = c_dict['a3']
        a1_ = c_dict['a1_']
        a2_ = c_dict['a2_']
        a3_ = c_dict['a3_']
        alpha = c_dict['alpha']
        beta = c_dict['beta']
        gamma = c_dict['gamma']

        print()
        print(
            f'  a1 = ( {a1[0]:8.4f}, {a1[1]:8.4f}, {a1[2]:8.4f} ), |a1| = {a1_:8.4f} Ang')
        print(
            f'  a2 = ( {a2[0]:8.4f}, {a2[1]:8.4f}, {a2[2]:8.4f} ), |a2| = {a2_:8.4f} Ang')
        print(
            f'  a3 = ( {a3[0]:8.4f}, {a3[1]:8.4f}, {a3[2]:8.4f} ), |a3| = {a3_:8.4f} Ang')
        print()
        print(
            f'  alpha, beta, gamma = (  {alpha:8.4f}, {beta:8.4f}, {gamma:8.4f}  )')
        print(f'  Volume = {np.linalg.det(cell):8.4f} Ang^3')
        print()

    def sort_atomic_pos(self):
        """
        Performs a lexicographical sort of the atomic positions.
        """
        el = np.array(self.slab.get_chemical_symbols())
        pos = self.slab.get_positions()
        idx = np.lexsort([pos[:, 1], pos[:, 0], pos[:, 2]])
        self.slab.set_positions(pos[idx])
        self.slab.set_chemical_symbols(el[idx])

    def square_cell(self, max_vol: int = 2, verbose: bool = False):
        """
        Attempts to generate a square-like surface cell while allowing the
        cell volume to change up to at most a factor of max_vol.
        """

        print('\nAttempting to find an optimally orthogonal square-like surface cell:\n')
        print(f'Max. volume = {max_vol}\n')

        A = np.array(self.slab.cell).T

        # Print the starting cell features
        print('Starting cell:')
        self.print_cell(A)

        # -- only modify first two vectors
        S = np.eye(3).astype(int)
        x = range(-max_vol, max_vol+1)
        y = range(1, max_vol+1)
        S_opt = np.zeros((3, 3)).astype(int)
        J_min = float_info.max  # set to a huge number
        for a, b, c in iter_prod(y, x, x):
            S[0, 0] = a
            S[1, 0] = b
            S[1, 1] = c
            det = int(np.rint(np.linalg.det(S)))
            if abs(det) > max_vol or det <= 0:
                continue

            # -- Generate the lattice vectors
            B = np.matmul(A, S)
            b1 = B[:, 0]
            b2 = B[:, 1]
            b1_ = np.linalg.norm(b1)
            b2_ = np.linalg.norm(b2)

            # -- Compute the orthogonality measure
            P = np.abs(np.dot(b1, b2) / (b1_ * b2_))

            # -- Compute the congruency measure
            R = 1 - min(b1_/b2_, b2_/b1_)

            # -- Test the objective function
            J = P + R
            if (J < J_min):
                S_opt[:, :] = S
                J_min = J

        super_slab = make_supercell(self.slab, S_opt.T)

        self.sort_atomic_pos()

        # -- Rotate cell + atoms so that a1 = (a, 0, 0)
        cell = np.array(super_slab.cell).T
        a1 = cell[:, 0]
        x = np.array([1, 0, 0])
        R = get_rot_mat(a1, x)

        rot_cell = np.round(R @ cell, 8)
        rot_cell[np.where(np.abs(rot_cell) < 1.E-8)] = 0
        super_slab.set_cell(rot_cell)

        pos = super_slab.get_positions()
        pos = (R @ pos.T).T
        super_slab.set_positions(pos)
        self.slab = super_slab
        self.sort_atomic_pos()

        # Print the final cell features
        print('Super cell:')
        B = self.slab.cell[:, :]
        self.print_cell(B)

        # if verbose:
        #     print('S_opt =')
        #     print(S_opt)
        #     print()

        det = int(np.rint(np.linalg.det(S_opt)))
        print(f'The relative volume of the supercell is {det}\n')

    def add_adsorbate(self, fname):
        """
        Uses ASE's add_adsorbate method to place a molecule on the surface.
        """
        cell = self.slab.get_cell()[:, :]
        center = 0.5 * (cell[0] + cell[1])
        center = tuple(center[0:2])
        adsorbate = read(fname, format='xyz')
        add_adsorbate(self.slab, adsorbate, 10., position=center)
        self.slab.center(vacuum=10.0, axis=2)


def get_rot_mat(a, b):
    """
    Returns the matrix that rotates a onto b.
    """

    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)

    v = np.cross(a, b)
    s = np.linalg.norm(v)
    if abs(s) < 1.E-4:
        return np.eye(3)
    s2 = s*s
    c = np.dot(a, b)

    v_x = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    v_x2 = v_x @ v_x

    R = np.eye(3) + v_x + v_x2 * (1-c)/s2
    return R

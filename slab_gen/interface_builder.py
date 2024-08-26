import numpy as np
from sys import float_info
from itertools import product, combinations
import matplotlib.pyplot as plt
import os
import json
import uuid
import numpy as np
import spglib as spg

from ase.build import make_supercell
# -- from xtb.ase.calculator import XTB
from ase import Atoms
from ase.visualize import view
from ase.io import write
from ase.build.tools import sort as ase_sort

# -- Local imports
from .npencoder import NpEncoder
from .ase_encoder import deserialize_atoms
from .slab import Slab


class InterfaceBuilder:
    """Implements methods for building solid-solid interphases.

    Based on the algorithm implemented in QuantumATK:
    https://docs.quantumatk.com/technicalnotes/interface_builder/interface_builder.html

    Can we use XTB to run structure minimization ?

    https://xtb-docs.readthedocs.io/en/latest/setup.html
    """
    
    def __init__(self, 
                slab1 : Atoms,
                slab2 : Atoms,
                ncand : int = 10,
                nmax : int = 5,
                mmax: int = 6,
                theta_min: float = 0,
                theta_max: float = 180,
                dtheta: float = 4,
                strain_cutoff: float = 0.1,
                theta_cutoff: float = 20,
                strain_weight: float = 1,
                area_weight: float = 0.1,
                angle_weight: float = 1,
                grid_spacing: float = 1.,
                z_shift: float = 2.5
        ):
        """Description: Create an interface by stacking two slabs on top of each other. Parameters below define the shape and orientation of the interface, controlling 
        rotations for generated translations, and controlling weights for strain calculations.

        Args:
            slab1 (Atoms): Chosen interface slab (fixed orientation)
            slab2 (Atoms): Chosen interface slab (rotated to find best interface)
            ncand (int, optional): Number of top candidates chosen. Defaults to 10.
            nmax (int, optional): --ASK PROF NICOLE--. Defaults to 5.
            mmax (int, optional): --ASK PROF NICOLE--. Defaults to 6.
            theta_min (float, optional): Mininum angle of rotation of slab2 on slab1 (Degree). Defaults to 0°.
            theta_max (float, optional): Maximum angle of rotation of slab2 on slab1(Degree). Defaults to 180°.
            dtheta (float, optional): Step size from theta_min to theta_max of rotations of slab2 on slab1 (Degree). Defaults to 4.
            strain_cutoff (float, optional): Max allowed strain calculated from lattice vectors for saved candidates.Rejected candidates are NOT saved. Defaults to 0.1.
            theta_cutoff (float, optional): --ASK STEVE-- (?angle of interface lattice vectors) (Degree). Defaults to 20°.
            strain_weight (float, optional): Weights to optimize for strain. Defaults to 1.
            area_weight (float, optional): Weights to optimize for area. Defaults to 0.1.
            angle_weight (float, optional): Weights to optimize for angle. Defaults to 1.
            grid_spacing (float, optional): Translation step size. Defaults to 1.
            z_shift (float, optional): Space between slabs in the interface (Angstrom). Defaults to 2.5Å.
        """

        self.slab_A = slab1
        self.slab_B = slab2
        self.nat1 = len(slab1.slab)
        self.nat2 = len(slab2.slab)
        
        # Candidate parameters
        self.ncand = ncand
        self.nmax = nmax # x
        self.mmax = mmax # y
        self.theta_min = theta_min
        self.theta_max = theta_max
        self.dtheta = dtheta
        self.strain_cutoff = strain_cutoff
        self.theta_cutoff = theta_cutoff
        self.strain_weight = strain_weight
        self.area_weight = area_weight
        self.angle_weight = angle_weight
        self.grid_spacing = grid_spacing
        self.z_shift = z_shift
        
        # Translation and interface building parameters
        self.super_slab_A = None
        self.super_slab_B = None
        self.composite_cell = None
        self.mesh = None
        self.shape = None
        self.itrans_list = []
        
        self.candidates = []
        self.opt_candidate = {}
        self.interfaces = []
        
        # Generate candidates
        self.gen_candidates()
        # Build single most optimal candidate
        self.build_interfaces()
        
    # -- Getters
    def get_slabA(self) -> Atoms: return self.slab_A
    def get_slabB(self) -> Atoms: return self.slab_B
    def get_NCAND(self) -> int : return self.ncand
    def get_nmax(self) -> int: return self.nmax
    def get_mmax(self) -> int: return self.mmax
    def get_theta_min(self) -> float: return self.theta_min
    def get_theta_max(self) -> float: return self.theta_max
    def get_dtheta(self) -> float: return self.dtheta
    def get_strain_cutoff(self) -> float: return self.strain_cutoff
    def get_theta_cutoff(self) -> float: return self.theta_cutoff
    def get_strain_weight(self) -> float: return self.strain_weight
    def get_area_weight(self) -> float: return self.area_weight
    def get_angle_weight(self) -> float: return self.angle_weight
    def get_grid_spacing(self) -> float: return self.grid_spacing
    def get_z_shift(self) -> float: return self.z_shift
    def get_interface(self, id: int) -> Atoms: return self.interfaces[id]
    def get_all_interfaces(self) -> list[Atoms]: return self.interfaces
    def get_candidate(self, id : int) -> dict: return self.candidates[id]
    def get_all_candidates(self) -> list: return self.candidates
    
    # -- Setters
    def set_slabA(self, slab : Atoms) -> None: self.slab_A = slab
    def set_slabB(self, slab : Atoms) -> None: self.slab_B = slab
    def set_NCAND(self, ncand : int) -> None: self.ncand = ncand
    def set_nmax(self, nmax : int) -> None: self.nmax = nmax
    def set_mmax(self, mmax : int) -> None: self.mmax = mmax
    def set_theta_min(self, theta_min : float) -> None: self.theta_min = theta_min
    def set_theta_max(self, theta_max : float) -> None: self.theta_max = theta_max
    def set_dtheta(self, dtheta: float) -> None: self.dtheta = dtheta
    def set_strain_cutoff(self, strain_cutoff : float) -> None: self.strain_cutoff = strain_cutoff
    def set_theta_cutoff(self, theta_cutoff : float) -> None: self.theta_cutoff = theta_cutoff
    def set_strain_weight(self, strain_weight : float) -> None: self.strain_weight = strain_weight
    def set_area_weight(self, area_weight : float) -> None: self.area_weight = area_weight
    def set_angle_weight(self, angle_weight : float) -> None: self.angle_weight = angle_weight
    def set_grid_spacing(self, grid_spacing : float) -> None: self.grid_spacing = grid_spacing
    def set_z_shift(self, z_shift : float) -> None: self.z_shift = z_shift

    @classmethod
    def from_json(self, path : str) -> Atoms:
        """Builds the Interface Builder class from JSON. Any generated candidates or interfaces are not preserved so they
        will need to be regenerated

        Args:
            path (str): Path to the json file that contains the class attributes

        Returns:
            Atoms: Returns an instance of InterfaceBuilder which is an extension of the ASE Atoms object. Returns None on any failure
        """
        
        if not os.path.exists(path):
            return None
        
        with open(path, 'r') as json_data:
            data = json.load(json_data)
            
            slab_A = deserialize_atoms(data['slab_A'])
            slab_B = deserialize_atoms(data['slab_B'])
            params = {
                'ncand' : data['ncand'],
                'nmax' : data['nmax'],
                'mmax' : data['mmax'],
                'theta_min' : data['theta_min'],
                'theta_max' : data['theta_max'],
                'dtheta' : data['dtheta'],
                'strain_cutoff' : data['strain_cutoff'],
                'theta_cutoff' : data['theta_cutoff'],
                'strain_weight' : data['strain_weight'],
                'area_weight' : data['area_weight'],
                'angle_weight' : data['angle_weight'],
                'grid_spacing' : data['grid_spacing'],
                'z_shift' : data['z_shift']
            }
            
            return InterfaceBuilder(slab_A, slab_B, **params)

    def to_json(self, name : str = None, out : str = None) -> None:
        """Serializes the class as json 
        """
        
        obj_dict = self.serialize()
                
        if out is not None:
            if os.path.exists(out):
                os.chdir(out)
        
        default = "InterfaceBuilder"
        with open(f'{name if name is not None else default}.json', 'w') as out:
            json.dump(obj_dict, out, indent=4, cls=NpEncoder)
            
    def serialize(self) -> dict[str, any]:
        return {
            'slab_A' : self.slab_A.serialize(),
            'slab_B' : self.slab_B.serialize(),
            'ncand' : self.ncand,
            'nmax' : self.nmax,
            'mmax' : self.mmax,
            'theta_min' : self.theta_min,
            'theta_max' : self.theta_max,
            'dtheta' : self.dtheta,
            'strain_cutoff' : self.strain_cutoff,
            'theta_cutoff' : self.theta_cutoff,
            'strain_weight' : self.strain_weight,
            'area_weight' : self.area_weight,
            'angle_weight' : self.angle_weight,
            'grid_spacing' : self.grid_spacing,
            'z_shift' : self.z_shift
        }
    
    @classmethod
    def deserialize(self, data : dict[str, any]) -> Atoms:
        if data is None: return None
        
        slab_A = Slab.deserialize(data['slab_A'])
        slab_B = Slab.deserialize(data['slab_B'])
        
        return InterfaceBuilder(slab_A, slab_B, **data)
        
    def export(self, ftype : str, out : str | None = None, name : str | None = None, id : int = 0) -> bool:
        """Uses ase.write() to return the specified interface as the chosen filetype
        Args:
            ftype (str): Specify the desire file output type. Accepted types include most image formats and JSON
            out (str | None, optional): Output location for the object. Defaults to current working directory.
            name (str | None, optional): Desired file name. Defaults to 'img'.
            id (int): Specifies which interface to export

        Returns:
            bool: True on successful export and False otherwise
        """
        default = 'img'
        if out is not None:
            if os.path.exists(out): # Default to cwd if no path
                os.chdir(out)
        
        if self.interfaces is None:
            return False

        if fname == "json":
            self.as_json(name, out)
            return True
        else:
            fname = f'{name if name is not None else default}.{ftype}'
            write(fname, self.interfaces[id])
            return True
        
    def view(self, id : int = 0) -> None:
        """Uses the ase gui to create a popout viewer for the Atoms interface object

        Args:
            id (int, optional): Use to target a specific interface. Defaults to 0 (most optimal interface). If the id is out of bounds of the generatd interfaces
            dictionary the function will throw an IndexError
        """
        if (id > len(self.interfaces)):
            raise IndexError('Selected ID out of bounds of generated interafaces')
        
        view(self.interfaces[str(id)])

    def write_interfaces(self, output : str = None) -> bool:
        """Will write all generated interfaces to POSCAR files labeled by their translation vector 

        Args:
            ouput (str | None): Path to the output directory that will contain the poscars

        Returns:
            bool: True on successful wirte and False otherwise
        """
        if self.candidates is None:
            return False
        
        if output is not None:
            if os.path.exists(output):
                os.chdir(output)
        
        count = 0
        print(len(self.interfaces))
        for composite_slab in self.interfaces:
            poscar_fname = f'POSCAR_itrans_{self.itrans_list[count]:03d}'
            write(poscar_fname, composite_slab, format="vasp")
            count += 1
        
        return True

    def write_candidates(self, output : str = None) -> bool:
        """Write the top NCAND number of interfaces generated after an optimization
        has been performed. The interfaces will ouput to a Candidates folder in the supplied
        directory. If no directory is supplied the folder is generated in the current working directory
        and if there is no Candidates folder present one will be created.
        Args:
            output (str | None): Location of output directory

        Returns:
            bool: returns True on successful write and False otherwise
        """
        if self.candidates is None or self.opt_candidate is None:
            return False
    
        if output is not None:
            if os.path.exists(output):
                os.chdir(output)
        
        with open('candidates.json', 'w') as cdjson:
            json.dump(self.candidates, cdjson, indent=4, cls=NpEncoder)

    def write_translations(self, output : str = None) -> bool:
        """Gets the surface primitive cell vectors and generate the unique set of fractional translations.
       writen to the desired output directory or to the current working directory
        
        Args:
            output (str | None): Desired output directory

        Returns:
            bool: True on successful write False otherwise
        """
        if output is not None:
            if os.path.exists(output):
                os.chdir(output)
        
        p1, p2 = self.get_prim_cell_vectors(self.super_slab_A, self.super_slab_B)
        
        with open('mesh.dat', 'w') as f:
                f.write(f'{self.shape[0]} {self.shape[1]} -1\n')
                for t in self.mesh:
                    f.write(f'{t[0]:12.6f} {t[1]:12.6f} {t[2]:12.6f}\n')

        with open('prim_cell.dat', 'w') as f:
            f.write(f'{p1[0]} {p1[1]} {p1[2]}\n')
            f.write(f'{p2[0]} {p2[1]} {p2[2]}\n')

        with open('super_cell.dat', 'w') as f:
            c1 = self.composite_cell[:, 0]
            c2 = self.composite_cell[:, 1]
            f.write(f'{c1[0]} {c1[1]} {c1[2]}\n')
            f.write(f'{c2[0]} {c2[1]} {c2[2]}\n')
            
    def generate_translations(self) -> None:
        """Generate translations for interface and populate the self.intrerfaces
        """
        p1, p2 = self.get_prim_cell_vectors(self.super_slab_A, self.super_slab_B)
        self.mesh, _ = self.generate_2d_mesh(np.vstack([p1, p2]))
        
        for itrans, trans_vec in enumerate(self.mesh):
            composite_slab = self.super_slab_A.copy()
            composite_slab.set_cell(self.composite_cell.T)
            
            # -- Get translated replica of super cell B
            temp_B = self.super_slab_B.copy()
            temp_B.translate(trans_vec)
            
            # -- Build the composite slab
            for atom in temp_B:
                composite_slab.append(atom)
            composite_slab = ase_sort(composite_slab)
            
            # -- Rottate the position to fit within the rotated cell
            composite_slab.set_pbc((True, True, True))
            composite_slab.wrap()
            
            self.interfaces.append(composite_slab)
            self.itrans_list.append(itrans)
              
    def build_interfaces(self, id : int = None) -> None:
        """Builds interfaces by combining the two slabs using either a selected candidate
        or the candidate with the highest strain score. CAUTION : Calling this function will overwrite
        the current saved array of interfaces and regenerate and sort based on the current list
        of candidates.

        Args:
            id (int | None): The index location of a chosen interface
        """
        if id is not None:
            if id <= len(self.candidates) and id >= 0:
                candidate = self.candidates[id]
        else:
            candidate = self.candidates[0]
        
        # Isolate candidate information
        R = np.asarray(candidate['R'])
        N_B = np.asarray(candidate['N_B'])
        N_A = np.asarray(candidate['N_A'])
        N_A[2,2] = 1
        
        # Rotate slab B - copy made because slab is rotated from original orientation
        copy_slab_B = self.slab_B.slab.copy()
        cell = (R @ copy_slab_B.cell[:, :].T).T
        pos = (R @ copy_slab_B.positions.T).T
        copy_slab_B.set_cell(cell)
        copy_slab_B.set_positions(pos)
        
        # Perform transformations
        self.super_slab_A = make_supercell(self.slab_A.slab, N_A.T)
        self.super_slab_B = make_supercell(copy_slab_B, N_B.T)
        
        # Translate A down & B up
        self.super_slab_A.center(vacuum=0, axis=2)
        z_min_A = min(self.super_slab_A.positions[:, 2])
        z_max_A = max(self.super_slab_A.positions[:, 2])
        height_A = z_max_A - z_min_A
        
        self.super_slab_B.center(vacuum=0, axis=2)
        z_min_B = min(self.super_slab_B.positions[:, 2])
        z_max_B = max(self.super_slab_B.positions[:, 2])
        height_B = z_max_B - z_min_B
        
        self.super_slab_A.translate([0, 0, -z_min_A])
        self.super_slab_B.translate([0, 0, z_max_A - z_min_B + self.z_shift])

        super_cell_A = self.super_slab_A.cell[:, :].T
        super_cell_B = self.super_slab_B.cell[:, :].T

        c1 = 0.5 * (super_cell_A[:, 0] + super_cell_B[:, 0])
        c2 = 0.5 * (super_cell_A[:, 1] + super_cell_B[:, 1])
        c3 = np.array([0, 0, height_A + height_B + 2*self.z_shift])
        self.composite_cell = np.vstack([c1, c2, c3]).T  # column-order

        # -- Rotate the composite cell so c1 is parallel to [1, 0, 0]
        R = self._get_rot_a_onto_b(self.composite_cell[:, 0], [1, 0, 0])
        tmp = self.composite_cell.copy()
        self.composite_cell = np.round(R @ tmp, 8)
        c1_rot = self.composite_cell[:, 0]
        c2_rot = self.composite_cell[:, 1]
        if c1_rot[0] < 0 and c2_rot[1] < 0:
            R_180 = np.array([[-1, 0, 0], [0, -1, 0], [0, 0, 1]]).astype(int)
            R = R_180 @ R
            self.composite_cell = np.round(R @ tmp, 8)
            c1_rot = self.composite_cell[:, 0]
            c2_rot = self.composite_cell[:, 1]

        self.super_slab_A.set_cell(np.vstack([c1, c2, super_cell_A[2]]),
                              scale_atoms=True)
        pos = self.super_slab_A.get_positions()
        pos = (R @ pos.T).T
        self.super_slab_A.set_positions(pos)
        self.super_slab_A.set_cell(np.vstack([c1_rot, c2_rot, super_cell_A[2]]))

        self.super_slab_B.set_cell(np.vstack([c1, c2, super_cell_B[2]]),
                              scale_atoms=True)
        pos = self.super_slab_B.get_positions()
        pos = (R @ pos.T).T
        self.super_slab_B.set_positions(pos)
        self.super_slab_B.set_cell(np.vstack([c1_rot, c2_rot, super_cell_B[2]]))
        
        # Set translation parameters
        self.super_slab_A = self.super_slab_A
        self.super_slab_B = self.super_slab_B
        
        composite_slab = self.super_slab_A.copy()
        composite_slab.set_cell(self.composite_cell.T)
        
        # Generate Mesh
        p1, p2 = self.get_prim_cell_vectors(self.super_slab_A, self.super_slab_B)
        self.mesh, self.shape = self.generate_2d_mesh(np.vstack([p1, p2]))
        
        self.mesh = np.zeros(3).reshape(1, 3)
        self.interfaces.clear()
        for itrans, trans_vec in enumerate(self.mesh):
            
            # -- Get a translated replica of super cell A
            composite_slab = self.super_slab_A.copy()
            composite_slab.set_cell(self.composite_cell.T)

            # -- Get a translated replica of super cell B
            temp_B = self.super_slab_B.copy()
            temp_B.translate(trans_vec)

            # -- Build the composite slab
            for atom in temp_B:
                composite_slab.append(atom)
            composite_slab = ase_sort(composite_slab)

            # -- Rotate the positions to fit within the rotated cell
            composite_slab.set_pbc((True, True, True))
            composite_slab.wrap()
            self.interfaces.append(composite_slab)
            self.itrans_list.append(itrans)
                
    def gen_candidates(self) -> None:
        """
        TODO:
        – Add a docstring to this function.
        – Create unique labels/code/index for each candidate structure.
        """
        # -- Store the lattice vectors of each cell column-wise
        A_slab1 = np.array(self.slab_A.slab.cell).T
        B_slab2 = np.array(self.slab_B.slab.cell).T

        # -- Generate  list of angles (in radians)
        theta_list = np.arange(self.theta_min, self.theta_max +
                               self.dtheta, self.dtheta).astype(np.float64)
        theta_list *= np.pi / 180
        # -- Generate the complete set of integer transformation matrices
        N = self._get_int_trans_mats()

        # -- Start the structure enumeration
        # NOTE: Picks out candidates to be built later populating the candidates and opt_candidates arrays. Optimizing for each theta in theta_list
        for theta in theta_list:
            # print(f'Testing theta = {theta * 180/np.pi:.1f}°...')
            self._get_candidates(A_slab1, B_slab2, theta, N)
        self._sort_candidates()
       
    def _compute_strain(self, A_t, B_t):
        """
        Returns a dictionary containing the strain elements and average strain.
        """

        # -- Unpack the tuples
        A_s, NA = A_t
        B_s, NB = B_t

        # -- Get the first two cell vectors of both super cells
        a1, a2 = A_s[:, 0], A_s[:, 1]
        b1, b2 = B_s[:, 0], B_s[:, 1]

        eps_11 = np.abs(b1[0] / a1[0]) - 1
        eps_22 = np.abs(b2[1] / a2[1]) - 1
        eps_12 = 0.5 * (b2[0] - (b1[0] / a1[0]) * a2[0]) / a2[1]
        eps_ave = (abs(eps_11) + abs(eps_22) + abs(eps_12)) / 3.

        strain_dict = {'eps_11': eps_11,
                       'eps_22': eps_22,
                       'eps_12': eps_12,
                       'eps_ave': eps_ave
                       }

        return strain_dict

    def _compute_lattice_match_score(
        self,
        A_t: tuple,
        eps_ave: float
    ) -> float:
        """
        Computes an empirical score to measure how closely two cells match.
        """
        A, N = A_t
        a1 = A[:, 0]
        a2 = A[:, 1]
        a1_, a2_, _ = np.linalg.norm(A, axis=0)
        area = np.linalg.norm(np.cross(a1, a2))
        phi = np.arccos(np.dot(a1, a2)/(a1_ * a2_))

        alpha = np.array([15, 30, 45, 60, 90, 120, 150]) * np.pi / 180
        # alpha = np.array([15, 30, 45, 60, 90, 120, 150]) * np.pi / 180
        # Try to keep the angles to 30, 45, 60, 90, 120
        S_phi = np.sum(1 - 7 * np.exp(-(phi - alpha)**2 / 15))
        score = np.exp(-self.strain_weight * eps_ave
                       - self.area_weight * area - self.angle_weight * S_phi)
        # eps_ave, area, angle in json, make little python script to load json file
        # have a couple of candidates, so that we can tune the score
        # write out supercell lattice vectors
        return score

    def _get_int_trans_mats(self):
        """
        Generates a set of integer transformation matrices for building supercells.

        Cell lattice vectors are stored column-wise in a matrix, A = (a_1, a_2, a_3).
        The supercell with lattice vectors S = (s_1, s_2, s_3) are generated via
        S = A N, where N is an integer transformation matrix of the form:

             ( a, 0, 0 )
        N =  ( b, c, 0 )
             ( 0, 0, 1 )

        Parameters
        ----------
        mmax : integer
            Max value for the N_1,1 element.
        nmax : integer
            Max value for the N_2,1 and N_2,2 elements.

        Returns
        -------
        N : np.ndarray, shape = (# of arrays, 3, 3)
            The set of 3x3 integer transformation matrices.

        TODO
        ----
        - Add hint typing for the return type here (np.ndarray)
        """
        x = range(1, self.mmax+1)
        y = range(-self.nmax, self.nmax+1)
        N = np.array([[a, 0, 0, b, c, 0, 0, 0, 1] for a, b, c in product(x, y, y)
                      if 0 not in (a, c)]).astype(int)
        N.shape = (-1, 3, 3)
        return N

    def _get_candidates(
            self,
            A: np.ndarray,
            B: np.ndarray,
            theta: float,
            N_B: np.ndarray
    ) -> None:
        """
        Finds the transformation matrix for B that best aligns with the A lattice.

        Parameters
        ----------
        A : np.ndarray, shape = (3, 3)
            The lattice vectors of cell A stored as columns.
        B : np.ndarray, shape = (3, 3)
            The lattice vectors of cell B stored as columns.
        theta : float
            The angle (in radians) used to rotate cell B.
        N_B : np.ndarray, shape = (N_trans, 3, 3)
            The list of integer transformation matrices used to build
            supercells from the rotated versions of cell B.
        strain_cutoff : float
            A cutoff value to set the largest tolerable average strain.
        theta_cutoff : float
            A cutoff value (in degrees) to set the smallest tolerable angle
            between a1 and a2 vectors.
        area_weight : float
            Weighting parameter for cell area in the interface scoring function.
        angle_weight : float
            Weighting parameter for cell angle in the interface scoring function.

        Returns
        -------
        None
        """
        N_super, _, _ = N_B.shape
        A_inv = np.linalg.inv(A)

        # -- Rotate the B cell vectors about the b_3 axis by theta rad
        R = self._get_rot_z(theta)
        B_rot = R @ B

        # -- Generate a list of supercell lattice vectors
        B_super = B_rot @ N_B
        B_tuple = [(B_super[i], N_B[i]) for i in range(N_super)]

        # -- Compute the integer transformation matrices for the A cell
        #    that approximates the shapes of the corresponding B supercells
        #    and generate a set of A supercells
        N_A = np.rint(A_inv @ B_super).astype(int)
        N_A[2, 2] = 1  # SW: Explicitly set this to 1
        A_super = A @ N_A
        A_tuple = [(A_super[i], N_A[i]) for i in range(N_super)]

        # -- Test all possible pairs of lattice vectors
        tmp_candidates = []
        for idx, (A_t, B_t) in enumerate(zip(A_tuple, B_tuple)):

            # -- Check that the candidate satisfies target geometric criteria
            c_dict = self._verify_interface_geom(A_t, B_t, R)
            if c_dict is None:
                continue

            # -- Generate the strain dictionary
            eps = self._compute_strain(A_t, B_t)

            if abs(eps['eps_ave']) <= self.strain_cutoff:
                score = self._compute_lattice_match_score(A_t, eps['eps_ave'])
                c_dict['eps'] = eps
                c_dict['score'] = score
                c_dict['theta'] = theta * 180 / np.pi
                c_dict['uuid'] = str(uuid.uuid4())
                tmp_candidates.append(c_dict)

        self.candidates += tmp_candidates

        # -- Identify the highest scoring candidate for this specific angle theta
        scores = [x['score'] for x in self.candidates]
        try:
            idx = np.argmax(scores)
        except ValueError:
            raise ValueError('Scores is an empty list')

        new_candidate = self.candidates[idx]
        old_score = self.opt_candidate.get('score')
        new_score = new_candidate.get('score')

        if old_score is None or new_score > old_score:
            self.opt_candidate = new_candidate
            self._print_candidate(new_candidate)

    def _get_rot_z(self, theta: float) -> np.ndarray:
        """Returns a 3D rotation matrix about the z-axis."""
        R_z = np.array([[np.cos(theta), -np.sin(theta), 0],
                        [np.sin(theta), np.cos(theta), 0],
                        [0, 0, 1]]
                       )
        return R_z

    def _get_rot_a_onto_b(self, a: np.ndarray, b: np.ndarray) -> np.ndarray:
        """
        Returns the matrix that rotates a onto b.
        """

        a_hat = a / np.linalg.norm(a)
        b_hat = b / np.linalg.norm(b)

        v = np.cross(a_hat, b_hat)
        s = np.linalg.norm(v)
        if abs(s) < 1.E-4:
            return np.eye(3)
        s2 = s*s
        c = np.dot(a_hat, b_hat)

        v_x = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        v_x2 = v_x @ v_x

        R = np.eye(3) + v_x + v_x2 * (1-c)/s2
        return R

    def _list_top_candidates(self, ):
        """
        Prints a list of the top `ncand` interfacial candidates.
        """
        header = f'\nTop {self.ncand} candidates\n'
        header += '-' * (len(header) - 2)
        print(header)
        for cnt, i in enumerate(range(self.ncand)):
            print(f'Candidate # {cnt+1:4d}:')
            self._print_candidate(self.candidates[i])

    def _print_candidate(self, candidate):
        """
        Prints a formatted summary of a candidate.
        """
        line = f"  score = {candidate['score']:12.4e},    "
        line += f"nat = {candidate['nat']:6d},    "
        line += f"strain = {candidate['eps']['eps_ave']:12.4e},    "
        line += f"area = {candidate['area']:6.2f},    "
        line += f"gamma = {candidate['angle']:6.2f}°,   "
        line += f"theta = {candidate['theta']:6.2f}°"
        print(line)

    def _sort_candidates(self):
        """
        Sorts the candidates in descending order of their cell match
        scores.
        """
        scores = [x['score'] for x in self.candidates]
        idx = np.argsort(scores)
        sorted_candidates = [self.candidates[i] for i in idx]
        self.candidates = sorted_candidates[::-1]

    def _verify_interface_geom(self, A_t, B_t, R, prec=1.E-3):
        """
        Confirms that the proposed interface meets certain geometric criteria.

        If the interface satisfies the criteria, a dictionary containing results
        is returned. Otherwise `None` is returned.
        """

        A_s, N_A = A_t
        B_s, N_B = B_t

        # -- Get the first two cell vectors of both super cells
        a1, a2 = A_s[:, 0], A_s[:, 1]
        b1, b2 = B_s[:, 0], B_s[:, 1]

        # -- Estimate the number of atoms in the interfacial cell
        A_rel_vol = np.abs(np.linalg.det(N_A))
        B_rel_vol = np.abs(np.linalg.det(N_B))
        # The number of atoms reported here might be an estimate, rather than the actual value.
        nat_interface = int(A_rel_vol * self.nat1 + B_rel_vol * self.nat2)

        # -- Compute the cell vector norms
        a1_, a2_, _ = np.linalg.norm(A_s, axis=0)
        b1_, b2_, _ = np.linalg.norm(B_s, axis=0)

        # -- Compute the cross-sectional areas
        A_area = np.linalg.norm(np.cross(a1, a2))
        B_area = np.linalg.norm(np.cross(b1, b2))

        # -- Compute the cell volumes
        A_vol = np.linalg.det(A_s)
        B_vol = np.linalg.det(B_s)

        # -- Compute the angles between cell vectors
        A_theta = np.arccos(np.dot(a1, a2) / (a1_ * a2_))
        A_theta_deg = A_theta * 180 / np.pi

        B_theta = np.arccos(np.dot(b1, b2) / (b1_ * b2_))
        B_theta_deg = B_theta * 180 / np.pi

        # -- Check signs and magnitudes
        tmp = np.array([nat_interface, abs(a1[0]), abs(
            a2[1]), A_area, B_area, A_vol, B_vol])
        if np.any(tmp < prec):
            # then these either have no area or they are left-handed cells
            return None

        # -- Check cell angles
        if B_theta_deg < self.theta_cutoff or B_theta_deg > 180 - self.theta_cutoff:
            return None

        # -- Compile results into a dictionary
        c_dict = {'nat': nat_interface,
                  'score': None,
                  'eps': None,
                  'area': A_area,
                  'angle': A_theta_deg,
                  'R': R,
                  'N_A': N_A,
                  'N_B': N_B
                  }

        return c_dict

    def generate_2d_mesh(self, cell):
        """
        Returns a uniform mesh in the surface cell basis.
        """

        a1, a2 = cell[0], cell[1]

        N = np.rint(np.linalg.norm(cell, axis=1) /
                    self.grid_spacing).astype(int)[:2]
        a1_step = a1 / (N[0])
        a2_step = a2 / (N[1])

        # if verbose:
        #     print(f'Target grid spacing = {grid_spacing}')
        #     print(f'Mesh shape = {N}')
        #     print(f'a1_step length = {np.linalg.norm(a1_step)}')
        #     print(f'a2_step length = {np.linalg.norm(a2_step)}')

        # -- Generate the 2D mesh
        v_mesh = [i * a1_step for i in range(N[0])]
        w_mesh = [i * a2_step for i in range(N[1])]
        mesh = np.vstack([v + w for v in v_mesh for w in w_mesh])
        return mesh, N

    def get_prim_cell_vectors(self, slab_A, slab_B, prec=1.E-4):
        """
        Returns a set of surface primitive cell vectors that corresponds to the
        smallest set of surface fractional translations between slab_A and slab_B.
        """

        # -- Get the supercell lattice vectors
        cell = np.array(slab_B.cell)
        a1 = cell[0]
        a2 = cell[1]

        # -- Extract all fractional translations parallel to interface
        parallel_translations_list = []
        for slab in [slab_A, slab_B]:

            spg_cell = self.get_spg_cell(slab)
            ase_cell = np.array(slab.cell).T

            # if verbose:
            #     spacegroup = spg.get_spacegroup(spg_cell, symprec=1e-4)
            #     line = '-' * 30
            #     print('\n' + line)
            #     print(f'Spacegroup = {spacegroup}')
            #     print(line + '\n')

            # -- Get spacegroup symmetry operations
            sym_dict = spg.get_symmetry(spg_cell, symprec=1e-4)
            rot = sym_dict['rotations']
            trans = sym_dict['translations']  # Note: in crystal coordinates
            # set floating point values close to 0 and 1 to be exactly 0
            trans[np.where(np.abs(trans) < prec)] = 0
            trans[np.where(np.abs(trans - 1) < prec)] = 0

            # -- Get subset of fractional translations parallel to interface
            parallel_translations_cryst = self.get_parallel_trans_vecs(
                rot, trans)
            parallel_translations_cart = (
                ase_cell @ parallel_translations_cryst.T).T
            parallel_translations_list.append(parallel_translations_cart)
            # if verbose:
            #     print()
            #     print('Fractional translations in cartesian units:')
            #     print(parallel_translations_cart)
            #     print()

        parallel_translations_list = np.vstack(parallel_translations_list)

        # -- Sort the translation vectors according to their projections along
        # the a1 and a2 axes
        a1_list = []
        a2_list = []
        for vec in parallel_translations_list:

            # -- Skip to the next vector if current vec is a null vector
            norm = np.linalg.norm(vec)
            if norm < 1.E-4:
                continue

            # -- Sort vec according to which super cell axis it projects onto the most
            proj_1 = self.scalar_projection(vec, a1)
            proj_2 = self.scalar_projection(vec, a2)
            if abs(proj_1 - proj_2) < prec or proj_1 > proj_2:
                a1_list.append(vec)
            elif proj_2 > proj_1:
                a2_list.append(vec)

        # -- Get the primitive cell vector along a1 and a2
        p1 = self.get_shortest_frac_trans_vector(a1_list, a1)
        p2 = self.get_shortest_frac_trans_vector(a2_list, a2)
        # if verbose:
        #     print('The detected primitive surface cell vectors:')
        #     print('\np1 = ', p1)
        #     print('p2 = ', p2, '\n')

        return p1, p2

    def get_shortest_frac_trans_vector(self, ft_vec_list, default):
        """
        Returns the shortest fractional translation vector.

        If the ft_list is empty, returns the default vector (which should be the
        supercell lattice vector).
        """

        if not ft_vec_list:
            return default

        min_norm = 1E6
        min_vec = None
        for vec in ft_vec_list:
            norm = np.linalg.norm(vec)
            if norm < min_norm:
                min_norm = norm
                min_vec = vec

        return min_vec

    def get_parallel_trans_vecs(self, rot, trans, prec=1.E-4):
        """
        Returns an N x 3 array of fractional translations vectors parallel
        to the interface.

        The fractional translations returned here are intended to be associated
        with purely translational spacegroup symmetry operations.
        """

        # -- Get the fractional translations parallel to the a1 - a2 plane
        idx = np.where(np.abs(trans[:, 2]) < prec)
        trans_para = trans[idx]
        rot_para = rot[idx]

        # -- Generate a unique list of parallel translation vectors
        trans_para_list = [trans_para[0]]
        rot_para_list = [rot_para[0]]
        eye3 = np.eye(3).astype(int)
        for r_test, t_test in zip(rot_para[1:, :, :], trans_para[1:]):
            add_row = True

            if not np.all(r_test - eye3 == 0):
                # only want pure translational operations
                continue

            for t in trans_para_list:  # check to see if vector is already present
                if np.all(np.abs(t - t_test) < prec):
                    add_row = False
                    break
            if add_row:
                rot_para_list.append(r_test)
                trans_para_list.append(t_test)

        trans_para_list = np.vstack(trans_para_list)

        # -- Print a summary of the detected parallel translational operations
        # if verbose:
        #     print(
        #         '\nThe unique set of translations (in fractional coordinates) parallel to the a1 - a2 plane:')
        #     eye3 = np.eye(3).astype(int)
        #     for r, t in zip(rot_para_list, trans_para_list):
        #         if np.all(r - eye3 != 0):
        #             print('Warning: Non-identity rotation operation detected.')
        #         print(r)
        #         print(t)
        #         print()

        return trans_para_list

    def scalar_projection(self, a, b):
        """
        Computes the projection of vector a onto vector b
        """
        return np.dot(a, b) / np.linalg.norm(b)

    def plot_surface_cell(self, ax, a1, a2, c='r', ls='-'):
        """
        Plots the outline of the surface cell
        """

        ax.plot([0, a1[0]], [0, a1[1]], c=c, ls=ls, zorder=0)
        ax.plot([0, a2[0]], [0, a2[1]], c=c, ls=ls, zorder=0)
        ax.plot([a2[0], a1[0] + a2[0]],
                [a2[1], a1[1] + a2[1]], c=c, ls=ls, zorder=0)
        ax.plot([a1[0], a1[0] + a2[0]],
                [a1[1], a1[1] + a2[1]], c=c, ls=ls, zorder=0)

        return

    def get_spg_cell(self, slab):
        lattice = np.array(slab.cell)
        frac_pos = slab.get_scaled_positions()
        numbers = slab.numbers
        cell = (lattice, frac_pos, numbers)
        return cell

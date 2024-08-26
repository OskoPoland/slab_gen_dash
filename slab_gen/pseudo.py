import yaml
import os
from slab_gen import PACKAGEDIR


class Pseudo:

    def __init__(self, pseudo_lib, el_list):

        self.lib_name = pseudo_lib.upper()
        assert self.lib_name in ['SSSP_EFF', 'GBRV', 'ONCV', 'CUSTOM', 'DOJO']

        if isinstance(el_list, str):
            self.el_tuple = tuple(el_list)
        elif isinstance(el_list, (list, tuple, set)):
            self.el_tuple = tuple(el_list)

        self.data_path = os.path.join(PACKAGEDIR, 'pseudo_data')

        # deffered attributes
        self.remote_path = None
        self.pseudo_tuple = None

        self.get_pseudo_path()
        self.get_pseudo_list()

    def get_pseudo_dir(self):
        pseudo_path = {
            'SSSP_EFF': '/usr/WS1/woodgrp/corrosion/Pseudopotentials/SSSP_efficiency_pseudos',
            'GBRV': '/usr/workspace/woodgrp/corrosion/Pseudopotentials/GBRV_PBE_UPF_v1.5',
            'ONCV': '/usr/workspace/woodgrp/corrosion/Pseudopotentials/sg15_oncv_upf_2015-10-07',
            'DOJO': '/g/g13/weitzner/Pseudos/Pseudo_DOJO_PAW',
            'CUSTOM': '/usr/workspace/woodgrp/corrosion/Pseudopotentials/Custom_pseudo_set'
        }
        self.path = pseudo_path[self.lib_name]

    def get_pseudo_tuple(self):
        """
        Returns the name of a pseudopotential from a particular library.
        """

        element = element.capitalize()
        el_list = [el.capitalize() for el in self.el_tuple]


# def get_lat_const(pseudo_name, lib_name='SSSP_EFF', XC='PBE'):
#     """
#     Retrieves the converged lattice constant for a metal
#     according to the cutoffs provided in the complementary
#     cutoffs dictionary (below).

#     Convention:

#     degauss = 0.007 Ry
#     smearing = 'mv'
#     k_grid = 18 18 18 0 0 0
#     """

#     XC = XC.upper()
#     lib_name = lib_name.upper()
#     assert XC.lower() in ['pbe', 'rpbe', 'pbesol']
#     assert lib_name.lower() in ['sssp_eff', 'gbrv', 'oncv', 'custom', 'dojo']

#     lat_const_pbe = {
#         'SSSP_EFF': {
#             #
#             # -- SSSP Efficiency v.1.1
#             #
#             'Al.pbe-n-kjpaw_psl.1.0.0.UPF': 4.0355,
#             'Cu_pbe_v1.2.uspp.F.UPF': 3.6285,
#             'Au_ONCV_PBE-1.0.oncvpsp.upf': 4.1501
#         },
#         'GBRV': {
#             'al_pbe_v1.uspp.F.UPF': 4.0447,
#             'au_pbe_v1.uspp.F.UPF': 4.1543,
#             'cu_pbe_v1.2.uspp.F.UPF': 3.6285,
#             'pt_pbe_v1.4.uspp.F.UPF': 3.9671
#         },
#         'ONCV': {
#             'Al_ONCV_PBE-1.0.upf': 4.0449
#         },
#         'CUSTOM': {
#             #
#             # -- Marzari custom generated
#             #    https://doi.org/10.1103/PhysRevB.77.172102
#             #
#             'Al.pbe-sp-van.UPF': 4.0385,
#             'Cu.temp.UPF': 3.67688450446
#         }
#     }

#     lat_const_rpbe = {
#         #
#         # XC = 'RPBE'
#         #
#         'SSSP_EFF': {
#             #
#             # -- SSSP Efficiency v.1.1
#             #
#             'Al.pbe-n-kjpaw_psl.1.0.0.UPF': 4.0650,
#             'Cu_pbe_v1.2.uspp.F.UPF': 3.6725
#         },
#         'GBRV': {
#             'au_pbe_v1.uspp.F.UPF': 4.1926,
#             'cu_pbe_v1.2.uspp.F.UPF': 3.6354,
#             'pt_pbe_v1.4.uspp.F.UPF': 3.9884
#         },
#         'ONCV': {
#             'Al_ONCV_PBE-1.0.upf': None
#         },
#         'DOJO': {
#             'Al.upf': 4.0576
#         },

#         'CUSTOM': {
#             #
#             # -- Marzari custom generated
#             #    https://doi.org/10.1103/PhysRevB.77.172102
#             #
#             'Al.pbe-sp-van.UPF': 4.0621
#         }
#     }

#     lat_const_pbesol = {
#         'SSSP_EFF': {
#             #
#             # -- SSSP Efficiency v.1.1
#             #
#             'Al.pbe-n-kjpaw_psl.1.0.0.UPF': 4.0163,
#             'Cu_pbe_v1.2.uspp.F.UPF': None
#         },
#         'GBRV': {
#             'au_pbe_v1.uspp.F.UPF': 4.0825
#         },
#         'ONCV': {
#             'Al_ONCV_PBE-1.0.upf': None
#         },
#         'CUSTOM': {
#             #
#             # -- Marzari custom generated
#             #    https://doi.org/10.1103/PhysRevB.77.172102
#             #
#             'Al.pbe-sp-van.UPF': None
#         }
#     }

#     if XC.lower() == "pbe":
#         return lat_const_pbe[lib_name][pseudo_name]

#     elif XC.lower() == "rpbe":
#         return lat_const_rpbe[lib_name][pseudo_name]

#     elif XC.lower() == "pbesol":
#         return lat_const_pbesol[lib_name][pseudo_name]


# def get_ecuts(pseudo_name):
#     """
#     Returns a tuple of cutoffs (ecutwfc, ecutrho) for
#     a given pseudopotential.
#     """

#     cutoffs = {
#         #
#         # -- SSSP Efficiency v.1.1
#         #
#         'H.pbe-rrkjus_psl.1.0.0.UPF': (60,  480),
#         'He_ONCV_PBE-1.0.oncvpsp.upf': (50,  400),
#         'li_pbe_v1.4.uspp.F.UPF': (40,  320),
#         'be_pbe_v1.4.uspp.F.UPF': (40,  320),
#         'b_pbe_v1.4.uspp.F.UPF': (35,  280),
#         'C.pbe-n-kjpaw_psl.1.0.0.UPF': (45,  360),
#         'N.pbe-n-radius_5.UPF': (60,  480),
#         'O.pbe-n-kjpaw_psl.0.1.UPF': (50,  400),
#         'f_pbe_v1.4.uspp.F.UPF': (45,  360),
#         'Ne_ONCV_PBE-1.0.oncvpsp.upf': (50,  200),
#         'na_pbe_v1.5.uspp.F.UPF': (40,  320),
#         'Mg.pbe-n-kjpaw_psl.0.3.0.UPF': (30,  240),
#         'Al.pbe-n-kjpaw_psl.1.0.0.UPF': (30,  240),
#         'Si.pbe-n-rrkjus_psl.1.0.0.UPF': (30,  240),
#         'P.pbe-n-rrkjus_psl.1.0.0.UPF': (30,  240),
#         's_pbe_v1.4.uspp.F.UPF': (35,  280),
#         'cl_pbe_v1.4.uspp.F.UPF': (40,  320),
#         'Ar_ONCV_PBE-1.1.oncvpsp.upf': (60,  240),
#         'K.pbe-spn-kjpaw_psl.1.0.0.UPF': (60,  480),
#         'Ca_pbe_v1.uspp.F.UPF': (30,  240),
#         'Sc_ONCV_PBE-1.0.oncvpsp.upf': (40,  160),
#         'ti_pbe_v1.4.uspp.F.UPF': (35,  280),
#         'v_pbe_v1.4.uspp.F.UPF': (35,  280),
#         'cr_pbe_v1.5.uspp.F.UPF': (40,  320),
#         'mn_pbe_v1.5.uspp.F.UPF': (65,  780),
#         'Fe.pbe-spn-kjpaw_psl.0.2.1.UPF': (90, 1080),
#         'Co_pbe_v1.2.uspp.F.UPF': (45,  360),
#         'ni_pbe_v1.4.uspp.F.UPF': (45,  360),
#         'Cu_pbe_v1.2.uspp.F.UPF': (55,  440),
#         'Zn_pbe_v1.uspp.F.UPF': (40,  320),
#         'Ga.pbe-dn-kjpaw_psl.1.0.0.UPF': (70,  560),
#         'ge_pbe_v1.4.uspp.F.UPF': (40,  320),
#         'As.pbe-n-rrkjus_psl.0.2.UPF': (35,  280),
#         'Se_pbe_v1.uspp.F.UPF': (30,  240),
#         'br_pbe_v1.4.uspp.F.UPF': (30,  240),
#         'Kr_ONCV_PBE-1.0.oncvpsp.upf': (45,  180),
#         'Rb_ONCV_PBE-1.0.oncvpsp.upf': (30,  120),
#         'Sr_pbe_v1.uspp.F.UPF': (30,  240),
#         'Y_pbe_v1.uspp.F.UPF': (35,  280),
#         'Zr_pbe_v1.uspp.F.UPF': (30,  240),
#         'Nb.pbe-spn-kjpaw_psl.0.3.0.UPF': (40,  320),
#         'Mo_ONCV_PBE-1.0.oncvpsp.upf': (35,  140),
#         'Tc_ONCV_PBE-1.0.oncvpsp.upf': (30,  120),
#         'Ru_ONCV_PBE-1.0.oncvpsp.upf': (35,  140),
#         'Rh_ONCV_PBE-1.0.oncvpsp.upf': (35,  140),
#         'Pd_ONCV_PBE-1.0.oncvpsp.upf': (45,  180),
#         'Ag_ONCV_PBE-1.0.oncvpsp.upf': (50,  200),
#         'Cd.pbe-dn-rrkjus_psl.0.3.1.UPF': (60,  480),
#         'In.pbe-dn-rrkjus_psl.0.2.2.UPF': (50,  400),
#         'Sn_pbe_v1.uspp.F.UPF': (60,  480),
#         'sb_pbe_v1.4.uspp.F.UPF': (40,  320),
#         'Te_pbe_v1.uspp.F.UPF': (30,  240),
#         'I.pbe-n-kjpaw_psl.0.2.UPF': (35,  280),
#         'Xe_ONCV_PBE-1.1.oncvpsp.upf': (60,  240),
#         'Cs_pbe_v1.uspp.F.UPF': (30,  240),
#         'Ba.pbe-spn-kjpaw_psl.1.0.0.UPF': (30,  240),
#         # -- skipping lanthanides
#         'Hf-sp.oncvpsp.upf': (50,  200),
#         'Ta_pbe_v1.uspp.F.UPF': (45,  360),
#         'W_pbe_v1.2.uspp.F.UPF': (30,  240),
#         'Re_pbe_v1.2.uspp.F.UPF': (30,  240),
#         'Os_pbe_v1.2.uspp.F.UPF': (40,  320),
#         'Ir_pbe_v1.2.uspp.F.UPF': (55,  440),
#         'pt_pbe_v1.4.uspp.F.UPF': (35,  280),
#         'Au_ONCV_PBE-1.0.oncvpsp.upf': (45,  180),
#         'Hg_ONCV_PBE-1.0.oncvpsp.upf': (50,  200),
#         'Tl_pbe_v1.2.uspp.F.UPF': (50,  400),
#         'Pb.pbe-dn-kjpaw_psl.0.2.2.UPF': (40,  320),
#         'Bi_pbe_v1.uspp.F.UPF': (45,  360),
#         'Po.pbe-dn-rrkjus_psl.1.0.0.UPF': (75,  600),
#         'Rn.pbe-dn-kjpaw_psl.1.0.0.UPF': (120,  960),
#         #
#         # -- ONCV norm conserving pseudos
#         #
#         'Al_ONCV_PBE-1.0.upf': (80, 320),
#         #
#         # -- Pseudo-Dojo PAW potentials
#         #
#         'Al.upf': (50, 500),
#         'O.upf': (50, 500),
#         'H.upf': (50, 500),
#         #
#         # -- Marzari custom generated
#         #    https://doi.org/10.1103/PhysRevB.77.172102
#         #
#         'Al.pbe-sp-van.UPF': (100, 800),
#         #
#         # -- GBRV
#         #
#         #
#         'al_pbe_v1.uspp.F.UPF': (40, 400),
#         'au_pbe_v1.uspp.F.UPF': (60, 480),
#         'cu_pbe_v1.2.uspp.F.UPF': (55, 440),
#         'Cu.temp.UPF': (60, 600)
#     }

#     return cutoffs[pseudo_name]

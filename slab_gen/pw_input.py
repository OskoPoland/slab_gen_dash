#!/usr/bin/env python

import sys

from pseudo_properties import get_pseudo_name, get_pseudo_dir, get_lat_const, get_ecuts


class PWInput:

    def __init__(self, element, pseudo_lib):
        self.el = element
        self.pslib = pseudo_lib


def get_pseudo_dict(element, pseudo_lib):
    """
    Returns a dictionary populated with element names
    and corresponding pseudopotential names.
    """
    if isinstance(element, str):
        pseudo_dict = {element: get_pseudo_name(element, lib_name=pseudo_lib)}
    elif isinstance(element, list) or isinstance(element, tuple):
        pseudo_dict = {el: get_pseudo_name(
            el, lib_name=pseudo_lib) for el in element}
    return pseudo_dict


def get_slab(element, hkl, size, vac_height, XC, pseudo_lib, fix_tag=2,
             square_cell=False, prec=1E-4, ortho_111=False,
             return_lat_const=False, highlight_free=False, extract_atom=False,
             adatom=False, ad_molecule=False, adsorbate=None, ads_site=None,
             z_shift=0.):
    """
    Returns a slab centered on z = 0. Layers with tags less than
    fix_tag are fixed completely. (Largest layer tag number has the
    lowest z position.)
    """
    if extract_atom and adatom:
        raise Exception("Cannot currently extract atoms with adatoms present.")
        sys.exit()

    pseudo_name = get_pseudo_name(element, lib_name=pseudo_lib)
    lat_const = get_lat_const(pseudo_name, lib_name=pseudo_lib, XC=XC)
    ecutwfc, ecutrho = get_ecuts(pseudo_name)

    if hkl == (1, 0, 0):
        slab = fcc100(element, size=size, a=lat_const)
    elif hkl == (1, 1, 0):
        slab = fcc110(element, size=size, a=lat_const)
    elif hkl == (1, 1, 1):
        slab = fcc111(element, size=size, a=lat_const, orthogonal=ortho_111)
    elif hkl == (2, 1, 1):
        assert size[0] % 3 == 0, "size[0] must be divisible by 3 for the step"
        slab = fcc211(element, size=size, a=lat_const)
    else:
        metal = bulk(element, crystalstructure='fcc', a=lat_const, cubic=True)
        slab = surface(metal, hkl, size[2])

    # -- Add vacuum in the z direction
    slab.center(vacuum=vac_height, axis=2)

    # -- Center the slab at z=0 for RISM to work
    C = slab.cell[2][2]
    slab.translate([0., 0., -0.5*C])

    if square_cell:
        print('Attempting to square the cell')
        # -- Attempt to make a1 and a2 cell vectors similar in length
        cell = slab.get_cell()
        # -- Determine longer basis vector
        a1 = cell[0]
        a2 = cell[1]
        ratio = np.linalg.norm(a2) / np.linalg.norm(a1)
        if abs(ratio - 1) < prec:
            # lengths are equal
            print('cell is approximately square already')
            pass
        elif ratio > 1:
            # a2 > a1
            ratio = int(round(ratio))
            slab = slab * (ratio, 1, 1)
            print('ratio > 1: ratio = ', ratio)
        elif ratio < 1:
            # a2 < a1
            ratio = int(round(1 / ratio))
            slab = slab * (1, ratio, 1)
            print('ratio < 1: ratio = ', ratio)

    # -- Scale general hkl slab if requested
    if hkl not in [(1, 0, 0), (1, 1, 0), (1, 1, 1), (2, 1, 1)]:
        slab = slab * (size[0], size[1], 1)

    # -- sort atom positions
    pos = slab.get_positions()
    idx = np.lexsort([pos[:, 1], pos[:, 0], pos[:, 2]])
    slab.set_positions(pos[idx])

    # -- Extract an atom if requested for PMF of metal dissolution
    if extract_atom:
        # -- get the element index to displace
        n_surf_atoms = size[0] * size[1]
        tags = slab.get_tags()
        surf_idx = np.where(tags == 1)[0]
        shift_idx = n_surf_atoms // 2 + 2
        shift_idx = surf_idx[shift_idx]

        # -- Displace the surface atom
        pos = slab.get_positions()
        pos[shift_idx, 2] += z_shift

        # -- Update the cell vectors so slab bottom is vac_height Ang.
        #    from -C/2 and extracted atom is 10 Ang. from +C/2
        max_z = np.max(pos[:, 2])
        min_z = np.min(pos[:, 2])
        new_C = max_z - min_z + 2 * vac_height
        slab.cell[2, 2] = new_C

        # -- Shift the atoms so they have vac_height Ang. of vacuum
        #    padding
        shift = new_C / 2. - vac_height - max_z
        pos[:, 2] = pos[:, 2] + shift
        slab.positions = pos

        # -- fix the extracted atom
        extract_fix = FixAtoms(indices=[shift_idx])

    # == Adjust slab / cell if an adatom is present
    if adatom:
        add_adsorbate(slab, adsorbate, z_shift, ads_site)

    elif ad_molecule:
        mol = molecule(adsorbate)
        mol_index = 0
        if adsorbate == 'CO':
            mol_index = 1
        add_adsorbate(slab, mol, z_shift, ads_site, mol_index=mol_index)
        pass

    if adatom or ad_molecule or extract_atom:
        pos = slab.get_positions()

        # -- Update the cell vectors so slab bottom is vac_height Ang.
        #    from -C/2 and extracted atom is 10 Ang. from +C/2
        max_z = np.max(pos[:, 2])
        min_z = np.min(pos[:, 2])
        new_C = max_z - min_z + 2 * vac_height
        slab.cell[2, 2] = new_C

        # -- Shift the atoms so they have vac_height Ang. of vacuum
        #    padding
        shift = new_C / 2. - vac_height - max_z
        pos[:, 2] = pos[:, 2] + shift
        slab.positions = pos

    # -- Fix the last layer
    if hkl not in [(1, 0, 0), (1, 1, 0), (1, 1, 1)]:
        hkl = np.array(hkl)
        slab = detect_planes(slab, lat_const, hkl, fix_tag,
                             highlight_free=highlight_free)
#        fix = FixAtoms(indices=[atom.index for atom in slab if atom.position[2] < 0 ])
    indices = [atom.index for atom in slab if atom.tag > fix_tag]
    fix_bottom_layers = FixAtoms(indices=indices)

    if extract_atom:
        slab.set_constraint([fix_bottom_layers, extract_fix])
    else:
        fix_bottom_layers = FixAtoms(
            indices=[atom.index for atom in slab if atom.position[2] < 0])  # temporary
        slab.set_constraint(fix_bottom_layers)

    if return_lat_const:
        return slab, lat_const
    else:
        return slab


def get_pwscf_input_dict(slab, element, pseudo_lib, XC, max_seconds,
                         calculation='relax', cutoffs=None, bulk=False,
                         degauss=0.007, tot_charge=0):
    """
    Returns a python dictionary of PWscf input parameters.
    """

    nat = len(slab.get_chemical_symbols())
    conv_thr = nat * 2E-9
    pseudo_dir = get_pseudo_dir(pseudo_lib)
    if not cutoffs:
        if isinstance(element, str):
            cutoffs = get_ecuts(get_pseudo_name(element, lib_name=pseudo_lib))
        elif isinstance(element, list) or isinstance(element, tuple):
            for el in element:
                temp = get_ecuts(get_pseudo_name(el, lib_name=pseudo_lib))
                if cutoffs == None:
                    cutoffs = temp
                    continue
                if temp[0] > cutoffs[0]:
                    cutoffs = temp
                elif temp[0] == cutoffs[0]:
                    if temp[1] > cutoffs[1]:
                        cutoffs = temp

    input_dict = {'calculation': calculation,
                  'nstep': 1000,
                  'outdir': './SCRATCH',
                  'prefix': 'calc',
                  'etot_conv_thr': 1E-5,
                  'forc_conv_thr': 1E-4,
                  'pseudo_dir': pseudo_dir,
                  'max_seconds': max_seconds,
                  'ecutwfc': cutoffs[0],
                  'ecutrho': cutoffs[1],
                  'nosym': True,
                  'occupations': 'smearing',
                  'degauss': degauss,
                  'smearing': 'mv',
                  'assume_isolated': 'esm',
                  'esm_bc': 'bc1',
                  'conv_thr': float(f'{conv_thr:.3e}'),
                  'mixing_mode': 'local-TF'
                  }
    if XC != 'PBE':
        input_dict['input_dft'] = XC

    if tot_charge != 0:
        input_dict['tot_charge'] = tot_charge

    if bulk:
        del input_dict['assume_isolated']
        del input_dict['esm_bc']
        del input_dict['nosym']

    return input_dict


def get_rism_input_dict(slab, solvent_dict, explicit_lj_dict, n_debye_lengths=10,
                        expand_right=None, temperature=300., closure='kh',
                        rism3d_conv_thr=1.E-6, rism1d_conv_thr=1.E-8):
    """
    This function generates a RISM namelist and attempts to set reasonable values for the
    following parameters:

    * laue_starting_right
      - Based on the position of the second layer of slab in the DFT unit cell.

    * laue_expand_right
      - Based on the ionic strength and Debye length of the electrolyte, scaled by
        a user-specified parameter n_debye_lengths. For PMF calculations, this may
        be influenced by the position of the ion being pulled away from a surface.

    Notes:
    -----
    * We assume the dielectric constant of the solution to be 78 (pure water)
    * We assume that the extracted atom is renamed to have a suffix of 2 and the
      slab atoms to have a suffix of 1 (e.g., Al1, Al2, with Al2 shifted)
    """

    # -- Compute the ionic strength of the solution
    ionic_strength = 0.  # [ mol / L ]
    for k, v in solvent_dict.items():
        mol = k.lower()
        if mol == 'h2o':
            continue
        conc = v['conc']
        charge = k.count('+')
        if charge == 0:
            charge = -k.count('-')
        ionic_strength += conc*charge**2
    ionic_strength *= 0.5

    # -- Compute the Debye length of the solution
    eps_r = 78
    eps_0 = 8.854188E-12    # [ F / m   ]
    k_B = 1.380649E-23    # [ J / K   ]
    N_A = 6.02214076E23   # [ 1 / mol ]
    e = 1.602176634E-19  # [ C       ]
    kappa_inv = eps_r * eps_0 * k_B * temperature
    kappa_inv /= 2000 * N_A * e**2 * ionic_strength
    kappa_inv = np.sqrt(kappa_inv)
    # -- convert Debye length to Bohr
    kappa_inv *= 1E10        # Angstroms
    print(f"Debye length = {kappa_inv:.2f}")
    print(f"laue expand right = {n_debye_lengths*kappa_inv:.2f} Ang.")
    kappa_inv *= 1.889725989  # Bohr
    print(f"laue expand right = {n_debye_lengths*kappa_inv:.2f} Bohr")

    # determine parameters from slab and solvents info
    idx = np.where(slab.get_tags() == 2)[0][0]
    laue_starting_right = round(
        slab.positions[idx, 2] * 1.889725989, 1)  # in Bohr
    nsolv = len(solvent_dict)

    laue_expand_right = None
    if expand_right is not None:
        laue_expand_right = expand_right

    elif n_debye_lenghts is not None:
        laue_expand_right = round(n_debye_lengths * kappa_inv, 1)

    # -- populate the rism input dictionary
    rism_input_dict = {
        "closure": closure,
        "tempv": temperature,
        "laue_expand_right": laue_expand_right,
        "laue_starting_right": laue_starting_right,
        "nsolv": nsolv,
        "rism3d_conv_level": 0.5,
        "mdiis1d_step": 0.1,
        "mdiis3d_step": 0.5,
        "rism1d_conv_thr": rism1d_conv_thr,
        "rism3d_conv_thr": rism3d_conv_thr,
        "rism1d_maxstep": 100000,
        "rism3d_maxstep": 100000
    }

    # -- Update the dictionary to have explicit atom LJ parameters
    i = 1
    for sym, ff_params in explicit_lj_dict.items():
        print(sym, ff_params)
        if isinstance(ff_params, dict):
            key = f"solute_sigma({i})"
            val = round(ff_params['sigma'], 4)
            rism_input_dict[key] = val

            key = f"solute_epsilon({i})"
            val = round(ff_params['epsilon'], 4)
            rism_input_dict[key] = val

        elif isinstance(ff_params, str):
            key = f"solute_lj({i})"
            val = ff_params
            if not val.startswith("'") and not val.endswith("'"):
                val = "'" + val + "'"
            rism_input_dict[key] = val
        i += 1

    return rism_input_dict


def update_with_rism_input(rism_input_dict, solvent_dict, infile='pw.in',
                           starting_charge=0, extract_atom=False):
    """
    Since ASE doesn't natively support ESM-RISM input, read the file back in
    and update the input with ESM-RISM parameters.
    """

    with open(infile, 'r') as f:
        contents = f.readlines()

    updated_contents = ""
    do_control = False
    do_system = False
    do_electrons = False
    do_rism = False
    do_cell = False
    do_species = False
    do_labels = False

    for line in contents:
        # -- Update the control namelist
        if '&CONTROL' in line:
            do_control = True
        if do_control and line.strip() == '/':
            updated_contents += "   trism            = .TRUE.\n"
            updated_contents += "   lfcp             = .FALSE.\n"
            do_control = False

        # -- Update the system namelist
        if '&SYSTEM' in line:
            do_system = True

        if do_system and "ntyp" in line:
            if extract_atom:
                line = line.split()
                line[2] = str(int(line[2].strip()) + 1) + "\n"
                line[0] = " "*3 + line[0]
                line[1] = " "*12 + line[1]
                line = " ".join(line)
            if starting_charge != 0:
                line += f"   starting_charge(2) = {starting_charge}\n"
            do_system = False

        # -- Update the electrons namelist
        if '&ELECTRONS' in line:
            do_electrons = True
        if do_electrons and "mixing_mode" in line:
            updated_contents += "   diagonalization  = 'rmm-davidson'\n"
            updated_contents += "   diago_rmm_conv   = .false.\n"
            updated_contents += "   mixing_beta      = 0.2\n"
            updated_contents += "   electron_maxstep = 200\n"
            do_electrons = False

        # -- Update the RISM namelist
        if "&IONS" in line:
            do_rism = True
        if do_rism and "/" in line:
            updated_contents += "/\n"
            updated_contents += "&RISM\n"
            for k, v in rism_input_dict.items():
                if v in ["kh", "hnc"]:
                    v = "'" + v + "'"
                updated_contents += f"   {k} = {v}\n"
            do_rism = False

        # -- Skip the cell namelist in the input
        if "&CELL" in line:
            do_cell = True

        if do_cell and "/" in line:
            do_cell = False
            continue
        elif do_cell:
            continue

        # -- Update the atomic species
        if extract_atom:
            if "ATOMIC_SPECIES" in line:
                do_species = True
                i_spc = 0
            if do_species and "upf" in line.lower():
                if i_spc == 0:
                    line1 = line.split()
                    line1[0] = line1[0] + "1"
                    line1 = " ".join(line1)
                    line1 += "\n"

                    line2 = line.split()
                    line2[0] = line2[0] + "2"
                    line2 = " ".join(line2)
                    line2 += "\n"
                    updated_contents += line1 + line2
                    i_spc += 1
                    continue
                else:
                    line = line.split()
                    line[0] = line[0] + "1"
                    line = " ".join(line)
                    line += "\n"
                if "K_POINTS" in line:
                    do_species = False
                    continue

        # -- Update the atomic labels
            if "ATOMIC_POSITIONS" in line:
                do_labels = True
                updated_contents += line
                i_free = 0
                continue

            if do_labels:
                lsplit = line.split()
                if len(lsplit) == 1:
                    do_labels = False
                    continue

                if len(lsplit) == 4:
                    lsplit[0] += "1"
                    i_free += 1

                if len(lsplit) == 7 and i_free == 0:
                    lsplit[0] += "1"
                elif len(lsplit) == 7 and i_free > 0:
                    lsplit[0] += "2"

                line = " ".join(lsplit) + "\n"

        # -- By default, print the (possibly modified) line
        updated_contents += line

    # -- Add the solvents card
    updated_contents += "SOLVENTS { mol/L }\n"
    for species, ff_dict in solvent_dict.items():
        mol_name = ff_dict['ff_name']
        conc = ff_dict['conc']
        updated_contents += f"  {species:<8s} {conc:>6.3f}  {mol_name:<10s}\n"

    updated_contents += "\n"
    with open(infile, 'w') as f:
        f.write(updated_contents)

    return


def fix_k_grid(fname):
    """
    Reduces the number of k-points in the z direction to 1.
    """
    with open(fname, 'r') as f:
        lines = f.readlines()

    update_k_points = False
    with open(fname, 'w') as f:
        for line in lines:
            if "K_POINTS" in line:
                update_k_points = True
                f.write(line)
                continue
            elif update_k_points:
                line = line.split()
                line[2] = "1"
                line = " ".join(line)
                f.write(line + "\n")
                update_k_points = False
            else:
                f.write(line)


def write_input(slab, ksp, input_dict, pseudo_dict, fname='pw.in', bulk=False):
    """
    Writes the PWscf input file with custom parameters.
    """

    ksp = ksp / (2 * pi)
    write(fname, slab, format='espresso-in', input_data=input_dict,
          pseudopotentials=pseudo_dict, kspacing=ksp
          )
    if not bulk:
        fix_k_grid(fname)
    return

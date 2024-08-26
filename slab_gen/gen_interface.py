import argparse, os, sys, yaml
from slab_gen import Slab, InterfaceBuilder
from ase.visualize import view

#TODO: Remove override variable. Determine if solution is needed or feature to be scrapped
#TODO: Add all other parameters to YAML dump

def main():
    parse_args(sys.argv[1:])
    return

def parse_args(args: list) -> bool: 
    parser = argparse.ArgumentParser(
        prog='gen_interface.py',
        description="Search for optimal interface between two slabs")

    #-- Config file will specify all arguments by default BUT if you specify something here then you 
    parser.add_argument('-c', '--config',
                        default=None,
                        nargs='?',
                        help='Configuration used to populate candidate parameters.')
    
    parser.add_argument('-sc', '--save_config',
                        action='store_true',
                        help='Toggle to save current arguments to a config to the specified output location')
    
    parser.add_argument('-s1', '--structure1',
                        default=None,
                        nargs='?',
                        help='Location of structure used to create slab.')     #POSCAR, OUTCAR, LAMPS, CIF
    
    parser.add_argument('-s2', '--structure2',
                        default=None,
                        nargs='?',
                        help='Location of structures used to create slab.')
    
    parser.add_argument('-o', '--output',
                        default=None,
                        nargs='?',
                        help='Location of output directory.')
    
    parser.add_argument('-v', '--view',
                        action='store_true',
                        help='Toggle if you would like to view the unit cells of provided structures.')
    
    parser.add_argument('-ov', '--override',
                        action='store_true',
                        help='If you want to override arguments in the config use this flag' +
                        ' and follow with fields that you want overriden')
    
    parser.add_argument('--hkl1', nargs='+', default='001', help='Miller indices of surface 1')                                    
    parser.add_argument('--hkl2', nargs='+', default='001', help='Miller indices of surface 2')                                    
    parser.add_argument('--nlay',type=int, nargs='?', default=2, help='Number of layres in the generated slab')   
    parser.add_argument('--vac_height', default=10, type=float, nargs='?', help='Aangstrom height of the vaccuum')
    parser.add_argument('--ncand', default=10, type=int,nargs='?', help='Number of candidates to seek')
    parser.add_argument('--nmax', default=6, type=int, nargs='?', help='Constrains the size of the interface')
    parser.add_argument('--mmax', default=6, type=int, nargs='?', help='Constrains the size of the interface')
    parser.add_argument('--theta_min', default=0, type=float, nargs='?')
    parser.add_argument('--theta_max', default=180, type=float, nargs='?')
    parser.add_argument('--dtheta', default=4, type=float, nargs='?')
    parser.add_argument('--strain_cutoff', default=0.1, type=float, nargs='?')
    parser.add_argument('--theta_cutoff', default=20, type=float, nargs='?')
    parser.add_argument('--strain_weight', default=1, type=float, nargs='?')
    parser.add_argument('--area_weight', default=0.1, type=float, nargs='?')
    parser.add_argument('--angle_weight', default=1, type=float, nargs='?')
    parser.add_argument('--grid_spacing', default=1, type=float, nargs='?')
    parser.add_argument('--gen_trans', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--z_shift', default=2.5, type=float, nargs='?')

    args = parser.parse_args(args)

    print(args)

    #-- If YAML provided load it and update values
    if args.config is not None:
        with open(args.config, 'r') as cfg:
            config = yaml.safe_load(cfg)
    else:
        args.override = True
        print('--CONFIG NOT PROVIDED: USING DEFAULTS--')

    #-- Set structure paths and output path
    structure1 = args.structure1 if args.config is None else config.structure1
    structure2 = args.structure2 if args.config is None else config.structure2

    #-- Check structure paths to be valid
    if structure1 is None or structure2 is None:
        print('No provided structure paths')
        return False
    elif not os.path.exists(structure1) or not os.path.exists(structure2):
        print('Invalid structure paths')
        return False

    #-- Set slab generation variables
    hkl1 = tuple(int(i) for i in args.hkl1) if args.config is None else tuple(int(i) for i in config.hkl1)
    hkl2 = tuple(int(i) for i in args.hkl2) if args.config is None else tuple(int(i) for i in config.hkl2)
    nlay = args.nlay if args.config is None else config.nlay
    vac_height = args.vac_height if args.config is None else config.vac_height

    #-- Generate SLABS
    surface1 = Slab(structure1, hkl1)
    surface1.gen_slab(nlay, vac_height)

    surface2 = Slab(structure2, hkl2)
    surface2.gen_slab(nlay, vac_height)

    if args.view:
        view(surface1.slab)
        view(surface2.slab)

    #-- On success populate the rest of the arguments
    #-- TODO: Determine if args.config is None is necessary for this to work
    ncand = args.ncand if args.config is None or args.override else config.ncand 
    nmax = args.nmax if args.config is None or args.override else config.nmax
    mmax = args.mmax if args.config is None or args.override else config.mmax
    theta_min = args.theta_min if args.config is None or args.override else config.theta_min
    theta_max = args.theta_max if args.config is None or args.override else config.theta_max
    dtheta = args.dtheta if args.config is None or args.override else config.dtheta
    strain_cutoff = args.strain_cutoff if args.config is None or args.override else config.strain_cutoff
    theta_cutoff = args.theta_cutoff if args.config is None or args.override else config.theta_cutoff
    strain_weight = args.strain_weight if args.config is None or args.override else config.strain_weight
    area_weight = args.area_weight if args.config is None or args.override else config.area_weight
    angle_weight = args.angle_weight if args.config is None or args.override else config.angle_weight
    grid_spacing = args.grid_spacing if args.config is None or args.override else config.grid_spacing
    gen_trans = args.gen_trans if args.config is None or args.override else config.gen_trans
    z_shift = args.z_shift if args.config is None or args.override else config.z_shift
    verbose = args.verbose if args.config is None else config.verbose


    #-- Generate candidates
    #-- Requires modifying gen_candidates to respect a target working directory.
    #   - Possible to force working directory to target location before calling gen_candidates since
    #     it just puts it into whatever the current working directory of the program
        
    #-- Save config if flag is triggered
    if args.save_config:
        print('Writing YAML')
        writeYAML(args)

    ib = InterfaceBuilder(surface1, surface2, verbose)

    if os.path.exists(args.output):
        os.chdir(args.output)

    ib.gen_candidates(ncand,
                        nmax,
                        mmax,
                        theta_min,
                        theta_max,
                        dtheta,
                        strain_cutoff,
                        theta_cutoff,
                        strain_weight,
                        area_weight,
                        angle_weight,
                        grid_spacing,
                        gen_trans,
                        z_shift)
    ib._list_top_candidates(ncand)
    return True

#TODO: Ensure 
def writeYAML(args: list) -> None:
    """Function to take parsed arguments and write to YAML"""   
    to_write = {
        'config': args.config if args.config is not None else "",
        'structure1': args.structure1,
        'structure2': args.structure2,
        'output': args.output if args.output is not None else "",
        'ncand': str(args.ncand),
        'hkl1' : str(args.hkl1),
        'hkl2' : str(args.hkl2),
        'nlay': str(args.nlay),
        'vac_height' : str(args.vac_height),
        'nmax': str(args.nmax),
        'mmax': str(args.mmax),
        'theta_min': str(args.theta_min),
        'theta_max': str(args.theta_max),
        'dtheta': str(args.dtheta),
        'strain_cutoff': str(args.strain_cutoff),
        'theta_cutoff': str(args.theta_cutoff),
        'strain_weight': str(args.strain_weight),
        'area_weight': str(args.area_weight),
        'angle_weight': str(args.angle_weight),
        'grid_spacing': str(args.grid_spacing),
        'gen_trans': str(args.gen_trans),
        'z_shift': str(args.z_shift),
        'view' : str(args.view),
        'verbose': str(args.verbose)
    }

    # data = yaml.safe_load(to_write)
    hkl1 = ''.join(args.hkl1)
    hkl2 = ''.join(args.hkl2)
    struct1 = os.path.basename(args.structure1).split('.')[0] + f'({hkl1})'
    struct2 = os.path.basename(args.structure2).split('.')[0] + f'({hkl2})'

    config_name = f'{struct1}_{struct2}.yaml'
    if os.path.exists(args.output):
        os.chdir(args.output)
    with open(config_name, 'w') as config:
        yaml.dump(to_write, config)
        
if __name__ == "__main__":
    main()

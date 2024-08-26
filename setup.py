from setuptools import setup, find_packages

setup(name='slab_gen',
      version='0.1',
      description='Utilities for generating slab surface models with ASE.',
#      url='https://github.com/sweitzner/cube_viskit',
      author='Stephen Weitzner',
      author_email='weitzner1@llnl.gov',
      license='MIT',
      packages=['slab_gen'],
#      packages=find_packages(),
#      package_data={'my_data_pack': ['my_data/*']}
      include_package_data=True,
#      scripts=['bin/convert_cube.py', 
#               'bin/extract_pw_coord.py', 
#               'bin/extract_xyz_traj.py'
#               'bin/planar_average.py'],
      install_requires=['numpy', 'scipy', 'ase', 'spglib'],
      zip_safe=False)


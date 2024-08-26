import pytest
from ase import Atoms
from ase.build import bulk, fcc211
from slab_gen import Graph


@pytest.fixture
def crystal() -> Atoms:
    return bulk('LiF', crystalstructure='rocksalt', a=4.08)


@pytest.fixture
def surface() -> Atoms:
    return fcc211('Al', (3, 3, 3), a=4.05, vacuum=10.)


@pytest.fixture
def molecule() -> Atoms:
    return Atoms('H2', positions=[(0, 0, 0), (0, 0, 1.2)])


def test_graph_crystal(crystal: Atoms) -> None:
    g = Graph(crystal)
    assert g.atoms_type == 'crystal'


def test_graph_surface(surface: Atoms) -> None:
    g = Graph(surface)
    assert g.atoms_type == 'surface'


def test_graph_molecule(molecule: Atoms) -> None:
    g = Graph(molecule)
    assert g.atoms_type == 'molecule'

def test_graph_coord_num_crystal(crystal: Atoms) -> None:
    g = Graph(crystal)
    cn = g.coord_num
    print(cn)
import sys
from os.path import *
import numpy as np
#from simtk.openmm.app.element import Element
sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
from AmberFD import DispersionPauli, VectorI, VectorD, VectorPairII, PairII
import time
from AssertEqual import AssertEqual

int_to_exp_frz = {1: 2.50, 6: 2.45, 7: 2.17, 8: 2.42}
int_to_exp_del = {1: 2.10, 6: 2.05, 7: 1.80, 8: 1.90}

def get_bonds(coords, atoms):
    bonds = VectorPairII()
    for i, atom_i in enumerate(atoms):
        for j, atom_j in enumerate(atoms):
            if j <= i : continue
            dist = np.linalg.norm(coords[i] - coords[j])

            if atom_i == 'H' or atom_j == 'H':
                check_dist = 1.2
            else:
                check_dist = 1.6
            if dist <= check_dist:
                pair = PairII(i, j)
                bonds.push_back(pair)
    return bonds


atom_to_nuc = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
if __name__ == "__main__":
    data = np.loadtxt('data/u_u_data_1.txt', dtype=object).T
    atom_names = data[0]
    nuclei = data[1].astype('int32')
    elms = data[2]
    x, y, z, frz_chg, exp_frz, exp_dyn, pauli_exp, pauli_radii =  data[3:].astype(float)
    coords = np.array([x, y, z]).T
    bonds = get_bonds(coords, elms)
    n_atoms = int(len(coords)/2)

    disp = DispersionPauli(nuclei, pauli_exp, pauli_radii)
    disp.create_exclusions_from_fragment(np.arange(0, n_atoms))
    disp.create_exclusions_from_fragment(np.arange(n_atoms, n_atoms*2))

    energy = disp.calc_energy(coords.flatten())
    pauli_energy = disp.get_pauli_energy()*2625.5009
    disp_energy = disp.get_disp_energy()*2625.5009

    AssertEqual(pauli_energy, 0.4749888820543172741, 1e-14)
    AssertEqual(disp_energy, -8.4318942845087168081, 1e-14)
    


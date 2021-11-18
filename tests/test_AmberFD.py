import sys
from os.path import *
import numpy as np
#from simtk.openmm.app.element import Element
sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
from AmberFD import FlucDens, VectorI, VectorD, VectorPairII
import time

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
    #   import and assign parameters
    xyz_data = np.loadtxt('data/u_u.xyz', skiprows=2, dtype=str)
    frz_chg = np.loadtxt('data/u_u.chg')
    n_atoms = int(len(xyz_data)/2)
    coords = xyz_data[:, [1,2,3]].astype(float)*1.88973
    atoms = xyz_data[:, 0]
    nuclei = np.array([atom_to_nuc[x] for x in atoms])
    exp_frz = np.array([int_to_exp_frz[x] for x in nuclei])
    exp_dyn = exp_frz*0.5
    exp_dyn[n_atoms:] += 0.25
    bonds = get_bonds(coords, atoms)

    start = time.time()
    for i in range(1):
        #   total system
        fluc = FlucDens(frz_chg, nuclei, exp_frz, exp_dyn)
        fluc.create_frz_exclusions_from_bonds(bonds, 3)
        fluc.add_fragment(VectorI(range(n_atoms)))
        fluc.add_fragment(VectorI(range(n_atoms, n_atoms*2)))
        total_fluc = fluc.calc_energy(coords.flatten())
        frz_energy = fluc.get_frozen_energy()

        #   individual fragments
        fluc_1 = FlucDens(frz_chg[:n_atoms], nuclei[:n_atoms], exp_frz[:n_atoms], exp_dyn[:n_atoms])
        fluc_2 = FlucDens(frz_chg[n_atoms:], nuclei[n_atoms:], exp_frz[n_atoms:], exp_dyn[n_atoms:])
        bonds_1 = get_bonds(coords[:n_atoms], atoms[:n_atoms])
        bonds_2 = get_bonds(coords[n_atoms:], atoms[n_atoms:])
        fluc_1.create_frz_exclusions_from_bonds(bonds_1, 3)
        fluc_2.create_frz_exclusions_from_bonds(bonds_2, 3)
        fluc_1.add_fragment(VectorI(range(n_atoms)))
        fluc_2.add_fragment(VectorI(range(n_atoms)))
        frz_energy_1 = fluc_1.calc_energy(VectorD(coords[:n_atoms].flatten()))
        frz_energy_2 = fluc_2.calc_energy(VectorD(coords[n_atoms:].flatten()))


        eng_diff = frz_energy - frz_energy_1 - frz_energy_2
        print("Frz Energy_12: {:15.8f}".format(frz_energy))
        print("Frz Energy_1:  {:15.8f}".format(frz_energy_1))
        print("Frz Energy_2:  {:15.8f}".format(frz_energy_2))
        print("Diff:          {:15.8f}".format(eng_diff))
        print("Diff (kJ/mol): {:15.8f}".format(eng_diff*2625.5009))
        print("Pol. (kJ/mol): {:15.8f}".format((total_fluc - frz_energy)*2625.5009))
        -84.19430925179343
        assert(abs(eng_diff*2625.5009 - -84.19430925179343) < 1e-13)
    end = time.time()

    #print(end - start)

    
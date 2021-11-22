import sys
from os.path import *
import numpy as np
#from simtk.openmm.app.element import Element
sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
from AmberFD import FlucDens, VectorI, VectorD, VectorPairII, PairII
import time
from AssertEqual import AssertEqual

int_to_exp_frz = {1: 2.50, 6: 2.45, 7: 2.17, 8: 2.42}
int_to_exp_dyn = {1: 2.10, 6: 2.05, 7: 1.80, 8: 1.90}

def get_bonds(coords, atoms):
    bonds = VectorPairII()
    for i, atom_i in enumerate(atoms):
        for j, atom_j in enumerate(atoms):
            if j <= i : continue
            dist = np.linalg.norm(coords[i] - coords[j])

            if atom_i == 'H' or atom_j == 'H':
                check_dist = 1.2*1.8897259886
            else:
                check_dist = 1.6*1.8897259886
            if dist <= check_dist:
                pair = PairII(i, j)
                bonds.push_back(pair)
    return bonds


atom_to_nuc = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
ANG2BOHR = 1.8897259886
AU_2_KJ_PER_MOL = 2625.5009
if __name__ == "__main__":
    data = np.loadtxt('data/u_u_data.txt', dtype=object).T
    atom_names = data[0]
    nuclei = data[1].astype(int)
    elms = data[2]
    x, y, z, frz_chg, exp_frz, exp_dyn, pauli_exp, pauli_radii =  data[3:].astype(float)
    coords = np.array([x, y, z]).T
    bonds = get_bonds(coords, elms)
    n_atoms = int(len(coords)/2)


    start = time.time()
    for i in range(1):
        #   total system
        fluc = FlucDens(frz_chg, nuclei, exp_frz, exp_dyn)
        fluc.create_frz_exclusions_from_bonds(bonds, 3)
        fluc.add_fragment(VectorI(range(n_atoms)))
        fluc.add_fragment(VectorI(range(n_atoms, n_atoms*2)))
        fluc.set_dampening(1.5467, 1.4364)
        total_fluc = fluc.calc_energy(coords.flatten())*AU_2_KJ_PER_MOL
        frz_energy = fluc.get_frozen_energy()*AU_2_KJ_PER_MOL
        

        #   individual fragments
        fluc_1 = FlucDens(frz_chg[:n_atoms], nuclei[:n_atoms], exp_frz[:n_atoms], exp_dyn[:n_atoms])
        fluc_2 = FlucDens(frz_chg[n_atoms:], nuclei[n_atoms:], exp_frz[n_atoms:], exp_dyn[n_atoms:])
        bonds_1 = get_bonds(coords[:n_atoms], elms[:n_atoms])
        bonds_2 = get_bonds(coords[n_atoms:], elms[n_atoms:])
        fluc_1.create_frz_exclusions_from_bonds(bonds_1, 3)
        fluc_2.create_frz_exclusions_from_bonds(bonds_2, 3)
        fluc_1.add_fragment(VectorI(range(n_atoms)))
        fluc_2.add_fragment(VectorI(range(n_atoms)))
        fluc_1.set_dampening(1.5467, 1.4364)
        fluc_2.set_dampening(1.5467, 1.4364)
        frz_energy_1 = fluc_1.calc_energy(VectorD(coords[:n_atoms].flatten()))*AU_2_KJ_PER_MOL
        frz_energy_2 = fluc_2.calc_energy(VectorD(coords[n_atoms:].flatten()))*AU_2_KJ_PER_MOL


        eng_diff = frz_energy - frz_energy_1 - frz_energy_2
        pol_energy = total_fluc - frz_energy
        print("Frz Energy_12: {:15.16f}".format(frz_energy))
        print("Frz Energy_1:  {:15.16f}".format(frz_energy_1))
        print("Frz Energy_2:  {:15.16f}".format(frz_energy_2))
        print("Diff:          {:15.16f}".format(eng_diff))
        print("Pol. (kJ/mol): {:15.16f}".format((total_fluc - frz_energy)))

        AssertEqual(frz_energy,  -153.9164545735984291, 1e-14)
        AssertEqual(frz_energy_1, -57.9401318692427196, 1e-14)
        AssertEqual(frz_energy_2, -58.5616039140674047, 1e-14)
        AssertEqual(eng_diff,     -37.4147187902883047, 1e-14)
        AssertEqual(pol_energy,   -10.0445370280187092, 1e-14)

        

    end = time.time()

    #print(end - start)

    
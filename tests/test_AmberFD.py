import sys
from os.path import *
import numpy as np
#from simtk.openmm.app.element import Element
sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
from AmberFD import AmberFD, FlucDens, VectorI, VectorD, VectorPairII, ParticleInfo, PairII
import time
from AssertEqual import AssertEqual

ANG2BOHR = 1.8897259886
AU_2_KJ_PER_MOL = 2625.5009

def get_bonds(coords, atoms):
    bonds = VectorPairII()
    for i, atom_i in enumerate(atoms):
        for j, atom_j in enumerate(atoms):
            if j <= i : continue
            dist = np.linalg.norm(coords[i] - coords[j])

            if atom_i == 'H' or atom_j == 'H':
                check_dist = 1.2*ANG2BOHR
            else:
                check_dist = 1.6*ANG2BOHR
            if dist <= check_dist:
                pair = PairII(i, j)
                bonds.push_back(pair)
    return bonds


atom_to_nuc = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
if __name__ == "__main__":
    #   import and assign parameters
    data = np.loadtxt('data/u_u_data.txt', dtype=object).T
    atom_names = data[0]
    nuclei = data[1].astype('int32')
    elms = data[2]
    x, y, z, frz_chg, exp_frz, exp_dyn, pauli_exp, pauli_radii =  data[3:].astype(float)
    coords = np.array([x, y, z]).T
    bonds = get_bonds(coords, elms)
    n_atoms = int(len(coords)/2)
    

    amber = AmberFD(len(coords))
    for n, nuc in enumerate(nuclei):
        particle = ParticleInfo(nuc)
        particle.dyn_exp = exp_dyn[n]
        particle.frz_exp = exp_frz[n]
        particle.frz_chg = frz_chg[n]
        particle.pauli_exp = pauli_exp[n]
        particle.pauli_radii = pauli_radii[n]
        amber.add_particle(particle)

    #   create dispersion-pauli force
    disp = amber.create_disp_pauli_force()
    disp.create_exclusions_from_fragment(np.arange(0, n_atoms))
    disp.create_exclusions_from_fragment(np.arange(n_atoms, n_atoms*2))

    #   create electrostatics force
    fluc = amber.create_fluc_dens_force()
    fluc.create_frz_exclusions_from_bonds(bonds, 3)
    fluc.add_fragment(np.arange(0, n_atoms))
    fluc.add_fragment(np.arange(n_atoms, n_atoms*2))
    fluc.set_dampening(1.5467, 1.4364)
    
    energies = amber.calc_energy_forces(coords.flatten())

    disp.calc_energy(coords.flatten())
    fluc.calc_energy(coords.flatten())

    AssertEqual(energies.E_frz*AU_2_KJ_PER_MOL,   -153.9164545735984291,    1e-14)
    AssertEqual(energies.E_pol*AU_2_KJ_PER_MOL,    -10.0445370280187092,    1e-14)
    AssertEqual(energies.E_pauli*AU_2_KJ_PER_MOL,    0.4749888820543172741, 1e-14)
    AssertEqual(energies.E_disp*AU_2_KJ_PER_MOL,    -8.4318942845087168081, 1e-14)
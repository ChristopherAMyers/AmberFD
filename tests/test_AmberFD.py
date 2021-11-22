import sys
from os.path import *
import numpy as np
#from simtk.openmm.app.element import Element
sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
from AmberFD import AmberFD, FlucDens, VectorI, VectorD, VectorPairII, ParticleInfo, PairII
import time

int_to_exp_frz = {1: 2.50, 6: 2.45, 7: 2.17, 8: 2.42}
int_to_exp_dyn = {1: 2.10, 6: 2.05, 7: 1.80, 8: 1.90}
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
    xyz_data = np.loadtxt('data/u_u.xyz', skiprows=2, dtype=str)
    frz_chg = np.loadtxt('data/u_u.chg')
    param_data = np.loadtxt('data/u_disp_pauli_params.txt', dtype=str)
    param_data = param_data.T

    #   pdb info is used for dispersion and pauli energies
    pdb_data = np.loadtxt('tmp/u_u.pdb', dtype=str)
    atom_names = pdb_data[:, 2]
    frz_exp_dict = dict(zip(param_data[0], param_data[1].astype(float)))
    frz_dyn_dict = dict(zip(param_data[0], param_data[2].astype(float)))
    frz_chg_dict = dict(zip(param_data[0], param_data[3].astype(float)))
    pauli_exp_dict = dict(zip(param_data[0], param_data[4].astype(float)))
    pauli_radii_dict = dict(zip(param_data[0], param_data[5].astype(float)))
    

    n_atoms = int(len(xyz_data)/2)
    #coords = xyz_data[:, [1,2,3]].astype(float)*ANG2BOHR
    #atoms = xyz_data[:, 0]
    coords = pdb_data[:, [6,7,8]].astype(float)*ANG2BOHR
    atoms = pdb_data[:, -1]
    nuclei = np.array([atom_to_nuc[x] for x in atoms], dtype='intc')
    exp_frz = np.array([int_to_exp_frz[x] for x in nuclei])
    exp_dyn = np.array([int_to_exp_dyn[x] for x in nuclei])
    #exp_dyn = exp_frz*0.5
    exp_dyn[n_atoms:] += 0.25
    bonds = get_bonds(coords, atoms)
    pauli_exp = [pauli_exp_dict[x] for x in atom_names]
    pauli_exp = [pauli_exp_dict[x] for x in atom_names]
    pauli_radii = [pauli_radii_dict[x] for x in atom_names]

    amber = AmberFD(len(atoms))
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

    for n, name in enumerate(atom_names):
        mol_str = ' {:4s} {:1d} {:8.3f}  {:8.3f}  {:8.3f}'.format(name, nuclei[n], *tuple(coords[n]))
        param_str = ('  {:8.3f}'*5).format(frz_chg[n], exp_frz[n], exp_dyn[n], pauli_exp[n], pauli_radii[n])
        print(mol_str + param_str)
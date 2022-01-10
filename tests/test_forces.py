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

    for data_set in (1, 2):
        print("     DATA SET %d" % data_set)
        #   import and assign parameters
        data = np.loadtxt('data/u_u_data_%d.txt' % data_set, dtype=object).T
        goal_fluc_forces = np.loadtxt('data/u_u_fluc_forces_%d.txt' % data_set)
        goal_disp_forces = np.loadtxt('data/u_u_disp_forces_%d.txt' % data_set)
        goal_total_forces = np.loadtxt('data/u_u_total_forces_%d.txt' % data_set)
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
            amber.add_particle(n, particle)

        #   create dispersion-pauli force
        disp = amber.create_disp_pauli_force()
        disp.set_use_secondary_radii(True)
        disp.create_exclusions_from_fragment(np.arange(0, n_atoms))
        disp.create_exclusions_from_fragment(np.arange(n_atoms, n_atoms*2))

        #   create electrostatics force
        fluc = amber.create_fluc_dens_force()
        fluc.create_frz_exclusions_from_bonds(bonds, 3)
        fluc.add_fragment(np.arange(0, n_atoms))
        fluc.add_fragment(np.arange(n_atoms, n_atoms*2))
        fluc.set_dampening(1.5467, 1.4364)

        #   calculate forces and energies
        disp.calc_energy(coords.flatten())
        fluc.calc_energy(coords.flatten())
        energies = amber.calc_energy_forces(coords.flatten())
        fluc_forces = np.array(fluc.get_forces())*2625.5009*ANG2BOHR
        disp_forces = np.array(disp.get_forces())*2625.5009*ANG2BOHR
        forces = np.array(amber.get_forces())*2625.5009*ANG2BOHR
        # np.savetxt('data/u_u_fluc_forces_%d.txt' % data_set, fluc_forces)
        # np.savetxt('data/u_u_disp_forces_%d.txt' % data_set, disp_forces)
        # np.savetxt('data/u_u_total_forces_%d.txt' % data_set, forces)

        
        #   numerical derivatives
        if True:
            eps = (np.finfo(np.float64).eps)**(1/3)
            for n, coord in enumerate(coords):
                numerical_force = np.zeros(3)
                numerical_disp_force = np.zeros(3)
                numerical_fluc_force = np.zeros(3)
                for x in [0, 1, 2]:
                    new_coords_p = coords.copy()
                    new_coords_m = coords.copy()
                    new_coords_p[n][x] += eps
                    new_coords_m[n][x] -= eps

                    energy_p = amber.calc_energy_forces(new_coords_p.flatten()).total()
                    energy_m = amber.calc_energy_forces(new_coords_m.flatten()).total()
                    energy_disp_p = disp.calc_energy(new_coords_p.flatten())
                    energy_disp_m = disp.calc_energy(new_coords_m.flatten())
                    energy_fluc_p = fluc.calc_energy(new_coords_p.flatten(), True, True)
                    energy_fluc_m = fluc.calc_energy(new_coords_m.flatten(), True, True)

                    numerical_force[x] = -(energy_p - energy_m)/(2*eps)*2625.5009*ANG2BOHR
                    numerical_disp_force[x] = -(energy_disp_p - energy_disp_m)/(2*eps)*2625.5009*ANG2BOHR
                    numerical_fluc_force[x] = -(energy_fluc_p - energy_fluc_m)/(2*eps)*2625.5009*ANG2BOHR

                print("Data set: {:d} force {:d}".format(data_set, n))
                np.testing.assert_allclose(forces[n],      numerical_force,      rtol=eps, atol=eps)
                np.testing.assert_allclose(disp_forces[n], numerical_disp_force, rtol=eps, atol=eps)
                np.testing.assert_allclose(fluc_forces[n], numerical_fluc_force, rtol=eps, atol=eps)
                
                # print(('{:15.12f} '*3 + ' | ' + '{:15.12f} '*3).format(*tuple(forces[n]), *tuple(numerical_force)))
                # print(('{:15.12f} '*3 + ' | ' + '{:15.12f} '*3).format(*tuple(disp_forces[n]), *tuple(numerical_disp_force)))
                #print(('{:15.12f} '*3 + ' | ' + '{:15.12f} '*3).format(*tuple(fluc_forces[n]), *tuple(numerical_fluc_force)))

                #print(('{:15.12f} '*3 + ' | ' + '{:15.12f} '*3).format(*tuple(forces[n]), *tuple(goal_total_forces[n])))
                #print(('{:15.12f} '*3 + ' | ' + '{:15.12f} '*3).format(*tuple(disp_forces[n]), *tuple(goal_disp_forces[n])))
                #print(('{:15.12f} '*3 + ' | ' + '{:15.12f} '*3).format(*tuple(fluc_forces[n]), *tuple(goal_fluc_forces[n])))
        
        
        # np.testing.assert_allclose(fluc_forces, goal_fluc_forces, atol=1e-12)
        # np.testing.assert_allclose(disp_forces, goal_disp_forces, atol=1e-12)
        # np.testing.assert_allclose(forces, goal_total_forces, atol=1e-12)
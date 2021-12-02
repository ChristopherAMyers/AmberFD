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

fluc_forces = np.array([
        [ 0.079750489342, -0.069744915447, -0.000000000001],
        [-0.075088793915,  0.030804851820,  0.000000000001],
        [-0.033267622919,  0.040494701632, -0.004912226011],
        [-0.033267622919,  0.040494701632,  0.004912226011],
        [ 1.269335471289, -0.028395335405,  0.000000000001],
        [-1.179613229037,  0.742241057582, -0.000000000001],
        [-0.094227161616,  0.135536665721, -0.000000000001],
        [-0.017184554324, -0.334609858714,  0.000000000001],
        [-0.257065660607, -0.371138257698,  0.000000000001],
        [ 1.206937042877,  2.035235001071, -0.000000000001],
        [ 1.321105599085, -2.000579736134,  0.000000000001],
        [-4.048609173259, -0.842483134260,  0.000000000001],
        [ 0.244732307170,  0.126262528018, -0.000000000001],
        [-7.497943939252,  5.750258159352, -0.000000000001],
        [ 0.027387697021,  0.153833284046, -0.000000000001],
        [-0.000351795945, -0.023831834407, -0.000000000001],
        [ 0.101716504696,  0.004278524065,  0.000000000001],
        [-0.249922542197, -0.056428826959,  0.000000000001],
        [-5.535463723347,  0.734408429604,  0.000000000001],
        [ 0.154732695821,  0.068616351440, -0.000000000001],
        [-0.183514868450,  0.001453171279, -0.000000000001],
        [-0.075969782539, -0.029374454493, -0.013275691411],
        [-0.075969782539, -0.029374454493,  0.013275691411],
        [ 3.403750700478,  1.043907803441,  0.000000000001],
        [-0.808559023151, -1.887774853761, -0.000000000001],
        [-0.098059766176, -0.182480057740, -0.000000000001],
        [-0.028070457879,  0.314631023502,  0.000000000001],
        [-0.208658698755,  0.244606484063,  0.000000000001],
        [ 1.380288629563, -0.699256020886, -0.000000000001],
        [-1.495994785066,  1.890420944506,  0.000000000001],
        [-1.224194734664, -0.070734166734,  0.000000000001],
        [ 0.116917443719, -0.100399744319, -0.000000000001],
        [ 8.278099337768, -5.830128245555, -0.000000000001],
        [-0.019850972733, -0.144780109662,  0.000000000001],
        [ 5.377191309807, -0.653547776751,  0.000000000001],
        [ 0.280125019394,  0.048756901223,  0.000000000001],
        [-0.089745950066, -0.010343471416,  0.000000000001],
        [ 0.088524393324, -0.040835329163, -0.000000000001]
])

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
    goal_forces = np.loadtxt('data/u_u_forces.txt')
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

    forces = np.array(fluc.get_forces())*2625.5009*ANG2BOHR
    np.testing.assert_allclose(forces, goal_forces, atol=1e-12)
    exit()

    eps = 1e-5
    for n, coord in enumerate(coords):
        numerical_force = np.zeros(3)
        for x in [0, 1, 2]:
            new_coords_p = coords.copy()
            new_coords_m = coords.copy()
            new_coords_p[n][x] += eps
            new_coords_m[n][x] -= eps

            frz_energy_p = fluc.calc_energy(new_coords_p.flatten(), False, True)
            frz_energy_m = fluc.calc_energy(new_coords_m.flatten(), False, True)

            numerical_force[x] = -(frz_energy_p - frz_energy_m)/(2*eps)*2625.5009*ANG2BOHR

        diff = (forces[n] - numerical_force)/numerical_force
        #print(('{:15.12f} '*3).format(*tuple(diff)))
        print(('{:15.12f} '*3 + ' | ' + '{:15.12f} '*3).format(*tuple(forces[n]), *tuple(numerical_force)))


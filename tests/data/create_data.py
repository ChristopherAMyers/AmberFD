from os.path import *
import numpy as np

ANG2BOHR = 1.8897259886

def get_bonds(coords, atoms):
    bonds = []
    for i, atom_i in enumerate(atoms):
        for j, atom_j in enumerate(atoms):
            if j <= i : continue
            dist = np.linalg.norm(coords[i] - coords[j])

            if atom_i == 'H' or atom_j == 'H':
                check_dist = 1.2*ANG2BOHR
            else:
                check_dist = 1.6*ANG2BOHR
            if dist <= check_dist:
                pair = (i, j)
                bonds.append(pair)
    return bonds


atom_to_nuc = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
if __name__ == "__main__":
    #   import and assign parameters
    param_data = np.loadtxt('u_disp_pauli_params.txt', dtype=str)
    param_data = param_data.T
    frz_exp_dict = dict(zip(param_data[0], param_data[1].astype(float)))
    dyn_exp_dict = dict(zip(param_data[0], param_data[2].astype(float)))
    frz_chg_dict = dict(zip(param_data[0], param_data[3].astype(float)))
    pauli_exp_dict = dict(zip(param_data[0], param_data[4].astype(float)))
    pauli_radii_dict = dict(zip(param_data[0], param_data[5].astype(float)))


    pdb_data = np.loadtxt('u_u_2.pdb', dtype=str)
    atom_names = pdb_data[:, 2]
    frz_exp = [frz_exp_dict[x] for x in atom_names]
    dyn_exp = [dyn_exp_dict[x] for x in atom_names]
    frz_chg = [frz_chg_dict[x] for x in atom_names]
    pauli_exp = [pauli_exp_dict[x] for x in atom_names]
    pauli_radii = [pauli_radii_dict[x] for x in atom_names]
    

    n_atoms = int(len(pdb_data)/2)
    coords = pdb_data[:, [6,7,8]].astype(float)*ANG2BOHR
    atoms = pdb_data[:, -1]
    nuclei = np.array([atom_to_nuc.get(x, 0) for x in atoms], dtype='intc')
    bonds = get_bonds(coords, atoms)
    
    print('# name nuclei elm   x(au)   y(au)      z(au)   frz_chg    exp_frz   exp_dyn  pauli_exp  pauli_radii')
    for n, name in enumerate(atom_names):
        mol_str = ' {:4s} {:4d} {:>4s} {:8.3f}  {:8.3f}  {:8.3f}'.format(name, nuclei[n], atoms[n], *tuple(coords[n]))
        param_str = ('  {:8.3f}'*5).format(frz_chg[n], frz_exp[n], dyn_exp[n], pauli_exp[n], pauli_radii[n])
        print(mol_str + param_str)
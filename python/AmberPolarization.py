from copy import deepcopy, copy
from openmm.app import forcefield as ff, Simulation, Element
from openmm.openmm import CustomExternalForce, NonbondedForce
import openmm.unit as uu
import sys
from os.path import *
import numpy as np
from openmm.app.gromacsgrofile import GromacsGroFile
from openmm.app.gromacstopfile import GromacsTopFile
from openmm.app.internal.pdbstructure import Model
from openmm.app.modeller import Modeller
from openmm.app.pdbfile import PDBFile
#sys.path.insert(1, '/network/rit/lab/ChenRNALab/awesomeSauce/code/fluctuating_density/build')
sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
from AmberFD import AmberFD, FlucDens, VectorI, VectorD, VectorPairII, ParticleInfo, MapID

class AmberFDGenerator(object):

    def __init__(self, forcefield) -> None:
        self.dispersion_params = {}
        self.dispersion_radii = {}
        self.dispersion_C6 = {}
        self.polPairwise_atoms = []
        self.fragments = {}
        self.res_to_fragment = {}
        self.damp_coeff = 0.0
        self.damp_exp = 1.0

        self.residues = {}
        
    def registerDispersion(self, element):
        for param in ['s6', 'a1', 'a2']:
            self.dispersion_params[param] = float(element.attrib[param])
        for atom in element.findall('Atom'):
            atom_type = atom.attrib['type']
            nuclei = Element.getBySymbol(atom_type).atomic_number
            self.dispersion_C6[nuclei]    = float(atom.attrib['C6'])
            self.dispersion_radii[nuclei] = float(atom.attrib['vdw-radii'])


    def registerPolPairwise(self, element):
        for atom in element.findall('Atom'):
            pol_atom = {}
            pol_atom['type'] = atom.attrib['type']
            for attrib in ['exp', 'radii', 'type']:
                pol_atom[attrib] = atom.attrib[attrib]
            self.polPairwise_atoms.append(pol_atom)

    def add_atom(self, atom_element):
        res_atom = {}
        for attrib in ['frz_chg', 'frz_exp', 'dyn_exp', 'pauli_exp', 'pauli_radii']:
            res_atom[attrib] = float(atom_element.attrib[attrib])
        return res_atom

    @staticmethod
    def parseElement(element, ff):
        print("IN PARSE ELEMENT")
        existing = [f for f in ff._forces if isinstance(f, AmberFDGenerator)]
        if len(existing) == 0:
            generator = AmberFDGenerator(ff)
            ff.registerGenerator(generator)
        else:
            generator = existing[0]

        generator.xml = element

        #   parse dispersion parameters
        disp_sections = element.findall('Dispersion')
        if len(disp_sections) > 1:
            raise ValueError("More than one <Dispersion> section specified in force field")
        if len(disp_sections) == 0:
            raise ValueError(" <Dispersion> section not found in <AmberFDForce>")
        generator.registerDispersion(disp_sections[0])

        #   pairwise polarization parameters
        pol_pairwise_sections = element.findall('PolarizationPairwise')
        if len(pol_pairwise_sections) > 1:
            raise ValueError("More than one <PolarizationPairwise> section specified in force field")
        generator.registerPolPairwise(pol_pairwise_sections[0])

        generator.damp_coeff = element.findall('PolarizationNonPairwise')[0].attrib['damp_coeff']
        generator.damp_coeff = element.findall('PolarizationNonPairwise')[0].attrib['damp_coeff']


    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        print("IN CREATE FORCE")
        force = AmberFD()
        pol_residues = set()
        pauli_residues = set()
        atom_indices = {}
        for atom, params in data.atomParameters.items():
            add_atom = False
            nuclei = atom.element.atomic_number if isinstance(atom.element, Element) else 0
            particle = ParticleInfo(nuclei)
            polar_params_present = [x in params for x in ['frz_chg', 'frz_exp', 'dyn_exp']]
            pauli_params_present = [x in params for x in ['pauli_exp', 'pauli_radii']]
            
            #   fluctuating density force
            if any(polar_params_present):
                if not all(polar_params_present):
                    raise ValueError("Missing polarization parameters in atom {:s} in res {:s}".format(atom.name, atom.residue.name))
                particle.frz_chg = params['frz_chg']
                particle.dyn_exp = params['dyn_exp']
                particle.frz_exp = params['frz_exp']
                pol_residues.add(atom.residue)
                add_atom = True

            #   pauli repulsion and dispersion force
            if any(pauli_params_present):
                if not all(pauli_params_present):
                    raise ValueError("Missing Pauli parameters in atom {:s} in res {:s}".format(atom.name, atom.residue.name))
                particle.pauli_exp = params['pauli_exp']
                particle.pauli_radii = params['pauli_radii']
                pauli_residues.add(atom.residue)
                add_atom = True
            if add_atom:
                force.add_particle(atom.index, particle)
                atom_indices[atom] = len(atom_indices)

        #   create polarization force
        if len(pol_residues) != 0:
            pol_force = force.create_fluc_dens_force()
            for res in pol_residues:
                fragment_idx = []
                for atom in res.atoms():
                    if atom in atom_indices:
                        fragment_idx.append(atom_indices[atom])
                pol_force.add_fragment(fragment_idx)

            bonds = []
            index_mapping = force.get_index_mapping()
            for bond in data.bonds:
                index_1 = bond.atom1
                index_2 = bond.atom2
                if index_1 in index_mapping and index_2 in index_mapping:
                    bonds.append((index_mapping[index_1], index_mapping[index_2]))
            pol_force.create_frz_exclusions_from_bonds(bonds, 3)

            

        #   create pauli repulsion and dispersion force
        if len(pauli_residues) != 0:
            pauli_force = force.create_disp_pauli_force()
            for res in pauli_residues:
                fragment_idx = []
                for atom in res.atoms():
                    if atom in atom_indices:
                        fragment_idx.append(atom_indices[atom])
                pauli_force.create_exclusions_from_fragment(fragment_idx)
            pauli_force.set_vdw_radii(MapID(self.dispersion_radii))
            pauli_force.set_C6_map(MapID(self.dispersion_C6))
            params = self.dispersion_params
            pauli_force.set_dispersion_params(params['s6'], params['a1'], params['a2'])
            

        external_force = CustomExternalForce('-Fx*x - Fy*y - Fz*z + fd_energy')
        external_force.addPerParticleParameter('Fx')
        external_force.addPerParticleParameter('Fy')
        external_force.addPerParticleParameter('Fz')
        external_force.addGlobalParameter('fd_energy', 0.0)
        index_mapping = force.get_index_mapping()
        for index in index_mapping:
            external_force.addParticle(index, (0, 0, 0))
        sys.addForce(external_force)

        self.force = force
        sys._amberFDData = {'ext_force': external_force, 'force': force, 'data': data}
        

    def postprocessSystem(self, sys, data, args):
        print("IN POSTPROCESS")
        #   exclude base-base interactions already accounted for in AmberFD
        for force in sys.getForces():
            if isinstance(force, NonbondedForce):
                index_mapping = self.force.get_index_mapping()
                for idx_1 in index_mapping:
                    for idx_2 in index_mapping:
                        force.addException(idx_1, idx_2, 0.0, 1.0, 0.0, replace=True)
                break


ff.parsers['AmberFDForce'] = AmberFDGenerator.parseElement

class MoleculeImporter():
    def __init__(self, ff_files, structure_files, solvate=False, **solvate_kwargs) -> None:
        ''' Loads multiple molecule files at once and stores the positions, topology, and forcefield information
        
        Parameters
        ----------
        files:  tuple of strings
            The location of the files to load. 
            Valid extensions are .pdb, .gro, .top, and .xml. 
            This mist MUST include a .xml forcefield file
        '''

        self.topology = None
        self.positions = None
        self.forcefield = None

        if isinstance(structure_files, str):
            structure_files = tuple([structure_files])
        elif isinstance(structure_files, list):
            structure_files = tuple(structure_files)
        if isinstance(ff_files, str):
            ff_files = tuple([ff_files])
        elif isinstance(ff_files, list):
            ff_files = tuple(ff_files)

        #   import force fields
        forcefield = ff.ForceField(*ff_files)

        #   determine tyoe of structure files provided
        extensions = [splitext(x)[-1] for x in structure_files]
        box_dims = None
        if set(extensions) == {'.pdb'}:
            pdb = PDBFile(*structure_files)
            topology  = pdb.getTopology()
            positions = pdb.getPositions(True)
        elif set(extensions) == {'.gro', '.top'}:
            gro = GromacsGroFile(structure_files[extensions.index('.gro')])
            top = GromacsTopFile(structure_files[extensions.index('.top')])
            positions = gro.getPositions(True)
            topology = top.topology
            topology.setPeriodicBoxVectors(gro.getPeriodicBoxVectors())
            box_dims = gro.getPeriodicBoxVectors()
        else:
            raise ValueError("Structure files must be .pdb or (.gro, .top)")

        #   add virtual sites to model
        modeller = Modeller(topology, positions)
        modeller.addExtraParticles(forcefield)

        #   solvate system
        if solvate:
            solvate_options = {'model': 'tip4pew', 'positiveIon': 'Na+', 'negativeIon': 'Cl-', 'ionicStrength': 0*uu.molar, 'neutralize': True}
            if box_dims is not None:
                solvate_options['boxVectors'] = box_dims
            else:
                solvate_options['padding'] = 1*uu.nanometers
            if len(solvate_kwargs) != 0:
                #   if solvation method is specified, remove defaults and use this instead
                solv_methods = ['boxSize', 'boxVectors', 'numAdded', 'padding']
                solv_methods_given = [x for x in solv_methods if x in solv_methods]
                for method in solv_methods_given:
                    solvate_options.pop(method, None)

                solvate_options.update(solvate_kwargs)
            modeller.addSolvent(forcefield, **solvate_options)

        #modeller.addExtraParticles(forcefield)

        self.positions = modeller.getPositions()
        self.topology = modeller.getTopology()
        self.forcefield = forcefield
        

class AmberFDSimulation(Simulation):
    def __init__(self, topology, system, integrator, platform=None, platformProperties=None, state=None):
        super().__init__(topology, system, integrator, platform=platform, platformProperties=platformProperties, state=state)

        print(system.__dict__)
        self._FD_solver = system._amberFDData['force']
        self._external_force = system._amberFDData['ext_force']
        self.index_mapping = dict(self._FD_solver.get_index_mapping())
        self.topology = topology

        self._dispPauliForce = self._FD_solver.get_disp_pauli_force()
        self._flucDensForce = self._FD_solver.get_fluc_dens_force()


    def _update_force(self):
        state = self.context.getState(getPositions=True)
        pos = state.getPositions(True)/uu.bohr
        fd_pos = np.array([pos[x] for x in self.index_mapping]).flatten()
        dim = len(self.index_mapping)

        # total_E = self._dispPauliForce.calc_energy(fd_pos)
        # total_E = self._dispPauliForce.get_pauli_energy()
        # total_E = self._dispPauliForce.get_disp_energy()
        # fd_forces = np.array(self._dispPauliForce.get_forces()).reshape((dim, 3))*49614.77640958472

        total_E = self._flucDensForce.calc_energy(fd_pos, True, False)
        fd_forces = np.array(self._flucDensForce.get_forces()).reshape((dim, 3))*49614.77640958472


        # energies = self._FD_solver.calc_energy_forces(fd_pos)
        # total_E = energies.total()
        # fd_forces = np.array(self._FD_solver.get_forces()).reshape((dim, 3))*49614.77640958472 # atomic to mm units

        atoms = list(self.topology.atoms())
        for (omm_particle, index), force in zip(self.index_mapping.items(), fd_forces):
            self._external_force.setParticleParameters(index, omm_particle, force)
            #print(atoms[omm_particle], omm_particle, index, force)


        self.context.setParameter('fd_energy', total_E*2625.5009)
        self._external_force.updateParametersInContext(self.context)
        #self.energies = energies


class BFGS(object):
    def __init__(self, context, constraints=None, out_pdb=None, topology=None):
        raise NotImplementedError("Still working on BFGS Minimizer")
        self._out_file = None
        self._topology = topology
        self._step_num = 0
        self._constraints = constraints
        self._context = context
        if out_pdb is not None and topology is not None:
            self._out_file = open(out_pdb, 'w')


    def _callback(self, pos):
        if self._out_file is not None:
            PDBFile.writeModel(self._topology, pos.reshape(-1,3)*nanometer, file=self._out_file, modelIndex=self._step_num)
        self._step_num += 1


    def minimize(self):
        #constraints = dict(zip(np.arange(64), ['Z']*64))

        init_state = self._context.getState(getForces=True, getEnergy=True, getPositions=True)
        init_pos = init_state.getPositions(True).value_in_unit(nanometer)
        init_energy, init_forces = self._target_func(init_pos, self._context, self._constraints)
        force_norms = [np.linalg.norm(f) for f in init_forces]
        print(" Initial max. force: {:15.3f} kJ/mol".format(np.max(force_norms)))
        print(" Initial energy:     {:15.3f} kJ/mol/nm".format(init_energy))


        self._step_num = 0
        args = (self._context, self._constraints)
        self._callback(init_pos)
        res = optimize.minimize(self._target_func, init_pos, args=args, method='L-BFGS-B', jac=True, callback=self._callback,
        options=dict(maxiter=500, disp=False, gtol=5))
        final_pos = res.x.reshape(-1,3)

        final_energy, final_forces = self._target_func(final_pos, self._context, self._constraints)
        force_norms = [np.linalg.norm(f) for f in final_forces]
        print(" Final max. force:   {:15.3f} kJ/mol".format(np.max(force_norms)))
        print(" Final energy:       {:15.3f} kJ/mol/nm".format(final_energy))


    def _target_func(self, pos, context, constraints=None):
        context.setPositions(pos.reshape(-1,3))
        state = context.getState(getEnergy=True, getForces=True)
        forces = state.getForces(asNumpy=True)
        energy = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
        forces = forces.value_in_unit(kilojoule_per_mole/nanometer)

        if constraints is not None:
            for n, constr in constraints.items():
                for idx in constr_2_idx[constr.upper()]:
                    forces[n][idx] *= 0

        return energy, -forces.flatten()
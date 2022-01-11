from copy import deepcopy, copy
import enum
from openmm.app import forcefield as ff, Simulation, Element, GromacsTopFile, GromacsGroFile, PDBFile, Modeller
from openmm.openmm import CustomExternalForce, Force, NonbondedForce, _openmm, Context as CT
import openmm.unit as uu
import sys
from os.path import *
import numpy as np
import openmm as mm
from openmm.app.simulation import string_types
from scipy.optimize import minimize
from itertools import combinations

try:
    sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
    from AmberFD import ANG2BOHR, AmberFD, FlucDens, VectorI, VectorD, VectorPairII, ParticleInfo, MapID
except:
    sys.path.insert(1, join(dirname(realpath(__file__)), '../build_minerva/'))
    from AmberFD import ANG2BOHR, AmberFD, FlucDens, VectorI, VectorD, VectorPairII, ParticleInfo, MapID

HARTREE_TO_KJ_MOL = 2625.5009
ANG_TO_BOHR = uu.angstrom.conversion_factor_to(uu.bohr)
FORCE_ATOMIC_TO_MM = ANG_TO_BOHR *10 * HARTREE_TO_KJ_MOL

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
        self.ct_coeff = 0.0

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

        generator.ct_coeff   = float(element.findall('PolarizationNonPairwise')[0].attrib['CT_coeff'])
        generator.damp_coeff = float(element.findall('PolarizationNonPairwise')[0].attrib['damp_coeff'])
        generator.damp_exp   = float(element.findall('PolarizationNonPairwise')[0].attrib['damp_exp'])


    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
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
            
            if 'LP' in atom.name and False:
               polar_params_present = [False]
               pauli_params_present = [False]

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
        #if len(pol_residues) != 0:
        if True:
            pol_force = force.create_fluc_dens_force()
            pol_force.set_dampening(self.damp_coeff, self.damp_exp)
            pol_force.set_ct_coeff(self.ct_coeff)
            for res in pol_residues:
                fragment_idx = []
                for atom in res.atoms():
                    if atom in atom_indices:
                        fragment_idx.append(atom_indices[atom])
                pol_force.add_fragment(fragment_idx)
            
            #   create frozen-frozen exclusions
            bonds = set()
            index_mapping = force.get_index_mapping()
            for bond in data.bonds:
                index_1 = bond.atom1
                index_2 = bond.atom2
                if index_1 in index_mapping and index_2 in index_mapping:
                    bond_tuple = tuple(sorted((index_mapping[index_1], index_mapping[index_2])))
                    bonds.add(bond_tuple)
            
            #  exclude the virtual sites with the same bonds as it's host atom
            vs_per_atom = {}
            for atom in data.atoms:
                if atom in data.virtualSites and atom.index in index_mapping:
                    site_data = data.virtualSites[atom]
                    excluded = site_data[2]
                    if excluded not in vs_per_atom:
                        vs_per_atom[excluded] = []
                    vs_per_atom[excluded].append(atom.index)
                    for bonded_to in data.bondedToAtom[excluded]:
                        atom1 = index_mapping[atom.index]
                        atom2 = index_mapping[bonded_to]
                        bond_tuple = tuple(sorted((atom1, atom2)))
                        #bond_tuple = tuple(sorted((atom.index, bonded_to)))
                        bonds.add(bond_tuple)

            # for atoms in vs_per_atom.values():
            #     if len(atoms) > 1:
            #         for combo in combinations(atoms, 2):
            #             bonds.add(tuple(sorted(combo)))
            
            pol_force.create_frz_exclusions_from_bonds(tuple(bonds), 3)
        
        #   create pauli repulsion and dispersion force
        #if len(pauli_residues) != 0:
        if True:
            pauli_force = force.create_disp_pauli_force()
            for res in pauli_residues:
                fragment_idx = []
                for atom in res.atoms():
                    if atom in atom_indices:
                        fragment_idx.append(atom_indices[atom])
                #pauli_force.create_exclusions_from_fragment(fragment_idx)
            pauli_force.create_exclusions_from_bonds(tuple(bonds), 3)
            pauli_force.set_vdw_radii(MapID(self.dispersion_radii))
            pauli_force.set_C6_map(MapID(self.dispersion_C6))
            params = self.dispersion_params
            pauli_force.set_dispersion_params(params['s6'], params['a1'], params['a2'])
            
        #   An "external" force is used to update forces on FD atoms.
        #   fd_offset is updated to cancel the distance dependant terms.
        #   fd_energy is as the actual energy generated by the force.
        #   This way, the total energy reported by contest.getState is correct
        external_force = CustomExternalForce('-Fx*x - Fy*y - Fz*z + fd_offset + fd_energy/fd_num_sites')
        external_force.addPerParticleParameter('Fx')
        external_force.addPerParticleParameter('Fy')
        external_force.addPerParticleParameter('Fz')
        external_force.addGlobalParameter('fd_offset', 0.0)
        external_force.addGlobalParameter('fd_energy', 0.0)
        index_mapping = force.get_index_mapping()
        for index in index_mapping:
            external_force.addParticle(index, (0, 0, 0))
        external_force.addGlobalParameter('fd_num_sites', max(1.0, external_force.getNumParticles()))
        sys.addForce(external_force)

        self.force = force
        sys._amberFDData = {'ext_force': external_force, 'force': force, 'data': data}

    def postprocessSystem(self, sys, data, args):
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
        

class Context(mm.Context):
    '''Construct a new Context in which to run a simulation. ONLY THE CPU PLATFORM IS CURRENTLY SUPPORTED

        Parameters
        ----------
        system : System
            the System which will be simulated
        integrator : Integrator
            the Integrator which will be used to simulate the System
    '''
    def __init__(self, *args):
        super().__init__(*args)

        self._system = args[0]
        self.integrator = args[1]
        self._FD_solver = self._system._amberFDData['force']
        self._external_force = self._system._amberFDData['ext_force']
        self._omm_index_to_FD = dict(self._FD_solver.get_index_mapping())
        self._n_sites = len(self._omm_index_to_FD)

        self._dispPauliForce = self._FD_solver.get_disp_pauli_force()
        self._flucDensForce = self._FD_solver.get_fluc_dens_force()

        self._current_energies = None
        self._current_forces = None

        self._dont_update = False
        self._currentStep = 0

    def _update_force(self):
        #return
        if self._dont_update: return
        state = self.getState(getPositions=True)
        all_pos = state.getPositions(True)
        pos_bohr = np.array([all_pos[x]/uu.bohr for x in self._omm_index_to_FD])
        pos_nm = np.array([all_pos[x]/uu.nanometers for x in self._omm_index_to_FD])
        
        self._current_energies = self._FD_solver.calc_energy_forces(pos_bohr.flatten())

        # total_E = self._dispPauliForce.calc_energy(pos_bohr.flatten())
        # total_E = self._dispPauliForce.get_pauli_energy()
        # total_E = self._dispPauliForce.get_disp_energy()
        # fd_forces = np.array(self._dispPauliForce.get_forces()).reshape((self._n_sites, 3))*49614.77640958472

        # total_E = self._flucDensForce.calc_energy(pos_bohr.flatten(), True, True)
        # fd_forces = np.array(self._flucDensForce.get_forces()).reshape((self._n_sites, 3))*49614.77640958472

        if False:
            ENG = HARTREE_TO_KJ_MOL
            self._dispPauliForce.calc_energy(pos_bohr.flatten())
            E_pauli = self._dispPauliForce.get_pauli_energy()*HARTREE_TO_KJ_MOL
            E_frz = self._flucDensForce.calc_energy(pos_bohr.flatten(), True, False)*HARTREE_TO_KJ_MOL
            E_pol = self._flucDensForce.calc_energy(pos_bohr.flatten(), False, True)*HARTREE_TO_KJ_MOL
            dist = np.linalg.norm(pos_nm[5] - pos_nm[22])*10
            #print(f'{dist=:10.5f}  {E_pauli=:14.5f}  {E_frz=:14.5f}  {E_pol=:14.5f}')

            eng_O_N =  self._flucDensForce.calc_one_frozen(pos_bohr.flatten(), 22, 5)
            eng_LP_N = self._flucDensForce.calc_one_frozen(pos_bohr.flatten(), 31, 5)
            eng_O_H =  self._flucDensForce.calc_one_frozen(pos_bohr.flatten(), 22, 9)
            eng_LP_H = self._flucDensForce.calc_one_frozen(pos_bohr.flatten(), 31, 9)

            eng_O_N_dp =  self._dispPauliForce.calc_one_pair(pos_bohr.flatten(), 22, 5)
            eng_O_H_dp =  self._dispPauliForce.calc_one_pair(pos_bohr.flatten(), 22, 9)

            # self._flucDensForce.nuclei[22] -=2
            # self._flucDensForce.nuclei[5] -=2
            # self._flucDensForce.frozen_pop[22] -=2
            # self._flucDensForce.frozen_pop[5] -=2
            # eng_O_N_2 =  self._flucDensForce.calc_one_frozen(pos_bohr.flatten(), 22, 5)
            # self._flucDensForce.nuclei[22]     +=2
            # self._flucDensForce.nuclei[5]      +=2
            # self._flucDensForce.frozen_pop[22] +=2
            # self._flucDensForce.frozen_pop[5]  +=2

            #print("PARTS: {:5d}  {:10.5f}  {:10.5f}  {:10.5f}  {:10.5f}  {:10.5f}  {:10.5f}  {:10.5f}".format(self._currentStep, dist, eng_O_N.frz*ENG, eng_LP_N.frz*ENG, eng_O_H.frz*ENG, eng_LP_H.frz*ENG, eng_O_N_dp.pauli*ENG, eng_O_N_2.frz*ENG))

        



            # self._dispPauliForce.calc_energy(pos_bohr.flatten())
            # print("PAULI ENERGY: ", self._dispPauliForce.get_pauli_energy()*HARTREE_TO_KJ_MOL)
            # print("FROZEN ENERGY: ", self._flucDensForce.calc_energy(pos_bohr.flatten(), True, False)*HARTREE_TO_KJ_MOL)
            # print("POL ENERGY: ", self._flucDensForce.calc_energy(pos_bohr.flatten(), False, True)*HARTREE_TO_KJ_MOL)
            # print("ELEC_ELEC", self._current_energies.elec_elec*HARTREE_TO_KJ_MOL)
            # print("ELEC_NUC", self._current_energies.elec_nuc*HARTREE_TO_KJ_MOL)
            # print("NUC_NUC", self._current_energies.nuc_nuc*HARTREE_TO_KJ_MOL)


        total_E = self._current_energies.total()
        fd_forces = np.array(self._FD_solver.get_forces()).reshape((self._n_sites, 3))*FORCE_ATOMIC_TO_MM

        if len(fd_forces) != 0:
            fd_offset = np.mean(fd_forces * pos_nm)*3
            fd_energy = total_E*HARTREE_TO_KJ_MOL
            self._current_forces = fd_forces
            # print(f'{fd_forces=}')
            # for stuff in zip(self._omm_index_to_FD.values(), pos_nm*10):
            #     print(stuff)
            
            #   update external force and global parameters
            #print(fd_forces[np.argmax(np.linalg.norm(fd_forces, axis=1))])
            for (omm_particle, index), force in zip(self._omm_index_to_FD.items(), fd_forces):
                self._external_force.setParticleParameters(index, omm_particle, force)
                #print(index, omm_particle, force)
            self.setParameter('fd_offset', fd_offset)
            self.setParameter('fd_energy', fd_energy)
            self._external_force.updateParametersInContext(self)
        
    def getState(self, getPositions=False, getVelocities=False,
                 getForces=False, getEnergy=False, getParameters=False,
                 getParameterDerivatives=False, getIntegratorParameters=False,
                 enforcePeriodicBox=False, groups=-1):

        if not self._dont_update:
            if getEnergy or getForces or getParameterDerivatives:
                if groups == -1:
                    self._update_force()
                elif isinstance(groups, int):
                    mask = (1 << self._external_force.getForceGroup())
                    if mask == groups:
                        self._update_force()
                elif self._external_force.getForceGroup() in groups:
                    self._update_force()

        return super().getState(getPositions=getPositions, getVelocities=getVelocities,
                 getForces=getForces, getEnergy=getEnergy, getParameters=getParameters,
                 getParameterDerivatives=getParameterDerivatives, getIntegratorParameters=getIntegratorParameters,
                 enforcePeriodicBox=enforcePeriodicBox, groups=groups)



class AmberFDSimulation(Simulation):
    def __init__(self, topology, system, integrator, state=None):
        """Create a Simulation.

        Parameters
        ----------
        topology : Topology
            A Topology describing the the system to simulate
        system : System or XML file name
            The OpenMM System object to simulate (or the name of an XML file
            with a serialized System)
        integrator : Integrator or XML file name
            The OpenMM Integrator to use for simulating the System (or the name
            of an XML file with a serialized System)
        state : XML file name=None
            The name of an XML file containing a serialized State. If not None,
            the information stored in state will be transferred to the generated
            Simulation object.
        """
        ######### Copied from openmm.py 7.6.0 #########
        self.topology = topology
        ## The System being simulated
        if isinstance(system, string_types):
            with open(system, 'r') as f:
                self.system = mm.XmlSerializer.deserialize(f.read())
        else:
            self.system = system
        ## The Integrator used to advance the simulation
        if isinstance(integrator, string_types):
            with open(integrator, 'r') as f:
                self.integrator = mm.XmlSerializer.deserialize(f.read())
        else:
            self.integrator = integrator
        ## The index of the current time step
        # self.currentStep = 0
        ## A list of reporters to invoke during the simulation
        self.reporters = []

        #######   only the CPU platform is currently supported, so use our overrided Context object  ######
        self.context = Context(self.system, self.integrator, mm.Platform_getPlatformByName('CPU'))
 
        if state is not None:
            with open(state, 'r') as f:
                self.context.setState(mm.XmlSerializer.deserialize(f.read()))
        ## Determines whether or not we are using PBC. Try from the System first,
        ## fall back to Topology if that doesn't work
        try:
            self._usesPBC = self.system.usesPeriodicBoundaryConditions()
        except Exception: # OpenMM just raises Exception if it's not implemented everywhere
            self._usesPBC = topology.getUnitCellDimensions() is not None
        ######### End Copied  #########

        self._integrator_step_old = self.integrator.step
        self.integrator.step = self._integrator_step_override

    def _integrator_step_override(self, n_steps):
        for n in range(n_steps):
            self.context._currentStep = self.currentStep
            self.context._update_force()
            self._integrator_step_old(1)


    def minimizeEnergy(self, tolerance=10*uu.kilojoules_per_mole/uu.nanometer, maxIterations=500, PDBOutFile=None):
        if isinstance(PDBOutFile, str):
            solver = BFGS(self.context, out_pdb=PDBOutFile, topology=self.topology)
        else:
            solver = BFGS(self.context, topology=self.topology)

        solver.minimize(tolerance, maxIterations)


class BFGS(object):
    def __init__(self, context, out_pdb=None, topology=None):
        #raise NotImplementedError("Still working on BFGS Minimizer")
        self._out_file = None
        self._topology = topology
        self._step_num = 0
        self._context = context
        self._system = context.getSystem()
        if out_pdb is not None and topology is not None:
            self._out_file = open(out_pdb, 'w')

        #   only use atoms that are not lone pairs (element == None) or have non-zero mass
        self._force_index_list = []
        for atom in topology.atoms():
            if atom.element is None: continue
            mass = atom.element.mass
            mass = mass/mass.unit
            if mass == 0: continue
            self._force_index_list.append(atom.index)

        self._current_all_pos = np.empty(len(self._force_index_list))
        self._current_energy = None
        self._current_forces = None

        #   params for harmonic cost functions used to add constraints
        self._constraint_tol = self._context.integrator.getConstraintTolerance()
        self._working_constraint_tol = np.max((1e-4, self._constraint_tol))
        self._k = 100/self._working_constraint_tol

        self._silent = False
        self._print_header = True

    def _callback(self, pos):
        
        if not self._silent:
            if self._print_header:
                print("\n {:>5s}  {:>15s}  {:>15s}".format("Step", "Func (kJ/mol)", "Max Force"))
                self._print_header = False

            max_force_norm = np.linalg.norm(self._current_forces, axis=1).max()
            print(" {:5d}  {:15.3f}  {:15.3f}".format(self._step_num, self._current_energy, max_force_norm))
        
        if self._out_file is not None:
            PDBFile.writeModel(self._topology, pos.reshape(-1,3)*uu.nanometer, file=self._out_file, modelIndex=self._step_num)
        self._step_num += 1

    def _test(self, x0):
        energy, grad = self._target_func(x0)

        eps = 1e-8
        for n in range(len(x0)):
            xp = x0.copy()
            xp[n] += eps
            xn = x0.copy()
            xn[n] -= eps

            energy_p = self._target_func(xp)[0]
            energy_n = self._target_func(xn)[0]
            num_grad = (energy_p - energy_n)/(2*eps)
            print(n, num_grad, grad[n])

        max_step = 0.0000001
        for step in np.arange(0, max_step, max_step/20):
            new_pos = x0 - grad*step
            new_energy, new_grad = self._target_func(new_pos)
            print("{:10.3e}  {:15.8f}".format(step, new_energy))

    def minimize(self, tolerance=10*uu.kilojoules_per_mole/uu.nanometer, maxIterations=500):
        self._context.applyConstraints(self._working_constraint_tol)
        gtol = tolerance
        init_state = self._context.getState(getForces=True, getEnergy=True, getPositions=True)
        init_pos_all = init_state.getPositions(True).value_in_unit(uu.nanometer)
        self._current_all_pos = init_pos_all
        init_pos = init_pos_all.flatten()
        init_energy, init_forces = self._target_func(init_pos)
        force_norms = [np.linalg.norm(f) for f in init_forces]
        print(" Initial max. force: {:15.3f} kJ/mol".format(np.max(force_norms)))
        print(" Initial energy:     {:15.3f} kJ/mol/nm".format(init_energy))
        dim = len(init_pos)

        # Repeatedly minimize, steadily increasing the strength of the springs until all constraints are satisfied.

        self._step_num = 0
        self._callback(init_pos)
        prevMaxError = 1e10

        while(True):
            self._print_header = True
            res = minimize(self._target_func, init_pos, 
                method='l-bfgs-b', 
                jac=True, 
                callback=self._callback,
                options={'maxiter': maxIterations, 'disp': False, 'gtol': gtol, 'maxfun': maxIterations*3*dim})

            print(" L-BFGS-B Solver says: ", res.message)
            #   Check whether all constraints are satisfied.

            positions = self._context.getState(getPositions=True).getPositions()
            maxError = 0.0
            for n in range(self._system.getNumConstraints()):
                p1, p2, distance = self._system.getConstraintParameters(n)
                distance = distance/uu.nanometer
                delta = (positions[p2] - positions[p1])/uu.nanometer
                r = np.sqrt(np.dot(delta, delta))
                error = np.abs(r-distance)
                if error > maxError:
                    maxError = error
            
            if self._step_num > maxIterations:
                print(" Maximum number of iterations exceded")
                break

            if maxError <= self._working_constraint_tol:
                print(" All constraints are satisified. Current Max error = {:15.10f} nm".format(maxError))
                break # All constraints are satisfied.
            #context.setPositions(initialPos)
            if maxError >= prevMaxError:
                print(" Constraint cost increase didn't help. Current Max error = {:15.10f} nm".format(maxError))
                break # Further tightening the springs doesn't seem to be helping, so just give up.
            prevMaxError = maxError

            self._k *= 10
            if maxError > 100*self._working_constraint_tol:
                print(" Too far from valid soln. Resetting positions. Current Max error = {:15.10f} nm".format(maxError))
                # We've gotten far enough from a valid state that we might have trouble getting
                # back, so reset to the original positions.
                #init_pos = init_pos_all[self._force_index_list].flatten()
                init_pos = np.copy(init_pos_all).flatten()
            else:
                print(" Increasing constraint cost. Current Max error = {:15.10f} nm".format(maxError))
                init_pos = res.x

        final_pos = res.x.reshape(-1,3)
        final_energy, final_forces = self._target_func(final_pos)
        force_norms = [np.linalg.norm(f) for f in final_forces]
        print(" Final max. force:   {:15.3f} kJ/mol".format(np.max(force_norms)))
        print(" Final energy:       {:15.3f} kJ/mol/nm".format(final_energy))
        if self._out_file is not None:
            self._out_file.close()

        # If necessary, do a final constraint projection to make sure they are satisfied
        # to the full precision requested by the user.
        if (self._constraint_tol < self._working_constraint_tol):
            self._context.applyConstraints(self._working_constraint_tol)


    def _target_func(self, pos):

        # pos_to_set = np.copy(self._current_all_pos).reshape(-1, 3)

        # curr_pos = context.getState(getPositions=True).getPositions(True)/uu.nanometer
        # updated_pos = np.copy(curr_pos)
        # pos = np.reshape(pos, (-1, 3))
        # for idx, new_pos in zip(self._force_index_list, pos):
        #     updated_pos[idx] = new_pos
        # context.setPositions(updated_pos)
        # context.applyConstraints(self._working_constraint_tol)

        new_pos = np.reshape(pos, (-1, 3))
        self._context.setPositions(new_pos)
        self._context.computeVirtualSites()
        #self._context.applyConstraints(self._working_constraint_tol)

        state = self._context.getState(getEnergy=True, getForces=True, getPositions=True)
        forces_all = state.getForces(asNumpy=True).value_in_unit(uu.kilojoule_per_mole/uu.nanometer)
        pos_all = state.getPositions(True).value_in_unit(uu.nanometer)
        self._current_all_pos = pos_all.flatten()
        energy = state.getPotentialEnergy().value_in_unit(uu.kilojoule_per_mole)
        
        #   zero out all forces not being updated
        for n in range(len(forces_all)):
            if n not in self._force_index_list:
                forces_all[n] *= 0

        #   Add harmonic forces for any constraints 
        for n in range(self._system.getNumConstraints()):
            p1, p2, distance = self._system.getConstraintParameters(n)
            distance = distance/uu.nanometer
            delta = (pos_all[p2] - pos_all[p1])
            r2 = np.dot(delta, delta)
            r = np.sqrt(r2)
            delta *= 1/r
            dr = r-distance
            kdr = self._k*dr
            energy += 0.5*kdr*dr

            forces_all[p1] += kdr*delta
            forces_all[p2] -= kdr*delta


    
        forces = np.array([forces_all[idx] for idx in self._force_index_list])
        atoms = list(self._topology.atoms())
        force_norms = np.linalg.norm(forces, axis=1)
        pos_diff = np.linalg.norm(new_pos - pos_all)
        max_diff = np.max(pos_diff)
        #print("MAX DIFF: ", np.argmax(max_diff), max_diff, atoms[np.argmax(max_diff)])
        # for n, idx in enumerate(self._force_index_list):
        #     print(idx, atoms[idx], force_norms[n])
        max_force_idx = np.argmax(force_norms)
        max_force_atom = atoms[self._force_index_list[max_force_idx]]
        #print(max_force_idx, np.max(force_norms), max_force_atom, max_force_atom.residue.id, max_force_atom.id)
        # exit()

        #print(len(forces), len(forces_all), forces)

        self._current_energy = energy
        self._current_forces = forces
        return energy, -forces_all.flatten()
 
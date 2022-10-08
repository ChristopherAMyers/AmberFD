from copy import deepcopy, copy
import enum
import os
from statistics import harmonic_mean
from openmm.app import forcefield as ff, Simulation, Element, GromacsTopFile, GromacsGroFile, PDBFile, Modeller
from openmm.openmm import CustomExternalForce, Vec3, NonbondedForce, CustomNonbondedForce, Context as CT
import openmm.unit as uu
from openmm.app import Topology
import sys
from os.path import *
import numpy as np
import openmm as mm
from openmm.app.simulation import string_types
from scipy.optimize import minimize, OptimizeResult
import time
import warnings

#from . import LJOverrides

try:
    sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
    from _AmberFD import ANG2BOHR, AmberFD, FlucDens, ParticleInfo, MapID
except:
    print(" WARNING: Using MINERVA Build")
    sys.path.insert(1, join(dirname(realpath(__file__)), '../build_minerva/'))
    from _AmberFD import ANG2BOHR, AmberFD, FlucDens, ParticleInfo, MapID


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
        self.damp_type = None
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

        pol_opts = element.findall('PolarizationNonPairwise')[0].attrib
        generator.ct_coeff   = float(pol_opts['CT_coeff'])
        generator.damp_coeff = float(pol_opts['damp_coeff'])
        generator.damp_exp   = float(pol_opts['damp_exp'])
        if 'damp_type' in pol_opts:
            generator.damp_type  = pol_opts['damp_type']


    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        force = AmberFD()
        pol_residues = set()
        pauli_residues = set()
        atom_indices = {}
        optional_hardness = []
        found_hardness = False
        for atom, params in data.atomParameters.items():
            add_atom = False
            nuclei = atom.element.atomic_number if isinstance(atom.element, Element) else 0
            particle = ParticleInfo(nuclei)
            polar_params_present =   [x in params for x in ['frz_chg', 'frz_exp', 'dyn_exp']]
            pauli_params_present =   [x in params for x in ['pauli_exp', 'pauli_radii']]
            harness_params_present = [x in params for x in ['hard_inc']]
            

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
                #   optional hardness increase
                if any(harness_params_present):
                    optional_hardness.append(params['hard_inc'])
                    found_hardness = True
                else:
                    optional_hardness.append(0.0)
                force.add_particle(atom.index, particle)
                atom_indices[atom] = len(atom_indices)
        
        #   if no atoms were found to have parameters, then don't use the custom force
        if not atom_indices:
            print("Ignoring AmberFD")
            sys._amberFDData = {'ext_force': None, 'force': None, 'data': None}
            return

        if sys.usesPeriodicBoundaryConditions():
            boxVecs = sys.getDefaultPeriodicBoxVectors()
            boxVecs = [x.value_in_unit(uu.bohr) for x in boxVecs]
            force.set_use_PBC(True, boxVecs[0][0], boxVecs[1][1], boxVecs[2][2])

        #   create polarization force
        #if len(pol_residues) != 0:
        if True:
            pol_force = force.create_fluc_dens_force()
            if found_hardness:
                pol_force.set_additional_hardness(optional_hardness)
            
            damp_type = pol_force.Quadratic
            if self.damp_type is not None:
                if self.damp_type.lower() == 'linear':
                    damp_type = pol_force.Linear
                elif self.damp_type.lower() == 'Quadratic':
                    damp_type = pol_force.Quadratic
                else:
                    raise ValueError("Polarization Dampening must be 'linear' or 'quadratic'")
            pol_force.set_dampening(self.damp_coeff, self.damp_exp, damp_type)
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
        if sys._amberFDData['force'] is None: return
        for force in sys.getForces():
            if isinstance(force, NonbondedForce):
                index_mapping = self.force.get_index_mapping()
                for idx_1 in index_mapping:
                    for idx_2 in index_mapping:
                        force.addException(idx_1, idx_2, 0.0, 1.0, 0.0, replace=True)
            elif isinstance(force, CustomNonbondedForce):
                index_mapping = self.force.get_index_mapping()
                current_exclusions = set([tuple(force.getExclusionParticles(n)) for n in range(force.getNumExclusions())])
                for idx_1 in index_mapping:
                    for idx_2 in index_mapping:
                        if (idx_1, idx_2) not in current_exclusions and (idx_2, idx_1) not in current_exclusions:
                            force.addExclusion(idx_1, idx_2)
                            current_exclusions.add((idx_1, idx_2))
                            current_exclusions.add((idx_2, idx_1))




ff.parsers['AmberFDForce'] = AmberFDGenerator.parseElement

class MoleculeImporter():
    def __init__(self, ff_files, structure_files, solvate=False, **solvate_kwargs) -> None:
        ''' 
        Loads multiple molecule files at once and stores the positions, topology, and forcefield information into one object. MoleculeImporter can also solvate the system, using all valid solvation options that are used in openmm.app.modeller.addSolvent()
        
        Parameters
        ----------
        ff_files: tuple or string
            The .xml OpenMM forcefield files to use 
        structure_files:  tuple or string
            Molecule files to load, depending on the type. 
            If using a .pdb file, only one file needs to be supplied.
            If using GROMACS files, then a tuple of (.gro, .top) files should be used
        solvate: bool
            Whether or not to solvate the system. This is controlled using the following solvent_kwargs, which take on the same options used to solvate the system that are used with  openmm.app.modeller.addSolvent(). However, the kward 'forcefield' will be ignored: instead, force field files used with the 'ff_files' option will be used.
            
        Keyword Args
        ------------
        model : str='tip3p'
            the water model to use.  Supported values are 'tip3p', 'spce', 'tip4pew', and 'tip5p'.
        boxSize : Vec3=None
            the size of the box to fill with water
        boxVectors : tuple of Vec3=None
            the vectors defining the periodic box to fill with water
        padding : distance=None
            the padding distance to use
        numAdded : int=None
            the total number of molecules (waters and ions) to add
        positiveIon : string='Na+'
            the type of positive ion to add.  Allowed values are 'Cs+', 'K+', 'Li+', 'Na+', and 'Rb+'
        negativeIon : string='Cl-'
            the type of negative ion to add.  Allowed values are 'Cl-', 'Br-', 'F-', and 'I-'. Be aware
            that not all force fields support all ion types.
        ionicStrength : concentration=0*molar
            the total concentration of ions (both positive and negative) to add.  This
            does not include ions that are added to neutralize the system.
            Note that only monovalent ions are currently supported.
        neutralize : bool=True
            whether to add ions to neutralize the system
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

        #   add additional bond deffinitions
        #   first create a blank topology to force openmm.app.Topology
        #   to load their standard bonds
        _blank_top = Topology()
        _blank_top.createStandardBonds()
        #   then add our own deffinitions
        res_file = join(dirname(__file__), 'residues.xml')
        if isfile(res_file):
            Topology.loadBondDefinitions(res_file)
        for file in ff_files:
            try:
                # this handles either filenames or open file-like objects
                Topology.loadBondDefinitions(file)
            except FileNotFoundError:
                for dataDir in ff._getDataDirectories(): 
                    f = join(dataDir, file)
                    if isfile(f):
                        Topology.loadBondDefinitions(f)
                        break

        #   import force fields
        forcefield = ff.ForceField(*ff_files)

        #   determine type of structure files provided
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

        #   set box size if it wasn't already
        if self.topology.getPeriodicBoxVectors() is None:
            if 'boxSize' in solvate_kwargs:
                boxSize = solvate_kwargs['boxSize']
                if uu.is_quantity(boxSize):
                    boxSize = boxSize.value_in_unit(uu.nanometer)
                vectors = (Vec3(boxSize[0], 0, 0), Vec3(0, boxSize[1], 0), Vec3(0, 0, boxSize[2]))
                self.topology.setPeriodicBoxVectors(vectors)
            elif 'boxVectors' in solvate_kwargs:
                boxVectors = solvate_kwargs['boxVectors']
                if uu.is_quantity(boxVectors[0]):
                    boxVectors = (boxVectors[0].value_in_unit(uu.nanometer), boxVectors[1].value_in_unit(uu.nanometer), boxVectors[2].value_in_unit(uu.nanometer))
                vectors = boxVectors
                self.topology.setPeriodicBoxVectors(vectors)
            elif topology.getPeriodicBoxVectors() is not None:
                self.topology.setPeriodicBoxVectors(topology.getPeriodicBoxVectors())
            else:
                min_pos = np.min(positions, axis=0)
                max_pos = np.max(positions, axis=0)
                boxSize = max_pos - min_pos
                if uu.is_quantity(boxSize):
                    boxSize = boxSize.value_in_unit(uu.nanometer)
                vectors = (Vec3(boxSize[0], 0, 0), Vec3(0, boxSize[1], 0), Vec3(0, 0, boxSize[2]))
                self.topology.setPeriodicBoxVectors(vectors)       

class Context(mm.Context):
    '''Construct a new Context in which to run a simulation. At the moment, only the CPU platform is implimented. 

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
        self._FD_forces = None
        self._external_force = self._system._amberFDData['ext_force']
        self._omm_index_to_FD = dict(self._FD_solver.get_index_mapping())
        self._omm_indicies = list(self._omm_index_to_FD.keys())
        self._n_sites = len(self._omm_index_to_FD)

        self._dispPauliForce = self._FD_solver.get_disp_pauli_force()
        self._flucDensForce = self._FD_solver.get_fluc_dens_force()

        self._current_energies = None
        self._current_forces = None

        self._dont_update = False
        self._currentStep = 0

        self._self_pos = None
        self._wall_time = 0

    def _get_fd_pos(self):
        state = self.getState(getPositions=True)
        all_pos = state.getPositions(True)
        pos_bohr = (all_pos[self._omm_indicies])/uu.bohr

        return pos_bohr.flatten()

    def get_interaction_energy(self):
        if self._dont_update: return

        state = self.getState(getPositions=True)
        all_pos = state.getPositions(True)
        pos_bohr = (all_pos[self._omm_indicies])/uu.bohr    
        return self._FD_solver.calc_energy_forces(pos_bohr.flatten(), True)


    def _update_force(self):     
        if self._dont_update: return

        state = self.getState(getPositions=True)
        all_pos = state.getPositions(True)
        pos_bohr = (all_pos[self._omm_indicies])/uu.bohr
        pos_nm = pos_bohr * uu.bohr.conversion_factor_to(uu.nanometer)        
        self._current_energies = self._FD_solver.calc_energy_forces(pos_bohr.flatten())

        total_E = self._current_energies.total()
        fd_forces = np.array(self._FD_solver.get_forces()).reshape((self._n_sites, 3))*FORCE_ATOMIC_TO_MM

        if False:

            forces_fluc_all = np.array(self._flucDensForce.get_forces() ).reshape((self._n_sites, 3))*FORCE_ATOMIC_TO_MM
            forces_dp_all =   np.array(self._dispPauliForce.get_forces()).reshape((self._n_sites, 3))*FORCE_ATOMIC_TO_MM

            self._flucDensForce.initialize_calculation()
            self._dispPauliForce.initialize()
            ENG = HARTREE_TO_KJ_MOL
            pos_flat = pos_bohr.flatten()

            ID_1 = self._omm_index_to_FD[650 - 1]
            ID_2 = self._omm_index_to_FD[791 - 1]

            print(np.linalg.norm(pos_nm[ID_1] - pos_nm[ID_2])*10)
            print(self._FD_solver.get_use_PBC())
            print(self._FD_solver.getDeltaR(pos_flat, ID_1, ID_2).r)
            print()
            eng = self._FD_solver.calc_one_pair(pos_flat, ID_1, ID_2)
            forces_dp =   np.array(self._dispPauliForce.get_forces()).reshape((self._n_sites, 3))*FORCE_ATOMIC_TO_MM
            forces_fluc = np.array(self._flucDensForce.get_forces() ).reshape((self._n_sites, 3))*FORCE_ATOMIC_TO_MM
            # print(eng.pauli*ENG)
            # print(eng.pauli_wall*ENG)
            # print(eng.frz*ENG)
            # print(eng.total()*ENG)
            print(forces_dp[ID_2], forces_dp[ID_1])
            print(forces_fluc[ID_2], forces_fluc[ID_1])
            print()
            print(forces_dp_all[ID_2], forces_dp_all[ID_1])
            print(forces_fluc_all[ID_2], forces_fluc_all[ID_1])


        if len(fd_forces) != 0:
            fd_offset = np.mean(fd_forces * pos_nm)*3
            fd_energy = total_E*HARTREE_TO_KJ_MOL
            self._current_forces = fd_forces
            
            #   update external force and global parameters
            for (omm_particle, index), force in zip(self._omm_index_to_FD.items(), fd_forces):
                #print(index, omm_particle, force)
                self._external_force.setParticleParameters(index, omm_particle, force)
                
            self.setParameter('fd_offset', fd_offset)
            self.setParameter('fd_energy', fd_energy)
            self._external_force.updateParametersInContext(self)
            self._FD_forces = fd_forces
        
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

    def create_cube_file(self, cube_file, grid_range=None, grid_spacing=0.3, grid_padding=3.0, density_type='delta'):

        #   update density and get current positions
        state = self.getState(getPositions=True, getEnergy=True)
        all_pos = state.getPositions(True)
        pos_bohr = (all_pos[self._omm_indicies])/uu.bohr

        #   determine bounds
        if grid_range is None:
            x_min, y_min, z_min = np.min(pos_bohr, axis=0) - 3*ANG_TO_BOHR
            x_max, y_max, z_max = np.max(pos_bohr, axis=0) + 3*ANG_TO_BOHR
        else:
            ((x_min, x_max), (y_min, y_max), (z_min, z_max)) = np.array(grid_range)*ANG_TO_BOHR

        #   create grid points
        x = np.arange(x_min, x_max, grid_spacing*ANG_TO_BOHR)
        y = np.arange(y_min, y_max, grid_spacing*ANG_TO_BOHR)
        z = np.arange(z_min, z_max, grid_spacing*ANG_TO_BOHR)
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')

        nx, ny, nz = len(x), len(y), len(z)
        dim = nx*ny*nz
        a = np.reshape(xx, dim)
        b = np.reshape(yy, dim)
        c = np.reshape(zz, dim)
        coords = np.array([a, b, c]).T

        #   determine type of density to print
        if not isinstance(density_type, str):
            raise TypeError("create_cube_file(): density_dype must be 'delta', 'frozen', or 'all'")
        if density_type.lower() == 'delta':
            density_type = self._flucDensForce.Delta
        elif density_type.lower() == 'frozen':
            density_type = self._flucDensForce.Frozen
        elif density_type.lower() == 'all':
            density_type = self._flucDensForce.All

        #   compute the actual density
        density = self._flucDensForce.calc_density(coords.flatten(), pos_bohr.flatten(), density_type)
        density = np.array(density)

        if isinstance(cube_file, str):
            file = open(cube_file, 'w')
        else:
            file = cube_file

        all_nuclei = np.array(self._flucDensForce.get_params_by_name('nuclei'), dtype=int)
        keep_idx = np.where(all_nuclei > 0)
        out_pos = pos_bohr[keep_idx]
        out_nuclei = all_nuclei[keep_idx]
        out_nuclei[np.where(out_nuclei > 1)] += 2

        #   write cube file header
        file.write("Cube file generated with AmberFD\n")
        file.write("Cube file generated with AmberFD\n")
        file.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n"
        .format(len(out_pos), x_min, y_min, z_min))
        file.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n".format(nx, grid_spacing*ANG_TO_BOHR, 0.0, 0.0))
        file.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n".format(ny, 0.0, grid_spacing*ANG_TO_BOHR, 0.0))
        file.write("{:5d} {:12.6f} {:12.6f} {:12.6f}\n".format(nz, 0.0, 0.0, grid_spacing*ANG_TO_BOHR))
        for n, nuclei in enumerate(out_nuclei):
            file.write("{:5d} {:12.6f} {:12.6f} {:12.6f} {:12.6f}\n"
            .format(int(nuclei), nuclei, out_pos[n][0], out_pos[n][1], out_pos[n][2]))

        #   write cube file data
        remainder = density.shape[0] % 6
        first_part = density[0:-remainder].reshape((-1, 6))
        second_part = density[-remainder:]
        np.savetxt(file, first_part)
        np.savetxt(file, [second_part])

        file.close()


class AmberFDSimulation(Simulation):
    '''
    This is a subsitute to openmm.simulation objects and enables the use of AmberFD polarization.
    
    To use it, create a simulation just like you normally would with OpenMM, but instead call AmberFDSimulation. If a system is provided that was not created using the AmberFD forcefield, then this will initialize using openmm.simulation instead. All of the same members that belong to openmm.simulation are also available, such as minimizeEnergy() and step(). Some of these members have additional options that are otherwise unavailable with vanilla openmm.simulation. These options are available regardless if an AmberFD force field is used or not.
    '''
    def __init__(self, topology, system, integrator, platform=None, platformProperties=None, state=None, FDProperties=None):
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
        platform : Platform=None
            If not None, the OpenMM Platform to use. When using AmberFD, only the CPU platform is available.
        platformProperties : map=None
            If not None, a set of platform-specific properties to pass to the
            Context's constructor. NOT IMPLIMENTED IN AmberFD
        state : XML file name=None
            The name of an XML file containing a serialized State. If not None,
            the information stored in state will be transferred to the generated
            Simulation object. NOT YET IMPLIMENTED IN AmberFD
        """
        if hasattr(system, '_amberFDData'):
            n_threads = 1
            if system._amberFDData['force'] is not None:
                if platform is not None:
                    if platform.getName() not in ['Reference', 'CPU']:
                        warnings.warn('Only the Reference or CPU Platform is implimented in AmberFD. AmberFD Forces will be calculated on the CPU')
                        #raise NotImplementedError("Only the Reference or CPU Platform is implimented in AmberFD")
                if FDProperties is not None:
                    if "num_threads" in FDProperties:
                        n_threads = int(FDProperties['num_threads'])
                
                # if platform is not None and FDProperties is not None:
                #     if platform.getName() == 'CPU' and "num_threads" in FDProperties:
                #         system._amberFDData['force'].set_threads(int(FDProperties['num_threads']))
                if state is not None:
                    raise NotImplementedError("States are not yet implimented in AmberFD")
            else:
                print("INFO: Using Standard OpenMM System")
                super().__init__(topology, system, integrator, platform=platform, platformProperties=platformProperties, state=state)
                return

            #   set multithreading
            if 'OMP_NUM_THREADS' in os.environ:
                try:
                    n_threads = int(os.environ['OMP_NUM_THREADS'])
                except ValueError:
                    n_threads = 1
            if n_threads > 1:
                system._amberFDData['force'].set_threads(n_threads)

        else:
            super().__init__(topology, system, integrator, platform=platform, platformProperties=platformProperties, state=state)
            return

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

        ## The index of the current time step
        #if hasattr(self, 'currentStep'):
        self.currentStep = 0

    def _integrator_step_override(self, n_steps):
        ''' Overrides the step function for the integrator. This is needed to that
            AmberFD can calculate it's force and updates the CustomExternalForce
            parameters before the actual integration step is taken in OpenMM.'''
        for n in range(n_steps):
            self.context._currentStep = self.currentStep
            
            self._integrator_step_old(1)
            self.context._update_force()

    def minimizeEnergy(self, tolerance=10*uu.kilojoules_per_mole/uu.nanometer, maxIterations=500, PDBOutFile=None, output_interval=1, excludeResidueNames=None):
        '''
        Minimize the energy of the system w.r.t. the molecular coordinates.

        Parameters
        ----------
        tolerance : force
            This specifies how precisely the energy minimum must be located.  Minimization
            is halted once themaximum value of all force components reaches
            this tolerance.
        maxIterations : int
            The maximum number of iterations to perform.  If this is 0,
            minimization is continued until the results converge without regard
            to how many iterations it takes.
        PDBOutFile: string
            .pdb file to print minimized coordinates to as the minimization procedes.
        output_interval: int
            Number of steps in between minimization to print coordinates to PDBOutFile
        excludeResidueNames: iterable
            A list of residue names to exclude in PDBOutFile. For example, one may want to exclude all waters and ions from KCl, so ('HOH', 'CL', 'K') would be used
        '''
        if isinstance(PDBOutFile, str):
            minimizer = Minimizer(self.context, self.topology, out_pdb=PDBOutFile, excludeResidueNames=excludeResidueNames)
        else:
            minimizer = Minimizer(self.context, self.topology, excludeResidueNames=excludeResidueNames)

        minimizer.minimize(tolerance, maxIterations, output_interval=output_interval)


class Minimizer(object):
    def __init__(self, context, topology, out_pdb=None, excludeResidueNames=None):
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
        self._zero_force_index_list = []
        for atom in topology.atoms():
            if atom.element is None:
                self._zero_force_index_list.append(atom.index)
                continue
            #mass = atom.element.mass
            mass = self._system.getParticleMass(atom.index)
            mass = mass/mass.unit
            if mass == 0:
                self._zero_force_index_list.append(atom.index)
                continue
            self._force_index_list.append(atom.index)

        self._current_all_pos = np.empty(len(self._force_index_list))
        self._current_energy = None
        self._current_forces = None

        #   params for harmonic cost functions used to add constraints
        integrator = self._context.getIntegrator()
        self._constraint_tol = integrator.getConstraintTolerance()
        self._working_constraint_tol = np.max((1e-4, self._constraint_tol))
        self._k = 100/self._working_constraint_tol

        self._silent = False
        self._print_header = True
        self._total_time = 0
        self._force_time = 0

        #   get constraint information
        self._constraint_idx_1 = []
        self._constraint_idx_2 = []
        self._constraint_distances = []
        for n in range(self._system.getNumConstraints()):
            p1, p2, distance = self._system.getConstraintParameters(n)
            self._constraint_idx_1.append(p1)
            self._constraint_idx_2.append(p2)
            self._constraint_distances.append(distance/uu.nanometer)

        self._keep_idx = self._force_index_list.copy()
        self._other_to_report = None
        self._print_interval = 1

        #   residues to exclude from output
        self._exclude_residues = None
        if excludeResidueNames is not None:
            if isinstance(excludeResidueNames, str):
                self._exclude_residues = tuple([excludeResidueNames])
            else:
                self._exclude_residues = tuple(excludeResidueNames)
        self._keep_atom_idx = []
        self._excluded_top = None

    def _create_excluded_top(self, orig_top, positions):
        modeller = Modeller(orig_top, positions)
        res_to_delete = []
        for res in orig_top.residues():
            if res.name in self._exclude_residues:
                res_to_delete.append(res)
            else:
                for atom in res.atoms():
                    self._keep_atom_idx.append(atom.index)
        modeller.delete(res_to_delete)
        return modeller.topology


    def _callback(self, pos, inc_step_count=True, override=False):
        if (self._step_num % self._print_interval == 0) or override:
            if not self._silent:
                if self._print_header:
                    print("\n {:>5s}  {:>15s}  {:>15s}".format("Step", "Func (kJ/mol)", "Max Force"))
                    self._print_header = False

                force_norms = np.linalg.norm(self._current_forces, axis=1)
                max_force_norm = force_norms.max()
                output_line = " {:5d}  {:15.3f}  {:15.3f}".format(self._step_num, self._current_energy, max_force_norm)
                if self._other_to_report is not None:
                    for val in self._other_to_report:
                        output_line += " {:15.5E}".format(val)
                #output_line += " {:15d}".format(np.argmax(force_norms))
                print(output_line)
            
            if self._out_file is not None:

                all_pos = self._current_all_pos.reshape(-1,3)*uu.nanometer
                if self._exclude_residues is not None:
                    if self._excluded_top is None:
                        self._excluded_top = self._create_excluded_top(self._topology, all_pos)
                    out_pos = all_pos[self._keep_atom_idx]
                    out_top = self._excluded_top
                else:
                    out_pos = all_pos
                    out_top = self._topology

                PDBFile.writeModel(out_top, out_pos, file=self._out_file, modelIndex=self._step_num)
                if hasattr(self._out_file, 'flush') and callable(self._out_file.flush):
                    self._out_file.flush()

        if inc_step_count:
            self._step_num += 1


    def minimize(self, tolerance=10*uu.kilojoules_per_mole/uu.nanometer, maxIterations=500, output_interval=1):
        self._context.applyConstraints(self._working_constraint_tol)
        gtol = tolerance
        if uu.is_quantity(gtol):
            gtol /= gtol.unit
        init_state = self._context.getState(getForces=True, getEnergy=True, getPositions=True)
        init_pos_all = init_state.getPositions(True).value_in_unit(uu.nanometer)
        self._current_all_pos = init_pos_all
        #self._keep_idx = list(range(len(init_pos_all)))
        init_pos = init_pos_all[self._keep_idx].flatten()
        init_energy, init_forces = self._target_func(init_pos)
        force_norms = [np.linalg.norm(f) for f in init_forces]


        print(" Number of total coordinates:        %d" % len(init_pos_all))
        print(" Number of coordinates to optimize:  %d" % len(self._keep_idx))
        print(" Number of Degrees of Freedom:       %d" % (len(self._keep_idx)*3))
        print(" Initial max. force: {:15.3f} kJ/mol".format(np.max(force_norms)))
        print(" Initial energy:     {:15.3f} kJ/mol/nm".format(init_energy))


        # Repeatedly minimize, steadily increasing the strength of the springs until all constraints are satisfied.

        self._step_num = 0
        prevMaxError = 1e10
        if output_interval:
            self._print_interval = output_interval
            
        start_all = time.time()
        while(True):
            self._print_header = True
            res = minimize(self._target_func, init_pos, 
                method='l-bfgs-b', 
                jac=True, 
                callback=self._callback,
                options={'maxiter': maxIterations, 'disp': False, 'gtol': gtol})

            print(" SciPy L-BFGS-B Solver says: \n         ", res.message)
            # if res.success == False or np.max(np.abs(res.jac)) > tolerance:
            #     print(" Solver failed fo find minimum")
            #     print(" Falling back to Steepest descent algorithm")

            self._callback(res.x, inc_step_count=False, override=True)

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
                init_pos = np.copy(init_pos_all[self._keep_idx]).flatten()
            else:
                print(" Increasing constraint cost. Current Max error = {:15.10f} nm".format(maxError))
                init_pos = res.x
        end_all = time.time()
        self._total_time += (end_all - start_all)

        final_pos = res.x.reshape(-1,3)
        final_energy, final_forces = self._target_func(final_pos)
        force_norms = [np.linalg.norm(f) for f in final_forces]
        print(" Final max. force:   {:15.3f} kJ/mol/nm".format(np.max(force_norms)))
        print(" Final energy:       {:15.3f} kJ/mol".format(final_energy))
        print(" Total time:                     {:10.2f} s".format(self._total_time))
        print(" Total Force and Energy time:    {:10.2f} s".format(self._force_time))
        if self._out_file is not None:
            self._out_file.close()

        # If necessary, do a final constraint projection to make sure they are satisfied
        # to the full precision requested by the user.
        if (self._constraint_tol < self._working_constraint_tol):
            self._context.applyConstraints(self._working_constraint_tol)

    def _test_forces(self, fun, x0):

        #self.k = 0
        f0, g0 = fun(x0)
        grads = []
        h = (np.finfo(np.float64).eps)**(1/3)*10
        h = 0.1*0.0001
        for n in range(len(x0)):
            xp = x0.copy()
            xp[n] += h
            xn = x0.copy()
            xn[n] -= h

            fp = fun(xp)[0]
            fn = fun(xn)[0]
            num_deriv = (fp - fn)/(2*h)
            grads.append(num_deriv)
            print(" {:6d}  {:15.10f}  {:15.10f}".format(n, num_deriv, g0[n]))
        grads = np.reshape(grads, (-1, 3))
        exit()

        self._context.setPositions(np.reshape(x0, (-1, 3)))
        pos = self._context.getState(getPositions=True).getPositions(True).value_in_unit(uu.nanometers)[self._keep_idx]
        g0 = -self._context.getState(getForces=True).getForces(True).value_in_unit(uu.kilojoule_per_mole/uu.nanometers)
        h = (np.finfo(np.float64).eps)**(1/3)*10
        h = 0.1*0.0001
        pos0 = np.array(pos)
        for n in range(len(pos)):
            num_deriv = [0, 0, 0]
            for i in range(3):
                xp = pos0.copy()
                xp[n, i] += h
                xn = pos0.copy()
                xn[n, i] -= h

                self._context.setPositions(xp)
                energy_p = self._context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(uu.kilojoule_per_mole)
                
                self._context.setPositions(xn)
                energy_n = self._context.getState(getEnergy=True).getPotentialEnergy().value_in_unit(uu.kilojoule_per_mole)


                num_deriv[i] = (energy_p - energy_n)/(2*h)
            print(" {:6d}  [ {:10.5} {:10.5} {:10.5} ]  [ {:10.5} {:10.5} {:10.5} ]".format(n, *tuple(num_deriv), *tuple(g0[n])))

        exit()

    def _grad_descent(self, fun, init_x, tolerance, callback, maxiters=500):
        
        max_line_search = 15
        line_searches = 0
        step_size = 1e-6
        message = ""
        n_iters = 1

        f0, g0 = fun(init_x)
        run = True

        x0 = np.copy(init_x)
        while True:
            if not run:
                break
            line_searches = 0
            step_size0 = step_size
            while True:
                x = x0 - g0*step_size
                f, g = fun(x)
                if f < f0:
                    f0 = f
                    x0 = x.copy()
                    g0 = g.copy()
                    break
                else:
                    step_size *= 0.1
                    line_searches += 1
                    if step_size < 1e-10 and (f - f0) < 10:
                        print("Proceeding ahead anyway ", f - f0, step_size0)
                        f0 = f
                        x0 = x.copy()
                        g0 = g.copy()
                        step_size = 1e-6
                        break
                    if line_searches == max_line_search:
                        run = False
                        message = "Maximum line searches reached"
                        break
            
            step_size *= 2
            n_iters += 1
            self._other_to_report = [step_size]
            callback(x0)
            if np.max(np.abs(g)) < tolerance:
                message = "Tolerance reached"
                run = False
            elif n_iters == maxiters:
                message = "Maximum iterations achieved"
                run = False

        success = not("Maximum" in message)
        res = OptimizeResult(x=x, success=success, fun=f, jac=g, nit=n_iters,
                            message=message)
        return res
            

    def _profile(self, pos):
        for n in range(5):
            self._target_func(pos)

    def _target_func(self, pos):

        new_pos = self._context.getState(getPositions=True).getPositions(True).value_in_unit(uu.nanometer)
        new_pos[self._keep_idx] = np.reshape(pos, (-1, 3))

        # calculate new virtual site positions
        #new_pos = np.reshape(pos, (-1, 3))
        self._context.setPositions(new_pos)
        self._context.computeVirtualSites()

        #   calculates forces and energies
        start_f = time.time()
        state = self._context.getState(getEnergy=True, getForces=True, getPositions=True)
        forces_all = state.getForces(asNumpy=True).value_in_unit(uu.kilojoule_per_mole/uu.nanometer)
        pos_all = state.getPositions(True).value_in_unit(uu.nanometer)
        self._current_all_pos = pos_all.flatten()
        energy = state.getPotentialEnergy().value_in_unit(uu.kilojoule_per_mole)
        end_f = time.time()
        self._force_time += (end_f - start_f)

        #   apply harmonic cost functions for constraints
        if len(self._constraint_distances) != 0:
            delta = pos_all[self._constraint_idx_2] - pos_all[self._constraint_idx_1]
            r = np.linalg.norm(delta, axis=1)
            delta /= r[:, None]
            dr = r - self._constraint_distances
            kdr = self._k*dr
            energy += 0.5*np.sum(kdr*dr)
            k_dr_delta = kdr[:, None]*delta
            np.add.at(forces_all, self._constraint_idx_1,  k_dr_delta)
            np.add.at(forces_all, self._constraint_idx_2, -k_dr_delta)

        #   zero out all forces not being updated
        forces_all[self._zero_force_index_list] *= 0
 
        #forces = forces_all[self._force_index_list]
        forces = forces_all[self._keep_idx]
        self._current_energy = energy
        self._current_forces = forces

        return energy, -forces.flatten()
        return energy, -forces_all.flatten()
 
from copy import deepcopy, copy
from openmm.app import forcefield as ff, Simulation, Element
import openmm as mm
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

    '''
    def registerResidue(self, element):
        res_name = element.attrib['name']
        self.residues[res_name] = {}
        for atom in element.findall('Atom'):
            res_atom = {}
            for attrib in ['frz_chg', 'frz_exp', 'dyn_exp', 'pauli_exp', 'pauli_radii']:
                res_atom[attrib] = float(atom.attrib[attrib])
            self.residues[res_atom][atom.attrib['name']] = res_atom

    def registerPatchResidue(self, element):
        res_name = element.attrib['name']
        self.residues[res_name] = {}
        residue = self.residues[res_name]
        base_res_elms = element.findall('BaseResidue')
        if len(base_res_elms) != 1:
            raise ValueError(" Exactly one BaseResidue must be specified in each PatchResidue")
        residue.update(self.residues[base_res_elms.attrib('name')])
        for atom in element.findall('AddAtom'):
            residue[atom.attrib['name']] = self.add_atom(atom)
        for atom in element.findall('RemoveAtom'):
            residue.pop(atom.attrib['name'])


    def registerPolFragment(self, element):
        frag_name = element.attrib['name']
        self.fragments[frag_name] = {}
        for res in element.findall('Residue'):
            self.res_to_fragment[res.attrib['resName']] = frag_name
        for atom in element.findall('Atom'):
            frag_atom = {}
            for attrib in ['frz_chg', 'frz_exp', 'dyn_exp', 'pauli_exp', 'pauli_radii']:
                frag_atom[attrib] = float(atom.attrib[attrib])
            self.fragments[frag_name][atom.attrib['name']] = frag_atom
    '''


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

        
        # frag_sections = element.findall('Fragment')
        # for frag in frag_sections:
        #     generator.registerPolFragment(frag)



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
                particle.dyn_chg = params['frz_exp']
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

            


        
        # for atom in []:
        #     res_name  = atom.residue.name
        #     frag_name = self.res_to_fragment[res_name]
        #     fragment = self.fragments[frag_name]
        #     if atom.name not in fragment: continue
        #     frag_atom = fragment[atom.name]
        #     particle = ParticleInfo(atom.element.atomic_number)
        #     particle.frz_chg = frag_atom['frz_chg']
        #     particle.frz_exp = frag_atom['frz_exp']
        #     particle.dyn_exp = frag_atom['dyn_exp']
        #     particle.pauli_exp= frag_atom['pauli_exp']
        #     particle.pauli_radii = frag_atom['pauli_radii']
        #     force.add_particle(particle)
            
        sys._amberFDForce = force
        sys._data = data
        sys._generator= self
        pass

    def postprocessSystem(self, sys, data, args):
        print("IN POSTPROCESS")
        # for force in sys.getForces():
        #     print(force)
        pass

ff.parsers['AmberFDForce'] = AmberFDGenerator.parseElement

class MoleculeImporter():
    def __init__(self, *files) -> None:
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

        if isinstance(files, tuple):
            pass
        else:
            files = tuple([files])

        topology = None
        positions = None
        forcefield = None
        file_extensions = []
        for file in files:
            if not isinstance(file, str):
                raise ValueError("files must be of type String or a tuple of Strings")
            ext = splitext(file)[-1]
            file_extensions.append(ext)
            if ext == '.pdb':
                pdb = PDBFile(file)
                topology  = copy(topology)
                positions = copy(pdb.positions)
            elif ext == '.gro':
                gro = GromacsGroFile(file)
                positions = gro.positions
            elif ext == '.top':
                topology = copy(GromacsTopFile(file).topology)
            elif ext == '.xml':
                forcefield = ff.ForceField(file)
            else:
                raise ValueError("Valid file extensions are .pdb, .gro, .top, and .xml")
            
        if forcefield is None:
            raise ValueError("No .xml forcefield was provided in *files list")
        if positions is None:
            raise ValueError("A molecule file with positions must be provided in *files list")
        if topology is None:
            raise ValueError("A file with a molecular topology (.pdb or .top) must be provided in the *files list")
                
        modeller = Modeller(topology, positions)
        modeller.addExtraParticles(forcefield)
        self.positions = modeller.getPositions()
        self.topology = modeller.getTopology()
        self.forcefield = forcefield
        

class AmberFDSimulation(Simulation):
    def __init__(self, topology, system, integrator, platform=None, platformProperties=None, state=None):
        super().__init__(topology, system, integrator, platform=platform, platformProperties=platformProperties, state=state)


class BFGS(object):
    def __init__(self, context, constraints=None, out_pdb=None, topology=None):
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
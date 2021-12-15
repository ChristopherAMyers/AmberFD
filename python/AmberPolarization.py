from openmm.app import forcefield as ff, Simulation
import openmm as mm
import sys
from os.path import *
import numpy as np
#sys.path.insert(1, '/network/rit/lab/ChenRNALab/awesomeSauce/code/fluctuating_density/build')
sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
from AmberFD import AmberFD, FlucDens, VectorI, VectorD, VectorPairII, ParticleInfo, PairII

class _DispersionParamAtom:
    C6 = 0
    type = None
    vdw_radii = 0

class _PolarizationPairwiseAtom:
    type = None
    exp = 1.0
    radii = 0.0

class _PolarizationFragment:
    residues = []
    atoms = []

class _AtomParams:
    frz_exp = 1.0
    dyn_exp = 1.0
    frz_chg = 0.0
    pauli_exp= 1.0
    pauli_radii = 0.0
    name = None

class AmberFDGenerator(object):

    def __init__(self, forcefield) -> None:
        self.dispersion_params = {}
        self.dispersion_atoms = []
        self.polPairwise_atoms = []
        self.fragments = {}
        self.res_to_fragment = {}
        
    def registerDispersion(self, element):
        for param in ['s6', 'a1', 'a2']:
            self.dispersion_params[param] = element.attrib[param]
        for atom in element.findall('Atom'):
            # disp_atom = _DispersionParamAtom
            # disp_atom.type = atom.attrib['type']
            # disp_atom.C6 = atom.attrib['C6']
            # disp_atom.vdw_radii = atom.attrib['vdw-radii']
            disp_atom = {}
            disp_atom['type'] = atom.attrib['type']
            for attrib in ['C6', 'vdw-radii']:
                disp_atom[attrib] = float(atom.attrib[attrib])

            self.dispersion_atoms.append(disp_atom)

    def registerPolPairwise(self, element):
        for atom in element.findall('Atom'):
            # pol_atom = _PolarizationPairwiseAtom()
            # pol_atom.exp = atom.attrib['exp']
            # pol_atom.radii = atom.attrib['radii']
            # pol_atom.type = atom.attrib['type']
            pol_atom = {}
            pol_atom['type'] = atom.attrib['type']
            for attrib in ['exp', 'radii', 'type']:
                pol_atom[attrib] = atom.attrib[attrib]
            self.polPairwise_atoms.append(pol_atom)

    def registerPolFragment(self, element):
        frag_name = element.attrib['name']
        self.fragments[frag_name] = {}
        for res in element.findall('Residue'):
            self.res_to_fragment[res.attrib['resName']] = frag_name
        for atom in element.findall('Atom'):
            # frag_atom = _AtomParams
            # frag_atom.frz_chg = atom.attrib['frz_chg']
            # frag_atom.dyn_exp = atom.attrib['dyn_exp']
            # frag_atom.frz_exp = atom.attrib['frz_exp']
            # frag_atom.pauli_exp = atom.attrib['pauli_exp']
            # frag_atom.pauli_radii = atom.attrib['pauli_radii']
            # frag_atom.name = atom.attrib['name']

            frag_atom = {}
            for attrib in ['frz_chg', 'frz_exp', 'dyn_exp', 'pauli_exp', 'pauli_radii']:
                frag_atom[attrib] = float(atom.attrib[attrib])
            self.fragments[frag_name][atom.attrib['name']] = frag_atom


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

        disp_sections = element.findall('Dispersion')
        if len(disp_sections) > 1:
            raise ValueError("More than one <Dispersion> section specified in force field")
        if len(disp_sections) == 0:
            raise ValueError(" <Dispersion> section not found in <AmberFDForce>")
        generator.registerDispersion(disp_sections[0])

        pol_pairwise_sections = element.findall('PolarizationPairwise')
        if len(pol_pairwise_sections) > 1:
            raise ValueError("More than one <PolarizationPairwise> section specified in force field")
        generator.registerPolPairwise(pol_pairwise_sections[0])

        frag_sections = element.findall('Fragment')
        for frag in frag_sections:
            generator.registerPolFragment(frag)



    def createForce(self, sys, data, nonbondedMethod, nonbondedCutoff, args):
        print("IN CREATE FORCE")
        force = AmberFD(5)
        for atom in data.atoms:
            res_name  = atom.residue.name
            frag_name = self.res_to_fragment[res_name]
            fragment = self.fragments[frag_name]
            if atom.name not in fragment: continue
            frag_atom = fragment[atom.name]
            particle = ParticleInfo(atom.element.atomic_number)
            particle.frz_chg = frag_atom['frz_chg']
            particle.frz_exp = frag_atom['frz_exp']
            particle.dyn_exp = frag_atom['dyn_exp']
            particle.pauli_exp= frag_atom['pauli_exp']
            particle.pauli_radii = frag_atom['pauli_radii']
            force.add_particle(particle)
            


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

class AmberFDSimulation(Simulation):
    def __init__(self, topology, system, integrator, platform=None, platformProperties=None, state=None):
        super().__init__(topology, system, integrator, platform=platform, platformProperties=platformProperties, state=state)

        
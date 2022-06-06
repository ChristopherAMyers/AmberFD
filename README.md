# AmberFD

This package incorporates a Fluctuating Density based model for MD simulations performed with OpenMM.
Similar to a fluctuating charge model, this model incorporates explic polarization into a simulation by
adjusting the density sites at each time step so that the electrostatic energy is minimized. Additionally, 
The forces incorporated here uses atom centered monopole densities, described by Slater functions for the 
electrons and point-particles for the nuclei, to encorporate effects of electron overlap between atomic sites.

This model replaces the nonbonded interactions for RNA nucleobases only. Everything else that is simulated in
an OpenMM system remains the same. In order to construct a simple, traditional openmm simulation, one may
write a simple python script as follows (from their user guide: http://docs.openmm.org/latest/userguide/application/02_running_sims.html)

```Python
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

pdb = PDBFile('input.pdb')
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(pdb.topology, system, integrator)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)
```

To enabme AmberFD forces, a similar python script can be use with very minor adjustments:

```Python
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from AmberFD.AmberFD import MoleculeImporter, AmberFDSimulation

### AmberFD Handles forcefield assignment and importing structres
pdb = MoleculeImporter('input.pdb', ('AmberFD.xml', 'amber14/tip3pfb.xml'), onbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds) 
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = AmberFDSimulation(pdb.topology, system, integrator)   ### Similar to Simulation(), but enables AmberFD forces 
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
simulation.step(10000)
```

Full documentation (work in progress) for all the methods, arguments, and tools available
in AmberFD can be found in it's separate documentation page: https://christopheramyers.github.io/AmberFD_Documentation/html/index.html

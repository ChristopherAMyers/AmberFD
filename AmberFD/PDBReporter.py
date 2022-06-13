from openmm.app import PDBFile, PDBReporter, Modeller, modeller
from openmm.openmm import NonbondedForce


class PDBReporterExtra(PDBReporter):
    ''' Extension to the PDB Reporter that can exclude residues to be exported '''
    def __init__(self, file, reportInterval, enforcePeriodicBox=None, excludeResidueNames=None):
        ''' Create a PDBReporter that can be added to an OpenMM Simulation.reporters list. This
        is simply an extension to the OpenMM PDBReporter with the added option to exclude
        particule residue names from the trajectory. For example, you do not wish to print solvent
        to the .pdb file, then excluded residue names might include ['HOH', Cl', K']. Any residue named
        in the provided excluded residue list will not be exported to the trajectory file.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        excludeResidueNames: list or string
            Residue names to exclude from the trajectory file. Any residue named
        in the provided excluded residue list will not be exported.
            
        '''
        super().__init__(file, reportInterval, enforcePeriodicBox=enforcePeriodicBox)

        self._exclude_residues = None
        if excludeResidueNames is not None:
            if isinstance(excludeResidueNames, str):
                self._exclude_residues = tuple([excludeResidueNames])
            else:
                self._exclude_residues = tuple(excludeResidueNames)
        self._keep_atom_idx = []

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

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation: Simulation
            The Simulation to generate a report for
        state: State
            The current state of the simulation
        """

        if self._nextModel == 0:
            if self._exclude_residues is not None:
                self._topology = self._create_excluded_top(simulation.topology, state.getPositions())
            else:
                self._topology = simulation.topology
            PDBFile.writeHeader(self._topology, self._out)
            self._nextModel += 1

        pos = state.getPositions(True)
        if self._exclude_residues is not None:
            pos = pos[self._keep_atom_idx]
        PDBFile.writeModel(self._topology, pos, self._out, self._nextModel)
        self._nextModel += 1
        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()
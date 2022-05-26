from openmm.app import PDBFile, PDBReporter, Modeller, modeller
from openmm.openmm import NonbondedForce


class PDBReporterExtra(PDBReporter):
    """ Extension to the PDB Reporter """
    def __init__(self, file, reportInterval, enforcePeriodicBox=None, excludeResidueNames=None):
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
        simulation : Simulation
            The Simulation to generate a report for
        state : State
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
from openmm.app import PDBFile, PDBReporter, Modeller, modeller
from openmm.openmm import NonbondedForce


class PDBReporterExtra(PDBReporter):
    def __init__(self, file, reportInterval, enforcePeriodicBox=None, excludeResidueNames=None):
        super().__init__(file, reportInterval, enforcePeriodicBox=enforcePeriodicBox)

        self._exclude_residues = None
        if self._exclude_residues is not None:
            self._exclude_residues = tuple(excludeResidueNames)

    def _create_excluded_top(self, orig_top, positions):
        modeller = Modeller(orig_top, positions)
        for res in orig_top.residues():
            if res.name in self._exclude_residues:
                modeller.delete(res)
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

        PDBFile.writeModel(self._topology, state.getPositions(), self._out, self._nextModel)
        self._nextModel += 1
        if hasattr(self._out, 'flush') and callable(self._out.flush):
            self._out.flush()
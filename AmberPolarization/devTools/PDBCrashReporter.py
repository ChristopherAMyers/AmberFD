import numpy as np
from openmm.app import PDBFile

class PDBCrashReporter(object):

    def __init__(self, file_loc, steps_to_save):

        self._n_steps_to_save = steps_to_save
        self._top = None
       
        self._file_loc = file_loc
        self._prev_index = 0
        self._cached_pos = None
        

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A five element tuple. The first element is the number of steps
            until the next report. The remaining elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.
        """
        #steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (1, True, False, False, False)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """

        if self._top is None:
            n_atoms = simulation.topology.getNumAtoms()
            self._cached_pos = [np.zeros((n_atoms, 3)) for x in range(self._n_steps_to_save)]
            self._top = simulation.topology

        state_index = state.getStepCount() % self._n_steps_to_save
        self._prev_index = state_index
        self._cached_pos[state_index] = state.getPositions(True)

    def create_crash_file(self):

        pos_to_save = []
        start = self._prev_index + 1
        end = start + self._n_steps_to_save + 1
        for n in range(start, end):
            pos_to_save.append(self._cached_pos[n % self._n_steps_to_save])

        with open(self._file_loc, 'w') as file:
            for n, pos in enumerate(pos_to_save):
                PDBFile.writeModel(self._top, pos, file, keepIds=True, modelIndex=n)
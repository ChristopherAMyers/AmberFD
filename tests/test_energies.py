import sys
from os.path import *
sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
from _AmberFD import Energies


eng1 = Energies()
eng2 = Energies()

eng1.pauli      = 1.0
eng1.disp       = 2.0
eng1.frz        = 3.0
eng1.pol        = 4.0
eng1.vct        = 5.0
eng1.elec_elec  = 6.0
eng1.elec_nuc   = 7.0
eng1.nuc_nuc    = 8.0
eng1.pauli_wall = 9.0

eng2.pauli      = 1.2
eng2.disp       = 2.2
eng2.frz        = 3.2
eng2.pol        = 4.2
eng2.vct        = 5.2
eng2.elec_elec  = 6.2
eng2.elec_nuc   = 7.2
eng2.nuc_nuc    = 8.2
eng2.pauli_wall = 9.2


import sys
from os.path import *
sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
from .AmberFD import *
from _AmberFD import FlucDens, DispersionPauli, AmberFD, Energies, ParticleInfo

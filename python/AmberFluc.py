import sys
from os.path import *
import numpy as np
sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
from AmberFD import AmberFD, FlucDens, VectorI, VectorD, VectorPairII, ParticleInfo, PairII

class AmberFFParams():
    def __init__(self, file_loc) -> None:
        with open(file_loc) as file:
            pass
import sys
from os.path import *
sys.path.insert(1, join(dirname(realpath(__file__)), '../build/'))
print("SWIG_IMPORT: ", abspath(curdir))
from FlucDens import FlucDens


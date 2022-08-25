# python version of the gpu code:
# *** requirements ***
# python3
# numba: https://numba.readthedocs.io/en/stable/index.html
# uproot: https://uproot.readthedocs.io/en/latest/


#!/usr/bin/python

import sys
import os
import datetime

from numba import jit
import numpy as np
import time
import uproot

if len(sys.argv) < 2:
    print("Supply input file name")
    sys.exit()

filename = sys.argv[1]
file = uproot.open(filename)

#I need the actual list of hit(s) though... and the event and hit class, and the geometry.
@jit
#(nopython=True)
def reco_tracklet_in_stations(stID, listID, pos_exp, window):
    


@jit
def reco_backtracks()



@jit
def reco_globaltracks()



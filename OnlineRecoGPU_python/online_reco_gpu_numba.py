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


# functions
#I need the actual list of hit(s) though... and the event and hit class, and the geometry.
@jit(nopython=True)
def reco_tracklet_in_stations(stID, listID, *pos_exp): #doesn't seem to take another array 
# - makes pairs of hits in xx', uu', vv', in selected station;
# - if a view doesn't have hits, stop here.
# - combination of hits to form tracklets: 
#   * loop on x hits: combine x with u: to each x can only corresponds a range in u_min<u_pos<u_max
#   * inside loop x, loop on u hits: reject all u hits which do not meet u_min<u_pos<u_max:
#     for those who do, calculate v_window, v_min, v_max
#   * inside loop u, loop on v hits: reject all v hits which do not meet v_min<v_pos<v_max:
#    for those who do, add a tracklet with the combination of hits, and fit it;
#   * if tracklet is "valid" (see below) it is kept, otherwise it isn't
# - Once the combinations have been made, tracklets are added intot the tracklet list
    print('reco_tracklet_in_station', stID)
#TODO: implement the whole function

@jit(nopython=True)
def reco_backtracks():
# - combination of tracklets from station 2 and 3 to form backtracks
#   * loop on station 3 tracklets; if not coarse mode, loop on the tracklet 3 hits to extract only the X hits
#   * inside loop 3, loop on station 2 tracklets; if not coarse mode, loop on the tracklet 2 hits to extract only the X hits...
#     then fit the backtrack in X; then check the proportional tubes: we want at least one hit there;
#     otherwise, add the two tracklets together to obtain tracklet 23 (aka backtrack), and fit it.
#     If the fit chi^2 is too high, reject tracklet; if not coarse mode, resolve left right for backtrack.
#     Then keep only the best backtrack (i.e. with best chi^2 or best proba)
    print("reco_backtracks")
#TODO: implement the whole function


@jit(nopython=True)
def reco_globaltracks():
# - Combinations of backtracks and hits in station 1 to form global tracks:
#   * Loop on backtracks: evaluation of windows with Sagitta method if KMag ON, with extrapolation otherwise;
#     then build tracklets in station 1 (using the windows obtained with the search window method);
#   * inside loop on backtracks, loop on station plane (2 stations)
#   * inside loop on station plane, loop on station  hits; multiply (?) tracklet 1 and backtrack, and fit; reject if no hodo hits;
#     if not coarse mode, resolve left right for backtrack, then remove bad hits (on what criteria?);
#     then keep the global track with the best fit; If Kalmann filter, reconstruct vertex, keep the track with the best vertex chi^2;
#   * after the loop on station1 hits, keep the very best track: The selection logic is, prefer the tracks with best p-value, as long as it's not low-pz;
#     otherwise select the one with best vertex chisq; then fall back to the default only choice
#   * After the loop on backtracks, if best track from each station have momentum less than a defined value, merge tracks;
#     if the merged track is is better than the separate ones, keep it, otherwise, keep the best one of the two (better = with best chi^2 or best proba)
    print("reco_globaltracks")
#TODO: implement the whole function





#main 
if len(sys.argv) < 2:
    print("Supply input file name")
    sys.exit()

#opening the input file 
filename = sys.argv[1]
file = uproot.open(filename)

# following the same order as in KalmanFastTracking.cxx
reco_tracklet_in_stations(3, 1)
reco_tracklet_in_stations(4, 1)
reco_tracklet_in_stations(5, 1)
reco_backtracks()
reco_globaltracks()

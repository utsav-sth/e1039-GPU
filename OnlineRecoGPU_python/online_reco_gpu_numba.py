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

# this code needs:
# hit containers;
hits = (1)#placeholder #TODO: define those
# track containers;
tracklets = (1)#placeholder #TODO: define those


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
    #it would look like all that stuff has to be defined here and cannot be defined globally???
    z_plane_x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]#placeholder: TODO: replace with actual values
    z_plane_u = [1.1, 2.1, 3.1, 4.1, 5.1, 6.1]#placeholder: TODO: replace with actual values
    z_plane_v = [1.2, 2.2, 3.2, 4.2, 5.2, 6.2]#placeholder: TODO: replace with actual values
    ucostheta = np.cos(10./180.*3.141592654)
    usintheta = np.sin(10./180.*3.141592654)
    uwin = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]#placeholder: TODO: replace with actual values
    vcostheta = np.cos(-10./180.*3.141592654)
    vsintheta = np.sin(-10./180.*3.141592654)
    vwin = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]#placeholder: TODO: replace with actual values
    spacingplane = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
    TXMAX = 1.
    TYMAX = 1.
    
    print('reco_tracklet_in_station', stID)
    hitpairs_in_x = ((1, 4), (2, -1))#placeholder
    hitpairs_in_u = ((2, 5), (4, -1))#placeholder
    hitpairs_in_v = ((3, 6), (6, -1))#placeholder
    nhitsx = len(hitpairs_in_x)
    nhitsu = len(hitpairs_in_u)
    nhitsv = len(hitpairs_in_v)
    detID = 10#placeholder
    if nhitsx==0 or nhitsu==0 or nhitsv==0:
        return 0
#TODO: implement the whole function
    for i in range (0, nhitsx):
        xpos = (hitpairs_in_x[i][0]+hitpairs_in_x[i][1])*0.5
        if hitpairs_in_x[i][1]<0:
            xpos = hitpairs_in_x[i][0]
        umin = xpos*ucostheta-uwin[stID]
        umax = umin+2*uwin[stID]
        for j in range (0, nhitsu):
            upos = (hitpairs_in_u[j][0]+hitpairs_in_u[j][1])*0.5
            if hitpairs_in_u[j][1]<0:
                upos = hitpairs_in_u[j][0]
            if upos<umin or upos>umax:
                continue
            z_x = z_plane_x[stID]
            z_u = z_plane_u[stID]
            z_v = z_plane_v[stID]
            vwin1 = spacingplane[detID]*2*ucostheta
            vwin2 = abs((z_u+z_v-2*z_x)*ucostheta*TXMAX)
            vwin3 = abs((z_v-z_u)*usintheta*TYMAX)
            vwin = vwin1+vwin2+vwin3+2*spacingplane[detID]
            vmin = 2*xpos*ucostheta-upos-vwin
            vmax = vmin+2*vwin
            for k in range (0, nhitsv):
                upos = (hitpairs_in_v[k][0]+hitpairs_in_v[k][1])*0.5
                if hitpairs_in_v[k][1]<0:
                    vpos = hitpairs_in_v[k][0]
                if vpos<vmin or vpos>vmax:
                    continue
                #build tracklet here... but I need to define it first... (also the hits)
           

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

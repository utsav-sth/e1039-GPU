# python version of the gpu code:
# *** requirements ***
# python3
# numba: https://numba.readthedocs.io/en/stable/index.html
# uproot: https://uproot.readthedocs.io/en/latest/


#!/usr/bin/python

import sys
import os
import datetime
import time

import math
import numpy as np
import uproot
import numba
#from numba import jit
from numba import cuda

# this code needs:
# hit containers;
# hits = (1)#placeholder #TODO: define those
# track containers;
# tracklets = (1)#placeholder #TODO: define those


# functions
@numba.jit(nopython=True)
def make_hitpairs_in_station(stID, projID, detectorid, pos):
    stID-=1
# makes pairs of hits: check all hits in both detector IDs corresponding to station ID 
    spacingplane = [0., 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 1.3, 1.3, 1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 4.0, 4.0, 7.0, 7.0, 8.0, 12.0, 12.0, 10.0, 3.0, 3.0, 3.0, 3.0]
    detsuperid = [[2, 1, 3], [5, 6, 4], [8, 9, 7], [11, 12, 10], [14, 15, 13], [25, 26], [24, 27]] 
    detid1 = detsuperid[stID][projID]*2
    detid2 = detsuperid[stID][projID]*2-1

    hitlist1_pos = np.zeros(len(detectorid), dtype='float32')
    hitlist2_pos = np.zeros(len(detectorid), dtype='float32')
    hitctr1 = 0
    hitctr2 = 0
    for i in range(0, len(detectorid)):
        if(detectorid[i]==detid1):
            hitlist1_pos[hitctr1] = pos[i]
            hitctr1+=1
        if(detectorid[i]==detid2):
            hitlist1_pos[hitctr2] = pos[i]
            hitctr2+=1
    maxsize = (hitctr1+1)*(hitctr2+1)
    hitpairs = np.zeros( (maxsize, 2), dtype='float32' )
    index1 = -1
    index2 = -1
    hitflag1 = np.zeros(hitctr1, dtype='int8')
    hitflag2 = np.zeros(hitctr2, dtype='int8')
    indexpair = 0
    for i in range(0, hitctr1):
        index1+=1
        index2 = -1
        for j in range(0, hitctr2):
            index2+=1
            if(hitlist1_pos[i]-hitlist2_pos[j]>spacingplane[detsuperid[stID][projID]]):
                continue
            hitpairs[indexpair][0] = hitlist1_pos[i]
            hitpairs[indexpair][1] = hitlist2_pos[j]
            indexpair+=1
            hitflag1[index1] = 1
            hitflag2[index2] = 1
    index1 = 0;
    for i in range(0, hitctr1):
        if(hitflag1[index1]<1):
            hitpairs[indexpair][0] = hitlist1_pos[i]
            hitpairs[indexpair][1] = -1
            indexpair+=1
        index1+=1
    index2 = 0;
    for j in range(0, hitctr2):
       if(hitflag2[index2]<1):
            hitpairs[indexpair][0] = -1
            hitpairs[indexpair][1] = hitlist2_pos[j]
            indexpair+=1
            index2+=1
    return hitpairs
    

@numba.jit(nopython=True)
def reco_tracklet_in_station(stID, detectorid, hitpairs_in_x, hitpairs_in_u, hitpairs_in_v, *pos_exp): #doesn't seem to take another array 
# - takes pairs of hits in xx', uu', vv', in selected station;
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
    stID-=1
    z_plane_x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]#placeholder: TODO: replace with actual values
    z_plane_u = [1.1, 2.1, 3.1, 4.1, 5.1, 6.1]#placeholder: TODO: replace with actual values
    z_plane_v = [1.2, 2.2, 3.2, 4.2, 5.2, 6.2]#placeholder: TODO: replace with actual values
    ucostheta = np.cos(10./180.*3.141592654)
    usintheta = np.sin(10./180.*3.141592654)
    uwin = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]#placeholder: TODO: replace with actual values
    vcostheta = np.cos(-10./180.*3.141592654)
    vsintheta = np.sin(-10./180.*3.141592654)
    vwin = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1]#placeholder: TODO: replace with actual values
    spacingplane = [0., 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 1.3, 1.3, 1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 4.0, 4.0, 7.0, 7.0, 8.0, 12.0, 12.0, 10.0, 3.0, 3.0, 3.0, 3.0]
    TXMAX = 1.
    TYMAX = 1.
    
    #print('reco_tracklet_in_station', stID)
    #hitpairs_in_x = make_hitpairs_in_station(stID, 0, detectorid, pos)
    #hitpairs_in_u = make_hitpairs_in_station(stID, 1, detectorid, pos)
    #hitpairs_in_v = make_hitpairs_in_station(stID, 2, detectorid, pos)
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
           

@numba.jit(nopython=True)
def reco_backtracks():
# - combination of tracklets from station 2 and 3 to form backtracks
#   * loop on station 3 tracklets; if not coarse mode, loop on the tracklet 3 hits to extract only the X hits
#   * inside loop 3, loop on station 2 tracklets; if not coarse mode, loop on the tracklet 2 hits to extract only the X hits...
#     then fit the backtrack in X; then check the proportional tubes: we want at least one hit there;
#     otherwise, add the two tracklets together to obtain tracklet 23 (aka backtrack), and fit it.
#     If the fit chi^2 is too high, reject tracklet; if not coarse mode, resolve left right for backtrack.
#     Then keep only the best backtrack (i.e. with best chi^2 or best proba)
    #print("reco_backtracks")
#TODO: implement the whole function
    return 0


@numba.jit(nopython=True)
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
    #print("reco_globaltracks")
#TODO: implement the whole function
    return 0




#main 
from timeit import default_timer as timer
start_all=timer()
if len(sys.argv) < 2:
    print("Supply input file name")
    sys.exit()


start_evload=timer()
#opening the input file 
filename = sys.argv[1]
inputtree = uproot.open(filename+':save')
#inputdata = inputtree.arrays(library="np")
detectorid=inputtree["fAllHits.detectorID"].arrays(library="np")["fAllHits.detectorID"]
elementid=inputtree["fAllHits.elementID"].arrays(library="np")["fAllHits.elementID"]
pos=inputtree["fAllHits.pos"].arrays(library="np")["fAllHits.pos"]
nevents=len(detectorid)
if(len(sys.argv)>2):
    if(int(sys.argv[2])<nevents):
        nevents = int(sys.argv[2])
print('number of events to process',nevents)
#detector_data=np.zeros((len(detectorid),30,200),dtype='int8')
end_evload=timer()

start_proc=timer()
# we will have to parallelize the thing so let's do the loop just to understand, but don't invest too much time in it.
for ev in range(0, nevents):
    #print(inputdata['fAllHits.detectorID'][ev])#for tests
    hitpairs_in_x = make_hitpairs_in_station(3, 0, detectorid[ev],pos[ev])
    hitpairs_in_u = make_hitpairs_in_station(3, 1, detectorid[ev],pos[ev])
    hitpairs_in_v = make_hitpairs_in_station(3, 2, detectorid[ev],pos[ev])
    #reco_tracklet_in_station(3, detectorid[ev], hitpairs_in_x, hitpairs_in_u, hitpairs_in_v)
# following the same order as in KalmanFastTracking.cxx
#    reco_tracklet_in_station(3, detectorid[ev],pos[ev])
#    reco_tracklet_in_station(4, detectorid[ev],pos[ev])
#    reco_tracklet_in_station(5, detectorid[ev],pos[ev])
    reco_backtracks()
    reco_globaltracks()

end_proc=timer()

end_all=timer()
print("data loading time: ", end_evload-start_evload)
print("data processing time: ", end_proc-start_proc)
print("total time of execution: ", end_all-start_all)

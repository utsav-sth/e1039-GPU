#include "reconstruction_classes.cuh"
#include "trackanalyticalminimizer.cuh"

// --------------------------------------------------------------- //
// functions to calculate bottom and top end wire points for a hit //
// --------------------------------------------------------------- //

__device__ float x_bep(const int detid, const int elid, const gPlane* plane)
{
	return plane->p1x_w1[detid]+plane->dp1x[detid]*(elid-1);
}

__device__ float x_tep(const int detid, const int elid, const gPlane* plane)
{
	return x_bep(detid, elid, plane)+plane->deltapx[detid];
}

__device__ float y_bep(const int detid, const int elid, const gPlane* plane)
{
	return plane->p1y_w1[detid]+plane->dp1y[detid]*(elid-1);
}

__device__ float y_tep(const int detid, const int elid, const gPlane* plane)
{
	return y_bep(detid, elid, plane)+plane->deltapy[detid];
}

__device__ float z_bep(const int detid, const int elid, const gPlane* plane)
{
	return plane->p1z_w1[detid]+plane->dp1z[detid]*(elid-1);
}

__device__ float z_tep(const int detid, const int elid, const gPlane* plane)
{
	return z_bep(detid, elid, plane)+plane->deltapz[detid];
}



__device__ float position(const float pos, const float drift, const short sign)
{
	return (pos+sign*drift);
}


// ---------------------------------------------------------------- //
// functions to return track state at a certain z
// ---------------------------------------------------------------- //

__device__ float x_trk(const float x0, const float tx, const float z)
{
	return z*tx+x0;
}

__device__ float err_x_trk(const float errx0, const float errtx, const float z)
{
	return z*errtx+errx0;
}

__device__ float y_trk(const float y0, const float ty, const float z)
{
	return z*ty+y0;
}

__device__ float err_y_trk(const float erry0, const float errty, const float z)
{
	return z*errty+erry0;
}


// --------------------------------------------------------------- //
// calculation of y for UV hits
// --------------------------------------------------------------- //

__device__ bool calculate_y_uvhit(const int detid, const int elid, const float drift, const short hitsign, const float x0, const float tx, const gPlane* planes, float &y, float &err_y){
	float p1x = x_bep(detid, elid, planes);
	float p2x = p1x+planes->deltapx[detid];//must win some time to no try to access the same plane variables twice//x_tep(detid, elid, planes);
	
	float x_trk = x0+planes->z[ detid ]*tx;

	if( x_trk-drift>p1x && x_trk-drift>p2x)return false;// if xtrk>p1x and >p2x, no overlap possible
	if( x_trk+drift<p1x && x_trk+drift<p2x)return false;// if xtrk<p1x and <p2x, no overlap possible

	float p1y = y_bep(detid, elid, planes);

	y = p1y + (x_trk-p1x) *  planes->deltapy[ detid ]/planes->deltapx[ detid ];
	if(hitsign!=0)//y += hitsign*drift*planes->sintheta[detid];
		y+= copysign(1.f, planes->sintheta[detid])*2.f*(drift-0.5f*hitsign);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf("det %d chan %d p1x %1.4f p1y %1.4f p2x %1.4f dpy %1.4f dpx% 1.4f x_trk %1.4f y %1.4f \n", detid, elid, p1x, p1y, p2x, planes->deltapy[ detid ], planes->deltapx[ detid ], x_trk, y );
#endif
	//if hitsign is zero, we don't want to toss a hit that could potentially be in range of the track accounting for the drift distance


	y = max(y, p1y);
	y = min(y, p1y+planes->deltapy[ detid ]);
	
	err_y = planes->resolution[ detid ] * fabs(planes->deltapy[ detid ]/planes->deltapx[ detid ]);
	return true;
}

__device__ bool calculate_y_uvhit(const int detid, const int elid, const float drift, const short hitsign, const float x0, const float tx, const gPlane* planes, float &z, float &y, float &err_y){
	float p1x = x_bep(detid, elid, planes);
	float p1y = y_bep(detid, elid, planes);
	float p1z = z_bep(detid, elid, planes);

	float dpx = planes->deltapx[ detid ];
	float dpz = planes->deltapz[ detid ];

	float p2x = p1x+dpx;//must win some time to no try to access the same plane variables twice//x_tep(detid, elid, planes);

	float x_trk = x0+planes->z[ detid ]*tx;
	z = p1z + (x_trk-p1x)*dpz/dpx;

	if( x_trk-drift>p1x && x_trk-drift>p2x)return false;// if xtrk>p1x and >p2x, no overlap possible
	if( x_trk+drift<p1x && x_trk+drift<p2x)return false;// if xtrk<p1x and <p2x, no overlap possible


	y = p1y + (x_trk-p1x) *  planes->deltapy[ detid ]/dpx;
	if(hitsign!=0)//y += hitsign*drift*planes->sintheta[detid];
		y+= copysign(1.f, planes->sintheta[detid])*2.f*(drift-0.5f*hitsign);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf("det %d chan %d p1x %1.4f p1y %1.4f p2x %1.4f dpy %1.4f dpx% 1.4f x_trk %1.4f y %1.4f \n", detid, elid, p1x, p1y, p2x, planes->deltapy[ detid ], planes->deltapx[ detid ], x_trk, y );
#endif
	//if hitsign is zero, we don't want to toss a hit that could potentially be in range of the track accounting for the drift distance


	y = max(y, p1y);
	y = min(y, p1y+planes->deltapy[ detid ]);
	
	err_y = planes->resolution[ detid ] * fabs(planes->deltapy[ detid ]/planes->deltapx[ detid ]);
	return true;
}


// --------------------------------------------------------------- //
// general calculation of x, z for y 
// --------------------------------------------------------------- //

__device__ bool calculate_xz_fy(const int detid, const int elid, const float drift, const short hitsign, const float y0, const float ty, const gPlane* planes, float &x, float &z){
	float p1x = x_bep(detid, elid, planes);
	float p1y = y_bep(detid, elid, planes);
	float p1z = z_bep(detid, elid, planes);

	float dpy =  planes->deltapy[ detid ];
	float dpz =  planes->deltapz[ detid ];
	
	z = p1z - p1y * dpz/dpy;//first we assume y = 0
	float y_trk = y0+z*ty;
	x = p1x + (y_trk - p1y) * planes->deltapx[ detid ]/dpy;
	z = p1z + (y_trk - p1y) * dpz/dpy;
	
	if(hitsign!=0)x+= hitsign*drift;
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf("det %d chan %d p1x %1.4f p1y %1.4f p1z %1.4f dpy %1.4f dpz %1.4f dpx %1.4f y_trk %1.4f x %1.4f z %1.4f \n", detid, elid, p1x, p1y, p1z, dpy, planes->deltapz[ detid ], planes->deltapx[ detid ], y_trk, x, z );
#endif
	return true;
}





// ---------------------------//
// event reducer
// -------------------------- //

__device__ int event_reduction(const gHits& hitcoll, short* hitflag, const int detid, const int nhits) {
	float w_max; // max drift distance of the hit furthest from the cluster avg position // current average position of cluster * 0.9
	
	float w_min; // max drift distance of the hit closest to the cluster avg position // current average position of cluster * 0.4
	float dt_mean; // tbd
	int cluster_iAH_arr_cur; // current cluster array
	int cluster_iAH_arr_size; // cluster size i.e number of hits in cluster
	int cluster_iAH_arr[ClusterSizeMax]; // global cluster array 
	int uniqueID; // hit unique element ID
	int uniqueID_curr; // current hit unique element ID
	float tdcTime_curr; // current hit TDC time
	int iAH; // hit index 
	int nAH_reduced; // number of hits after hit quality filtering
	// int nHitsPerDetector[nDetectors+1];
	
	//if(ic[index].EventID==0){
	//	printf("evt = %d, nAH = %d: \n", ic[index].EventID, ic[index].nAH);
	//	for(int i = 1; i<=nDetectors; i++)printf(" det %d: %d;  ", i, ic[index].NHits[i]);
	//	printf("\n");
	//}

	// initialization of array size
	cluster_iAH_arr_size = 0;
	nAH_reduced = 0;
	
	// event reducing/hit filtering
	for(iAH = 0; iAH<nhits; ++iAH) {
		
		hitflag[iAH] = 1;
		// if hit not good, set its flag to 0 and continue;
		if( (int(hitcoll.flag(iAH)) & hitFlagBit(1)) == 0) {
			//if(blockIdx.x==debug::EvRef && threadIdx.x==0)printf("detid %d chan %1.4f, flag %d \n", detid, hitcoll.chan(iAH), hitflag[iAH]);
			//if(ic[index].EventID==0)printf("hit det %d Skip out-of-time...\n", detid);
			hitflag[iAH] = -1;
			continue;
		}
		uniqueID = uniqueID_curr = -1;

		// hits in DCs or Prop tubes
		if(detid < 31 || detid > 46) {
			// evaluate "unique ID"
			uniqueID = detid*1000 + hitcoll.chan(iAH);
			// compare with current unique element ID; if different, update the unique element ID and time info 
			if(uniqueID != uniqueID_curr) {
				uniqueID_curr = uniqueID;
				tdcTime_curr = hitcoll.tdc(iAH);
			}
			// if next hit and current hit belong to the same element: if detID>36 => prop tubes (reminder that hodoscpes are out of the picture in this scope), 
			// we're suppose to have one signal (a second signal would be after-pulsing)
			// if time difference between new hit and current hit is less than 80ns for DCs, it's also considered after-pulsing
			else {
#ifdef REFINED_ER 
				if(detid > 36 || ((hitcoll.tdc(iAH) - tdcTime_curr >= 0.0) && (hitcoll.tdc(iAH) - tdcTime_curr < 80.0)) || ((hitcoll.tdc(iAH) - tdcTime_curr <= 0.0) && (hitcoll.tdc(iAH) - tdcTime_curr > -80.0))) {
#endif
					hitflag[iAH] = -1;
					continue;
#ifdef REFINED_ER 
				}
				else {
					tdcTime_curr = hitcoll.tdc(iAH);
				}
#endif
			}
		}
		// declustering of hits in DCs (from CPU code, I understand this one better)
		// if there are hits in the same plane and hitting to neighboring wires, they both give redundant information: 
		if(detid <= nChamberPlanes) {
			//if(ic[index].EventID==0)printf("%d\n", cluster_iAH_arr_size);
//			printf("Decluster...\n");
			if(cluster_iAH_arr_size == ClusterSizeMax) {
				printf("Oversized cluster...\n");
			}
			// if array size is zero, start storing the hit in the array
			if(cluster_iAH_arr_size == 0) {
				cluster_iAH_arr[0] = iAH;
				++cluster_iAH_arr_size;
			} else { // otherwise

				// current hit and previous hit are *not* in same detector plane OR next hit and current hit are *not* in neighbors cells
				// we "declusterize" i.e. we remove the hit/hits which information is redundant with other hits and/or useless
				//if(ic[index].EventID==0){
//printf("hit indices: %d %d %d\n", iAH, cluster_iAH_arr[cluster_iAH_arr_size-1], cluster_iAH_arr[0]);
//printf("hit det/elem: %d, %d; %d, %d\n", detid, ic->AllHits[cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID, 
//	    	      	      	  	 hitcoll.chan(iAH), hitcoll.chan(cluster_iAH_arr[cluster_iAH_arr_size-1]));
//printf("diffs: %d, %d\n", (detid - ic->AllHits[cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID),
//	       	   	  (hitcoll.chan(iAH) - hitcoll.chan(cluster_iAH_arr[cluster_iAH_arr_size-1])));
//printf("bools: %d, %d\n", (detid != ic->AllHits[cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID),
//	       	   	  (abs(hitcoll.chan(iAH) - hitcoll.chan(cluster_iAH_arr[cluster_iAH_arr_size-1])) > 1));
	    	     	   	//}
	    	     	   	//all hits are in the same plane now
				if(abs(hitcoll.chan(iAH) - hitcoll.chan(cluster_iAH_arr[cluster_iAH_arr_size-1])) > 1 || iAH==nhits-1) {
					if(iAH==nhits-1 && abs(hitcoll.chan(iAH) - hitcoll.chan(cluster_iAH_arr[cluster_iAH_arr_size-1])) <= 1){
						cluster_iAH_arr[cluster_iAH_arr_size] = iAH;
				  		++cluster_iAH_arr_size;
					}
					
					//if(blockIdx.x==debug::EvRef && detid>24 && detid<31)printf("%d %1.0f %1.0f %d \n", detid, hitcoll.chan(cluster_iAH_arr[0]), hitcoll.chan(cluster_iAH_arr[cluster_iAH_arr_size-1]), cluster_iAH_arr_size ); 

					// if 2 hits in cluster, evaluate w_max and w_min; drift distance has to be < w_min for one of the hits, while it has to be < w_max for the other hit 
					if(cluster_iAH_arr_size == 2) {
#ifdef REFINED_ER 
						w_max = 0.9*0.5*(hitcoll.pos(cluster_iAH_arr[cluster_iAH_arr_size-1]) - hitcoll.pos(cluster_iAH_arr[0]));
						w_min = 4.0/9.0*w_max;
						
						if((hitcoll.drift(cluster_iAH_arr[0]) > w_max && hitcoll.drift(cluster_iAH_arr[cluster_iAH_arr_size-1]) > w_min) || (hitcoll.drift(cluster_iAH_arr[0]) > w_min && hitcoll.drift(cluster_iAH_arr[cluster_iAH_arr_size-1]) > w_max)) {
#endif
						//if(ic[index].EventID==0)printf("hit indices: %d %d %d\n", iAH, cluster_iAH_arr[cluster_iAH_arr_size-1], cluster_iAH_arr[0]);
							//eliminating the existing hit with the lagest drift distance
							if(hitcoll.drift(cluster_iAH_arr[0]) > hitcoll.drift(cluster_iAH_arr[cluster_iAH_arr_size-1])) {
								//if(ic[index].EventID==0)printf("1 - hit det %d elem %d Skip cluster...\n", ic->AllHits[cluster_iAH_arr[0]].detectorID, ic->AllHits[cluster_iAH_arr[0]].elementID);
								hitflag[cluster_iAH_arr[0]] = -1;
							}
							else {
								//if(ic[index].EventID==0)printf("2 - hit det %d elem %d Skip cluster...\n", ic->AllHits[cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID, hitcoll.chan(cluster_iAH_arr[cluster_iAH_arr_size-1]));
								hitflag[cluster_iAH_arr[cluster_iAH_arr_size-1]] = -1;
							}
#ifdef REFINED_ER 

						}
						// if the time difference is less than 8 ns for detectors 19 to 24 (which btw are DC3p) we remove both
						else 
#endif
						if((((hitcoll.tdc(cluster_iAH_arr[0]) - hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_size-1])) >= 0.0 && (hitcoll.tdc(cluster_iAH_arr[0]) - hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_size-1])) < 8.0) || ((hitcoll.tdc(cluster_iAH_arr[0]) - hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_size-1])) <= 0.0 && (hitcoll.tdc(cluster_iAH_arr[0]) - hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_size-1])) > -8.0)) && (detid >= 19 && detid <= 24)) {
						        //if(ic[index].EventID==0)printf("3 - hit det %d elem %d Skip cluster...\n", detid, hitcoll.chan(iAH));
							hitflag[cluster_iAH_arr[0]] = -1;
							hitflag[cluster_iAH_arr[cluster_iAH_arr_size-1]] = -1;
						}
					}
					// if 3 hits or more in cluster: we essentially discard them all;
					if(cluster_iAH_arr_size >= 3) {
						// evaluate the mean time difference;
						dt_mean = 0.0;
						for(cluster_iAH_arr_cur = 1; cluster_iAH_arr_cur < cluster_iAH_arr_size; ++cluster_iAH_arr_cur) {
							dt_mean += ((hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_cur]) - hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_cur-1])) > 0.0 ? (hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_cur]) - hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_cur-1])) : (hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_cur-1]) - hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_cur])));
						}
						dt_mean = dt_mean/(cluster_iAH_arr_size - 1);
						// if mean time difference is less than 10, that's electronic noise, so we remove them all.
						if(dt_mean < 10.0) {
						        //if(ic[index].EventID==0)printf("4 - hit det %d elem %d Skip cluster...\n", detid, hitcoll.chan(iAH));
							for(cluster_iAH_arr_cur = 0; cluster_iAH_arr_cur < cluster_iAH_arr_size; ++cluster_iAH_arr_cur) {
								hitflag[cluster_iAH_arr[cluster_iAH_arr_cur]] = -1;
							}
						}
						// otherwise, we remove them all except first and last
						else {
						        //if(ic[index].EventID==0)printf("5 - hit det %d elem %d Skip cluster...\n", detid, hitcoll.chan(iAH));
							for(cluster_iAH_arr_cur = 1; cluster_iAH_arr_cur < cluster_iAH_arr_size; ++cluster_iAH_arr_cur) {
								hitflag[cluster_iAH_arr[cluster_iAH_arr_cur]] = -1;
							}
						}
					}
					cluster_iAH_arr_size = 0;
				} else { // if hits are in the same detector and in two neighboring cells!
				        // current hit and previous hit are in same detector plane and in neighbor wires: 
				  	// we count how many hits we have in this case, until we find a hit in a different detector or in a wire that is not a neighbor to the previous hit.
				  	//if(ic[index].EventID==0)printf("h1: det %d el %d; h2 det %d el %d \n", detid, hitcoll.chan(iAH), ic->AllHits[cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID, hitcoll.chan(cluster_iAH_arr[cluster_iAH_arr_size-1]));
				  	cluster_iAH_arr[cluster_iAH_arr_size] = iAH;
				  	++cluster_iAH_arr_size;
				}
			}
		}
	}
	//end of the hit loop
	
	//array to check for double hits: if one channel already has a hit, the next hit in the same channel is suppressed.
	// better declare and initialize this array just before doing the check 
	bool chanhashit[202];
	for(short k = 0; k<202; k++)chanhashit[k] = 0;
	
	for(iAH = 0; iAH<nhits; ++iAH) {
	//for(iAH = nhits-1; iAH>=0; --iAH) {
		if(hitflag[iAH]>0) {
			//if the chan already has a hit, then set the hit flag in the same chan as not good.
			if(chanhashit[(short)hitcoll.chan(iAH)]){
				hitflag[iAH] = -1;
			}else{
				//otherwise, set the "chanhashit" flag as true, and proceed as usual 
				chanhashit[(short)hitcoll.chan(iAH)] = true;
				++nAH_reduced;
			}
		}
	}
	
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){
		for(iAH = 0; iAH<nhits; ++iAH) {
			printf("detid %d iah %d chan %1.0f time %1.4f hitflag %d \n", detid, iAH, hitcoll.chan(iAH), hitcoll.tdc(iAH), hitflag[iAH]);
		}
	}	
#endif
	

	return nAH_reduced;
}

// --------------------------------------------------------------- //
// Hit pairing functions
// --------------------------------------------------------------- //

//new function
__device__ void make_hitpairs_in_station_bins(const gHits hitcoll1, const int nhits1, const gHits hitcoll2, const int nhits2, thrust::pair<int, int>* hitpairs, int* npairs, const short bin0, const short Nbins, short* hitflag1, short* hitflag2, const int stID){
	// I think we assume that by default we want to know where we are
	//printf("stID %d projID %d bin0 %d\n", stID, projID, bin0);
	short projID = 0;
	short bin;
	const short MaxHits = globalconsts::MaxHitsProj_X;

	for(bin = bin0; bin<bin0+Nbins; bin++){
		npairs[bin-bin0] = 0;
	}
	
	//declaring arrays for the hit lists
	for(int i = 0; i<100; i++){
		hitflag1[i] = hitflag2[i] = 0;
	}
	
	//building the lists of hits for each detector plane
	//const int detid1 = globalconsts::detsuperid[stID][projID]*2;
	//const int detid2 = globalconsts::detsuperid[stID][projID]*2-1;
	const int superdetid = globalconsts::detsuperid[stID][projID];
	
	// pair the hits by position:
	// if one hit on e.g. x and one hit on x' are closer than
	// the "spacing" defined for the planes, then the hits can be paired together.
	int idx1 = -1;
	int idx2 = -1;
	
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0)printf("nhits %d %d \n", nhits1, nhits2);
#endif
	for(int i = 0; i<nhits1; i++){
		idx1++;
		idx2 = -1;
		for(int j = 0; j<nhits2; j++){
			idx2++;
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef && threadIdx==0)printf("i %d j %d pos %1.4f %1.4f\n", i, j, );
#endif
			if( abs(hitcoll1.pos(idx1) - hitcoll2.pos(idx2)) > globalconsts::spacingplane[superdetid] ){
				continue;
			}
			
			for(bin = bin0; bin<bin0+Nbins; bin++){
				if( globalconsts::WCHitsBins_X[stID-2][0][bin] <= hitcoll1.chan(idx1) && 
				    hitcoll1.chan(idx1) <= globalconsts::WCHitsBins_X[stID-2][1][bin]){
#ifdef DEBUG
					printf("bin %d low %d high %d hit 1 elem %d hit 2 elem %d global bin %d \n", bin, globalconsts::WCHitsBins_X[stID-2][0][bin-bin0], globalconsts::WCHitsBins_X[stID-2][1][bin-bin0], ic[index].AllHits[ i ].elementID, ic[index].AllHits[ idx2 ].elementID, bin+npairs[bin]*Nbins);
#endif
					if(npairs[bin-bin0]<=MaxHits)hitpairs[bin-bin0+npairs[bin-bin0]*Nbins] = thrust::make_pair(idx1, idx2);
					npairs[bin-bin0]++;
				}
			}
			hitflag1[idx1] = 1;
			hitflag2[idx2] = 1;
		}
	}
	// here the hits that cannot be paired to another hit are paired to "nothing"
	// (but they still have to be paired to be used in the trackletteing)
	for(int i = 0; i<nhits1; i++){
		if(hitflag1[i]<1){
			for(bin = bin0; bin<bin0+Nbins; bin++){
#ifdef DEBUG
				printf("bin %d low %d high %d hit elem %d global bin %d \n", bin, globalconsts::WCHitsBins_X[stID-2][0][bin], globalconsts::WCHitsBins_X[stID-2][1][bin], ic[index].AllHits[ i ].elementID, bin+npairs[bin]*Nbins);
#endif
				if( globalconsts::WCHitsBins_X[stID-2][0][bin] <= hitcoll1.chan(i) && 
				    hitcoll1.chan(i) <= globalconsts::WCHitsBins_X[stID-2][1][bin]){
#ifdef DEBUG
					printf("bin %d low %d high %d hit elem %d global bin %d \n", bin, globalconsts::WCHitsBins_X[stID-2][0][bin-bin0], globalconsts::WCHitsBins_X[stID-2][1][bin-bin0], ic[index].AllHits[ i ].elementID, bin+npairs[bin]*Nbins);
#endif
					if(npairs[bin-bin0]<=MaxHits)hitpairs[bin-bin0+npairs[bin-bin0]*Nbins] = thrust::make_pair(i, -1);
					npairs[bin-bin0]++;
				}
			}
		}
	}
	for(int i = 0; i<nhits2; i++){
		if(hitflag2[i]<1){
			for(bin = bin0; bin<bin0+Nbins; bin++){
#ifdef DEBUG
				printf("bin %d low %d high %d hit elem %d global bin %d \n", bin, globalconsts::WCHitsBins_X[stID-2][0][bin], globalconsts::WCHitsBins_X[stID-2][1][bin], ic[index].AllHits[ i ].elementID, bin+npairs[bin]*Nbins);
#endif
				if( globalconsts::WCHitsBins_X[stID-2][0][bin] <= hitcoll2.chan(i) && 
				    hitcoll2.chan(i) <= globalconsts::WCHitsBins_X[stID-2][1][bin]){
#ifdef DEBUG
					printf("bin %d low %d high %d hit elem %d global bin %d \n", bin, globalconsts::WCHitsBins_X[stID-2][0][bin-bin0], globalconsts::WCHitsBins_X[stID-2][1][bin-bin0], ic[index].AllHits[ i ].elementID, bin+npairs[bin]*Nbins);
#endif
					if(npairs[bin-bin0]<=MaxHits)hitpairs[bin-bin0+npairs[bin-bin0]*Nbins] = thrust::make_pair(-1, i);
					npairs[bin-bin0]++;
				}
			}
		 }
	}
}



__device__ int make_hitpairs_in_station(const gHits hitcoll1, const int nhits1, const gHits hitcoll2, const int nhits2, thrust::pair<int, int>* hitpairs, short* hitidx1, short* hitidx2, short* hitflag1, short* hitflag2, const int stID, const int projID, const gPlane* planes, const float xmin, const float xmax){
	// I think we assume that by default we want to know where we are
	int npairs = 0;
	
	//declaring arrays for the hit lists
	for(int i = 0; i<100; i++){
		hitidx1[i] = hitidx2[i] = 0;
		hitflag1[i] = hitflag2[i] = 0;
	}
	
	//building the lists of hits for each detector plane
	const int detid1 = globalconsts::detsuperid[stID][projID]*2;
	const int detid2 = globalconsts::detsuperid[stID][projID]*2-1;
	const int superdetid = globalconsts::detsuperid[stID][projID];
	float p1x, p2x;
	int hitctr1 = 0, hitctr2 = 0;
	for(int i = 0; i<nhits1; i++){
		if(planes->deltapx[ detid1 ]>0){
			p1x = x_bep(detid1, hitcoll1.chan(i), planes);
			p2x = x_tep(detid1, hitcoll1.chan(i), planes);
		}else{
			p1x = x_tep(detid1, hitcoll1.chan(i), planes);
			p2x = x_bep(detid1, hitcoll1.chan(i), planes);
		}
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("detid %d xmin %1.6f xmax %1.6f deltapx %1.6f p1x %1.6f p2x %1.6f \n", detid1, xmin, xmax, planes->deltapx[ detid1 ], p1x, p2x);
#endif
		if( (p1x <= xmax) && (p2x >= xmin) ){ 
			hitidx1[hitctr1] = i;
			hitctr1++;
		}
	}
	for(int i = 0; i<nhits2; i++){
		if(planes->deltapx[ detid2 ]>0){
			p1x = x_bep(detid2, hitcoll2.chan(i), planes);
			p2x = x_tep(detid2, hitcoll2.chan(i), planes);
		}else{
			p1x = x_tep(detid2, hitcoll2.chan(i), planes);
			p2x = x_bep(detid2, hitcoll2.chan(i), planes);
		}
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("detid %d xmin %1.6f xmax %1.6f deltapx %1.6f p1x %1.6f p2x %1.6f \n", detid2, xmin, xmax, planes->deltapx[ detid2 ], p1x, p2x);
#endif
		if( (p1x <= xmax) && (p2x >= xmin) ){ 
			hitidx2[hitctr2] = i;
			hitctr2++;
		}
	}
	
	// pair the hits by position:
	// if one hit on e.g. x and one hit on x' are closer than
	// the "spacing" defined for the planes, then the hits can be paired together.
	int idx1 = -1;
	int idx2 = -1;
	for(int i = 0; i<hitctr1; i++){
		idx1++;
		idx2 = -1;
		for(int j = 0; j<hitctr2; j++){
			idx2++;
			if( abs( hitcoll1.pos(hitidx1[idx1]) - hitcoll2.pos(hitidx2[idx2]) ) > globalconsts::spacingplane[superdetid] ){
				continue;
			}
			
			hitpairs[npairs] = thrust::make_pair(hitidx1[idx1], hitidx2[idx2]);
			npairs++;
			hitflag1[idx1] = 1;
			hitflag2[idx2] = 1;
		}
	}
	// here the hits that cannot be paired to another hit are paired to "nothing"
	// (but they still have to be paired to be used in the trackletteing)
	for(int i = 0; i<hitctr1; i++){
		if(hitflag1[i]<1){
			hitpairs[npairs] = thrust::make_pair(hitidx1[i], -1);
			npairs++;
		}
	}
	for(int i = 0; i<hitctr2; i++){
		if(hitflag2[i]<1){
			hitpairs[npairs] = thrust::make_pair(-1, hitidx2[i]);
			npairs++;
		}
	}
	   	   
	return npairs;
}


// ----------------------------------------------------------------- //
// functions for selection of station 1 hits for back partial tracks //
// ----------------------------------------------------------------- // 

__device__ void SagittaRatioInStation1(const float x0, const float tx, const float y0, const float ty, int detid_lasthit, float* pos_exp, float* window, const float* z_, const float* costheta_, const float* sintheta_)
{
	float z_st3 = z_[detid_lasthit];
	float x_st3 = x0+tx*z_st3;
	float y_st3 = y0+ty*z_st3;

	//printf("det id %d z %1.3f x %1.3f y %1.3f \n", detid_lasthit, z_st3, x_st3, y_st3);
	
	float z_st1;
	float z_st2, x_st2, y_st2;
		
	short detid, detid_2, idx;
	float pos_st3, pos_st2;
	
	float s2_target, s2_dump;
	
	float pos_exp_target, pos_exp_dump;
	float win_target, win_dump;
	float p_min, p_max;
	
	for(int i = 0; i<3; i++){
		detid = 2*i+2;
		idx = globalconsts::planetype[detid/2-1];
		detid_2 = globalconsts::detsuperid[2][idx]*2;
		
		pos_st3 = x_st3*costheta_[detid_2] + y_st3*sintheta_[detid_2];

#ifdef DEBUG
		printf(" i %d idx %d  pos %1.3f \n", i, idx, pos_st3);
#endif	
		
		z_st1 = z_[detid];
		z_st2 = z_[detid_2];
		
		x_st2 = x0+tx*z_st2;
		y_st2 = y0+ty*z_st2;
		
		pos_st2 = x_st2*costheta_[detid_2] + y_st2*sintheta_[detid_2];
		
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("det id %d z %1.3f s_det id %d z %1.3f x %1.3f y %1.3f \n", detid, z_st1, detid_2, z_st2, x_st2, y_st2);
#endif	
        	s2_target = pos_st2 - pos_st3*(z_st2 - globalconsts::Z_TARGET)/(z_st3 - globalconsts::Z_TARGET);
        	s2_dump   = pos_st2 - pos_st3*(z_st2 - globalconsts::Z_DUMP)/(z_st3 - globalconsts::Z_DUMP);

		pos_exp_target = globalconsts::SAGITTA_TARGET_CENTER*s2_target + pos_st3*(z_st1 - globalconsts::Z_TARGET)/(z_st3 - globalconsts::Z_TARGET);
		pos_exp_dump   = globalconsts::SAGITTA_DUMP_CENTER*s2_dump + pos_st3*(z_st1 - globalconsts::Z_DUMP)/(z_st3 - globalconsts::Z_DUMP);
		win_target = fabs(s2_target*globalconsts::SAGITTA_TARGET_WIDTH);
		win_dump   = fabs(s2_dump*globalconsts::SAGITTA_DUMP_WIDTH);
		
		p_min = min(pos_exp_target - win_target, pos_exp_dump - win_dump);
		p_max = max(pos_exp_target + win_target, pos_exp_dump + win_dump);
		
		pos_exp[idx] = 0.5*(p_max + p_min);
		window[idx]  = 0.5*(p_max - p_min);
		//printf("idx %d pos_exp %1.3f window %1.3f \n", idx, pos_exp[idx], window[idx]);
	}
}



// --------------------------------------------- //
// functions to calculate x0 and tx in station 1 //
// and function to calculate inverse momentum    //
// --------------------------------------------- //

__device__ void calculate_x0_tx_st1(const float x0, const float tx, const float invP, const short charge, float &x0_st1, float &tx_st1)
{	
	tx_st1 = tx + globalconsts::PT_KICK_KMAG * invP * charge;
	x0_st1 = tx*globalconsts::Z_KMAG_BEND + x0 - tx_st1 * globalconsts::Z_KMAG_BEND;
}

__device__ void calculate_x0_tx_st1_with_errors(const float x0, const float tx, const float invP, const short charge, const float err_x0, const float err_tx, const float err_invP, float &x0_st1, float &tx_st1, float &err_x0_st1, float &err_tx_st1)
{	
	tx_st1 = tx + globalconsts::PT_KICK_KMAG * invP * charge;
	x0_st1 = tx*globalconsts::Z_KMAG_BEND + x0 - tx_st1 * globalconsts::Z_KMAG_BEND;
	
	err_tx_st1 = err_tx + fabs(err_invP*globalconsts::PT_KICK_KMAG);
	err_x0_st1 = err_x0 + fabs(err_invP*globalconsts::PT_KICK_KMAG)*globalconsts::Z_KMAG_BEND;
}

__device__ void calculate_x0_tx_st1(const gTracklet tkl, float &x0, float &tx)
{	
	tx = tkl.tx() + globalconsts::PT_KICK_KMAG * tkl.invP() * tkl.charge();
	x0 = tkl.tx()*globalconsts::Z_KMAG_BEND + tkl.x0() - tx * globalconsts::Z_KMAG_BEND;
}

__device__ void calculate_x0_tx_st1_with_errors(const gTracklet tkl, float &x0, float &tx, float &err_x0, float &err_tx)
{	
	tx = tkl.tx() + globalconsts::PT_KICK_KMAG * tkl.invP() * tkl.charge();
	x0 = tkl.tx()*globalconsts::Z_KMAG_BEND + tkl.x0() - tx * globalconsts::Z_KMAG_BEND;
	
	err_tx = tkl.err_tx() + fabs(tkl.err_invP()*globalconsts::PT_KICK_KMAG);
	err_x0 = tkl.err_x0() + fabs(tkl.err_invP()*globalconsts::PT_KICK_KMAG)*globalconsts::Z_KMAG_BEND;
}

__device__ float calculate_invP(float tx, float tx_st1, const short charge)
{
	return (tx_st1 - tx)*charge / globalconsts::PT_KICK_KMAG;
}

//EF: NB: prolly not. The function below is not as robust as I thought it was... should not be used and probably decommissioned 
__device__ float calculate_invP_charge(float tx, float tx_st1, short& charge)
{
	float invP = (tx_st1 - tx) / globalconsts::PT_KICK_KMAG;
	if(invP<0){
		charge = -1;
		invP*= charge;
	}else{
		charge = +1;
	}
	return invP;
}

__device__ float calculate_invP_error(float err_tx, float err_tx_st1)
{
	return ( err_tx - err_tx )/ globalconsts::PT_KICK_KMAG;
} 

/**
 * This function should be as simple as possible, in order to reduce the
 * computation time.  Therefore the condition of the charge determination
 * uses only "x0" and "tx" (at St. 2+3).  The formula was obtained 
 * practically by the study in DocDB 9505.  This function is valid for both 
 * the parallel and anti-parallel FMag+KMag polarity combination.  But it 
 * is _not_ guaranteed to be valid when the FMag and/or KMag field strength
 * is changed largely.
 */
__device__ short calculate_charge(float tx, float x0)
{
#ifdef E1039
	return -0.0033f * copysign(1.0, globalconsts::FMAGSTR) * x0 < tx  ?  +1  :  -1;
#else
	return x0*globalconsts::KMAGSTR > tx ? 1 : -1;
#endif
}

#ifdef TEST_MOMENTUM
/*
__device__ float calculate_x_fmag(const float tx_st1_tgt, const float tx, const short charge)
{
	float invP = calculate_invP(tx, tx_st1_tgt, charge);
	float tx_tgt = tx_st1_tgt + globalconsts::PT_KICK_FMAG * invP * charge;
	return (tx_tgt*(globalconsts::Z_FMAG_BEND-globalconsts::Z_TARGET));
}
*/

__device__ float calculate_invp_tgt(const float tx_st1, const float tx_tgt,  const short charge)
{
	return (tx_tgt - tx_st1)*charge / globalconsts::PT_KICK_FMAG;
}
#endif


// ----------------------------------------------------------- //
// function to resolve the left right ambiguities in the track //
// ----------------------------------------------------------- //

__device__ void resolve_leftright_newhits(const float x0, const float tx, const float y0, const float ty, const float err_x0, const float err_tx, const float err_y0, const float err_ty, const short nhits, const short* hits_detid, const float* hits_pos, const float* hits_drift, short* hits_sign, const gPlane* planes, const float thr)
{
	short i, j;
	int indexmin = -1;
	float pull_min = 1.e6;
	float pull;
	float slope_local, inter_local;
	float slope_exp, inter_exp;
	float err_slope, err_inter;
	short detID_i, detID_j;
	short k, n;
	for(n = 0; n<nhits; n+=2){
		i = n;
		j = i+1;
		detID_i = hits_detid[i];
		detID_j = hits_detid[j];
		//check this

		if( abs(detID_i-detID_j)!=1 ){
		    n--;//step back by 1 to move by 1 hit instead of 2
		    continue;		    
		}
		
		//that should be correct
		slope_exp = planes->costheta[detID_i]*tx + planes->sintheta[detID_i]*ty;
		err_slope = fabs(planes->costheta[detID_i]*err_tx) + fabs(planes->sintheta[detID_i]*err_ty);
		//that should be correct too
		inter_exp = planes->costheta[detID_i]*x0 + planes->sintheta[detID_i]*y0;
		err_inter = fabs(planes->costheta[detID_i]*err_x0) + fabs(planes->sintheta[detID_i]*err_y0);
		//position(const float pos, const float drift, const short sign)
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("hits dets %d %d; exp slope %1.4f +- %1.4f inter %1.4f +- %1.4f \n", detID_i, hits_detid[j], slope_exp, err_slope, inter_exp, err_inter);
		if(blockIdx.x==debug::EvRef)printf("hit 1 positions %1.4f, %1.4f hit 2 positions %1.4f, %1.4f \n", 
						position(hits_pos[i], hits_drift[i], +1), position(hits_pos[i], hits_drift[i], -1), 
						position(hits_pos[j], hits_drift[j], +1), position(hits_pos[i], hits_drift[j], -1));
#endif
		
		if(hits_sign[i]*hits_sign[j]==0){
			indexmin = -1;
			pull_min = 1.e6;
			for(k = 0; k<4; k++){
				slope_local = ( position( hits_pos[i], hits_drift[i],  globalconsts::lrpossibility[k][0]) - position(hits_pos[j], hits_drift[j],  globalconsts::lrpossibility[k][1]) ) / ( planes->z[detID_i]-planes->z[detID_j] );
				inter_local = position(hits_pos[i], hits_drift[i], globalconsts::lrpossibility[k][0]) - slope_local*planes->z[detID_i];
				
				if(fabs(slope_local) > planes->slope_max[detID_i] || fabs(inter_local) > planes->inter_max[detID_i])continue;
				
				pull = sqrtf( (slope_exp-slope_local)*(slope_exp-slope_local)/err_slope/err_slope + (inter_exp-inter_local)*(inter_exp-inter_local)/err_inter/err_inter );
				
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("lr %d %d, slope %1.4f inter %1.4f\n", globalconsts::lrpossibility[k][0], globalconsts::lrpossibility[k][1], slope_local, inter_local);
				if(blockIdx.x==debug::EvRef)printf("pull %1.4f\n", pull);
#endif		
				if(pull<pull_min){
					indexmin = k;
					pull_min = pull;
				}
			}
			
			if(indexmin>0 && pull_min<thr){
				hits_sign[i] = globalconsts::lrpossibility[indexmin][0];
				hits_sign[j] = globalconsts::lrpossibility[indexmin][1];
				//isUpdated = true;
			}
		}
	//	++nresolved;
	}
}



//This function just does the left right resolution for the newest hits (as a way to naturally handle st2-3 and st1 segments)
__device__ void resolve_single_leftright_newhits(const float x0, const float tx, const float y0, const float ty, const short nhits, const short* hits_detid, const float* hits_pos, short* hits_sign, const gPlane* planes)
{
	float pos_exp;
	short detID;
	
	for(short n = 0; n<nhits; n++){
		// don't do anything for hits whichs already have a sign...
		if(hits_sign[n])continue;
		detID = hits_detid[n];
		if(detID<0)continue;
				
		pos_exp = (planes->z[detID]*tx+x0)*planes->costheta[detID]+(planes->z[detID]*ty+y0)*planes->sintheta[detID];
		hits_sign[n] = pos_exp>hits_pos[n]? +1 : -1;

#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("pos_exp %1.4f hit_pos %1.4f sign %d \n", pos_exp, hits_pos[n], hits_sign[n]);
#endif
	}
	
}



__device__ void resolve_leftright_xhits(const float x0, const float tx, const float err_x0, const float err_tx, const short nhits, const short* hits_detid, const float* hits_pos, const float* hits_drift, short* hits_sign, const float* zarray, const float thr)
//const float* slope_max, const float* inter_max, 
{
	short i, j;
	int indexmin = -1;
	float pull_min = 1.e6;
	float pull;
	float slope_local, inter_local;
	float slope_exp, inter_exp;
	float err_slope, err_inter;
	short detID_i, detID_j;

	for(short n = 0; n<nhits; n+=2){
		i = n;
		j = i+1;
		detID_i = hits_detid[i];
		detID_j = hits_detid[j];
		if( abs(detID_i-detID_j)!=1 ){
		    n--;//step back by 1 to move by 1 hit instead of 2
		    continue;		    
		}
				
		slope_exp = tx;
		err_slope = err_tx;
		
		inter_exp = x0;
		err_inter = err_x0;
		//position(const float pos, const float drift, const short sign)
#ifdef DEBUG
		printf("hits dets %d %d; exp slope %1.4f +- %1.4f inter %1.4f +- %1.4f \n", detID_i, hits_detid[j], slope_exp, err_slope, inter_exp, err_inter);
		printf("hit 1 positions %1.4f, %1.4f hit 2 positions %1.4f, %1.4f \n", 
			position(hits_pos[i], hits_drift[i], +1), position(hits_pos[i], hits_drift[i], -1), 
			position(hits_pos[j], hits_drift[j], +1), position(hits_pos[i], hits_drift[j], -1));
#endif
		
		if(hits_sign[i]*hits_sign[j]==0){
			indexmin = -1;
			pull_min = 1.e6;
			for(int k = 0; k<4; k++){
				slope_local = ( position( hits_pos[i], hits_drift[i],  globalconsts::lrpossibility[k][0]) - position(hits_pos[j], hits_drift[j],  globalconsts::lrpossibility[k][1]) ) / ( zarray[detID_i]-zarray[detID_j] );
				inter_local = position(hits_pos[i], hits_drift[i], globalconsts::lrpossibility[k][0]) - slope_local*zarray[detID_i];
				
				//if(fabs(slope_local) > slope_max[detID_i] || fabs(inter_local) > inter_max[detID_i])continue;
				
				pull = sqrtf( (slope_exp-slope_local)*(slope_exp-slope_local)/err_slope/err_slope + (inter_exp-inter_local)*(inter_exp-inter_local)/err_inter/err_inter );
				
#ifdef DEBUG
				printf("lr %d %d, slope %1.4f inter %1.4f\n", globalconsts::lrpossibility[k][0], globalconsts::lrpossibility[k][1], slope_local, inter_local);
				printf("pull %1.4f\n", pull);
#endif		
				if(pull<pull_min){
					indexmin = k;
					pull_min = pull;
				}
			}
			
			if(indexmin>0 && pull_min<thr){
				hits_sign[i] = globalconsts::lrpossibility[indexmin][0];
				hits_sign[j] = globalconsts::lrpossibility[indexmin][1];
				//isUpdated = true;
			}
		}
	//	++nresolved;
	}
}


__device__ void resolve_single_leftright_xhits(const float x0, const float tx, const short nhits, const short* hits_detid, const float* hits_pos, short* hits_sign, const float* z_array)
{
	float pos_exp;
	short detID;
	
	for(short n = 0; n<nhits; n++){
		// don't do anything for hits whichs already have a sign...
		if(hits_sign[n])continue;
		detID = hits_detid[n];
		
		pos_exp = z_array[detID]*tx+x0;
		hits_sign[n] = pos_exp>hits_pos[n]? +1 : -1;
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("pos_exp %1.4f hit_pos %1.4f sign %d \n", pos_exp, hits_pos[n], hits_sign[n]);
#endif
	}
	
}


__device__ void resolve_leftright(const gTracklet tkl, float* hitsign, const gPlane* planes, const float thr)
{
	//bool isUpdated = false;
	short nhits = tkl.nHits();
	//short nresolved = 0;
	short i, j;
	int indexmin = -1;
	float pull_min = 1.e6;
	float pull;
	float slope_local, inter_local;
	float slope_exp, inter_exp;
	float err_slope, err_inter;
	float x0, tx;// x0 and tx are different for global track station 1 hits 
	float err_x0, err_tx;// x0 and tx are different for global track station 1 hits 
	short detID_i, detID_j;

	for(short n = 0; n<nhits; n+=2){
		i = n;
		j = i+1;
		detID_i = (short)tkl.hits_detid(i);
		detID_j = (short)tkl.hits_detid(j);
		if( abs(detID_i-detID_j)!=1 ){
		    n--;//step back by 1 to move by 1 hit instead of 2
		    continue;		    
		}
		
		if(tkl.stationID()>=6 && detID_i<=12){
			calculate_x0_tx_st1_with_errors(tkl, x0, tx, err_x0, err_tx);
		}else{
			tx = tkl.tx();
			x0 = tkl.x0();
			err_x0 = tkl.err_x0();
			err_tx = tkl.err_tx();
		}
		
		slope_exp = planes->costheta[detID_i]*tx + planes->sintheta[detID_i]*tkl.ty();
		err_slope = fabs(planes->costheta[detID_i]*err_tx) + fabs(planes->sintheta[detID_i]*tkl.err_ty());
		
		inter_exp = planes->costheta[detID_i]*x0 + planes->sintheta[detID_i]*tkl.y0();
		err_inter = fabs(planes->costheta[detID_i]*err_x0) + fabs(planes->sintheta[detID_i]*tkl.err_y0());
		//position(const float pos, const float drift, const short sign)
#ifdef DEBUG
		printf("hits dets %d %d; exp slope %1.4f +- %1.4f inter %1.4f +- %1.4f \n", detID_i, tkl.hits_detid(j), slope_exp, err_slope, inter_exp, err_inter);
		printf("hit 1 positions %1.4f, %1.4f hit 2 positions %1.4f, %1.4f \n", 
			position(tkl.hits_pos(i), tkl.hits_drift(i), +1), position(tkl.hits_pos(i), tkl.hits_drift(i), -1), 
			position(tkl.hits_pos(j), tkl.hits_drift(j), +1), position(tkl.hits_pos(i), tkl.hits_drift(j), -1));
#endif

		if(tkl.hits_sign(i)*tkl.hits_sign(j)==0){
			indexmin = -1;
			pull_min = 1.e6;
			for(int k = 0; k<4; k++){
				slope_local = ( position( tkl.hits_pos(i), tkl.hits_drift(i),  globalconsts::lrpossibility[k][0]) - position(tkl.hits_pos(j), tkl.hits_drift(j),  globalconsts::lrpossibility[k][1]) ) / ( planes->z[detID_i]-planes->z[detID_j] );
				inter_local = position(tkl.hits_pos(i), tkl.hits_drift(i), globalconsts::lrpossibility[k][0]) - slope_local*planes->z[detID_i];
				
				if(fabs(slope_local) > planes->slope_max[detID_i] || fabs(inter_local) > planes->inter_max[detID_i])continue;
				
				pull = sqrtf( (slope_exp-slope_local)*(slope_exp-slope_local)/err_slope/err_slope + (inter_exp-inter_local)*(inter_exp-inter_local)/err_inter/err_inter );
				
#ifdef DEBUG
				printf("lr %d %d, slope %1.4f inter %1.4f\n", globalconsts::lrpossibility[k][0], globalconsts::lrpossibility[k][1], slope_local, inter_local);
				printf("pull %1.4f\n", pull);
#endif		
				if(pull<pull_min){
					indexmin = k;
					pull_min = pull;
				}
			}
			
			if(indexmin>0 && pull_min<thr){
				hitsign[i] = globalconsts::lrpossibility[indexmin][0];
				hitsign[j] = globalconsts::lrpossibility[indexmin][1];
				//isUpdated = true;
			}
		}
	//	++nresolved;
	}
}


__device__ void resolve_single_leftright(const gTracklet tkl, float* hitsign, const gPlane* planes)
{
	short nhits = tkl.nHits();
	float pos_exp;
	short detID;
	float x0, tx;// x0 and tx are different for global track station 1 hits 
	
	for(short n = 0; n<nhits; n++){
		// don't do anything for hits whichs already have a sign...
		if(tkl.hits_sign(n)!=0 || hitsign[n])continue;
		detID = tkl.hits_detid(n);
				
		if(tkl.stationID()>=6 && detID<=12){
			calculate_x0_tx_st1(tkl, x0, tx);
		}else{
			tx = tkl.tx();
			x0 = tkl.x0();
		}
		
		pos_exp = (planes->z[detID]*tx+x0)*planes->costheta[detID]+(planes->z[detID]*tkl.ty()+tkl.y0())*planes->sintheta[detID];
		hitsign[n] = pos_exp>tkl.hits_pos(n)? +1 : -1;
	}
	
}


// functions to match a tracklet to a hodoscope hit
__device__ bool match_tracklet_to_hodo(int stid, const int detid, const int nhits, const gHits hits, const float x0, const float y0, const float tx, const float ty, const float err_x0, const float err_y0, const float err_tx, const float err_ty, const gPlane* planes)
{
	int masked = 0;//false
	if(nhits==0)return masked;
	
	if(detid>40)stid = 3;
	// first, define the search region, and foremost, the planes on which we define this search region, which depends on the station ID we're looking at
	// define the region in which we are supposed to have hits:

	const float fudgefac = globalconsts::hodofudgefac[stid];
	const float cellwidth = planes->cellwidth[detid];
	
	const float xhodo = planes->z[detid]*tx+x0;
	const float yhodo = planes->z[detid]*ty+y0;
	
	float err_x = 3.f*(fabs(planes->z[detid]*err_tx)+err_x0);
	float err_y = 3.f*(fabs(planes->z[detid]*err_ty)+err_y0);
	
	float xmin, xmax, ymin, ymax;
	
	//we only consider hits in the hodoscopes planes corresponding to the station where the tracklet is reconstructed 
	//calculate the track position at the hodoscope plane z
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf(" evt %d stid %d detid %d, xhodo %1.4f yhodo %1.4f zhodo %1.4f posErrX %1.4f posErrY %1.4f errx %1.4f erry %1.4f fudge factor %1.2f \n", blockIdx.x, stid, detid, xhodo, yhodo, planes->z[detid], fabs(planes->z[detid]*err_tx)+err_x0, fabs(planes->z[detid]*err_ty)+err_y0, err_x, err_y, globalconsts::hodofudgefac[stid] );
#endif
		
	// loop on the hits and select hodoscope hits corresponding to the station
	for(int i = 0; i<nhits; i++){
		//calculate "xmin, xmax, ymin, ymax" in which the track is supposed to pass through; 
		//these are basicially defined as the spatial coverage of the hit hodoscope element (plus some fudge factor for x)
#ifdef MATCH_Y_HODO
		if(planes->costheta[detid]>0.99){
#endif
			xmin = hits.pos(i)-cellwidth*0.5f;
			xmax = hits.pos(i)+cellwidth*0.5f;
			
			ymin = planes->y1[detid];
			ymax = planes->y2[detid];
			
			xmin-=(xmax-xmin)*fudgefac;
			xmax+=(xmax-xmin)*fudgefac;
			
			ymin-=(ymax-ymin)*fudgefac;
			ymax+=(ymax-ymin)*fudgefac;
#ifdef MATCH_Y_HODO
		}else{

			xmin = planes->x1[detid];
			xmax = planes->x2[detid];
			
			ymin = hits.pos(i)-cellwidth*0.5f;
			ymax = hits.pos(i)+cellwidth*0.5f;
			
			xmin-=(xmax-xmin)*fudgefac;
			xmax+=(xmax-xmin)*fudgefac;
			
			ymin-=(ymax-ymin)*fudgefac;
			ymax+=(ymax-ymin)*fudgefac;
		}
#endif

#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf(" evt %d detid %d elid %1.0f pos %1.4f xmin %1.4f xmax %1.4f, ymin %1.4f ymax %1.4f\n", blockIdx.x, detid, hits.chan(i), hits.pos(i), xmin, xmax, ymin, ymax );
#endif
		err_x+= (xmax-xmin)*0.15f;

		xmin-= err_x;
		xmax+= err_x;
		ymin-= err_y;
		ymax+= err_y;
		
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf(" errx %1.4f xmin %1.4f xmax %1.4f, erry %1.4f ymin %1.4f ymax %1.4f \n", err_x, xmin, xmax, err_y, ymin, ymax);
#endif
		if(xmin <= xhodo && xhodo <= xmax && ymin <= yhodo && yhodo <= ymax ){
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("evt %d track x0 %1.4f det %d match in element %d \n", blockIdx.x, x0, detid, (short)hits.chan(i) );
#endif
			masked++;
			break;
		}
		
	}
	// 
	return masked>0;
}


__device__ bool match_track_XZ_to_hodo(int stid, const int detid, const int nhits, const gHits hits, const float x0, const float tx, const float err_x0, const float err_tx, const float* z_array)
{
	int masked = 0;//false
	if(nhits==0)return masked;
	
	if(detid>40)stid = 3;
	// first, define the search region, and foremost, the planes on which we define this search region, which depends on the station ID we're looking at
	// define the region in which we are supposed to have hits:
	
	const float xhodo = z_array[detid]*tx+x0;
	const float fudgefac = globalconsts::hodofudgefac[stid];
	const float cellwidth = globalconsts::hodoXcellwidth[stid];
	
	float err_x = 3.f*(fabs(z_array[detid]*err_tx)+err_x0);
	float xmin, xmax;
	
	//we only consider hits in the hodoscopes planes corresponding to the station where the tracklet is reconstructed 
	//calculate the track position at the hodoscope plane z
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf(" evt %d stid %d detid %d, xhodo %1.4f yhodo %1.4f zhodo %1.4f posErrX %1.4f posErrY %1.4f errx %1.4f erry %1.4f fudge factor %1.2f \n", blockIdx.x, stid, detid, xhodo, yhodo, z_array[detid], fabs(z_array[detid]*err_tx)+err_x0, fabs(z_array[detid]*err_ty)+err_y0, err_x, err_y, globalconsts::hodofudgefac[stid] );
#endif
		
	// loop on the hits and select hodoscope hits corresponding to the station
	for(int i = 0; i<nhits; i++){
		//calculate "xmin, xmax, ymin, ymax" in which the track is supposed to pass through; 
		//these are basicially defined as the spatial coverage of the hit hodoscope element (plus some fudge factor for x)
		xmin = hits.pos(i)-cellwidth*0.5f;
		xmax = hits.pos(i)+cellwidth*0.5f;
			
		xmin-=(xmax-xmin)*fudgefac;
		xmax+=(xmax-xmin)*fudgefac;

#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf(" evt %d detid %d elid %1.0f pos %1.4f xmin %1.4f xmax %1.4f, ymin %1.4f ymax %1.4f\n", blockIdx.x, detid, hits.chan(i), hits.pos(i), xmin, xmax, ymin, ymax );
#endif
		err_x+= (xmax-xmin)*0.15f;

		xmin-= err_x;
		xmax+= err_x;
		
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf(" errx %1.4f xmin %1.4f xmax %1.4f \n", err_x, xmin, xmax);
#endif
		if(xmin <= xhodo && xhodo <= xmax ){
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("evt %d track x0 %1.4f det %d match in element %d \n", blockIdx.x, x0, detid, (short)hits.chan(i) );
#endif
			masked++;
			break;
		}
		
	}
	// 
	return masked>0;
}


//Functions to calculate residuals and chi2

__device__ float residual(const short detid, const short elid, const float drift, const short sign, const gPlane* planes, const float x0, const float y0, const float tx, const float ty)
{
	float p1x = x_bep(detid, elid, planes);
	float p1y = y_bep(detid, elid, planes);
	float p1z = z_bep(detid, elid, planes);
	float deltapx = planes->deltapx[detid];
	float deltapy = planes->deltapy[detid];
	float deltapz = planes->deltapz[detid];
		
	float den2 = deltapy*deltapy*(1+tx*tx) + deltapx*deltapx*(1+ty*ty) - 2*( ty*deltapx*deltapz + ty*deltapy*deltapz + tx*ty*deltapx*deltapy);
	float dca = ( (ty*deltapz-deltapy)*(p1x-x0) + (deltapx-tx*deltapz)*(p1y-y0) + p1z*(tx*deltapy-ty*deltapx) ) / sqrtf(den2);

#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf("x0 %1.4f y0 %1.4f tx %1.4f ty %1.4f p1x %1.4f p1y %1.4f p1z %1.4f deltapx %1.4f deltapy %1.4f deltapz %1.4f den2 %1.4f dca %1.4f drift %1.4f sign %d \n", x0, y0, tx, ty, p1x, p1y, p1z, deltapx, deltapy, deltapz, den2, dca, drift, sign);
#endif	
	return drift*sign - dca;
}

__device__ float residual(float const p1x, float const p1y, float const p1z,
			float const deltapx, float const deltapy, float const deltapz,
			const float drift, const short sign, 
			const float x0, const float y0, const float tx, const float ty)
{
	float den2 = deltapy*deltapy*(1+tx*tx) + deltapx*deltapx*(1+ty*ty) - 2*( ty*deltapx*deltapz + ty*deltapy*deltapz + tx*ty*deltapx*deltapy);
	float dca = ( (ty*deltapz-deltapy)*(p1x-x0) + (deltapx-tx*deltapz)*(p1y-y0) + p1z*(tx*deltapy-ty*deltapx) ) / sqrtf(den2);

#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf("x0 %1.4f y0 %1.4f tx %1.4f ty %1.4f p1x %1.4f p1y %1.4f p1z %1.4f deltapx %1.4f deltapy %1.4f deltapz %1.4f den2 %1.4f dca %1.4f drift %1.4f sign %d \n", x0, y0, tx, ty, p1x, p1y, p1z, deltapx, deltapy, deltapz, den2, dca, drift, sign);
#endif	
	return drift*sign - dca;
}


__device__ float chi2_track(size_t const n_points, float* residuals,
			float* const driftdist, short* const sign, float* const resolutions,
			float* const p1x, float* const p1y, float* const p1z,
			float* const deltapx, float* const deltapy, float* const deltapz,
			const float x0, const float y0, const float tx, const float ty)
{
	float dca;
	float chi2 = 0;
	float den2;
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf(" x0 %1.6f tx %1.6f y0 %1.6f ty %1.6f \n", x0, tx, y0, ty);
#endif
	for( size_t i=0; i<n_points; i++ ){
		den2 = deltapy[i]*deltapy[i]*(1+tx*tx) + deltapx[i]*deltapx[i]*(1+ty*ty) - 2*( ty*deltapx[i]*deltapz[i] + ty*deltapy[i]*deltapz[i] + tx*ty*deltapx[i]*deltapy[i]);
		dca = ( (ty*deltapz[i]-deltapy[i])*(p1x[i]-x0) + (deltapx[i]-tx*deltapz[i])*(p1y[i]-y0) + p1z[i]*(tx*deltapy[i]-ty*deltapx[i]) ) / sqrtf(den2);
		residuals[i] = driftdist[i]*sign[i] - dca;
		chi2+= residuals[i] * residuals[i] / resolutions[i] / resolutions[i];
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf(" thread %d p1x %1.6f p1y %1.6f p1z %1.6f dpx %1.6f dpy %1.6f dpz %1.6f dca %1.6f drift dist %1.6f * sign %d resid %1.6f resol %1.6f chi2 %1.6f \n", threadIdx.x, p1x[i], p1y[i], p1z[i], deltapx[i], deltapy[i], deltapz[i], dca, driftdist[i], sign[i], residuals[i], resolutions[i], chi2);
#endif
	}
	return chi2;
}


// --------------------------------------------- //
//
// !!! track cleaning kernels !!!
// 
// --------------------------------------------- //


__global__ void gKernel_TrackOutlierHitRemoval(
	gEventTrackCollection* tklcoll,
	const short track_stid_ref,
	const gPlane* planes,
	const bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x])return;
	
	//tracks
	unsigned int Ntracks;
	const gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);
	
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	int i;
	short ihit, nhitsi;
	//float res_curr, res_max, res_max2;
	//short hit_rm_idx; 
	short hit_rm_neighbor;
	float cut;
	
	float x0, errx0, y0, erry0, tx, errtx, ty, errty, x0_st1, tx_st1;
	float y, err_y;
	short charge;
	float invP, errinvP;
#ifdef TEST_MOMENTUM
	float tx_tgt, invp_tgt;
#endif
	const short ndof = track_stid_ref>5? 5 : 4;
	
	bool isupdated;
	//short signflip[18];
	//for(short k = 0; k<18; k++){
	//	signflip[k] = 0;
	//}
	float resid_[18];	
	
	short stid;
	short detid[18];
	short chan[18];
	short sign[18];
	float pos[18];
	float drift[18];
	
	float X_[4];
	float errX_[4];
	float Z_23[4];
	float X_st1_[3];
	float errX_st1_[3];
	float Z_1[3];
	float Y_[12];
	float errY_[12];
	float Z_[12];
	
	float A_[4];
	float Ainv_[4];
	float B_[2]; 
	float Par[2];
	float ParErr[2];
	float chi2;
	
	float tdc[18];
	short nhits_st23, nhits_st1, nhits_x, nhits_uv, nhits, nhits_x_st23, nhits_x_st1;
	
	for(i = 0; i<Ntracks; i++){
		if(Tracks.stationID(i)<track_stid_ref)continue;
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
#endif
		nhitsi = Tracks.nHits(i);
		nhits_st23 = 0;
		nhits_st1 = 0;
		
		x0 = Tracks.x0(i);
		y0 = Tracks.y0(i);
		tx = Tracks.tx(i);
		ty = Tracks.ty(i);

		//res_max = -1.;
		//res_max2 = -1.;
		//hit_rm_idx = -1;
		hit_rm_neighbor = -1;
		
		if(track_stid_ref>5){
			calculate_x0_tx_st1(x0, tx, Tracks.invP(i), Tracks.charge(i), x0_st1, tx_st1);
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("st1 var calculation: x0 %1.4f tx %1.4f x0_st1 %1.4f tx_st1 %1.4f  \n", x0, tx, x0_st1, tx_st1);
#endif
		}
			
		isupdated = false;
		//first fill in the arrays... that makes sense...
		for(ihit = 0; ihit<nhitsi; ihit++){
			detid[ihit] = Tracks.hits_detid(i, ihit);
			chan[ihit] = Tracks.hits_chan(i, ihit);
			sign[ihit] = Tracks.hits_sign(i, ihit);
			pos[ihit] = Tracks.hits_pos(i, ihit);
			drift[ihit] = Tracks.hits_drift(i, ihit);
			tdc[ihit] = Tracks.hits_tdc(i, ihit);
			resid_[ihit] = Tracks.hits_residual(i, ihit);

			if(detid[ihit]>12){
				nhits_st23++;
			//	resid_[ihit] = residual(detid[ihit], chan[ihit], drift[ihit], sign[ihit], planes, x0, y0, tx, ty);
			}else{
				nhits_st1++;
			//	resid_[ihit] = residual(detid[ihit], chan[ihit], drift[ihit], sign[ihit], planes, x0_st1, y0, tx_st1, ty);
			}
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("hit %d detid %d drift %1.4f resid %1.4f sign %d \n", ihit, detid[ihit], drift[ihit], resid_[ihit], sign[ihit]);
#endif
		}
		//then reloop *once* to apply the logic....
		for(ihit = 0; ihit<nhitsi; ihit++){
			
			hit_rm_neighbor = detid[ihit]%2==0? ihit+1: ihit-1;
			if( abs( detid[hit_rm_neighbor] - detid[ihit] )>1 ||  detid[hit_rm_neighbor]<0 ||  hit_rm_neighbor>=nhitsi )hit_rm_neighbor = -1;
			
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf(" hit %d det %d neighbor %d det %d \n", ihit, detid[ihit], hit_rm_neighbor, detid[hit_rm_neighbor]);
#endif			
			stid = (detid[ihit]-1)/6-1;
			if(stid<0 && detid[ihit]>=0)stid = 0;
			
			cut = sign[ihit]==0 ? selection::rejectwin[stid]+fabs(drift[ihit]) : selection::rejectwin[stid];

			if( fabs(resid_[ihit]) > cut ){
				//if the neighbor hit id good, is larger than current hit and has a larger residual, skip because we will do it next 
				if(hit_rm_neighbor>ihit && fabs(resid_[hit_rm_neighbor])>fabs(resid_[ihit])){
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("Check next hit instead\n");
#endif
					continue;
				}

#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("cut %1.4f hit %d detid %d resid %1.4f drift %1.4f sign %d \n", cut, ihit, detid[ihit], resid_[ihit], drift[ihit], sign[ihit]);
#endif
						
				isupdated = true;
				
				if(fabs(resid_[ihit]-2*drift[ihit]*sign[ihit]) < cut){
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("hit %d det %d sign changed \n", ihit, detid[ihit]);
#endif
					sign[ihit] = -sign[ihit];
					if(hit_rm_neighbor>=0)sign[hit_rm_neighbor] = 0;
				}else{
					detid[ihit] = -1;
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("removing hit %d detid %d resid %1.4f, %1.4f > cut %1.4f \n", ihit, detid[ihit], fabs(resid_[ihit]), fabs(resid_[ihit]-2*drift[ihit]*sign[ihit]), cut);
#endif
					if(hit_rm_neighbor<0){
						// if we have two hits in the same chamber that are not good, then the track cannot be kept 
						// we don't need to do the follow up either
#ifdef DEBUG
						if(blockIdx.x==debug::EvRef)
						printf("evt %d track removed \n", blockIdx.x, threadIdx.x);
#endif
						tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, Tracks.stationID(i)-2);
						isupdated = false;
						break;
					}
					if(detid[ihit]>12){
						nhits_st23--;
					}else{
						nhits_st1--;
					}
				}
				
				
			}
		}
		if(isupdated){
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("updating: evt %d thread %d x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f stid %1.0f \n", blockIdx.x, threadIdx.x, x0, tx, y0, ty, Tracks.stationID(i));
#endif			
			nhits_x = 0;
			nhits_uv = 0;
			nhits_x_st23 = 0;
			nhits_x_st1 = 0;
			
			resolve_single_leftright_newhits(x0, tx, y0, ty, nhits_st23, detid, pos, sign, planes);
			resolve_single_leftright_newhits(x0_st1, tx_st1, y0, ty, nhits_st1, detid+nhits_st23, pos+nhits_st23, sign+nhits_st23, planes);

			for(ihit = 0; ihit<nhitsi; ihit++){
				if(detid[ihit]>=0 && (detid[ihit]-1)%6==2 || (detid[ihit]-1)%6==3 ){
					if( detid[ihit]>12 ){
						Z_23[nhits_x_st23] = planes->z[detid[ihit]];
						X_[nhits_x_st23] = pos[ihit];
						errX_[nhits_x_st23] = planes->resolution[detid[ihit]];
						nhits_x_st23++;
					}else{
						Z_1[nhits_x_st1] = planes->z[detid[ihit]];
						X_st1_[nhits_x_st1] = pos[ihit];
						errX_st1_[nhits_x_st1] = planes->resolution[detid[ihit]];
						nhits_x_st1++;
					}
				}
			}
			
			fit_2D_track(nhits_x_st23, X_, Z_23, errX_, A_, Ainv_, B_, Par, ParErr, chi2);
			x0 = Par[0];
			errx0 = ParErr[0];
			tx = Par[1];
			errtx = ParErr[1];
			
			if(Tracks.stationID(i)>5){
				Z_1[nhits_x_st1] = globalconsts::Z_KMAG_BEND;
				X_st1_[nhits_x_st1] = x_trk(x0, tx, Z_1[nhits_x_st1]);
				errX_st1_[nhits_x_st1] = err_x_trk(errx0, errtx, Z_1[nhits_x_st1]);
				
				fit_2D_track(nhits_x_st1+1, X_st1_, Z_1, errX_st1_, A_, Ainv_, B_, Par, ParErr, chi2);
				x0_st1 = Par[0];
				tx_st1 = Par[1];
				
				charge = calculate_charge(tx, x0);
				invP = calculate_invP(tx, tx_st1, charge);
				//invP = calculate_invP_charge(tx, tx_st1, charge);
				errinvP = calculate_invP_error(errtx, ParErr[1]);
				
#ifdef TEST_MOMENTUM
				tx_tgt = x_trk(x0_st1, tx_st1, globalconsts::Z_FMAG_BEND)/(globalconsts::Z_FMAG_BEND-globalconsts::Z_TARGET);
				invp_tgt = calculate_invp_tgt(tx_st1, tx_tgt, charge);
#endif
#ifdef DEBUG
				if(invP>0.25)printf(" evt %d thread %d x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f x0_st1 %1.4f tx_st1 %1.4f charge %d \n", blockIdx.x, threadIdx.x, x0, tx, y0, ty, x0_st1, tx_st1, charge);
#endif
			}
			
			for(ihit = 0; ihit<nhitsi; ihit++){
				if(detid[ihit]>=0 && (detid[ihit]-1)%6!=2 && (detid[ihit]-1)%6!=3 ){
					if( detid[ihit]>12 ){
						calculate_y_uvhit(detid[ihit], chan[ihit], drift[ihit], sign[ihit], x0, tx, planes, y, err_y);
					}else{
						calculate_y_uvhit(detid[ihit], chan[ihit], drift[ihit], sign[ihit], x0_st1, tx_st1, planes, y, err_y);
					}
					Y_[nhits_uv] = y;
					errY_[nhits_uv] = err_y;
					Z_[nhits_uv] = planes->z[detid[ihit]];
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("hit %d detid %d z %1.4f y %1.4f erry %1.4f\n", nhits_uv, detid[ihit], Z_[nhits_uv], Y_[nhits_uv], errY_[nhits_uv] );
#endif
					nhits_uv++;
				}
			}
			
			fit_2D_track(nhits_uv, Y_, Z_, errY_, A_, Ainv_, B_, Par, ParErr, chi2);
			y0 = Par[0];
			erry0 = ParErr[0];
			ty = Par[1];
			errty = ParErr[1];
			//update the track parameters:
			
			tklcoll->setTx(tkl_coll_offset+array_thread_offset, i, tx);
			tklcoll->setX0(tkl_coll_offset+array_thread_offset, i, x0);
			tklcoll->setErrTx(tkl_coll_offset+array_thread_offset, i, errtx);
			tklcoll->setErrX0(tkl_coll_offset+array_thread_offset, i, errx0);
			tklcoll->setTy(tkl_coll_offset+array_thread_offset, i, ty);
			tklcoll->setY0(tkl_coll_offset+array_thread_offset, i, y0);
			tklcoll->setErrTy(tkl_coll_offset+array_thread_offset, i, errty);
			tklcoll->setErrY0(tkl_coll_offset+array_thread_offset, i, erry0);
			tklcoll->setinvP(tkl_coll_offset+array_thread_offset, i, invP);
			tklcoll->setErrinvP(tkl_coll_offset+array_thread_offset, i, errinvP);
#ifdef TEST_MOMENTUM
			tklcoll->setInvPTarget(tkl_coll_offset+array_thread_offset, i, invp_tgt);
#endif			
			chi2 = 0;
			nhits = 0;
			for(ihit = 0; ihit<nhitsi; ihit++){
				if(detid[ihit]>=0){
					resid_[nhits] = residual(detid[ihit], chan[ihit], drift[ihit], sign[ihit], planes, x0, y0, tx, ty);
					chi2+= resid_[nhits]*resid_[nhits]/planes->resolution[detid[ihit]]/planes->resolution[detid[ihit]];
					tklcoll->setHitDetID(tkl_coll_offset+array_thread_offset, i, nhits, detid[ihit]);
					tklcoll->setHitChan(tkl_coll_offset+array_thread_offset, i, nhits, chan[ihit]);
					tklcoll->setHitPos(tkl_coll_offset+array_thread_offset, i, nhits, pos[ihit]);
					tklcoll->setHitDrift(tkl_coll_offset+array_thread_offset, i, nhits, drift[ihit]);
					tklcoll->setHitSign(tkl_coll_offset+array_thread_offset, i, nhits, sign[ihit]);
					tklcoll->setHitTDC(tkl_coll_offset+array_thread_offset, i, nhits, tdc[ihit]);
					tklcoll->setHitResidual(tkl_coll_offset+array_thread_offset, i, nhits, resid_[ihit]);
					nhits++;
				}
			}
			tklcoll->setnHits(tkl_coll_offset+array_thread_offset, i, nhits);
			tklcoll->setChisq(tkl_coll_offset+array_thread_offset, i, chi2);
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("evt %d thread %d x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f stid %1.0f \n", blockIdx.x, threadIdx.x, x0, tx, y0, ty, Tracks.stationID(i));
#endif

		}
	}
}

// --------------------------------------------- //
//
// function to clean the tracks after XZ-YZ processing
// 
// --------------------------------------------- //

__global__ void gKernel_BackTrackCleaning(
	gEventTrackCollection* tklcoll,
	const bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x])return;
	
	//tracks
	unsigned int Ntracks;
	const gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);
	
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	int i, j; 
	short ihit, jhit;
	short nhitsi, nhitsj;
	
	int nhits_common = 0;	
	
	for(i = 0; i<Ntracks; i++){
		//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
		if(Tracks.stationID(i)<5)continue;
		nhitsi = Tracks.nHits(i);
		nhits_common = 0;
		
		for(j = i+1; j<Ntracks; j++){
		//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
			if(Tracks.stationID(j)<5)continue;
			nhitsj = Tracks.nHits(j);
			
			for(ihit = 0; ihit<nhitsi; ihit++){
				for(jhit = 0; jhit<nhitsj; jhit++){
					if( Tracks.hits_detid(j, jhit)==Tracks.hits_detid(i, ihit) && Tracks.hits_chan(j, jhit)==Tracks.hits_chan(i, ihit) ){
						nhits_common++;
						break;
					}
				}
			}
			
			if(Tracks.chisq(i)/(nhitsi-4) < Tracks.chisq(j)/(nhitsj-4)){
				if(nhits_common>nhitsi*0.3333f ){
					tklcoll->setStationID(tkl_coll_offset+array_thread_offset, j, 3);
				}
			
			}else{
				if(nhits_common>nhitsj*0.3333f ){
					tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 3);
				}
			
			} 

		}	

	}
}


__global__ void gKernel_BackTrackCleaning_crossthreads(
	gEventTrackCollection* tklcoll,
	const bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x])return;
	
	//tracks
	unsigned int Ntracks, Ntracks_trd;
	gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);
	unsigned int array_thread_offset[THREADS_PER_BLOCK];// = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	int i_trd;
	for(i_trd = 0; i_trd<THREADS_PER_BLOCK; i_trd++){
		//Tracks[i_trd] = tklcoll->tracks(blockIdx.x, i_trd, Ntracks[i_trd]);
		array_thread_offset[i_trd] = i_trd*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	}
	
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	
	int i, j; 
	short ihit, jhit;
	short nhitsi, nhitsj;
	
	int nhits_common = 0;	

	for(i_trd = 0; i_trd<THREADS_PER_BLOCK; i_trd++){
		if(i_trd==threadIdx.x)continue;
		gTracks Tracks_trd = tklcoll->tracks(blockIdx.x, i_trd, Ntracks_trd);
	
		for(i = 0; i<Ntracks; i++){
		//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
			if(Tracks.stationID(i)<5)continue;
			nhitsi = Tracks.nHits(i);
			nhits_common = 0;
		
			for(j = 0; j<Ntracks_trd; j++){
			//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
				if(Tracks_trd.stationID(j)<5)continue;
				nhitsj = Tracks_trd.nHits(j);
				
				for(ihit = 0; ihit<nhitsi; ihit++){
					for(jhit = 0; jhit<nhitsj; jhit++){
						if( Tracks_trd.hits_detid(j, jhit)==Tracks.hits_detid(i, ihit) && Tracks_trd.hits_chan(j, jhit)==Tracks.hits_chan(i, ihit) ){	
							nhits_common++;
							break;
						}
					}
				}
			
				if(Tracks.chisq(i)/(nhitsi-4) < Tracks_trd.chisq(j)/(nhitsj-4)){
					if(nhits_common>nhitsi*0.3333f ){
						tklcoll->setStationID(tkl_coll_offset+array_thread_offset[i_trd], j, 3);
					}
			
				}else{
					if(nhits_common>nhitsj*0.3333f ){
						tklcoll->setStationID(tkl_coll_offset+array_thread_offset[threadIdx.x], i, 3);
					}
				} 
			}
		}
	}
	
}


// --------------------------------------------- //
//
// proportional segment mtaching after global tracking
// 
// --------------------------------------------- //

__global__ void gKernel_PropSegmentMatching(
	gEventHitCollections* hitcolls,
	gEventTrackCollection* tklcoll,
	const float* z_array,
	const bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x])return;
	
	int i, j;

	float x0, tx;
	float y0, ty;
	
	float invP;
	
	//proportional tubes
	short stid, projid, detid;
	short nprop;
	float xExp, yExp, ipos;
	bool checknext;
	float prop_pos[2];
	float prop_z[2];
	float a, b;
	bool goodsegmentx;
	bool goodsegmenty;

	int n, m;
	
	projid = 0;
	stid = 6-1;
	detid = globalconsts::detsuperid[stid][projid]*2;
	int nhits_p1x1;
	const gHits hits_p1x1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1x1);
	const float z_p1x1 = z_array[detid];

	detid-= 1;
	int nhits_p1x2;
	const gHits hits_p1x2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1x2);
	const float z_p1x2 = z_array[detid];
	
	stid = 7-1;
	detid = globalconsts::detsuperid[stid][projid]*2;
	int nhits_p2x1;
	const gHits hits_p2x1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2x1);
	const float z_p2x1 = z_array[detid];

	detid-= 1;
	int nhits_p2x2;
	const gHits hits_p2x2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2x2);
	const float z_p2x2 = z_array[detid];
	
	// proportional tube hits
	projid = 1;
	stid = 6-1;
	detid = globalconsts::detsuperid[stid][projid]*2;
	int nhits_p1y1;
	const gHits hits_p1y1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1y1);
	const float z_p1y1 = z_array[detid];

	detid-= 1;
	int nhits_p1y2;
	const gHits hits_p1y2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1y2);
	const float z_p1y2 = z_array[detid];
	
	stid = 7-1;
	detid = globalconsts::detsuperid[stid][projid]*2;
	int nhits_p2y1;
	const gHits hits_p2y1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2y1);
	const float z_p2y1 = z_array[detid];

	detid-= 1;
	int nhits_p2y2;
	const gHits hits_p2y2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2y2);
	const float z_p2y2 = z_array[detid];
	
	//tracks
	unsigned int Ntracks;
	const gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);
	
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	float cut = 0.03;
	
	for(i = 0; i<Ntracks; i++){
		x0 = Tracks.x0(i);
		y0 = Tracks.y0(i);
		tx = Tracks.tx(i);
		ty = Tracks.ty(i);
		
		if(Tracks.stationID(i)>=6){	
			invP = Tracks.invP(i);
			//cut = max(invP*0.11825f, 0.00643f-0.00009f/invP+0.00000046f/invP/invP);
			cut = 0.03f;
		}
		
		goodsegmentx = false;
		// *very rough* tracklet segment
		// loop on first plane first
		for(n = 0; n<nhits_p1x1; n++){
			ipos = hits_p1x1.pos(n);
			xExp = tx*z_p1x1+x0;
			prop_pos[0] = ipos;
			prop_z[0] = z_p1x1;
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("evt %d p1x1, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
			if(fabs(ipos-xExp)<5.08f){
				nprop++;
				for(m = 0; m<nhits_p2x1; m++){
					ipos = hits_p2x1.pos(m);
					xExp = tx*z_p2x1+x0;
					prop_pos[1] = ipos;
					prop_z[1] = z_p2x1;
#ifdef DEBUG
					if(blockIdx.x<=debug::EvRef)printf("evt %d p2x1, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
					if(fabs(ipos-xExp)>7.62f)continue;

					a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
					//if(blockIdx.x<=debug::EvRef && st3==3 && hits_p1x1.chan(n)==21 && hits_p2x1.chan(m)==19)printf("a %1.4f tx %1.4f, chan p1x1 %1.0f pos %1.4f, chan p2x1 %1.0f  pos %1.4f \n", a, tx, hits_p1x1.chan(n), prop_pos[0], hits_p2x1.chan(m), prop_pos[1] );
					if(fabs(a)>TX_MAX)continue;
					if(fabs(a-tx)>cut)continue;
					b = prop_pos[0]-a*prop_z[0];
					//if(blockIdx.x<=debug::EvRef && st3==3 && hits_p1x1.chan(n)==21 && hits_p2x1.chan(m)==19)printf("a %1.4f tx %1.4f, b %1.4f, x0 %1.4f, chan p1x1 %1.0f, chan p2x1 %1.0f \n", a, tx, b, x0, hits_p1x1.chan(n), hits_p2x1.chan(m) );
					if(fabs(b)>X0_MAX)continue;
					if(fabs( (tx*2028.19f+x0)-(a*2028.19f+b) )>3.0f )continue;
					//if(blockIdx.x<=debug::EvRef && threadIdx.x%2==0 && hits_p1x1.chan(n)==21 && hits_p2x1.chan(m)==19)printf("thread %d a %1.4f tx %1.4f, b %1.4f, x0 %1.4f, pos_tkl %1.4f pos_seg %1.4f, detID[0] %d, detID[2] %d elID[0-3] %d %d %d %d \n", threadIdx.x, a, tx, b, x0, tx*2028.19f+x0, a*2028.19f+b, detID[0], detID[2], elID[0], elID[1], elID[2], elID[3] );
					
					nprop++;
					checknext = false;
					goodsegmentx = true;
					break;
				}
				if(!goodsegmentx){
					for(m = 0; m<nhits_p2x2; m++){
					ipos = hits_p2x2.pos(m);
					xExp = tx*z_p2x2+x0;
					prop_pos[1] = ipos;
					prop_z[1] = z_p2x2;
#ifdef DEBUG
					if(blockIdx.x<=debug::EvRef)printf("evt %d p2x2, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
					if(fabs(ipos-xExp)>7.62f)continue;

					a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
					if(fabs(a)>TX_MAX)continue;
					if(fabs(a-tx)>cut)continue;
					b = prop_pos[0]-a*prop_z[0];
					if(fabs(b)>X0_MAX)continue;
					if(fabs( (tx*2028.19f+x0)-(a*2028.19f+b) )>3.0f )continue;
						
					nprop++;
					checknext = false;
					goodsegmentx = true;
					break;
					}
				}
			}
		}
		
		if(!goodsegmentx){
			for(n = 0; n<nhits_p1x2; n++){
				ipos = hits_p1x2.pos(n);
				xExp = tx*z_p1x2+x0;
				prop_pos[0] = ipos;
				prop_z[0] = z_p1x2;
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("evt %d p1x1, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
				if(fabs(ipos-xExp)<5.08f){
					nprop++;
					for(m = 0; m<nhits_p2x1; m++){
						ipos = hits_p2x1.pos(m);
						xExp = tx*z_p2x1+x0;
						prop_pos[1] = ipos;
						prop_z[1] = z_p2x1;
#ifdef DEBUG
						if(blockIdx.x<=debug::EvRef)printf("evt %d p2x1, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
						if(fabs(ipos-xExp)>7.62f)continue;
	
						a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
						if(fabs(a)>TX_MAX)continue;
						if(fabs(a-tx)>cut)continue;
						b = prop_pos[0]-a*prop_z[0];
						if(fabs(b)>X0_MAX)continue;
						if(fabs( (tx*2028.19f+x0)-(a*2028.19f+b) )>3.0f )continue;
							
						nprop++;
						checknext = false;
						goodsegmentx = true;
						break;
					}
					if(!goodsegmentx){
						for(m = 0; m<nhits_p2x2; m++){
						ipos = hits_p2x2.pos(m);
						xExp = tx*z_p2x2+x0;
						prop_pos[1] = ipos;
						prop_z[1] = z_p2x2;
#ifdef DEBUG
						if(blockIdx.x<=debug::EvRef)printf("evt %d p2x2, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
						if(fabs(ipos-xExp)>7.62f)continue;
	
						a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
						if(fabs(a)>TX_MAX)continue;
						if(fabs(a-tx)>cut)continue;
						b = prop_pos[0]-a*prop_z[0];
						if(fabs(b)>X0_MAX)continue;
						if(fabs( (tx*2028.19f+x0)-(a*2028.19f+b) )>3.0f )continue;
							
						nprop++;
						checknext = false;
						goodsegmentx = true;
						break;
						}
					}
				}
			}
		}
		if(!goodsegmentx){
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 5);
			continue;
		}
		
		goodsegmenty = false;
		// *very rough* tracklet segment
		// loop on first plane first
		for(n = 0; n<nhits_p1y1; n++){
			ipos = hits_p1y1.pos(n);
			yExp = ty*z_p1y1+y0;
			prop_pos[0] = ipos;
			prop_z[0] = z_p1y1;
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("evt %d p1y1, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
			if(fabs(ipos-yExp)<5.08f){
				nprop++;
				for(m = 0; m<nhits_p2y1; m++){
					ipos = hits_p2y1.pos(m);
					yExp = ty*z_p2y1+y0;
					prop_pos[1] = ipos;
					prop_z[1] = z_p2y1;
#ifdef DEBUG
					if(blockIdx.x<=debug::EvRef)printf("evt %d p2y1, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
					if(fabs(ipos-yExp)>10.16f)continue;

					a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
					if(fabs(a)>TY_MAX)continue;
					if(fabs(a-ty)>cut)continue;
					b = prop_pos[0]-a*prop_z[0];
					if(fabs(b)>Y0_MAX)continue;
					if(fabs( (ty*2028.19f+y0)-(a*2028.19f+b) )>3.0f )continue;
				
					nprop++;
					checknext = false;
					goodsegmenty = true;
					break;
				}
				if(!goodsegmenty){
					for(m = 0; m<nhits_p2y2; m++){
						ipos = hits_p2y2.pos(m);
						yExp = ty*z_p2y2+y0;
						prop_pos[1] = ipos;
						prop_z[1] = z_p2y2;
#ifdef DEBUG
						if(blockIdx.x<=debug::EvRef)printf("evt %d p2y2, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
						if(fabs(ipos-yExp)>10.16f)continue;
		
						a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
						if(fabs(a)>TY_MAX)continue;
						if(fabs(a-ty)>cut)continue;
						b = prop_pos[0]-a*prop_z[0];
						if(fabs(b)>Y0_MAX)continue;
						if(fabs( (ty*2028.19f+y0)-(a*2028.19f+b) )>3.0f )continue;
							
						nprop++;
						checknext = false;
						goodsegmenty = true;
						break;
					}
				}
			}
		}
			
		if(!goodsegmenty){
			for(n = 0; n<nhits_p1y2; n++){
				ipos = hits_p1y2.pos(n);
				yExp = ty*z_p1y2+y0;
				prop_pos[0] = ipos;
				prop_z[0] = z_p1y2;
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("evt %d p1y1, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
				if(fabs(ipos-yExp)<5.08f){
					nprop++;
					for(m = 0; m<nhits_p2y1; m++){
						ipos = hits_p2y1.pos(m);
						yExp = ty*z_p2y1+y0;
						prop_pos[1] = ipos;
						prop_z[1] = z_p2y1;
#ifdef DEBUG
						if(blockIdx.x<=debug::EvRef)printf("evt %d p2y1, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
						if(fabs(ipos-yExp)>10.16f)continue;
		
						a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
						if(fabs(a)>TY_MAX)continue;
						if(fabs(a-ty)>cut)continue;
						b = prop_pos[0]-a*prop_z[0];
						if(fabs(b)>Y0_MAX)continue;
						if(fabs( (ty*2028.19f+y0)-(a*2028.19f+b) )>3.0f )continue;
							
						nprop++;
						checknext = false;
						goodsegmenty = true;
						break;
					}
					if(!goodsegmenty){
						for(m = 0; m<nhits_p2y2; m++){
							ipos = hits_p2y2.pos(m);
							yExp = ty*z_p2y2+y0;
							prop_pos[1] = ipos;
							prop_z[1] = z_p2y2;
#ifdef DEBUG
							if(blockIdx.x<=debug::EvRef)printf("evt %d p2y2, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
							if(fabs(ipos-yExp)>10.16f)continue;
	
							a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
							if(fabs(a)>TY_MAX)continue;
							if(fabs(a-ty)>cut)continue;
							b = prop_pos[0]-a*prop_z[0];
							if(fabs(b)>Y0_MAX)continue;
							if(fabs( (ty*2028.19f+y0)-(a*2028.19f+b) )>3.0f )continue;
								
							nprop++;
							checknext = false;
							goodsegmenty = true;
							break;
						}
					}
				}
			}
		}
		if(!goodsegmenty){
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 5);
			continue;
		}
	
	}
}



// --------------------------------------------- //
//
// function to clean the tracks after full processing
// 
// --------------------------------------------- //

__global__ void gKernel_GlobalTrackCleaning(
	gEventTrackCollection* tklcoll,
	const bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x])return;
	
	int i, j;

	float x0, tx;
	float y0, ty;
	
	float invP;
		
	//tracks
	unsigned int Ntracks;
	const gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);
	
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	float cut = 0.03;
	
	for(i = 0; i<Ntracks; i++){
		//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
		if(Tracks.stationID(i)<6)continue;
		
		if(Tracks.chisq(i)/(Tracks.nHits(i)-5)>selection::chi2dofmax){//  || Tracks.chisq(i)>100.f
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 5);
			continue;
		}
		
		for(j = i+1; j<Ntracks; j++){
		//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
			if(Tracks.stationID(j)<6)continue;
			
			
			if(Tracks.chisq(i)/(Tracks.nHits(i)-5) < Tracks.chisq(j)/(Tracks.nHits(j)-5)){
				//"merging" tracks with similar momentum, but only if their charge is of same sign
				if( (Tracks.invP(i)-Tracks.invP(j))/Tracks.invP(i) < selection::merge_thres && Tracks.charge(i)*Tracks.charge(j)>0 )  tklcoll->setStationID(tkl_coll_offset+array_thread_offset, j, 5);
			}else{
				if( (Tracks.invP(i)-Tracks.invP(i))/Tracks.invP(j) < selection::merge_thres && Tracks.charge(i)*Tracks.charge(j)>0 )  tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 5);
			} 

		}
	}
}

__global__ void gKernel_GlobalTrackCleaning_crossthreads(
	gEventTrackCollection* tklcoll,
	const bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x])return;
	
	//tracks
	unsigned int Ntracks, Ntracks_trd;
	gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);
	unsigned int array_thread_offset[THREADS_PER_BLOCK];// = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	int i_trd;
	for(i_trd = 0; i_trd<THREADS_PER_BLOCK; i_trd++){
		//Tracks[i_trd] = tklcoll->tracks(blockIdx.x, i_trd, Ntracks[i_trd]);
		array_thread_offset[i_trd] = i_trd*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	}
	
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	
	int i, j; 
	short ihit, jhit;
	short nhitsi, nhitsj;
	
	int nhits_common = 0;	

	for(i_trd = 0; i_trd<THREADS_PER_BLOCK; i_trd++){
		if(i_trd==threadIdx.x)continue;
		gTracks Tracks_trd = tklcoll->tracks(blockIdx.x, i_trd, Ntracks_trd);
	
		for(i = 0; i<Ntracks; i++){
		//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
			if(Tracks.stationID(i)<6)continue;
						
			nhitsi = Tracks.nHits(i);
			nhits_common = 0;
		
			for(j = 0; j<Ntracks_trd; j++){
			//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
				if(Tracks_trd.stationID(j)<6)continue;
				
				if(Tracks.chisq(i)/(Tracks.nHits(i)-5) < Tracks_trd.chisq(j)/(Tracks_trd.nHits(j)-5)){
					//"merging" tracks with similar momentum, but only if their charge is of same sign
					if( (Tracks.invP(i)-Tracks_trd.invP(j))/Tracks.invP(i) < selection::merge_thres && Tracks.charge(i)*Tracks_trd.charge(j)>0 )  tklcoll->setStationID(tkl_coll_offset+array_thread_offset[i_trd], j, 5);
				}else{
					if( (Tracks.invP(i)-Tracks_trd.invP(j))/Tracks_trd.invP(j) < selection::merge_thres && Tracks.charge(i)*Tracks_trd.charge(j)>0 )  tklcoll->setStationID(tkl_coll_offset+array_thread_offset[threadIdx.x], i, 5);
				}

			}
		}
	}
	
}


// --------------------------------------------- //
//
// simple track printing function for debugging. 
// 
// --------------------------------------------- //
__global__ void gKernel_check_tracks(gEventTrackCollection* tklcoll, const bool* hastoomanyhits, const int blockID)
{
	if(hastoomanyhits[blockIdx.x])return;
	
	unsigned int nTracks;
	
	if(blockIdx.x==blockID){
		const gTracks Tracks = tklcoll->tracks(blockID, threadIdx.x, nTracks);
		
		for(int i = 0; i<nTracks; i++){
			printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
		}
	}
	
}


__global__ void gkernel_checkHistos(gHistsArrays *HistsArrays, const short hist_num){
	const int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if(blockIdx.x==hist_num || hist_num==-1){
		if(threadIdx.x==0)printf(" hw %1.4f ", HistsArrays->pts_hw[blockIdx.x]);
		printf(" thread %d idx %d xpts %1.4f values %1.4f\n", threadIdx.x, idx, HistsArrays->xpts[idx], HistsArrays->values[idx]);
	}
}



#ifdef EXTRASTUFF

__device__ float delta_tx(const float delta_x0)
{
  return -1.41477e-05-0.00054203*delta_x0;
}
__device__ float delta_ty(const float delta_y0)
{
  return -1.71247e-05-0.00071988*delta_y0;
}
// float delta_ty(const float delta_y0)
// {
//   return -0.001*delta_y0;
// }
__device__ float invp_ratio(const float delta_x0, const short charge)
{
  if(charge>0){
    return 1.00317-0.1043*delta_x0;
  }else{
    return 0.996416+0.1032*delta_x0;
  }
}

__device__ void adjust_track_parameters(size_t const n_points, float* residuals,
				short* const detid, float* const drift, short* const sign, float* const resolutions,
				float* const p1x, float* const p1y, float* const p1z,
				float* const deltapx, float* const deltapy, float* const deltapz,
				const float x0, const float y0, const float tx, const float ty, const float invp, const short charge,
				float &x0_new, float &y0_new, float &tx_new, float &ty_new, float &invp_new)
{
	short nn, mm, l;
	float x0_mod, tx_mod, y0_mod, ty_mod, invp_mod, x0_st1_mod, tx_st1_mod;
 	float dx0, dy0;
	float chi2, chi2_min;
		
	for(nn = 0; nn<=100; nn++){
		dx0 = deltax0_+deltax0_sigma*(gauss_quantiles[nn]);
		x0_mod = x0+dx0;
		tx_mod = tx+delta_tx(dx0);
		invp_mod = invp/invp_ratio(dx0, charge);
		calculate_x0_tx_st1(x0_mod, tx_mod, invp_mod, charge, x0_st1_mod, tx_st1_mod);
		
		chi2 = 0;
		
		for(l = 0; l<n_points; l++){
		  if(detid[l]<12){
		    residuals[l] = residual(p1x[l], p1y[l], p1z[l], deltapx[l], deltapy[l], deltapz[l], drift[l], sign[l], x0_st1_mod, y0, tx_st1_mod, ty);
		     //if(fabs(residuals[l])>fabs(residuals[l]-2*sign[l]*drift[l]))residuals[l] = resid_mod[l]-2*sign[l]*drift[l];
		  }else{
		    residuals[l] = residual(p1x[l], p1y[l], p1z[l], deltapx[l], deltapy[l], deltapz[l], drift[l], sign[l], x0_mod, y0, tx_mod, ty);
		    //if(fabs(residuals[l])>fabs(residuals[l]-2*sign[l]*drift[l]))residuals[l] = resid_mod[l]-2*sign[l]*drift[l];
		  }
		  chi2+= residuals[l]*residuals[l]/resolutions[l]/resolutions[l];
		}
		
		if(chi2<chi2_min){
			chi2_min = chi2;
		    	x0_new = x0_mod;
			tx_new = tx_mod;
			invp_new = invp_mod;
		}
	}

	calculate_x0_tx_st1(x0_new, tx_new, invp_new, charge, x0_st1_mod, tx_st1_mod);
	
	for(nn = 0; nn<=100; nn++){
		dy0 = deltay0_+deltay0_sigma*(gauss_quantiles[nn]);
		y0_mod = x0+dy0;
		
		for(int mm = 0; mm<=10; mm++){
	    		ty_mod = ty+delta_ty(dy0)+gauss_quantiles[mm*10]*0.0005f;
			chi2 = 0;
			
			for(l = 0; l<n_points; l++){
				if(detid[l]<12){
					residuals[l] = residual(p1x[l], p1y[l], p1z[l], deltapx[l], deltapy[l], deltapz[l], drift[l], sign[l], x0_st1_mod, y0_mod, tx_st1_mod, ty_mod);
					//if(fabs(residuals[l])>fabs(residuals[l]-2*sign[l]*drift[l]))residuals[l] = resid_mod[l]-2*sign[l]*drift[l];
				}else{
		    			residuals[l] = residual(p1x[l], p1y[l], p1z[l], deltapx[l], deltapy[l], deltapz[l], drift[l], sign[l], x0_mod, y0_new, tx_new, ty_mod);
					//if(fabs(residuals[l])>fabs(residuals[l]-2*sign[l]*drift[l]))residuals[l] = resid_mod[l]-2*sign[l]*drift[l];
				}
				chi2+= residuals[l]*residuals[l]/resolutions[l]/resolutions[l];
			}
		
			if(chi2<chi2_min){
				chi2_min = chi2;
			    	y0_new = y0_mod;
				ty_new = ty_mod;
			}
		}
	}

	
	for(l = 0; l<n_points; l++){
		if(detid[l]<12){
			residuals[l] = residual(p1x[l], p1y[l], p1z[l], deltapx[l], deltapy[l], deltapz[l], drift[l], sign[l], x0_st1_mod, y0_new, tx_st1_mod, ty_mod);
		}else{
			residuals[l] = residual(p1x[l], p1y[l], p1z[l], deltapx[l], deltapy[l], deltapz[l], drift[l], sign[l], x0_new, y0_new, tx_new, ty_mod);
		}
	}
}


// -------------------------------------- //
// functions to check the track at target //
// -------------------------------------- //

__device__ float tx_tgt_(const float tx, const float invp, const short charge)
{
	return -charge*(0.0178446+2.56307*invp)+tx;
}

__device__ float x0_tgt_(const float x0, const float invp, const short charge)
{
	return charge*(3.78335+1010.71*invp)+x0;
}

__device__ float z_min_(const float x0_tgt, const float tx_tgt, const float y0, const float ty)
{	
	return -(x0_tgt*tx_tgt + y0*ty)/(tx_tgt*tx_tgt + ty*ty);
}


__device__ bool check_target_pointing_quick(const float x0, const float tx, const float y0, const float ty, const float invp, const short charge)
{
	float tx_tgt = tx_tgt_(tx, invp, charge);
	float x0_tgt = x0_tgt_(x0, invp, charge);
	
	float z_min = z_min_(x0_tgt, tx_tgt, y0, ty);
	float x_tgt = x_trk(x0_tgt, tx_tgt, z_min);
	float y_tgt = y_trk(y0, ty, z_min);
	
	if(fabs(x_tgt)<=selection::x_vtx_cut && fabs(y_tgt)<=selection::y_vtx_cut && fabs(tx_tgt)<=selection::tx_vtx_cut && fabs(ty)<=selection::ty_vtx_cut){
		return true;
	}
	return false;
}

#endif


#ifdef KALMAN_TRACKING

// --------------------------------------------- //
// functions to calculate x0 and tx in station 1 //
// and function to calculate inverse momentum    //
// --------------------------------------------- //


__device__ void extrapolate_track_position_st1(gTracklet& tkl, float* x_st1_mean, float* x_st1_width, const gPlane* planes, const short hyp)
{
	//1st order;
	tkl.invP = extrapolation_tools::invP_x0_[hyp][0]+fabs(tkl.x0)*extrapolation_tools::invP_x0_[hyp][1];
	tkl.err_invP = extrapolation_tools::err_invP_x0[hyp];

	//printf("x0 %1.6f, p0 %1.6f p1 %1.6f => invP %1.6f\n", tkl.x0, extrapolation_tools::invP_x0_[0][0], extrapolation_tools::invP_x0_[0][1], tkl.invP);
	
	float x_st1_trk_diff, dx_st1_trk_diff;
	bool xpos = tkl.x0>0;
	
	for(int i = 3; i<=4; i++){
		x_st1_trk_diff = extrapolation_tools::straight_st1_det_extrap[i-3][xpos][0]+tkl.invP*extrapolation_tools::straight_st1_det_extrap[i-3][xpos][1];
		dx_st1_trk_diff = tkl.err_invP*extrapolation_tools::straight_st1_det_extrap[i-3][xpos][1];
		//if(xpos)printf(" %d %1.6f %1.6f %1.6f \n", i, x_st1_trk_diff, extrapolation_tools::straight_st1_det_extrap[i-3][xpos][0], extrapolation_tools::straight_st1_det_extrap[i-3][xpos][1]);
		x_st1_mean[i-3] = x_st1_trk_diff+tkl.x0+tkl.tx*planes[i].z;
		x_st1_width[i-3] = dx_st1_trk_diff+tkl.err_x0+tkl.err_tx*planes[i].z+planes[i].spacing;
	}
}


// ----------------------------------- //
// matrix operations for Kalman filter //
// ----------------------------------- //

__device__ float similarity_1x2_S5x5_2x1(const float* M_2x1, const float* M_S5x5)
{
	return( M_2x1[0]*M_2x1[0]*M_S5x5[0] + M_2x1[0]*M_2x1[1]*M_S5x5[1] + M_2x1[1]*M_2x1[1]*M_S5x5[6]);
}

__device__ void multiply_S5x5_2x1(const float* M_2x1, const float* M_S5x5, float* M_5x1)
{
	for(int k = 0; k<5; k++){
		M_5x1[k] = M_S5x5[k]*M_2x1[0]+M_S5x5[5+k] * M_2x1[1];
		//printf("k = %d M5x5_k = %1.6f, M5x5_5+k = %1.6f \n", k, M_S5x5[k], M_S5x5[5+k]);
	}
}

__device__ void tensor_product(const float* M_5x1_1, const float* M_5x1_2, float* M_5x5, const float Cnorm)
{
	for(short i = 0; i<5; i++){
		for(short j = 0; j<5; j++){
			M_5x5[i+5*j] = M_5x1_1[i] * M_5x1_2[j] * Cnorm;
		}
	}
}


__device__ void scale_matrix(float* M1, const float C, const short size)
{
	for(short i = 0; i<size; i++){
		M1[i]*= C;
	}
}

__device__ void add_matrix(float* M1, const float* M2, const short size, const float C = 1.f)
{
	for(short i = 0; i<size; i++){
		M1[i] += M2[i]*C;
	}
}


// ------------------------------ //
// functions for Kalman filtering // 
// ------------------------------ //

__device__ void initialize_arrays(gTracklet& tkl, gKalmanFitArrays &fitarray)
{
	fitarray.state[0] = tkl.x0;
	fitarray.state[1] = tkl.y0;
	fitarray.state[2] = tkl.tx;
	fitarray.state[3] = tkl.ty;
	fitarray.state[4] = tkl.invP;
	
	//initialize covariance matrix:

	fitarray.Cov[0] = 100.f;                //c00
	fitarray.Cov[1] = fitarray.Cov[5] = 0.f; //c01
	fitarray.Cov[2] = fitarray.Cov[10] = 0.f; //c02
	fitarray.Cov[3] = fitarray.Cov[15] = 0.f; //c03
	fitarray.Cov[4] = fitarray.Cov[20] = 0.f; //c03
	
	fitarray.Cov[6] = 100.f;                 //c11
	fitarray.Cov[7] = fitarray.Cov[11] = 0.f; //c12
	fitarray.Cov[8] = fitarray.Cov[16] = 0.f; //c13
	fitarray.Cov[9] = fitarray.Cov[21] = 0.f; //c14
	
	fitarray.Cov[12] = 0.01f;                 //c22
	fitarray.Cov[13] = fitarray.Cov[17] = 0.f; //c23
	fitarray.Cov[14] = fitarray.Cov[22] = 0.f; //c24
	
	fitarray.Cov[18] = 0.01f;                  //c33
	fitarray.Cov[19] = fitarray.Cov[23] = 0.f; //c34
	
	fitarray.Cov[24] = 0.09f*tkl.invP*tkl.invP; //c44
}

__device__ void predict_state(const gTracklet tkl, const float z, gKalmanFitArrays& fitarray)//float* pred_state)
//if we do it this way we have to make sure we update the tracklet parameters at each step
{
	fitarray.state[1] = tkl.y0+z*tkl.ty;
	fitarray.state[3] = tkl.ty;
	fitarray.state[4] = tkl.invP;

	if(z<globalconsts::Z_KMAG_BEND){
		float x0_st1, tx_st1;
		calculate_x0_tx_st1(tkl, x0_st1, tx_st1);

		fitarray.state[0] = x0_st1+z*tx_st1;
		fitarray.state[2] = tx_st1;
	}else{
		fitarray.state[0] = tkl.x0+z*tkl.tx;
		fitarray.state[2] = tkl.tx;
	}
}

__device__ void update_tkl(gTracklet& tkl, const float z, gKalmanFitArrays& fitarray)//const float* updated_state)
{
	float x0, tx;
	float x0_st1, tx_st1;
		
	//update not only x but also invP with new info as well:
	//get the other part of the 
        if(z<globalconsts::Z_KMAG_BEND){
		x0 = tkl.x0;
		tx = tkl.tx;
		x0_st1 = fitarray.state[0]-z*fitarray.state[2];
		tx_st1 = fitarray.state[2];
	}else{
		calculate_x0_tx_st1(tkl, x0_st1, tx_st1);
		x0 = fitarray.state[0]-z*fitarray.state[2];
		tx = fitarray.state[2];
	}
	
	tkl.x0 = x0;
        tkl.tx = tx;
	
	tkl.y0 = fitarray.state[1]-z*fitarray.state[3];
	tkl.ty = fitarray.state[3];
	
	tkl.invP = calculate_invP(tx, tx_st1);	
	
}

__device__ void update_state(gTracklet& tkl, const gHit hit, gKalmanFitArrays& fitarray, const gPlane plane)//will optimize after if something is redundant
{
	const float dxdy = plane.deltapx/plane.deltapy;
	const float y0 = y_bep(hit, plane);
	const float y1 = y_tep(hit, plane);
	const float x0 = hit.pos*plane.costheta + y0*dxdy;
	const float x1 = x0 + (y1-y0) *dxdy;
	const float z = plane.z;

	const float x2y2 = sqrtf( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) );
	
	//printf("x0 %1.6f x1 %1.6f y0 %1.6f y1 %1.6f x2y2 %1.6f \n", x0, x1, y0, y1, x2y2);

	fitarray.H[0] = (y1-y0) / x2y2;
	fitarray.H[1] = (x1-x0) / x2y2;
	
	const float res = fitarray.H[0] * x0 + fitarray.H[1] * y0 - (fitarray.H[0] * fitarray.state[0] + fitarray.H[1] * fitarray.state[1]); 
	const float CRes = similarity_1x2_S5x5_2x1(fitarray.H, fitarray.Cov) + plane.resolution * plane.resolution;

	//printf("h0 %1.6f h1 %1.6f res %1.6f CRes %1.6f \n", fitarray.H[0], fitarray.H[1], res, CRes);
		
	multiply_S5x5_2x1(fitarray.H, fitarray.Cov, fitarray.K);
	
	scale_matrix(fitarray.K, 1./CRes, 5);
	for(int k = 0; k<5; k++){
		printf("%d, %1.6f, %1.6f\n", k, fitarray.state[k], fitarray.K[k]);
	}
	add_matrix(fitarray.state, fitarray.K, 5, res);
	
	tensor_product(fitarray.K, fitarray.K, fitarray.KCResKt, CRes);
	
	add_matrix(fitarray.Cov, fitarray.KCResKt, 25, -1.f);
	
	fitarray.chi2 = res*res / CRes;
	
	// now update the tracklet parameters!
	
	//update_tkl(tkl, z, fitarray.state);
	update_tkl(tkl, z, fitarray);
}

#endif


#ifdef LEGACYCODE

__device__ void make_hitpairs_in_station_bins(const gHits hitcoll1, const int nhits1, const gHits hitcoll2, const int nhits2, thrust::pair<int, int>* hitpairs, int* npairs, const short bin0, const short Nbins, short* hitflag1, short* hitflag2, const int stID, const int projID){
	// I think we assume that by default we want to know where we are
	//printf("stID %d projID %d bin0 %d\n", stID, projID, bin0);
	
	short bin;
	const short MaxHits = globalconsts::MaxHitsProj[projID];

	for(bin = bin0; bin<bin0+Nbins; bin++){
		npairs[bin-bin0] = 0;
	}
	
	//declaring arrays for the hit lists
	for(int i = 0; i<100; i++){
		hitflag1[i] = hitflag2[i] = 0;
	}
	
	//building the lists of hits for each detector plane
	//const int detid1 = globalconsts::detsuperid[stID][projID]*2;
	//const int detid2 = globalconsts::detsuperid[stID][projID]*2-1;
	const int superdetid = globalconsts::detsuperid[stID][projID];
	
	// pair the hits by position:
	// if one hit on e.g. x and one hit on x' are closer than
	// the "spacing" defined for the planes, then the hits can be paired together.
	int idx1 = -1;
	int idx2 = -1;
	
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0)printf("nhits %d %d \n", nhits1, nhits2);
#endif
	for(int i = 0; i<nhits1; i++){
		idx1++;
		idx2 = -1;
		for(int j = 0; j<nhits2; j++){
			idx2++;
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef && threadIdx==0)printf("i %d j %d pos %1.4f %1.4f\n", i, j, );
#endif
			if( abs(hitcoll1.pos(idx1) - hitcoll2.pos(idx2)) > globalconsts::spacingplane[superdetid] ){
				continue;
			}
			
			for(bin = bin0; bin<bin0+Nbins; bin++){
				if( globalconsts::WCHitsBins[stID-1][projID][0][bin] <= hitcoll1.chan(idx1) && 
				    hitcoll1.chan(idx1) <= globalconsts::WCHitsBins[stID-1][projID][1][bin]){
					//printf("bin %d low %d high %d hit 1 elem %d hit 2 elem %d global bin %d \n", bin, globalconsts::WCHitsBins[stID-1][projID][0][bin-bin0], globalconsts::WCHitsBins[stID-1][projID][1][bin-bin0], ic[index].AllHits[ i ].elementID, ic[index].AllHits[ idx2 ].elementID, bin+npairs[bin]*Nbins);
					if(npairs[bin-bin0]<=MaxHits)hitpairs[bin-bin0+npairs[bin-bin0]*Nbins] = thrust::make_pair(idx1, idx2);
					npairs[bin-bin0]++;
				}
			}
			hitflag1[idx1] = 1;
			hitflag2[idx2] = 1;
		}
	}
	// here the hits that cannot be paired to another hit are paired to "nothing"
	// (but they still have to be paired to be used in the trackletteing)
	for(int i = 0; i<nhits1; i++){
		if(hitflag1[i]<1){
			for(bin = bin0; bin<bin0+Nbins; bin++){
			//printf("bin %d low %d high %d hit elem %d global bin %d \n", bin, globalconsts::WCHitsBins[stID-1][projID][0][bin], globalconsts::WCHitsBins[stID-1][projID][1][bin], ic[index].AllHits[ i ].elementID, bin+npairs[bin]*Nbins);
				if( globalconsts::WCHitsBins[stID-1][projID][0][bin] <= hitcoll1.chan(i) && 
				    hitcoll1.chan(i) <= globalconsts::WCHitsBins[stID-1][projID][1][bin]){
					//printf("bin %d low %d high %d hit elem %d global bin %d \n", bin, globalconsts::WCHitsBins[stID-1][projID][0][bin-bin0], globalconsts::WCHitsBins[stID-1][projID][1][bin-bin0], ic[index].AllHits[ i ].elementID, bin+npairs[bin]*Nbins);
					if(npairs[bin-bin0]<=MaxHits)hitpairs[bin-bin0+npairs[bin-bin0]*Nbins] = thrust::make_pair(i, -1);
					npairs[bin-bin0]++;
				}
			}
		}
	}
	for(int i = 0; i<nhits2; i++){
		if(hitflag2[i]<1){
			for(bin = bin0; bin<bin0+Nbins; bin++){
			//printf("bin %d low %d high %d hit elem %d global bin %d \n", bin, globalconsts::WCHitsBins[stID-1][projID][0][bin], globalconsts::WCHitsBins[stID-1][projID][1][bin], ic[index].AllHits[ i ].elementID, bin+npairs[bin]*Nbins);
				if( globalconsts::WCHitsBins[stID-1][projID][0][bin] <= hitcoll2.chan(i) && 
				    hitcoll2.chan(i) <= globalconsts::WCHitsBins[stID-1][projID][1][bin]){
					//printf("bin %d low %d high %d hit elem %d global bin %d \n", bin, globalconsts::WCHitsBins[stID-1][projID][0][bin-bin0], globalconsts::WCHitsBins[stID-1][projID][1][bin-bin0], ic[index].AllHits[ i ].elementID, bin+npairs[bin]*Nbins);
					if(npairs[bin-bin0]<=MaxHits)hitpairs[bin-bin0+npairs[bin-bin0]*Nbins] = thrust::make_pair(-1, i);
					npairs[bin-bin0]++;
				}
			}
		 }
	}
}

// Hit comparison
struct lessthan {
	__host__ __device__ bool operator()(const gHit& lhs, const gHit& rhs)
	{
	//returns true if :
	// hit1.detID<hit2.detID;  
		if(lhs.detectorID < rhs.detectorID)
		{
			return true;
		}
		else if(lhs.detectorID > rhs.detectorID)
		{
			return false;
		}
	//hit1.detID=hit2.detID & hit1.elID<hit2.elID;
		if(lhs.elementID < rhs.elementID)
		{
			return true;
		}
		else if(lhs.elementID > rhs.elementID)
		{
			return false;
		}
	//hit1.detID=hit2.detID & hit1.elID=hit2.elID & hit1.time>hit2.time;
		if(lhs.tdcTime > rhs.tdcTime)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};


// ------------------------------------------------------------- //
// functions to evaluate the hit selection window for a 2D track //
// and to calculate y from u, v hits given x                     //
// ------------------------------------------------------------- //


__device__ void find_xmin_xmax_in_chamber(float &xmin, float &xmax, const gTrack2D track2d, const short stID, const short projID, const gPlane* planes)
{
	const int detid1 = globalconsts::detsuperid[stID][projID]*2;
	const int detid2 = globalconsts::detsuperid[stID][projID]*2-1;
	
	xmin = min(planes->z[detid1]*(track2d.tx_-track2d.err_tx_), planes->z[detid2]*(track2d.tx_-track2d.err_tx_));
	xmax = max(planes->z[detid1]*(track2d.tx_+track2d.err_tx_), planes->z[detid2]*(track2d.tx_+track2d.err_tx_));
	
	xmin = xmin + track2d.x_0-track2d.err_x_0-planes->spacing[detid1];
	xmax = xmax + track2d.x_0+track2d.err_x_0+planes->spacing[detid1];
}

__device__ void find_xmin_xmax_in_chamber(float &xmin, float &xmax, const gTracklet tkl, const short stID, const short projID, const gPlane* planes)
{
	const int detid1 = globalconsts::detsuperid[stID][projID]*2;
	const int detid2 = globalconsts::detsuperid[stID][projID]*2-1;
	
	xmin = min(planes->z[detid1]*(tkl.tx-tkl.err_tx), planes->z[detid2]*(tkl.tx-tkl.err_tx));
	xmax = max(planes->z[detid1]*(tkl.tx+tkl.err_tx), planes->z[detid2]*(tkl.tx+tkl.err_tx));
	
	xmin = xmin + tkl.x0-tkl.err_x0-planes->spacing[detid1];
	xmax = xmax + tkl.x0+tkl.err_x0+planes->spacing[detid1];
}




__device__ void FillChi2Arrays(const int n, const gHit hit, const short hitsign, gStraightFitArrays *fitarray, const gPlane* planes){
	const int index = threadIdx.x + blockIdx.x * blockDim.x;
	const int idxoff_global = index*MaxHitsPerTrack;	
	
	fitarray->drift_dist[n] = hit.driftDistance*hitsign;
	fitarray->resolution[n] = planes->resolution[ hit.detectorID ];
	if(hitsign==0){
		fitarray->resolution[n] = planes->spacing[ hit.detectorID ]*3.4641f;
	}else{
		fitarray->resolution[n] = planes->resolution[ hit.detectorID ];
	}	       
	fitarray->p1x[n] = x_bep( hit, planes);
	fitarray->p1y[n] = y_bep( hit, planes);
	fitarray->p1z[n] = z_bep( hit, planes);
	
	fitarray->deltapx[n] = planes->deltapx[ hit.detectorID ];
	fitarray->deltapy[n] = planes->deltapy[ hit.detectorID ];
	fitarray->deltapz[n] = planes->deltapz[ hit.detectorID ];
}

__device__ void FillChi2Arrays(const int n, const gHit hit, const short hitsign, float* drift_dist, float* resolution, float* p1x, float* p1y, float* p1z, float* deltapx, float* deltapy, float* deltapz, const gPlane* planes){
	const int index = threadIdx.x + blockIdx.x * blockDim.x;
	const int idxoff_global = index*MaxHitsPerTrack;	
	
	drift_dist[n] = hit.driftDistance*hitsign;
	resolution[n] = planes->resolution[ hit.detectorID ];
	if(hitsign==0){
		resolution[n] = planes->spacing[ hit.detectorID ]*3.4641f;
	}else{
		resolution[n] = planes->resolution[ hit.detectorID ];
	}	       
	p1x[n] = x_bep( hit, planes);
	p1y[n] = y_bep( hit, planes);
	p1z[n] = z_bep( hit, planes);
	
	deltapx[n] = planes->deltapx[ hit.detectorID ];
	deltapy[n] = planes->deltapy[ hit.detectorID ];
	deltapz[n] = planes->deltapz[ hit.detectorID ];
}


// ------------------------------------ //
// function to refit the track after    //
// left right ambiguity resolution      //
// ------------------------------------ //


__device__ float refit_backpartialtrack_with_drift(gTracklet& tkl, gStraightFitArrays* fitarray, const gPlane* planes){
	const int index = threadIdx.x + blockIdx.x * blockDim.x;
	const int idxoff_array = index*MaxHitsPerTrack;
	
	float A_[4];
	float Ainv_[4];
	float B_[2];
	float Par[4];
	float ParErr[4];
	float chi2;
	
	//X hits are stored first, so we fit them first;
	for(int i = 0; i<tkl.nXHits; i++){
		FillFitArrays_X(i, tkl.hits[i], tkl.hitsign[i], fitarray, planes);
		FillChi2Arrays(i, tkl.hits[i], tkl.hitsign[i], fitarray, planes);
	}
	fit_2D_track(tkl.nXHits, fitarray->x_array+idxoff_array, fitarray->z_array+idxoff_array, fitarray->dx_array+idxoff_array, A_, Ainv_, B_, Par, ParErr, chi2);

	tkl.x0 = Par[0];
	tkl.err_x0 = ParErr[0];
	
	tkl.tx = Par[1];
	tkl.err_tx = ParErr[1];
	
	float y, err_y;
	
	for(int i = tkl.nXHits; i<tkl.nXHits+tkl.nUHits+tkl.nVHits; i++){
		if( calculate_y_uvhit(y, err_y, tkl.hits[i], tkl.hitsign[i], tkl, planes) ){
			FillFitArrays_UV(i-tkl.nXHits, tkl.hits[i], fitarray, planes, y, err_y);
			FillChi2Arrays(i, tkl.hits[i], tkl.hitsign[i], fitarray, planes);
		}
	}
	fit_2D_track(tkl.nXHits, fitarray->x_array+idxoff_array, fitarray->z_array+idxoff_array, fitarray->dx_array+idxoff_array, A_, Ainv_, B_, Par, ParErr, chi2);
	
	tkl.y0 = Par[0];
	tkl.err_y0 = ParErr[0];
	
	tkl.ty = Par[1];
	tkl.err_ty = ParErr[1];
	
	Par[0] = tkl.x0;
	Par[1] = tkl.y0;
	Par[2] = tkl.tx;
	Par[3] = tkl.ty;

	chi2_straight(tkl.nXHits+tkl.nUHits+tkl.nVHits, fitarray->drift_dist+idxoff_array, fitarray->resolution+idxoff_array,
			fitarray->p1x+idxoff_array, fitarray->p1y+idxoff_array, fitarray->p1z+idxoff_array,
			fitarray->deltapx+idxoff_array, fitarray->deltapy+idxoff_array, fitarray->deltapz+idxoff_array,
			Par, chi2);
	
	tkl.chisq = chi2;
	return chi2;
}

#endif


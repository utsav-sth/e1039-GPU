#include "reconstruction_classes.cuh"


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

// ---------------------------------------------------------------- //
//
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
	float p1y = y_bep(detid, elid, planes);
	float p2x = x_tep(detid, elid, planes);
	
	float x_trk = x0+planes->z[ detid ]*tx;

	y = p1y + (x_trk-p1x) *  planes->deltapy[ detid ]/planes->deltapx[ detid ];
	
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf("det %d chan %d p1x %1.4f p1y %1.4f p2x %1.4f dpy %1.4f dpx% 1.4f x_trk %1.4f y %1.4f \n", detid, elid, p1x, p1y, p2x, planes->deltapy[ detid ], planes->deltapx[ detid ], x_trk, y );
#endif
	//if hitsign is zero, we don't want to toss a hit that could potentially be in range of the track accounting for the drift distance
	if(hitsign==0){
		if( x_trk-drift>p1x && x_trk-drift>p2x)return false;// if xtrk>p1x and >p2x, no overlap possible
		if( x_trk+drift<p1x && x_trk+drift<p2x)return false;// if xtrk<p1x and <p2x, no overlap possible

		y = max(y, p1y);
		y = min(y, p1y+planes->deltapy[ detid ]);
	}else{
		if( x_trk>p1x && x_trk>p2x)return false;// if xtrk>p1x and >p2x, no overlap possible
		if( x_trk<p1x && x_trk<p2x)return false;// if xtrk<p1x and <p2x, no overlap possible
	}

	err_y = planes->spacing[ detid ]/3.4641f * fabs(planes->deltapy[ detid ]/planes->deltapx[ detid ]);
	return true;
}


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
				if(detid > 36 || ((hitcoll.tdc(iAH) - tdcTime_curr >= 0.0) && (hitcoll.tdc(iAH) - tdcTime_curr < 80.0)) || ((hitcoll.tdc(iAH) - tdcTime_curr <= 0.0) && (hitcoll.tdc(iAH) - tdcTime_curr > -80.0))) {
					//if(ic[index].EventID==0)printf("hit det %d el %d Skip after-pulse...\n", detid, hitcoll.chan(iAH));
					hitflag[iAH] = -1;
					continue;
				}
				else {
					tdcTime_curr = hitcoll.tdc(iAH);
				}
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
						w_max = 0.9*0.5*(hitcoll.pos(cluster_iAH_arr[cluster_iAH_arr_size-1]) - hitcoll.pos(cluster_iAH_arr[0]));
						w_min = 4.0/9.0*w_max;
						
						if((hitcoll.drift(cluster_iAH_arr[0]) > w_max && hitcoll.drift(cluster_iAH_arr[cluster_iAH_arr_size-1]) > w_min) || (hitcoll.drift(cluster_iAH_arr[0]) > w_min && hitcoll.drift(cluster_iAH_arr[cluster_iAH_arr_size-1]) > w_max)) {
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
						}
						// if the time difference is less than 8 ns for detectors 19 to 24 (which btw are DC3p) we remove both
						else if((((hitcoll.tdc(cluster_iAH_arr[0]) - hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_size-1])) >= 0.0 && (hitcoll.tdc(cluster_iAH_arr[0]) - hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_size-1])) < 8.0) || ((hitcoll.tdc(cluster_iAH_arr[0]) - hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_size-1])) <= 0.0 && (hitcoll.tdc(cluster_iAH_arr[0]) - hitcoll.tdc(cluster_iAH_arr[cluster_iAH_arr_size-1])) > -8.0)) && (detid >= 19 && detid <= 24)) {
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
	
	for(iAH = 0; iAH<nhits; ++iAH) {
		if(hitflag[iAH]>0) {
			++nAH_reduced;
		}
	}
	return nAH_reduced;
}

// --------------------------------------------------------------- //
// Hit pairing functions
// --------------------------------------------------------------- //

__device__ void make_hitpairs_in_station_bins(const gHits hitcoll1, const int nhits1, const gHits hitcoll2, const int nhits2, thrust::pair<int, int>* hitpairs, int* npairs, const short bin0, const short Nbins, short* hitflag1, short* hitflag2, const int stID, const int projID){
	// I think we assume that by default we want to know where we are
	//printf("stID %d projID %d bin0 %d\n", stID, projID, bin0);
	
	short bin;
	const short MaxHits = geometry::MaxHitsProj[projID];

#ifdef DEBUG
#endif
	for(bin = bin0; bin<bin0+Nbins; bin++){
		npairs[bin-bin0] = 0;
	}
	
	//declaring arrays for the hit lists
	for(int i = 0; i<100; i++){
		hitflag1[i] = hitflag2[i] = 0;
	}
	
	//building the lists of hits for each detector plane
	//const int detid1 = geometry::detsuperid[stID][projID]*2;
	//const int detid2 = geometry::detsuperid[stID][projID]*2-1;
	const int superdetid = geometry::detsuperid[stID][projID];
	
	// pair the hits by position:
	// if one hit on e.g. x and one hit on x' are closer than
	// the "spacing" defined for the planes, then the hits can be paired together.
	int idx1 = -1;
	int idx2 = -1;
	
	//if(blockIdx.x==debug::EvRef && threadIdx.x==0)printf("nhits %d %d \n", nhits1, nhits2);

	for(int i = 0; i<nhits1; i++){
		idx1++;
		idx2 = -1;
		for(int j = 0; j<nhits2; j++){
			idx2++;
			//if(blockIdx.x==debug::EvRef && threadIdx==0)printf("i %d j %d pos %1.4f %1.4f\n", i, j, );

			if( abs(hitcoll1.pos(idx1) - hitcoll2.pos(idx2)) > geometry::spacingplane[superdetid] ){
				continue;
			}
			
			for(bin = bin0; bin<bin0+Nbins; bin++){
				if( geometry::WCHitsBins[stID-1][projID][0][bin] <= hitcoll1.chan(idx1) && 
				    hitcoll1.chan(idx1) <= geometry::WCHitsBins[stID-1][projID][1][bin]){
					//printf("bin %d low %d high %d hit 1 elem %d hit 2 elem %d global bin %d \n", bin, geometry::WCHitsBins[stID-1][projID][0][bin-bin0], geometry::WCHitsBins[stID-1][projID][1][bin-bin0], ic[index].AllHits[ i ].elementID, ic[index].AllHits[ idx2 ].elementID, bin+npairs[bin]*Nbins);
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
			//printf("bin %d low %d high %d hit elem %d global bin %d \n", bin, geometry::WCHitsBins[stID-1][projID][0][bin], geometry::WCHitsBins[stID-1][projID][1][bin], ic[index].AllHits[ i ].elementID, bin+npairs[bin]*Nbins);
				if( geometry::WCHitsBins[stID-1][projID][0][bin] <= hitcoll1.chan(i) && 
				    hitcoll1.chan(i) <= geometry::WCHitsBins[stID-1][projID][1][bin]){
					//printf("bin %d low %d high %d hit elem %d global bin %d \n", bin, geometry::WCHitsBins[stID-1][projID][0][bin-bin0], geometry::WCHitsBins[stID-1][projID][1][bin-bin0], ic[index].AllHits[ i ].elementID, bin+npairs[bin]*Nbins);
					if(npairs[bin-bin0]<=MaxHits)hitpairs[bin-bin0+npairs[bin-bin0]*Nbins] = thrust::make_pair(i, -1);
					npairs[bin-bin0]++;
				}
			}
		}
	}
	for(int i = 0; i<nhits2; i++){
		if(hitflag2[i]<1){
			for(bin = bin0; bin<bin0+Nbins; bin++){
			//printf("bin %d low %d high %d hit elem %d global bin %d \n", bin, geometry::WCHitsBins[stID-1][projID][0][bin], geometry::WCHitsBins[stID-1][projID][1][bin], ic[index].AllHits[ i ].elementID, bin+npairs[bin]*Nbins);
				if( geometry::WCHitsBins[stID-1][projID][0][bin] <= hitcoll2.chan(i) && 
				    hitcoll2.chan(i) <= geometry::WCHitsBins[stID-1][projID][1][bin]){
					//printf("bin %d low %d high %d hit elem %d global bin %d \n", bin, geometry::WCHitsBins[stID-1][projID][0][bin-bin0], geometry::WCHitsBins[stID-1][projID][1][bin-bin0], ic[index].AllHits[ i ].elementID, bin+npairs[bin]*Nbins);
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
	const int detid1 = geometry::detsuperid[stID][projID]*2;
	const int detid2 = geometry::detsuperid[stID][projID]*2-1;
	const int superdetid = geometry::detsuperid[stID][projID];
	float p1x, p2x;
	int hitctr1 = 0, hitctr2 = 0;
	for(int i = 0; i<nhits1; i++){
		//p1x = planes[ ic->AllHits[i].detectorID ].p1x_w1 + planes[ ic->AllHits[i].detectorID ].dp1x * (ic->AllHits[i].elementID-1);
		if(planes->deltapx[ detid1 ]>0){
			p1x = x_bep(detid1, hitcoll1.chan(i), planes);
			p2x = x_tep(detid1, hitcoll1.chan(i), planes);
		}else{
			p1x = x_tep(detid1, hitcoll1.chan(i), planes);
			p2x = x_bep(detid1, hitcoll1.chan(i), planes);
		}
		//printf("%d %d %1.6f p1-2x %1.6f %1.6f xmin-max %1.6f %1.6f \n", ic->AllHits[i].detectorID, ic->AllHits[i].elementID, ic->AllHits[i].pos, p1x, p2x, xmin, xmax);
		//if(xmin>-999.)printf("xmin %1.6f xmax %1.6f p1x %1.6f p2x %1.6f \n", xmin, xmax, p1x, p2x);
		if( (p1x <= xmax) && (p2x >= xmin) ){ 
			hitidx1[hitctr1] = i;
			hitctr1++;
		}
	}
	for(int i = 0; i<nhits2; i++){
		//p1x = planes[ ic->AllHits[i].detectorID ].p1x_w1 + planes[ ic->AllHits[i].detectorID ].dp1x * (ic->AllHits[i].elementID-1);
		if(planes[ detid2 ].deltapx>0){
			p1x = x_bep(detid2, hitcoll2.chan(i), planes);
			p2x = x_tep(detid2, hitcoll2.chan(i), planes);
		}else{
			p1x = x_tep(detid2, hitcoll2.chan(i), planes);
			p2x = x_bep(detid2, hitcoll2.chan(i), planes);
		}
		//printf("%d %d %1.6f p1-2x %1.6f %1.6f xmin-max %1.6f %1.6f \n", ic->AllHits[i].detectorID, ic->AllHits[i].elementID, ic->AllHits[i].pos, p1x, p2x, xmin, xmax);
		//if(xmin>-999.)printf("xmin %1.6f xmax %1.6f p1x %1.6f p2x %1.6f \n", xmin, xmax, p1x, p2x);
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
			if( abs( hitcoll1.pos(hitidx1[idx1]) - hitcoll2.pos(hitidx2[idx2]) ) > geometry::spacingplane[superdetid] ){
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
		idx = geometry::planetype[detid/2-1];
		detid_2 = geometry::detsuperid[2][idx]*2;
		
		pos_st3 = x_st3*costheta_[detid_2] + y_st3*sintheta_[detid_2];

		//printf(" i %d idx %d  pos %1.3f \n", i, idx, pos_st3);
		
		z_st1 = z_[detid];
		z_st2 = z_[detid_2];
		
		x_st2 = x0+tx*z_st2;
		y_st2 = y0+ty*z_st2;
		
		pos_st2 = x_st2*costheta_[detid_2] + y_st2*sintheta_[detid_2];
		
		//printf("det id %d z %1.3f s_det id %d z %1.3f x %1.3f y %1.3f \n", detid, z_st1, detid_2, z_st2, x_st2, y_st2);
		
        	s2_target = pos_st2 - pos_st3*(z_st2 - geometry::Z_TARGET)/(z_st3 - geometry::Z_TARGET);
        	s2_dump   = pos_st2 - pos_st3*(z_st2 - geometry::Z_DUMP)/(z_st3 - geometry::Z_DUMP);

		pos_exp_target = geometry::SAGITTA_TARGET_CENTER*s2_target + pos_st3*(z_st1 - geometry::Z_TARGET)/(z_st3 - geometry::Z_TARGET);
		pos_exp_dump   = geometry::SAGITTA_DUMP_CENTER*s2_dump + pos_st3*(z_st1 - geometry::Z_DUMP)/(z_st3 - geometry::Z_DUMP);
		win_target = fabs(s2_target*geometry::SAGITTA_TARGET_WIDTH);
		win_dump   = fabs(s2_dump*geometry::SAGITTA_DUMP_WIDTH);
		
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


__device__ void calculate_x0_tx_st1(const gTracklet tkl, float &x0, float &tx)
{	
	tx = tkl.tx + geometry::PT_KICK_KMAG * tkl.invP * tkl.charge;
	x0 = tkl.tx*geometry::Z_KMAG_BEND + tkl.x0 - tx * geometry::Z_KMAG_BEND;
}

__device__ void calculate_x0_tx_st1_with_errors(const gTracklet tkl, float &x0, float &tx, float &err_x0, float &err_tx)
{	
	tx = tkl.tx + geometry::PT_KICK_KMAG * tkl.invP * tkl.charge;
	x0 = tkl.tx*geometry::Z_KMAG_BEND + tkl.x0 - tx * geometry::Z_KMAG_BEND;
	
	err_tx = tkl.err_tx + fabs(tkl.err_invP*geometry::PT_KICK_KMAG);
        err_x0 = tkl.err_x0 + fabs(tkl.err_invP*geometry::PT_KICK_KMAG)*geometry::Z_KMAG_BEND;
}

__device__ float calculate_invP(float tx, float tx_st1, const short charge)
{
	return (tx_st1 - tx)*charge / geometry::PT_KICK_KMAG;
}

__device__ float calculate_invP_charge(float tx, float tx_st1, short& charge)
{
	float invP = (tx_st1 - tx) / geometry::PT_KICK_KMAG;
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
	return ( err_tx - err_tx )/ geometry::PT_KICK_KMAG;
} 



#ifdef OLDCODE


// ------------------------------------------------------------- //
// functions to evaluate the hit selection window for a 2D track //
// and to calculate y from u, v hits given x                     //
// ------------------------------------------------------------- //


__device__ void find_xmin_xmax_in_chamber(float &xmin, float &xmax, const gTrack2D track2d, const short stID, const short projID, const gPlane* planes)
{
	const int detid1 = geometry::detsuperid[stID][projID]*2;
	const int detid2 = geometry::detsuperid[stID][projID]*2-1;
	
	xmin = min(planes->z[detid1]*(track2d.tx_-track2d.err_tx_), planes->z[detid2]*(track2d.tx_-track2d.err_tx_));
	xmax = max(planes->z[detid1]*(track2d.tx_+track2d.err_tx_), planes->z[detid2]*(track2d.tx_+track2d.err_tx_));
	
	xmin = xmin + track2d.x_0-track2d.err_x_0-planes->spacing[detid1];
	xmax = xmax + track2d.x_0+track2d.err_x_0+planes->spacing[detid1];
}

__device__ void find_xmin_xmax_in_chamber(float &xmin, float &xmax, const gTracklet tkl, const short stID, const short projID, const gPlane* planes)
{
	const int detid1 = geometry::detsuperid[stID][projID]*2;
	const int detid2 = geometry::detsuperid[stID][projID]*2-1;
	
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


// --------------------------------------------- //
// functions to calculate x0 and tx in station 1 //
// and function to calculate inverse momentum    //
// --------------------------------------------- //



// ----------------------------------------------------------- //
// function to resolve the left right ambiguities in the track //
// ----------------------------------------------------------- //

__device__ void resolve_leftright(gTracklet &tkl, const gPlane* planes, const float thr)
{
	//bool isUpdated = false;
	short nhits = tkl.nXHits+tkl.nUHits+tkl.nVHits;
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
	for(short n = 0; n<nhits; n+=2){
		i = n;
		j = i+1;
		if( abs(tkl.hits[i].detectorID-tkl.hits[j].detectorID)!=1 ){
		    n--;//step back by 1 to move by 1 hit instead of 2
		    continue;		    
		}
		
		if(tkl.stationID>=6 && tkl.hits[i].detectorID<=6){
			calculate_x0_tx_st1_with_errors(tkl, x0, tx, err_x0, err_tx);
		}else{
			tx = tkl.tx;
			x0 = tkl.x0;
			err_x0 = tkl.err_x0;
			err_tx = tkl.err_tx;
		}
		
		slope_exp = planes->costheta[tkl.hits[i].detectorID]*tx + planes->sintheta[tkl.hits[i].detectorID]*tkl.ty;
		err_slope = fabs(planes->costheta[tkl.hits[i].detectorID]*err_tx) + fabs(planes->sintheta[tkl.hits[i].detectorID]*tkl.err_ty);
		
		inter_exp = planes->costheta[tkl.hits[i].detectorID]*x0 + planes->sintheta[tkl.hits[i].detectorID]*tkl.y0;
		err_inter = fabs(planes->costheta[tkl.hits[i].detectorID]*err_x0) + fabs(planes->sintheta[tkl.hits[i].detectorID]*tkl.err_y0);
		
#ifdef DEBUG
		printf("hits dets %d %d; exp slope %1.4f +- %1.4f inter %1.4f +- %1.4f \n", tkl.hits[i].detectorID, tkl.hits[j].detectorID, slope_exp, err_slope, inter_exp, err_inter);
		printf("hit 1 positions %1.4f, %1.4f hit 2 positions %1.4f, %1.4f \n", 
			position(tkl.hits[i], +1), position(tkl.hits[i], -1), 
			position(tkl.hits[j], +1), position(tkl.hits[j], -1));
#endif

		if(tkl.hitsign[i]*tkl.hitsign[j]==0){
			indexmin = -1;
			pull_min = 1.e6;
			for(int k = 0; k<4; k++){
				slope_local = ( position(tkl.hits[i], geometry::lrpossibility[k][0]) - position(tkl.hits[j], geometry::lrpossibility[k][1]) )/(planes->z[tkl.hits[i].detectorID]-planes->z[tkl.hits[j].detectorID]);
				inter_local = position(tkl.hits[i], geometry::lrpossibility[k][0]) - slope_local*planes->z[tkl.hits[i].detectorID];
				
				if(fabs(slope_local) > planes->slope_max[tkl.hits[i].detectorID] || fabs(inter_local) > planes->inter_max[tkl.hits[i].detectorID])continue;
				
				pull = sqrtf( (slope_exp-slope_local)*(slope_exp-slope_local)/err_slope/err_slope + (inter_exp-inter_local)*(inter_exp-inter_local)/err_inter/err_inter );
				
#ifdef DEBUG
				printf("lr %d %d, slope %1.4f inter %1.4f\n", geometry::lrpossibility[k][0], geometry::lrpossibility[k][1], slope_local, inter_local);
				printf("pull %1.4f\n", pull);
#endif		
				if(pull<pull_min){
					indexmin = k;
					pull_min = pull;
				}
			}
			
			if(indexmin>0 && pull_min<thr){
				tkl.hitsign[i] = geometry::lrpossibility[indexmin][0];
				tkl.hitsign[j] = geometry::lrpossibility[indexmin][1];
				//isUpdated = true;
			}
		}
	//	++nresolved;
	}
}

__device__ void resolve_single_leftright(gTracklet &tkl, const gPlane* planes)
{
	short nhits = tkl.nXHits+tkl.nUHits+tkl.nVHits;
	float pos_exp;
	short detID;
	float x0, tx;// x0 and tx are different for global track station 1 hits 
	
	for(short n = 0; n<nhits; n++){

		if(tkl.stationID>=6 && tkl.hits[n].detectorID<=6){
			calculate_x0_tx_st1(tkl, x0, tx);
		}else{
			tx = tkl.tx;
			x0 = tkl.x0;
		}
		
		// don't do anything for hits whichs already have a sign...
		if(tkl.hitsign[n]!=0)continue;
		
		detID = tkl.hits[n].detectorID;
		pos_exp = (planes->z[detID]*tx+x0)*planes->costheta[detID]+(planes->z[detID]*tkl.ty+tkl.y0)*planes->sintheta[detID];
		tkl.hitsign[n] = pos_exp>tkl.hits[n].pos? +1 : -1;
	}
	
}


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

	if(z<geometry::Z_KMAG_BEND){
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
        if(z<geometry::Z_KMAG_BEND){
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


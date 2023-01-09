#include "reconstruction_classes.cuh"
//#include "tracknumericalminimizer.cuh"
#include "trackanalyticalminimizer.cuh"

// kernel functions: 
// CUDA C++ extends C++ by allowing the programmer to define C++ functions, called kernels, that, when called, 
// are executed N times in parallel by N different CUDA threads, as opposed to only once like regular C++ functions. 

// event reducer: 
__global__ void gkernel_eR(gEvent* ic) {
	//printf("Running the kernel function...\n");
	// retrieve global thread index
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	double w_max[EstnEvtMax]; // max drift distance of the hit furthest from the cluster avg position // current average position of cluster * 0.9
	
	double w_min[EstnEvtMax]; // max drift distance of the hit closest to the cluster avg position // current average position of cluster * 0.4
	double dt_mean[EstnEvtMax]; // tbd
	int cluster_iAH_arr_cur[EstnEvtMax]; // current cluster array
	int cluster_iAH_arr_size[EstnEvtMax]; // cluster size i.e number of hits in cluster
	static int cluster_iAH_arr[EstnEvtMax][ClusterSizeMax]; // global cluster array 
	int uniqueID[EstnEvtMax]; // hit unique element ID
	int uniqueID_curr[EstnEvtMax]; // current hit unique element ID
	double tdcTime_curr[EstnEvtMax]; // current hit TDC time
	int iAH[EstnEvtMax]; // hit index 
	int nAH_reduced[EstnEvtMax]; // number of hits after hit quality filtering
	// int nHitsPerDetector[nDetectors+1];
	

	//if(ic[index].EventID==0){
	//	printf("evt = %d, nAH = %d: \n", ic[index].EventID, ic[index].nAH);
	//	for(int i = 1; i<=nDetectors; i++)printf(" det %d: %d;  ", i, ic[index].NHits[i]);
	//	printf("\n");
	//}

	// initialization of array size
	cluster_iAH_arr_size[index] = 0;
	nAH_reduced[index] = 0;
	
	// event reducing/hit filtering
	for(iAH[index] = 0; iAH[index]<ic[index].nAH; ++iAH[index]) {
		// if hit not good, set its detID to 0 and continue;
		if((ic[index].AllHits[iAH[index]].flag & hitFlagBit(1)) == 0) {
			//if(ic[index].EventID==0)printf("hit det %d Skip out-of-time...\n", ic[index].AllHits[iAH[index]].detectorID);
			ic[index].AllHits[iAH[index]].detectorID = 0;
			continue;
		}
		uniqueID[index] = uniqueID_curr[index] = -1;

		// hits in DCs or Prop tubes
		if(ic[index].AllHits[iAH[index]].detectorID < 31 || ic[index].AllHits[iAH[index]].detectorID > 46) {
			// evaluate "unique ID"
			uniqueID[index] = ic[index].AllHits[iAH[index]].detectorID*1000 + ic[index].AllHits[iAH[index]].elementID;
			// compare with current unique element ID; if different, update the unique element ID and time info 
			if(uniqueID[index] != uniqueID_curr[index]) {
				uniqueID_curr[index] = uniqueID[index];
				tdcTime_curr[index] = ic[index].AllHits[iAH[index]].tdcTime;
			}
			// if next hit and current hit belong to the same element: if detID>36 => prop tubes (reminder that hodoscpes are out of the picture in this scope), 
			// we're suppose to have one signal (a second signal would be after-pulsing)
			// if time difference between new hit and current hit is less than 80ns for DCs, it's also considered after-pulsing
			else {
				if(ic[index].AllHits[iAH[index]].detectorID > 36 || ((ic[index].AllHits[iAH[index]].tdcTime - tdcTime_curr[index] >= 0.0) && (ic[index].AllHits[iAH[index]].tdcTime - tdcTime_curr[index] < 80.0)) || ((ic[index].AllHits[iAH[index]].tdcTime - tdcTime_curr[index] <= 0.0) && (ic[index].AllHits[iAH[index]].tdcTime - tdcTime_curr[index] > -80.0))) {
					//if(ic[index].EventID==0)printf("hit det %d el %d Skip after-pulse...\n", ic[index].AllHits[iAH[index]].detectorID, ic[index].AllHits[iAH[index]].elementID);
					ic[index].AllHits[iAH[index]].detectorID = 0;
					continue;
				}
				else {
					tdcTime_curr[index] = ic[index].AllHits[iAH[index]].tdcTime;
				}
			}
		}
		// declustering of hits in DCs (from CPU code, I understand this one better)
		// if there are hits in the same plane and hitting to neighboring wires, they both give redundant information: 
		if(ic[index].AllHits[iAH[index]].detectorID <= nChamberPlanes) {
			//if(ic[index].EventID==0)printf("%d\n", cluster_iAH_arr_size[index]);
//			printf("Decluster...\n");
			if(cluster_iAH_arr_size[index] == ClusterSizeMax) {
//				printf("Oversized cluster...\n");
			}
			// if array size is zero, start storing the hit in the array
			if(cluster_iAH_arr_size[index] == 0) {
				cluster_iAH_arr[index][0] = iAH[index];
				++cluster_iAH_arr_size[index];
			} else { // otherwise
				// current hit and previous hit are *not* in same detector plane OR next hit and current hit are *not* in neighbors cells
				// we "declusterize" i.e. we remove the hit/hits which information is redundant with other hits and/or useless
				//if(ic[index].EventID==0){
//printf("hit indices: %d %d %d\n", iAH[index], cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1], cluster_iAH_arr[index][0]);
//printf("hit det/elem: %d, %d; %d, %d\n", ic[index].AllHits[iAH[index]].detectorID, ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID, 
//	    	      	      	  	 ic[index].AllHits[iAH[index]].elementID, ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID);
//printf("diffs: %d, %d\n", (ic[index].AllHits[iAH[index]].detectorID - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID),
//	       	   	  (ic[index].AllHits[iAH[index]].elementID - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID));
//printf("bools: %d, %d\n", (ic[index].AllHits[iAH[index]].detectorID != ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID),
//	       	   	  (abs(ic[index].AllHits[iAH[index]].elementID - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID) > 1));
	    	     	   	//}
				if((ic[index].AllHits[iAH[index]].detectorID != ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID) || (abs(ic[index].AllHits[iAH[index]].elementID - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID) > 1)) {
					// if 2 hits in cluster, evaluate w_max and w_min; drift distance has to be < w_min for one of the hits, while it has to be < w_max for the other hit 
					if(cluster_iAH_arr_size[index] == 2) {
						w_max[index] = 0.9*0.5*(ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].pos - ic[index].AllHits[cluster_iAH_arr[index][0]].pos);
						w_min[index] = 4.0/9.0*w_max[index];
						if((ic[index].AllHits[cluster_iAH_arr[index][0]].driftDistance > w_max[index] && ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].driftDistance > w_min[index]) || (ic[index].AllHits[cluster_iAH_arr[index][0]].driftDistance > w_min[index] && ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].driftDistance > w_max[index])) {
						//if(ic[index].EventID==0)printf("hit indices: %d %d %d\n", iAH[index], cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1], cluster_iAH_arr[index][0]);
							//eliminating the existing hit with the lagest drift distance
							if(ic[index].AllHits[cluster_iAH_arr[index][0]].driftDistance > ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].driftDistance) {
								//if(ic[index].EventID==0)printf("1 - hit det %d elem %d Skip cluster...\n", ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID, ic[index].AllHits[cluster_iAH_arr[index][0]].elementID);
								ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID = 0;
							}
							else {
								//if(ic[index].EventID==0)printf("2 - hit det %d elem %d Skip cluster...\n", ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID, ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID);
								ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID = 0;
							}
						}
						// if the time difference is less than 8 ns for detectors 19 to 24 (which btw are DC3p) we remove both
						else if((((ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) >= 0.0 && (ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) < 8.0) || ((ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) <= 0.0 && (ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) > -8.0)) && (ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID >= 19 && ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID <= 24)) {
						        //if(ic[index].EventID==0)printf("3 - hit det %d elem %d Skip cluster...\n", ic[index].AllHits[iAH[index]].detectorID, ic[index].AllHits[iAH[index]].elementID);
							ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID = 0;
							ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID = 0;
						}
					}
					// if 3 hits or more in cluster: we essentially discard them all;
					if(cluster_iAH_arr_size[index] >= 3) {
						// evaluate the mean time difference;
						dt_mean[index] = 0.0;
						for(cluster_iAH_arr_cur[index] = 1; cluster_iAH_arr_cur[index] < cluster_iAH_arr_size[index]; ++cluster_iAH_arr_cur[index]) {
							dt_mean[index] += ((ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]-1]].tdcTime) > 0.0 ? (ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]-1]].tdcTime) : (ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]-1]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].tdcTime));
						}
						dt_mean[index] = dt_mean[index]/(cluster_iAH_arr_size[index] - 1);
						// if mean time difference is less than 10, that's electronic noise, so we remove them all.
						if(dt_mean[index] < 10.0) {
						        //if(ic[index].EventID==0)printf("4 - hit det %d elem %d Skip cluster...\n", ic[index].AllHits[iAH[index]].detectorID, ic[index].AllHits[iAH[index]].elementID);
							for(cluster_iAH_arr_cur[index] = 0; cluster_iAH_arr_cur[index] < cluster_iAH_arr_size[index]; ++cluster_iAH_arr_cur[index]) {
								ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].detectorID = 0;
							}
						}
						// otherwise, we remove them all except first and last
						else {
						        //if(ic[index].EventID==0)printf("5 - hit det %d elem %d Skip cluster...\n", ic[index].AllHits[iAH[index]].detectorID, ic[index].AllHits[iAH[index]].elementID);
							for(cluster_iAH_arr_cur[index] = 1; cluster_iAH_arr_cur[index] < cluster_iAH_arr_size[index]; ++cluster_iAH_arr_cur[index]) {
								ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].detectorID = 0;
							}
						}
					}
					cluster_iAH_arr_size[index] = 0;
				} else { // if hits are in the same detector and in two neighboring cells!
				        // current hit and previous hit are in same detector plane and in neighbor wires: 
				  	// we count how many hits we have in this case, until we find a hit in a different detector or in a wire that is not a neighbor to the previous hit.
				  	//if(ic[index].EventID==0)printf("h1: det %d el %d; h2 det %d el %d \n", ic[index].AllHits[iAH[index]].detectorID, ic[index].AllHits[iAH[index]].elementID, ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID, ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID);
				  	cluster_iAH_arr[index][cluster_iAH_arr_size[index]] = iAH[index];
				  	++cluster_iAH_arr_size[index];
				}
			}
		}
	}
	//end of the hit loop

	// Hit reduction: 
	// store in "AllHits" containers only hits with non-zero detectorID and couting those with nAH_reduced
	for(iAH[index] = 0; iAH[index]<ic[index].nAH; ++iAH[index]) {
		if(ic[index].AllHits[iAH[index]].detectorID != 0) {
			ic[index].AllHits[nAH_reduced[index]] = ic[index].AllHits[iAH[index]];
			++nAH_reduced[index];

						
		}
	}

	// compute hits per detector
	int nEventHits = nAH_reduced[index];
	// reinitialize number of hits per detector
	for(auto iDetector = 1; iDetector <= nDetectors; ++iDetector) {
		ic[index].NHits[iDetector] = 0;
	}
	// loop on reduced hits and counting number of hits per detector
	for(auto iHit = 0; iHit < nEventHits; ++iHit) {
		auto detectorId = ic[index].AllHits[iHit].detectorID;
		if(detectorId != 0) {
			++ic[index].NHits[detectorId];
		}
	}


	ic[index].nAH = nAH_reduced[index];
}

// function to match a tracklet to a hodoscope hit
__device__ int match_tracklet_to_hodo(const gTracklet tkl, const int stID, gEvent* ic, const gPlane* planes)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int masked = 0;//false
	// first, define the search region, and foremost, the planes on which we define this search region, which depends on the station ID we're looking at
	// define the region in which we are supposed to have hits:

	//printf(" stID %d hodo plane[0] %d [1] %d \n", stID, geometry::hodoplanerange[stID][0], geometry::hodoplanerange[stID][1]);
	//printf(" x0 %1.4f +- %1.4f, y0 %1.4f +- %1.4f, tx %1.4f +- %1.4f, ty %1.4f +- %1.4f \n", tkl.x0, tkl.err_x0, tkl.y0, tkl.err_y0, tkl.tx, tkl.err_tx, tkl.ty, tkl.err_ty);		
	
	REAL xhodo, yhodo, err_x, err_y, xmin, xmax, ymin, ymax;
		
	// loop on the hits and select hodoscope hits corresponding to the station
	for(int i = 0; i<ic[index].nAH; i++){
		//we only consider hits in the hodoscopes planes corresponding to the station where the tracklet is reconstructed 
		if(geometry::hodoplanerange[stID][0]>ic[index].AllHits[i].detectorID || geometry::hodoplanerange[stID][1]<ic[index].AllHits[i].detectorID)continue;
				
		//calculate the track position at the hodoscope plane z
		xhodo = planes[ic[index].AllHits[i].detectorID].z*tkl.tx+tkl.x0;
		yhodo = planes[ic[index].AllHits[i].detectorID].z*tkl.ty+tkl.y0;
		
		err_x = 3.f*(fabs(planes[ic[index].AllHits[i].detectorID].z*tkl.err_tx)+tkl.err_x0);
		err_y = 3.f*(fabs(planes[ic[index].AllHits[i].detectorID].z*tkl.err_ty)+tkl.err_y0);

		//printf(" det %d elem %d z_hodo %1.4f x_hodo %1.4f y_hodo %1.4f err x %1.4f err y %1.4f \n", ic[index].AllHits[i].detectorID, ic[index].AllHits[i].elementID, planes[ic[index].AllHits[i].detectorID].z, xhodo, yhodo, err_x, err_y);

		//calculate "xmin, xmax, ymin, ymax" in which the track is supposed to pass through; 
		//these are basicially defined as the spatial coverage of the hit hodoscope element (plus some fudge factor for x)
		if(planes[ic[index].AllHits[i].detectorID].costheta>0.99){
			xmin = ic[index].AllHits[i].pos-planes[ic[index].AllHits[i].detectorID].cellwidth*0.5f;
			xmax = ic[index].AllHits[i].pos+planes[ic[index].AllHits[i].detectorID].cellwidth*0.5f;
			
			ymin = planes[ic[index].AllHits[i].detectorID].y1;
			ymax = planes[ic[index].AllHits[i].detectorID].y2;
			
			xmin-=(xmax-xmin)*geometry::hodofudgefac[stID];
			xmax+=(xmax-xmin)*geometry::hodofudgefac[stID];
			
			ymin-=(ymax-ymin)*geometry::hodofudgefac[stID];
			ymax+=(ymax-ymin)*geometry::hodofudgefac[stID];
			
			ymin+=planes[ic[index].AllHits[i].detectorID].y0;
			ymax+=planes[ic[index].AllHits[i].detectorID].y0;
		}else{
			xmin = planes[ic[index].AllHits[i].detectorID].x1;
			xmax = planes[ic[index].AllHits[i].detectorID].x2;
			
			ymin = ic[index].AllHits[i].pos-planes[ic[index].AllHits[i].detectorID].cellwidth*0.5f;
			ymax = ic[index].AllHits[i].pos+planes[ic[index].AllHits[i].detectorID].cellwidth*0.5f;
			
			xmin-=(xmax-xmin)*geometry::hodofudgefac[stID];
			xmax+=(xmax-xmin)*geometry::hodofudgefac[stID];
			
			ymin-=(ymax-ymin)*geometry::hodofudgefac[stID];
			ymax+=(ymax-ymin)*geometry::hodofudgefac[stID];
		}

		//printf(" xmin %1.4f xmax %1.4f, ymin %1.4f ymax %1.4f xfudge %1.4f\n", xmin, xmax, ymin, ymax, (xmax-xmin)*0.15 );
		err_x+= (xmax-xmin)*0.15;

		xmin-= err_x;
		xmax+= err_x;
		ymin-= err_y;
		ymax+= err_y;
		
		//printf(" xmin %1.4f xmax %1.4f, ymin %1.4f ymax %1.4f \n", xmin, xmax, ymin, ymax);
		if(xmin <= xhodo && xhodo <= xmax && ymin <= yhodo && yhodo <= ymax ){
			masked++;
			break;
		}
		
	}
	// 
	return masked>0;
}




// function to make the hit pairs in station;
// I assume it will only be called by the tracklet builder
// (not by the main function), so I can make it a "device" function. 
__device__ int make_hitpairs_in_station(gEvent* ic, thrust::pair<int, int>* hitpairs, int* hitidx1, int* hitidx2, short* hitflag1, short* hitflag2, const int stID, const int projID){
	// I think we assume that by default we want to know where we are
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	int npairs = 0;
	
	//declaring arrays for the hit lists
	for(int i = 0; i<100; i++){
		hitidx1[i] = hitidx2[i] = 0;
		hitflag1[i] = hitflag2[i] = 0;
	}
	
	//building the lists of hits for each detector plane
	int detid1 = geometry::detsuperid[stID][projID]*2;
	int detid2 = geometry::detsuperid[stID][projID]*2-1;
	int superdetid = geometry::detsuperid[stID][projID];
	int hitctr1 = 0, hitctr2 = 0;
	for(int i = 0; i<ic[index].nAH; i++){
		if(ic[index].AllHits[i].detectorID==detid1){
			hitidx1[hitctr1] = i;
			hitctr1++;
		}
		if(ic[index].AllHits[i].detectorID==detid2){
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
			if( abs(ic[index].AllHits[ hitidx1[idx1] ].pos - ic[index].AllHits[ hitidx2[idx2] ].pos) > geometry::spacingplane[superdetid] ){
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
			hitpairs[npairs] = thrust::make_pair(hitidx2[i], -1);
			npairs++;
		   }
	}
	   	   
	return npairs;
}



__device__ int make_hitpairs_in_station(gEvent* ic, thrust::pair<int, int>* hitpairs, int* hitidx1, int* hitidx2, short* hitflag1, short* hitflag2, const int stID, const int projID, const gPlane* planes, const REAL xmin, const REAL xmax){
	// I think we assume that by default we want to know where we are
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	int npairs = 0;
	
	//declaring arrays for the hit lists
	for(int i = 0; i<100; i++){
		hitidx1[i] = hitidx2[i] = 0;
		hitflag1[i] = hitflag2[i] = 0;
	}
	
	//building the lists of hits for each detector plane
	int detid1 = geometry::detsuperid[stID][projID]*2;
	int detid2 = geometry::detsuperid[stID][projID]*2-1;
	int superdetid = geometry::detsuperid[stID][projID];
	int hitctr1 = 0, hitctr2 = 0;
	float p1x, p2x;
	for(int i = 0; i<ic[index].nAH; i++){
		if(ic[index].AllHits[i].detectorID==detid1){
			p1x = planes[ ic[index].AllHits[i].detectorID ].p1x_w1 + planes[ ic[index].AllHits[i].detectorID ].dp1x * (ic[index].AllHits[i].elementID-1);
			if(planes[ ic[index].AllHits[i].detectorID ].deltapx>0){
				p2x = p1x + planes[ ic[index].AllHits[i].detectorID ].deltapx;
			}else{
				p2x = p1x;
				p1x+= planes[ ic[index].AllHits[i].detectorID ].deltapx;
			}
			//printf("%d %d %1.6f p1-2x %1.6f %1.6f xmin-max %1.6f %1.6f \n", ic[index].AllHits[i].detectorID, ic[index].AllHits[i].elementID, ic[index].AllHits[i].pos, p1x, p2x, xmin, xmax);
			//if(xmin>-999.)printf("xmin %1.6f xmax %1.6f p1x %1.6f p2x %1.6f \n", xmin, xmax, p1x, p2x);
			if( (p1x <= xmax) && (p2x >= xmin) ){ 
				hitidx1[hitctr1] = i;
				hitctr1++;
			}
		}
		if(ic[index].AllHits[i].detectorID==detid2){
			p1x = planes[ ic[index].AllHits[i].detectorID ].p1x_w1 + planes[ ic[index].AllHits[i].detectorID ].dp1x * (ic[index].AllHits[i].elementID-1);
			if(planes[ ic[index].AllHits[i].detectorID ].deltapx>0){
				p2x = p1x + planes[ ic[index].AllHits[i].detectorID ].deltapx;
			}else{
				p2x = p1x;
				p1x+= planes[ ic[index].AllHits[i].detectorID ].deltapx;
			}
			//printf("%d %d %1.6f p1-2x %1.6f %1.6f xmin-max %1.6f %1.6f \n", ic[index].AllHits[i].detectorID, ic[index].AllHits[i].elementID, ic[index].AllHits[i].pos, p1x, p2x, xmin, xmax);
			//if(xmin>-999.)printf("xmin %1.6f xmax %1.6f p1x %1.6f p2x %1.6f \n", xmin, xmax, p1x, p2x);
			if( (p1x <= xmax) && (p2x >= xmin) ){ 
				hitidx2[hitctr2] = i;
				hitctr2++;
			}
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
			if( abs(ic[index].AllHits[ hitidx1[idx1] ].pos - ic[index].AllHits[ hitidx2[idx2] ].pos) > geometry::spacingplane[superdetid] ){
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
			hitpairs[npairs] = thrust::make_pair(hitidx2[i], -1);
			npairs++;
		   }
	}
	   	   
	return npairs;
}




__device__ void calculate_x0_tx_st1(const gTracklet tkl, float &x0, float &tx, float &err_x0, float &err_tx)
{	
	tx = tkl.tx + geometry::PT_KICK_KMAG * tkl.invP * tkl.charge;
	x0 = tkl.tx*geometry::Z_KMAG_BEND + tkl.x0 - tx * geometry::Z_KMAG_BEND;
	
	err_tx = tkl.err_tx + fabs(tkl.err_invP*geometry::PT_KICK_KMAG);
        err_x0 = tkl.err_x0 + fabs(tkl.err_invP*geometry::PT_KICK_KMAG)*geometry::Z_KMAG_BEND;
}


__device__ void resolve_leftright(gTracklet &tkl, const gPlane* planes, const float thr)
{
	//bool isUpdated = false;
	short nhits = tkl.nXHits+tkl.nUHits+tkl.nVHits;
	short nresolved = 0;
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
		
		if(tkl.hitsign[i]*tkl.hitsign[j]==0){
			indexmin = -1;
			pull_min = 1.e6;
			for(int k = 0; k<4; k++){
				slope_local = ( tkl.hits[i].pos+tkl.hits[i].driftDistance*geometry::lrpossibility[k][0] - tkl.hits[j].pos+tkl.hits[j].driftDistance*geometry::lrpossibility[k][1] )/(planes[tkl.hits[i].detectorID].z-planes[tkl.hits[j].detectorID].z);
				inter_local = tkl.hits[i].pos+tkl.hits[i].driftDistance*geometry::lrpossibility[k][0] - slope_local*planes[tkl.hits[i].detectorID].z;
				
				if(fabs(slope_local) > planes[tkl.hits[i].detectorID].slope_max || fabs(inter_local) > planes[tkl.hits[i].detectorID].inter_max)continue;
				
				if(tkl.stationID>=6 && tkl.hits[i].detectorID<=6){
					calculate_x0_tx_st1(tkl, x0, tx, err_x0, err_tx);
				}else{
					tx = tkl.tx;
					x0 = tkl.x0;
					err_x0 = tkl.err_x0;
					err_tx = tkl.err_tx;
				}
				
				slope_exp = planes[tkl.hits[i].detectorID].costheta*tx + planes[tkl.hits[i].detectorID].sintheta*tkl.ty;
				err_slope = fabs(planes[tkl.hits[i].detectorID].costheta*err_tx) + fabs(planes[tkl.hits[i].detectorID].sintheta*tkl.err_ty);
				
				inter_exp = planes[tkl.hits[i].detectorID].costheta*x0 + planes[tkl.hits[i].detectorID].sintheta*tkl.y0;
				err_inter = fabs(planes[tkl.hits[i].detectorID].costheta*err_x0) + fabs(planes[tkl.hits[i].detectorID].sintheta*tkl.err_y0);
				
				pull = sqrtf( (slope_exp-slope_local)*(slope_exp-slope_local)/err_slope/err_slope + (inter_exp-inter_local)*(inter_exp-inter_local)/err_inter/err_inter );
				
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
		++nresolved;
	}
}




__device__ void find_xmin_xmax_in_chamber(float &xmin, float &xmax, const gTrackXZ trackxz, const short stID, const short projID, const gPlane* planes)
{
	int detid1 = geometry::detsuperid[stID][projID]*2;
	int detid2 = geometry::detsuperid[stID][projID]*2-1;
	
	xmin = min(planes[detid1].z*(trackxz.tx-trackxz.err_tx), planes[detid2].z*(trackxz.tx-trackxz.err_tx));
	xmax = max(planes[detid1].z*(trackxz.tx+trackxz.err_tx), planes[detid2].z*(trackxz.tx+trackxz.err_tx));
	
	xmin = xmin + trackxz.x0-trackxz.err_x0-planes[detid1].spacing;
	xmax = xmax + trackxz.x0+trackxz.err_x0+planes[detid1].spacing;
}

__device__ void FillFitArrays(const int n, const gHit hit, const short hitsign, gStraightFitArrays &fitarray, const gPlane* planes){
	fitarray.drift_dist[n] = hit.driftDistance*hitsign;
	fitarray.resolution[n] = planes[ hit.detectorID ].resolution;
	if(hitsign==0)fitarray.resolution[n] = planes[ hit.detectorID ].spacing*3.4641f;
	
	fitarray.p1x[n] = planes[ hit.detectorID ].p1x_w1 + planes[ hit.detectorID ].dp1x * (hit.elementID-1);
	fitarray.p1y[n] = planes[ hit.detectorID ].p1y_w1 + planes[ hit.detectorID ].dp1y * (hit.elementID-1);
	fitarray.p1z[n] = planes[ hit.detectorID ].p1z_w1 + planes[ hit.detectorID ].dp1z * (hit.elementID-1);
	
	fitarray.deltapx[n] = planes[ hit.detectorID ].deltapx;
	fitarray.deltapy[n] = planes[ hit.detectorID ].deltapy;
	fitarray.deltapz[n] = planes[ hit.detectorID ].deltapz;
}

__device__ bool calculate_y_uvhit(float &y, float &err_y, const gHit hit, const short hitsign, const gTrackXZ trackxz, const gPlane* planes){
	float p1x = planes[ hit.detectorID ].p1x_w1 + planes[ hit.detectorID ].dp1x * (hit.elementID-1);
	float p1y = planes[ hit.detectorID ].p1y_w1 + planes[ hit.detectorID ].dp1y * (hit.elementID-1);
	float p2x = p1x+planes[ hit.detectorID ].deltapx;
	
	float x_trk = trackxz.x0+planes[ hit.detectorID ].z*trackxz.tx - hit.driftDistance*hitsign;

	y = p1y + (x_trk-hit.driftDistance*hitsign-p1x) *  planes[ hit.detectorID ].deltapy/planes[ hit.detectorID ].deltapx;
	
	//if hitsign is zero, we don't want to toss a hit that could potentially be in range of the track accounting for the drift distance
	if(hitsign==0){
		if( x_trk-hit.driftDistance>p1x && x_trk-hit.driftDistance>p2x)return false;// if xtrk>p1x and >p2x, no overlap possible
		if( x_trk+hit.driftDistance<p1x && x_trk+hit.driftDistance<p2x)return false;// if xtrk<p1x and <p2x, no overlap possible

		y = max(y, p1y);
		y = min(y, p1y+planes[ hit.detectorID ].deltapy);
	}else{
		if( x_trk>p1x && x_trk>p2x)return false;// if xtrk>p1x and >p2x, no overlap possible
		if( x_trk<p1x && x_trk<p2x)return false;// if xtrk<p1x and <p2x, no overlap possible
	}

	err_y = planes[ hit.detectorID ].spacing/3.4641f * fabs(planes[ hit.detectorID ].deltapy/planes[ hit.detectorID ].deltapx);
	return true;
}



__global__ void gKernel_XZ_YZ_tracking(gEvent* ic, gOutputEvent* oc, gStraightTrackBuilder* straighttrackbuilder, gStraightFitArrays* fitarrays, const gPlane* planes)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	straighttrackbuilder[index].nTracksXZ = 0;
	straighttrackbuilder[index].nTracksYZ = 0;
	
	oc[index].EventID = ic[index].EventID;
	oc[index].nAH = ic[index].nAH;
	oc[index].nTracklets = 0;
	
	//for the time being, declare in function;
	//REAL x_pos[4];
	//REAL x_res[4];
	//REAL z_pos[4];
	//REAL a, b, sum, det, sx, sy, sxx, syy, sxy;
	short nprop, iprop;
	REAL xExp;
	short nhits_X2, nhits_X3;
	short nhits_U2, nhits_U3;
	short nhits_V2, nhits_V3;
	
	short stid, projid;
	
        // 1- get X pairs in st2, st3:
	projid = 0;
	//D2: stid = 3-1
	stid = 2;
	int nx2 = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_x2, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid);
	
	//D3p: stid = 4-1
	stid = 3;
	int nx3p = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_x3p, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid);

	//D3p: stid = 5-1
	stid = 4;
	int nx3m = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_x3m, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid);

        // 2- loop on X hit pairs; calculate slope between the hit X pairs (i.e. XZ tracking):
	for(int i = 0; i<nx2; i++){
		for(int n = 0; n<nChamberPlanes; n++){
			fitarrays[index].x_array[n] = 0;
			fitarrays[index].z_array[n] = 0;
			fitarrays[index].dx_array[n] = 0;
		}
		
		nhits_X2 = 0;
		if(straighttrackbuilder[index].hitpairs_x2[i].first>=0){
			fitarrays[index].x_array[nhits_X2] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x2[i].first].pos;
			fitarrays[index].dx_array[nhits_X2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x2[i].first].detectorID].spacing/3.4641f;
			fitarrays[index].z_array[nhits_X2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x2[i].first].detectorID].z;
			straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].hitlist[nhits_X2] = straighttrackbuilder[index].hitpairs_x2[i].first;
			nhits_X2++;
		}
		if(straighttrackbuilder[index].hitpairs_x2[i].second>=0){
			fitarrays[index].x_array[nhits_X2] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x2[i].second].pos;
			fitarrays[index].dx_array[nhits_X2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x2[i].second].detectorID].spacing/3.4641f;
			fitarrays[index].z_array[nhits_X2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x2[i].second].detectorID].z;
			straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].hitlist[nhits_X2] = straighttrackbuilder[index].hitpairs_x2[i].second;
			nhits_X2++;
		}
			
		for(int j = 0; j<nx3p+nx3m; j++){
			nhits_X3 = nhits_X2;
			if(j<nx3p){
				if(straighttrackbuilder[index].hitpairs_x3p[j].first>=0){
					fitarrays[index].x_array[nhits_X3] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3p[j].first].pos;
					fitarrays[index].dx_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3p[j].first].detectorID].spacing/3.4641f;
					fitarrays[index].z_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3p[j].first].detectorID].z;
					straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].hitlist[nhits_X3] = straighttrackbuilder[index].hitpairs_x3p[j].first;
					nhits_X3++;
				}
				if(straighttrackbuilder[index].hitpairs_x3p[j].second>=0){
					fitarrays[index].x_array[nhits_X3] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3p[j].second].pos;
					fitarrays[index].dx_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3p[j].second].detectorID].spacing/3.4641f;
					fitarrays[index].z_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3p[j].second].detectorID].z;
					straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].hitlist[nhits_X3] = straighttrackbuilder[index].hitpairs_x3p[j].second;
					nhits_X3++;
				}
			}else{
				if(straighttrackbuilder[index].hitpairs_x3m[j].first>=0){
					fitarrays[index].x_array[nhits_X3] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[j].first].pos;
					fitarrays[index].dx_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[j].first].detectorID].spacing/3.4641f;
					fitarrays[index].z_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[j].first].detectorID].z;
					straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].hitlist[nhits_X3] = straighttrackbuilder[index].hitpairs_x3m[j].first;
					nhits_X3++;
				}
				if(straighttrackbuilder[index].hitpairs_x3m[j].second>=0){
					fitarrays[index].x_array[nhits_X3] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[j].second].pos;
					fitarrays[index].dx_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[j].second].detectorID].spacing/3.4641f;
					fitarrays[index].z_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[j].second].detectorID].z;
					straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].hitlist[nhits_X3] = straighttrackbuilder[index].hitpairs_x3m[j].second;
					nhits_X3++;
				}
			}
			//chi2_simplefit(nhits_X3, fitarrays[index].z_array, fitarrays[index].x_array, a, b, sum, det, sx, sy, sxx, syy, sxy);
			//if(fabs(a)>2*TX_MAX || fabs(b)>2*X0_MAX)continue;
			
			fit_2D_track(nhits_X3, fitarrays[index].x_array, fitarrays[index].z_array, fitarrays[index].dx_array, fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B, fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, fitarrays[index].chi2_2d);
			if(fabs(fitarrays[index].output_parameters[0])>2*X0_MAX || fabs(fitarrays[index].output_parameters[1])>2*TX_MAX)continue;
			
			//if(fitarrays[index].output_parameters[0]+fitarrays[index].output_parameters_errors[0]<-X0_MAX || 
			//   fitarrays[index].output_parameters[0]-fitarrays[index].output_parameters_errors[0]>X0_MAX)continue;
			//if(fitarrays[index].output_parameters[1]+fitarrays[index].output_parameters_errors[1]<-TX_MAX || 
			//   fitarrays[index].output_parameters[1]-fitarrays[index].output_parameters_errors[1]>TX_MAX)continue;
					
			//prop matching
			nprop = 0;
			for(short ip = 0; ip<4; ip++){
				iprop = 48+ip;
				xExp = fitarrays[index].output_parameters[1]*planes[iprop].z+fitarrays[index].output_parameters[0];
				//loop on hits to find prop
				for(int k = 0; k<ic[index].nAH; k++){
					if(ic[index].AllHits[k].detectorID==iprop){
						if(fabs(ic[index].AllHits[k].pos-xExp)<5.08f){
							nprop++;
							break;
						}
					}
				}
				if(nprop>0)break;
			}
			if(nprop==0)continue;
			
			//fill in the XZ track
			straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].nXHits = nhits_X3;
			straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].x0 = fitarrays[index].output_parameters[0];
			straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].err_x0 = fitarrays[index].output_parameters_errors[0];
			straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].tx = fitarrays[index].output_parameters[1];
			straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].err_tx = fitarrays[index].output_parameters_errors[1];

			straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].chisq = fitarrays[index].chi2_2d;
			
			straighttrackbuilder[index].nTracksXZ++;
			
			
		}// loops on st3 X hits
	
	}//loop on st2 X hits

		
	// 3- get U, V pairs compatible with the X pair in each station: already implemented, tested;
        // 4- U, V pairs also have to be compatible with the X pair slope => this last point needs more thought, and will probably need adjustments (2-3 days to implement and test);
	// don't we just want U, V pairs compatible with the track slope?
	float xmin, xmax;

#ifdef KTRACKERHITSELECTION
	float umin2, umax2;
	float umin3p, umax3p;
	float umin3m, umax3m;
	float vmin2, vmax2;
	float vmin3p, vmax3p;
	float vmin3m, vmax3m;
	float xpos2, xpos3p, xpos3m;
	float upos2, upos3p, upos3m;
	float vpos2, vpos3p, vpos3m;
	
	float z_x2, z_u2, z_v2;
	float z_x3, z_u3, z_v3;
	float v_win1, v_win2, v_win3, v_win;
#endif

	int nu2, nu3p, nu3m;
	int nv2, nv3p, nv3m;
	int nxhits;
#ifdef BEST_TRACKYZ_ONLY
	int bestYZcand = -1;
#endif
	float y, err_y;
	// loop on XZ tracks
	float chi2min = 100000.;//-1.0f;
	bool yz_in_acceptance;
	
	for(int i = 0; i<straighttrackbuilder[index].nTracksXZ; i++){
		nxhits = straighttrackbuilder[index].TrackXZ[i].nXHits;
#ifdef BEST_TRACKYZ_ONLY
		bestYZcand = -1;
#endif
		chi2min = 100000.;//-1.0f;
		straighttrackbuilder[index].nTracksYZ = 0;

#ifdef KTRACKERHITSELECTION
		xpos2 = xpos3p = xpos3m = 0;
		short n2 = 0, n3p = 0, n3m = 0;
		for(int n = 0; n<nxhits; n++){
			FillFitArrays(n, ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]], 0, fitarrays[index], planes);

			if(ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]].detectorID==geometry::detsuperid[2][0]*2 ||
			   ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]].detectorID==geometry::detsuperid[2][0]*2-1){
			   xpos2+=ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]].pos;
			   n2++;
			}
			if(ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]].detectorID==geometry::detsuperid[3][0]*2 ||
			   ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]].detectorID==geometry::detsuperid[3][0]*2-1){
			   xpos3p+=ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]].pos;
			   n3p++;
			}
			if(ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]].detectorID==geometry::detsuperid[4][0]*2 ||
			   ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]].detectorID==geometry::detsuperid[4][0]*2-1){
			   xpos3m+=ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]].pos;
			   n3m++;
			}
		}
		
		xpos2 = xpos2/n2;
		if(n3p)xpos3p = xpos3p/n3p;
		if(n3m)xpos3m = xpos3m/n3m;

		umin2 = xpos2*planes[geometry::detsuperid[2][1]*2].costheta-planes[geometry::detsuperid[2][1]*2].u_win;
		umax2 = umin2+2*planes[geometry::detsuperid[2][1]*2].u_win;
		
		umin3p = xpos3p*planes[geometry::detsuperid[3][1]*2].costheta-planes[geometry::detsuperid[3][1]*2].u_win;
		umax3p = umin3p+2*planes[geometry::detsuperid[3][1]*2].u_win;
		
		umin3m = xpos3m*planes[geometry::detsuperid[4][1]*2].costheta-planes[geometry::detsuperid[4][1]*2].u_win;
		umax3m = umin3m+2*planes[geometry::detsuperid[4][1]*2].u_win;
#endif

		//straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].nUHits = 0;
		//straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].nVHits = 0;

		//if(ic[index].EventID!=13)return;
		
		projid = 1;
		stid = 2;
		find_xmin_xmax_in_chamber(xmin, xmax, straighttrackbuilder[index].TrackXZ[i], stid, projid, planes);
		nu2 = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_u2, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);

		stid = 3;
		find_xmin_xmax_in_chamber(xmin, xmax, straighttrackbuilder[index].TrackXZ[i], stid, projid, planes);
		nu3p = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_u3p, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);

		stid = 4;
		find_xmin_xmax_in_chamber(xmin, xmax, straighttrackbuilder[index].TrackXZ[i], stid, projid, planes);
		nu3m = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_u3m, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);
		
		projid = 2;
		stid = 2;
		find_xmin_xmax_in_chamber(xmin, xmax, straighttrackbuilder[index].TrackXZ[i], stid, projid, planes);
		nv2 = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_v2, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);
		
		stid = 3;
		find_xmin_xmax_in_chamber(xmin, xmax, straighttrackbuilder[index].TrackXZ[i], stid, projid, planes);
		nv3p = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_v3p, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);
		
		stid = 4;
		find_xmin_xmax_in_chamber(xmin, xmax, straighttrackbuilder[index].TrackXZ[i], stid, projid, planes);
		nv3m = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_v3m, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);
		
		if(nu2==0 || nv2==0 || nu3p+nu3m==0 || nv3p+nv3m==0) continue;
		
		// build tracks with UHits only first?
		for(int j = 0; j<nu2; j++){
			for(int n = 0; n<nChamberPlanes; n++){
				fitarrays[index].y_array[n] = 0;
				fitarrays[index].z_array[n] = 0;
				fitarrays[index].dy_array[n] = 0;
			}

#ifdef KTRACKERHITSELECTION
			upos2 = straighttrackbuilder[index].hitpairs_u2[j].second>=0 ? 0.5f*(ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].first].pos+ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].second].pos) : ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].first].pos;
			
			if(upos2<umin2 || upos2>umax2)continue;
			
			z_x2 = n2==2 ? planes[geometry::detsuperid[2][0]*2].z_mean : planes[geometry::detsuperid[2][0]*2].z; 
			z_u2 = straighttrackbuilder[index].hitpairs_u2[j].second>=0 ? planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].first].detectorID].z_mean : planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].first].detectorID].z; 
			z_v2 = planes[geometry::detsuperid[2][2]*2].z_mean;
			
			v_win1 = planes[ ic[index].AllHits[ straighttrackbuilder[index].hitpairs_u2[j].first ].detectorID].v_win_fac1;
			v_win2 = fabs(z_u2+z_v2-2*z_x2)*planes[ geometry::detsuperid[2][2]*2 ].v_win_fac2;
			v_win3 = fabs((z_v2-z_u2)*planes[ geometry::detsuperid[2][2]*2 ].v_win_fac3);
			v_win = v_win1+v_win2+v_win3+2*planes[ ic[index].AllHits[ straighttrackbuilder[index].hitpairs_u2[j].first ].detectorID].spacing;
			
			vmin2 = 2*xpos2*planes[geometry::detsuperid[2][1]*2].costheta-upos2-v_win;
			vmax2 = vmin2+2*v_win;
#endif
			
			nhits_U2 = 0;
			//if(j<nu2){
			if(straighttrackbuilder[index].hitpairs_u2[j].first>=0){
				// 5- calculate Y from U, V hits 
				if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].first], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
					FillFitArrays(nxhits+nhits_U2, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].first], 0, fitarrays[index], planes);
					fitarrays[index].y_array[nhits_U2] = y;
					fitarrays[index].dy_array[nhits_U2] = err_y;
					fitarrays[index].z_array[nhits_U2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].first].detectorID].z;
					straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_U2] = straighttrackbuilder[index].hitpairs_u2[j].first;
					nhits_U2++;
				}
			}
			if(straighttrackbuilder[index].hitpairs_u2[j].second>=0){
				if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].second], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
					FillFitArrays(nxhits+nhits_U2, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].second], 0, fitarrays[index], planes);
					fitarrays[index].y_array[nhits_U2] = y;
					fitarrays[index].dy_array[nhits_U2] = err_y;
					fitarrays[index].z_array[nhits_U2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].second].detectorID].z;
					straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_U2] = straighttrackbuilder[index].hitpairs_u2[j].second;
					nhits_U2++;
				}
			}
			//}
			
			for(int k = 0; k<nu3p+nu3m; k++){
				nhits_U3 = nhits_U2;
				//if(k<nu3p+nu3m){
				if(k<nu3p){

#ifdef KTRACKERHITSELECTION
					upos3p = straighttrackbuilder[index].hitpairs_u3p[k].second>=0 ? 0.5f*(ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].first].pos+ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].second].pos) : ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].first].pos;
			
					if(upos3p<umin3p || upos3p>umax3p)continue;

					z_x3 = n3p==2 ? planes[geometry::detsuperid[3][0]*2].z_mean : planes[geometry::detsuperid[3][0]*2].z; 
					z_u3 = straighttrackbuilder[index].hitpairs_u3p[k].second>=0 ? planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].first].detectorID].z_mean : planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].first].detectorID].z; 
					z_v3 = planes[geometry::detsuperid[3][2]*2].z_mean;
			
					v_win1 = planes[ ic[index].AllHits[ straighttrackbuilder[index].hitpairs_u3p[k].first ].detectorID].v_win_fac1;
					v_win2 = fabs(z_u3+z_v3-2*z_x3)*planes[ geometry::detsuperid[3][2]*2 ].v_win_fac2;
					v_win3 = fabs((z_v3-z_u3)*planes[ geometry::detsuperid[3][2]*2 ].v_win_fac3);
					v_win = v_win1+v_win2+v_win3+2*planes[ ic[index].AllHits[ straighttrackbuilder[index].hitpairs_u3p[0].first ].detectorID].spacing;
			
					vmin3p = 2*xpos3p*planes[geometry::detsuperid[3][1]*2].costheta-upos3p-v_win;
 					vmax3p = vmin3p+2*v_win;
#endif
							
					if(straighttrackbuilder[index].hitpairs_u3p[k].first>=0){
						if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].first], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
							FillFitArrays(nxhits+nhits_U3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].first], 0, fitarrays[index], planes);
							fitarrays[index].y_array[nhits_U3] = y;
							fitarrays[index].dy_array[nhits_U3] = err_y;
							fitarrays[index].z_array[nhits_U3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].first].detectorID].z;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_U3] = straighttrackbuilder[index].hitpairs_u3p[k].first;
							nhits_U3++;
						}
					}
					if(straighttrackbuilder[index].hitpairs_u3p[k].second>=0){
						if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].second], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
							FillFitArrays(nxhits+nhits_U3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].second], 0, fitarrays[index], planes);
							fitarrays[index].y_array[nhits_U3] = y;
							fitarrays[index].dy_array[nhits_U3] = err_y;
							fitarrays[index].z_array[nhits_U3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].second].detectorID].z;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_U3] = straighttrackbuilder[index].hitpairs_u3p[k].second;
							nhits_U3++;
						}
					}
				}else{

#ifdef KTRACKERHITSELECTION
					upos3m = straighttrackbuilder[index].hitpairs_u3m[k].second>=0 ? 0.5f*(ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k].first].pos+ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k].second].pos) : ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k].first].pos;
			
					if(upos3m<umin3m || upos3m>umax3m)continue;

					z_x3 = n3m==2 ? planes[geometry::detsuperid[4][0]*2].z_mean : planes[geometry::detsuperid[4][0]*2].z; 
					z_u3 = straighttrackbuilder[index].hitpairs_u3m[k].second>=0 ? planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k].first].detectorID].z_mean : planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k].first].detectorID].z; 
					z_v3 = planes[geometry::detsuperid[4][2]*2].z_mean;
			
					v_win1 = planes[ ic[index].AllHits[ straighttrackbuilder[index].hitpairs_u3m[k].first ].detectorID].v_win_fac1;
					v_win2 = fabs(z_u3+z_v3-2*z_x3)*planes[ geometry::detsuperid[4][2]*2 ].v_win_fac2;
					v_win3 = fabs((z_v3-z_u3)*planes[ geometry::detsuperid[4][2]*2 ].v_win_fac3);
					v_win = v_win1+v_win2+v_win3+2*planes[ ic[index].AllHits[ straighttrackbuilder[index].hitpairs_u3m[0].first ].detectorID].spacing;
			
					vmin3m = 2*xpos3m*planes[geometry::detsuperid[4][1]*2].costheta-upos3m-v_win;
 					vmax3m = vmin3m+2*v_win;

					//if(ic[index].EventID<=20){
					//	printf("%d %1.6f %1.6f %1.6f \n", ic[index].EventID, z_x3, z_u3, z_v3);
					//	printf("%d %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f \n", ic[index].EventID, vmin3m, vmax3m, v_win, v_win1, v_win2, v_win3);
					//}

					//if(ic[index].EventID<=20){
					//printf("%d %1.6f %1.6f %1.6f %1.6f %1.6f %1.6f \n", ic[index].EventID, vmin3m, vmax3m, v_win, v_win1, v_win2, v_win3);
					//}
#endif
					
					if(straighttrackbuilder[index].hitpairs_u3m[k].first>=0){
						if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k].first], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
							FillFitArrays(nxhits+nhits_U3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k].first], 0, fitarrays[index], planes);
							fitarrays[index].y_array[nhits_U3] = y;
							fitarrays[index].dy_array[nhits_U3] = err_y;
							fitarrays[index].z_array[nhits_U3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k].first].detectorID].z;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_U3] = straighttrackbuilder[index].hitpairs_u3m[k].first;
							nhits_U3++;
						}
					}
					if(straighttrackbuilder[index].hitpairs_u3m[k].second>=0){
						if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k].second], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
							FillFitArrays(nxhits+nhits_U3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k].second], 0, fitarrays[index], planes);
							fitarrays[index].y_array[nhits_U3] = y;
							fitarrays[index].dy_array[nhits_U3] = err_y;
							fitarrays[index].z_array[nhits_U3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k].second].detectorID].z;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].nUHits++;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_U3] = straighttrackbuilder[index].hitpairs_u3m[k].second;
							nhits_U3++;
						}
					}
				}	
				//}
				
				for(int l = 0; l<nv2; l++){
					nhits_V2 = nhits_U3;
					
#ifdef KTRACKERHITSELECTION
					vpos2 = straighttrackbuilder[index].hitpairs_v2[l].second>=0 ? 0.5f*(ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].first].pos+ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].second].pos) : ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].first].pos;
					
					if(vpos2<vmin2 || vpos2>vmax2)continue;
#endif

					//if(l<nv2){
					if(straighttrackbuilder[index].hitpairs_v2[l].first>=0){
						if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].first], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
							FillFitArrays(nxhits+nhits_V2, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].first], 0, fitarrays[index], planes);
							fitarrays[index].y_array[nhits_V2] = y;
							fitarrays[index].dy_array[nhits_V2] = err_y;
							fitarrays[index].z_array[nhits_V2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].first].detectorID].z;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_V2] = straighttrackbuilder[index].hitpairs_v2[l].first;
							nhits_V2++;
						}
					}

					if(straighttrackbuilder[index].hitpairs_v2[l].second>=0){
						if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].second], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
							FillFitArrays(nxhits+nhits_V2, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].second], 0, fitarrays[index], planes);
							fitarrays[index].y_array[nhits_V2] = y;
							fitarrays[index].dy_array[nhits_V2] = err_y;
							fitarrays[index].z_array[nhits_V2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].second].detectorID].z;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_V2] = straighttrackbuilder[index].hitpairs_v2[l].second;
							nhits_V2++;
						}
					}
					//}
					for(int m = 0; m<nv3p+nv3m; m++){
						nhits_V3 = nhits_V2;
						if(m<nv3p){

#ifdef KTRACKERHITSELECTION
							vpos3p = straighttrackbuilder[index].hitpairs_v3p[m].second>=0 ? 0.5f*(ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].first].pos+ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].second].pos) : ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].first].pos;

							if(vpos3p<vmin3p || vpos3p>vmax3p)continue;
#endif

							if(straighttrackbuilder[index].hitpairs_v3p[m].first>=0){
								if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].first], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
									FillFitArrays(nxhits+nhits_V3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].first], 0, fitarrays[index], planes);
									fitarrays[index].y_array[nhits_V3] = y;
									fitarrays[index].dy_array[nhits_V3] = err_y;
									fitarrays[index].z_array[nhits_V3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].first].detectorID].z;
									straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_V3] = straighttrackbuilder[index].hitpairs_v3p[m].first;
									nhits_V3++;
								}
							}
						if(straighttrackbuilder[index].hitpairs_v3p[m].second>=0){
								if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].second], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
									FillFitArrays(nxhits+nhits_V3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].second], 0, fitarrays[index], planes);
									fitarrays[index].y_array[nhits_V3] = y;
									fitarrays[index].dy_array[nhits_V3] = err_y;
									fitarrays[index].z_array[nhits_V3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].second].detectorID].z;
									straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_V3] = straighttrackbuilder[index].hitpairs_v3p[m].second;
									nhits_V3++;
								}
							}
						}else{

#ifdef KTRACKERHITSELECTION
							vpos3m = straighttrackbuilder[index].hitpairs_v3m[j].second>=0 ? 0.5f*(ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[j].first].pos+ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[j].second].pos) : ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[j].first].pos;
							
							if(vpos3m<vmin3m || vpos3m>vmax3m)continue;
#endif

							if(straighttrackbuilder[index].hitpairs_v3m[m].first>=0){
								if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m].first], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
									FillFitArrays(nxhits+nhits_V3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m].first], 0, fitarrays[index], planes);
									fitarrays[index].y_array[nhits_V3] = y;
									fitarrays[index].dy_array[nhits_V3] = err_y;
									fitarrays[index].z_array[nhits_V3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m].first].detectorID].z;
									straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_V3] = straighttrackbuilder[index].hitpairs_v3m[m].first;
									nhits_V3++;
								}
							}
							if(straighttrackbuilder[index].hitpairs_v3m[m].second>=0){
								if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m].second], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
									FillFitArrays(nxhits+nhits_V3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m].second], 0, fitarrays[index], planes);
									fitarrays[index].y_array[nhits_V3] = y;
									fitarrays[index].dy_array[nhits_V3] = err_y;
									fitarrays[index].z_array[nhits_V3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m].second].detectorID].z;
									straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_V3] = straighttrackbuilder[index].hitpairs_v3m[m].second;
									nhits_V3++;
								}
							}
						}
						
						// 6- fit the X-Y slope
						
						//fit here
						fit_2D_track(nhits_V3, fitarrays[index].y_array, fitarrays[index].z_array, fitarrays[index].dy_array, fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B, fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, fitarrays[index].chi2_2d);
						
						yz_in_acceptance = true;
						if(fitarrays[index].output_parameters[0]+fitarrays[index].output_parameters_errors[0]<-Y0_MAX || 
						     fitarrays[index].output_parameters[0]-fitarrays[index].output_parameters_errors[0]>Y0_MAX)yz_in_acceptance = false;
						if(fitarrays[index].output_parameters[1]+fitarrays[index].output_parameters_errors[1]<-TY_MAX || 
						   fitarrays[index].output_parameters[1]-fitarrays[index].output_parameters_errors[1]>TY_MAX)yz_in_acceptance = false;
						
						straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].nUHits = nhits_U3;
						straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].nVHits = nhits_V3-nhits_U3;
						
						straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].y0 = fitarrays[index].output_parameters[0];
						straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].ty = fitarrays[index].output_parameters[1];
						straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].err_y0 = fitarrays[index].output_parameters_errors[0];
						straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].err_ty = fitarrays[index].output_parameters_errors[1];
						
						fitarrays[index].output_parameters[0] = straighttrackbuilder[index].TrackXZ[i].x0;
						fitarrays[index].output_parameters[1] = straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].y0;
						fitarrays[index].output_parameters[2] = straighttrackbuilder[index].TrackXZ[i].tx;
						fitarrays[index].output_parameters[3] = straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].ty;
						
						fitarrays[index].output_parameters_errors[0] = straighttrackbuilder[index].TrackXZ[i].err_x0;
						fitarrays[index].output_parameters_errors[1] = straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].err_y0;
						fitarrays[index].output_parameters_errors[2] = straighttrackbuilder[index].TrackXZ[i].err_tx;
						fitarrays[index].output_parameters_errors[3] = straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].err_ty;
						
						nhits_V3+= nxhits;
						
						// then, calculate the "ktracker" chi2 for the back partial track candidate
						chi2_straight(nhits_V3, fitarrays[index].drift_dist, fitarrays[index].resolution,
									fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z,
									fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz,
									fitarrays[index].output_parameters, fitarrays[index].chi2);
						
						
						//what happens if we keep all YZ tracks???
						//if(fitarrays[index].chi2>9000)continue;// allowing for even looser track parameters
						
#ifdef BEST_TRACKYZ_ONLY
						bestYZcand = straighttrackbuilder[index].nTracksYZ;
#endif
						oc[index].AllTracklets[oc[index].nTracklets].stationID = 5;
						oc[index].AllTracklets[oc[index].nTracklets].nXHits = straighttrackbuilder[index].TrackXZ[i].nXHits;
						oc[index].AllTracklets[oc[index].nTracklets].nUHits = straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].nUHits;
						oc[index].AllTracklets[oc[index].nTracklets].nVHits = straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].nVHits;
						
						oc[index].AllTracklets[oc[index].nTracklets].x0 = straighttrackbuilder[index].TrackXZ[i].x0;
						oc[index].AllTracklets[oc[index].nTracklets].y0 = straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].y0;
						oc[index].AllTracklets[oc[index].nTracklets].tx = straighttrackbuilder[index].TrackXZ[i].tx;
						oc[index].AllTracklets[oc[index].nTracklets].ty = straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].ty;
							
						oc[index].AllTracklets[oc[index].nTracklets].err_x0 = straighttrackbuilder[index].TrackXZ[i].err_x0;
						oc[index].AllTracklets[oc[index].nTracklets].err_y0 = straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].err_y0;
						oc[index].AllTracklets[oc[index].nTracklets].err_tx = straighttrackbuilder[index].TrackXZ[i].err_tx;
						oc[index].AllTracklets[oc[index].nTracklets].err_ty = straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].err_ty;
						
						oc[index].AllTracklets[oc[index].nTracklets].chisq = fitarrays[index].chi2;
						
						if(!match_tracklet_to_hodo(oc[index].AllTracklets[oc[index].nTracklets], 2, ic, planes))continue;
						if(!match_tracklet_to_hodo(oc[index].AllTracklets[oc[index].nTracklets], 3, ic, planes) &&
						   !match_tracklet_to_hodo(oc[index].AllTracklets[oc[index].nTracklets], 4, ic, planes))continue;
						
						for(int n = 0; n<nxhits; n++){
							oc[index].AllTracklets[oc[index].nTracklets].hits[n] = ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]];
						}
						for(int n = nxhits; n<nxhits+oc[index].AllTracklets[oc[index].nTracklets].nUHits+oc[index].AllTracklets[oc[index].nTracklets].nVHits; n++){
							oc[index].AllTracklets[oc[index].nTracklets].hits[n] = ic[index].AllHits[straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[n-nxhits]];
						}
						
						resolve_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes, 40.);
						resolve_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes, 150.);
						
						//refit after???

						oc[index].nTracklets++;
						

						if(fitarrays[index].chi2<=chi2min){
							chi2min = fitarrays[index].chi2;
#ifdef BEST_TRACKYZ_ONLY
							bestYZcand = straighttrackbuilder[index].nTracksYZ;
#endif
						}
						straighttrackbuilder[index].nTracksYZ++;
						//
					}// end loop on st3 v hits
				}// end loop on st2 v hits
			}// end loop on st3 u hits
		}// end loop on st2 u hits
		
		//if(straighttrackbuilder[index].nTracksYZ==0)printf("evt %d, %d pairs u2, %d pairs v2, %d pairs v2, %d pairs v3\n", ic[index].EventID, nu2, nu3p+nu3m, nv2, nv3p+nv3m);
		
#ifdef BEST_TRACKYZ_ONLY
		if(straighttrackbuilder[index].nTracksYZ>0 && bestYZcand>=0){
			oc[index].AllTracklets[oc[index].nTracklets].stationID = 5;
			oc[index].AllTracklets[oc[index].nTracklets].nXHits = straighttrackbuilder[index].TrackXZ[i].nXHits;
			oc[index].AllTracklets[oc[index].nTracklets].nUHits = straighttrackbuilder[index].TrackYZ[bestYZcand].nUHits;
			oc[index].AllTracklets[oc[index].nTracklets].nVHits = straighttrackbuilder[index].TrackYZ[bestYZcand].nVHits;

			oc[index].AllTracklets[oc[index].nTracklets].x0 = straighttrackbuilder[index].TrackXZ[i].x0;
			oc[index].AllTracklets[oc[index].nTracklets].y0 = straighttrackbuilder[index].TrackYZ[bestYZcand].y0;
			oc[index].AllTracklets[oc[index].nTracklets].tx = straighttrackbuilder[index].TrackXZ[i].tx;
			oc[index].AllTracklets[oc[index].nTracklets].ty = straighttrackbuilder[index].TrackYZ[bestYZcand].ty;
			
			oc[index].AllTracklets[oc[index].nTracklets].err_x0 = straighttrackbuilder[index].TrackXZ[i].err_x0;
			oc[index].AllTracklets[oc[index].nTracklets].err_y0 = straighttrackbuilder[index].TrackYZ[bestYZcand].err_y0;
			oc[index].AllTracklets[oc[index].nTracklets].err_tx = straighttrackbuilder[index].TrackXZ[i].err_tx;
			oc[index].AllTracklets[oc[index].nTracklets].err_ty = straighttrackbuilder[index].TrackYZ[bestYZcand].err_ty;
			
			oc[index].AllTracklets[oc[index].nTracklets].chisq = chi2min;
			
			if(!match_tracklet_to_hodo(oc[index].AllTracklets[oc[index].nTracklets], 2, ic, planes))continue;
			//if(!match_tracklet_to_hodo(oc[index].AllTracklets[oc[index].nTracklets], 3, ic, planes) &&
			//   !match_tracklet_to_hodo(oc[index].AllTracklets[oc[index].nTracklets], 4, ic, planes))continue;
			
			match_tracklet_to_hodo(oc[index].AllTracklets[oc[index].nTracklets], 3, ic, planes);
			match_tracklet_to_hodo(oc[index].AllTracklets[oc[index].nTracklets], 4, ic, planes);
			
			for(int n = 0; n<nxhits; n++){
				oc[index].AllTracklets[oc[index].nTracklets].hits[n] = ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]];
			}
			if(oc[index].AllTracklets[oc[index].nTracklets].ty==0 || oc[index].AllTracklets[oc[index].nTracklets].y0==0){
				printf("1- evt %d, bestYZcand %d, nYZ %d, chi2 %1.6f, nuhits = %d (%d+%d+%d), nvhits = %d (%d+%d+%d), x0 = %1.6f, y0 = %1.6f, tx = %1.6f, ty = %1.6f \n",
					ic[index].EventID, bestYZcand, straighttrackbuilder[index].nTracksYZ, chi2min,
					oc[index].AllTracklets[oc[index].nTracklets].nUHits, nu2, nu3p, nu3m,
					oc[index].AllTracklets[oc[index].nTracklets].nVHits, nv2, nv3p, nv3m,
					oc[index].AllTracklets[oc[index].nTracklets].x0, 
					oc[index].AllTracklets[oc[index].nTracklets].y0, 
					oc[index].AllTracklets[oc[index].nTracklets].tx, 
					oc[index].AllTracklets[oc[index].nTracklets].ty 
					);

			}

			if(nxhits+oc[index].AllTracklets[oc[index].nTracklets].nUHits+oc[index].AllTracklets[oc[index].nTracklets].nVHits>nChamberPlanes)
			printf("%d + %d + %d  > %d \n", nxhits, oc[index].AllTracklets[oc[index].nTracklets].nUHits, oc[index].AllTracklets[oc[index].nTracklets].nVHits, nChamberPlanes);
			for(int n = nxhits; n<nxhits+oc[index].AllTracklets[oc[index].nTracklets].nUHits+oc[index].AllTracklets[oc[index].nTracklets].nVHits; n++){
				oc[index].AllTracklets[oc[index].nTracklets].hits[n] = ic[index].AllHits[straighttrackbuilder[index].TrackYZ[bestYZcand].hitlist[n-nxhits]];
			}
			if(oc[index].AllTracklets[oc[index].nTracklets].ty==0 || oc[index].AllTracklets[oc[index].nTracklets].y0==0){
				printf("2- evt %d, bestYZcand %d, nYZ %d, chi2 %1.6f, nuhits = %d (%d+%d+%d), nvhits = %d (%d+%d+%d), x0 = %1.6f, y0 = %1.6f, tx = %1.6f, ty = %1.6f \n",
					ic[index].EventID, bestYZcand, straighttrackbuilder[index].nTracksYZ, chi2min,
					oc[index].AllTracklets[oc[index].nTracklets].nUHits, nu2, nu3p, nu3m,
					oc[index].AllTracklets[oc[index].nTracklets].nVHits, nv2, nv3p, nv3m,
					oc[index].AllTracklets[oc[index].nTracklets].x0, 
					oc[index].AllTracklets[oc[index].nTracklets].y0, 
					oc[index].AllTracklets[oc[index].nTracklets].tx, 
					oc[index].AllTracklets[oc[index].nTracklets].ty 
					);

			}
			//if(ic[index].EventID==13){
			//	printf("x0 = %1.6f, y0 = %1.6f, tx = %1.6f, ty = %1.6f \n", 
			//		oc[index].AllTracklets[oc[index].nTracklets].x0, 
			//		oc[index].AllTracklets[oc[index].nTracklets].y0, 
			//		oc[index].AllTracklets[oc[index].nTracklets].tx, 
			//		oc[index].AllTracklets[oc[index].nTracklets].ty 
			//		);
			//}
			oc[index].nTracklets++;
		}else{
			//printf("evt %d: track XZ %d without YZ track: nuhits = (%d+%d+%d), nvhits = (%d+%d+%d)\n", ic[index].EventID, i, nu2, nu3p, nu3m, nv2, nv3p, nv3m);
			continue;
		}
#endif
	}// end loop on XZ tracks
	
	

}


__device__ void calculate_invP(gTracklet& tkl, const gTrackXZ tkl_st1)
{
// invP = ( tx_st1 - tx ) / ( PT_KICK_KMAG*Charge );
//      = ( ( tx*Z_KMAG_BEND + x0 - x0_st1 )/Z_KMAG_BEND - tx ) / ( PT_KICK_KMAG*Charge );
	
	tkl.invP = ( tkl_st1.tx - tkl.tx )/( geometry::PT_KICK_KMAG );
	if(tkl.invP<0){
		tkl.charge = -1;
		tkl.invP*= tkl.charge;
	}else{
		tkl.charge = +1;
	}

//Error: err_invP = err_kick/PT_KICK_KMAG
//                = (err_tx_st1 - err_tx)/PT_KICK_KMAG
	
	tkl.err_invP = ( tkl_st1.err_tx - tkl.err_tx )/( geometry::PT_KICK_KMAG );
}


__device__ void SagittaRatioInStation1(const gTracklet tkl, float* pos_exp, float* window, const gPlane* planes)
{
	float z_st3 = planes[tkl.hits[tkl.nXHits+tkl.nUHits+tkl.nVHits-1].detectorID].z;
	float x_st3 = tkl.x0+tkl.tx*z_st3;
	float y_st3 = tkl.y0+tkl.ty*z_st3;

	//printf("det id %d z %1.3f x %1.3f y %1.3f \n", tkl.hits[tkl.nXHits+tkl.nUHits+tkl.nVHits-1].detectorID, z_st3, x_st3, y_st3);
	
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
		
		pos_st3 = x_st3*planes[detid_2].costheta + y_st3*planes[detid_2].sintheta;

		//printf(" i %d idx %d  pos %1.3f \n", i, idx, pos_st3);
		
		z_st1 = planes[detid].z;
		z_st2 = planes[detid_2].z;
		
		x_st2 = tkl.x0+tkl.tx*z_st2;
		y_st2 = tkl.y0+tkl.ty*z_st2;
		
		pos_st2 = x_st2*planes[detid_2].costheta + y_st2*planes[detid_2].sintheta;
		
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



__global__ void gKernel_GlobalTracking(gEvent* ic, gOutputEvent* oc, gFullTrackBuilder* fulltrackbuilder, gStraightFitArrays* fitarrays, const gPlane* planes, bool ktracker)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;

	int N_tkl = oc[index].nTracklets;
	
	int nx1, nu1, nv1;
	float pos_exp[3];
	float window[3];
	short projid;
	short stid = 0;
	short nhits_X1, nhits_U1, nhits_V1, nhits_1;
	//float y, err_y;
	
	float umin1, umax1;
	float vmin1, vmax1;
	float xpos1, upos1, vpos1;
	
	float z_x1, z_u1, z_v1;
	float v_win1, v_win2, v_win3, v_win;
	
	if(ktracker){
	
	for(int i = 0; i<oc[index].nTracklets; i++){
		SagittaRatioInStation1(oc[index].AllTracklets[i], pos_exp, window, planes);
		
		projid = 0;
		nx1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_x1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);
		projid = 1;
		nu1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_u1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);
		projid = 2;
		nv1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_v1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);
		
		if(nx1==0 || nu1==0 || nv1==0)continue;
		
		//triple loop on hits
		for(int j = 0; j<nx1; j++){
			nhits_X1 = 0;
			//in x, we can already make a first evaluation of tx_st1, and a first evaluation of Pinv
			xpos1 = fulltrackbuilder[index].hitpairs_x1[j].second>=0 ? 0.5f*(ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].first].pos+ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].second].pos) : ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].first].pos;
			
			umin1 = xpos1*planes[geometry::detsuperid[0][1]*2].costheta-planes[geometry::detsuperid[0][1]*2].u_win;
			umax1 = umin1+2*planes[geometry::detsuperid[0][1]*2].u_win;
			
			if(fulltrackbuilder[index].hitpairs_x1[j].first>=0){
				fitarrays[index].x_array[nhits_X1] = ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].first].pos;
				fitarrays[index].dx_array[nhits_X1] = planes[ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].first].detectorID].spacing/3.4641f;
				fitarrays[index].z_array[nhits_X1] = planes[ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].first].detectorID].z;
				fulltrackbuilder[index].hitlist[nhits_X1] = fulltrackbuilder[index].hitpairs_x1[j].first;
				nhits_X1++;
			}
			if(fulltrackbuilder[index].hitpairs_x1[j].second>=0){
				fitarrays[index].x_array[nhits_X1] = ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].second].pos;
				fitarrays[index].dx_array[nhits_X1] = planes[ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].second].detectorID].spacing/3.4641f;
				fitarrays[index].z_array[nhits_X1] = planes[ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].second].detectorID].z;
				fulltrackbuilder[index].hitlist[nhits_X1] = fulltrackbuilder[index].hitpairs_x1[j].second;
				nhits_X1++;
			}
			
			//fit to obtain the st1 track segment: this time fitting from KMag, not from the origin:
			fitarrays[index].x_array[nhits_X1] = oc[index].AllTracklets[i].x0+oc[index].AllTracklets[i].tx*geometry::Z_KMAG_BEND;
			fitarrays[index].dx_array[nhits_X1] = oc[index].AllTracklets[i].err_x0+oc[index].AllTracklets[i].err_tx*geometry::Z_KMAG_BEND;
			fitarrays[index].z_array[nhits_X1] = geometry::Z_KMAG_BEND;
			
			fit_2D_track(nhits_X1+1, fitarrays[index].x_array, fitarrays[index].z_array, fitarrays[index].dx_array, fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B, fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, fitarrays[index].chi2_2d);
			
			fulltrackbuilder[index].TrackXZ_st1.x0 = fitarrays[index].output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.err_x0 = fitarrays[index].output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.tx = fitarrays[index].output_parameters[1];
			fulltrackbuilder[index].TrackXZ_st1.err_tx = fitarrays[index].output_parameters[1];
			
			for(int k = 0; k<nu1; k++){
				nhits_U1 = nhits_X1;
				upos1 = fulltrackbuilder[index].hitpairs_u1[k].second>=0 ? 0.5f*(ic[index].AllHits[fulltrackbuilder[index].hitpairs_u1[k].first].pos+ic[index].AllHits[fulltrackbuilder[index].hitpairs_u1[k].second].pos) : ic[index].AllHits[fulltrackbuilder[index].hitpairs_u1[k].first].pos;
				
				if(upos1<umin1 || upos1>umax1)continue;
				
				z_x1 = fulltrackbuilder[index].hitpairs_x1[j].second>=0 ? planes[geometry::detsuperid[0][0]*2].z_mean : planes[geometry::detsuperid[0][0]*2].z; 
				z_u1 = fulltrackbuilder[index].hitpairs_u1[k].second>=0 ? planes[ic[index].AllHits[fulltrackbuilder[index].hitpairs_u1[k].first].detectorID].z_mean : planes[ic[index].AllHits[fulltrackbuilder[index].hitpairs_u1[k].first].detectorID].z; 
				z_v1 = planes[geometry::detsuperid[0][2]*2].z_mean;
				
				v_win1 = planes[ ic[index].AllHits[ fulltrackbuilder[index].hitpairs_u1[k].first ].detectorID].v_win_fac1;
				v_win2 = fabs(z_u1+z_v1-2*z_x1)*planes[ geometry::detsuperid[0][2]*2 ].v_win_fac2;
				v_win3 = fabs((z_v1-z_u1)*planes[ geometry::detsuperid[0][2]*2 ].v_win_fac3);
				v_win = v_win1+v_win2+v_win3+2*planes[ ic[index].AllHits[ fulltrackbuilder[index].hitpairs_u1[k].first ].detectorID].spacing;
			
				vmin1 = 2*xpos1*planes[geometry::detsuperid[0][1]*2].costheta-upos1-v_win;
				vmax1 = vmin1+2*v_win;
				
				if(fulltrackbuilder[index].hitpairs_u1[k].first>=0){
					fulltrackbuilder[index].hitlist[nhits_U1] = fulltrackbuilder[index].hitpairs_u1[k].first;
					nhits_U1++;
				}
				if(fulltrackbuilder[index].hitpairs_u1[k].second>=0){
					fulltrackbuilder[index].hitlist[nhits_U1] = fulltrackbuilder[index].hitpairs_u1[k].second;
					nhits_U1++;
				}
				
				for(int l = 0; l<nv1; l++){
					nhits_V1 = nhits_U1;
					vpos1 = fulltrackbuilder[index].hitpairs_v1[l].second>=0 ? 0.5f*(ic[index].AllHits[fulltrackbuilder[index].hitpairs_v1[l].first].pos+ic[index].AllHits[fulltrackbuilder[index].hitpairs_v1[l].second].pos) : ic[index].AllHits[fulltrackbuilder[index].hitpairs_v1[l].first].pos;
					
					if(vpos1<vmin1 || vpos1>vmax1)continue;
					
					if(fulltrackbuilder[index].hitpairs_v1[l].first>=0){
						fulltrackbuilder[index].hitlist[nhits_V1] = fulltrackbuilder[index].hitpairs_v1[l].first;
						nhits_V1++;
					}
					if(fulltrackbuilder[index].hitpairs_v1[l].second>=0){
						fulltrackbuilder[index].hitlist[nhits_V1] = fulltrackbuilder[index].hitpairs_v1[l].second;
						nhits_V1++;
					}
					
					//building here the global track candidate
					
					oc[index].AllTracklets[N_tkl].stationID = 6;
					oc[index].AllTracklets[N_tkl].nXHits = oc[index].AllTracklets[i].nXHits + nhits_X1;
					oc[index].AllTracklets[N_tkl].nUHits = oc[index].AllTracklets[i].nUHits + nhits_U1-nhits_X1;
					oc[index].AllTracklets[N_tkl].nVHits = oc[index].AllTracklets[i].nVHits + nhits_V1-nhits_U1;
					
					oc[index].AllTracklets[N_tkl].x0 = oc[index].AllTracklets[i].x0;
					oc[index].AllTracklets[N_tkl].err_x0 = oc[index].AllTracklets[i].err_x0;
					oc[index].AllTracklets[N_tkl].y0 = oc[index].AllTracklets[i].y0;
					oc[index].AllTracklets[N_tkl].err_y0 = oc[index].AllTracklets[i].err_y0;
					oc[index].AllTracklets[N_tkl].tx = oc[index].AllTracklets[i].tx;
					oc[index].AllTracklets[N_tkl].err_tx = oc[index].AllTracklets[i].err_tx;
					oc[index].AllTracklets[N_tkl].ty = oc[index].AllTracklets[i].ty;
					oc[index].AllTracklets[N_tkl].err_ty = oc[index].AllTracklets[i].err_ty;
					
					calculate_invP(oc[index].AllTracklets[N_tkl], fulltrackbuilder[index].TrackXZ_st1);
					
					nhits_1 = 0;
					for(int n = 0; n<oc[index].AllTracklets[i].nXHits+oc[index].AllTracklets[i].nUHits+oc[index].AllTracklets[i].nVHits; n++){
						oc[index].AllTracklets[N_tkl].hits[nhits_1] = oc[index].AllTracklets[i].hits[n];
						nhits_1++;
					}
					for(int n = 0; n<nhits_V1; n++){
						oc[index].AllTracklets[N_tkl].hits[nhits_1] = ic[index].AllHits[fulltrackbuilder[index].hitlist[n]];
						nhits_1++;
					}
					resolve_leftright(oc[index].AllTracklets[N_tkl], planes, 75.);
					resolve_leftright(oc[index].AllTracklets[N_tkl], planes, 150.);
					
					//refit after???
					
					N_tkl++;
				}
			}//end loop on u hits
		}//end loop on x hits
	}//end loop on straight tracks
	
	}else{
	
	float x_st1_mean[2];
	float x_st1_width[2];
	
	float xmin, xmax;
	
	short hyp;
	
	projid = 1;
	projid = 2;

	for(int i = 0; i<oc[index].nTracklets; i++){
		hyp = 0;
		extrapolate_track_position_st1(oc[index].AllTracklets[i], x_st1_mean, x_st1_width, planes, hyp);
		xmin = min(x_st1_mean[0]-x_st1_width[0], x_st1_mean[1]-x_st1_width[1]);
		xmax = max(x_st1_mean[0]+x_st1_width[0], x_st1_mean[1]+x_st1_width[1]);
		
		projid = 0;
		nx1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_x1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);
		
		if(nx1==0){
			hyp = 1;
			extrapolate_track_position_st1(oc[index].AllTracklets[i], x_st1_mean, x_st1_width, planes, hyp);
			xmin = min(x_st1_mean[0]-x_st1_width[0], x_st1_mean[1]-x_st1_width[1]);
			xmax = max(x_st1_mean[0]+x_st1_width[0], x_st1_mean[1]+x_st1_width[1]);
			projid = 0;
			nx1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_x1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);
		}
		
		if(nx1==0)continue;
		//triple loop on hits
		for(int j = 0; j<nx1; j++){
			nhits_X1 = 0;
			
			if(fulltrackbuilder[index].hitpairs_x1[j].first>=0){
				fitarrays[index].x_array[nhits_X1] = ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].first].pos;
				fitarrays[index].dx_array[nhits_X1] = planes[ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].first].detectorID].spacing/3.4641f;
				fitarrays[index].z_array[nhits_X1] = planes[ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].first].detectorID].z;
				fulltrackbuilder[index].hitlist[nhits_X1] = fulltrackbuilder[index].hitpairs_x1[j].first;
				nhits_X1++;
			}
			if(fulltrackbuilder[index].hitpairs_x1[j].second>=0){
				fitarrays[index].x_array[nhits_X1] = ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].second].pos;
				fitarrays[index].dx_array[nhits_X1] = planes[ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].second].detectorID].spacing/3.4641f;
				fitarrays[index].z_array[nhits_X1] = planes[ic[index].AllHits[fulltrackbuilder[index].hitpairs_x1[j].second].detectorID].z;
				fulltrackbuilder[index].hitlist[nhits_X1] = fulltrackbuilder[index].hitpairs_x1[j].second;
				nhits_X1++;
			}
			
			//fit to obtain the st1 track segment: this time fitting from KMag, not from the origin:
			fitarrays[index].x_array[nhits_X1] = oc[index].AllTracklets[i].x0+oc[index].AllTracklets[i].tx*geometry::Z_KMAG_BEND;
			fitarrays[index].dx_array[nhits_X1] = oc[index].AllTracklets[i].err_x0+oc[index].AllTracklets[i].err_tx*geometry::Z_KMAG_BEND;
			fitarrays[index].z_array[nhits_X1] = geometry::Z_KMAG_BEND;
			
			fit_2D_track(nhits_X1+1, fitarrays[index].x_array, fitarrays[index].z_array, fitarrays[index].dx_array, fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B, fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, fitarrays[index].chi2_2d);
			
			fulltrackbuilder[index].TrackXZ_st1.x0 = fitarrays[index].output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.err_x0 = fitarrays[index].output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.tx = fitarrays[index].output_parameters[1];
			fulltrackbuilder[index].TrackXZ_st1.err_tx = fitarrays[index].output_parameters[1];

			projid = 1;
			find_xmin_xmax_in_chamber(xmin, xmax, fulltrackbuilder[index].TrackXZ_st1, stid, projid, planes);
			nu1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_u1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);
			
			projid = 2;
			find_xmin_xmax_in_chamber(xmin, xmax, fulltrackbuilder[index].TrackXZ_st1, stid, projid, planes);
			nv1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_v1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);
			
			
			//
			for(int k = 0; k<nu1; k++){
				nhits_U1 = nhits_X1;
				if(fulltrackbuilder[index].hitpairs_u1[k].first>=0){
					fulltrackbuilder[index].hitlist[nhits_U1] = fulltrackbuilder[index].hitpairs_u1[k].first;
					nhits_U1++;
				}
				if(fulltrackbuilder[index].hitpairs_u1[k].second>=0){
					fulltrackbuilder[index].hitlist[nhits_U1] = fulltrackbuilder[index].hitpairs_u1[k].second;
					nhits_U1++;
				}
				
				for(int l = 0; l<nv1; l++){
					nhits_V1 = nhits_U1;

					if(fulltrackbuilder[index].hitpairs_v1[l].first>=0){
						fulltrackbuilder[index].hitlist[nhits_V1] = fulltrackbuilder[index].hitpairs_v1[l].first;
						nhits_V1++;
					}
					if(fulltrackbuilder[index].hitpairs_v1[l].second>=0){
						fulltrackbuilder[index].hitlist[nhits_V1] = fulltrackbuilder[index].hitpairs_v1[l].second;
						nhits_V1++;
					}
					//printf("nhits V1 = %d\n", nhits_V1-nhits_U1);

					//building here the global track candidate
					oc[index].AllTracklets[N_tkl].stationID = 6;
					oc[index].AllTracklets[N_tkl].nXHits = oc[index].AllTracklets[i].nXHits + nhits_X1;
					oc[index].AllTracklets[N_tkl].nUHits = oc[index].AllTracklets[i].nUHits + nhits_U1-nhits_X1;
					oc[index].AllTracklets[N_tkl].nVHits = oc[index].AllTracklets[i].nVHits + nhits_V1-nhits_U1;
					
					oc[index].AllTracklets[N_tkl].x0 = oc[index].AllTracklets[i].x0;
					oc[index].AllTracklets[N_tkl].err_x0 = oc[index].AllTracklets[i].err_x0;
					oc[index].AllTracklets[N_tkl].y0 = oc[index].AllTracklets[i].y0;
					oc[index].AllTracklets[N_tkl].err_y0 = oc[index].AllTracklets[i].err_y0;
					oc[index].AllTracklets[N_tkl].tx = oc[index].AllTracklets[i].tx;
					oc[index].AllTracklets[N_tkl].err_tx = oc[index].AllTracklets[i].err_tx;
					oc[index].AllTracklets[N_tkl].ty = oc[index].AllTracklets[i].ty;
					oc[index].AllTracklets[N_tkl].err_ty = oc[index].AllTracklets[i].err_ty;
					
					calculate_invP(oc[index].AllTracklets[N_tkl], fulltrackbuilder[index].TrackXZ_st1);
					
					nhits_1 = 0;
					for(int n = 0; n<oc[index].AllTracklets[i].nXHits+oc[index].AllTracklets[i].nUHits+oc[index].AllTracklets[i].nVHits; n++){
						oc[index].AllTracklets[N_tkl].hits[nhits_1] = oc[index].AllTracklets[i].hits[n];
						nhits_1++;
					}
					for(int n = 0; n<nhits_V1; n++){
						oc[index].AllTracklets[N_tkl].hits[nhits_1] = ic[index].AllHits[fulltrackbuilder[index].hitlist[n]];
						nhits_1++;
					}

					resolve_leftright(oc[index].AllTracklets[N_tkl], planes, 75.);
					resolve_leftright(oc[index].AllTracklets[N_tkl], planes, 150.);
					
					//refit after???
					
					N_tkl++;
				}//end loop on v hits
			}//end loop on u hits
		}//end loop on x_hits
		//printf("nTracklets = %d \n", N_tkl);
	}//end loop on tracklets
	}	

	oc[index].nTracklets = N_tkl;
	
}





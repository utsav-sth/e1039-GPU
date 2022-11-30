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
	//ic[index].NHits[0] = NHits_reduced[index];
	
	//if(ic[index].EventID==0){
	//	printf("evt = %d, nAH = %d: \n", ic[index].EventID, ic[index].nAH);
	//	for(int i = 1; i<=nDetectors; i++)printf(" det %d: %d;  ", i, ic[index].NHits[i]);
	//	printf("\n");
	//}

	//if(((ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6])>0) || ((ic[index].NHits[7]+ic[index].NHits[8]+ic[index].NHits[9]+ic[index].NHits[10]+ic[index].NHits[11]+ic[index].NHits[12])>0) || ((ic[index].NHits[13]+ic[index].NHits[14]+ic[index].NHits[15]+ic[index].NHits[16]+ic[index].NHits[17]+ic[index].NHits[18])>0) || ((ic[index].NHits[19]+ic[index].NHits[20]+ic[index].NHits[21]+ic[index].NHits[22]+ic[index].NHits[23]+ic[index].NHits[24])>0) || ((ic[index].NHits[25]+ic[index].NHits[26]+ic[index].NHits[27]+ic[index].NHits[28]+ic[index].NHits[29]+ic[index].NHits[30])>0)){	
	//if(((ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6])<270) || ((ic[index].NHits[7]+ic[index].NHits[8]+ic[index].NHits[9]+ic[index].NHits[10]+ic[index].NHits[11]+ic[index].NHits[12])>350) || ((ic[index].NHits[13]+ic[index].NHits[14]+ic[index].NHits[15]+ic[index].NHits[16]+ic[index].NHits[17]+ic[index].NHits[18])>170) || ((ic[index].NHits[19]+ic[index].NHits[20]+ic[index].NHits[21]+ic[index].NHits[22]+ic[index].NHits[23]+ic[index].NHits[24])>140) || ((ic[index].NHits[25]+ic[index].NHits[26]+ic[index].NHits[27]+ic[index].NHits[28]+ic[index].NHits[29]+ic[index].NHits[30])>140))

	//we do not accept the event unless there is at least one hit in the first DC

	/*
	if( (ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6])<1){
		//printf("Event rejected...\n");
		}
		else {
			//counting total hit number, for all events < 6668? why? because she wanted just a subset!
			if( (ic[index].EventID)>10000 && (ic[index].EventID)<10050 ){//just look at a subset with something in it
				int totalDetectorHits = 0;
				for(int i = 1; i <= nDetectors; ++i) {
					totalDetectorHits += ic[index].NHits[i];
					//printf("%d ", ic[index].NHits[i]);
				}//printf("\n");
				
				//int nFirstRegionHits = 0;
				//for(int i = 1; i < 6; ++i) {
					//nFirstRegionHits += ic[index].NHits[i];
					//printf("nHits[%d] = %d\n", i, ic[index].NHits[i]);
				//}
			
				//printf("AllHits value : %d\n", (ic[index].NHits[0]));
				//printf("event : %d\n", (ic[index].EventID));
				//printf("all detector hit sum : %d; reduced AllHits value : %d\n", totalDetectorHits, (nAH_reduced[index]));
				//printf("sum of detectors : %d (%d)\n", totalDetectorHits, nFirstRegionHits);
				//}
			}
		    
	//		Process the accepted events (tracking) here.
			// where is the tracking though?
		}
	*/
}


// function to match a tracklet to a hodoscope hit
__device__ int match_tracklet_to_hodo(gTracklet tkl, int stID, gEvent* ic, gPlane* planes)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int masked = 0;//false
	// first, define the search region, and foremost, the planes on which we define this search region, which depends on the station ID we're looking at
	int hodoplanerange[5][2] = {{31, 34}, {31, 34}, {35, 38}, {39, 42}, {39, 42}};// range of planes to look for hits
	REAL fudgefac[5] = {0.25f, 0.25f, 0.2f, 0.15f, 0.15f};
	// define the region in which we are supposed to have hits:
	
	REAL xhodo, yhodo, xmin, xmax, ymin, ymax;
	
	// loop on the hits and select hodoscope hits corresponding to the station
	for(int i = 0; i<ic[index].nAH; i++){
		//we only consider hits in the hodoscopes planes corresponding to the station where the tracklet is reconstructed 
		if(hodoplanerange[stID][0]>ic[index].AllHits[i].detectorID || hodoplanerange[stID][1]<ic[index].AllHits[i].detectorID)continue;
		
		//calculate the track position at the hodoscope plane z
		xhodo = planes[ic[index].AllHits[i].detectorID-1].z*tkl.tx+tkl.x0;
		yhodo = planes[ic[index].AllHits[i].detectorID-1].z*tkl.ty+tkl.y0;
		
		//calculate "xmin, xmax, ymin, ymax" in which the track is supposed to pass through; 
		//these are basicially defined as the spatial coverage of the hit hodoscope element (plus some fudge factor for x)
		if(planes[ic[index].AllHits[i].detectorID-1].costheta>0.99){
			xmin = ic[index].AllHits[i].pos-planes[ic[index].AllHits[i].detectorID-1].spacing*(0.5f+fudgefac[stID])-planes[ic[index].AllHits[i].detectorID-1].z*tkl.err_tx+tkl.err_x0;
			xmax = ic[index].AllHits[i].pos+planes[ic[index].AllHits[i].detectorID-1].spacing*(0.5f+fudgefac[stID])+planes[ic[index].AllHits[i].detectorID-1].z*tkl.err_tx+tkl.err_x0;
			ymax = 0.5f*planes[ic[index].AllHits[i].detectorID-1].scaley+planes[ic[index].AllHits[i].detectorID-1].z*tkl.err_ty+tkl.err_y0;
			ymin = -ymax;
		}else{
			xmax = 0.5f*planes[ic[index].AllHits[i].detectorID-1].scalex+planes[ic[index].AllHits[i].detectorID-1].z*tkl.err_tx+tkl.err_x0;
			xmin = -xmax;
			ymin = ic[index].AllHits[i].pos-planes[ic[index].AllHits[i].detectorID-1].spacing-planes[ic[index].AllHits[i].detectorID-1].z*tkl.err_ty+tkl.err_y0;
			ymax = ic[index].AllHits[i].pos+planes[ic[index].AllHits[i].detectorID-1].spacing+planes[ic[index].AllHits[i].detectorID-1].z*tkl.err_ty+tkl.err_y0;
		}

		//fudge = planes[ic[index].AllHits[i].detectorID-1].spacing*fudgefac[stID];
		if(xmin <= xhodo && xhodo <= xmax && ymin <= yhodo && yhodo <= ymax ){
			masked++;
		}
		
	}
	// 
	return masked>0;
}




// function to make the hit pairs in station;
// I assume it will only be called by the tracklet builder
// (not by the main function), so I can make it a "device" function. 
__device__ int make_hitpairs_in_station(gEvent* ic, thrust::pair<int, int>* hitpairs, int* hitidx1, int* hitidx2, short* hitflag1, short* hitflag2, int stID, int projID, gPlane* planes){
	   //TODO: fix memory management: there's a lot of stuff declared here (though it's probably destroyed after so it *might* be fine)
	   int npairs = 0;
	   
	   // I think we assume that by default we want to know where we are
	   int index = threadIdx.x + blockIdx.x * blockDim.x;

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


// tracklet in station builder: 
__global__ void gkernel_TrackletinStation(gEvent* ic, gSW* oc, gFitArrays* fitarrays, int stID, gPlane* planes, gFitParams* fitparams) {
	//int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
	//int index = threadIdx.x + blockId * blockDim.x;
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	stID--;
	oc[index].EventID = ic[index].EventID;
	oc[index].nAH = ic[index].nAH;
	
	//if(10000<ic[index].EventID && ic[index].EventID<10050){
	//	for(int m = 0; m<30; m++){
	//		if(planes[m].u_win!=0)printf("index= %d, m = %d, u_win = %1.6f, costheta = %1.6f\n", index, m, planes[m].u_win, planes[m].costheta);
	//	}
	//}
	// loop on hits
	//if( (ic[index].EventID)>10000 && (ic[index].EventID)<10100 ){//just look at a subset with something in it
	//	printf("core idx %d, evt %d: reduced AllHits value : %d\n", (index), ic[index].EventID, (Nhits));
	//}
	// answer is yes, we still have the info from the previous function i.e. running after eR we still benefit from hit reduction;
	// was worth checking, just in case...

	//we don't need pairs of *HITS* necessarily, we just need pairs of indices...
	//thrust::pair<int, int> hitpairs_x[100];
	//thrust::pair<int, int> hitpairs_u[100];
	//thrust::pair<int, int> hitpairs_v[100];

	int nx = make_hitpairs_in_station(ic, oc[index].hitpairs_x, oc[index].hitidx1, oc[index].hitidx2, oc[index].hitflag1, oc[index].hitflag2, stID, 0, planes);
	int nu = make_hitpairs_in_station(ic, oc[index].hitpairs_u, oc[index].hitidx1, oc[index].hitidx2, oc[index].hitflag1, oc[index].hitflag2, stID, 1, planes);
	int nv = make_hitpairs_in_station(ic, oc[index].hitpairs_v, oc[index].hitidx1, oc[index].hitidx2, oc[index].hitflag1, oc[index].hitflag2, stID, 2, planes);
	
	short uidx = stID==0? stID*6 : stID*6+4;
	short vidx = stID==0? stID*6+4 : stID*6;
	
	//bool print = false;
	//if(0<=ic[index].EventID && ic[index].EventID<20){
	//	print = true;
	//printf("evt %d, nx = %d, nu = %d, nv = %d, ucostheta(plane %d) = %1.6f, uwin(plane %d) = %1.6f\n", ic[index].EventID, nx, nu, nv, uidx, planes[uidx].costheta, uidx, planes[uidx].u_win);
	//	printf("evt %d, nx = %d, nu = %d, nv = %d\n", ic[index].EventID, nx, nu, nv);
	//}
	//one has to have at least one hit in x, u, v
	if(nx==0 || nu==0 || nv==0)return;
	int n_tkl = oc[index].nTracklets;

	int nhits_tkl = 0;
	int npts = 0;

	REAL fixedpoint[3] = {0, 0, 0};
		
	//X-U combinations first
	for(int i = 0; i< nx; i++){
		double xpos = oc[index].hitpairs_x[i].second>=0 ? 0.5*(ic[index].AllHits[ oc[index].hitpairs_x[i].first ].pos+ic[index].AllHits[ oc[index].hitpairs_x[i].second ].pos): ic[index].AllHits[ oc[index].hitpairs_x[i].first ].pos;
		double umin = xpos*planes[uidx].costheta-planes[uidx].u_win;
		double umax = umin+2*planes[uidx].u_win;
		//if(print){
		//	printf("evt %d, xpos = %1.6f, umin = %1.6f, umax = %1.6f\n", ic[index].EventID, xpos, umin, umax);
		//	printf("evt %d, x1 pos = %1.6f, x2 pos =%1.6f\n", ic[index].EventID, ic[index].AllHits[ oc[index].hitpairs_x[i].first ].pos, 
		//		oc[index].hitpairs_x[i].second >=0 ? ic[index].AllHits[ oc[index].hitpairs_x[i].second ].pos : -1000000);
		//}
		for(int j = 0; j< nu; j++){
			double upos = oc[index].hitpairs_u[j].second>=0 ? 0.5*(ic[index].AllHits[ oc[index].hitpairs_u[j].first ].pos+ic[index].AllHits[ oc[index].hitpairs_u[j].second ].pos): ic[index].AllHits[ oc[index].hitpairs_u[j].first ].pos;

			if(upos<umin || upos>umax)continue;
			//if(print)printf("evt %d, %1.6f <? upos = %1.6f <? %1.6f \n", ic[index].EventID, umin, upos, umax);
			//we chose by convention to start the numbering of the "planes" object arrays to 0 instead of 1, this is why we have to subtract 1 to detectorID to get the correct information
			double z_x = oc[index].hitpairs_x[i].second>=0 ? planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].z_mean : planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].z;
			double z_u = oc[index].hitpairs_u[j].second>=0 ? planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].z_mean : planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].z;
			double z_v = planes[vidx].z_mean;
			//if(ic[index].EventID==0)printf("detid x = %d, detid u = %d, z_x = %1.6f, z_u = %1.6f, z_v = %1.6f\n", ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID, ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID, z_x, z_u, z_v);
			double v_win1 = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].v_win_fac1;
			double v_win2 = fabs(z_u+z_v-2*z_x)*planes[ vidx ].v_win_fac2;
			double v_win3 = fabs(z_v-z_u)*planes[ vidx ].v_win_fac3;
			double v_win = v_win1+v_win2+v_win3+2*planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].spacing;

			double vmin = 2*xpos*planes[uidx].costheta-upos-v_win;
			double vmax = vmin+2*v_win;
			//if(ic[index].EventID==0)printf("vmin = %1.6f, vmax = %1.6f, vwin = %1.6f, vwin1 = %1.6f, vwin2 = %1.6f, vwin3 = %1.6f\n", vmin, vmax, v_win, v_win1, v_win2, v_win3);
			for(int k = 0; k< nv; k++){
				double vpos = oc[index].hitpairs_v[k].second>=0 ? 0.5*(ic[index].AllHits[ oc[index].hitpairs_v[k].first ].pos+ic[index].AllHits[ oc[index].hitpairs_v[k].second ].pos): ic[index].AllHits[ oc[index].hitpairs_v[k].first ].pos;
				//if(ic[index].EventID<20)printf("evt %d: vmin = %1.6f <? vpos = %1.6f <? vmax = %1.6f\n", ic[index].EventID, vmin, vpos, vmax);
				if(vpos<vmin || vpos>vmax)continue;

				oc[index].AllTracklets[n_tkl].stationID = stID;
				nhits_tkl = 0;
				npts = 0;

				if(oc[index].hitpairs_x[i].first>=0){
					if(ic[index].EventID==0)printf("x first hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID, ic[index].AllHits[ oc[index].hitpairs_x[i].first ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_x[i].first ];
					oc[index].AllTracklets[n_tkl].nXHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;
					
					//fitarrays[index].x_array[npts] = ic[index].AllHits[ oc[index].hitpairs_x[i].first ].pos; // costheta = 1.
					//fitarrays[index].dx_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].spacing*0.5f; // costheta = 1.
					//fitarrays[index].y_array[npts] = 0.; // sintheta = 0.
					//fitarrays[index].dy_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].scaley*0.5f; // sintheta = 0.
					//fitarrays[index].z_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].z;
					
					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_x[i].first ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_x[i].first ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_x[i].first ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_x[i].first ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_x[i].second>=0){
					if(ic[index].EventID==0)printf("x second hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID, ic[index].AllHits[ oc[index].hitpairs_x[i].second ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_x[i].second ];
					oc[index].AllTracklets[n_tkl].nXHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;

					//fitarrays[index].x_array[npts] = ic[index].AllHits[ oc[index].hitpairs_x[i].second ].pos; // costheta = 1.
					//fitarrays[index].dx_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].spacing*0.5f; // costheta = 1.
					//fitarrays[index].y_array[npts] = 0.0; // sintheta = 0.
					//fitarrays[index].dy_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].scaley*0.5f; // sintheta = 0.
					//fitarrays[index].z_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].z;

					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_x[i].second ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_x[i].second ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_x[i].second ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_x[i].second ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_u[j].first>=0){
					if(ic[index].EventID==0)printf("u first hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID, ic[index].AllHits[ oc[index].hitpairs_u[j].first ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_u[j].first ];
					oc[index].AllTracklets[n_tkl].nUHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;
					
					//fitarrays[index].x_array[npts] = ic[index].AllHits[ oc[index].hitpairs_u[j].first ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].costheta;
					//fitarrays[index].dx_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].costheta);
					//fitarrays[index].y_array[npts] = ic[index].AllHits[ oc[index].hitpairs_u[j].first ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].sintheta;
					//fitarrays[index].dy_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].sintheta);
					//fitarrays[index].z_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].z;
					
					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_u[j].first ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_u[j].first ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[i].first ].detectorID-1 ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[i].first ].detectorID-1 ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_u[j].first ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_u[j].first ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_u[j].second>=0){
					if(ic[index].EventID==0)printf("u second hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID, ic[index].AllHits[ oc[index].hitpairs_u[j].second ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_u[j].second ];
					oc[index].AllTracklets[n_tkl].nUHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;

					//fitarrays[index].x_array[npts] = ic[index].AllHits[ oc[index].hitpairs_u[j].second ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].costheta;
					//fitarrays[index].dx_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].costheta);
					//fitarrays[index].y_array[npts] = ic[index].AllHits[ oc[index].hitpairs_u[j].second ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].sintheta;
					//fitarrays[index].dy_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].sintheta);
					//fitarrays[index].z_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].z;
					
					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_u[j].second ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_u[j].second ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[i].second ].detectorID-1 ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[i].second ].detectorID-1 ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_u[j].second ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_u[j].second ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_v[k].first>=0){
					if(ic[index].EventID==0)printf("v first hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID, ic[index].AllHits[ oc[index].hitpairs_v[k].first ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_v[k].first ];
					oc[index].AllTracklets[n_tkl].nVHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;
					
					//fitarrays[index].x_array[npts] = ic[index].AllHits[ oc[index].hitpairs_v[k].first ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].costheta;
					//fitarrays[index].dx_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].costheta);
					//fitarrays[index].y_array[npts] = ic[index].AllHits[ oc[index].hitpairs_v[k].first ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].sintheta;
					//fitarrays[index].dy_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].sintheta);
					//fitarrays[index].z_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].z;
					
					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_v[k].first ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_v[k].first ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[i].first ].detectorID-1 ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[i].first ].detectorID-1 ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_v[k].first ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_v[k].first ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_v[k].second>=0){
					if(ic[index].EventID==0)printf("v second hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID, ic[index].AllHits[ oc[index].hitpairs_v[k].second ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_v[k].second ];
					oc[index].AllTracklets[n_tkl].nVHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;
					
					//fitarrays[index].x_array[npts] = ic[index].AllHits[ oc[index].hitpairs_v[k].second ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].costheta;
					//fitarrays[index].dx_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].costheta);
					//fitarrays[index].y_array[npts] = ic[index].AllHits[ oc[index].hitpairs_v[k].second ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].sintheta;
					//fitarrays[index].dy_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].sintheta);
					//fitarrays[index].z_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].z;
					
					fitarrays[index].drift_dist[npts] = 0;//ic[index].AllHits[ oc[index].hitpairs_v[k].second ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_v[k].second ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[i].second ].detectorID-1 ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[i].second ].detectorID-1 ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_v[k].second ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_v[k].second ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(npts<4)continue;
				if(oc[index].AllTracklets[n_tkl].nXHits<1 || oc[index].AllTracklets[n_tkl].nUHits<1 || oc[index].AllTracklets[n_tkl].nVHits<1)continue;
				

				fitarrays[index].output_parameters[0] = 0.;
				fitarrays[index].output_parameters[1] = 0.;
				fitarrays[index].output_parameters[2] = 0.;
				fitarrays[index].output_parameters[3] = 0.;
				/*
				get_straighttrack_fixedpoint(npts, fitarrays[index].drift_dist, fitarrays[index].resolution, 
								fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z, 
								fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz,
								fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B,
								fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, 
								fixedpoint, fitarrays[index].chi2);

				fitarrays[index].lambda = 0.01f;
				fitarrays[index].chi2prev = 1.e20;
				
				gpufit_algorithm_fitter(fitparams[0].max_iterations,
										fitarrays[index].n_iter, fitarrays[index].iter_failed, fitarrays[index].finished,
										fitarrays[index].state, fitarrays[index].skip, fitarrays[index].singular,
										fitparams[0].tolerance, fitparams[0].nparam,
										fitparams[0].parameter_limits_min, fitparams[0].parameter_limits_max,  
										fitarrays[index].output_parameters, fitarrays[index].prev_parameters, 
										fitarrays[index].chi2, fitarrays[index].chi2prev,
										fitarrays[index].values, fitarrays[index].derivatives, fitarrays[index].gradients,
										fitarrays[index].hessians, fitarrays[index].scaling_vector, 
										fitarrays[index].deltas, fitarrays[index].lambda,
										fitarrays[index].Ainv,
										//fitarrays[index].calc_matrix, fitarrays[index].abs_row, fitarrays[index].abs_row_index,
										npts, fitarrays[index].drift_dist, fitarrays[index].resolution, 
										fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z, 
										fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz);					
				*/
				oc[index].AllTracklets[n_tkl].x0 = fitarrays[index].output_parameters[0];
				oc[index].AllTracklets[n_tkl].y0 = fitarrays[index].output_parameters[1];
				oc[index].AllTracklets[n_tkl].tx = fitarrays[index].output_parameters[2];
				oc[index].AllTracklets[n_tkl].ty = fitarrays[index].output_parameters[3];
				oc[index].AllTracklets[n_tkl].err_x0 = fitarrays[index].output_parameters_errors[0];
				oc[index].AllTracklets[n_tkl].err_y0 = fitarrays[index].output_parameters_errors[1];
				oc[index].AllTracklets[n_tkl].err_tx = fitarrays[index].output_parameters_errors[2];
				oc[index].AllTracklets[n_tkl].err_ty = fitarrays[index].output_parameters_errors[3];
				oc[index].AllTracklets[n_tkl].chisq = fitarrays[index].chi2;
				
				/*
				if(ic[index].EventID==0)printf("track: x0 = %1.6f +- %1.6f, y0 = %1.6f +- %1.6f, tx = %1.6f +- %1.6f, ty = %1.6f +- %1.6f; chi2 = %1.6f\n", 
					oc[index].AllTracklets[n_tkl].x0, oc[index].AllTracklets[n_tkl].err_x0, 
					oc[index].AllTracklets[n_tkl].y0, oc[index].AllTracklets[n_tkl].err_y0, 
					oc[index].AllTracklets[n_tkl].tx, oc[index].AllTracklets[n_tkl].err_tx, 
					oc[index].AllTracklets[n_tkl].ty, oc[index].AllTracklets[n_tkl].err_ty,
					oc[index].AllTracklets[n_tkl].chisq);
				*/

				if(!match_tracklet_to_hodo(oc[index].AllTracklets[n_tkl], stID, ic, planes))continue;

				
				if(n_tkl>TrackletSizeMax){
					continue;
					//printf("evt %d: n_tkl = %d > %d\n", oc[index].EventID, n_tkl, TrackletSizeMax);
				}
				n_tkl++;
				oc[index].nTKL_stID[stID]++;

			}
		}
	}
	
	//if(print)printf("evt %d number of tracklets %d\n", oc[index].EventID, n_tkl);
	//printf("%d ", n_tkl);
	oc[index].nTracklets = n_tkl;
}



//device int resolve_leftright(gTracklet tkl)
//{
//}


// kernel to combine tracklets into back partial tracks
__global__ void gkernel_BackPartialTracks(gEvent* ic, gSW* oc, gFitArrays* fitarrays, gPlane* planes, gFitParams* fitparams){
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	short stID = 5;
	//for the time being, declare in function;
	REAL x_pos[4];
	REAL z_pos[4];
	REAL a, b, sum, det, sx, sy, sxx, syy, sxy;
	short nprop, iprop;
	REAL xExp;
	short nhits_X2, nhits_X3;
	int n_tkl = oc[index].nTracklets;
	int nhits_2 = 0, nhits_3 = 0;
	
	REAL fixedpoint[3] = {0, 0, 0};

	//tracklets in station 3  (D3p, stID 4-1, D3m, stID 5-1)
	for(int i = oc[index].nTKL_stID[2]; i<oc[index].nTKL_stID[2]+oc[index].nTKL_stID[3]+oc[index].nTKL_stID[4]; i++){
		nhits_2 = oc[index].AllTracklets[i].nXHits+oc[index].AllTracklets[i].nUHits+oc[index].AllTracklets[i].nVHits;
		nhits_X2 = 0;
		for(int k = 0; k<nhits_2; k++){
			if(oc[index].AllTracklets[i].hits[k].detectorID==geometry::dets_x[0] || oc[index].AllTracklets[i].hits[k].detectorID==geometry::dets_x[1])
			{
				x_pos[nhits_X2] = oc[index].AllTracklets[i].hits[k].pos;
				z_pos[nhits_X2] = planes[oc[index].AllTracklets[i].hits[k].detectorID].z;
				nhits_X2++;
			}
		}
		
		//tracklets in station 2 (D2, stID 3-1)
		for(int j = 0; j<oc[index].nTKL_stID[2]; j++){
			if(fabs(oc[index].AllTracklets[i].tx - oc[index].AllTracklets[j].tx) > TX_MAX || fabs(oc[index].AllTracklets[i].ty - oc[index].AllTracklets[i].ty) > TY_MAX)continue;
			
			nhits_3 = oc[index].AllTracklets[j].nXHits+oc[index].AllTracklets[j].nUHits+oc[index].AllTracklets[j].nVHits;
			nhits_X3 = nhits_X2;
			for(int k = 0; k<nhits_3; k++){
				if(oc[index].AllTracklets[j].hits[k].detectorID==geometry::dets_x[2] || oc[index].AllTracklets[j].hits[k].detectorID==geometry::dets_x[3] ||
			   	oc[index].AllTracklets[j].hits[k].detectorID==geometry::dets_x[4] || oc[index].AllTracklets[j].hits[k].detectorID==geometry::dets_x[5])
				{
					x_pos[nhits_X3] = oc[index].AllTracklets[i].hits[k].pos;
					z_pos[nhits_X3] = planes[oc[index].AllTracklets[i].hits[k].detectorID].z;
					nhits_X3++;
				}
			}
			
			chi2_simplefit(nhits_X2+nhits_X3, z_pos, x_pos, a, b, sum, det, sx, sy, sxx, syy, sxy);
			if(fabs(a)>X0_MAX || fabs(b)>TX_MAX)continue;
			//prop matching
			nprop = 0;
			for(short ip = 0; ip<4; ip++){
				iprop = 48+ip;
				xExp = a*planes[iprop-1].z+b;
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
				
			//Add tracklets: first, get number of hits in each
			oc[index].AllTracklets[n_tkl].stationID = stID;
			oc[index].AllTracklets[n_tkl].nXHits = oc[index].AllTracklets[i].nXHits+oc[index].AllTracklets[j].nXHits;
			oc[index].AllTracklets[n_tkl].nUHits = oc[index].AllTracklets[i].nUHits+oc[index].AllTracklets[j].nUHits;
			oc[index].AllTracklets[n_tkl].nVHits = oc[index].AllTracklets[i].nVHits+oc[index].AllTracklets[j].nVHits;
			
			//add in the new tracklet the hits from the two existing tracklets;
			//also fill the fit arrays
			for(int k = 0; k<nhits_2; k++){
				oc[index].AllTracklets[n_tkl].hits[k]=oc[index].AllTracklets[i].hits[k];
				oc[index].AllTracklets[n_tkl].hitsign[k]=oc[index].AllTracklets[i].hitsign[k];
				
				fitarrays[index].drift_dist[k] = 0.;//oc[index].AllTracklets[n_tkl].hits[k].driftDistance;
				fitarrays[index].resolution[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].spacing/3.4641f;//.resolution;
				fitarrays[index].p1x[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].p1x_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].dp1x * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1);
				fitarrays[index].p1y[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].p1y_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].dp1y * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1);
				fitarrays[index].p1z[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].p1z_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].dp1z * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1);
				fitarrays[index].deltapx[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].deltapx;
				fitarrays[index].deltapy[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].deltapy;
				fitarrays[index].deltapz[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].deltapz;
			}
				
			for(int k = nhits_2; k<nhits_2+nhits_3; k++){
				oc[index].AllTracklets[n_tkl].hits[k]=oc[index].AllTracklets[j].hits[k-nhits_2];
				oc[index].AllTracklets[n_tkl].hitsign[k]=oc[index].AllTracklets[j].hitsign[k-nhits_2];
					fitarrays[index].drift_dist[k] = 0.;//oc[index].AllTracklets[n_tkl].hits[k].driftDistance;
				fitarrays[index].resolution[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].spacing/3.4641f;//.resolution;
				fitarrays[index].p1x[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].p1x_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].dp1x * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1) ;
				fitarrays[index].p1y[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].p1y_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].dp1y * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1) ;
				fitarrays[index].p1z[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].p1z_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].dp1z * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1) ;
				fitarrays[index].deltapx[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].deltapx;
				fitarrays[index].deltapy[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].deltapy;
				fitarrays[index].deltapz[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].deltapz;
			}
			
			//then fit:
				
			get_straighttrack_fixedpoint(nhits_2+nhits_3, fitarrays[index].drift_dist, fitarrays[index].resolution, 
							fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z, 
							fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz,
							fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B,
							fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, 
							fixedpoint, fitarrays[index].chi2);
			
			fitarrays[index].output_parameters[0] = b;
			fitarrays[index].output_parameters[2] = a;
			
			// do a rough estimation of the starting fit
			//fitarrays[index].output_parameters[0] = (oc[index].AllTracklets[i].x0+oc[index].AllTracklets[j].x0)*0.5f;
			//fitarrays[index].output_parameters[1] = (oc[index].AllTracklets[i].y0+oc[index].AllTracklets[j].y0)*0.5f;
			//fitarrays[index].output_parameters[2] = (oc[index].AllTracklets[i].tx+oc[index].AllTracklets[j].tx)*0.5f;
			//fitarrays[index].output_parameters[3] = (oc[index].AllTracklets[i].ty+oc[index].AllTracklets[j].ty)*0.5f;
			fitarrays[index].lambda = 0.01f;
			fitarrays[index].chi2prev = 1.e20;
			
			//include fit here:
			/*
			gpufit_algorithm_fitter(fitparams[0].max_iterations,
								fitarrays[index].n_iter, fitarrays[index].iter_failed, fitarrays[index].finished,
								fitarrays[index].state, fitarrays[index].skip, fitarrays[index].singular,
								fitparams[0].tolerance, fitparams[0].nparam,
								fitparams[0].parameter_limits_min, fitparams[0].parameter_limits_max,  
								fitarrays[index].output_parameters, fitarrays[index].prev_parameters, 
								fitarrays[index].chi2, fitarrays[index].chi2prev,
								fitarrays[index].values, fitarrays[index].derivatives, fitarrays[index].gradients,
								fitarrays[index].hessians, fitarrays[index].scaling_vector, 
								fitarrays[index].deltas, fitarrays[index].lambda,
								fitarrays[index].Ainv, 
								nhits_2+nhits_3, fitarrays[index].drift_dist, fitarrays[index].resolution, 
								fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z, 
								fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz);					
			
			*/
			
			if(fitarrays[index].chi2>9000)continue;
			
			oc[index].AllTracklets[n_tkl].x0 = fitarrays[index].output_parameters[0];
			oc[index].AllTracklets[n_tkl].y0 = fitarrays[index].output_parameters[1];
			oc[index].AllTracklets[n_tkl].tx = fitarrays[index].output_parameters[2];
			oc[index].AllTracklets[n_tkl].ty = fitarrays[index].output_parameters[3];
			oc[index].AllTracklets[n_tkl].err_x0 = fitarrays[index].output_parameters_errors[0];
			oc[index].AllTracklets[n_tkl].err_y0 = fitarrays[index].output_parameters_errors[1];
			oc[index].AllTracklets[n_tkl].err_tx = fitarrays[index].output_parameters_errors[2];
			oc[index].AllTracklets[n_tkl].err_ty = fitarrays[index].output_parameters_errors[3];
			oc[index].AllTracklets[n_tkl].chisq = fitarrays[index].chi2;
			
			n_tkl++;
			oc[index].nTKL_stID[stID]++;
		}
	}
	
	oc[index].nTracklets = n_tkl;
	
}



__global__ void gKernel_XZ_YZ_tracking(gEvent* ic, gOutputEvent* oc, gStraightTrackBuilder* straighttrackbuilder, gStraightFitArrays* fitarrays, gPlane* planes)
{
        int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	oc[index].EventID = ic[index].EventID;
	oc[index].nAH = ic[index].nAH;
	
	//for the time being, declare in function;
	//REAL x_pos[4];
	//REAL x_res[4];
	//REAL z_pos[4];
	REAL a, b, sum, det, sx, sy, sxx, syy, sxy;
	short nprop, iprop;
	REAL xExp;
	short nhits_X2, nhits_X3;

        // 1- get X pairs in st2, st3:
	//D2: stid = 3-1
	int nx2 = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_x2, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, 2, 0, planes);
	int nu2 = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_u2, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, 2, 1, planes);
	int nv2 = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_v2, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, 2, 2, planes);

	//D3p: stid = 4-1
	int nx3p = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_x3p, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, 3, 0, planes);
	int nu3p = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_u3p, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, 3, 1, planes);
	int nv3p = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_v3p, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, 3, 2, planes);
	
	//D3p: stid = 5-1
	int nx3m = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_x3m, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, 4, 0, planes);
	int nu3m = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_u3m, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, 4, 1, planes);
	int nv3m = make_hitpairs_in_station(ic, straighttrackbuilder[index].hitpairs_v3m, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, 4, 2, planes);
	
	//if(ic[index].EventID==13)printf("evt %d, %d pairs in st2, %d pairs in st3+, %d pairs in st3-\n", ic[index].EventID, nx2, nx3p, nx3m);
	
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
			fitarrays[index].dx_array[nhits_X2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x2[i].first].detectorID-1].spacing/3.4641;
			fitarrays[index].z_array[nhits_X2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x2[i].first].detectorID-1].z;
			nhits_X2++;
		}
		if(straighttrackbuilder[index].hitpairs_x2[i].second>=0){
			fitarrays[index].x_array[nhits_X2] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x2[i].second].pos;
			fitarrays[index].dx_array[nhits_X2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x2[i].second].detectorID-1].spacing/3.4641;
			fitarrays[index].z_array[nhits_X2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x2[i].second].detectorID-1].z;
			nhits_X2++;
		}
				
		for(int j = 0; j<nx3p+nx3m; j++){
			nhits_X3 = nhits_X2;
			if(j<nx3p){
				if(straighttrackbuilder[index].hitpairs_x3p[i].first>=0){
					fitarrays[index].x_array[nhits_X3] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3p[i].first].pos;
					fitarrays[index].dx_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3p[i].first].detectorID-1].spacing/3.4641;
					fitarrays[index].z_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3p[i].first].detectorID-1].z;
					nhits_X3++;
				}
				if(straighttrackbuilder[index].hitpairs_x3p[i].second>=0){
					fitarrays[index].x_array[nhits_X3] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3p[i].second].pos;
					fitarrays[index].dx_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3p[i].second].detectorID-1].spacing/3.4641;
					fitarrays[index].z_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3p[i].second].detectorID-1].z;
					nhits_X3++;
				}
			}else{
				if(straighttrackbuilder[index].hitpairs_x3m[i].first>=0){
					fitarrays[index].x_array[nhits_X3] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[i].first].pos;
					fitarrays[index].dx_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[i].first].detectorID-1].spacing/3.4641;
					fitarrays[index].z_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[i].first].detectorID-1].z;
					nhits_X3++;
				}
				if(straighttrackbuilder[index].hitpairs_x3m[i].second>=0){
					fitarrays[index].x_array[nhits_X3] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[i].second].pos;
					fitarrays[index].dx_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[i].second].detectorID-1].spacing/3.4641;
					fitarrays[index].z_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[i].second].detectorID-1].z;
					nhits_X3++;
				}
			}			
			//don't we want something with at least the errors???
			
			//if(ic[index].EventID==13)for(int k = 0; k<nhits_X3; k++)printf("k %d, x %1.6f, y %1.6f +- %1.6f \n", k, fitarrays[index].z_array[k], fitarrays[index].x_array[k], fitarrays[index].dx_array[k]);
			//chi2_simplefit(nhits_X3, fitarrays[index].z_array, fitarrays[index].x_array, a, b, sum, det, sx, sy, sxx, syy, sxy);
			//if(ic[index].EventID==13)printf("evt %d: %1.6f %1.6f\n", ic[index].EventID, a, b); 
			
			fit_2D_track(nhits_X3, fitarrays[index].x_array, fitarrays[index].z_array, fitarrays[index].dx_array, fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B, fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, fitarrays[index].chi2_xz);
			//if(0<=ic[index].EventID && ic[index].EventID<=20)printf("evt %d: (%1.6f+-%1.6f) + z * (%1.6f+-%1.6f), %1.6f \n", ic[index].EventID, 
			//				fitarrays[index].output_parameters[0], fitarrays[index].output_parameters_errors[0], 
			//				fitarrays[index].output_parameters[1], fitarrays[index].output_parameters_errors[1], 
			//				fitarrays[index].chi2_xz); 
			
			if(fabs(fitarrays[index].output_parameters[0])>X0_MAX || fabs(fitarrays[index].output_parameters[1])>TX_MAX)continue;
			
			//prop matching
			nprop = 0;
			for(short ip = 0; ip<4; ip++){
				iprop = 48+ip;
				xExp = a*planes[iprop-1].z+b;
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

			
		}// loops on st3 X hits
	
	}//loop on st2 X hits
	
        // 3- get U, V pairs compatible with the X pair in each station: already implemented, tested;


        // 4- U, V pairs also have to be compatible with the X pair slope => this last point needs more thought, and will probably need adjustments (2-3 days to implement and test);


        // 5- calculate Y from U, V hits (0.5-1 day to implement and test);


        // 6- fit the X-Y slope (0.5-1 day to implement and test);


        // 7- combine the XZ and XY track and combine hodoscope matching (already implemented, needs test: 1-1.5 day)

}


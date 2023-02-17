#include "reconstruction_helper.cuh"
//#include "tracknumericalminimizer.cuh"
#include "trackanalyticalminimizer.cuh"

// kernel functions: 
// CUDA C++ extends C++ by allowing the programmer to define C++ functions, called kernels, that, when called, 
// are executed N times in parallel by N different CUDA threads, as opposed to only once like regular C++ functions. 


// event reducer: 
__global__ void gkernel_eR(gEvent* ic) {
	//printf("Running the kernel function...\n");
	// retrieve global thread index
	const int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	double w_max; // max drift distance of the hit furthest from the cluster avg position // current average position of cluster * 0.9
	
	double w_min; // max drift distance of the hit closest to the cluster avg position // current average position of cluster * 0.4
	double dt_mean; // tbd
	int cluster_iAH_arr_cur; // current cluster array
	int cluster_iAH_arr_size; // cluster size i.e number of hits in cluster
	static int cluster_iAH_arr[ClusterSizeMax]; // global cluster array 
	int uniqueID; // hit unique element ID
	int uniqueID_curr; // current hit unique element ID
	double tdcTime_curr; // current hit TDC time
	int iAH; // hit index 
	int nAH_reduced; // number of hits after hit quality filtering
	// int nHitsPerDetector[nDetectors+1];
	
	const int idxoff_global = index*EstnAHMax;
	int iAH_global;
	//if(ic[index].EventID==0){
	//	printf("evt = %d, nAH = %d: \n", ic[index].EventID, ic[index].nAH);
	//	for(int i = 1; i<=nDetectors; i++)printf(" det %d: %d;  ", i, ic[index].NHits[i]);
	//	printf("\n");
	//}

	// initialization of array size
	cluster_iAH_arr_size = 0;
	nAH_reduced = 0;
	
	// event reducing/hit filtering
	for(iAH = 0; iAH<ic->nAH[index]; ++iAH) {
		iAH_global = idxoff_global+iAH;
		// if hit not good, set its detID to 0 and continue;
		if((ic->AllHits[iAH_global].flag & hitFlagBit(1)) == 0) {
			//if(ic[index].EventID==0)printf("hit det %d Skip out-of-time...\n", ic->AllHits[iAH_global].detectorID);
			ic->AllHits[iAH_global].detectorID = 0;
			continue;
		}
		uniqueID = uniqueID_curr = -1;

		// hits in DCs or Prop tubes
		if(ic->AllHits[iAH_global].detectorID < 31 || ic->AllHits[iAH_global].detectorID > 46) {
			// evaluate "unique ID"
			uniqueID = ic->AllHits[iAH_global].detectorID*1000 + ic->AllHits[iAH_global].elementID;
			// compare with current unique element ID; if different, update the unique element ID and time info 
			if(uniqueID != uniqueID_curr) {
				uniqueID_curr = uniqueID;
				tdcTime_curr = ic->AllHits[iAH_global].tdcTime;
			}
			// if next hit and current hit belong to the same element: if detID>36 => prop tubes (reminder that hodoscpes are out of the picture in this scope), 
			// we're suppose to have one signal (a second signal would be after-pulsing)
			// if time difference between new hit and current hit is less than 80ns for DCs, it's also considered after-pulsing
			else {
				if(ic->AllHits[iAH_global].detectorID > 36 || ((ic->AllHits[iAH_global].tdcTime - tdcTime_curr >= 0.0) && (ic->AllHits[iAH_global].tdcTime - tdcTime_curr < 80.0)) || ((ic->AllHits[iAH_global].tdcTime - tdcTime_curr <= 0.0) && (ic->AllHits[iAH_global].tdcTime - tdcTime_curr > -80.0))) {
					//if(ic[index].EventID==0)printf("hit det %d el %d Skip after-pulse...\n", ic->AllHits[iAH_global].detectorID, ic->AllHits[iAH_global].elementID);
					ic->AllHits[iAH_global].detectorID = 0;
					continue;
				}
				else {
					tdcTime_curr = ic->AllHits[iAH_global].tdcTime;
				}
			}
		}
		// declustering of hits in DCs (from CPU code, I understand this one better)
		// if there are hits in the same plane and hitting to neighboring wires, they both give redundant information: 
		if(ic->AllHits[iAH_global].detectorID <= nChamberPlanes) {
			//if(ic[index].EventID==0)printf("%d\n", cluster_iAH_arr_size);
//			printf("Decluster...\n");
			if(cluster_iAH_arr_size == ClusterSizeMax) {
//				printf("Oversized cluster...\n");
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
//printf("hit det/elem: %d, %d; %d, %d\n", ic->AllHits[iAH_global].detectorID, ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID, 
//	    	      	      	  	 ic->AllHits[iAH_global].elementID, ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].elementID);
//printf("diffs: %d, %d\n", (ic->AllHits[iAH_global].detectorID - ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID),
//	       	   	  (ic->AllHits[iAH_global].elementID - ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].elementID));
//printf("bools: %d, %d\n", (ic->AllHits[iAH_global].detectorID != ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID),
//	       	   	  (abs(ic->AllHits[iAH_global].elementID - ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].elementID) > 1));
	    	     	   	//}
				if((ic->AllHits[iAH_global].detectorID != ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID) || (abs(ic->AllHits[iAH_global].elementID - ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].elementID) > 1)) {
					// if 2 hits in cluster, evaluate w_max and w_min; drift distance has to be < w_min for one of the hits, while it has to be < w_max for the other hit 
					if(cluster_iAH_arr_size == 2) {
						w_max = 0.9*0.5*(ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].pos - ic->AllHits[idxoff_global+cluster_iAH_arr[0]].pos);
						w_min = 4.0/9.0*w_max;
						if((ic->AllHits[idxoff_global+cluster_iAH_arr[0]].driftDistance > w_max && ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].driftDistance > w_min) || (ic->AllHits[idxoff_global+cluster_iAH_arr[0]].driftDistance > w_min && ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].driftDistance > w_max)) {
						//if(ic[index].EventID==0)printf("hit indices: %d %d %d\n", iAH, cluster_iAH_arr[cluster_iAH_arr_size-1], cluster_iAH_arr[0]);
							//eliminating the existing hit with the lagest drift distance
							if(ic->AllHits[idxoff_global+cluster_iAH_arr[0]].driftDistance > ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].driftDistance) {
								//if(ic[index].EventID==0)printf("1 - hit det %d elem %d Skip cluster...\n", ic->AllHits[idxoff_global+cluster_iAH_arr[0]].detectorID, ic->AllHits[idxoff_global+cluster_iAH_arr[0]].elementID);
								ic->AllHits[idxoff_global+cluster_iAH_arr[0]].detectorID = 0;
							}
							else {
								//if(ic[index].EventID==0)printf("2 - hit det %d elem %d Skip cluster...\n", ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID, ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].elementID);
								ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID = 0;
							}
						}
						// if the time difference is less than 8 ns for detectors 19 to 24 (which btw are DC3p) we remove both
						else if((((ic->AllHits[idxoff_global+cluster_iAH_arr[0]].tdcTime - ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].tdcTime) >= 0.0 && (ic->AllHits[idxoff_global+cluster_iAH_arr[0]].tdcTime - ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].tdcTime) < 8.0) || ((ic->AllHits[idxoff_global+cluster_iAH_arr[0]].tdcTime - ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].tdcTime) <= 0.0 && (ic->AllHits[idxoff_global+cluster_iAH_arr[0]].tdcTime - ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].tdcTime) > -8.0)) && (ic->AllHits[idxoff_global+cluster_iAH_arr[0]].detectorID >= 19 && ic->AllHits[idxoff_global+cluster_iAH_arr[0]].detectorID <= 24)) {
						        //if(ic[index].EventID==0)printf("3 - hit det %d elem %d Skip cluster...\n", ic->AllHits[iAH_global].detectorID, ic->AllHits[iAH_global].elementID);
							ic->AllHits[idxoff_global+cluster_iAH_arr[0]].detectorID = 0;
							ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID = 0;
						}
					}
					// if 3 hits or more in cluster: we essentially discard them all;
					if(cluster_iAH_arr_size >= 3) {
						// evaluate the mean time difference;
						dt_mean = 0.0;
						for(cluster_iAH_arr_cur = 1; cluster_iAH_arr_cur < cluster_iAH_arr_size; ++cluster_iAH_arr_cur) {
							dt_mean += ((ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_cur]].tdcTime - ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_cur-1]].tdcTime) > 0.0 ? (ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_cur]].tdcTime - ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_cur-1]].tdcTime) : (ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_cur-1]].tdcTime - ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_cur]].tdcTime));
						}
						dt_mean = dt_mean/(cluster_iAH_arr_size - 1);
						// if mean time difference is less than 10, that's electronic noise, so we remove them all.
						if(dt_mean < 10.0) {
						        //if(ic[index].EventID==0)printf("4 - hit det %d elem %d Skip cluster...\n", ic->AllHits[iAH_global].detectorID, ic->AllHits[iAH_global].elementID);
							for(cluster_iAH_arr_cur = 0; cluster_iAH_arr_cur < cluster_iAH_arr_size; ++cluster_iAH_arr_cur) {
								ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_cur]].detectorID = 0;
							}
						}
						// otherwise, we remove them all except first and last
						else {
						        //if(ic[index].EventID==0)printf("5 - hit det %d elem %d Skip cluster...\n", ic->AllHits[iAH_global].detectorID, ic->AllHits[iAH_global].elementID);
							for(cluster_iAH_arr_cur = 1; cluster_iAH_arr_cur < cluster_iAH_arr_size; ++cluster_iAH_arr_cur) {
								ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_cur]].detectorID = 0;
							}
						}
					}
					cluster_iAH_arr_size = 0;
				} else { // if hits are in the same detector and in two neighboring cells!
				        // current hit and previous hit are in same detector plane and in neighbor wires: 
				  	// we count how many hits we have in this case, until we find a hit in a different detector or in a wire that is not a neighbor to the previous hit.
				  	//if(ic[index].EventID==0)printf("h1: det %d el %d; h2 det %d el %d \n", ic->AllHits[iAH_global].detectorID, ic->AllHits[iAH_global].elementID, ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].detectorID, ic->AllHits[idxoff_global+cluster_iAH_arr[cluster_iAH_arr_size-1]].elementID);
				  	cluster_iAH_arr[cluster_iAH_arr_size] = iAH;
				  	++cluster_iAH_arr_size;
				}
			}
		}
	}
	//end of the hit loop

	// Hit reduction: 
	// store in "AllHits" containers only hits with non-zero detectorID and couting those with nAH_reduced
	for(iAH = 0; iAH<ic->nAH[index]; ++iAH) {
		iAH_global = idxoff_global+iAH;
		if(ic->AllHits[iAH_global].detectorID != 0) {
			ic->AllHits[idxoff_global+nAH_reduced] = ic->AllHits[iAH_global];
			++nAH_reduced;
		}
	}

	// compute hits per detector
	int nEventHits = nAH_reduced;
	// reinitialize number of hits per detector
	for(auto iDetector = 1; iDetector <= nDetectors; ++iDetector) {
		ic->NHits[index*nDetectors+iDetector] = 0;
	}
	// loop on reduced hits and counting number of hits per detector
	for(auto iHit = 0; iHit < nEventHits; ++iHit) {
		auto detectorId = ic->AllHits[iHit].detectorID;
		if(detectorId != 0) {
			++ic->NHits[index*nDetectors+detectorId];
		}
	}

	ic->nAH[index] = nAH_reduced;
	ic->HasTooManyHits[index] = false;

	int hitmult_[5] = {0, 0, 0, 0, 0};
	for(auto iDetector = 1; iDetector <= 12; ++iDetector)hitmult_[0]+= ic->NHits[index*nDetectors+iDetector];
	if(hitmult_[0]>250)ic->HasTooManyHits[index] = true;
	for(auto iDetector = 13; iDetector <= 18; ++iDetector)hitmult_[1]+= ic->NHits[index*nDetectors+iDetector];
	if(hitmult_[1]>200)ic->HasTooManyHits[index] = true;
	for(auto iDetector = 19; iDetector <= 24; ++iDetector)hitmult_[2]+= ic->NHits[index*nDetectors+iDetector];
	if(hitmult_[2]>150)ic->HasTooManyHits[index] = true;
	for(auto iDetector = 25; iDetector <= 30; ++iDetector)hitmult_[3]+= ic->NHits[index*nDetectors+iDetector];
	if(hitmult_[3]>120)ic->HasTooManyHits[index] = true;
	for(auto iDetector = 47; iDetector <= 54; ++iDetector)hitmult_[4]+= ic->NHits[index*nDetectors+iDetector];
	if(hitmult_[4]>250)ic->HasTooManyHits[index] = true;
}


// function to match a tracklet to a hodoscope hit
__device__ int match_tracklet_to_hodo(const gTracklet tkl, const int stID, gEvent* ic, const gPlane* planes)
{
	const int index = threadIdx.x + blockIdx.x * blockDim.x;
	const int idxoff_global = index*EstnAHMax;
	int masked = 0;//false
	// first, define the search region, and foremost, the planes on which we define this search region, which depends on the station ID we're looking at
	// define the region in which we are supposed to have hits:

	//printf(" stID %d hodo plane[0] %d [1] %d \n", stID, geometry::hodoplanerange[stID][0], geometry::hodoplanerange[stID][1]);
	//printf(" x0 %1.4f +- %1.4f, y0 %1.4f +- %1.4f, tx %1.4f +- %1.4f, ty %1.4f +- %1.4f \n", tkl.x0, tkl.err_x0, tkl.y0, tkl.err_y0, tkl.tx, tkl.err_tx, tkl.ty, tkl.err_ty);		
	
	REAL xhodo, yhodo, err_x, err_y, xmin, xmax, ymin, ymax;
		
	// loop on the hits and select hodoscope hits corresponding to the station
	for(int i = 0; i<ic->nAH[index]; i++){
		//we only consider hits in the hodoscopes planes corresponding to the station where the tracklet is reconstructed 
		if(geometry::hodoplanerange[stID][0]>ic->AllHits[idxoff_global+i].detectorID || geometry::hodoplanerange[stID][1]<ic->AllHits[idxoff_global+i].detectorID)continue;
		
		//calculate the track position at the hodoscope plane z
		xhodo = planes->z[ic->AllHits[idxoff_global+i].detectorID]*tkl.tx+tkl.x0;
		yhodo = planes->z[ic->AllHits[idxoff_global+i].detectorID]*tkl.ty+tkl.y0;
		
		err_x = 3.f*(fabs(planes->z[ic->AllHits[idxoff_global+i].detectorID]*tkl.err_tx)+tkl.err_x0);
		err_y = 3.f*(fabs(planes->z[ic->AllHits[idxoff_global+i].detectorID]*tkl.err_ty)+tkl.err_y0);

		//printf(" det %d elem %d z_hodo %1.4f x_hodo %1.4f y_hodo %1.4f err x %1.4f err y %1.4f \n", ic->AllHits[idxoff_global+i].detectorID, ic->AllHits[idxoff_global+i].elementID, planes[ic->AllHits[idxoff_global+i].detectorID].z, xhodo, yhodo, err_x, err_y);

		//calculate "xmin, xmax, ymin, ymax" in which the track is supposed to pass through; 
		//these are basicially defined as the spatial coverage of the hit hodoscope element (plus some fudge factor for x)
		if(planes->costheta[ic->AllHits[idxoff_global+i].detectorID]>0.99){
			xmin = ic->AllHits[idxoff_global+i].pos-planes->cellwidth[ic->AllHits[idxoff_global+i].detectorID]*0.5f;
			xmax = ic->AllHits[idxoff_global+i].pos+planes->cellwidth[ic->AllHits[idxoff_global+i].detectorID]*0.5f;
			
			ymin = planes->y1[ic->AllHits[idxoff_global+i].detectorID];
			ymax = planes->y2[ic->AllHits[idxoff_global+i].detectorID];
			
			xmin-=(xmax-xmin)*geometry::hodofudgefac[stID];
			xmax+=(xmax-xmin)*geometry::hodofudgefac[stID];
			
			ymin-=(ymax-ymin)*geometry::hodofudgefac[stID];
			ymax+=(ymax-ymin)*geometry::hodofudgefac[stID];
			
			ymin+=planes->y0[ic->AllHits[idxoff_global+i].detectorID];
			ymax+=planes->y0[ic->AllHits[idxoff_global+i].detectorID];
		}else{
			xmin = planes->x1[ic->AllHits[idxoff_global+i].detectorID];
			xmax = planes->x2[ic->AllHits[idxoff_global+i].detectorID];
			
			ymin = ic->AllHits[idxoff_global+i].pos-planes->cellwidth[ic->AllHits[idxoff_global+i].detectorID]*0.5f;
			ymax = ic->AllHits[idxoff_global+i].pos+planes->cellwidth[ic->AllHits[idxoff_global+i].detectorID]*0.5f;
			
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


__global__ void gKernel_XZ_YZ_tracking_new(gEvent* ic, gOutputEvent* oc, gStraightTrackBuilder* straighttrackbuilder, gStraightFitArrays* fitarrays, const gPlane* planes, bool best_candyz_only = true)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	//printf("index %d, thread %d, block %d \n", index, threadIdx.x, blockIdx.x);

	oc->EventID[index] = ic->EventID[index];
	oc->nAH[index] = ic->nAH[index];
	oc->HasTooManyHits[index] = ic->HasTooManyHits[index];
	oc->nTracklets[index] = 0;
	
	if(oc->HasTooManyHits[index])return;
	
	const int idxoff_global = index*EstnAHMax;
	const int idxoff_tkl = index*TrackletSizeMax;
	
	const short nbins_st2 = geometry::N_WCHitsBins[1];
	const short nbins_st3 = geometry::N_WCHitsBins[2]+geometry::N_WCHitsBins[3];
	
	const int st2x_offset = nbins_st2*geometry::MaxHitsProj[0]*index;
	const int st2uv_offset = nbins_st2*geometry::MaxHitsProj[1]*index;
	
	const int st3x_offset = nbins_st3*geometry::MaxHitsProj[0]*index;
	const int st3uv_offset = nbins_st3*geometry::MaxHitsProj[1]*index;
	
	short stid, projid;
	
	unsigned int hitidx1[100];
	unsigned int hitidx2[100];
	short hitflag1[100];
	short hitflag2[100];

	int nhitpairs_x2[28];
	int nhitpairs_u2[28];
	int nhitpairs_v2[28];
	
	int nhitpairs_x3[58];
	int nhitpairs_u3[58];
	int nhitpairs_v3[58];
	
	float A_[4];
	float Ainv_[4];
	float B_[2];
	float Par[4];
	float ParErr[4];
	float chi2;
	
        // 1- get hit pairs in bins in st2, st3:
	// the pairs for st3p and st3m are now in the same array. to allow for parallelization of d3p/d3m if needed;
	projid = 0;
	//D2: stid = 3-1
	stid = 2;
	//make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_x2, straighttrackbuilder->nhitpairs_x2, straighttrackbuilder->hitidx1, straighttrackbuilder->hitidx2, straighttrackbuilder->hitflag1, straighttrackbuilder->hitflag2, stid, projid);
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_x2, nhitpairs_x2, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);

	//D3p: stid = 4-1
	stid = 3;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_x3, nhitpairs_x3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);

	//D3p: stid = 5-1
	stid = 4;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_x3, nhitpairs_x3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
		
#ifdef DEBUG
	for(int bin = 0; bin<58; bin++){
	if(nhitpairs_x3[bin])printf("evt %d bin %d nhits x3 = %d\n", ic[index].EventID, bin, nhitpairs_x3[bin]);
	for(int i = 0; i<nhitpairs_x3[bin]; i++){
		if(hitpairs_x3[bin+28*i].first>-1)printf(" evt %d bin %d  i %d detid1 %d elid1 %d\n", ic[index].EventID, bin, i, 
				ic->AllHits[idxoff_global+hitpairs_x3[bin+29*i].first].detectorID, ic->AllHits[idxoff_global+hitpairs_x3[bin+29*i].first].elementID);
		if(hitpairs_x3[bin+28*i].second>-1)printf(" evt %d bin %d i %d detid2 %d elid2 %d\n", ic[index].EventID, bin, i, 
				ic->AllHits[idxoff_global+hitpairs_x3[bin+29*i].second].detectorID, ic->AllHits[idxoff_global+hitpairs_x3[bin+29*i].second].elementID);
	}
	}
#endif
	
	projid = 1;
	stid = 2;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_u2, nhitpairs_u2, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
	
	stid = 3;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_u3, nhitpairs_u3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
	
	stid = 4;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_u3, nhitpairs_u3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
		
	projid = 2;
	stid = 2;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_v2, nhitpairs_v2, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
		
	stid = 3;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_v3, nhitpairs_v3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
		
	stid = 4;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_v3, nhitpairs_v3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
	
	bool bin_overflows = false;
	for(int bin = 0; bin<28; bin++){
		if(nhitpairs_x2[bin]>geometry::MaxHitsProj[0])bin_overflows = true;
		//if(nhitpairs_x2[bin])printf("evt %d bin %d nhits x2 = %d\n", ic[index].EventID, bin, nhitpairs_x2[bin]);
		if(nhitpairs_u2[bin]>geometry::MaxHitsProj[1])bin_overflows = true;
		//if(nhitpairs_u2[bin])printf("evt %d bin %d nhits u2 = %d\n", ic[index].EventID, bin, nhitpairs_u2[bin]);
		if(nhitpairs_v2[bin]>geometry::MaxHitsProj[2])bin_overflows = true;
		//if(nhitpairs_v2[bin])printf("evt %d bin %d nhits v2 = %d\n", ic[index].EventID, bin, nhitpairs_v2[bin]);
	}
	for(int bin = 0; bin<58; bin++){
		if(nhitpairs_x3[bin]>geometry::MaxHitsProj[0])bin_overflows = true;
		//if(nhitpairs_x3[bin])printf("evt %d bin %d nhits x3 = %d\n", ic[index].EventID, bin, nhitpairs_x3[bin]);
		if(nhitpairs_u3[bin]>geometry::MaxHitsProj[1])bin_overflows = true;
		//if(nhitpairs_u3[bin])printf("evt %d bin %d nhits u3 = %d\n", ic[index].EventID, bin, nhitpairs_u3[bin]);
		if(nhitpairs_v3[bin]>geometry::MaxHitsProj[2])bin_overflows = true;
		//if(nhitpairs_v3[bin])printf("evt %d bin %d nhits v3 = %d\n", ic[index].EventID, bin, nhitpairs_v3[bin]);
	}
	if(bin_overflows){
		printf("Evt %d discarded, too many hits\n", ic[index].EventID);
		return;
	}
	
	// declare variables useful for the loop 
	int ntkl = 0;
	
	//for the time being, declare in function;
	short nprop, iprop, idet;
	REAL xExp, ipos;
	
	short nhits_x;
	short nhits_uv;
	short nhits_v;
	
	short nhits_x2, nhits_x3;
	short nhits_u2, nhits_u3;
	short nhits_v2, nhits_v3;
	
	// if parallel:
	// const short nbins_st3 = geometry::N_WCHitsBins[2];
	const short nbins_total = nbins_st2*nbins_st3;
	short bin2, bin3;
	
	int nx2, nu2, nv2;
	int nx3, nu3, nv3;
	
	int ncomb_x, ncomb_uv;
	
	short i_x2, i_x3;
	short i_u2, i_u3, i_v2, i_v3;
	
	int i_x, i_uv, i_hit;

	float y, err_y;
	float ty;
	
	float chi2min = 10000.1f;
	
	int n_goodxz;
	//int hitidx;
	//gHit hit;
	gTrack2D trackXZ, trackYZ, besttrackYZ;
	
	//loop on bins FIRST
	for(short i = 0; i<nbins_total; i++){
		bin2 = i%nbins_st2;
		bin3 = (i-bin2)/nbins_st2;
		//printf("bin %d %d\n", bin2, bin3);
		// if parallel:
		// bin3+=nbins_st3*i_thread;

		nx2 = nhitpairs_x2[bin2];
		nx3 = nhitpairs_x3[bin3];
		
		//printf("nx %d %d \n", nx2, nx3);
		if(nx2 == 0 || nx3==0) continue;
		
		nu2 = nhitpairs_u2[bin2];
		nu3 = nhitpairs_u3[bin3];

		//printf("nu %d %d \n", nu2, nu3);
		if(nu2 == 0 || nu3==0) continue;
		
		nv2 = nhitpairs_v2[bin2];
		nv3 = nhitpairs_v3[bin3];

		//printf("nv %d %d \n", nv2, nv3);
		if(nv2 == 0 || nv3==0) continue;
		
		// evaluating the number of combinations;
		ncomb_x = nx2*nx3;
		ncomb_uv = nu2*nu3*nv2*nv3;
		
		n_goodxz = 0;
				
		for(i_x = 0; i_x<ncomb_x; i_x++){
			i_x2 = i_x%nx3;
			i_x3 = (i_x-i_x2)/nx3;
			
			nhits_v = -1;
			nhits_x = 0;

			if(straighttrackbuilder->hitpairs_x2[st2x_offset+bin2+nbins_st2*i_x2].first>=0){
				trackXZ.hitlist[nhits_x] = straighttrackbuilder->hitpairs_x2[st2x_offset+bin2+nbins_st2*i_x2].first;
				FillFitArrays_X(nhits_x, ic->AllHits[trackXZ.hitlist[nhits_x]], 0, fitarrays, planes);
				nhits_x++;
			}
			if(straighttrackbuilder->hitpairs_x2[st2x_offset+bin2+nbins_st2*i_x2].second>=0){
				trackXZ.hitlist[nhits_x] = straighttrackbuilder->hitpairs_x2[st2x_offset+bin2+nbins_st2*i_x2].second;
				FillFitArrays_X(nhits_x, ic->AllHits[trackXZ.hitlist[nhits_x]], 0, fitarrays, planes);
				nhits_x++;
			}
			
			nhits_x2 = nhits_x;
			if(nhits_x2==0) continue;
			
			if(straighttrackbuilder->hitpairs_x3[st3x_offset+bin3+nbins_st3*i_x3].first>=0){
				trackXZ.hitlist[nhits_x] = straighttrackbuilder->hitpairs_x3[st3x_offset+bin3+nbins_st3*i_x3].first;
				FillFitArrays_X(nhits_x, ic->AllHits[trackXZ.hitlist[nhits_x]], 0, fitarrays, planes);
				nhits_x++;
			}
			if(straighttrackbuilder->hitpairs_x3[st3x_offset+bin3+nbins_st3*i_x3].second>=0){
				trackXZ.hitlist[nhits_x] = straighttrackbuilder->hitpairs_x3[st3x_offset+bin3+nbins_st3*i_x3].second;
				FillFitArrays_X(nhits_x, ic->AllHits[trackXZ.hitlist[nhits_x]], 0, fitarrays, planes);
				nhits_x++;
			}
			
			nhits_x3 = nhits_x-nhits_x2;
			if(nhits_x3==0) continue;
			
			fit_2D_track(nhits_x, fitarrays->x_array, fitarrays->z_array, fitarrays->dx_array, A_, Ainv_, B_, Par, ParErr, chi2);
			if(fabs(Par[0])>1.05*X0_MAX || fabs(Par[1])>1.05*TX_MAX)continue;
			
			trackXZ.x_0 = Par[0];
			trackXZ.err_x_0 = ParErr[0];
			trackXZ.tx_ = Par[1];
			trackXZ.err_tx_ = ParErr[1];
			trackXZ.nhits = nhits_x;
			
			n_goodxz++;
			
			//prop matching
			nprop = 0;
			for(int n = 0; n<ic->nAH[index]; n++){
				idet = ic->AllHits[idxoff_global+n].detectorID;
				ipos = ic->AllHits[idxoff_global+n].pos;
				for(short ip = 0; ip<4; ip++){
					iprop = 48+ip;
					xExp = trackXZ.tx_*planes->z[iprop]+trackXZ.x_0;
					//loop on hits to find prop
					if(idet==iprop){
						if(fabs(ipos-xExp)<5.08f){
							nprop++;
							break;
						}
					}
				}
				if(nprop>0)break;
			}
			if(nprop==0)continue;

			//partially fill the tracklet xz info:
			oc->AllTracklets[idxoff_tkl+ntkl].stationID = 5;
			oc->AllTracklets[idxoff_tkl+ntkl].x0 = trackXZ.x_0;
			oc->AllTracklets[idxoff_tkl+ntkl].err_x0 = trackXZ.err_x_0;
			oc->AllTracklets[idxoff_tkl+ntkl].tx = trackXZ.tx_;
			oc->AllTracklets[idxoff_tkl+ntkl].err_tx = trackXZ.err_tx_;

			for(i_hit = 0; i_hit<nhits_x; i_hit++){
				FillChi2Arrays(i_hit, ic->AllHits[idxoff_global+trackXZ.hitlist[i_hit]], 0, fitarrays, planes);
				oc->AllTracklets[idxoff_tkl+ntkl].hits[i_hit] = ic->AllHits[idxoff_global+trackXZ.hitlist[i_hit]];
			}
			oc->AllTracklets[idxoff_tkl+ntkl].nXHits = nhits_x;
			
			for(i_uv = 0; i_uv<ncomb_uv; i_uv++){
				nhits_uv = 0;
				i_u2 = i_uv%nu2;
				i_v2 = ((i_uv-i_u2)/nu2)%nv2;
				i_u3 = (((i_uv-i_u2)/nu2-i_v2)/nv2)%nu3;
				i_v3 = ((((i_uv-i_u2)/nu2)-i_v2)/nv2-i_u3)/nu3;
				
				if(straighttrackbuilder->hitpairs_u2[st2uv_offset+bin2+nbins_st2*i_u2].first>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_u2[st2uv_offset+bin2+nbins_st2*i_u2].first], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_u2[st2uv_offset+bin2+nbins_st2*i_u2].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder->hitpairs_u2[st2uv_offset+bin2+nbins_st2*i_u2].second>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_u2[st2uv_offset+bin2+nbins_st2*i_u2].second], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_u2[st2uv_offset+bin2+nbins_st2*i_u2].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}
				
				nhits_u2 = nhits_uv;
				if(nhits_u2==0) continue;
				
				if(straighttrackbuilder->hitpairs_v2[st2uv_offset+bin2+nbins_st2*i_v2].first>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_v2[st2uv_offset+bin2+nbins_st2*i_v2].first], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_v2[st2uv_offset+bin2+nbins_st2*i_v2].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder->hitpairs_v2[st2uv_offset+bin2+nbins_st2*i_v2].second>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_v2[st2uv_offset+bin2+nbins_st2*i_v2].second], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_v2[st2uv_offset+bin2+nbins_st2*i_v2].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}

				nhits_v2 = nhits_uv-nhits_u2;
				if(nhits_v2==0) continue;
				
				if(straighttrackbuilder->hitpairs_u3[st3uv_offset+bin3+nbins_st3*i_u3].first>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_u3[st3uv_offset+bin3+nbins_st3*i_u3].first], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_u3[st3uv_offset+bin3+nbins_st3*i_u3].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder->hitpairs_u3[st3uv_offset+bin3+nbins_st3*i_u3].second>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_u3[st3uv_offset+bin3+nbins_st3*i_u3].second], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_u3[st3uv_offset+bin3+nbins_st3*i_u3].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}

				nhits_u3 = nhits_uv-nhits_u2-nhits_v2;
				if(nhits_u3==0) continue;

				if(straighttrackbuilder->hitpairs_v3[st3uv_offset+bin3+nbins_st3*i_v3].first>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_v3[st3uv_offset+bin3+nbins_st3*i_v3].first], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_v3[st3uv_offset+bin3+nbins_st3*i_v3].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder->hitpairs_v3[st3uv_offset+bin3+nbins_st3*i_v3].second>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_v3[st3uv_offset+bin3+nbins_st3*i_v3].second], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_v3[st3uv_offset+bin3+nbins_st3*i_v3].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}

				nhits_v3 = nhits_uv-nhits_u2-nhits_v2-nhits_u3;
				if(nhits_v3==0) continue;

				fit_2D_track(nhits_uv, fitarrays->y_array, fitarrays->z_array, fitarrays->dy_array, A_, Ainv_, B_, Par, ParErr, chi2);
				
				trackYZ.x_0 = Par[0];
				trackYZ.err_x_0 = ParErr[0];
				trackYZ.tx_ = Par[1];
				trackYZ.err_tx_ = ParErr[1];
				trackYZ.nhits = nhits_uv;
				
				// now evaluate the back track candidate.
				oc->AllTracklets[idxoff_tkl+ntkl].y0 = trackYZ.x_0;
				oc->AllTracklets[idxoff_tkl+ntkl].err_y0 = trackYZ.err_x_0;
				oc->AllTracklets[idxoff_tkl+ntkl].ty = trackYZ.tx_;
				oc->AllTracklets[idxoff_tkl+ntkl].err_ty = trackYZ.err_tx_;

				
				// filter with hodoscope matching before evaluating chi2.
				if(!match_tracklet_to_hodo(oc[index].AllTracklets[ntkl], 2, ic, planes))continue;
				if(!match_tracklet_to_hodo(oc[index].AllTracklets[ntkl], 3, ic, planes) &&
				 !match_tracklet_to_hodo(oc[index].AllTracklets[ntkl], 4, ic, planes))continue;
				
				for(i_hit = nhits_x; i_hit<nhits_x+nhits_uv; i_hit++){
					oc->AllTracklets[idxoff_tkl+ntkl].hits[i_hit] = ic->AllHits[idxoff_global+trackYZ.hitlist[i_hit-nhits_x]];
				}
				oc->AllTracklets[idxoff_tkl+ntkl].nUHits = nhits_u2+nhits_u3;
				oc->AllTracklets[idxoff_tkl+ntkl].nVHits = nhits_v2+nhits_v3;
				
				resolve_leftright(oc[index].AllTracklets[ntkl], planes, 40.);
				resolve_leftright(oc[index].AllTracklets[ntkl], planes, 150.);
				resolve_single_leftright(oc[index].AllTracklets[ntkl], planes);
			
				chi2 = refit_backpartialtrack_with_drift(oc[index].AllTracklets[ntkl], fitarrays, planes);

				resolve_leftright(oc[index].AllTracklets[ntkl], planes, 40.);
				resolve_leftright(oc[index].AllTracklets[ntkl], planes, 150.);
				resolve_single_leftright(oc[index].AllTracklets[ntkl], planes);
			
				for(i_hit = 0; i_hit<nhits_x; i_hit++){
					trackXZ.hitsign[i_hit] = oc->AllTracklets[idxoff_tkl+ntkl].hitsign[i_hit];
				}
				for(i_hit = nhits_x; i_hit<nhits_x+nhits_uv; i_hit++){
					trackYZ.hitsign[i_hit-nhits_x] = oc->AllTracklets[idxoff_tkl+ntkl].hitsign[i_hit];
				}
				
				if(best_candyz_only){
					if(chi2<chi2min){
						chi2min = chi2;
						besttrackYZ.x_0 = oc->AllTracklets[idxoff_tkl+ntkl].y0;
						besttrackYZ.err_x_0 = oc->AllTracklets[idxoff_tkl+ntkl].err_y0;
						besttrackYZ.tx_ = oc->AllTracklets[idxoff_tkl+ntkl].ty;
						besttrackYZ.err_tx_ = oc->AllTracklets[idxoff_tkl+ntkl].err_ty;
						besttrackYZ.nhits = nhits_uv;
						nhits_v = nhits_v2+nhits_v3;
						for(i_hit = 0; i_hit<nhits_uv; i_hit++){
							besttrackYZ.hitlist[i_hit] = trackYZ.hitlist[i_hit];
							besttrackYZ.hitsign[i_hit] = trackYZ.hitsign[i_hit];
						}
					}
				}else{
					if(ntkl<TrackletSizeMax){
						ntkl++;
						oc->AllTracklets[idxoff_tkl+ntkl].chisq = chi2;
						oc->AllTracklets[idxoff_tkl+ntkl].stationID = 5;
						oc->AllTracklets[idxoff_tkl+ntkl].x0 = trackXZ.x_0;
						oc->AllTracklets[idxoff_tkl+ntkl].err_x0 = trackXZ.err_x_0;
						oc->AllTracklets[idxoff_tkl+ntkl].tx = trackXZ.tx_;
						oc->AllTracklets[idxoff_tkl+ntkl].err_tx = trackXZ.err_tx_;
						oc->AllTracklets[idxoff_tkl+ntkl].y0 = besttrackYZ.x_0;
						oc->AllTracklets[idxoff_tkl+ntkl].err_y0= besttrackYZ.err_x_0;
						oc->AllTracklets[idxoff_tkl+ntkl].ty = besttrackYZ.tx_;
						oc->AllTracklets[idxoff_tkl+ntkl].err_ty = besttrackYZ.err_tx_;
					}
				}
			}// end loop on uv hits

			if(best_candyz_only){
				if(nhits_v<0 || chi2min>=10000.f)continue;
			
				oc->AllTracklets[idxoff_tkl+ntkl].chisq = chi2min;
				oc->AllTracklets[idxoff_tkl+ntkl].y0 = besttrackYZ.x_0;
				oc->AllTracklets[idxoff_tkl+ntkl].err_y0= besttrackYZ.err_x_0;
				oc->AllTracklets[idxoff_tkl+ntkl].ty = besttrackYZ.tx_;
				oc->AllTracklets[idxoff_tkl+ntkl].err_ty = besttrackYZ.err_tx_;
			
				nhits_uv = trackYZ.nhits;
			
				for(i_hit = 0; i_hit<nhits_uv; i_hit++){
					oc->AllTracklets[idxoff_tkl+ntkl].hits[i_hit+nhits_x] = ic->AllHits[idxoff_global+trackYZ.hitlist[i_hit]];
				}
				if(ntkl<TrackletSizeMax)ntkl++;	
			}
		}// end loop on x hits
		//if(n_goodxz==0)printf("bin2 %d bin3 %d\n", bin2, bin3);
	}//end loop on bins
		
	oc->nTracklets[index] = ntkl;
}

#ifdef OLDCODE

// ------------------------------
// Global Track candidates
// ------------------------------


__global__ void gKernel_GlobalTrack_building(gEvent* ic, gOutputEvent* oc, gFullTrackBuilder* fulltrackbuilder, gStraightFitArrays* fitarrays, const gPlane* planes, bool ktracker)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;

	//if(oc[index].HasTooManyHits)return;
	
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
			xpos1 = fulltrackbuilder[index].hitpairs_x1[j].second>=0 ? 0.5f*(ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].pos+ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].pos) : ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].pos;
			
			umin1 = xpos1*planes[geometry::detsuperid[0][1]*2].costheta-planes[geometry::detsuperid[0][1]*2].u_win;
			umax1 = umin1+2*planes[geometry::detsuperid[0][1]*2].u_win;
			
			if(fulltrackbuilder[index].hitpairs_x1[j].first>=0){
				fitarrays->x_array[nhits_X1] = ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].pos;
				fitarrays->dx_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].detectorID].spacing/3.4641f;
				fitarrays->z_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].detectorID].z;
				fulltrackbuilder[index].hitlist[nhits_X1] = fulltrackbuilder[index].hitpairs_x1[j].first;
				nhits_X1++;
			}
			if(fulltrackbuilder[index].hitpairs_x1[j].second>=0){
				fitarrays->x_array[nhits_X1] = ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].pos;
				fitarrays->dx_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].detectorID].spacing/3.4641f;
				fitarrays->z_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].detectorID].z;
				fulltrackbuilder[index].hitlist[nhits_X1] = fulltrackbuilder[index].hitpairs_x1[j].second;
				nhits_X1++;
			}
			
			//fit to obtain the st1 track segment: this time fitting from KMag, not from the origin:
			fitarrays->x_array[nhits_X1] = oc[index].AllTracklets[i].x0+oc[index].AllTracklets[i].tx*geometry::Z_KMAG_BEND;
			fitarrays->dx_array[nhits_X1] = oc[index].AllTracklets[i].err_x0+oc[index].AllTracklets[i].err_tx*geometry::Z_KMAG_BEND;
			fitarrays->z_array[nhits_X1] = geometry::Z_KMAG_BEND;
			
			fit_2D_track(nhits_X1+1, fitarrays->x_array, fitarrays->z_array, fitarrays->dx_array, fitarrays->A, fitarrays->Ainv, fitarrays->B, fitarrays->output_parameters, fitarrays->output_parameters_errors, fitarrays->chi2_2d);
			
			fulltrackbuilder[index].TrackXZ_st1.x_0 = fitarrays->output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.err_x_0 = fitarrays->output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.tx_ = fitarrays->output_parameters[1];
			fulltrackbuilder[index].TrackXZ_st1.err_tx_ = fitarrays->output_parameters[1];
			
			for(int k = 0; k<nu1; k++){
				nhits_U1 = nhits_X1;
				upos1 = fulltrackbuilder[index].hitpairs_u1[k].second>=0 ? 0.5f*(ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_u1[k].first].pos+ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_u1[k].second].pos) : ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_u1[k].first].pos;
				
				if(upos1<umin1 || upos1>umax1)continue;
				
				z_x1 = fulltrackbuilder[index].hitpairs_x1[j].second>=0 ? planes[geometry::detsuperid[0][0]*2].z_mean : planes[geometry::detsuperid[0][0]*2].z; 
				z_u1 = fulltrackbuilder[index].hitpairs_u1[k].second>=0 ? planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_u1[k].first].detectorID].z_mean : planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_u1[k].first].detectorID].z; 
				z_v1 = planes[geometry::detsuperid[0][2]*2].z_mean;
				
				v_win1 = planes[ ic->AllHits[idxoff_global+ fulltrackbuilder[index].hitpairs_u1[k].first ].detectorID].v_win_fac1;
				v_win2 = fabs(z_u1+z_v1-2*z_x1)*planes[ geometry::detsuperid[0][2]*2 ].v_win_fac2;
				v_win3 = fabs((z_v1-z_u1)*planes[ geometry::detsuperid[0][2]*2 ].v_win_fac3);
				v_win = v_win1+v_win2+v_win3+2*planes[ ic->AllHits[idxoff_global+ fulltrackbuilder[index].hitpairs_u1[k].first ].detectorID].spacing;
			
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
					vpos1 = fulltrackbuilder[index].hitpairs_v1[l].second>=0 ? 0.5f*(ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_v1[l].first].pos+ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_v1[l].second].pos) : ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_v1[l].first].pos;
					
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
						oc[index].AllTracklets[N_tkl].hits[nhits_1] = ic->AllHits[idxoff_global+fulltrackbuilder[index].hitlist[n]];
						nhits_1++;
					}
					
					resolve_leftright(oc[index].AllTracklets[N_tkl], planes, 75.);
					resolve_leftright(oc[index].AllTracklets[N_tkl], planes, 150.);
					resolve_single_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes);
					
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
				fitarrays->x_array[nhits_X1] = ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].pos;
				fitarrays->dx_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].detectorID].spacing/3.4641f;
				fitarrays->z_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].detectorID].z;
				fulltrackbuilder[index].hitlist[nhits_X1] = fulltrackbuilder[index].hitpairs_x1[j].first;
				nhits_X1++;
			}
			if(fulltrackbuilder[index].hitpairs_x1[j].second>=0){
				fitarrays->x_array[nhits_X1] = ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].pos;
				fitarrays->dx_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].detectorID].spacing/3.4641f;
				fitarrays->z_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].detectorID].z;
				fulltrackbuilder[index].hitlist[nhits_X1] = fulltrackbuilder[index].hitpairs_x1[j].second;
				nhits_X1++;
			}
			
			//fit to obtain the st1 track segment: this time fitting from KMag, not from the origin:
			fitarrays->x_array[nhits_X1] = oc[index].AllTracklets[i].x0+oc[index].AllTracklets[i].tx*geometry::Z_KMAG_BEND;
			fitarrays->dx_array[nhits_X1] = oc[index].AllTracklets[i].err_x0+oc[index].AllTracklets[i].err_tx*geometry::Z_KMAG_BEND;
			fitarrays->z_array[nhits_X1] = geometry::Z_KMAG_BEND;
			
			fit_2D_track(nhits_X1+1, fitarrays->x_array, fitarrays->z_array, fitarrays->dx_array, fitarrays->A, fitarrays->Ainv, fitarrays->B, fitarrays->output_parameters, fitarrays->output_parameters_errors, fitarrays->chi2_2d);
			
			fulltrackbuilder[index].TrackXZ_st1.x_0 = fitarrays->output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.err_x_0 = fitarrays->output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.tx_ = fitarrays->output_parameters[1];
			fulltrackbuilder[index].TrackXZ_st1.err_tx_ = fitarrays->output_parameters[1];

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
						oc[index].AllTracklets[N_tkl].hits[nhits_1] = ic->AllHits[idxoff_global+fulltrackbuilder[index].hitlist[n]];
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


// ------------------------------
// Track Kalman fitting 
// ------------------------------

__global__ void gKernel_GlobalTrack_KalmanFitting(gOutputEvent* oc, gKalmanFitArrays* fitarrays, const gPlane* planes)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	if(oc[index].HasTooManyHits)return;
	
	int nhits_tkl;
	
	for(int k = 0; k<oc[index].nTracklets; k++){
		if(oc[index].AllTracklets[k].stationID<6)continue;
		
		nhits_tkl = oc[index].AllTracklets[k].nXHits+oc[index].AllTracklets[k].nUHits+oc[index].AllTracklets[k].nVHits;
		
		//initialize state
		//initialize_arrays(oc[index].AllTracklets[k], fitarrays[index]);
		
		
			if(oc[index].EventID<20)
		for(int i = 0; i<nhits_tkl; i++){
printf("hit %d, z = %1.3f \n", i, planes[oc[index].AllTracklets[k].hits[i].detectorID].z);
			// first predict the state at the hit z
			predict_state(oc[index].AllTracklets[k], planes[oc[index].AllTracklets[k].hits[i].detectorID].z, fitarrays[index]);
			update_state(oc[index].AllTracklets[k], oc[index].AllTracklets[k].hits[i], fitarrays[index], planes[oc[index].AllTracklets[k].hits[i].detectorID]);
		}
	}
	
	//
}


#endif

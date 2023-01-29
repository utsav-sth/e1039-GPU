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
	ic[index].HasTooManyHits = false;

	int hitmult_[5] = {0, 0, 0, 0, 0};
	for(auto iDetector = 1; iDetector <= 12; ++iDetector)hitmult_[0]+= ic[index].NHits[iDetector];
	if(hitmult_[0]>250)ic[index].HasTooManyHits = true;
	for(auto iDetector = 13; iDetector <= 18; ++iDetector)hitmult_[1]+= ic[index].NHits[iDetector];
	if(hitmult_[1]>200)ic[index].HasTooManyHits = true;
	for(auto iDetector = 19; iDetector <= 24; ++iDetector)hitmult_[2]+= ic[index].NHits[iDetector];
	if(hitmult_[2]>150)ic[index].HasTooManyHits = true;
	for(auto iDetector = 25; iDetector <= 30; ++iDetector)hitmult_[3]+= ic[index].NHits[iDetector];
	if(hitmult_[3]>120)ic[index].HasTooManyHits = true;
	for(auto iDetector = 47; iDetector <= 54; ++iDetector)hitmult_[4]+= ic[index].NHits[iDetector];
	if(hitmult_[4]>250)ic[index].HasTooManyHits = true;
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



__device__ void refit_backpartialtrack_with_drift(gTracklet& tkl, gStraightFitArrays& fitarray, const gPlane* planes){

	//X hits are stored first, so we fit them first;
	for(int i = 0; i<tkl.nXHits; i++){
		FillFitArrays_X(i, tkl.hits[i], tkl.hitsign[i], fitarray, planes);
		FillChi2Arrays(i, tkl.hits[i], tkl.hitsign[i], fitarray, planes);
	}
	fit_2D_track(tkl.nXHits, fitarray.x_array, fitarray.z_array, fitarray.dx_array, fitarray.A, fitarray.Ainv, fitarray.B, fitarray.output_parameters, fitarray.output_parameters_errors, fitarray.chi2_2d);

	tkl.x0 = fitarray.output_parameters[0];
	tkl.err_x0 = fitarray.output_parameters_errors[0];
	
	tkl.tx = fitarray.output_parameters[1];
	tkl.err_tx = fitarray.output_parameters_errors[1];
	
	float y, err_y;
	
	for(int i = tkl.nXHits; i<tkl.nXHits+tkl.nUHits+tkl.nVHits; i++){
		if( calculate_y_uvhit(y, err_y, tkl.hits[i], tkl.hitsign[i], tkl, planes) ){
			FillFitArrays_UV(i-tkl.nXHits, tkl.hits[i], fitarray, planes, y, err_y);
			FillChi2Arrays(i, tkl.hits[i], tkl.hitsign[i], fitarray, planes);
		}
	}
	fit_2D_track(tkl.nUHits+tkl.nVHits, fitarray.y_array, fitarray.z_array, fitarray.dy_array, fitarray.A, fitarray.Ainv, fitarray.B, fitarray.output_parameters, fitarray.output_parameters_errors, fitarray.chi2_2d);
		
	tkl.y0 = fitarray.output_parameters[0];
	tkl.err_y0 = fitarray.output_parameters_errors[0];
	
	tkl.ty = fitarray.output_parameters[1];
	tkl.err_ty = fitarray.output_parameters_errors[1];

	fitarray.output_parameters[0] = tkl.x0;
	fitarray.output_parameters[1] = tkl.y0;
	fitarray.output_parameters[2] = tkl.tx;
	fitarray.output_parameters[3] = tkl.ty;

	chi2_straight(tkl.nXHits+tkl.nUHits+tkl.nVHits, fitarray.drift_dist, fitarray.resolution,
			fitarray.p1x, fitarray.p1y, fitarray.p1z,
			fitarray.deltapx, fitarray.deltapy, fitarray.deltapz,
			fitarray.output_parameters, fitarray.chi2);
	
	tkl.chisq = fitarray.chi2;
}



__global__ void gKernel_XZ_YZ_tracking_new(gEvent* ic, gOutputEvent* oc, gStraightTrackBuilder* straighttrackbuilder, gStraightFitArrays* fitarrays, const gPlane* planes)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	oc[index].EventID = ic[index].EventID;
	oc[index].nAH = ic[index].nAH;
	oc[index].HasTooManyHits = ic[index].HasTooManyHits;
	oc[index].nTracklets = 0;
	
	//if(oc[index].HasTooManyHits)return;
	//if(oc[index].EventID>20)return;
	
	const short nbins_st2 = geometry::N_WCHitsBins[1];
	const short nbins_st3 = geometry::N_WCHitsBins[2]+geometry::N_WCHitsBins[3];
	
	short stid, projid;

        // 1- get hit pairs in bins in st2, st3:
	// the pairs for st3p and st3m are now in the same array. to allow for parallelization of d3p/d3m if needed;
	projid = 0;
	//D2: stid = 3-1
	stid = 2;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder[index].hitpairs_x2, straighttrackbuilder[index].nhitpairs_x2, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid);
	
	//D3p: stid = 4-1
	stid = 3;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder[index].hitpairs_x3, straighttrackbuilder[index].nhitpairs_x3, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid);

	//D3p: stid = 5-1
	stid = 4;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder[index].hitpairs_x3, straighttrackbuilder[index].nhitpairs_x3, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid);
	
#ifdef DEBUG
	for(int bin = 0; bin<58; bin++){
	if(straighttrackbuilder[index].nhitpairs_x3[bin])printf("evt %d bin %d nhits x3 = %d\n", ic[index].EventID, bin, straighttrackbuilder[index].nhitpairs_x3[bin]);
	for(int i = 0; i<straighttrackbuilder[index].nhitpairs_x3[bin]; i++){
		if(straighttrackbuilder[index].hitpairs_x3[bin+28*i].first>-1)printf(" evt %d bin %d  i %d detid1 %d elid1 %d\n", ic[index].EventID, bin, i, 
				ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3[bin+29*i].first].detectorID, ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3[bin+29*i].first].elementID);
		if(straighttrackbuilder[index].hitpairs_x3[bin+28*i].second>-1)printf(" evt %d bin %d i %d detid2 %d elid2 %d\n", ic[index].EventID, bin, i, 
				ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3[bin+29*i].second].detectorID, ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3[bin+29*i].second].elementID);
	}
	}
#endif
	
	projid = 1;
	stid = 2;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder[index].hitpairs_u2, straighttrackbuilder[index].nhitpairs_u2, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid);
	
	stid = 3;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder[index].hitpairs_u3, straighttrackbuilder[index].nhitpairs_u3, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid);
	
	stid = 4;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder[index].hitpairs_u3, straighttrackbuilder[index].nhitpairs_u3, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid);
	
	projid = 2;
	stid = 2;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder[index].hitpairs_v2, straighttrackbuilder[index].nhitpairs_v2, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid);
	
	stid = 3;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder[index].hitpairs_v3, straighttrackbuilder[index].nhitpairs_v3, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid);
		
	stid = 4;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder[index].hitpairs_v3, straighttrackbuilder[index].nhitpairs_v3, straighttrackbuilder[index].hitidx1, straighttrackbuilder[index].hitidx2, straighttrackbuilder[index].hitflag1, straighttrackbuilder[index].hitflag2, stid, projid);
	
	bool bin_overflows = false;
	for(int bin = 0; bin<28; bin++){
		if(straighttrackbuilder[index].nhitpairs_x2[bin]>geometry::MaxHitsProj[0])bin_overflows = true;
		//if(straighttrackbuilder[index].nhitpairs_x2[bin])printf("evt %d bin %d nhits x2 = %d\n", ic[index].EventID, bin, straighttrackbuilder[index].nhitpairs_x2[bin]);
		if(straighttrackbuilder[index].nhitpairs_u2[bin]>geometry::MaxHitsProj[1])bin_overflows = true;
		//if(straighttrackbuilder[index].nhitpairs_u2[bin])printf("evt %d bin %d nhits u2 = %d\n", ic[index].EventID, bin, straighttrackbuilder[index].nhitpairs_u2[bin]);
		if(straighttrackbuilder[index].nhitpairs_v2[bin]>geometry::MaxHitsProj[2])bin_overflows = true;
		//if(straighttrackbuilder[index].nhitpairs_v2[bin])printf("evt %d bin %d nhits v2 = %d\n", ic[index].EventID, bin, straighttrackbuilder[index].nhitpairs_v2[bin]);
	}
	for(int bin = 0; bin<58; bin++){
		if(straighttrackbuilder[index].nhitpairs_x3[bin]>geometry::MaxHitsProj[0])bin_overflows = true;
	//if(straighttrackbuilder[index].nhitpairs_x3[bin])printf("evt %d bin %d nhits x3 = %d\n", ic[index].EventID, bin, straighttrackbuilder[index].nhitpairs_x3[bin]);
		if(straighttrackbuilder[index].nhitpairs_u3[bin]>geometry::MaxHitsProj[1])bin_overflows = true;
		//if(straighttrackbuilder[index].nhitpairs_u3[bin])printf("evt %d bin %d nhits u3 = %d\n", ic[index].EventID, bin, straighttrackbuilder[index].nhitpairs_u3[bin]);
		if(straighttrackbuilder[index].nhitpairs_v3[bin]>geometry::MaxHitsProj[2])bin_overflows = true;
		//if(straighttrackbuilder[index].nhitpairs_v3[bin])printf("evt %d bin %d nhits v3 = %d\n", ic[index].EventID, bin, straighttrackbuilder[index].nhitpairs_v3[bin]);
	}
	if(bin_overflows){
		printf("Evt %d discarded, too many hits\n", ic[index].EventID);
		return;
	}
	
	// declare variables useful for the loop 
	int ntkl = 0;
	
	//for the time being, declare in function;
	short nprop, iprop;
	REAL xExp;
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
	
	short i_x, i_u, i_v;
	
	float y, err_y;
	float ty;
	
	float chi2min = 10000.1f;
	
	int n_goodxz;
	//int hitidx;
	//gHit hit;
	
	//loop on bins FIRST
	for(short i = 0; i<nbins_total; i++){
		bin2 = i%nbins_st2;
		bin3 = (i-bin2)/nbins_st2;
		//printf("bin %d %d\n", bin2, bin3);
		// if parallel:
		// bin3+=nbins_st3*i_thread;

		nx2 = straighttrackbuilder[index].nhitpairs_x2[bin2];
		nx3 = straighttrackbuilder[index].nhitpairs_x3[bin3];
		
		//printf("nx %d %d \n", nx2, nx3);
		if(nx2 == 0 || nx3==0) continue;
		
		nu2 = straighttrackbuilder[index].nhitpairs_u2[bin2];
		nu3 = straighttrackbuilder[index].nhitpairs_u3[bin3];

		//printf("nu %d %d \n", nu2, nu3);
		if(nu2 == 0 || nu3==0) continue;
		
		nv2 = straighttrackbuilder[index].nhitpairs_v2[bin2];
		nv3 = straighttrackbuilder[index].nhitpairs_v3[bin3];

		//printf("nv %d %d \n", nv2, nv3);
		if(nv2 == 0 || nv3==0) continue;
		
		// evaluating the number of combinations;
		ncomb_x = nx2*nx3;
		ncomb_uv = nu2*nu3*nv2*nv3;
		
		n_goodxz = 0;
				
		int i_x, i_uv, i_hit;

		for(i_x = 0; i_x<ncomb_x; i_x++){
			i_x2 = i_x%nx3;
			i_x3 = (i_x-i_x2)/nx3;
			
			nhits_v = -1;
			nhits_x = 0;
			
			if(straighttrackbuilder[index].hitpairs_x2[bin2+nbins_st2*i_x2].first>=0){
				straighttrackbuilder[index].trackXZ.hitlist[nhits_x] = straighttrackbuilder[index].hitpairs_x2[bin2+nbins_st2*i_x2].first;
				FillFitArrays_X(nhits_x, ic[index].AllHits[straighttrackbuilder[index].trackXZ.hitlist[nhits_x]], 0, fitarrays[index], planes);
				nhits_x++;
			}
			if(straighttrackbuilder[index].hitpairs_x2[bin2+nbins_st2*i_x2].second>=0){
				straighttrackbuilder[index].trackXZ.hitlist[nhits_x] = straighttrackbuilder[index].hitpairs_x2[bin2+nbins_st2*i_x2].second;
				FillFitArrays_X(nhits_x, ic[index].AllHits[straighttrackbuilder[index].trackXZ.hitlist[nhits_x]], 0, fitarrays[index], planes);
				nhits_x++;
			}
			
			nhits_x2 = nhits_x;
			if(nhits_x2==0) continue;
			
			if(straighttrackbuilder[index].hitpairs_x3[bin3+nbins_st3*i_x3].first>=0){
				straighttrackbuilder[index].trackXZ.hitlist[nhits_x] = straighttrackbuilder[index].hitpairs_x3[bin3+nbins_st3*i_x3].first;
				FillFitArrays_X(nhits_x, ic[index].AllHits[straighttrackbuilder[index].trackXZ.hitlist[nhits_x]], 0, fitarrays[index], planes);
				nhits_x++;
			}
			if(straighttrackbuilder[index].hitpairs_x3[bin3+nbins_st3*i_x3].second>=0){
				straighttrackbuilder[index].trackXZ.hitlist[nhits_x] = straighttrackbuilder[index].hitpairs_x3[bin3+nbins_st3*i_x3].second;
				FillFitArrays_X(nhits_x, ic[index].AllHits[straighttrackbuilder[index].trackXZ.hitlist[nhits_x]], 0, fitarrays[index], planes);
				nhits_x++;
			}
			
			nhits_x3 = nhits_x-nhits_x2;
			if(nhits_x3==0) continue;
			
			fit_2D_track(nhits_x, fitarrays[index].x_array, fitarrays[index].z_array, fitarrays[index].dx_array, fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B, fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, fitarrays[index].chi2_2d);
			if(fabs(fitarrays[index].output_parameters[0])>1.05*X0_MAX || fabs(fitarrays[index].output_parameters[1])>1.05*TX_MAX)continue;
			
			straighttrackbuilder[index].trackXZ.x_0 = fitarrays[index].output_parameters[0];
			straighttrackbuilder[index].trackXZ.err_x_0 = fitarrays[index].output_parameters_errors[0];
			straighttrackbuilder[index].trackXZ.tx_ = fitarrays[index].output_parameters[1];
			straighttrackbuilder[index].trackXZ.err_tx_ = fitarrays[index].output_parameters_errors[1];
			straighttrackbuilder[index].trackXZ.nhits = nhits_x;
			
			n_goodxz++;
			
			//prop matching
			nprop = 0;
			for(short ip = 0; ip<4; ip++){
				iprop = 48+ip;
				xExp = straighttrackbuilder[index].trackXZ.tx_*planes[iprop].z+straighttrackbuilder[index].trackXZ.x_0;
				//loop on hits to find prop
				for(int n = 0; n<ic[index].nAH; n++){
					if(ic[index].AllHits[n].detectorID==iprop){
						if(fabs(ic[index].AllHits[n].pos-xExp)<5.08f){
							nprop++;
							break;
						}
					}
				}
				if(nprop>0)break;
			}
			if(nprop==0)continue;

			//partially fill the tracklet xz info:
			oc[index].AllTracklets[ntkl].x0 = straighttrackbuilder[index].trackXZ.x_0;
			oc[index].AllTracklets[ntkl].err_x0 = straighttrackbuilder[index].trackXZ.err_x_0;
			oc[index].AllTracklets[ntkl].tx = straighttrackbuilder[index].trackXZ.tx_;
			oc[index].AllTracklets[ntkl].err_tx = straighttrackbuilder[index].trackXZ.err_tx_;

			for(i_hit = 0; i_hit<nhits_x; i_hit++){
				FillChi2Arrays(i_hit, ic[index].AllHits[straighttrackbuilder[index].trackXZ.hitlist[i_hit]], 0, fitarrays[index], planes);
				oc[index].AllTracklets[ntkl].hits[i_hit] = ic[index].AllHits[straighttrackbuilder[index].trackXZ.hitlist[i_hit]];
			}
			oc[index].AllTracklets[ntkl].nXHits = nhits_x;
			
			for(i_uv = 0; i_uv<ncomb_uv; i_uv++){
				nhits_uv = 0;
				i_u2 = i_uv%nu2;
				i_v2 = ((i_uv-i_u2)/nu2)%nv2;
				i_u3 = (((i_uv-i_u2)/nu2-i_v2)/nv2)%nu3;
				i_v3 = ((((i_uv-i_u2)/nu2)-i_v2)/nv2-i_u3)/nu3;
				
#ifdef DEBUG
				hitidx = straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].first;
				hit = ic[index].AllHits[hitidx];

				if(hitidx>=0){
					if(calculate_y_uvhit(y, err_y, hit, 0, straighttrackbuilder[index].trackXZ, planes)){
					//straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = hitidx;
					//FillFitArrays_UV(nhits_x+nhits_uv, hit, fitarrays[index], planes, y, err_y);
					fitarrays[index].y_array[nhits_x+nhits_uv] = y;
					nhits_uv++;
					}
				}
#endif
				
				if(straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].first>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].first], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].first;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}

				if(straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].second>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].second], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].second;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}
				
				nhits_u2 = nhits_uv;
				if(nhits_u2==0) continue;
				
				if(straighttrackbuilder[index].hitpairs_v2[bin2+nbins_st2*i_v2].first>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[bin2+nbins_st2*i_v2].first], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_v2[bin2+nbins_st2*i_v2].first;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder[index].hitpairs_v2[bin2+nbins_st2*i_v2].second>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[bin2+nbins_st2*i_v2].second], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_v2[bin2+nbins_st2*i_v2].second;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}

				nhits_v2 = nhits_uv-nhits_u2;
				if(nhits_v2==0) continue;
				
				if(straighttrackbuilder[index].hitpairs_u3[bin3+nbins_st3*i_u3].first>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3[bin3+nbins_st3*i_u3].first], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_u3[bin3+nbins_st3*i_u3].first;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder[index].hitpairs_u3[bin3+nbins_st3*i_u3].second>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3[bin3+nbins_st3*i_u3].second], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_u3[bin3+nbins_st3*i_u3].second;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}

				nhits_u3 = nhits_uv-nhits_u2-nhits_v2;
				if(nhits_u3==0) continue;

				if(straighttrackbuilder[index].hitpairs_v3[bin3+nbins_st3*i_v3].first>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3[bin3+nbins_st3*i_v3].first], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_v3[bin3+nbins_st3*i_v3].first;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder[index].hitpairs_v3[bin3+nbins_st3*i_v3].second>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3[bin3+nbins_st3*i_v3].second], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_v3[bin3+nbins_st3*i_v3].second;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}

				nhits_v3 = nhits_uv-nhits_u2-nhits_v2-nhits_u3;
				if(nhits_v3==0) continue;
				
				fit_2D_track(nhits_uv, fitarrays[index].y_array, fitarrays[index].z_array, fitarrays[index].dy_array, fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B, fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, fitarrays[index].chi2_2d);
				
				straighttrackbuilder[index].trackYZ.x_0 = fitarrays[index].output_parameters[0];
				straighttrackbuilder[index].trackYZ.err_x_0 = fitarrays[index].output_parameters_errors[0];
				straighttrackbuilder[index].trackYZ.tx_ = fitarrays[index].output_parameters[1];
				straighttrackbuilder[index].trackYZ.err_tx_ = fitarrays[index].output_parameters_errors[1];
				straighttrackbuilder[index].trackYZ.nhits = nhits_uv;
				
				// now evaluate the back track candidate.
				oc[index].AllTracklets[ntkl].y0 = straighttrackbuilder[index].trackYZ.x_0;
				oc[index].AllTracklets[ntkl].err_y0 = straighttrackbuilder[index].trackYZ.err_x_0;
				oc[index].AllTracklets[ntkl].ty = straighttrackbuilder[index].trackYZ.tx_;
				oc[index].AllTracklets[ntkl].err_ty = straighttrackbuilder[index].trackYZ.err_tx_;

				
				// filter with hodoscope matching before evaluating chi2.
				if(!match_tracklet_to_hodo(oc[index].AllTracklets[ntkl], 2, ic, planes))continue;
				if(!match_tracklet_to_hodo(oc[index].AllTracklets[ntkl], 3, ic, planes) &&
				 !match_tracklet_to_hodo(oc[index].AllTracklets[ntkl], 4, ic, planes))continue;

				chi2_straight(nhits_x+nhits_uv, fitarrays[index].drift_dist, fitarrays[index].resolution,
					fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z,
					fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz,
					fitarrays[index].output_parameters, fitarrays[index].chi2);
				
				for(i_hit = nhits_x; i_hit<nhits_x+nhits_uv; i_hit++){
					oc[index].AllTracklets[ntkl].hits[i_hit] = ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[i_hit-nhits_x]];
				}
				oc[index].AllTracklets[ntkl].nUHits = nhits_u2+nhits_u3;
				oc[index].AllTracklets[ntkl].nVHits = nhits_v2+nhits_v3;
				
				resolve_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes, 40.);
				resolve_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes, 150.);
				resolve_single_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes);
			
				refit_backpartialtrack_with_drift(oc[index].AllTracklets[ntkl], fitarrays[index], planes);
				for(i_hit = 0; i_hit<nhits_x; i_hit++){
					straighttrackbuilder[index].trackXZ.hitsign[i_hit] = oc[index].AllTracklets[ntkl].hitsign[i_hit];
				}
				for(i_hit = nhits_x; i_hit<nhits_x+nhits_uv; i_hit++){
					straighttrackbuilder[index].trackYZ.hitsign[i_hit-nhits_x] = oc[index].AllTracklets[ntkl].hitsign[i_hit];
				}

				if(fitarrays[index].chi2<chi2min){
					chi2min = fitarrays[index].chi2;
					straighttrackbuilder[index].besttrackYZ.x_0 = oc[index].AllTracklets[ntkl].y0;
					straighttrackbuilder[index].besttrackYZ.err_x_0 = oc[index].AllTracklets[ntkl].err_y0;
					straighttrackbuilder[index].besttrackYZ.tx_ = oc[index].AllTracklets[ntkl].ty;
					straighttrackbuilder[index].besttrackYZ.err_tx_ = oc[index].AllTracklets[ntkl].err_ty;
					straighttrackbuilder[index].besttrackYZ.nhits = nhits_uv;
					nhits_v = nhits_v2+nhits_v3;
					for(i_hit = 0; i_hit<nhits_uv; i_hit++){
						straighttrackbuilder[index].besttrackYZ.hitlist[i_hit] = straighttrackbuilder[index].trackYZ.hitlist[i_hit];
						straighttrackbuilder[index].besttrackYZ.hitsign[i_hit] = straighttrackbuilder[index].trackYZ.hitsign[i_hit];
					}
				}
			}// end loop on uv hits

			if(nhits_v<0 || chi2min>=10000.f)continue;
			
			oc[index].AllTracklets[ntkl].chisq = chi2min;
			oc[index].AllTracklets[ntkl].y0 = straighttrackbuilder[index].besttrackYZ.x_0;
			oc[index].AllTracklets[ntkl].err_y0= straighttrackbuilder[index].besttrackYZ.err_x_0;
			oc[index].AllTracklets[ntkl].ty = straighttrackbuilder[index].besttrackYZ.tx_;
			oc[index].AllTracklets[ntkl].err_ty = straighttrackbuilder[index].besttrackYZ.err_tx_;
			
			nhits_uv = straighttrackbuilder[index].trackYZ.nhits;
			
			for(i_hit = 0; i_hit<nhits_uv; i_hit++){
				oc[index].AllTracklets[ntkl].hits[i_hit+nhits_x] = ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[i_hit]];
			}
			if(ntkl<TrackletSizeMax)ntkl++;	

		}// end loop on x hits
		//if(n_goodxz==0)printf("bin2 %d bin3 %d\n", bin2, bin3);
	}//end loop on bins
		
	oc[index].nTracklets = ntkl;
}


#ifdef OLDCODE
//
__global__ void gKernel_XZ_YZ_tracking(gEvent* ic, gOutputEvent* oc, gStraightTrackBuilder* straighttrackbuilder, gStraightFitArrays* fitarrays, const gPlane* planes, const bool best_cand_yz_only = false)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	straighttrackbuilder[index].nTracksXZ = 0;
	straighttrackbuilder[index].nTracksYZ = 0;
	
	oc[index].EventID = ic[index].EventID;
	oc[index].nAH = ic[index].nAH;
	oc[index].nTracklets = 0;
	
	//for the time being, declare in function;
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
				if(straighttrackbuilder[index].hitpairs_x3m[j-nx3p].first>=0){
					fitarrays[index].x_array[nhits_X3] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[j-nx3p].first].pos;
					fitarrays[index].dx_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[j-nx3p].first].detectorID].spacing/3.4641f;
					fitarrays[index].z_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[j-nx3p].first].detectorID].z;
					straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].hitlist[nhits_X3] = straighttrackbuilder[index].hitpairs_x3m[j-nx3p].first;
					nhits_X3++;
				}
				if(straighttrackbuilder[index].hitpairs_x3m[j-nx3p].second>=0){
					fitarrays[index].x_array[nhits_X3] = ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[j-nx3p].second].pos;
					fitarrays[index].dx_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[j-nx3p].second].detectorID].spacing/3.4641f;
					fitarrays[index].z_array[nhits_X3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_x3m[j-nx3p].second].detectorID].z;
					straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].hitlist[nhits_X3] = straighttrackbuilder[index].hitpairs_x3m[j-nx3p].second;
					nhits_X3++;
				}
			}

			fit_2D_track(nhits_X3, fitarrays[index].x_array, fitarrays[index].z_array, fitarrays[index].dx_array, fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B, fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, fitarrays[index].chi2_2d);
			if(fabs(fitarrays[index].output_parameters[0])>1.05*X0_MAX || fabs(fitarrays[index].output_parameters[1])>1.05*TX_MAX)continue;
			
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
			if(j<nx3p){
				straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].d3p = true;
			}else{
				straighttrackbuilder[index].TrackXZ[straighttrackbuilder[index].nTracksXZ].d3p = false;
			}
			if(straighttrackbuilder[index].nTracksXZ>=Track2DSizeMax)break;//continue;
			
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
	int bestYZcand = -1;
	bool d3p;

	float y, err_y;
	// loop on XZ tracks
	float chi2min = 100000.;//-1.0f;
	bool yz_in_acceptance;
	
	for(int i = 0; i<straighttrackbuilder[index].nTracksXZ; i++){
		nxhits = straighttrackbuilder[index].TrackXZ[i].nXHits;
		bestYZcand = -1;
		chi2min = 100000.;//-1.0f;
		straighttrackbuilder[index].nTracksYZ = 0;
		d3p = straighttrackbuilder[index].TrackXZ[i].d3p;

#ifdef KTRACKERHITSELECTION
		xpos2 = xpos3p = xpos3m = 0;
		short n2 = 0, n3p = 0, n3m = 0;
		for(int n = 0; n<nxhits; n++){
			FillChi2Arrays(n, ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]], 0, fitarrays[index], planes);

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

		straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].nUHits = 0;
		straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].nVHits = 0;

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
					FillChi2Arrays(nxhits+nhits_U2, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].first], 0, fitarrays[index], planes);
					fitarrays[index].y_array[nhits_U2] = y;
					fitarrays[index].dy_array[nhits_U2] = err_y;
					fitarrays[index].z_array[nhits_U2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].first].detectorID].z;
					straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_U2] = straighttrackbuilder[index].hitpairs_u2[j].first;
					nhits_U2++;
				}
			}
			if(straighttrackbuilder[index].hitpairs_u2[j].second>=0){
				if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].second], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
					FillChi2Arrays(nxhits+nhits_U2, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[j].second], 0, fitarrays[index], planes);
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
					if(!d3p)continue;
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
							FillChi2Arrays(nxhits+nhits_U3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].first], 0, fitarrays[index], planes);
							fitarrays[index].y_array[nhits_U3] = y;
							fitarrays[index].dy_array[nhits_U3] = err_y;
							fitarrays[index].z_array[nhits_U3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].first].detectorID].z;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_U3] = straighttrackbuilder[index].hitpairs_u3p[k].first;
							nhits_U3++;
						}
					}
					if(straighttrackbuilder[index].hitpairs_u3p[k].second>=0){
						if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].second], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
							FillChi2Arrays(nxhits+nhits_U3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].second], 0, fitarrays[index], planes);
							fitarrays[index].y_array[nhits_U3] = y;
							fitarrays[index].dy_array[nhits_U3] = err_y;
							fitarrays[index].z_array[nhits_U3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3p[k].second].detectorID].z;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_U3] = straighttrackbuilder[index].hitpairs_u3p[k].second;
							nhits_U3++;
						}
					}
				}else{
					if(d3p)continue;
#ifdef KTRACKERHITSELECTION
					upos3m = straighttrackbuilder[index].hitpairs_u3m[k-nu3p].second>=0 ? 0.5f*(ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k-nu3p].first].pos+ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k-nu3p].second].pos) : ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k-nu3p].first].pos;
			
					if(upos3m<umin3m || upos3m>umax3m)continue;

					z_x3 = n3m==2 ? planes[geometry::detsuperid[4][0]*2].z_mean : planes[geometry::detsuperid[4][0]*2].z; 
					z_u3 = straighttrackbuilder[index].hitpairs_u3m[k-nu3p].second>=0 ? planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k-nu3p].first].detectorID].z_mean : planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k-nu3p].first].detectorID].z; 
					z_v3 = planes[geometry::detsuperid[4][2]*2].z_mean;
			
					v_win1 = planes[ ic[index].AllHits[ straighttrackbuilder[index].hitpairs_u3m[k-nu3p].first ].detectorID].v_win_fac1;
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
					
					if(straighttrackbuilder[index].hitpairs_u3m[k-nu3p].first>=0){
						if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k-nu3p].first], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
							FillChi2Arrays(nxhits+nhits_U3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k-nu3p].first], 0, fitarrays[index], planes);
							fitarrays[index].y_array[nhits_U3] = y;
							fitarrays[index].dy_array[nhits_U3] = err_y;
							fitarrays[index].z_array[nhits_U3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k-nu3p].first].detectorID].z;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_U3] = straighttrackbuilder[index].hitpairs_u3m[k-nu3p].first;
							nhits_U3++;
						}
					}
					if(straighttrackbuilder[index].hitpairs_u3m[k-nu3p].second>=0){
						if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k-nu3p].second], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
							FillChi2Arrays(nxhits+nhits_U3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k-nu3p].second], 0, fitarrays[index], planes);
							fitarrays[index].y_array[nhits_U3] = y;
							fitarrays[index].dy_array[nhits_U3] = err_y;
							fitarrays[index].z_array[nhits_U3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3m[k-nu3p].second].detectorID].z;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].nUHits++;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_U3] = straighttrackbuilder[index].hitpairs_u3m[k-nu3p].second;
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
							FillChi2Arrays(nxhits+nhits_V2, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].first], 0, fitarrays[index], planes);
							fitarrays[index].y_array[nhits_V2] = y;
							fitarrays[index].dy_array[nhits_V2] = err_y;
							fitarrays[index].z_array[nhits_V2] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].first].detectorID].z;
							straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_V2] = straighttrackbuilder[index].hitpairs_v2[l].first;
							nhits_V2++;
						}
					}

					if(straighttrackbuilder[index].hitpairs_v2[l].second>=0){
						if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].second], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
							FillChi2Arrays(nxhits+nhits_V2, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[l].second], 0, fitarrays[index], planes);
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
							if(!d3p)continue;
#ifdef KTRACKERHITSELECTION
							vpos3p = straighttrackbuilder[index].hitpairs_v3p[m].second>=0 ? 0.5f*(ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].first].pos+ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].second].pos) : ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].first].pos;

							if(vpos3p<vmin3p || vpos3p>vmax3p)continue;
#endif

							if(straighttrackbuilder[index].hitpairs_v3p[m].first>=0){
								if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].first], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
									FillChi2Arrays(nxhits+nhits_V3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].first], 0, fitarrays[index], planes);
									fitarrays[index].y_array[nhits_V3] = y;
									fitarrays[index].dy_array[nhits_V3] = err_y;
									fitarrays[index].z_array[nhits_V3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].first].detectorID].z;
									straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_V3] = straighttrackbuilder[index].hitpairs_v3p[m].first;
									nhits_V3++;
								}
							}
						if(straighttrackbuilder[index].hitpairs_v3p[m].second>=0){
								if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].second], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
									FillChi2Arrays(nxhits+nhits_V3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].second], 0, fitarrays[index], planes);
									fitarrays[index].y_array[nhits_V3] = y;
									fitarrays[index].dy_array[nhits_V3] = err_y;
									fitarrays[index].z_array[nhits_V3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3p[m].second].detectorID].z;
									straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_V3] = straighttrackbuilder[index].hitpairs_v3p[m].second;
									nhits_V3++;
								}
							}
						}else{
							if(d3p)continue;
#ifdef KTRACKERHITSELECTION
							vpos3m = straighttrackbuilder[index].hitpairs_v3m[m-nv3p].second>=0 ? 0.5f*(ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m-nv3p].first].pos+ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m-nv3p].second].pos) : ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m-nv3p].first].pos;
							
							if(vpos3m<vmin3m || vpos3m>vmax3m)continue;
#endif

							if(straighttrackbuilder[index].hitpairs_v3m[m-nv3p].first>=0){
								if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m-nv3p].first], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
									FillChi2Arrays(nxhits+nhits_V3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m-nv3p].first], 0, fitarrays[index], planes);
									fitarrays[index].y_array[nhits_V3] = y;
									fitarrays[index].dy_array[nhits_V3] = err_y;
									fitarrays[index].z_array[nhits_V3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m-nv3p].first].detectorID].z;
									straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_V3] = straighttrackbuilder[index].hitpairs_v3m[m-nv3p].first;
									nhits_V3++;
								}
							}
							if(straighttrackbuilder[index].hitpairs_v3m[m-nv3p].second>=0){
								if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m-nv3p].second], 0, straighttrackbuilder[index].TrackXZ[i], planes)){
									FillChi2Arrays(nxhits+nhits_V3, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m-nv3p].second], 0, fitarrays[index], planes);
									fitarrays[index].y_array[nhits_V3] = y;
									fitarrays[index].dy_array[nhits_V3] = err_y;
									fitarrays[index].z_array[nhits_V3] = planes[ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3m[m-nv3p].second].detectorID].z;
									straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[nhits_V3] = straighttrackbuilder[index].hitpairs_v3m[m-nv3p].second;
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
						
						bestYZcand = straighttrackbuilder[index].nTracksYZ;
						
						if(!best_cand_yz_only){
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
						
							//geometry::detmap_zincr[31] 
							//= {-1,  0,  1,  3,  2,  4,  5,  0,  1,  3,  2,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 12, 13, 14, 15, 16, 17};		
							// geometry::detmap_zincr[ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]].detectorID]
						
							for(int n = 0; n<nxhits; n++){
								oc[index].AllTracklets[oc[index].nTracklets].hits[ n ] = ic[index].AllHits[straighttrackbuilder[index].TrackXZ[i].hitlist[n]];
							}
							for(int n = nxhits; n<nxhits+oc[index].AllTracklets[oc[index].nTracklets].nUHits+oc[index].AllTracklets[oc[index].nTracklets].nVHits; n++){
								oc[index].AllTracklets[oc[index].nTracklets].hits[ n ] = ic[index].AllHits[straighttrackbuilder[index].TrackYZ[straighttrackbuilder[index].nTracksYZ].hitlist[n-nxhits]];
							}
						
							resolve_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes, 40.);
							resolve_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes, 150.);
							resolve_single_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes);
							
							if(oc[index].nTracklets>=TrackletSizeMax)break;
							
							//refit after???
							oc[index].nTracklets++;
							
							if(fitarrays[index].chi2<=chi2min){
								chi2min = fitarrays[index].chi2;
								bestYZcand = straighttrackbuilder[index].nTracksYZ;
							}
						}
						straighttrackbuilder[index].nTracksYZ++;
						//
					}// end loop on st3 v hits
				}// end loop on st2 v hits
			}// end loop on st3 u hits
		}// end loop on st2 u hits
		
		//if(straighttrackbuilder[index].nTracksYZ==0)printf("evt %d, %d pairs u2, %d pairs v2, %d pairs v2, %d pairs v3\n", ic[index].EventID, nu2, nu3p+nu3m, nv2, nv3p+nv3m);
		
		if(best_cand_yz_only){
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
				for(int n = nxhits; n<nxhits+oc[index].AllTracklets[oc[index].nTracklets].nUHits+oc[index].AllTracklets[oc[index].nTracklets].nVHits; n++){
					oc[index].AllTracklets[oc[index].nTracklets].hits[n] = ic[index].AllHits[straighttrackbuilder[index].TrackYZ[bestYZcand].hitlist[n-nxhits]];
				}

				resolve_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes, 40.);
				resolve_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes, 150.);
				resolve_single_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes);

#ifdef DEBUG
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
#endif
				
				oc[index].nTracklets++;
			}else{
				//printf("evt %d: track XZ %d without YZ track: nuhits = (%d+%d+%d), nvhits = (%d+%d+%d)\n", ic[index].EventID, i, nu2, nu3p, nu3m, nv2, nv3p, nv3m);
				continue;
			}
		}
	}// end loop on XZ tracks
	
	//reloop in tracklets to refit???
	for(int n = 0; n<oc[index].nTracklets; n++){
		refit_backpartialtrack_with_drift(oc[index].AllTracklets[n], fitarrays[index], planes);
	}

}
#endif


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
			
			fulltrackbuilder[index].TrackXZ_st1.x_0 = fitarrays[index].output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.err_x_0 = fitarrays[index].output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.tx_ = fitarrays[index].output_parameters[1];
			fulltrackbuilder[index].TrackXZ_st1.err_tx_ = fitarrays[index].output_parameters[1];
			
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
			
			fulltrackbuilder[index].TrackXZ_st1.x_0 = fitarrays[index].output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.err_x_0 = fitarrays[index].output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.tx_ = fitarrays[index].output_parameters[1];
			fulltrackbuilder[index].TrackXZ_st1.err_tx_ = fitarrays[index].output_parameters[1];

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

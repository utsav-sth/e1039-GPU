#include "reconstruction_classes.cuh"


// --------------------------------------------------------------- //
// functions to calculate bottom and top end wire points for a hit //
// --------------------------------------------------------------- //

__device__ float x_bep(const gHit hit, const gPlane plane)
{
	return plane.p1x_w1+plane.dp1x*(hit.elementID-1);
}

__device__ float x_tep(const gHit hit, const gPlane plane)
{
	return x_bep(hit, plane)+plane.deltapx;
}

__device__ float y_bep(const gHit hit, const gPlane plane)
{
	return plane.p1y_w1+plane.dp1y*(hit.elementID-1);
}

__device__ float y_tep(const gHit hit, const gPlane plane)
{
	return y_bep(hit, plane)+plane.deltapy;
}

__device__ float z_bep(const gHit hit, const gPlane plane)
{
	return plane.p1z_w1+plane.dp1z*(hit.elementID-1);
}

__device__ float z_tep(const gHit hit, const gPlane plane)
{
	return z_bep(hit, plane)+plane.deltapz;
}


__device__ float position(const gHit hit, const short sign)
{
	return (hit.pos+sign*hit.driftDistance);
}


// ---------------------------------------------------

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
			//p1x = planes[ ic[index].AllHits[i].detectorID ].p1x_w1 + planes[ ic[index].AllHits[i].detectorID ].dp1x * (ic[index].AllHits[i].elementID-1);
			if(planes[ ic[index].AllHits[i].detectorID ].deltapx>0){
				p1x = x_bep(ic[index].AllHits[i], planes[ ic[index].AllHits[i].detectorID ]);
				p2x = x_tep(ic[index].AllHits[i], planes[ ic[index].AllHits[i].detectorID ]);
				//p2x = p1x + planes[ ic[index].AllHits[i].detectorID ].deltapx;
			}else{
				p1x = x_tep(ic[index].AllHits[i], planes[ ic[index].AllHits[i].detectorID ]);
				p2x = x_bep(ic[index].AllHits[i], planes[ ic[index].AllHits[i].detectorID ]);
				//p2x = p1x;
				//p1x+= planes[ ic[index].AllHits[i].detectorID ].deltapx;
			}
			//printf("%d %d %1.6f p1-2x %1.6f %1.6f xmin-max %1.6f %1.6f \n", ic[index].AllHits[i].detectorID, ic[index].AllHits[i].elementID, ic[index].AllHits[i].pos, p1x, p2x, xmin, xmax);
			//if(xmin>-999.)printf("xmin %1.6f xmax %1.6f p1x %1.6f p2x %1.6f \n", xmin, xmax, p1x, p2x);
			if( (p1x <= xmax) && (p2x >= xmin) ){ 
				hitidx1[hitctr1] = i;
				hitctr1++;
			}
		}
		if(ic[index].AllHits[i].detectorID==detid2){
			//p1x = planes[ ic[index].AllHits[i].detectorID ].p1x_w1 + planes[ ic[index].AllHits[i].detectorID ].dp1x * (ic[index].AllHits[i].elementID-1);
			if(planes[ ic[index].AllHits[i].detectorID ].deltapx>0){
				p1x = x_bep(ic[index].AllHits[i], planes[ ic[index].AllHits[i].detectorID ]);
				p2x = x_tep(ic[index].AllHits[i], planes[ ic[index].AllHits[i].detectorID ]);
				//p2x = p1x + planes[ ic[index].AllHits[i].detectorID ].deltapx;
			}else{
				p1x = x_tep(ic[index].AllHits[i], planes[ ic[index].AllHits[i].detectorID ]);
				p2x = x_bep(ic[index].AllHits[i], planes[ ic[index].AllHits[i].detectorID ]);
				//p2x = p1x;
				//p1x+= planes[ ic[index].AllHits[i].detectorID ].deltapx;
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


// ------------------------------------------------------------- //
// functions to evaluate the hit selection window for a 2D track //
// and to calculate y from u, v hits given x                     //
// ------------------------------------------------------------- //


__device__ void find_xmin_xmax_in_chamber(float &xmin, float &xmax, const gTrackXZ trackxz, const short stID, const short projID, const gPlane* planes)
{
	int detid1 = geometry::detsuperid[stID][projID]*2;
	int detid2 = geometry::detsuperid[stID][projID]*2-1;
	
	xmin = min(planes[detid1].z*(trackxz.tx-trackxz.err_tx), planes[detid2].z*(trackxz.tx-trackxz.err_tx));
	xmax = max(planes[detid1].z*(trackxz.tx+trackxz.err_tx), planes[detid2].z*(trackxz.tx+trackxz.err_tx));
	
	xmin = xmin + trackxz.x0-trackxz.err_x0-planes[detid1].spacing;
	xmax = xmax + trackxz.x0+trackxz.err_x0+planes[detid1].spacing;
}

__device__ bool calculate_y_uvhit(float &y, float &err_y, const gHit hit, const short hitsign, const gTrackXZ trackxz, const gPlane* planes){
	float p1x = x_bep( hit, planes[ hit.detectorID ]);
	float p1y = y_bep( hit, planes[ hit.detectorID ]);
	float p2x = x_tep( hit, planes[ hit.detectorID ]);
	//float p1x = planes[ hit.detectorID ].p1x_w1 + planes[ hit.detectorID ].dp1x * (hit.elementID-1);
	//float p1y = planes[ hit.detectorID ].p1y_w1 + planes[ hit.detectorID ].dp1y * (hit.elementID-1);
	//float p2x = p1x+planes[ hit.detectorID ].deltapx;
	
	float x_trk = trackxz.x0+planes[ hit.detectorID ].z*trackxz.tx;

	y = p1y + (x_trk-p1x) *  planes[ hit.detectorID ].deltapy/planes[ hit.detectorID ].deltapx;
	
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

__device__ bool calculate_y_uvhit(float &y, float &err_y, const gHit hit, const short hitsign, const float x0, const float tx, const gPlane* planes){
	float p1x = x_bep( hit, planes[ hit.detectorID ]);
	float p1y = y_bep( hit, planes[ hit.detectorID ]);
	float p2x = x_tep( hit, planes[ hit.detectorID ]);
	//float p1x = planes[ hit.detectorID ].p1x_w1 + planes[ hit.detectorID ].dp1x * (hit.elementID-1);
	//float p1y = planes[ hit.detectorID ].p1y_w1 + planes[ hit.detectorID ].dp1y * (hit.elementID-1);
	//float p2x = p1x+planes[ hit.detectorID ].deltapx;
	
	float x_trk = x0+planes[ hit.detectorID ].z*tx;//x_trk is x_trk
	
	y = p1y + (x_trk-p1x) *  planes[ hit.detectorID ].deltapy/planes[ hit.detectorID ].deltapx;
	
	// we build a virtual wire, parallel to the wire hit, at a distance driftDistance*hitsign from the wire
	// the track will intersect this wire...
	//y = p1y + (x_trk-p1x+hit.driftDistance*hitsign*planes[ hit.detectorID ].costheta) * planes[ hit.detectorID ].deltapy/planes[ hit.detectorID ].deltapx;
	    
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

	err_y = hitsign==0? planes[ hit.detectorID ].spacing/3.4641f : planes[ hit.detectorID ].resolution; 
	err_y*= fabs(planes[ hit.detectorID ].deltapy/planes[ hit.detectorID ].deltapx);

	return true;
}




__device__ void FillFitArrays(const int n, const gHit hit, const short hitsign, gStraightFitArrays &fitarray, const gPlane* planes){
	fitarray.drift_dist[n] = hit.driftDistance*hitsign;
	fitarray.resolution[n] = planes[ hit.detectorID ].resolution;
	if(hitsign==0){
		fitarray.resolution[n] = planes[ hit.detectorID ].spacing*3.4641f;
	}else{
		fitarray.resolution[n] = planes[ hit.detectorID ].resolution;
	}	       
	fitarray.p1x[n] = x_bep( hit, planes[ hit.detectorID ]);
	fitarray.p1y[n] = y_bep( hit, planes[ hit.detectorID ]);
	fitarray.p1z[n] = z_bep( hit, planes[ hit.detectorID ]);
	//fitarray.p1x[n] = planes[ hit.detectorID ].p1x_w1 + planes[ hit.detectorID ].dp1x * (hit.elementID-1);
	//fitarray.p1y[n] = planes[ hit.detectorID ].p1y_w1 + planes[ hit.detectorID ].dp1y * (hit.elementID-1);
	//fitarray.p1z[n] = planes[ hit.detectorID ].p1z_w1 + planes[ hit.detectorID ].dp1z * (hit.elementID-1);
	
	fitarray.deltapx[n] = planes[ hit.detectorID ].deltapx;
	fitarray.deltapy[n] = planes[ hit.detectorID ].deltapy;
	fitarray.deltapz[n] = planes[ hit.detectorID ].deltapz;
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

__device__ float calculate_invP(float tx, float tx_st1)
{
	return (tx_st1 - tx) / geometry::PT_KICK_KMAG;
}

__device__ float calculate_invP_error(float err_tx, float err_tx_st1)
{
	return ( err_tx - err_tx )/ geometry::PT_KICK_KMAG;
} 



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
		
		slope_exp = planes[tkl.hits[i].detectorID].costheta*tx + planes[tkl.hits[i].detectorID].sintheta*tkl.ty;
		err_slope = fabs(planes[tkl.hits[i].detectorID].costheta*err_tx) + fabs(planes[tkl.hits[j].detectorID].sintheta*tkl.err_ty);
		
		inter_exp = planes[tkl.hits[i].detectorID].costheta*x0 + planes[tkl.hits[i].detectorID].sintheta*tkl.y0;
		err_inter = fabs(planes[tkl.hits[i].detectorID].costheta*err_x0) + fabs(planes[tkl.hits[j].detectorID].sintheta*tkl.err_y0);
		
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
				slope_local = ( position(tkl.hits[i], geometry::lrpossibility[k][0]) - position(tkl.hits[j], geometry::lrpossibility[k][1]) )/(planes[tkl.hits[i].detectorID].z-planes[tkl.hits[j].detectorID].z);
				inter_local = position(tkl.hits[i], geometry::lrpossibility[k][0]) - slope_local*planes[tkl.hits[i].detectorID].z;
				
				if(fabs(slope_local) > planes[tkl.hits[i].detectorID].slope_max || fabs(inter_local) > planes[tkl.hits[i].detectorID].inter_max)continue;
				
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
		pos_exp = (planes[detID].z*tx+x0)*planes[detID].costheta+(planes[detID].z*tkl.ty+tkl.y0)*planes[detID].sintheta;
		tkl.hitsign[n] = pos_exp>tkl.hits[n].pos? +1 : -1;
	}
	
}


// ----------------------------------------------------------------- //
// functions for selection of station 1 hits for back partial tracks //
// ----------------------------------------------------------------- // 

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


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
	
	// we build a virtual wire, parallel to the wire hit, at a distance driftDistance*hitsign from the wire
	// the track will intersect this wire...
	//y = p1y + (x_trk-p1x) *  planes[ hit.detectorID ].deltapy/planes[ hit.detectorID ].deltapx;
	
	y = p1y + (x_trk-p1x+hit.driftDistance*hitsign*planes[ hit.detectorID ].costheta) * planes[ hit.detectorID ].deltapy/planes[ hit.detectorID ].deltapx;
	
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
	M_5x1[0] = M_S5x5[0] * M_2x1[0] + M_S5x5[5] * M_2x1[1];
	M_5x1[1] = M_S5x5[1] * M_2x1[0] + M_S5x5[6] * M_2x1[1];
	M_5x1[2] = M_S5x5[2] * M_2x1[0] + M_S5x5[7] * M_2x1[1];
	M_5x1[3] = M_S5x5[3] * M_2x1[0] + M_S5x5[8] * M_2x1[1];
	M_5x1[4] = M_S5x5[4] * M_2x1[0] + M_S5x5[9] * M_2x1[1];
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

__device__ void initialize_arrays(gTracklet& tkl, gKalmanFitArrays fitarray)
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
	fitarray.state[1] = tkl.y0+z*tkl.tx;
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

__device__ void update_state(gTracklet& tkl, const gHit hit, gKalmanFitArrays fitarray, const gPlane plane)//will optimize after if something is redundant
{
	const float dxdy = plane.deltapx/plane.deltapy;
	const float y0 = y_bep(hit, plane);
	const float y1 = y_tep(hit, plane);
	const float x0 = hit.pos*plane.costheta + y0*dxdy;
	const float x1 = x0 + (y1-y0) *dxdy;
	const float z = plane.z;
	
	const float x2y2 = sqrtf( (x1-x0)*(x1-x0) + (y1-y0)*(y1-y0) );
	
	fitarray.H[0] = (y1-y0) / x2y2;
	fitarray.H[1] = (x1-x0) / x2y2;
	
	const float res = fitarray.H[0] * x0 + fitarray.H[1] * y0 - (fitarray.H[0] * fitarray.state[0] + fitarray.H[1] * fitarray.state[1]); 
	const float CRes = similarity_1x2_S5x5_2x1(fitarray.H, fitarray.Cov) + plane.resolution * plane.resolution;
	
	multiply_S5x5_2x1(fitarray.H, fitarray.Cov, fitarray.K);
	
	scale_matrix(fitarray.K, 1./CRes, 5);	
	add_matrix(fitarray.state, fitarray.K, 5, res);
	
	tensor_product(fitarray.K, fitarray.K, fitarray.KCResKt, CRes);
	
	add_matrix(fitarray.Cov, fitarray.KCResKt, 25, -1.f);
	
	fitarray.chi2 = res*res / CRes;
	
	// now update the tracklet parameters!
	
	//update_tkl(tkl, z, fitarray.state);
	update_tkl(tkl, z, fitarray);
}


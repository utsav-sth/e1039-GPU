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


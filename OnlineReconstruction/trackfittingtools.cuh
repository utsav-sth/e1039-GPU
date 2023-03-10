#include "operations.h"

__device__ void matinv_2x2_matrix_per_thread (const float *A, float *Ainv)
{
    //const int blkNum = blockIdx.y * gridDim.x + blockIdx.x;
    //const int thrdNum = blkNum * blockDim.x + threadIdx.x;
    //const int N = 2;
    int perm0, perm1;
    int icol0, icol1;
    float AA00, AA01; 
    float AA10, AA11;
    float tmp;
//#if USE_PIVOTING
    //typename config<T,arch>::absValType t;
    //typename config<T,arch>::absValType p;
    float t;
    float p;
    int i, pvt;
//#endif

//    A    += thrdNum * N * N;
//    Ainv += thrdNum * N * N;

    //if (thrdNum < batch) {

        AA00 = A[0];
        AA10 = A[1];
        AA01 = A[2];
        AA11 = A[3];

        perm0 = 0;
        perm1 = 1;
        
        /****************** iteration 0 ***********/

//#if USE_PIVOTING
        /* search pivot row */
        p = absOp (AA00);
        pvt = 0;
        t = absOp (AA10);
        if (t > p) { p = t;  pvt = 1; }
        
        /* swap pivot row with row 0 */
        if (pvt == 1) {
            tmp = AA00;  AA00 = AA10;  AA10 = tmp;
            tmp = AA01;  AA01 = AA11;  AA11 = tmp;
            /* update permutation vector based on row swap */
            i = perm0;  perm0 = perm1;  perm1 = i;
        }
//#endif // USE_PIVOTING

        /* scale current row */
        tmp = rcpOp (AA00);
        icol0 = perm0;
        AA00 = tmp;
        AA01 = mulOp (tmp, AA01);

        /* eliminate above and below current row */
        tmp = AA10;
        AA10 = mulOp (negOp(tmp), AA00);
        AA11 = fmnaOp (tmp, AA01, AA11);
        
        /****************** iteration 1 ***********/

        /* scale current row */
        tmp = rcpOp (AA11);
        icol1 = perm1;
        AA10 = mulOp (tmp, AA10);
        AA11 = tmp;

        /* eliminate above and below current row */
        tmp = AA01;
        AA00 = fmnaOp (tmp, AA10, AA00);
        AA01 = mulOp (negOp(tmp), AA11);

        /* sort columns into the correct order */
        Ainv[0*2+icol0] = AA00;
        Ainv[1*2+icol0] = AA10;
        Ainv[0*2+icol1] = AA01;
        Ainv[1*2+icol1] = AA11;
    //}
}

__device__ void matinv_4x4_matrix_per_thread (const float *A, float *Ainv)
{
    //const int blkNum = blockIdx.y * gridDim.x + blockIdx.x;
    //const int thrdNum = blkNum * blockDim.x + threadIdx.x;
    //const int N = 4;
    int perm0, perm1, perm2, perm3;
    int icol0, icol1, icol2, icol3;
    float AA00, AA01, AA02, AA03; 
    float AA10, AA11, AA12, AA13;
    float AA20, AA21, AA22, AA23;
    float AA30, AA31, AA32, AA33;
    float tmp;

//#if USE_PIVOTING
    float t;
    float p;
    int i, pvt;
//#endif

    //A    += thrdNum * N * N;
    //Ainv += thrdNum * N * N;

    //if (thrdNum < batch) {

        AA00 = A[0*4+0];
        AA10 = A[1];
        AA20 = A[2];
        AA30 = A[3];
        AA01 = A[4];
        AA11 = A[5];
        AA21 = A[6];
        AA31 = A[7];
        AA02 = A[8];
        AA12 = A[9];
        AA22 = A[10];
        AA32 = A[11];
        AA03 = A[12];
        AA13 = A[13];
        AA23 = A[14];
        AA33 = A[15];

        perm0 = 0;
        perm1 = 1;
        perm2 = 2;
        perm3 = 3;
        
        /****************** iteration 0 ***********/

//#if USE_PIVOTING
        /* search pivot row */
        p = absOp (AA00);
        pvt = 0;
        t = absOp (AA10);
        if (t > p) { p = t;  pvt = 1; }
        t = absOp (AA20);
        if (t > p) { p = t;  pvt = 2; }
        t = absOp (AA30);
        if (t > p) { p = t;  pvt = 3; }
        
        /* swap pivot row with row 0 */
        if (pvt == 1) {
            tmp = AA00;  AA00 = AA10;  AA10 = tmp;
            tmp = AA01;  AA01 = AA11;  AA11 = tmp;
            tmp = AA02;  AA02 = AA12;  AA12 = tmp;
            tmp = AA03;  AA03 = AA13;  AA13 = tmp;
            /* update permutation vector based on row swap */
            i = perm0;  perm0 = perm1;  perm1 = i;
        }
        if (pvt == 2) {
            tmp = AA00;  AA00 = AA20;  AA20 = tmp;
            tmp = AA01;  AA01 = AA21;  AA21 = tmp;
            tmp = AA02;  AA02 = AA22;  AA22 = tmp;
            tmp = AA03;  AA03 = AA23;  AA23 = tmp;
            /* update permutation vector based on row swap */
            i = perm0;  perm0 = perm2;  perm2 = i;
        }
        if (pvt == 3) {
            tmp = AA00;  AA00 = AA30;  AA30 = tmp;
            tmp = AA01;  AA01 = AA31;  AA31 = tmp;            
            tmp = AA02;  AA02 = AA32;  AA32 = tmp;
            tmp = AA03;  AA03 = AA33;  AA33 = tmp;
            /* update permutation vector based on row swap */
            i = perm0;  perm0 = perm3;  perm3 = i;
        }
//#endif // USE_PIVOTING

        /* scale current row */
        tmp = rcpOp (AA00);
        icol0 = perm0;
        AA00 = tmp;
        AA01 = mulOp (tmp, AA01);
        AA02 = mulOp (tmp, AA02);
        AA03 = mulOp (tmp, AA03);

        /* eliminate above and below current row */
        tmp = AA10;
        AA10 = mulOp (negOp(tmp), AA00);
        AA11 = fmnaOp (tmp, AA01, AA11);
        AA12 = fmnaOp (tmp, AA02, AA12);
        AA13 = fmnaOp (tmp, AA03, AA13);

        tmp = AA20;
        AA20 = mulOp (negOp(tmp), AA00);
        AA21 = fmnaOp (tmp, AA01, AA21);
        AA22 = fmnaOp (tmp, AA02, AA22);
        AA23 = fmnaOp (tmp, AA03, AA23);

        tmp = AA30;
        AA30 = mulOp (negOp(tmp), AA00);
        AA31 = fmnaOp (tmp, AA01, AA31);
        AA32 = fmnaOp (tmp, AA02, AA32);
        AA33 = fmnaOp (tmp, AA03, AA33);

        
        /****************** iteration 1 ***********/

//#if USE_PIVOTING
        /* search pivot row */
        p = absOp (AA11);
        pvt = 1;
        t = absOp (AA21);
        if (t > p) { p = t;  pvt = 2; }
        t = absOp (AA31);
        if (t > p) { p = t;  pvt = 3; }

        /* swap pivot row with row 1 */
        if (pvt == 2) {
            tmp = AA10;   AA10 = AA20;   AA20 = tmp;
            tmp = AA11;   AA11 = AA21;   AA21 = tmp;
            tmp = AA12;   AA12 = AA22;   AA22 = tmp;
            tmp = AA13;   AA13 = AA23;   AA23 = tmp;
            /* update permutation vector based on row swap */
            i = perm1;   perm1 = perm2;   perm2 = i;
        }
        if (pvt == 3) {
            tmp = AA10;   AA10 = AA30;   AA30 = tmp;
            tmp = AA11;   AA11 = AA31;   AA31 = tmp;
            tmp = AA12;   AA12 = AA32;   AA32 = tmp;
            tmp = AA13;   AA13 = AA33;   AA33 = tmp;
            /* update permutation vector based on row swap */
            i = perm1;   perm1 = perm3;   perm3 = i;
        }
//#endif // USE_PIVOTING

        /* scale current row */
        tmp = rcpOp (AA11);
        icol1 = perm1;
        AA10 = mulOp (tmp, AA10);
        AA11 = tmp;
        AA12 = mulOp (tmp, AA12);
        AA13 = mulOp (tmp, AA13);

        /* eliminate above and below current row */
        tmp = AA01;
        AA00 = fmnaOp (tmp, AA10, AA00);
        AA01 = mulOp (negOp(tmp), AA11);
        AA02 = fmnaOp (tmp, AA12, AA02);
        AA03 = fmnaOp (tmp, AA13, AA03);
        
        tmp = AA21;
        AA20 = fmnaOp (tmp, AA10, AA20);
        AA21 = mulOp (negOp(tmp), AA11);
        AA22 = fmnaOp (tmp, AA12, AA22);
        AA23 = fmnaOp (tmp, AA13, AA23);
        
        tmp = AA31;
        AA30 = fmnaOp (tmp, AA10, AA30);
        AA31 = mulOp (negOp(tmp), AA11);
        AA32 = fmnaOp (tmp, AA12, AA32);
        AA33 = fmnaOp (tmp, AA13, AA33);
        
        /****************** iteration 2 ****************/

//#if USE_PIVOTING
        /* search pivot row */
        p = absOp (AA22);
        pvt = 2;
        t = absOp (AA32);
        if (t > p) { p = t;  pvt = 3; }

        /* swap pivot row with row 2 */
        if (pvt == 3) {
            tmp = AA20;   AA20 = AA30;   AA30 = tmp;
            tmp = AA21;   AA21 = AA31;   AA31 = tmp;
            tmp = AA22;   AA22 = AA32;   AA32 = tmp;
            tmp = AA23;   AA23 = AA33;   AA33 = tmp;
            /* update permutation vector based on row swap */
            i = perm2;   perm2 = perm3;   perm3 = i;
        }
//#endif // USE_PIVOTING

        /* scale current row */
        tmp = rcpOp (AA22);
        icol2 = perm2;
        AA20 = mulOp (tmp, AA20);
        AA21 = mulOp (tmp, AA21);
        AA22 = tmp;
        AA23 = mulOp (tmp, AA23);

        /* eliminate above and below current row */
        tmp = AA02;
        AA00 = fmnaOp (tmp, AA20, AA00);
        AA01 = fmnaOp (tmp, AA21, AA01); 
        AA02 = mulOp (negOp(tmp), AA22);
        AA03 = fmnaOp (tmp, AA23, AA03);

        tmp = AA12;
        AA10 = fmnaOp (tmp, AA20, AA10);
        AA11 = fmnaOp (tmp, AA21, AA11);
        AA12 = mulOp (negOp(tmp), AA22);
        AA13 = fmnaOp (tmp, AA23, AA13);

        tmp = AA32;
        AA30 = fmnaOp (tmp, AA20, AA30);
        AA31 = fmnaOp (tmp, AA21, AA31);
        AA32 = mulOp (negOp(tmp), AA22);
        AA33 = fmnaOp (tmp, AA23, AA33);

        /****************** iteration 3 ****************/

        /* scale current row */
        tmp = rcpOp (AA33);
        icol3 = perm3;
        AA30 = mulOp (tmp, AA30);
        AA31 = mulOp (tmp, AA31);
        AA32 = mulOp (tmp, AA32);
        AA33 = tmp;

        /* eliminate above and below current row */
        tmp = AA03;
        AA00 = fmnaOp (tmp, AA30, AA00);
        AA01 = fmnaOp (tmp, AA31, AA01);
        AA02 = fmnaOp (tmp, AA32, AA02);
        AA03 = mulOp (negOp(tmp), AA33);

        tmp = AA13;
        AA10 = fmnaOp (tmp, AA30, AA10);
        AA11 = fmnaOp (tmp, AA31, AA11);
        AA12 = fmnaOp (tmp, AA32, AA12);
        AA13 = mulOp (negOp(tmp), AA33);

        tmp = AA23;
        AA20 = fmnaOp (tmp, AA30, AA20);
        AA21 = fmnaOp (tmp, AA31, AA21);
        AA22 = fmnaOp (tmp, AA32, AA22);
        AA23 = mulOp (negOp(tmp), AA33);

        /* sort columns into the correct order */
        Ainv[0*4+icol0] = AA00;
        Ainv[1*4+icol0] = AA10;
        Ainv[2*4+icol0] = AA20;
        Ainv[3*4+icol0] = AA30;
        Ainv[0*4+icol1] = AA01;
        Ainv[1*4+icol1] = AA11;
        Ainv[2*4+icol1] = AA21;
        Ainv[3*4+icol1] = AA31;
        Ainv[0*4+icol2] = AA02;
        Ainv[1*4+icol2] = AA12;
        Ainv[2*4+icol2] = AA22;
        Ainv[3*4+icol2] = AA32;
        Ainv[0*4+icol3] = AA03;
        Ainv[1*4+icol3] = AA13;
        Ainv[2*4+icol3] = AA23;
        Ainv[3*4+icol3] = AA33;
    //}
}


__device__ float chi2_track(size_t const n_points, 
			float* const driftdist, float* const resolutions,
			float* const p1x, float* const p1y, float* const p1z,
			float* const deltapx, float* const deltapy, float* const deltapz,
			const float x0, const float y0, const float tx, const float ty)
{
	float dca;
	float chi2 = 0;
	for( size_t i=0; i<n_points; i++ ){
		dca = ( -deltapy[i]*(p1x[i]-x0) + deltapx[i]*(p1y[i]-y0) + p1z[i]*(tx*deltapy[i]-ty*deltapx[i]) ) / sqrtf( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*tx*ty*deltapy[i]*deltapx[i] );
		chi2+= ( driftdist[i] - dca ) * ( driftdist[i] - dca ) / resolutions[i] / resolutions[i];
	}
	return chi2;
}





#ifdef OLDCODE

__device__ void chi2_straight(size_t const n_points, 
                              float* const driftdist, float* const resolutions,
                              float* const p1x, float* const p1y, float* const p1z,
                              float* const deltapx, float* const deltapy, float* const deltapz,
			      float* output_parameters, float& chi2)
{
	float dca;
	chi2 = 0;
	for( int i=0; i<n_points; i++ ){
		//printf("%1.6f %1.6f, %1.6f %1.6f %1.6f, %1.6f %1.6f %1.6f \n", driftdist[i], resolutions[i], p1x[i], p1y[i], p1z[i], deltapx[i], deltapy[i], deltapz[i]);
		//printf("%1.6f %1.6f, %1.6f ; %1.6f %1.6f %1.6f (%1.6f) ; \n", p1x[i]-output_parameters[0], p1y[i]-output_parameters[1], p1z[i], output_parameters[3]*deltapz[i] - deltapy[i], deltapx[i] - output_parameters[2]*deltapz[i], output_parameters[2]*deltapy[i] - output_parameters[3]*deltapx[i], deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i]);
		//this is the simplified expression of the chi2 where we neglect deltapz
		//TODO: plug in the float expression
		dca = ( -deltapy[i]*(p1x[i]-output_parameters[0]) + deltapx[i]*(p1y[i]-output_parameters[1]) + p1z[i]*(output_parameters[2]*deltapy[i]-output_parameters[3]*deltapx[i]) ) / sqrtf( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i] );
		//printf(" dca = %1.6f \n", dca);
		chi2+= ( driftdist[i] - dca ) * ( driftdist[i] - dca ) / resolutions[i] / resolutions[i];
	}
}

__device__ void chi2_global(size_t const n_points, 
                            float* const driftdist, float* const resolutions,
                            float* const p1x, float* const p1y, float* const p1z,
                            float* const deltapx, float* const deltapy, float* const deltapz,
			    float* output_parameters, float& chi2)
{
	float dca;
	chi2 = 0;
	for( int i=0; i<n_points; i++ ){
	     	//this is the simplified expression of the chi2 where we neglect deltapz
		//TODO: plug in the float expression
		dca = ( -deltapy[i]*(p1x[i]-output_parameters[0]) + deltapx[i]*(p1y[i]-output_parameters[1]) + p1z[i]*(output_parameters[2]*deltapy[i]-output_parameters[3]*deltapx[i]) ) / sqrtf( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i] );
		chi2+= ( driftdist[i] - dca ) * ( driftdist[i] - dca ) / resolutions[i] / resolutions[i];
	}
}

__device__ void chisquare(size_t const n_points, 
                            float* const driftdist, float* const resolutions,
                            float* const p1x, float* const p1y, float* const p1z,
                            float* const deltapx, float* const deltapy, float* const deltapz,
			    size_t const nparam, float* output_parameters, float& chi2)
{
	if(nparam==4)chi2_straight(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, output_parameters, chi2);
	if(nparam==5)chi2_global(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, output_parameters, chi2);
}

__device__ void chi2_straight_res(size_t const n_points, 
                              float* const driftdist, float* const resolutions,
                              float* const p1x, float* const p1y, float* const p1z,
                              float* const deltapx, float* const deltapy, float* const deltapz,
			      float* output_parameters, float& chi2, float* values)
{
	chi2 = 0;
	for( int i=0; i<n_points; i++ ){
		//this is the simplified expression of the chi2 where we neglect deltapz
		//TODO: plug in the float expression
	}
}

__device__ void chi2_global_res(size_t const n_points, 
                            float* const driftdist, float* const resolutions,
                            float* const p1x, float* const p1y, float* const p1z,
                            float* const deltapx, float* const deltapy, float* const deltapz,
			    float* output_parameters, float& chi2, float* values)
{
	chi2 = 0;
	for( int i=0; i<n_points; i++ ){
	     	//this is the simplified expression of the chi2 where we neglect deltapz
		//TODO: plug in the float expression
		values[i] = driftdist[i] - ( -deltapy[i]*(p1x[i]-output_parameters[0]) + deltapx[i]*(p1y[i]-output_parameters[1]) + p1z[i]*(output_parameters[2]*deltapy[i]-output_parameters[3]*deltapx[i]) ) / sqrtf( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i] );
		chi2+= ( values[i] ) * ( values[i] ) / resolutions[i] / resolutions[i];
	}
}

__device__ void chisquare_res(size_t const n_points, 
                            float* const driftdist, float* const resolutions,
                            float* const p1x, float* const p1y, float* const p1z,
                            float* const deltapx, float* const deltapy, float* const deltapz,
			    size_t const nparam, float* output_parameters, float& chi2, float *values)
{
	if(nparam==4)chi2_straight_res(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, output_parameters, chi2, values);
	if(nparam==5)chi2_global_res(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, output_parameters, chi2, values);
}



// numerical calculation of derivative
__device__ void calc_derivatives_num(size_t const n_points,
				     float* const driftdist, float* const resolutions,
				     float* const p1x, float* const p1y, float* const p1z,
				     float* const deltapx, float* const deltapy, float* const deltapz,
				     size_t const nparam, float* output_parameters, float* output_parameters_steps, float const lambda, 
				     float* derivatives, //float* doublederivatives, 
				     float& chi2, float &chi2prev, float dca)
{
	//calculation of chi2 at the nominal x, 
	chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2prev);
	
	for(int i = 0; i<nparam; i++)
	{
		//lambda scales the existing step
		output_parameters[i]+=output_parameters_steps[i]*lambda;
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2);
		derivatives[i] = (chi2-chi2prev)/(output_parameters_steps[i]*lambda);   // (chi2(X+dX)-chi2(X))/dX
		
		//output_parameters[i]-=2*output_parameters_steps[i]*lambda;
		//chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2);
		//doublederivatives[i] = ((chi2prev-chi2)/(output_parameters_steps[i]*lambda)-derivatives[i])/output_parameters_steps[i]*lambda;
		output_parameters[i]-=output_parameters_steps[i]*lambda;
		
		if(derivatives[i]>0){
			output_parameters_steps[i] = -fabs(output_parameters_steps[i]);
		}else{
			output_parameters_steps[i] = fabs(output_parameters_steps[i]);
		}
		
		output_parameters[i]+=output_parameters_steps[i];
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2);
		//if at longer range the chi2 is higher, we shrink the step by a factor 2
		if(chi2-chi2prev>0){
			output_parameters_steps[i]*=0.5f;
		}
		output_parameters[i]-=output_parameters_steps[i]*lambda;
		
//output_parameters_steps[i] = -derivatives[i]/doublederivatives[i];
	}

}


__device__ void calc_corr_derivatives_num(size_t const n_points,
					  float* const driftdist, float* const resolutions,
					  float* const p1x, float* const p1y, float* const p1z,
					  float* const deltapx, float* const deltapy, float* const deltapz,
					  size_t const nparam, float* output_parameters, float* output_parameters_steps, float const lambda, 
					  float* derivatives, 
					  float& chi2, float &chi2prev, float dca)
{
	chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2prev);
	
	for(int i = 0; i<2; i++)
	{
		//lambda scales the existing step
		output_parameters[i]+=output_parameters_steps[i]*lambda;
		output_parameters[i+2]+=output_parameters_steps[i+2]*lambda;
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2);
		derivatives[i] = (chi2-chi2prev)/(fabs(output_parameters_steps[i])*lambda);   // (chi2(X+dX)-chi2(X))/dX
		output_parameters[i]-=output_parameters_steps[i]*lambda;
		output_parameters[i+2]-=output_parameters_steps[i+2]*lambda;


		if(derivatives[i]>0){
			output_parameters_steps[i] = -fabs(output_parameters_steps[i]);
			output_parameters_steps[i+2] = -fabs(output_parameters_steps[i+2]);
		}else{
			output_parameters_steps[i] = fabs(output_parameters_steps[i]);
			output_parameters_steps[i+2] = fabs(output_parameters_steps[i+2]);
			
		}

		output_parameters[i]+=output_parameters_steps[i];
		output_parameters[i+2]+=output_parameters_steps[i+2];
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2);
		//if at longer range the chi2 is higher, we shrink the step by a factor 2
		if(chi2-chi2prev>0){
			output_parameters_steps[i]*=0.5f;
			output_parameters_steps[i+2]*=0.5f;
		}
		output_parameters[i]-=output_parameters_steps[i];
		output_parameters[i+2]-=output_parameters_steps[i+2];
	}
}

#endif


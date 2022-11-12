#include "operations.h"

__device__ void matinv_2x2_matrix_per_thread (const REAL *A, REAL *Ainv)
{
    const int blkNum = blockIdx.y * gridDim.x + blockIdx.x;
    const int thrdNum = blkNum * blockDim.x + threadIdx.x;
    const int N = 2;
    int perm0, perm1;
    int icol0, icol1;
    REAL AA00, AA01; 
    REAL AA10, AA11;
    REAL tmp;
//#if USE_PIVOTING
    //typename config<T,arch>::absValType t;
    //typename config<T,arch>::absValType p;
    REAL t;
    REAL p;
    int i, pvt;
//#endif

    A    += thrdNum * N * N;
    Ainv += thrdNum * N * N;

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

__device__ void matinv_4x4_matrix_per_thread (const REAL *A, REAL *Ainv)
{
    //const int blkNum = blockIdx.y * gridDim.x + blockIdx.x;
    //const int thrdNum = blkNum * blockDim.x + threadIdx.x;
    //const int N = 4;
    int perm0, perm1, perm2, perm3;
    int icol0, icol1, icol2, icol3;
    REAL AA00, AA01, AA02, AA03; 
    REAL AA10, AA11, AA12, AA13;
    REAL AA20, AA21, AA22, AA23;
    REAL AA30, AA31, AA32, AA33;
    REAL tmp;

//#if USE_PIVOTING
    REAL t;
    REAL p;
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

__device__ void chi2_simplefit(size_t const n_points, REAL* const z_array, REAL* const x_array, REAL& tx, REAL& x0, REAL sum, REAL det, REAL sz, REAL sx, REAL szz, REAL sxx, REAL szx)
{
	sx = 0;
	sz = 0;
	szz = 0;
	szx = 0;
	sxx = 0;
	
	for(int i = 0; i<n_points; i++){
		sum+=1.;
		sz+=z_array[i];
		sx+=x_array[i];
		szz+=z_array[i]*z_array[i];
		sxx+=x_array[i]*x_array[i];
		szx+=x_array[i]*z_array[i];
	}
	det = sum*szz - sz*sz; 	
	if(fabs(det)<1.e-20){
		tx = 0.;
		x0 = 0.;
		return;
	}
	tx = (sum*szx - sz*sx)/det;
	x0 = (sum*sxx - szx*sx)/det;

}

__device__ void chi2_straight(size_t const n_points, 
                              REAL* const driftdist, REAL* const resolutions,
                              REAL* const p1x, REAL* const p1y, REAL* const p1z,
                              REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
			      REAL* output_parameters, REAL& chi2, REAL dca)
{
	chi2 = 0;
	for( int i=0; i<n_points; i++ ){
		//printf("%1.6f %1.6f, %1.6f %1.6f %1.6f, %1.6f %1.6f %1.6f \n", driftdist[i], resolutions[i], p1x[i], p1y[i], p1z[i], deltapx[i], deltapy[i], deltapz[i]);
		//printf("%1.6f %1.6f, %1.6f ; %1.6f %1.6f %1.6f (%1.6f) ; \n", p1x[i]-output_parameters[0], p1y[i]-output_parameters[1], p1z[i], output_parameters[3]*deltapz[i] - deltapy[i], deltapx[i] - output_parameters[2]*deltapz[i], output_parameters[2]*deltapy[i] - output_parameters[3]*deltapx[i], deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i]);
		//this is the simplified expression of the chi2 where we neglect deltapz
		//TODO: plug in the real expression
		dca = ( -deltapy[i]*(p1x[i]-output_parameters[0]) + deltapx[i]*(p1y[i]-output_parameters[1]) + p1z[i]*(output_parameters[2]*deltapy[i]-output_parameters[3]*deltapx[i]) ) / sqrtf( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i] );
		//printf(" dca = %1.6f \n", dca);
		chi2+= ( driftdist[i] - dca ) * ( driftdist[i] - dca ) / resolutions[i] / resolutions[i];
	}
}

__device__ void chi2_global(size_t const n_points, 
                            REAL* const driftdist, REAL* const resolutions,
                            REAL* const p1x, REAL* const p1y, REAL* const p1z,
                            REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
			    REAL* output_parameters, REAL& chi2, REAL dca)
{
	chi2 = 0;
	for( int i=0; i<n_points; i++ ){
	     	//this is the simplified expression of the chi2 where we neglect deltapz
		//TODO: plug in the real expression
		dca = ( -deltapy[i]*(p1x[i]-output_parameters[0]) + deltapx[i]*(p1y[i]-output_parameters[1]) + p1z[i]*(output_parameters[2]*deltapy[i]-output_parameters[3]*deltapx[i]) ) / sqrtf( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i] );
		chi2+= ( driftdist[i] - dca ) * ( driftdist[i] - dca ) / resolutions[i] / resolutions[i];
	}
}

__device__ void chisquare(size_t const n_points, 
                            REAL* const driftdist, REAL* const resolutions,
                            REAL* const p1x, REAL* const p1y, REAL* const p1z,
                            REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
			    size_t const nparam, REAL* output_parameters, REAL& chi2, REAL dca)
{
	if(nparam==4)chi2_straight(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, output_parameters, chi2, dca);
	if(nparam==5)chi2_global(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, output_parameters, chi2, dca);
}

__device__ void chi2_straight_res(size_t const n_points, 
                              REAL* const driftdist, REAL* const resolutions,
                              REAL* const p1x, REAL* const p1y, REAL* const p1z,
                              REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
			      REAL* output_parameters, REAL& chi2, REAL* values)
{
	chi2 = 0;
	for( int i=0; i<n_points; i++ ){
		//this is the simplified expression of the chi2 where we neglect deltapz
		//TODO: plug in the real expression
	}
}

__device__ void chi2_global_res(size_t const n_points, 
                            REAL* const driftdist, REAL* const resolutions,
                            REAL* const p1x, REAL* const p1y, REAL* const p1z,
                            REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
			    REAL* output_parameters, REAL& chi2, REAL* values)
{
	chi2 = 0;
	for( int i=0; i<n_points; i++ ){
	     	//this is the simplified expression of the chi2 where we neglect deltapz
		//TODO: plug in the real expression
		values[i] = driftdist[i] - ( -deltapy[i]*(p1x[i]-output_parameters[0]) + deltapx[i]*(p1y[i]-output_parameters[1]) + p1z[i]*(output_parameters[2]*deltapy[i]-output_parameters[3]*deltapx[i]) ) / sqrtf( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i] );
		chi2+= ( values[i] ) * ( values[i] ) / resolutions[i] / resolutions[i];
	}
}

__device__ void chisquare_res(size_t const n_points, 
                            REAL* const driftdist, REAL* const resolutions,
                            REAL* const p1x, REAL* const p1y, REAL* const p1z,
                            REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
			    size_t const nparam, REAL* output_parameters, REAL& chi2, REAL *values)
{
	if(nparam==4)chi2_straight_res(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, output_parameters, chi2, values);
	if(nparam==5)chi2_global_res(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, output_parameters, chi2, values);
}



// numerical calculation of derivative
__device__ void calc_derivatives_num(size_t const n_points,
				     REAL* const driftdist, REAL* const resolutions,
				     REAL* const p1x, REAL* const p1y, REAL* const p1z,
				     REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
				     size_t const nparam, REAL* output_parameters, REAL* output_parameters_steps, REAL const lambda, 
				     REAL* derivatives, //REAL* doublederivatives, 
				     REAL& chi2, REAL &chi2prev, REAL dca)
{
	//calculation of chi2 at the nominal x, 
	chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2prev, dca);
	
	for(int i = 0; i<nparam; i++)
	{
		//lambda scales the existing step
		output_parameters[i]+=output_parameters_steps[i]*lambda;
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2, dca);
		derivatives[i] = (chi2-chi2prev)/(output_parameters_steps[i]*lambda);   // (chi2(X+dX)-chi2(X))/dX
		
		//output_parameters[i]-=2*output_parameters_steps[i]*lambda;
		//chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2, dca);
		//doublederivatives[i] = ((chi2prev-chi2)/(output_parameters_steps[i]*lambda)-derivatives[i])/output_parameters_steps[i]*lambda;
		output_parameters[i]-=output_parameters_steps[i]*lambda;
		
		if(derivatives[i]>0){
			output_parameters_steps[i] = -fabs(output_parameters_steps[i]);
		}else{
			output_parameters_steps[i] = fabs(output_parameters_steps[i]);
		}
		
		output_parameters[i]+=output_parameters_steps[i];
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2, dca);
		//if at longer range the chi2 is higher, we shrink the step by a factor 2
		if(chi2-chi2prev>0){
			output_parameters_steps[i]*=0.5f;
		}
		output_parameters[i]-=output_parameters_steps[i]*lambda;
		
//output_parameters_steps[i] = -derivatives[i]/doublederivatives[i];
	}

}


__device__ void calc_corr_derivatives_num(size_t const n_points,
					  REAL* const driftdist, REAL* const resolutions,
					  REAL* const p1x, REAL* const p1y, REAL* const p1z,
					  REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
					  size_t const nparam, REAL* output_parameters, REAL* output_parameters_steps, REAL const lambda, 
					  REAL* derivatives, 
					  REAL& chi2, REAL &chi2prev, REAL dca)
{
	chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2prev, dca);
	
	for(int i = 0; i<2; i++)
	{
		//lambda scales the existing step
		output_parameters[i]+=output_parameters_steps[i]*lambda;
		output_parameters[i+2]+=output_parameters_steps[i+2]*lambda;
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2, dca);
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
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2, dca);
		//if at longer range the chi2 is higher, we shrink the step by a factor 2
		if(chi2-chi2prev>0){
			output_parameters_steps[i]*=0.5f;
			output_parameters_steps[i+2]*=0.5f;
		}
		output_parameters[i]-=output_parameters_steps[i];
		output_parameters[i+2]-=output_parameters_steps[i+2];
	}
}


__device__ void calc_val_derivatives(size_t const n_points,
				     REAL* const driftdist, REAL* const resolutions,
				     REAL* const p1x, REAL* const p1y, REAL* const p1z,
				     REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
				     size_t const nparam, REAL* const output_parameters, 
				     REAL* values, REAL* derivatives, REAL& chi2)
{
	REAL Den2, Den;
	chi2 = 0;
	for(int i = 0; i<n_points; i++){
		values[i] = driftdist[i] - ( -deltapy[i]*(p1x[i]-output_parameters[0]) + deltapx[i]*(p1y[i]-output_parameters[1]) + p1z[i]*(output_parameters[2]*deltapy[i]-output_parameters[3]*deltapx[i]) ) / sqrtf( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i] );
		chi2+= ( values[i] ) * ( values[i] ) / resolutions[i] / resolutions[i];

	     	Den2 = deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2 * ( deltapx[i]*deltapy[i]*output_parameters[2]*output_parameters[3]);
	     	Den = sqrtf(Den2);
		
		//dchi2/dx0:
		derivatives[0*n_points] = -deltapy[i]/Den;
		//dchi2/dy0:
		derivatives[1*n_points] = deltapx[i]/Den;
		//dchi2/dtx:
		derivatives[2*n_points] = -deltapy[i]*p1z[i]/Den+values[i]*output_parameters[3]*deltapx[i]*deltapy[i]/Den2;
		//dchi2/dty:
		derivatives[3*n_points] = deltapx[i]*p1z[i]/Den+values[i]*output_parameters[2]*deltapx[i]*deltapy[i]/Den2;
	}
}


__device__ void calc_gradients(size_t const n_points, 
			       REAL* const resolutions, REAL* const driftdist,
			       int const n_parameters,
			       REAL* values, REAL* derivatives, REAL* gradients)
{
	for(int i = 0; i<n_points; i++){
		for(int j = 0; j< n_parameters; j++){
			gradients[j*n_points] = derivatives[j*n_points]*(0.0-values[i])*resolutions[i];
		}
	}
}


__device__ void calc_hessians(size_t const n_points,
                              REAL* const resolutions, 
			      int const n_parameters,
			      REAL* values, REAL* derivatives, REAL* hessians)
{
	for(int i = 0; i<n_points; i++){
		//j: row; k: col; 
		for(int j = 0; j< n_parameters; j++){
			for(int k = 0; k<n_parameters; k++){
                		hessians[j*n_parameters+k]+= derivatives[j*n_points]*derivatives[k*n_points]*resolutions[i];
			}
		}
        }
}




__device__ void get_straighttrack_fixedpoint(size_t const n_points,
					     REAL* driftdist, REAL* resolutions,
					     REAL* p1x, REAL* p1y, REAL* p1z,
					     REAL* deltapx, REAL* deltapy, REAL* deltapz,
					     REAL *A, REAL* Ainv, REAL* B,
					     REAL* output_parameters, REAL* output_parameters_errors, 
					     REAL* fixed_point, REAL& chi2, REAL dca)
{
	REAL Den2, Den;
	// slopes are set to ZERO on purpose
	for( int i=0; i<n_points; i++ ){
		Den2 = deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i];
		Den = sqrtf(Den2);
		
		//tx, ty, treated as constants...
		B[0] += (driftdist[i]*deltapy[i]*Den + deltapy[i]*deltapy[i]*p1x[i] - deltapx[i]*deltapy[i]*p1y[i])/(resolutions[i]*resolutions[i]*Den2);
		B[1] += (-driftdist[i]*deltapx[i]*Den - deltapx[i]*deltapy[i]*p1x[i] + deltapx[i]*deltapx[i]*p1y[i])/(resolutions[i]*resolutions[i]*Den2);
		
		A[0*2+0] += deltapy[i]*deltapy[i]/(resolutions[i]*resolutions[i]*Den2);
		A[0*2+1] += -deltapx[i]*deltapy[i]/(resolutions[i]*resolutions[i]*Den2);
		
		A[1*2+0] += -deltapx[i]*deltapy[i]/(resolutions[i]*resolutions[i]*Den2);
		A[1*2+1] += deltapx[i]*deltapx[i]/(resolutions[i]*resolutions[i]*Den2);
		fixed_point[2] += p1z[i]/n_points;//the z of the fixed point is just the average of the tracks
	}

	matinv_2x2_matrix_per_thread(A, Ainv);
	
	for(int j = 0; j<2; j++){//row
		output_parameters[j] = 0.0;
		output_parameters_errors[j] = sqrtf(fabs(A[j*2+j]));
		for(int k = 0; k<2; k++){//column
			output_parameters[j]+= Ainv[j*2+k]*B[k];
		}
		fixed_point[j] = output_parameters[j];
	}
	chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, 4, output_parameters, chi2, dca);
}


#include "operations.h"

const double TX_MAX = 0.15;
const double TY_MAX = 0.1;
const double X0_MAX = 150;
const double Y0_MAX = 50;
const double INVP_MAX = 0.2;
const double INVP_MIN = 0.01;


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


/*
__device__ void prefit_straight_tracks_parameters(size_t const n_points,
						  REAL* driftdist, REAL* resolutions,
						  REAL* p1x, REAL* p1y, REAL* p1z,
						  REAL* deltapx, REAL* deltapy, REAL* deltapz,
						  REAL *A, REAL* Ainv, REAL* B,
						  REAL* output_parameters, REAL* output_parameters_errors, REAL& chi2)
{
	REAL Den2, Den;

	for( int i=0; i<n_points; i++ ){
		Den2 = deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2 * ( deltapx[i]*deltapy[i]*output_parameters[2]*output_parameters[3]);
		Den = sqrtf(Den2);
		
		//tx, ty, treated as constants...
		B[0] += (driftdist[i]*deltapy[i]*Den + deltapy[i]*deltapy[i]*p1x[i] - deltapx[i]*deltapy[i]*p1y[i] +
		         deltapy[i]*deltapy[i]*p1z[i]*output_parameters[2] - deltapx[i]*deltapy[i]*p1z[i]*output_parameters[3])/(resolutions[i]*resolutions[i]*Den2);
		B[1] += (-driftdist[i]*deltapx[i]*Den - deltapx[i]*deltapy[i]*p1x[i] + deltapx[i]*deltapx[i]*p1y[i] - 
		     	 deltapx[i]*deltapy[i]*p1z[i]*output_parameters[2] + deltapx[i]*deltapx[i]*p1z[i]*output_parameters[3])/(resolutions[i]*resolutions[i]*Den2);
		
		A[0*2+0] += deltapy[i]*deltapy[i]/(resolutions[i]*resolutions[i]*Den2);
		A[0*2+1] += -deltapx[i]*deltapy[i]/(resolutions[i]*resolutions[i]*Den2);
		
		A[1*2+0] += -deltapx[i]*deltapy[i]/(resolutions[i]*resolutions[i]*Den2);
		A[1*2+1] += deltapx[i]*deltapx[i]/(resolutions[i]*resolutions[i]*Den2);
	}	

	matinv_2x2_matrix_per_thread(A, Ainv);
	
	for(int j = 0; j<2; j++){//row
		output_parameters[j] = 0.0;
		output_parameters_errors[j] = sqrtf(fabs(A[j*2+j]));
		for(int k = 0; k<2; k++){//column
			output_parameters[j]+= Ainv[j*2+k]*B[k];
		}
	}
	
	printf("x0 = %1.6f +- %1.6f; y0 = %1.6f +- %1.6f; \n", output_parameters[0], output_parameters_errors[0], output_parameters[1], output_parameters_errors[1]);
	
	//chi2 = 0;
	//REAL dca;
	//for( int i=0; i<n_points; i++ ){
	//	dca = ( -deltapy[i]*(p1x[i]-output_parameters[0]) + deltapx[i]*(p1y[i]-output_parameters[1]) + p1z[i]*(output_parameters[2]*deltapy[i]-output_parameters[3]*deltapx[i]) ) / ( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i] );
	//	chi2+= ( driftdist[i] - dca ) * ( driftdist[i] - dca ) / resolutions[i] / resolutions[i];
	//}
	
	// now we should have "starting" values for x0, y0, let's evaluate the slopes
	for(int j = 0; j<5; j++){
		B[j] = 0;
        	for(int k = 0; k<5; k++){
        		A[j*5+k] = 0;
			Ainv[j*5+k] = 0;
		}
	}
	
	for( int i=0; i<n_points; i++ ){
		Den2 = deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2 * ( deltapx[i]*deltapy[i]*output_parameters[2]*output_parameters[3]);
		Den = sqrtf(Den2);
		
		//tx, ty, treated as constants...
		B[0] += (driftdist[i]*deltapy[i]*Den + deltapy[i]*deltapy[i]*p1x[i] - deltapx[i]*deltapy[i]*p1y[i] + 
		         deltapy[i]*deltapy[i]*p1z[i]*output_parameters[0] - deltapx[i]*deltapy[i]*p1z[i]*output_parameters[1])/(resolutions[i]*resolutions[i]*Den2);
		B[1] += (-driftdist[i]*deltapx[i]*Den - deltapx[i]*deltapy[i]*p1x[i] + deltapx[i]*deltapx[i]*p1y[i] - 
		     	 deltapx[i]*deltapy[i]*p1z[i]*output_parameters[0] + deltapx[i]*deltapx[i]*p1z[i]*output_parameters[1])/(resolutions[i]*resolutions[i]*Den2);
		
		A[0*2+0] += deltapy[i]*deltapy[i]*p1z[i]*p1z[i]/(resolutions[i]*resolutions[i]*Den2);
		A[0*2+1] += -(deltapx[i]*deltapy[i]*p1z[i]*p1z[i] + driftdist[i]*(deltapx[i]*deltapy[i]*deltapy[i]*p1x[i] - deltapx[i]*deltapx[i]*deltapy[i]*p1y[i])/Den -
			     deltapx[i]*deltapy[i]*(deltapy[i]*deltapy[i]*p1x[i]*p1x[i] - deltapx[i]*deltapx[i]*p1y[i]*p1y[i] + 2*deltapx[i]*deltapy[i]*p1x[i]*p1y[i] )/Den2 )/(resolutions[i]*resolutions[i]*Den2);
		
		A[1*2+0] += -(deltapx[i]*deltapy[i]*p1z[i]*p1z[i] + driftdist[i]*(deltapx[i]*deltapy[i]*deltapy[i]*p1x[i] - deltapx[i]*deltapx[i]*deltapy[i]*p1y[i])/Den -
			     deltapx[i]*deltapy[i]*(deltapy[i]*deltapy[i]*p1x[i]*p1x[i] - deltapx[i]*deltapx[i]*p1y[i]*p1y[i] + 2*deltapx[i]*deltapy[i]*p1x[i]*p1y[i] )/Den2 )/(resolutions[i]*resolutions[i]*Den2);
		A[1*2+1] += deltapx[i]*deltapx[i]*p1z[i]*p1z[i]/(resolutions[i]*resolutions[i]*Den2);
	}	
	
	matinv_2x2_matrix_per_thread(A, Ainv);

	for(int j = 0; j<2; j++){//row
		output_parameters[j+2] = 0.0;
		output_parameters_errors[j+2] = sqrtf(fabs(A[j*2+j]));
		for(int k = 0; k<2; k++){//column
			output_parameters[j+2]+= Ainv[j*2+k]*B[k];
		}
	}
	
        for(int j = 0; j<5; j++){
                B[j] = 0;
                for(int k = 0; k<5; k++){
                        A[j*5+k] = 0;
                        Ainv[j*5+k] = 0;
                }
        }
	
	printf("tx = %1.6f +- %1.6f; ty = %1.6f +- %1.6f; \n", output_parameters[2], output_parameters_errors[2], output_parameters[3], output_parameters_errors[3]);
	
}
*/

//This function allows to calculate the chi2 for a set of given parameters
__device__ REAL chisquare(size_t const n_points, 
			  REAL* driftdist, REAL* resolutions, 
			  REAL* p1x, REAL* p1y, REAL* p1z, 
			  REAL* deltapx, REAL* deltapy, REAL* deltapz, 
			  REAL* output_parameters)
{
	REAL dca, chi2 = 0;
	for( int i=0; i<n_points; i++ ){
		dca = ( -deltapy[i]*(p1x[i]-output_parameters[0]) + deltapx[i]*(p1y[i]-output_parameters[1]) + p1z[i]*(output_parameters[2]*deltapy[i]-output_parameters[3]*deltapx[i]) ) / ( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i] );
		chi2+= ( driftdist[i] - dca ) * ( driftdist[i] - dca ) / resolutions[i] / resolutions[i];
	}
	return chi2;
}


// This function gets a first estimate of the parameters, by calculating the chi^2 for ( {-X0_MAX, 0, X0_MAX}, {-Y0_MAX, 0, Y0_MAX}, {-TX_MAX, 0, YX_MAX}, {-TY_MAX, 0, TY_MAX} ) 
// and keep the best chi^2 for those conditions... That is very rough, but hopefully good enough.
__device__ void get_track_initial_parameters(size_t const n_points,
					     REAL* driftdist, REAL* resolutions,
					     REAL* p1x, REAL* p1y, REAL* p1z,
					     REAL* deltapx, REAL* deltapy, REAL* deltapz,
					     REAL* output_parameters)
{
	REAL chi2min = 10000000.0f;
	REAL chi2;
	int i_min, j_min, k_min, l_min;//, m_min;
	for(int i = 0; i<4; i++){
		output_parameters[0] = X0_MAX*(i*2-3)/3.f;
	printf("%1.1f \n", output_parameters[0]);
		for(int j = 0; j<4; j++){
			output_parameters[1] = Y0_MAX*(j*2-3)/3.f;
			for(int k = 0; k<4; k++){
				output_parameters[2] = TX_MAX*(k*2-3)/3.f;
				for(int l = 0; l<4; l++){
					//for(int m = 0; m<4; m++){}
					output_parameters[3] = TY_MAX*(l*2-3)/3.f;
					chi2 = chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, output_parameters);
					if(chi2<chi2min){
						chi2min = chi2;
						i_min = i;
						j_min = j;
						k_min = k;
						l_min = l;
					}
				}
			}
		}
	}
	//printf("%1.6f \n", chi2min);

	output_parameters[0] = X0_MAX*(i_min*2-3)/3.f;
	output_parameters[1] = Y0_MAX*(j_min*2-3)/3.f;
	output_parameters[2] = TX_MAX*(k_min*2-3)/3.f;
	output_parameters[3] = TY_MAX*(l_min*2-3)/3.f;
}



__device__ void straight_track_residual_minimizer(size_t const n_points, 
	   					  REAL* driftdist, REAL* resolutions, 
						  REAL* p1x, REAL* p1y, REAL* p1z, 
						  REAL* deltapx, REAL* deltapy, REAL* deltapz, 
						  REAL *A, REAL* Ainv, REAL* B, 
						  REAL* output_parameters, REAL* output_parameters_errors, REAL& chi2)
{
	// start minimizing chi2 for x0, y0 only 
	for(int j = 0; j<5; j++){
		B[j] = 0;
		for(int k = 0; k<5; k++){
			A[j*5+k] = 0;
			Ainv[j*5+k] = 0;
		}
	}
	
	//bool do_prefit = true;	
	//for(int k = 0; k<4; k++){
	//	//printf("%1.6f +- %1.6f \n", output_parameters[k], output_parameters_errors[k]);
	//	// if one of the parameters has non zero value, then we probably don't need to do the "prefit"
	//	if(output_parameters[k]!=0 || output_parameters_errors[k]!=0)do_prefit = false;
	//}
	
	// "prefit" function... 
	//if(do_prefit){
	//	prefit_straight_tracks_parameters(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, A, Ainv, B, output_parameters, output_parameters_errors, chi2);
	//}
	
	//printf("x0 %1.6f, y0 %1.6f, tx %1.6f, ty %1.6f\n", output_parameters[0], output_parameters[1], output_parameters[2], output_parameters[3]);
	
	REAL Den2, Den;
	
	for( int i=0; i<n_points; i++ ){
		//printf("%1.6f %1.6f, %1.6f %1.6f %1.6f, %1.6f %1.6f %1.6f \n", driftdist[i], resolutions[i], p1x[i], p1y[i], p1z[i], deltapx[i], deltapy[i], deltapz[i]);
		
	     	Den2 = deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2 * ( deltapx[i]*deltapy[i]*output_parameters[2]*output_parameters[3]);
	     	//Den2 = deltapy[i]*deltapy[i] + output_parameters[2]*output_parameters[2]*(deltapy[i]*deltapy[i]+deltapz[i]*deltapz[i]) 
		//     + deltapx[i]*deltapx[i] + output_parameters[3]*output_parameters[3]*(deltapx[i]*deltapx[i]+deltapz[i]*deltapz[i]) 
		//     - 2 * ( deltapx[i]*deltapz[i]*output_parameters[2] + deltapy[i]*deltapz[i]*output_parameters[3] + deltapx[i]*deltapy[i]*output_parameters[2]*output_parameters[3]);
	     	Den = sqrtf(Den2);
		
		B[0] += (driftdist[i]*deltapy[i]*Den + deltapy[i]*deltapy[i]*p1x[i] - deltapx[i]*deltapy[i]*p1y[i])/(resolutions[i]*resolutions[i]*Den2);
		B[1] += (-driftdist[i]*deltapx[i]*Den - deltapx[i]*deltapy[i]*p1x[i] + deltapx[i]*deltapx[i]*p1y[i])/(resolutions[i]*resolutions[i]*Den2);
		B[2] += (driftdist[i]*deltapy[i]*p1z[i]*Den + deltapy[i]*deltapy[i]*p1x[i]*p1z[i] - deltapx[i]*deltapy[i]*p1x[i]*p1y[i])/(resolutions[i]*resolutions[i]*Den2);
		B[3] += (-driftdist[i]*deltapx[i]*p1z[i]*Den - deltapx[i]*deltapy[i]*p1x[i]*p1z[i] + deltapx[i]*deltapx[i]*p1y[i]*p1z[i])/(resolutions[i]*resolutions[i]*Den2);
		
		// most complicated parameters.
		//A[0*4+0] += (deltapy[i]*deltapy[i] - 2*deltapy[i]*deltapz[i]*output_parameters[3] + deltapz[i]*deltapz[i]*output_parameters[3]*output_parameters[3])/(resolutions[i]*resolutions[i]*Den2);
		//A[0*4+1] += (-deltapx[i]*deltapy[i] + deltapy[i]*deltapz[i]*output_parameters[2] + deltapx[i]*deltapz[i]*output_parameters[3] - deltapz[i]*deltapz[i]*output_parameters[2]*output_parameters[3] )/(resolutions[i]*resolutions[i]*Den2);
		//A[0*4+2] += (-deltapy[i]*deltapz[i]*p1y[i] + deltapy[i]*deltapy[i]*p1z[i] + (deltapz[i]*deltapz[i]*p1y[i]-deltapy[i]*deltapz[i]*p1z[i])*output_parameters[3])/(resolutions[i]*resolutions[i]*Den2);
		//A[1*4+3] += (driftdist[i]*deltapz[i]*Den +2*deltapy[i]*deltapz[i]*p1x[i] - deltapx[i]*deltapz[i]*p1y[i] - deltapx[i]*deltapy[i]*p1z[i] - (deltapz[i]*deltapz[i]*p1x[i]-deltapx[i]*deltapz[i]*p1z[i])*output_parameters[3])/(resolutions[i]*resolutions[i]*Den2);
				
		//A[1*4+0] += (-deltapx[i]*deltapy[i] + deltapy[i]*deltapz[i]*output_parameters[2] + deltapx[i]*deltapz[i]*output_parameters[3] - deltapz[i]*deltapz[i]*output_parameters[2]*output_parameters[3] )/(resolutions[i]*resolutions[i]*Den2);
		//A[1*4+1] += (deltapx[i]*deltapy[i] - 2*deltapx[i]*deltapz[i]*output_parameters[2] + deltapz[i]*deltapz[i]*output_parameters[2]*output_parameters[2])/(resolutions[i]*resolutions[i]*Den2);
		//A[1*4+2] += (-driftdist[i]*deltapx[i]*p1z[i]*Den - deltapy[i]*deltapz[i]*p1z[i] + 2*deltapx[i]*deltapz[i]*p1y[i] - deltapx[i]*deltapy[i]*p1z[i] + (- deltapz[i]*deltapz[i]*p1y[i] + deltapy[i]*deltapz[i]*p1z[i])*output_parameters[2] + (deltapz[i]*deltapz[i]*p1x[i] + deltapx[i]*deltapz[i]*p1z[i])*output_parameters[3])/(resolutions[i]*resolutions[i]*Den2);
		//A[1*4+3] += (- deltapx[i]*deltapz[i]*p1x[i] + 2*deltapx[i]*deltapx[i]*p1z[i])/(resolutions[i]*resolutions[i]*Den2);
		
		// first index: row; second index: col;
		// (consistent with what is used in the matrix inversion routine) 
		A[0*4+0] += deltapy[i]*deltapy[i]/(resolutions[i]*resolutions[i]*Den2);
		A[0*4+1] += -deltapx[i]*deltapy[i]/(resolutions[i]*resolutions[i]*Den2);
		A[0*4+2] += p1z[i]*deltapy[i]*deltapy[i]/(resolutions[i]*resolutions[i]*Den2);
		A[0*4+3] += -p1z[i]*deltapx[i]*deltapy[i]/(resolutions[i]*resolutions[i]*Den2);
		
		A[1*4+0] += -deltapx[i]*deltapy[i]/(resolutions[i]*resolutions[i]*Den2);
		A[1*4+1] += deltapx[i]*deltapx[i]/(resolutions[i]*resolutions[i]*Den2);
		A[1*4+2] += -p1z[i]*deltapx[i]*deltapy[i]/(resolutions[i]*resolutions[i]*Den2);
		A[1*4+3] += p1z[i]*deltapx[i]*deltapx[i]/(resolutions[i]*resolutions[i]*Den2);
		
		// bilinear terms for dchi2/dt(x, y)  
		// +driftdist[i]*deltapx[i]*deltapy[i]*deltapy[i]/Den - 2*deltapx[i]*deltapy[i]*(deltapy[i]*deltapy[i]*p1x[i]-deltapx[i]*p1y[i])/Den2 //x0 t(y, x)
		// -driftdist[i]*deltapx[i]*deltapx[i]*deltapy[i]/Den - 2*deltapx[i]*deltapy[i]*(deltapx[i]*deltapx[i]*p1y[i]-deltapy[i]*p1x[i])/Den2 //y0 t(y, x)
		// +driftdist[i]*deltapx[i]*deltapy[i]*deltapy[i]*p1z[i]/Den - 2*deltapx[i]*deltapy[i]*deltapy[i]*p1z[i]*(deltapx[i]*p1y[i]-deltapy[i]*p1x[i])/Den2 //tx t(y, x)
		// -driftdist[i]*deltapx[i]*deltapx[i]*deltapy[i]*p1z[i]/Den + 2*deltapx[i]*deltapx[i]*deltapy[i]*p1z[i]*(deltapx[i]*p1y[i]-deltapy[i]*p1x[i])/Den2 //ty t(y, x)
		// trilinear terms
		// +deltapx[i]*deltapy[i]*deltapy[i]*deltapy[i]/Den2 //x0^2 t(y, x)
		// +deltapx[i]*deltapx[i]*deltapx[i]*deltapy[i]/Den2 //y0^2 t(y, x)
		// +deltapx[i]*deltapy[i]*deltapy[i]*deltapy[i]*p1z[i]*p1z[i]/Den2  // tx^2 t(y, x)
		// +deltapx[i]*deltapx[i]*deltapx[i]*deltapy[i]*p1z[i]*p1z[i]/Den2  // ty^2 t(y, x)
		// -2*deltapx[i]*deltapy[i]*deltapy[i]*deltapy[i]/Den2 //x0 y0 t(y, x)
		// +2*deltapx[i]*deltapy[i]*deltapy[i]*deltapy[i]*p1z[i]/Den2 //x0 tx t(y, x)
		// -2*deltapx[i]*deltapx[i]*deltapy[i]*deltapy[i]*p1z[i]/Den2 //x0 ty t(y, x)
		// +2*deltapx[i]*deltapx[i]*deltapy[i]*deltapy[i]*p1z[i]/Den2 //y0 tx t(y, x)
		// -2*deltapx[i]*deltapx[i]*deltapx[i]*deltapy[i]*p1z[i]/Den2 //y0 ty t(y, x)
		// -2*deltapx[i]*deltapx[i]*deltapy[i]*deltapy[i]*p1z[i]*p1z[i]/Den2 //tx ty t(y, x)
			
		A[2*4+0] += (deltapy[i]*deltapy[i]*p1z[i])/(resolutions[i]*resolutions[i]*Den2);
		A[2*4+1] += (-deltapx[i]*deltapy[i]*p1z[i])/(resolutions[i]*resolutions[i]*Den2);
		A[2*4+2] += (deltapy[i]*deltapy[i]*p1z[i]*p1z[i])/(resolutions[i]*resolutions[i]*Den2);
		A[2*4+3] += (deltapx[i]*deltapy[i]*p1z[i]*p1z[i] + driftdist[i]*(deltapx[i]*deltapy[i]*deltapy[i]*p1x[i] - deltapx[i]*deltapx[i]*deltapy[i]*p1y[i])/Den -
			     deltapx[i]*deltapy[i]*(deltapy[i]*deltapy[i]*p1x[i]*p1x[i] - deltapx[i]*deltapx[i]*p1y[i]*p1y[i] + 2*deltapx[i]*deltapy[i]*p1x[i]*p1y[i])/Den2
//			     +output_parameters[0]*(+driftdist[i]*deltapx[i]*deltapy[i]*deltapy[i]/Den - 2*deltapx[i]*deltapy[i]*(deltapy[i]*deltapy[i]*p1x[i]-deltapx[i]*p1y[i])/Den2)
//			     +output_parameters[1]*(-driftdist[i]*deltapx[i]*deltapx[i]*deltapy[i]/Den - 2*deltapx[i]*deltapy[i]*(deltapx[i]*deltapx[i]*p1y[i]-deltapy[i]*p1x[i])/Den2) 
//			     +output_parameters[2]*(+driftdist[i]*deltapx[i]*deltapy[i]*deltapy[i]*p1z[i]/Den - 2*deltapx[i]*deltapy[i]*deltapy[i]*p1z[i]*(deltapx[i]*p1y[i]-deltapy[i]*p1x[i])/Den2)
//			     +output_parameters[3]*(-driftdist[i]*deltapx[i]*deltapx[i]*deltapy[i]*p1z[i]/Den + 2*deltapx[i]*deltapx[i]*deltapy[i]*p1z[i]*(deltapx[i]*p1y[i]-deltapy[i]*p1x[i])/Den2)
			     )/(resolutions[i]*resolutions[i]*Den2);  
				
		A[2*4+0] += (-deltapx[i]*deltapy[i]*p1z[i])/(resolutions[i]*resolutions[i]*Den2);
		A[2*4+1] += (deltapx[i]*deltapx[i]*p1z[i])/(resolutions[i]*resolutions[i]*Den2);
		A[3*4+2] += (-deltapx[i]*deltapy[i]*p1z[i]*p1z[i] + driftdist[i]*(deltapx[i]*deltapy[i]*deltapy[i]*p1x[i] - deltapx[i]*deltapx[i]*deltapy[i]*p1y[i])/Den -
			     deltapx[i]*deltapy[i]*(deltapy[i]*deltapy[i]*p1x[i]*p1x[i] - deltapx[i]*deltapx[i]*p1y[i]*p1y[i] + 2*deltapx[i]*deltapy[i]*p1x[i]*p1y[i] )/Den2 
//			     +output_parameters[0]*(+driftdist[i]*deltapx[i]*deltapy[i]*deltapy[i]/Den - 2*deltapx[i]*deltapy[i]*(deltapy[i]*deltapy[i]*p1x[i]-deltapx[i]*p1y[i])/Den2) 
//			     +output_parameters[1]*(-driftdist[i]*deltapx[i]*deltapx[i]*deltapy[i]/Den - 2*deltapx[i]*deltapy[i]*(deltapx[i]*deltapx[i]*p1y[i]-deltapy[i]*p1x[i])/Den2) 
//			     +output_parameters[2]*(+driftdist[i]*deltapx[i]*deltapy[i]*deltapy[i]*p1z[i]/Den - 2*deltapx[i]*deltapy[i]*deltapy[i]*p1z[i]*(deltapx[i]*p1y[i]-deltapy[i]*p1x[i])/Den2)
//			     +output_parameters[3]*(-driftdist[i]*deltapx[i]*deltapx[i]*deltapy[i]*p1z[i]/Den + 2*deltapx[i]*deltapx[i]*deltapy[i]*p1z[i]*(deltapx[i]*p1y[i]-deltapy[i]*p1x[i])/Den2)
			     )/(resolutions[i]*resolutions[i]*Den2);
			     
		
		A[3*4+3] += (-deltapx[i]*deltapx[i]*p1z[i]*p1z[i])/(resolutions[i]*resolutions[i]*Den2);
	}
	
	matinv_4x4_matrix_per_thread(A, Ainv);

	printf("\n");
	for(int j = 0; j<4; j++){//row
                for(int k = 0; k<4; k++){//column
			printf("%1.6f ,", A[j*4+k]);
		}printf("\n");
	}printf("\n");

	printf("\n");
	for(int j = 0; j<4; j++){//row
                for(int k = 0; k<4; k++){//column
			printf("%1.6f ,", Ainv[j*4+k]);
		}printf("\n");
	}printf("\n");

	for(int j = 0; j<4; j++){//row
		output_parameters[j] = 0.0;
		output_parameters_errors[j] = sqrtf(fabs(A[j*4+j]));
		for(int k = 0; k<4; k++){//column
			output_parameters[j]+= Ainv[j*4+k]*B[k];
		}
	}
	
	//calculate chi2;
	chi2 = chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, output_parameters);
	
}



__device__ void linear_regression_3D(size_t const n_points, REAL *x_points, REAL *y_points, REAL *z_points, REAL *x_weights, REAL *y_weights, REAL *A, REAL* Ainv, REAL* B, REAL* output_parameters, REAL* output_parameters_errors, REAL& chi2)
{
	//For a 3D fit to a straight-line:
	// chi^2 = sum_i wxi * (xi- (X + Xp*zi))^2 + wyi*(y - (Y+Yp*zi))^2
	// dchi^2/dX = -2 * (xi - (X+Xp*zi))* wxi = 0
	// dchi^2/dY = -2 * (yi - (Y+Yp*zi))* wyi = 0
	// dchi^2/dXp = -2 * (xi - (X+Xp*zi))*zi * wxi = 0
	// dchi^2/dYp = -2 * (yi - (Y+Yp*zi))*zi * wyi = 0

	for(int j = 0; j<4; j++){
		B[j] = 0;
		for(int k = 0; k<4; k++){
			A[j*4+k] = 0;
			Ainv[j*4+k] = 0;
		}
	}
	
	//printf("\n");
	//for(int j = 0; j<4; j++){//row
        //        for(int k = 0; k<4; k++){//column
	//		printf("%1.6f ,", A[j*4+k]);
	//	}printf("\n");
	//}printf("\n");
	
	for( int i=0; i<n_points; i++ ){
		B[0] += x_weights[i]*x_points[i];
		B[1] += y_weights[i]*y_points[i];
		B[2] += x_weights[i]*x_points[i]*z_points[i];
		B[3] += y_weights[i]*y_points[i]*z_points[i];
		
		// first index: row; second index: col;
		// (consistent with what is used in the matrix inversion routine) 
		A[0*4+0] += x_weights[i];
		A[0*4+1] += 0.0;
		A[0*4+2] += x_weights[i]*z_points[i];
		A[0*4+3] += 0.0;

		A[1*4+0] += 0.0;
		A[1*4+1] += y_weights[i];
		A[1*4+2] += 0.0;
		A[1*4+3] += y_weights[i]*z_points[i];

		A[2*4+0] += x_weights[i]*z_points[i];
    		A[2*4+1] += 0.0;
		A[2*4+2] += x_weights[i]*z_points[i]*z_points[i];
		A[2*4+3] += 0.0;

    		A[3*4+0] += 0.0;
    		A[3*4+1] += y_weights[i]*z_points[i];
    		A[3*4+2] += 0.0;
    		A[3*4+3] += y_weights[i]*z_points[i]*z_points[i];

		printf("pt %d: x: %1.6f +- %1.6f, y: %1.6f +- %1.6f, z = %1.6f\n", i, x_points[i], x_weights[i], y_points[i], y_weights[i], z_points[i]);
    	}
	
	//printf("\n");
	//for(int j = 0; j<4; j++){//row
        //        for(int k = 0; k<4; k++){//column
	//		printf("%1.6f ,", A[j*4+k]);
	//	}printf("\n");
	//}printf("\n");
	
	matinv_4x4_matrix_per_thread(A, Ainv);
	
	//printf("\n");
	//for(int j = 0; j<4; j++){//row
        //        for(int k = 0; k<4; k++){//column
	//		printf("%1.6f ,", Ainv[j*4+k]);
	//	}printf("\n");
	//}printf("\n");
	
	//printf("\n");
	//REAL  I_test[16];
	//for(int j = 0; j<4; j++){//row
        //        for(int k = 0; k<4; k++){//column
	//		I_test[j*4+k] = 0.0;
	//		for(int l = 0; l<4; l++){
	//			I_test[j*4+k]+= A[j*4+l]*Ainv[l*4+k];
	//		}
	//		printf("%1.6f ,", I_test[j*4+k]);
	//	}printf("\n");
	//}printf("\n");
	
	for(int j = 0; j<4; j++){//row
		output_parameters[j] = 0.0;
		output_parameters_errors[j] = sqrtf(fabs(A[j*4+j]));
		for(int k = 0; k<4; k++){//column
			output_parameters[j]+= Ainv[j*4+k]*B[k];
		}
	}
	
	chi2 = 0;
	for( int i=0; i<n_points; i++ ){
		chi2+= x_weights[i]*(output_parameters[0]+z_points[i]*output_parameters[2])*(output_parameters[0]+z_points[i]*output_parameters[2])
		+y_weights[i]*(output_parameters[1]+z_points[i]*output_parameters[3])*(output_parameters[1]+z_points[i]*output_parameters[3]);
	}

}

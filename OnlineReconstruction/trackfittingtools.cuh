#include "operations.h"
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

        AA00 = A[0];
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



__device__ void track_residual_minimizer(size_t const n_points, 
	   				 REAL* driftdist, REAL* resolutions, 
					 REAL* p1x, REAL* p1y, REAL* p1z, 
					 REAL* deltapx, REAL* deltapy, REAL* deltapz, 
					 REAL *A, REAL* Ainv, REAL* B, 
					 REAL* output_parameters, REAL* output_parameters_errors, REAL& chi2)
{
	for(int j = 0; j<4; j++){
		B[j] = 0;
		for(int k = 0; k<4; k++){
			A[j*4+k] = 0;
			Ainv[j*4+k] = 0;
		}
	}
	
	for( int i=0; i<n_points; i++ ){
		B[0] += 0.0;
		B[1] += 0.0;
		B[2] += 0.0;
		B[3] += 0.0;
		
		// first index: row; second index: col;
		// (consistent with what is used in the matrix inversion routine) 
		A[0*4+0] += 0.0;
		A[0*4+1] += 0.0;
		A[0*4+2] += 0.0;
		A[0*4+3] += 0.0;

		A[1*4+0] += 0.0;
		A[1*4+1] += 0.0;
		A[1*4+2] += 0.0;
		A[1*4+3] += 0.0;

		A[2*4+0] += 0.0;
    		A[2*4+1] += 0.0;
		A[2*4+2] += 0.0;
		A[2*4+3] += 0.0;

    		A[3*4+0] += 0.0;
    		A[3*4+1] += 0.0;
    		A[3*4+2] += 0.0;
    		A[3*4+3] += 0.0;
	}
	
	matinv_4x4_matrix_per_thread(A, Ainv);
	
	for(int j = 0; j<4; j++){//row
		output_parameters[j] = 0.0;
		output_parameters_errors[j] = sqrtf(fabs(A[j*4+j]));
		for(int k = 0; k<4; k++){//column
			output_parameters[j]+= Ainv[j*4+k]*B[k];
		}
	}
	
	chi2 = 0;
	REAL dca;
	for( int i=0; i<n_points; i++ ){
		dca = ( -deltapy[i]*(p1x[i]-output_parameters[0]) + deltapx[i]*(p1y[i]-output_parameters[1]) + p1z[i]*(output_parameters[2]*deltapy[i]-output_parameters[3]*deltapx[i]) ) / ( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i] );
		chi2+= ( driftdist[i] - dca ) * ( driftdist[i] - dca ) / resolutions[i] / resolutions[i];
	}
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
	
	printf("\n");
	for(int j = 0; j<4; j++){//row
                for(int k = 0; k<4; k++){//column
			printf("%1.6f ,", A[j*4+k]);
		}printf("\n");
	}printf("\n");
	
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
	
	printf("\n");
	for(int j = 0; j<4; j++){//row
                for(int k = 0; k<4; k++){//column
			printf("%1.6f ,", A[j*4+k]);
		}printf("\n");
	}printf("\n");
	
	matinv_4x4_matrix_per_thread(A, Ainv);
	
	printf("\n");
	for(int j = 0; j<4; j++){//row
                for(int k = 0; k<4; k++){//column
			printf("%1.6f ,", Ainv[j*4+k]);
		}printf("\n");
	}printf("\n");
	
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

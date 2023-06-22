#include "trackfittingtools.cuh"


__device__ void chi2_simplefit(size_t const n_points, float* const x_array, float* const y_array, float& a, float& b, float sum, float det, float sx, float sy, float sxx, float syy, float sxy)
{
	sx = 0;
	sy = 0;
	sxx = 0;
	syy = 0;
	sxy = 0;
	
	for(int i = 0; i<n_points; i++){
		sum+=1.;
		sx+=x_array[i];
		sy+=y_array[i];
		sxx+=x_array[i]*x_array[i];
		syy+=y_array[i]*y_array[i];
		sxy+=x_array[i]*y_array[i];
	}
	det = sum*sxx - sx*sx; 	
	if(fabs(det)<1.e-20){
		a = 0.;
		b = 0.;
		return;
	}
	a = (sum*sxy - sx*sy)/det;
	b = (sy*sxx - sxy*sx)/det;

}

//Function to fit 2D tracks with errors
__device__ void fit_2D_track(size_t const n_points, const float *x_points, const float *z_points, const float *x_weights, float *A, float* Ainv, float* B, float* output_parameters, float* output_parameters_errors, float& chi2){
	//For a 2D fit to a straight-line:
	// chi^2 = sum_i wxi * (xi- (X + Xp*zi))^2
	// dchi^2/dX = -2 * (xi - (X+Xp*zi))* wxi = 0
	// dchi^2/dXp = -2 * (xi - (X+Xp*zi))*zi * wxi = 0
	
	for(int j = 0; j<2; j++){
		B[j] = 0;
		for(int k = 0; k<2; k++){
			A[j*2+k] = 0;
			Ainv[j*2+k] = 0;
		}
	}

	for( int i=0; i<n_points; i++ ){
	//if(blockIdx.x==11)printf("thread %d i %d x %1.4f z %1.4f \n", threadIdx.x, i, x_points[i], z_points[i]);
	//if(isnan(x_points[i]) || isnan(z_points[i]) || isnan(x_weights[i]))printf("%1.4f %1.4f %1.4f \n", x_points[i], z_points[i], x_weights[i]);
		B[0] += x_weights[i]*x_points[i];
		B[1] += x_weights[i]*x_points[i]*z_points[i];
		
		// first index: row; second index: col;
		// (consistent with what is used in the matrix inversion routine) 
		A[0] += x_weights[i];//0*2+0
		A[1] += x_weights[i]*z_points[i];//0*2+1
		
		A[2] += x_weights[i]*z_points[i];//1*2+0
		A[3] += x_weights[i]*z_points[i]*z_points[i];//1*2+1
    	}
	matinv_2x2_matrix_per_thread(A, Ainv);

	for(int j = 0; j<2; j++){//row
		output_parameters[j] = 0.0;
		output_parameters_errors[j] = sqrtf(fabs(Ainv[j*2+j]));
		for(int k = 0; k<2; k++){//column
			output_parameters[j]+= Ainv[j*2+k]*B[k];
		}
	}
	
	chi2 = 0;
	for( int i=0; i<n_points; i++ ){
		chi2+= x_weights[i]*x_weights[i]*(output_parameters[0]+z_points[i]*output_parameters[1])*(output_parameters[0]+z_points[i]*output_parameters[1]);
	}
}


//Same as previous, but the position is "corrected" for the drift
__device__ void fit_2D_track_drift(size_t const n_points, const float *x_points, const float *drift, const short *sign, const float *z_points, const float *x_weights, float *A, float* Ainv, float* B, float* output_parameters, float* output_parameters_errors, float& chi2){
	//For a 2D fit to a straight-line:
	// chi^2 = sum_i wxi * (xi- (X + Xp*zi))^2
	// dchi^2/dX = -2 * (xi - (X+Xp*zi))* wxi = 0
	// dchi^2/dXp = -2 * (xi - (X+Xp*zi))*zi * wxi = 0
	
	for(int j = 0; j<2; j++){
		B[j] = 0;
		for(int k = 0; k<2; k++){
			A[j*2+k] = 0;
			Ainv[j*2+k] = 0;
		}
	}

	for( int i=0; i<n_points; i++ ){
	//if(isnan(x_points[i]) || isnan(z_points[i]) || isnan(x_weights[i]))printf("%1.4f %1.4f %1.4f \n", x_points[i], z_points[i], x_weights[i]);
		B[0] += x_weights[i]*(x_points[i]+drift[i]*sign[i]);
		B[1] += x_weights[i]*(x_points[i]+drift[i]*sign[i])*z_points[i];
		
		// first index: row; second index: col;
		// (consistent with what is used in the matrix inversion routine) 
		A[0] += x_weights[i];//0*2+0
		A[1] += x_weights[i]*z_points[i];//0*2+1
		
		A[2] += x_weights[i]*z_points[i];//1*2+0
		A[3] += x_weights[i]*z_points[i]*z_points[i];//1*2+1
    	}
	matinv_2x2_matrix_per_thread(A, Ainv);

	for(int j = 0; j<2; j++){//row
		output_parameters[j] = 0.0;
		output_parameters_errors[j] = sqrtf(fabs(Ainv[j*2+j]));
		for(int k = 0; k<2; k++){//column
			output_parameters[j]+= Ainv[j*2+k]*B[k];
		}
	}
	
	chi2 = 0;
	for( int i=0; i<n_points; i++ ){
		chi2+= x_weights[i]*x_weights[i]*(output_parameters[0]+z_points[i]*output_parameters[1])*(output_parameters[0]+z_points[i]*output_parameters[1]);
	}
}






__device__ void fit_3D_track(size_t const n_points, float *x_points, float *y_points, float *z_points, float *x_weights, float *y_weights, float *A, float* Ainv, float* B, float* output_parameters, float* output_parameters_errors, float& chi2)
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
	
	for( int i=0; i<n_points; i++ ){
		B[0] += x_weights[i]*x_points[i];
		B[1] += y_weights[i]*y_points[i];
		B[2] += x_weights[i]*x_points[i]*z_points[i];
		B[3] += y_weights[i]*y_points[i]*z_points[i];
		
		// first index: row; second index: col;
		// (consistent with what is used in the matrix inversion routine) 
		A[0] += x_weights[i];//0*4+0
		A[1] += 0.0;//0*4+1
		A[2] += x_weights[i]*z_points[i];//0*4+2
		A[3] += 0.0;//0*4+3

		A[4] += 0.0;//1*4+0
		A[5] += y_weights[i];//1*4+1
		A[6] += 0.0;//1*4+2
		A[7] += y_weights[i]*z_points[i];//1*4+3

		A[8] += x_weights[i]*z_points[i];//2*4+0
    		A[9] += 0.0;//2*4+1
		A[10] += x_weights[i]*z_points[i]*z_points[i];//2*4+2
		A[11] += 0.0;//2*4+3

    		A[12] += 0.0;//3*4+0
    		A[13] += y_weights[i]*z_points[i];//3*4+1
    		A[14] += 0.0;//3*4+2
    		A[15] += y_weights[i]*z_points[i]*z_points[i];//3*4+3

		printf("pt %d: x: %1.6f +- %1.6f, y: %1.6f +- %1.6f, z = %1.6f\n", i, x_points[i], x_weights[i], y_points[i], y_weights[i], z_points[i]);
    	}
	
	for(int j = 0; j<4; j++){//row
		output_parameters[j] = 0.0;
		output_parameters_errors[j] = sqrtf(fabs(Ainv[j*4+j]));
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


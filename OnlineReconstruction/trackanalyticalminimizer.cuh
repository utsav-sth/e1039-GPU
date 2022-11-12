#include "trackfittingtools.cuh"

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
	
	REAL Den2, Den, dca;
	
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
	chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, 4, output_parameters, chi2, dca);
	
}

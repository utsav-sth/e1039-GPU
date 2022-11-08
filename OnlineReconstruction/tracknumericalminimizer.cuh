// numerical minimizer for the GPU tracks
// miminimzes chi^2 = sum (drift_dist - dca)^2/sigma^2
// 
// 1) from a starting point in the X = {x0, y0, tx, ty (, pinv)} and with starting step size (e.g 0.1 or 0.05 * the size of the acceptance in each variable) 
// evaluate (chi2(X+dX) - chi2(X))/dX and (chi2(X-dX) - chi2(X))/dX => go in the direction which gives the lowest dchi2; 
// 2) while (chi2(X+-dX) - chi2(X))/dX < 0 continue; 
// 3) if (chi2(X+dX) - chi2(X))/dX > 0, reduce the step in half and check both directions
// 4) if |(chi2(X+dX) - chi2(X))| < tolerance, stop;

// X+-dX is {x0+-dx0, y0, tx, ty (, pinv)}, {x0, y0+-dy0, tx, ty (, pinv)}, ...
// 

__device__ void chi2_straight(size_t const n_points, 
                              REAL* const driftdist, REAL* const resolutions,
                              REAL* const p1x, REAL* const p1y, REAL* const p1z,
                              REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
			      REAL* output_parameters, REAL& chi2)
{
	chi2 = 0;
	REAL dca;
	for( int i=0; i<n_points; i++ ){
	     	//this is the simplified expression of the chi2 where we neglect deltapz
		//TODO: plug in the real expression
		dca = ( -deltapy[i]*(p1x[i]-output_parameters[0]) + deltapx[i]*(p1y[i]-output_parameters[1]) + p1z[i]*(output_parameters[2]*deltapy[i]-output_parameters[3]*deltapx[i]) ) / ( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i] );
		chi2+= ( driftdist[i] - dca ) * ( driftdist[i] - dca ) / resolutions[i] / resolutions[i];
	}
}

__device__ void chi2_global(size_t const n_points, 
                            REAL* const driftdist, REAL* const resolutions,
                            REAL* const p1x, REAL* const p1y, REAL* const p1z,
                            REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
			    REAL* output_parameters, REAL& chi2)
{
	chi2 = 0;
	REAL dca;
	for( int i=0; i<n_points; i++ ){
	     	//this is the simplified expression of the chi2 where we neglect deltapz
		//TODO: plug in the real expression
		dca = ( -deltapy[i]*(p1x[i]-output_parameters[0]) + deltapx[i]*(p1y[i]-output_parameters[1]) + p1z[i]*(output_parameters[2]*deltapy[i]-output_parameters[3]*deltapx[i]) ) / ( deltapy[i]*deltapy[i] + deltapx[i]*deltapx[i] - 2*output_parameters[2]*output_parameters[3]*deltapy[i]*deltapx[i] );
		chi2+= ( driftdist[i] - dca ) * ( driftdist[i] - dca ) / resolutions[i] / resolutions[i];
	}
}

__device__ void chisquare(size_t const n_points, 
                            REAL* const driftdist, REAL* const resolutions,
                            REAL* const p1x, REAL* const p1y, REAL* const p1z,
                            REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
			    size_t const nparam, REAL* output_parameters, REAL& chi2)
{
	if(nparam==4)chi2_straight(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, output_parameters, chi2);
	if(nparam==5)chi2_global(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, output_parameters, chi2);
}


// numerical calculation of derivative
__device__ void calc_derivatives_num(size_t const n_points,
				     REAL* const driftdist, REAL* const resolutions,
				     REAL* const p1x, REAL* const p1y, REAL* const p1z,
				     REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
				     size_t const nparam, REAL* output_parameters, REAL* output_parameters_steps, REAL const lambda, 
				     REAL* derivatives, REAL* doublederivatives, REAL& chi2, REAL chi2prev)
{
	//calculation of chi2 at the nominal x, 
	chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2prev);
	
	for(int i = 0; i<nparam; i++)
	{
		//lambda scales the existing step
		output_parameters[i]+=output_parameters_steps[i]*lambda*0.01f;
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2);
		derivatives[i] = (chi2-chi2prev)/(output_parameters_steps[i]*lambda);   // (chi2(X+dX)-chi2(X))/dX
		
		output_parameters[i]-=2*output_parameters_steps[i]*lambda*0.01f;
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2);
		doublederivatives[i] = ((chi2prev-chi2)/(output_parameters_steps[i]*lambda)-derivatives[i])/output_parameters_steps[i]*lambda;
		output_parameters[i]+=output_parameters_steps[i]*lambda;
		
		output_parameters_steps[i] = -derivatives[i]/doublederivatives[i];
		//printf("%d, %1.3f, %1.3f, %1.3f\n", i, chi2prev, output_parameters[i], output_parameters_steps[i]);
	}

}



__device__ void calc_derivatives_ana(size_t const n_points,
				     REAL* const driftdist, REAL* const resolutions,
				     REAL* const p1x, REAL* const p1y, REAL* const p1z,
				     REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
				     size_t const nparam, REAL* output_parameters, REAL* output_parameters_steps, REAL const lambda,
				     REAL* derivatives, REAL& chi2, REAL chi2prev)
{
}


__device__ void minimize_chi2(size_t const max_iterations, REAL const tolerance,
	   		      REAL* const parameter_limits_min, REAL* const parameter_limits_max, 
			      size_t const n_points, 
                              REAL* const driftdist, REAL* const resolutions,
                              REAL* const p1x, REAL* const p1y, REAL* const p1z,
                              REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
			      size_t const nparam, REAL* output_parameters, REAL* prev_parameters, 
			      REAL* output_parameters_error, REAL* output_parameters_steps, REAL lambda,
			      REAL* derivatives, REAL* doublederivatives, REAL chi2, REAL chi2prev)
{
	//first, calculate the derivatives
	for(int n = 0; n<max_iterations; n++){
		calc_derivatives_num(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, output_parameters_steps, lambda, derivatives, doublederivatives, chi2, chi2prev);
		REAL deriv_prod = 1.;
		
		for(int i = 0; i<nparam; i++){
			prev_parameters[i] = output_parameters[i];
			output_parameters[i]+= output_parameters_steps[i];
			if(output_parameters[i]>parameter_limits_max[i])output_parameters[i]=parameter_limits_max[i];
			if(output_parameters[i]<parameter_limits_min[i])output_parameters[i]=parameter_limits_min[i];
			deriv_prod*= derivatives[i];
		}
		
		//calculate chi2 with new parameters:
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2);

		// "win condition": the delta chi2 is smaller than the allowed tolerance
		// TODO: add conditions to avoid accidental stops
		if(fabs(chi2prev-chi2)<tolerance || deriv_prod<tolerance){
			return;
		}
		
		// if the chi2 happens to be higher than the previous one, roll back the parameters:
		// and change the step scale factor 
		if(chi2>chi2prev){
			for(int i = 0; i<nparam; i++){
				output_parameters[i] = prev_parameters[i];
				lambda*0.5;
				//output_parameters_steps[i]*=0.5f;
			}
		}//do we want to reset lambda if we go in the right direction again?
	}
}
//#include "trackfittingtools.cuh"
#include "gpufit_core_funcs.cuh"

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

__device__ void minimize_chi2(size_t const max_iterations, REAL const tolerance,
	   		      REAL* const parameter_limits_min, REAL* const parameter_limits_max, 
			      size_t const n_points, 
                              REAL* const driftdist, REAL* const resolutions,
                              REAL* const p1x, REAL* const p1y, REAL* const p1z,
                              REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
			      size_t const nparam, REAL* output_parameters, REAL* prev_parameters, 
			      REAL* output_parameters_errors, REAL* output_parameters_steps, REAL lambda,
			      REAL* derivatives, //REAL* doublederivatives, 
			      REAL &chi2, REAL chi2prev, REAL dca)
{
	//first, calculate the derivatives
	for(int n = 0; n<max_iterations; n++){
		//printf("before deriv %d  %1.4f  %1.4f \n", n, chi2, chi2prev);
		calc_derivatives_num(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, output_parameters_steps, lambda, derivatives, 
		//doublederivatives, 
		chi2, chi2prev, dca);
		//deriv_prod = 1.;
		//printf("after deriv %d  %1.4f  %1.4f \n", n, chi2, chi2prev);
		for(int i = 0; i<nparam; i++){
		//printf("before: %d  %d  %1.4e  %1.4e  %1.4e  %1.4e \n", n, i, output_parameters[i], output_parameters_steps[i], derivatives[i], doublederivatives[i]);
			prev_parameters[i] = output_parameters[i];
			output_parameters[i]+= output_parameters_steps[i]*lambda;
			if(output_parameters[i]>parameter_limits_max[i])output_parameters[i]=parameter_limits_max[i];
			if(output_parameters[i]<parameter_limits_min[i])output_parameters[i]=parameter_limits_min[i];
			//deriv_prod*= derivatives[i];
		//printf("after: %d  %d  %1.4e  %1.4e  %1.4e  %1.4e \n", n, i, output_parameters[i], output_parameters_steps[i], derivatives[i], doublederivatives[i]);
		}
		
		//calculate chi2 with new parameters:
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2);
		//printf("before chi2 checks %d  %1.4f  %1.4f \n", n, chi2, chi2prev);

		// "win condition": the delta chi2 is smaller than the allowed tolerance
		// TODO: add conditions to avoid accidental stops
		//printf("%1.4f <? %1.4f \n", fabs(chi2prev-chi2), tolerance);
		if(fabs(chi2prev-chi2)<tolerance){
			chi2 = chi2>chi2prev ? chi2prev : chi2;
			return;
		}
				
		// if the chi2 happens to be higher than the previous one, roll back the parameters:
		// and change the step scale factor 
		if(chi2>chi2prev){
			chi2 = chi2prev;
			for(int i = 0; i<nparam; i++){
				output_parameters[i] = prev_parameters[i];
				//lambda*=0.5f;
				// 1: are we sure we want to do that here? other wise we get s
				// 
				//output_parameters_steps[i]*=0.5f;
			}
		}else{//do we want to reset lambda if we go in the right direction again?
			//update chi2prev
			chi2prev = chi2;
			
		}

		//printf("after chi2 checks %d  %1.4f  %1.4f \n", n, chi2, chi2prev);

	}
}



__device__ void brute_force_search_area(REAL* const parameter_limits_min, REAL* const parameter_limits_max, 
					size_t const n_points, 
					REAL* const driftdist, REAL* const resolutions,
					REAL* const p1x, REAL* const p1y, REAL* const p1z,
					REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
					REAL* output_parameters, REAL* output_parameters_error, 
					REAL* const fixed_point, REAL &chi2, REAL chi2min, REAL dca)
{
	chi2min = 10000000.0f;
	
	int i_min, j_min;
	for(int i = 0; i<29; i++){
		output_parameters[0] = parameter_limits_min[0]+5.+10.*i;
		output_parameters[2] = (fixed_point[0]-output_parameters[0])/fixed_point[2];
		
		if(output_parameters[2]<parameter_limits_min[2]||output_parameters[2]>parameter_limits_max[2])continue;
				
		for(int j = 0; j<9; j++){
			output_parameters[1] = parameter_limits_min[1]+5.+10.*i;
			output_parameters[3] = (fixed_point[1]-output_parameters[1])/fixed_point[2];
			
			if(output_parameters[3]<parameter_limits_min[3]||output_parameters[3]>parameter_limits_max[3])continue;
			
			chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, 4, output_parameters, chi2);
			if(chi2<chi2min){
				chi2min = chi2;
				i_min = i;
				j_min = j;
			}
		}
	}
	//printf("%1.6f \n", chi2min);
	
	output_parameters[0] = parameter_limits_min[0]+5.+10.*i_min;
	output_parameters[1] = parameter_limits_min[1]+5.+10.*j_min;
	output_parameters[2] = (fixed_point[0]-output_parameters[0])/fixed_point[2];
	output_parameters[3] = (fixed_point[1]-output_parameters[1])/fixed_point[2];
}





//"fixed-point" minimizer: we vary x0, tx proportionally to each other, as well as y0, ty, so that the track points to the fixed point
__device__ void fixedpoint_straighttrack_chi2minimizer(size_t const max_iterations, REAL const tolerance,
	   						REAL* const parameter_limits_min, REAL* const parameter_limits_max, 
							size_t const n_points, 
							REAL* const driftdist, REAL* const resolutions,
							REAL* const p1x, REAL* const p1y, REAL* const p1z,
							REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
							REAL* A, REAL* Ainv, REAL* B, 
							size_t const nparam, REAL* output_parameters, REAL* output_parameters_errors, 
							REAL* prev_parameters, 
							REAL* output_parameters_steps, REAL lambda, REAL* derivatives, 
							REAL* const fixed_point, REAL &chi2, REAL chi2prev, REAL dca)
{
	dca = 0;
	//first, get the fixed point
	get_straighttrack_fixedpoint(n_points,
				     driftdist, resolutions,
				     p1x, p1y, p1z,
				     deltapx, deltapy, deltapz,
				     A, Ainv, B,
				     output_parameters, output_parameters_errors, 
				     fixed_point, chi2);

	// correlate the steps with each other:
	output_parameters_steps[0] = fixed_point[0]-output_parameters_steps[2]*fixed_point[2];
	output_parameters_steps[1] = fixed_point[0]-output_parameters_steps[2]*fixed_point[3];

	/*
	// for tests: brute force:
	brute_force_search_area(parameter_limits_min, parameter_limits_max, 
				n_points, 
				driftdist, resolutions,
				p1x, p1y, p1z,
				deltapx, deltapy, deltapz,
				output_parameters, output_parameters_errors, 
				fixed_point, chi2, chi2prev, dca);
	
	for(int k = 0; k<nparam; k++){
		output_parameters_steps[k] = (parameter_limits_max[k]-parameter_limits_min[k])*0.025f;
	}
	minimize_chi2(max_iterations, tolerance,
		      parameter_limits_min, parameter_limits_max, 
		      n_points, 
                      driftdist, resolutions,
                      p1x, p1y, p1z,
                      deltapx, deltapy, deltapz,
		      nparam, output_parameters, prev_parameters, 
		      output_parameters_errors, output_parameters_steps, lambda,
		      derivatives, //doublederivatives, 
		      chi2, chi2prev, dca);
	*/
	
	// then get in the loop:
	for(int i = 0; i<max_iterations; i++){
		//first get the derivatives
		calc_corr_derivatives_num(n_points, 
					  driftdist, resolutions, 
			    	  	  p1x, p1y, p1z, 
			     	  	  deltapx, deltapy, deltapz, 
			     	  	  nparam, output_parameters, 
			     		  output_parameters_steps, lambda, 
			     		  derivatives, 
			     		  chi2, chi2prev, dca);

		for(int i = 0; i<nparam; i++){
			prev_parameters[i] = output_parameters[i];
			output_parameters[i]+= output_parameters_steps[i];
			if(output_parameters[i]>parameter_limits_max[i])output_parameters[i]=parameter_limits_max[i];
			if(output_parameters[i]<parameter_limits_min[i])output_parameters[i]=parameter_limits_min[i];
		}
		
		//calculate chi2 with new parameters:
		chisquare(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, nparam, output_parameters, chi2);
		
		// "win condition": the delta chi2 is smaller than the allowed tolerance
		if(fabs(chi2prev-chi2)<tolerance){
			//how do we make sure we're not stuck in a local minimum?
			chi2 = chi2>chi2prev ? chi2prev : chi2;
			return;
		}
	}
	
}


//device?
//global? // as is, I might not have the choice, since I need to "reconfigure" the thread/block for the gauss_jordan_method
__device__ void gpufit_algorithm_fitter(int const max_n_iterations, 
					// status indices
					int& n_iterations, short& iteration_failed, short& finished, 
					short& state, short& skip, short& singular, 
					// fit configuration parameters
					REAL const tolerance,
					int const n_parameters,
					REAL* const lower_bounds, REAL* const upper_bounds,
					// variables to optimize
					REAL* parameters, REAL* prev_parameters, 
					REAL& chi_square, REAL& prev_chi_square, 
					REAL* values, REAL* derivatives, REAL* gradients, 
					REAL* hessians, REAL* scaling_vector, 
					REAL* deltas, REAL& lambda,
					REAL* Ainv,
					//REAL* calc_matrix, REAL* abs_row, int* abs_row_index,
					// exp data arrays 
					int const n_points,
					REAL* const driftdist, REAL* const resolutions,
					REAL* const p1x, REAL* const p1y, REAL* const p1z,
					REAL* const deltapx, REAL* const deltapy, REAL* const deltapz)
{
	// call: 0, 10, ...
	project_parameter_to_box(n_parameters,
				 parameters,
				 lower_bounds, 
				 upper_bounds);

	// call: 1, 11, ...
 	_calc_curve_values(parameters, 
			   n_points,
			   n_parameters,
			   finished,
			   values, derivatives,
			   driftdist, resolutions,
			   p1x, p1y, p1z,
			   deltapx, deltapy, deltapz,
			   chi_square);
	// call: 3, 13, ...
	_check_fit_improvement(iteration_failed,
			       chi_square,
   			       prev_chi_square,
			       finished);
	//		       printf("%d ?\n", iteration_failed);
	// call: 4, 14, ...
	_calculate_gradients(gradients,
			     values,
			     derivatives,
			     driftdist, resolutions,
			     n_points,
			     n_parameters,
			     finished,
			     skip);
	//printf("%d %1.4f %1.4f %1.4f %1.4f \n", finished, derivatives[0], gradients[0], parameters[0], chi_square);
	//call: 5, 15, ...
	_calculate_hessians(hessians,
			    values,
			    derivatives,
			    n_parameters,
			    n_points,
			    resolutions,
			    skip,
			    finished);
	
	for(int n = 0; n<max_n_iterations; n++){
	//printf("%d %d, x0 = %1.4f, chi2 %1.4f, prev_chi2 = %1.4f \n", n, finished, parameters[0], chi_square, prev_chi_square);
	//printf("%d %d, derivative[0,0] %1.8f, gradients[0] %1.4f, x0 = %1.4f \n", n, finished, derivatives[0*6], gradients[0], parameters[0]);
	//printf("%d %d, derivative[1,0] %1.8f, gradients[1] %1.4f, y0 = %1.4f \n", n, finished, derivatives[1*6], gradients[1], parameters[1]);
	//printf("%d %d, derivative[2,0] %1.8f, gradients[2] %1.4f, tx = %1.4f \n", n, finished, derivatives[2*6], gradients[2], parameters[2]);
	//printf("%d %d, derivative[3,0] %1.8f, gradients[3] %1.4f, ty = %1.4f \n", n, finished, derivatives[3*6], gradients[3], parameters[3]);
		//call: 6, 19,...
		_modify_step_widths(hessians,
				    lambda,
				    scaling_vector,
				    n_parameters,
				    iteration_failed,
				    finished);

		//need a lighter matrix inversion!
		solve_equations_systems(n_parameters,
					hessians, gradients,
					Ainv, deltas,
					skip);		
		/*
		//call: 7, 20, ... 
		// configure first
		dim3  threads(1, 1, 1);
		dim3  blocks(1, 1, 1);

		threads.x = n_parameters + 1;
		threads.y = n_parameters;
		blocks.x = 1;
		int nparam_pow2 = 1;
		for(int i = 0; i<n_parameters; i++){
			nparam_pow2*=2;
		}
		int const shared_size = sizeof(REAL) * ((threads.x * threads.y) + nparam_pow2 + nparam_pow2);
		
		_gaussjordan<<<blocks, threads, shared_size>>>(deltas,
				gradients,
				hessians,
				finished,
				singular,
				//calc_matrix, 
				//abs_row,
				//abs_row_index,
				n_parameters,
				nparam_pow2);
		*/
		//call: 8, 21, ...
	//printf("%d %d, deltas[0] = %1.4f, deltas[1] = %1.4f, deltas[2] =  %1.4f, deltas[3] = %1.4f\n", n, finished, deltas[0], deltas[1], deltas[2], deltas[3]);

		_update_state_after_solving(singular,
					    finished,
					    state);
		//call: 9, 22, ...
		_update_parameters(
			parameters,
			prev_parameters,
			deltas,
			n_parameters,
			finished);
		// call: 0, 10, ...
		project_parameter_to_box(n_parameters,
				parameters,
				lower_bounds, 
				upper_bounds);

		// call: 1, 11, ...
		_calc_curve_values(parameters, 
			   n_points,
			   n_parameters,
			   finished,
			   values, derivatives,
			   driftdist, resolutions,
			   p1x, p1y, p1z,
			   deltapx, deltapy, deltapz,
			   chi_square);
		// call: 3, 13, ...
		_check_fit_improvement(iteration_failed,
			       chi_square,
   			       prev_chi_square,
			       finished);
		//       printf("%d ?\n", iteration_failed);
		// call: 4, 14, ...
		_calculate_gradients(gradients,
			     values,
			     derivatives,
			     driftdist, resolutions,
			     n_points,
			     n_parameters,
			     finished,
			     skip);
		//call: 5, 15, ...
		_calculate_hessians(hessians,
			    values,
			    derivatives,
			    n_parameters,
			    n_points,
			    resolutions,
			    skip,
			    finished);
		//call: 16, ...
		_check_for_convergence(finished,
			tolerance,
    			state,
    			chi_square,
    			prev_chi_square,
    			n,
    			max_n_iterations);
		//call: 17
		_evaluate_iteration(n_iterations,
				finished,
    				n,
    				state);
		// call: 18
		_prepare_next_iteration(lambda,
			chi_square,
			prev_chi_square,
			parameters,
			prev_parameters,
			n_parameters);
	//printf("%d %d, x0 = %1.4f, prev_chi2 = %1.4f, chi2 = %1.4f \n", n, finished, parameters[0], prev_chi_square, chi_square);

	}
	//printf("%d, %1.6f \n", n_iterations, chi_square);

}
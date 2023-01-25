#ifndef DEVICE_FIT_INTERFACE_CUH_INCLUDED
#define DEVICE_FIT_INTERFACE_CUH_INCLUDED

#include "constants.h"
#include "cuda_kernels.cu"
#include "models/models.cuh"
#include "estimators/estimators.cuh"


__device__ void _gpu_configure_model(int const model_id, int& n_parameters_, int& n_dimensions_)
{
    switch (model_id)
    {
    case GAUSS_1D:              n_parameters_ = 4; n_dimensions_ = 1; break;
    case GAUSS_2D:              n_parameters_ = 5; n_dimensions_ = 2; break;
    case GAUSS_2D_ELLIPTIC:     n_parameters_ = 6; n_dimensions_ = 2; break;
    case GAUSS_2D_ROTATED:      n_parameters_ = 7; n_dimensions_ = 2; break;
    case CAUCHY_2D_ELLIPTIC:    n_parameters_ = 6; n_dimensions_ = 2; break;
    case LINEAR_1D:             n_parameters_ = 2; n_dimensions_ = 1; break;
    case LINEAR_3D:             n_parameters_ = 4; n_dimensions_ = 2; break;
    case FLETCHER_POWELL_HELIX: n_parameters_ = 3; n_dimensions_ = 1; break;
    case BROWN_DENNIS:          n_parameters_ = 4; n_dimensions_ = 1; break;
    case SPLINE_1D:             n_parameters_ = 3; n_dimensions_ = 1; break;
    case SPLINE_2D:             n_parameters_ = 4; n_dimensions_ = 2; break;
    case SPLINE_3D:             n_parameters_ = 5; n_dimensions_ = 3; break;
    case SPLINE_3D_MULTICHANNEL:         n_parameters_ = 5; n_dimensions_ = 4; break;
    case SPLINE_3D_PHASE_MULTICHANNEL:   n_parameters_ = 6; n_dimensions_ = 4; break;
    default: printf("unknown model ID"); n_parameters_ = 0; n_dimensions_ = 0; break;
    }
}

__device__ int _gpu_run_fit(
	std::size_t n_fits,
	int max_n_iterations,
	REAL tolerance,
	int model_id,
	int n_parameters, 
	int n_dimensions, 
	int estimator_id,
	std::size_t n_points,
	REAL* data,
	REAL* weights,
	std::size_t user_info_size,
	char* user_info,
	REAL* initial_parameters,
	int* parameters_to_fit,
	REAL* output_parameters,
	int* output_states,
	REAL* output_chi_squares,
	int* output_n_iterations,
	REAL* constraints,
	int* constraint_types
)
{
	/*
	for(int i = 0; i<n_parameters; i++){
		output_parameters[i] = initial_parameters[i];
	}
	REAL* values;// same as number of points
	REAL* derivatives;
	REAL* chi_squares;
	REAL* gradients;
	int* states;
	
	//for(int k = 0; k<n+parameters; k++)derivatives[k] = 0;
		
	if(constraints!=NULL){
		for(int i = 0; i<n_parameters; i++){
			REAL lower_bound = constraints[i*2+LOWER_BOUND];
			REAL upper_bound = constraints[i*2+UPPER_BOUND];
			project_parameter_to_box(output_parameters[i], lower_bound, upper_bound, constraint_types[i]);
		}
	}
	
	//_calc_curve_values();
	calculate_model(
		(ModelID)model_id, 
		output_parameters, 
		n_fits, 
		n_points, 
		values, 
		derivatives, 
		0, 
		0, 
		0,
		user_info,
		user_info_size
	);
	
	calculate_chi_square(
		estimator_id,
		chi_squares,
		0,
		data,
		values,
		weights,
		states,
		user_info,
		user_info_size
	);
	
	for (int parameter_index = 0; parameter_index < n_parameters; parameter_index++){
		int derivative_index = parameters_to_fit_indices[parameter_index]*n_points;
		calculate_gradient(
			estimator_id,
                    	gradients,
                    	0,
                    	derivative_index,
                    	data,
                    	values,
                    	derivatives,
                    	weights,
                    	user_info,
			user_info_size
		);
	}
	*/	
	

	return 1;
}


__device__ int device_fit_interface(
	std::size_t n_fits,
	int max_n_iterations,
	REAL tolerance,
	int model_id,
	int estimator_id,
	std::size_t n_points,
	REAL* data,
	REAL* weights,
	std::size_t user_info_size,
	char* user_info,
	REAL* initial_parameters, 
	int* parameters_to_fit,
	REAL* fit_parameters,
	int* output_states,
	REAL* output_chi_squares,
	int* output_n_iterations,
	REAL* constraints,
	int* constraint_types
)
{
	int n_parameters = 0, n_dimensions = 0;
	_gpu_configure_model(model_id, n_parameters, n_dimensions);
	if(n_parameters==0 || n_dimensions==0)return (0);
	
	int status = _gpu_run_fit(
		n_fits,
		max_n_iterations,
		tolerance,
		model_id,
		n_parameters,
		n_dimensions,
		estimator_id,
		n_points,
		data,
		weights,
		user_info_size,
		user_info, 
		initial_parameters, 
		parameters_to_fit,
        	fit_parameters,
		output_states,
        	output_chi_squares,
        	output_n_iterations,
		constraints,
		constraint_types
	);
		
	return status;
}





#endif

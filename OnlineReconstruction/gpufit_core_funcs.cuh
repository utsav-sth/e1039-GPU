//#include "definitions.h"
#include "trackfittingtools.cuh"

//calculate calculates residual values, derivatives and chi2 at the same time
// call: 1, 11, ...
__device__ void _calc_curve_values(
    REAL* const parameters,
    int const n_points,
    int const n_parameters,
    short const finished,
    REAL* values,
    REAL* derivatives,
    //data input array
    REAL* const driftdist, REAL* const resolutions,
    REAL* const p1x, REAL* const p1y, REAL* const p1z,
    REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
    REAL& chi2)
{
    if (finished)
        return;
	
    calc_val_derivatives(n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, n_parameters, parameters, values, derivatives, chi2);
        
}

__device__ void _calc_curve_corr_values(
    REAL* const parameters,
    int const n_points,
    int const n_parameters,
    short const finished,
    REAL* values,
    REAL* derivatives,
    //data input array
    REAL* fixedpoint,
    REAL* const driftdist, REAL* const resolutions,
    REAL* const p1x, REAL* const p1y, REAL* const p1z,
    REAL* const deltapx, REAL* const deltapy, REAL* const deltapz,
    REAL& chi2)
{
    if (finished)
        return;
	
    calc_corr_val_derivatives(fixedpoint, n_points, driftdist, resolutions, p1x, p1y, p1z, deltapx, deltapy, deltapz, n_parameters, parameters, values, derivatives, chi2);
        
}

// call: 2, 12, ... chi2 (merged with 1)


// call: 3, 13, ...
__device__ void _check_fit_improvement(
    short& iteration_failed,
    REAL const chi_square,
    REAL const prev_chi_square,
    short const finished)
{
    bool const prev_chi_square_initialized = prev_chi_square != 0.;
    bool const chi_square_decreased = (chi_square < prev_chi_square);
    if (prev_chi_square_initialized && !chi_square_decreased)
    {
        iteration_failed = 1;
    }
    else
    {
        iteration_failed = 0;
    }
}


// call: 4, 14, ...
__device__ void _calculate_gradients(
    REAL* gradients,
    REAL const * values,
    REAL const * derivatives,
    REAL* const driftdist, REAL* const resolutions,
    int const n_points,
    int const n_parameters,
    short const finished,
    short const skip)
{

    if (finished || skip)
    {
        return;
    }
    for(int j = 0; j< n_parameters; j++){
    	    gradients[j] = 0;
    	    for(int i = 0; i<n_points; i++){
		gradients[j]+= derivatives[j*n_points]*(0.0-values[i])/resolutions[i]/resolutions[i];
	}
    }
//    calc_gradients(n_points, 
//	           resolutions, driftdist,
//		   n_parameters,
//		   values, derivatives, gradients);
}


//call: 5, 15, ...
__device__ void _calculate_hessians(
    REAL * hessians,
    REAL const * values,
    REAL const * derivatives,
    int const n_parameters,
    int const n_points,
    REAL* const resolutions,
    short const skip,
    short const finished)
{
    if (finished || skip)
    {
        return;
    }

    for(int i = 0; i<n_points; i++){
    	//j: row; k: col; 
	for(int j = 0; j< n_parameters; j++){
		for(int k = 0; k<n_parameters; k++){
               		hessians[j*n_parameters+k]+= derivatives[j*n_points]*derivatives[k*n_points]/resolutions[i]/resolutions[i];
		}
	}
    }
}


//call: 6
__device__ void _modify_step_widths(
    REAL* hessians,
    REAL const lambda,
    REAL * scaling_vector,
    unsigned int const n_parameters,
    short const iteration_failed,
    short const finished)
{
    for(int parameter_index = 0; parameter_index<n_parameters; parameter_index++){

        if (finished)
    	{
	return;
    	}

        int const diagonal_index = parameter_index * n_parameters + parameter_index;

    	if (iteration_failed)
    	{
    	    hessians[diagonal_index] -= scaling_vector[parameter_index] * lambda / 10.;
    	}

	// adaptive scaling
    	scaling_vector[parameter_index]
		= max(scaling_vector[parameter_index], hessians[diagonal_index]);

    	// continuous scaling
    	//scaling_vector[parameter_index] = hessian[diagonal_index];
    	
    	// initial scaling
    	//if (scaling_vector[parameter_index] == 0.)
    	//    scaling_vector[parameter_index] = hessian[diagonal_index];

    	hessians[diagonal_index] += scaling_vector[parameter_index] * lambda;
    }
}

//call: 7 
__device__ void solve_equations_systems(std::size_t const n_param,
	  				REAL* A, REAL* const B, //A: hessians; B: gradients; 
					REAL* Ainv, REAL* output,     //output: deltas;
					short const skip_calculation)
{
    if (skip_calculation)
        return;

    for(int j = 0; j<n_param; j++){
            for(int k = 0; k<n_param; k++){
                    Ainv[j*n_param+k] = 0;
            }
    }

    matinv_4x4_matrix_per_thread(A, Ainv);
    
    for(int j = 0; j<n_param; j++){
    	    output[j] = 0;
            for(int k = 0; k<n_param; k++){
                    output[j]+= Ainv[j*n_param+k]*B[k];
            }
    }
}

//call: 7 solve_equation_systems_gj()
__global__ void _gaussjordan(
    REAL * delta,
    REAL const * beta,
    REAL const * alpha,
    short const skip_calculation,
    short& singular,
    //REAL* calculation_matrix,
    //REAL* abs_row,
    //int* abs_row_index,
    std::size_t const n_equations,
    std::size_t const n_equations_pow2)
{
    extern __shared__ REAL extern_array[];     //shared memory between threads of a single block, 
    //used for storing the calculation_matrix, the 
    //abs_row vector, and the abs_row_index vector

    // In this routine we will store the augmented matrix (A|B), referred to here
    // as the calculation matrix in a shared memory space which is visible to all
    // threads within a block.  Also stored in shared memory are two vectors which 
    // are used to find the largest element in each row (the pivot).  These vectors 
    // are called abs_row and abs_row_index.
    //
    // Sizes of data stored in shared memory:
    //
    //      calculation_matrix: n_equations * (n_equations+1)
    //      abs_row:            n_equations_pow2
    //      abs_row_index:      n_equations_pow2
    //  
    // Note that each thread represents an element of the augmented matrix, with
    // the column and row indicated by the x and y index of the thread.  Each 
    // solution is calculated within one block, and the solution index is the 
    // block index x value.

    int const col_index = threadIdx.x;                  //column index in the calculation_matrix
    int const row_index = threadIdx.y;                  //row index in the calculation_matrix
    int const solution_index = blockIdx.x;

    int const n_col = blockDim.x;                       //number of columns in calculation matrix (=threads.x)
    int const n_row = blockDim.y;                       //number of rows in calculation matrix (=threads.y)
    int const alpha_size = blockDim.y * blockDim.y;     //number of entries in alpha matrix for one solution (NxN)

    if (skip_calculation)
        return;

    REAL p;                                            //local variable used in pivot calculation

    REAL * calculation_matrix = extern_array;                          //point to the shared memory

    REAL * abs_row = extern_array + n_equations * (n_equations + 1);     //abs_row is located after the calculation_matrix
    //within the shared memory

    int * abs_row_index = (int *)(abs_row + n_equations_pow2);            //abs_row_index is located after abs_row
    //
    //note that although the shared memory is defined as
    //REAL, we are storing data of type int in this
    //part of the shared memory

    //initialize the singular vector
    if (col_index == 0 && row_index == 0)
    {
        singular = 0;
    }

    //initialize abs_row and abs_row_index, using only the threads on the diagonal
    if (col_index == row_index)
    {
        abs_row[col_index + (n_equations_pow2 - n_equations)] = 0.0;
        abs_row_index[col_index + (n_equations_pow2 - n_equations)] = col_index + (n_equations_pow2 - n_equations);
    }

    //initialize the calculation_matrix (alpha and beta, concatenated, for one solution)
    if (col_index != n_equations)
        calculation_matrix[row_index*n_col + col_index] = alpha[solution_index * alpha_size + row_index * n_equations + col_index];
    else
        calculation_matrix[row_index*n_col + col_index] = beta[solution_index * n_equations + row_index];

    //wait for thread synchronization

    __syncthreads();

    //start of main outer loop over the rows of the calculation matrix

    for (int current_row = 0; current_row < n_equations; current_row++)
    {

        // work in only one row, skipping the last column
        if (row_index == current_row && col_index != n_equations)
        {

            //save the absolute values of the current row
            abs_row[col_index] = abs(calculation_matrix[row_index * n_col + col_index]);

            //save the column indices
            abs_row_index[col_index] = col_index;

            __threadfence();

            //find the largest absolute value in the current row and write its index in abs_row_index[0]
            for (int n = 2; n <= n_equations_pow2; n = n * 2)
            {
                if (col_index < (n_equations_pow2 / n))
                {
                    if (abs_row[abs_row_index[col_index]] < abs_row[abs_row_index[col_index + (n_equations_pow2 / n)]])
                    {
                        abs_row_index[col_index] = abs_row_index[col_index + (n_equations_pow2 / n)];
                    }
                }
            }
        }

        __syncthreads();

        //singularity check - if all values in the row are zero, no solution exists
        if (row_index == current_row && col_index != n_equations)
        {
            if (abs_row[abs_row_index[0]] == 0.0)
            {
                singular = 1;
            }
        }

        //devide the row by the biggest value in the row
        if (row_index == current_row)
        {
            calculation_matrix[row_index * n_col + col_index]
                = calculation_matrix[row_index * n_col + col_index] / calculation_matrix[row_index * n_col + abs_row_index[0]];
        }

        __syncthreads();

        //The value of the largest element of the current row was found, and then current
        //row was divided by this value such that the largest value of the current row 
        //is equal to one.  
        //
        //Next, the matrix is manipulated to reduce to zero all other entries in the column 
        //in which the largest value was found.   To do this, the values in the current row
        //are scaled appropriately and substracted from the other rows of the matrix. 
        //
        //For each element of the matrix that is not in the current row, calculate the value
        //to be subtracted and let each thread store this value in the scalar variable p.

        p = calculation_matrix[current_row * n_col + col_index] * calculation_matrix[row_index * n_col + abs_row_index[0]];
        __syncthreads();

        if (row_index != current_row)
        {
            calculation_matrix[row_index * n_col + col_index] = calculation_matrix[row_index * n_col + col_index] - p;
        }
        __syncthreads();

    }

    //At this point, if the solution exists, the calculation matrix has been reduced to the 
    //identity matrix on the left side, and the solution vector on the right side.  However
    //we have not swapped rows during the procedure, so the identity matrix is out of order.
    //
    //For example, starting with the following augmented matrix as input:
    //
    //  [  3  2 -4 |  4 ]
    //  [  2  3  3 | 15 ]
    //  [  5 -3  1 | 14 ]
    //
    //we will obtain:
    //
    //  [  0  0  1 |  2 ]
    //  [  0  1  0 |  1 ]
    //  [  1  0  0 |  3 ]
    //
    //Which needs to be re-arranged to obtain the correct solution vector.  In the final
    //step, each thread checks to see if its value equals 1, and if so it assigns the value
    //in its rightmost column to the appropriate entry in the beta vector.  The solution is
    //stored in beta upon completetion.

    if (col_index != n_equations && calculation_matrix[row_index * n_col + col_index] == 1)
        delta[n_row * solution_index + col_index] = calculation_matrix[row_index * n_col + n_equations];

    __syncthreads();
}



//call: 8
__device__ void _update_state_after_solving(
    short const cublas_info,
    short const finished,
    short& state)
{
    if (finished)
        return;

    if (cublas_info!= 0)
        state = SINGULAR_HESSIAN;
}


//call: 9
__device__ void _update_parameters(
    REAL * parameters,
    REAL * prev_parameters,
    REAL const * deltas,
    int const n_parameters,
    short const finished)
{
    for(int parameter_index = 0; parameter_index<n_parameters; parameter_index++){
        prev_parameters[parameter_index] = parameters[parameter_index];

    	if (finished)
	    {
		return;
	    }

    	parameters[parameter_index] += deltas[parameter_index];
    }
    
}

__device__ void _update_corr_parameters(
    REAL * parameters,
    REAL * prev_parameters,
    REAL* fixedpoint,
    REAL const * deltas,
    int const n_parameters,
    int const n_total_parameters,
    short const finished)
{
    if (finished)
    {
	return;
    }

    for(int parameter_index = 0; parameter_index<n_total_parameters; parameter_index++){
        prev_parameters[parameter_index] = parameters[parameter_index];
	
    	//if()parameters[parameter_index] += deltas[parameter_index];
    }
    parameters[2]+= deltas[0];
    parameters[3]+= deltas[1];
    parameters[0]+= -deltas[0]*fixedpoint[2];
    parameters[1]+= -deltas[1]*fixedpoint[2];
    
}

//call: 0, 10, ...
__device__ void project_parameter_to_box(
	   	int const n_parameters,
	   	REAL* parameters, 
		REAL* const lower_bounds, 
		REAL* const upper_bounds)
{
	for(int i = 0; i< n_parameters; i++){
		if(parameters[i]<lower_bounds[i])parameters[i]=lower_bounds[i];
		if(parameters[i]>upper_bounds[i])parameters[i]=upper_bounds[i];
	}
}

__device__ void project_parameter_to_box_fixedpoint(
	   	int const n_parameters,
	   	REAL* parameters, 
		REAL* fixedpoint,
		REAL* const lower_bounds, 
		REAL* const upper_bounds)
{
	for(int i = 0; i< n_parameters; i++){
		if(i<2){
			if(parameters[i]<lower_bounds[i]){
				parameters[2+i]=parameters[2+i]-(parameters[i]-lower_bounds[i])/fixedpoint[2];
				parameters[i]=lower_bounds[i];
			}
			if(parameters[i]>upper_bounds[i]){
				parameters[2+i]=parameters[2+i]-(parameters[i]-upper_bounds[i])/fixedpoint[2];
				parameters[i]=upper_bounds[i];
			}
		}else{
			if(parameters[i]<lower_bounds[i]){
				parameters[i-2]=parameters[i-2]-(parameters[i]-lower_bounds[i])*fixedpoint[2];
				parameters[i]=lower_bounds[i];
			}
			if(parameters[i]>upper_bounds[i]){
				parameters[i-2]=parameters[i-2]-(parameters[i]-upper_bounds[i])*fixedpoint[2];
				parameters[i]=upper_bounds[i];
			}
		}
	}
}


//call: 16
__device__ void _check_for_convergence(
    short& finished,
    REAL const tolerance,
    short& state,
    REAL const chi_square,
    REAL const prev_chi_square,
    int const iteration,
    int const max_n_iterations)
{
    if (finished)
    {
        return;
    }

    int const fit_found
        = abs(chi_square - prev_chi_square)
        < tolerance * max(1., chi_square);

    int const max_n_iterations_reached = iteration == max_n_iterations - 1;

    if (fit_found)
    {
        finished = 1;
    }
    else if (max_n_iterations_reached)
    {
        state = MAX_ITERATION;
    }
}

//call: 17
__device__ void _evaluate_iteration(
    int& n_iterations,
    short& finished,
    int const iteration,
    short const state)
{
    if (state != CONVERGED)
    {
        finished = 1;
    }

    if (finished && n_iterations == 0)
    {
        n_iterations = iteration + 1;
    }
}

// call: 18
__device__ void _prepare_next_iteration(
    REAL& lambda,
    REAL& chi_square,
    REAL& prev_chi_square,
    REAL* parameters,
    REAL const* prev_parameters,
    int const n_parameters)
{
    if (chi_square < prev_chi_square)
    {
        lambda *= 0.1f;
        prev_chi_square = chi_square;
    }
    else
    {
        lambda *= 10.;
        chi_square = prev_chi_square;
        for (int iparameter = 0; iparameter < n_parameters; iparameter++)
        {
            parameters[iparameter] = prev_parameters[iparameter];
        }
    }
}

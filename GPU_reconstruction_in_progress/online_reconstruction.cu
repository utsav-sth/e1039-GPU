#include "gpufit.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <string>
#include <ctime>
#include <chrono>

// CUDA runtime
// #include <cuda_runtime.h>
#include <cublas_v2.h>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>
#include <thrust/transform.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <thrust/sort.h>

// #include <cublasLt.h>


#include <TObject.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TRandom.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include "LoadInput.h"

#define nChamberPlanes 30
#define nHodoPlanes 16
#define nPropPlanes 8
#define nDetectors (nChamberPlanes+nHodoPlanes+nPropPlanes)
#define Epsilon 0.00001f

#define triggerBit(n) (1 << (n))
#define hitFlagBit(n) (1 << (n))

using namespace std;

const int EstnEvtMax = 10240;
const int THREADS_PER_BLOCK = 512;
int BLOCKS_NUM = EstnEvtMax/THREADS_PER_BLOCK;
const int EstnAHMax = 5000;
const int EstnTHMax = 200;
const int ClusterSizeMax = 100;

// function to check GPU status
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

// looping on all GPUs to check their status
void printDeviceStatus() {
	int nDevices;

	
	gpuErrchk( cudaGetDeviceCount(&nDevices) );
	for (int i = 0; i < nDevices; i++) {
	  cudaDeviceProp prop;
	  cudaGetDeviceProperties(&prop, i);
	  printf("Device Number: %d\n", i);
	  printf("  Device name: %s\n", prop.name);
	  printf("  Memory Clock Rate (KHz): %d\n",
			 prop.memoryClockRate);
	  printf("  Memory Bus Width (bits): %d\n",
			 prop.memoryBusWidth);
	  printf("  Peak Memory Bandwidth (GB/s): %f\n\n",
			 2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
	}
}





//clone of LoadEvent::Hit:
class gHit {
	public:
	int index; // global hit index in the hit array
	short detectorID; // ID of the detector: one ID for each DC wire plane (30 total), hodoscope plane (16 total), proportional tube plane (8 total).
	short elementID; // ID of the element in the detector: wire/slat/tube number
	float tdcTime; // raw TDC time from the DAQ 
	float driftDistance; // calculated drift distance from RT profile (supplied in database) IF tdcTime between tmin and tmax defined for detector; 
	float pos; // position in the projection of the detector (e.g. X in a X plane, etc)
	short flag; // 1: in time; 2: hodo mask; 3: trigger mask
};

class gEvent {
	public:
	int RunID; // Run Number
	int EventID; // Event number
	int SpillID; // Spill number
	int TriggerBits; // hash of the trigger bits: 0-4: MATRIX1-5; 5-9: NIM1-5;
	short TargetPos; // target position: proxy for target ID?
	int TurnID; // => related to beam intensity
	int RFID; // => related to beam intensity
	int Intensity[33]; //  16 before, one onset, and 16 after
	short TriggerEmu; // 1 if MC event
	short NRoads[4]; // 0, positive top; 1, positive bottom; 2, negative top; 3, negative bottom
	int NHits[nDetectors+1]; // number of hits in each detector plane
	int nAH; // size of AllHits
	int nTH; // size of TriggerHits
	gHit AllHits[EstnAHMax]; // array of all hits
	gHit TriggerHits[EstnTHMax]; // array of trigger hits
};

// SW = SearchWindows?
class gSW {
public:
	int EventID;
	int nAH;
};

// Hit comparison
struct lessthan {
	__host__ __device__ bool operator()(const gHit& lhs, const gHit& rhs)
	{
	//returns true if :
	// hit1.detID<hit2.detID;  
		if(lhs.detectorID < rhs.detectorID)
		{
			return true;
		}
		else if(lhs.detectorID > rhs.detectorID)
		{
			return false;
		}
	//hit1.detID=hit2.detID & hit1.elID<hit2.elID;
		if(lhs.elementID < rhs.elementID)
		{
			return true;
		}
		else if(lhs.elementID > rhs.elementID)
		{
			return false;
		}
	//hit1.detID=hit2.detID & hit1.elID=hit2.elID & hit1.time>hit2.time;
		if(lhs.tdcTime > rhs.tdcTime)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
};


// position of first hit 
struct get_first_event_hit_pos
{
  __host__ __device__
  float operator()(gEvent& evt)
  {
    return evt.AllHits[0].pos;
  }
};

// Linear regression: fit of tracks
void linear_regression_example(int n_points_per_fit, REAL *device_input, thrust::device_vector<REAL> &d_parameters)
{
	// number of fits, fit points and parameters
	size_t const n_fits = 1;
	size_t const n_model_parameters = 2;


	// custom z positions, stored in user info
	// NOTE: this is the way to initialize a device_vector with many constant values
	vector< REAL > user_info_values {-130.43f, -131.437f, -116.782f, -115.712f, -128.114f, -129.111f };
	thrust::device_vector<REAL> d_user_info(user_info_values.size());
	thrust::copy(user_info_values.begin(), user_info_values.end(), d_user_info.begin());

	// size of user info in bytes
	size_t const user_info_size = d_user_info.size() * sizeof(REAL); 

	// // initial parameters (randomized)
	// std::vector< REAL > initial_parameters(n_fits * n_model_parameters);
	// for (size_t i = 0; i != n_fits; i++)
	// {
	// 	// random offset
	// 	initial_parameters[i * n_model_parameters + 0] = true_parameters[0] * (0.8f + 0.4f * uniform_dist(rng));
	// 	// random slope
	// 	initial_parameters[i * n_model_parameters + 1] = true_parameters[0] * (0.8f + 0.4f * uniform_dist(rng));
	// }

	// // generate data
	// std::vector< REAL > data(n_points_per_fit * n_fits);
	// for (size_t i = 0; i != data.size(); i++)
	// {

	// 	size_t j = i / n_points_per_fit; // the fit
	// 	size_t k = i % n_points_per_fit; // the position within a fit

	// 	REAL x = user_info[k];
	// 	REAL y = true_parameters[0] + x * true_parameters[1];
	// 	data[i] = y;



	// }

	// tolerance
	REAL const tolerance = 0.001f;

	// maximum number of iterations
	int const max_number_iterations = 20;

	// estimator ID
	int const estimator_id = LSE;

	// model ID
	int const model_id = LINEAR_1D;

	// parameters to fit (all of them)
	std::vector< int > parameters_to_fit(n_model_parameters, 1);

	
	thrust::device_vector< int > d_states(n_fits);
	thrust::device_vector< REAL > d_chi_square(n_fits);
	thrust::device_vector< int > d_number_iterations(n_fits);

	//call to gpufit (C interface)
	// can be found in https://github.com/gpufit/Gpufit/blob/master/Gpufit/gpufit.cpp
	// parameters:
	// size_t n_fits,
    	// size_t n_points,
	// float * gpu_data,
	// float * gpu_weights,
    	// int model_id,
    	// float tolerance,
    	// int max_n_iterations,
    	// int * parameters_to_fit,
    	// int estimator_id,
	// size_t user_info_size,
    	// char * gpu_user_info,
    	// float * gpu_fit_parameters,
    	// int * gpu_output_states,
    	// float * gpu_output_chi_squares,
    	// int * gpu_output_n_iterations
	// size_t n_fits,

	int const status = gpufit_cuda_interface
       (
        	n_fits,
        	n_points_per_fit,
            device_input,
            0,
            model_id,
            // initial_parameters.data(),
            tolerance,
            max_number_iterations,
            parameters_to_fit.data(),
			// true_parameters.data(),
            estimator_id,
            user_info_size,
            reinterpret_cast< char * >( thrust::raw_pointer_cast(d_user_info.data()) ),
            thrust::raw_pointer_cast(d_parameters.data()),
            thrust::raw_pointer_cast(d_states.data()),
            thrust::raw_pointer_cast(d_chi_square.data()),
            thrust::raw_pointer_cast(d_number_iterations.data())
        );

		

	// check status
	if (status != ReturnState::OK)
	{
		throw std::runtime_error(gpufit_get_last_error());
	}

	// get fit states
	std::vector< int > output_states_histogram(5, 0);
	for (auto it = d_states.begin(); it != d_states.end(); ++it)
	{
		output_states_histogram[*it]++;
	}

	std::cout << "ratio converged              " << (REAL) output_states_histogram[0] / n_fits << "\n";
	std::cout << "ratio max iteration exceeded " << (REAL) output_states_histogram[1] / n_fits << "\n";
	std::cout << "ratio singular hessian       " << (REAL) output_states_histogram[2] / n_fits << "\n";
	std::cout << "ratio neg curvature MLE      " << (REAL) output_states_histogram[3] / n_fits << "\n";
	std::cout << "ratio gpu not read           " << (REAL) output_states_histogram[4] / n_fits << "\n";

	// compute mean fitted parameters for converged fits
	std::vector< REAL > output_parameters_mean(n_model_parameters, 0);
	for (size_t i = 0; i != n_fits; i++)
	{
		if (d_states[i] == FitState::CONVERGED)
		{
			// add offset
			output_parameters_mean[0] += d_parameters[i * n_model_parameters + 0];
			// add slope
			output_parameters_mean[1] += d_parameters[i * n_model_parameters + 1];
		}
	}
	output_parameters_mean[0] /= output_states_histogram[0];
	output_parameters_mean[1] /= output_states_histogram[0];

	// compute std of fitted parameters for converged fits
	std::vector< REAL > output_parameters_std(n_model_parameters, 0);
	for (size_t i = 0; i != n_fits; i++)
	{
		if (d_states[i] == FitState::CONVERGED)
		{
			// add squared deviation for offset
			output_parameters_std[0] += (d_parameters[i * n_model_parameters + 0] - output_parameters_mean[0]) * (d_parameters[i * n_model_parameters + 0] - output_parameters_mean[0]);
			// add squared deviation for slope
			output_parameters_std[1] += (d_parameters[i * n_model_parameters + 1] - output_parameters_mean[1]) * (d_parameters[i * n_model_parameters + 1] - output_parameters_mean[1]);
		}
	}
	// divide and take square root
	output_parameters_std[0] = sqrt(output_parameters_std[0] / output_states_histogram[0]);
	output_parameters_std[1] = sqrt(output_parameters_std[1] / output_states_histogram[0]);

	// print mean and std
	std::cout << "offset  true " << d_parameters[0] << " mean " << output_parameters_mean[0] << " std " << output_parameters_std[0] << "\n";
	std::cout << "slope   true " << d_parameters[1] << " mean " << output_parameters_mean[1] << " std " << output_parameters_std[1] << "\n";

	// compute mean chi-square for those converged
	REAL  output_chi_square_mean = 0;
	for (size_t i = 0; i != n_fits; i++)
	{
		if (d_states[i] == FitState::CONVERGED)
		{
			output_chi_square_mean += d_chi_square[i];
		}
	}
	output_chi_square_mean /= static_cast<REAL>(output_states_histogram[0]);
	std::cout << "mean chi square " << output_chi_square_mean << "\n";

	// compute mean number of iterations for those converged
	REAL  output_number_iterations_mean = 0;
	for (size_t i = 0; i != n_fits; i++)
	{
		if (d_states[i] == FitState::CONVERGED)
		{
			output_number_iterations_mean += static_cast<REAL>(d_number_iterations[i]);
		}
	}

	// normalize
	output_number_iterations_mean /= static_cast<REAL>(output_states_histogram[0]);
	std::cout << "mean number of iterations " << output_number_iterations_mean << "\n";
}


// kernel functions: 
// CUDA C++ extends C++ by allowing the programmer to define C++ functions, called kernels, that, when called, 
// are executed N times in parallel by N different CUDA threads, as opposed to only once like regular C++ functions. 
// I guess TKL is for tracklet selection? 
__global__ void gkernel_TKL(gEvent* ic, gSW* sw) {
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	sw[index].nAH = ic[index].nAH;
	sw[index].EventID = ic[index].EventID;
}

// event reducer: 
__global__ void gkernel_eR(gEvent* ic) {
	//printf("Running the kernel function...\n");
	// retrieve global thread index
	int index = threadIdx.x + blockIdx.x * blockDim.x;

	double w_max[EstnEvtMax]; // max drift distance of the hit furthest from the cluster avg position // current average position of cluster * 0.9
	
	double w_min[EstnEvtMax]; // max drift distance of the hit closest to the cluster avg position // current average position of cluster * 0.4
	double dt_mean[EstnEvtMax]; // tbd
	int cluster_iAH_arr_cur[EstnEvtMax]; // current cluster array
	int cluster_iAH_arr_size[EstnEvtMax]; // cluster size i.e number of hits in cluster
	static int cluster_iAH_arr[EstnEvtMax][ClusterSizeMax]; // global cluster array 
	int uniqueID[EstnEvtMax]; // hit unique element ID
	int uniqueID_curr[EstnEvtMax]; // current hit unique element ID
	double tdcTime_curr[EstnEvtMax]; // current hit TDC time
	int iAH[EstnEvtMax]; // hit index 
	int nAH_reduced[EstnEvtMax]; // number of hits after hit quality filtering
	// int nHitsPerDetector[nDetectors+1];
	

	// initialization of array size
	cluster_iAH_arr_size[index] = 0;
	nAH_reduced[index] = 0;
		
	// event reducing/hit filtering
	for(iAH[index] = 0; iAH[index]<ic[index].nAH; ++iAH[index]) {
		// if hit not good, set its detID to 0 and continue;
		if((ic[index].AllHits[iAH[index]].flag & hitFlagBit(1)) == 0) {
//			printf("Skip out-of-time...\n");
			ic[index].AllHits[iAH[index]].detectorID = 0;
			continue;
		}
		// hits in DCs or Prop tubes
		if(ic[index].AllHits[iAH[index]].detectorID < 31 || ic[index].AllHits[iAH[index]].detectorID > 46) {
			// evaluate "unique ID"
			uniqueID[index] = ic[index].AllHits[iAH[index]].detectorID*1000 + ic[index].AllHits[iAH[index]].elementID;
			// compare with current unique element ID; if different, update the unique element ID and time info 
			if(uniqueID[index] != uniqueID_curr[index]) {
				uniqueID_curr[index] = uniqueID[index];
				tdcTime_curr[index] = ic[index].AllHits[iAH[index]].tdcTime;
			}
			// if next hit and current hit belong to the same element: if detID>36 => prop tubes (reminder that hodoscpes are out of the picture in this scope), 
			// we're suppose to have one signal (a second signal would be after-pulsing)
			// if time difference between new hit and current hit is less than 80ns for DCs, it's also considered after-pulsing
			else {
				if(ic[index].AllHits[iAH[index]].detectorID > 36 || ((ic[index].AllHits[iAH[index]].tdcTime - tdcTime_curr[index] >= 0.0) && (ic[index].AllHits[iAH[index]].tdcTime - tdcTime_curr[index] < 80.0)) || ((ic[index].AllHits[iAH[index]].tdcTime - tdcTime_curr[index] <= 0.0) && (ic[index].AllHits[iAH[index]].tdcTime - tdcTime_curr[index] > -80.0))) {
//					printf("Skip after-pulse...\n");
					ic[index].AllHits[iAH[index]].detectorID = 0;
					continue;
				}
				else {
					tdcTime_curr[index] = ic[index].AllHits[iAH[index]].tdcTime;
				}
			}
		}
		// declustering of hits in DCs (from CPU code, I understand this one better)
		// if there are hits in the same plane and hitting two neighboring wires, they both give redundant information: 
		if(ic[index].AllHits[iAH[index]].detectorID <= nChamberPlanes) {
//			printf("%d\n", cluster_iAH_arr_size[index]);
//			printf("Decluster...\n");
			if(cluster_iAH_arr_size[index] == ClusterSizeMax) {
//				printf("Oversized cluster...\n");
			}
			// if array size is zero, start storing the hit in the array
			if(cluster_iAH_arr_size[index] == 0) {
				cluster_iAH_arr[index][0] = iAH[index];
				++cluster_iAH_arr_size[index];
			}
			// otherwise
			else {
				// current hit and previous hit are *not* in same detector plane OR next hit and current hit are *not* in neighbors cells
				// we "declusterize" i.e. we remove the hit/hits which information is redundant with other hits and/or useless
				if((ic[index].AllHits[iAH[index]].detectorID != ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID) || (ic[index].AllHits[iAH[index]].elementID - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID > 1)) {
					// if 2 hits in cluster, evaluate w_max and w_min; drift distance has to be < w_min for one of the hits, while it has to be < w_max for the other hit 
					if(cluster_iAH_arr_size[index] == 2) {
						w_max[index] = 0.9*0.5*(ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].pos - ic[index].AllHits[cluster_iAH_arr[index][0]].pos);
						w_min[index] = 4.0/9.0*w_max[index];
						if((ic[index].AllHits[cluster_iAH_arr[index][0]].driftDistance > w_max[index] && ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].driftDistance > w_min[index]) || (ic[index].AllHits[cluster_iAH_arr[index][0]].driftDistance > w_min[index] && ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].driftDistance > w_max[index])) {
							//eliminating the existing hit with the lagest drift distance
							if(ic[index].AllHits[cluster_iAH_arr[index][0]].driftDistance > ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].driftDistance) {
//								printf("Skip cluster...\n");
								ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID = 0;
							}
							else {
//								printf("Skip cluster...\n");
								ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID = 0;
							}
						}
						// if the time difference is less than 8 ns for detectors 19 to 24 (which btw are DC3p) we remove both
						else if((((ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) >= 0.0 && (ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) < 8.0) || ((ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) <= 0.0 && (ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) > -8.0)) && (ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID >= 19 && ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID <= 24)) {
//							printf("Skip cluster...\n");
							ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID = 0;
							ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID = 0;
						}
					}
					// if 3 hits or more in cluster: we essentially discard them all;
					if(cluster_iAH_arr_size[index] >= 3) {
						// evaluate the mean time difference;
						dt_mean[index] = 0.0;
						for(cluster_iAH_arr_cur[index] = 1; cluster_iAH_arr_cur[index] < cluster_iAH_arr_size[index]; ++cluster_iAH_arr_cur[index]) {
							dt_mean[index] += ((ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]-1]].tdcTime) > 0.0 ? (ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]-1]].tdcTime) : (ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]-1]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].tdcTime));
						}
						dt_mean[index] = dt_mean[index]/(cluster_iAH_arr_size[index] - 1);
						// if mean time difference is less than 10, that's electronic noise, so we remove them all.
						if(dt_mean[index] < 10.0) {
//							printf("Skip cluster...\n");
							for(cluster_iAH_arr_cur[index] = 0; cluster_iAH_arr_cur[index] < cluster_iAH_arr_size[index]; ++cluster_iAH_arr_cur[index]) {
								ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].detectorID = 0;
							}
						}
						// otherwise, we remove them all except first and last
						else {
//							printf("Skip cluster...\n");
							for(cluster_iAH_arr_cur[index] = 1; cluster_iAH_arr_cur[index] < cluster_iAH_arr_size[index]; ++cluster_iAH_arr_cur[index]) {
								ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].detectorID = 0;
							}
						}
					}
					cluster_iAH_arr_size[index] = 0;
				}
				// current hit and previous hit are in same detector plane and in neighbor wires: 
				// we count how many hits we have in this case, until we find a hit in a different detector or in a wire that is not a neighbor to the previous hit.
				cluster_iAH_arr[index][cluster_iAH_arr_size[index]] = iAH[index];
				++cluster_iAH_arr_size[index];
			}
		}
	}
	//end of the hit loop

	// Hit reduction: 
	// store in "AllHits" containers only hits with non-zero detectorID and couting those with nAH_reduced
	for(iAH[index] = 0; iAH[index]<ic[index].nAH; ++iAH[index]) {
		if(ic[index].AllHits[iAH[index]].detectorID != 0) {
			ic[index].AllHits[nAH_reduced[index]] = ic[index].AllHits[iAH[index]];
			++nAH_reduced[index];

						
		}
	}

	// compute hits per detector
	int nEventHits = nAH_reduced[index];
	// reinitialize number of hits per detector
	for(auto iDetector = 1; iDetector <= nDetectors; ++iDetector) {
		ic[index].NHits[iDetector] = 0;
	}
	// loop on reduced hits and counting number of hits per detector
	for(auto iHit = 0; iHit < nEventHits; ++iHit) {
		auto detectorId = ic[index].AllHits[iHit].detectorID;
		if(detectorId != 0) {
			++ic[index].NHits[detectorId];
		}
	}


	ic[index].nAH = nAH_reduced[index];
	//ic[index].NHits[0] = NHits_reduced[index];
	


	//if(((ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6])>0) || ((ic[index].NHits[7]+ic[index].NHits[8]+ic[index].NHits[9]+ic[index].NHits[10]+ic[index].NHits[11]+ic[index].NHits[12])>0) || ((ic[index].NHits[13]+ic[index].NHits[14]+ic[index].NHits[15]+ic[index].NHits[16]+ic[index].NHits[17]+ic[index].NHits[18])>0) || ((ic[index].NHits[19]+ic[index].NHits[20]+ic[index].NHits[21]+ic[index].NHits[22]+ic[index].NHits[23]+ic[index].NHits[24])>0) || ((ic[index].NHits[25]+ic[index].NHits[26]+ic[index].NHits[27]+ic[index].NHits[28]+ic[index].NHits[29]+ic[index].NHits[30])>0)){	
	//if(((ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6])<270) || ((ic[index].NHits[7]+ic[index].NHits[8]+ic[index].NHits[9]+ic[index].NHits[10]+ic[index].NHits[11]+ic[index].NHits[12])>350) || ((ic[index].NHits[13]+ic[index].NHits[14]+ic[index].NHits[15]+ic[index].NHits[16]+ic[index].NHits[17]+ic[index].NHits[18])>170) || ((ic[index].NHits[19]+ic[index].NHits[20]+ic[index].NHits[21]+ic[index].NHits[22]+ic[index].NHits[23]+ic[index].NHits[24])>140) || ((ic[index].NHits[25]+ic[index].NHits[26]+ic[index].NHits[27]+ic[index].NHits[28]+ic[index].NHits[29]+ic[index].NHits[30])>140))

	//we do not accept the event unless there is at least one hit in the first DC

if( (ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6])<1)


	{
	
		//printf("Event rejected...\n");
	}
	else {
		//counting total hit number, for all events < 6668? why?
		if( (ic[index].EventID)<6668){ 	
			int totalDetectorHits = 0;
			for(int i = 1; i <= nDetectors; ++i) {
				totalDetectorHits += ic[index].NHits[i];
			}

			int nFirstRegionHits = 0;
			for(int i = 1; i < 6; ++i) {
				nFirstRegionHits += ic[index].NHits[i];
				printf("nHits[%d] = %d\n", i, ic[index].NHits[i]);
			}
		
		// printf("AllHits value : %d\n", (ic[index].NHits[0]));
		printf("reduced AllHits value : %d\n", (nAH_reduced[index]));
		printf("sum of detectors : %d (%d)\n", totalDetectorHits, nFirstRegionHits);
		//}
}
	    
//		Process the accepted events (tracking) here.
		// where is the tracking though?
	}
}






// test code

#include <curand.h>

// // Fill the array A(nr_rows_A, nr_cols_A) with random numbers on GPU
// void GPU_fill_rand(float *A, int nr_rows_A, int nr_cols_A) {
// 	// Create a pseudo-random number generator
// 	curandGenerator_t prng;
// 	curandCreateGenerator(&prng, CURAND_RNG_PSEUDO_DEFAULT);

// 	// Set the seed for the random number generator using the system clock
// 	curandSetPseudoRandomGeneratorSeed(prng, (unsigned long long) clock());

// 	// Fill the array with random numbers on the device
// 	curandGenerateUniform(prng, A, nr_rows_A * nr_cols_A);
// }

// Multiply the arrays A and B on GPU and save the result in C
// C(m,n) = A(m,k) * B(k,n)
void gpu_blas_mmul(const float *A, const float *B, float *C, const int m, const int k, const int n) {
	int lda=m,ldb=k,ldc=m;
	const float alf = 1;
	const float bet = 0;
	const float *alpha = &alf;
	const float *beta = &bet;

	// Create a handle for CUBLAS
	cublasHandle_t handle;
	cublasCreate(&handle);

	// Do the actual multiplication
	cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

	// Destroy the handle
	cublasDestroy(handle);
}

//Print matrix A(nr_rows_A, nr_cols_A) storage in column-major format
void print_matrix(const thrust::device_vector<float> &A, int nr_rows_A, int nr_cols_A) {

    for(int i = 0; i < nr_rows_A; ++i){
        for(int j = 0; j < nr_cols_A; ++j){
            std::cout << A[j * nr_rows_A + i] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// /**
//  * Online sample
//  * @see https://github.com/sol-prog/cuda_cublas_curand_thrust/blob/master/mmul_2.cu
//  */
// int main(int argc, char* argv[]) {
// 	// Allocate 3 arrays on CPU
// int nr_rows_A, nr_cols_A, nr_rows_B, nr_cols_B, nr_rows_C, nr_cols_C;
 
// // for simplicity we are going to use square arrays
// nr_rows_A = nr_cols_A = nr_rows_B = nr_cols_B = nr_rows_C = nr_cols_C = 3;
 
// thrust::device_vector<float> d_A(nr_rows_A * nr_cols_A), d_B(nr_rows_B * nr_cols_B), d_C(nr_rows_C * nr_cols_C);
 
// // Fill the arrays A and B on GPU with random numbers
// GPU_fill_rand(thrust::raw_pointer_cast(&d_A[0]), nr_rows_A, nr_cols_A);
// GPU_fill_rand(thrust::raw_pointer_cast(&d_B[0]), nr_rows_B, nr_cols_B);
 
// // Optionally we can print the data
// std::cout << "A =" << std::endl;
// print_matrix(d_A, nr_rows_A, nr_cols_A);
// std::cout << "B =" << std::endl;
// print_matrix(d_B, nr_rows_B, nr_cols_B);
 
// // Multiply A and B on GPU
// gpu_blas_mmul(thrust::raw_pointer_cast(&d_A[0]), thrust::raw_pointer_cast(&d_B[0]), thrust::raw_pointer_cast(&d_C[0]), nr_rows_A, nr_cols_A, nr_cols_B);
 
// //Print the result
// std::cout << "C =" << std::endl;
// print_matrix(d_C, nr_rows_C, nr_cols_C);
// }




int main(int argn, char * argv[]) {
	
	// initialization: declaration of SRaw event, opening file/tree, affecting rawEvent object to input tree
	// declaring array of gEvent;
	auto start = std::chrono::system_clock::now();
	clock_t cp1 = clock();

	TString inputFile;
	TString outputFile;
	inputFile = argv[1];
	outputFile = argv[2];

	cout<<"Running "<<argv[0]<<endl;
	cout<<"Loading "<<argv[1]<<endl;
	cout<<"Writing "<<argv[2]<<endl;

	SRawEvent* rawEvent = new SRawEvent();
	TFile* dataFile = new TFile(inputFile.Data(), "READ");
	TTree* dataTree = (TTree *)dataFile->Get("save");
	dataTree->SetBranchAddress("rawEvent", &rawEvent);
	int nEvtMax = dataTree->GetEntries();
	static gEvent host_gEvent[EstnEvtMax];


	// loop on event: get RawEvent information and load it into gEvent
	for(int i = 0; i < nEvtMax; ++i) {
		dataTree->GetEntry(i);
//		cout<<"Converting "<<i<<"/"<<nEvtMax<<endl;

		host_gEvent[i].RunID = rawEvent->fRunID;
		host_gEvent[i].EventID = rawEvent->fEventID;
		host_gEvent[i].SpillID = rawEvent->fSpillID;
		host_gEvent[i].TriggerBits = rawEvent->fTriggerBits;
		host_gEvent[i].TargetPos = rawEvent->fTargetPos;
		host_gEvent[i].TurnID = rawEvent->fTurnID;
		host_gEvent[i].RFID = rawEvent->fRFID;
		for(int j=0; j<33; j++) {
			host_gEvent[i].Intensity[j] = rawEvent->fIntensity[j];
		}
		host_gEvent[i].TriggerEmu = rawEvent->fTriggerEmu;
		for(int k=0; k<4; k++) {
			host_gEvent[i].NRoads[k] = rawEvent->fNRoads[k];
		}
		for(int l=0; l<(nChamberPlanes+nHodoPlanes+nPropPlanes+1); l++) {
			host_gEvent[i].NHits[l] = rawEvent->fNHits[l];
		}
		host_gEvent[i].nAH = rawEvent->fAllHits.size();
		host_gEvent[i].nTH = rawEvent->fTriggerHits.size();
		for(int m=0; m<rawEvent->fAllHits.size(); m++) {
			host_gEvent[i].AllHits[m].index=(rawEvent->fAllHits[m]).index;
			host_gEvent[i].AllHits[m].detectorID=(rawEvent->fAllHits[m]).detectorID;
			host_gEvent[i].AllHits[m].elementID=(rawEvent->fAllHits[m]).elementID;
			host_gEvent[i].AllHits[m].tdcTime=(rawEvent->fAllHits[m]).tdcTime;
			host_gEvent[i].AllHits[m].driftDistance=(rawEvent->fAllHits[m]).driftDistance;
			host_gEvent[i].AllHits[m].pos=(rawEvent->fAllHits[m]).pos;
			host_gEvent[i].AllHits[m].flag=(rawEvent->fAllHits[m]).flag;
		}
		for(int n=0; n<rawEvent->fTriggerHits.size(); n++) {
			host_gEvent[i].TriggerHits[n].index=(rawEvent->fTriggerHits[n]).index;
			host_gEvent[i].TriggerHits[n].detectorID=(rawEvent->fTriggerHits[n]).detectorID;
			host_gEvent[i].TriggerHits[n].elementID=(rawEvent->fTriggerHits[n]).elementID;
			host_gEvent[i].TriggerHits[n].tdcTime=(rawEvent->fTriggerHits[n]).tdcTime;
			host_gEvent[i].TriggerHits[n].driftDistance=(rawEvent->fTriggerHits[n]).driftDistance;
			host_gEvent[i].TriggerHits[n].pos=(rawEvent->fTriggerHits[n]).pos;
			host_gEvent[i].TriggerHits[n].flag=(rawEvent->fTriggerHits[n]).flag;
		}
	}



//If the decoded has NOT been sorted...
//	for(int i = 0; i < nEvtMax; ++i) {
//		thrust::stable_sort(host_gEvent[i].AllHits, host_gEvent[i].AllHits+host_gEvent[i].nAH, lessthan());
//	}


	// evaluate the total size of the gEvent array (and the SW array) for memory allocation 
	// (the memory cannot be dynamically allocated) 
	size_t NBytesAllEvent = EstnEvtMax * sizeof(gEvent);
	size_t NBytesAllSearchWindow = EstnEvtMax * sizeof(gSW);

	gEvent *host_output_eR = (gEvent*)malloc(NBytesAllEvent);
	gSW * host_output_TKL = (gSW*)malloc(NBytesAllSearchWindow);

	// declaring gEvent objects for the device (GPU) to use.
	gEvent *device_gEvent;
	// gEvent *device_output_eR;
	gEvent *device_input_TKL;
	gSW *device_output_TKL;
	
	
	// copy of data from host to device: evaluate operation time 
	clock_t cp2 = clock();
	auto start_kernel = std::chrono::system_clock::now();

	// printDeviceStatus();
	// Allocating memory for GPU (pointer to allocated device ); check for errors in the process; stops the program if issues encountered
	gpuErrchk( cudaMalloc((void**)&device_gEvent, NBytesAllEvent));
	gpuErrchk( cudaMalloc((void**)&device_input_TKL, NBytesAllEvent));
	gpuErrchk( cudaMalloc((void**)&device_output_TKL, NBytesAllSearchWindow));

	// cudaMemcpy(dst, src, count, kind): copies data between host and device:
	// dst: destination memory address; src: source memory address; count: size in bytes; kind: type of transfer
	// cudaMalloc((void**)&device_output_eR, sizeofoutput_eR);
	gpuErrchk( cudaMemcpy(device_gEvent, host_gEvent, NBytesAllEvent, cudaMemcpyHostToDevice));
	// cudaMemcpy(device_output_eR, host_output, sizeofoutput_eR, cudaMemcpyHostToDevice);
	auto end_kernel = std::chrono::system_clock::now();
	
	// now data is transfered in the device: kernel function for event reconstruction called;
	// note that the function call is made requesting a number of blocks and a number of threads per block
	// in practice we have as many threads total as number of events; 
	//auto start_kernel = std::chrono::system_clock::now();
	gkernel_eR<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent);
	//auto end_kernel = std::chrono::system_clock::now();
	
	// check status of device and synchronize;
	size_t nEvents = EstnEvtMax;
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	// copy result of event reconstruction from device_gEvent to device_input_TKL
	// this input_tkl should be the information that the device uses to reconstruct the tracklets
	gpuErrchk( cudaMemcpy(device_input_TKL, device_gEvent, NBytesAllEvent, cudaMemcpyDeviceToDevice));

	// shouldn't this function actually be called? should it be the function that puts together tracklets? and then call the fitting???
	// gkernel_TKL<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_input_TKL, device_output_TKL);

	// check status of device and synchronize again;
	
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	// data transfer from device to host
	gpuErrchk( cudaMemcpy(host_output_eR, device_gEvent, NBytesAllEvent, cudaMemcpyDeviceToHost));
	gpuErrchk( cudaMemcpy(host_output_TKL, device_output_TKL, NBytesAllSearchWindow, cudaMemcpyDeviceToHost));


	// thrust objects: C++ template library based on STL
	// convert raw pointer device_gEvent to device_vector
	// TODO: just don't use raw pointers to begin with
    thrust::device_ptr<gEvent> d_p_events(device_gEvent);
    thrust::device_vector<gEvent> d_events(d_p_events, d_p_events + nEvents);
	std::vector<gEvent> h_events(nEvents);
	std::copy(d_events.begin(), d_events.end(), h_events.begin());

	thrust::device_vector<float> d_hit_pos(nEvents);
	// std::vector<float> h_hit_pos;

	// copy hit pos from event vector to dedicated hit pos vector
	// TODO: do this on the GPU instead (possibly using zip_iterator)
	// for (auto j = h_events.begin(); j < h_events.begin() + 100; ++j) {
	// 	// cout << "e " << j->EventID << endl;
	// 	for (auto i = 0; i < EstnAHMax; ++i) {
	// 		// float pos = static_cast<gEvent>(*j).AllHits[i].pos;
	// 		float pos = static_cast<gEvent>(*j).AllHits[i].driftDistance;
	// 		if (abs(pos) > Epsilon) {
	// 			// h_hit_pos.push_back(pos);
	// 			d_hit_pos.push_back(pos);
	// 			// cout  << " " << pos << endl;
	// 		}
	// 	}
	// }
	// thrust::copy(h_hit_pos.begin(), h_hit_pos.end(), d_hit_pos.begin());


	// thrust::transform(device_gEvent, device_gEvent + nEvents, d_hit_pos.begin(), get_first_event_hit_pos());
	// cout << "First event positions (10 / " << d_hit_pos.size() << "):";
	// thrust::copy(d_hit_pos.begin(), d_hit_pos.begin()+10, std::ostream_iterator<int>(std::cout, ", "));
	// cout <<  endl;

	// int NCheck = 10;
	// for (int i = 0; i < NCheck; ++i) {
	// 	gEvent& evt = (host_output_eR)[i];
	// 	gSW& sw = (host_output_TKL)[i];
	// 	cout <<  i << ". " << evt.EventID <<  ". " <<  evt.nAH << ", " << sw.EventID << ", " << sw.nAH << endl;
	// }
	
	
	
	// ###############################################################
	// Gpufit
	// ###############################################################

	thrust::device_vector< REAL > d_parameters(2);
	
	//data is array of xz-positions of v,v',x,x',u,u' planes of each tracklet
	vector< REAL > _data_tkl_x {-2.48f,-2.50f, -0.824f, -0.826f, -0.473f,-0.474f};
	thrust::device_vector<REAL> d_tkl_x(_data_tkl_x.size());
	thrust::copy(_data_tkl_x.begin(), _data_tkl_x.end(), d_tkl_x.begin());


	// true parameters fo xz view
	std::vector< REAL > true_parameters_x { 150, 0.15f }; // offset, slope
	thrust::copy(true_parameters_x.begin(), true_parameters_x.end(), d_parameters.begin());

	// linear_regression_example(d_hit_pos.size(), d_hit_pos.data().get());
	// calling linear regression with 6 points... 
	// it looks like we're fitting a single tracklet. 
	// Where is the part where we're getting all tracklet candidates and fit them? 
	// what about the part where we are fitting a full track?
	linear_regression_example(d_tkl_x.size(), d_tkl_x.data().get(), d_parameters);

	//data is array of yz-positions of v,v',x,x',y,y' planes of each tracklet

	vector< REAL > _data_tkl_y {-0.761f, -0.764f, -0.067f, -0.069f, -0.742f, -0.75f};
	thrust::device_vector<REAL> d_tkl_y(_data_tkl_y.size());
	thrust::copy(_data_tkl_y.begin(), _data_tkl_y.end(), d_tkl_y.begin());


	// true parameters fo yz view
	std::vector< REAL > true_parameters_y { 150, 0.15f }; // offset, slope
	thrust::copy(true_parameters_y.begin(), true_parameters_y.end(), d_parameters.begin());

	linear_regression_example(d_tkl_y.size(), d_tkl_y.data().get(), d_parameters);



	// cudaMemcpy(host_gEvent, device_gEvent, NBytesAllEvent, cudaMemcpyDeviceToHost);
	// cudaMemcpy(host_output, device_output_eR, sizeofoutput_eR, cudaMemcpyDeviceToHost);
	// cudaFree(device_gEvent);
	// // cudaFree(device_output_eR);
	// cudaFree(device_input_TKL);
	// cudaFree(device_output_TKL);

	//auto end_kernel = std::chrono::system_clock::now();
	clock_t cp3 = clock();

	delete rawEvent;

	//cout<<"output: "<<(device_output_eR)<<endl
	
	//for(int i = 0; i < host_gEvent[0].nAH; ++i) {
		  //cout<<"D0_1st_wire:" << (host_gEvent[0].NHits[1])<<endl;
		//cout<<"output: "<<(host_gEvent[0].nAH)<<endl;
		//cout<<"output: "<<(device_output_eR)<<endl;
		//cout<<"output: "<<(sizeof(int))<<endl;
		//cout<<"size: "<<i<<endl;
	//}
	// printing the time required for all operations
	clock_t cp4 = clock();
	auto end = std::chrono::system_clock::now();

	//double cpu_secs = double(cp4-cp3+cp2-cp1) / CLOCKS_PER_SEC;
	double cpu_secs = double(cp2-cp1) / CLOCKS_PER_SEC;
	auto gpu_ns = end_kernel - start_kernel;
	auto overall = end - start;
	cout<<"Read/prepare events: "<<cpu_secs<<endl;
	cout<<"eR "<<(gpu_ns.count()/1000000000.0)<<endl;
	cout<<"Total time: "<<(overall.count()/1000000000.0)<<endl;

	return 0;
}

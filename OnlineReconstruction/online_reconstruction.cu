#include "gpufit.h"
#include "interface.h"
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <fstream>
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
//#include "LoadInput.h"
#include "OROutput.h"
//#include "trackfittingtools.cuh"
#include "tracknumericalminimizer.cuh"
//#include "cuda_fit_device_funcs.cuh"
//#include "gpufit_core_funcs.cuh"

#include "SQEvent_v1.h"
#include "SQHit_v1.h"
#include "SQHitVector_v1.h"

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
const int TrackletSizeMax = 200;

const double TX_MAX = 0.15;
const double TY_MAX = 0.1;
const double X0_MAX = 150;
const double Y0_MAX = 50;
const double INVP_MAX = 0.2;
const double INVP_MIN = 0.01;


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

class gTracklet {
      public:
      gTracklet(){
	nXHits = nUHits = nVHits = 0;
      }
      
      short stationID;
      short nXHits;
      short nUHits;
      short nVHits;

      float chisq;
      float chisq_vtx;

      //maybe we can be a bit more sober in memory, and just require hit "event" index?
      gHit hits[nDetectors];// array of all hits
      short hitsign[nDetectors];
      
      float tx;
      float ty;
      float x0;
      float y0;
      float invP;
      
      float err_tx;
      float err_ty;
      float err_x0;
      float err_y0;
      float err_invP;
      
      float residual[nChamberPlanes];
};

class gFitParams {
public:
      int max_iterations;
      short nparam;
      
      float parameter_limits_min[5];
      float parameter_limits_max[5];
      
      float tolerance;
};




class gFitArrays {
public:
      int npoints;
      float drift_dist[nChamberPlanes]; // hit drift distance
      float resolution[nChamberPlanes]; // detector resolution
      
      float p1x[nChamberPlanes];// x bottom end point of the wire hit 
      float p1y[nChamberPlanes];// y bottom end point of the wire hit 
      float p1z[nChamberPlanes];// z bottom end point of the wire hit 
      
      float deltapx[nChamberPlanes];// x distance between bottom and top end points of the wire hit 
      float deltapy[nChamberPlanes];// y distance between bottom and top end points of the wire hit 
      float deltapz[nChamberPlanes];// z distance between bottom and top end points of the wire hit 
      
      float prev_parameters[5];
      float output_parameters[5];
      float output_parameters_errors[5];
      float chi2prev;
      float chi2;

      float values[nChamberPlanes];
      float derivatives[5*nChamberPlanes];
      float gradients[5];
      float hessians[25];
      
      float scaling_vector[5];
      float deltas[5];
      
      float lambda;//scale factor
      
      //float calc_matrix[25];
      //float abs_row[25];
      //int abs_row_index[25];

      int n_iter;
      short iter_failed;
      short finished;
      short state;
      short skip;
      short singular;
      
      //float x_array[nChamberPlanes];// x position arrays
      //float y_array[nChamberPlanes];// y position arrays
      //float z_array[nChamberPlanes];// z position arrays
      //float dx_array[nChamberPlanes];// x position uncertainty
      //float dy_array[nChamberPlanes];// x position uncertainty

      float A[25];// matrix: max size 5x5, but we can use this unique array for all possible sizes
      float Ainv[25];// matrix
      float B[5];// input vector

      //float output_parameters_steps[5];
      //float doublederivatives[5];
};



class gEvent {
	public:
	gEvent(){
		RunID = EventID = SpillID = -1;
	}
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

//Output class
class gSW {
public:
	gSW(){
		EventID = -1;
		nAH = 0;
		nTracklets = 0;
	}
	thrust::pair<int, int> hitpairs_x[100];
	thrust::pair<int, int> hitpairs_u[100];
	thrust::pair<int, int> hitpairs_v[100];
	int hitidx1[100];
	int hitidx2[100];
	short hitflag1[100];
	short hitflag2[100];
	int EventID;
	int nAH;
	int nTracklets;
	//really tempted to replace tracklet array with array of IDs
	gTracklet AllTracklets[TrackletSizeMax];
	short nTKL_stID[7];//0: D0; 1: D1; 2: D2; 3: D3p; 4: D3m; 5: back partial; 6: global
};

//geometry carrier
class gPlane {
      public:
      float z;
      int nelem;
      float spacing;
      float xoffset;
      float scalex;
      float x0;
      float costheta;
      float scaley;
      float y0;
      float sintheta;
      float resolution;
      float deltaW_[9];
      float z_mean;
      float u_win;
      float v_win_fac1;
      float v_win_fac2;
      float v_win_fac3;
      float p1x_w1;
      float p1y_w1;
      float p1z_w1;
      float deltapx;
      float deltapy;
      float deltapz;
      float dp1x;
      float dp1y;
      float dp1z;
      float spacingplane[28];
      short detsuperid[7][3];
      short dets_x[6];

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


// kernel functions: 
// CUDA C++ extends C++ by allowing the programmer to define C++ functions, called kernels, that, when called, 
// are executed N times in parallel by N different CUDA threads, as opposed to only once like regular C++ functions. 

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
	

	//if(ic[index].EventID==0){
	//	printf("evt = %d, nAH = %d: \n", ic[index].EventID, ic[index].nAH);
	//	for(int i = 1; i<=nDetectors; i++)printf(" det %d: %d;  ", i, ic[index].NHits[i]);
	//	printf("\n");
	//}

	// initialization of array size
	cluster_iAH_arr_size[index] = 0;
	nAH_reduced[index] = 0;
	
	// event reducing/hit filtering
	for(iAH[index] = 0; iAH[index]<ic[index].nAH; ++iAH[index]) {
		// if hit not good, set its detID to 0 and continue;
		if((ic[index].AllHits[iAH[index]].flag & hitFlagBit(1)) == 0) {
			//if(ic[index].EventID==0)printf("hit det %d Skip out-of-time...\n", ic[index].AllHits[iAH[index]].detectorID);
			ic[index].AllHits[iAH[index]].detectorID = 0;
			continue;
		}
		uniqueID[index] = uniqueID_curr[index] = -1;
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
					//if(ic[index].EventID==0)printf("hit det %d el %d Skip after-pulse...\n", ic[index].AllHits[iAH[index]].detectorID, ic[index].AllHits[iAH[index]].elementID);
					ic[index].AllHits[iAH[index]].detectorID = 0;
					continue;
				}
				else {
					tdcTime_curr[index] = ic[index].AllHits[iAH[index]].tdcTime;
				}
			}
		}
		// declustering of hits in DCs (from CPU code, I understand this one better)
		// if there are hits in the same plane and hitting to neighboring wires, they both give redundant information: 
		if(ic[index].AllHits[iAH[index]].detectorID <= nChamberPlanes) {
			//if(ic[index].EventID==0)printf("%d\n", cluster_iAH_arr_size[index]);
//			printf("Decluster...\n");
			if(cluster_iAH_arr_size[index] == ClusterSizeMax) {
//				printf("Oversized cluster...\n");
			}
			// if array size is zero, start storing the hit in the array
			if(cluster_iAH_arr_size[index] == 0) {
				cluster_iAH_arr[index][0] = iAH[index];
				++cluster_iAH_arr_size[index];
			} else { // otherwise
				// current hit and previous hit are *not* in same detector plane OR next hit and current hit are *not* in neighbors cells
				// we "declusterize" i.e. we remove the hit/hits which information is redundant with other hits and/or useless
				//if(ic[index].EventID==0){
//printf("hit indices: %d %d %d\n", iAH[index], cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1], cluster_iAH_arr[index][0]);
//printf("hit det/elem: %d, %d; %d, %d\n", ic[index].AllHits[iAH[index]].detectorID, ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID, 
//	    	      	      	  	 ic[index].AllHits[iAH[index]].elementID, ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID);
//printf("diffs: %d, %d\n", (ic[index].AllHits[iAH[index]].detectorID - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID),
//	       	   	  (ic[index].AllHits[iAH[index]].elementID - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID));
//printf("bools: %d, %d\n", (ic[index].AllHits[iAH[index]].detectorID != ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID),
//	       	   	  (abs(ic[index].AllHits[iAH[index]].elementID - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID) > 1));
	    	     	   	//}
				if((ic[index].AllHits[iAH[index]].detectorID != ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID) || (abs(ic[index].AllHits[iAH[index]].elementID - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID) > 1)) {
					// if 2 hits in cluster, evaluate w_max and w_min; drift distance has to be < w_min for one of the hits, while it has to be < w_max for the other hit 
					if(cluster_iAH_arr_size[index] == 2) {
						w_max[index] = 0.9*0.5*(ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].pos - ic[index].AllHits[cluster_iAH_arr[index][0]].pos);
						w_min[index] = 4.0/9.0*w_max[index];
						if((ic[index].AllHits[cluster_iAH_arr[index][0]].driftDistance > w_max[index] && ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].driftDistance > w_min[index]) || (ic[index].AllHits[cluster_iAH_arr[index][0]].driftDistance > w_min[index] && ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].driftDistance > w_max[index])) {
						//if(ic[index].EventID==0)printf("hit indices: %d %d %d\n", iAH[index], cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1], cluster_iAH_arr[index][0]);
							//eliminating the existing hit with the lagest drift distance
							if(ic[index].AllHits[cluster_iAH_arr[index][0]].driftDistance > ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].driftDistance) {
								//if(ic[index].EventID==0)printf("1 - hit det %d elem %d Skip cluster...\n", ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID, ic[index].AllHits[cluster_iAH_arr[index][0]].elementID);
								ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID = 0;
							}
							else {
								//if(ic[index].EventID==0)printf("2 - hit det %d elem %d Skip cluster...\n", ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID, ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID);
								ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID = 0;
							}
						}
						// if the time difference is less than 8 ns for detectors 19 to 24 (which btw are DC3p) we remove both
						else if((((ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) >= 0.0 && (ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) < 8.0) || ((ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) <= 0.0 && (ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) > -8.0)) && (ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID >= 19 && ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID <= 24)) {
						        //if(ic[index].EventID==0)printf("3 - hit det %d elem %d Skip cluster...\n", ic[index].AllHits[iAH[index]].detectorID, ic[index].AllHits[iAH[index]].elementID);
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
						        //if(ic[index].EventID==0)printf("4 - hit det %d elem %d Skip cluster...\n", ic[index].AllHits[iAH[index]].detectorID, ic[index].AllHits[iAH[index]].elementID);
							for(cluster_iAH_arr_cur[index] = 0; cluster_iAH_arr_cur[index] < cluster_iAH_arr_size[index]; ++cluster_iAH_arr_cur[index]) {
								ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].detectorID = 0;
							}
						}
						// otherwise, we remove them all except first and last
						else {
						        //if(ic[index].EventID==0)printf("5 - hit det %d elem %d Skip cluster...\n", ic[index].AllHits[iAH[index]].detectorID, ic[index].AllHits[iAH[index]].elementID);
							for(cluster_iAH_arr_cur[index] = 1; cluster_iAH_arr_cur[index] < cluster_iAH_arr_size[index]; ++cluster_iAH_arr_cur[index]) {
								ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].detectorID = 0;
							}
						}
					}
					cluster_iAH_arr_size[index] = 0;
				} else { // if hits are in the same detector and in two neighboring cells!
				        // current hit and previous hit are in same detector plane and in neighbor wires: 
				  	// we count how many hits we have in this case, until we find a hit in a different detector or in a wire that is not a neighbor to the previous hit.
				  	//if(ic[index].EventID==0)printf("h1: det %d el %d; h2 det %d el %d \n", ic[index].AllHits[iAH[index]].detectorID, ic[index].AllHits[iAH[index]].elementID, ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID, ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID);
				  	cluster_iAH_arr[index][cluster_iAH_arr_size[index]] = iAH[index];
				  	++cluster_iAH_arr_size[index];
				}
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
	
	//if(ic[index].EventID==0){
	//	printf("evt = %d, nAH = %d: \n", ic[index].EventID, ic[index].nAH);
	//	for(int i = 1; i<=nDetectors; i++)printf(" det %d: %d;  ", i, ic[index].NHits[i]);
	//	printf("\n");
	//}

	//if(((ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6])>0) || ((ic[index].NHits[7]+ic[index].NHits[8]+ic[index].NHits[9]+ic[index].NHits[10]+ic[index].NHits[11]+ic[index].NHits[12])>0) || ((ic[index].NHits[13]+ic[index].NHits[14]+ic[index].NHits[15]+ic[index].NHits[16]+ic[index].NHits[17]+ic[index].NHits[18])>0) || ((ic[index].NHits[19]+ic[index].NHits[20]+ic[index].NHits[21]+ic[index].NHits[22]+ic[index].NHits[23]+ic[index].NHits[24])>0) || ((ic[index].NHits[25]+ic[index].NHits[26]+ic[index].NHits[27]+ic[index].NHits[28]+ic[index].NHits[29]+ic[index].NHits[30])>0)){	
	//if(((ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6])<270) || ((ic[index].NHits[7]+ic[index].NHits[8]+ic[index].NHits[9]+ic[index].NHits[10]+ic[index].NHits[11]+ic[index].NHits[12])>350) || ((ic[index].NHits[13]+ic[index].NHits[14]+ic[index].NHits[15]+ic[index].NHits[16]+ic[index].NHits[17]+ic[index].NHits[18])>170) || ((ic[index].NHits[19]+ic[index].NHits[20]+ic[index].NHits[21]+ic[index].NHits[22]+ic[index].NHits[23]+ic[index].NHits[24])>140) || ((ic[index].NHits[25]+ic[index].NHits[26]+ic[index].NHits[27]+ic[index].NHits[28]+ic[index].NHits[29]+ic[index].NHits[30])>140))

	//we do not accept the event unless there is at least one hit in the first DC

	/*
	if( (ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6])<1){
		//printf("Event rejected...\n");
		}
		else {
			//counting total hit number, for all events < 6668? why? because she wanted just a subset!
			if( (ic[index].EventID)>10000 && (ic[index].EventID)<10050 ){//just look at a subset with something in it
				int totalDetectorHits = 0;
				for(int i = 1; i <= nDetectors; ++i) {
					totalDetectorHits += ic[index].NHits[i];
					//printf("%d ", ic[index].NHits[i]);
				}//printf("\n");
				
				//int nFirstRegionHits = 0;
				//for(int i = 1; i < 6; ++i) {
					//nFirstRegionHits += ic[index].NHits[i];
					//printf("nHits[%d] = %d\n", i, ic[index].NHits[i]);
				//}
			
				//printf("AllHits value : %d\n", (ic[index].NHits[0]));
				//printf("event : %d\n", (ic[index].EventID));
				//printf("all detector hit sum : %d; reduced AllHits value : %d\n", totalDetectorHits, (nAH_reduced[index]));
				//printf("sum of detectors : %d (%d)\n", totalDetectorHits, nFirstRegionHits);
				//}
			}
		    
	//		Process the accepted events (tracking) here.
			// where is the tracking though?
		}
	*/
}

// function to match a tracklet to a hodoscope hit
__device__ int match_tracklet_to_hodo(gTracklet tkl, int stID, gEvent* ic, gPlane* planes)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	int masked = 0;//false
	// first, define the search region, and foremost, the planes on which we define this search region, which depends on the station ID we're looking at
	int hodoplanerange[5][2] = {{31, 34}, {31, 34}, {35, 38}, {39, 42}, {39, 42}};// range of planes to look for hits
	REAL fudgefac[5] = {0.25f, 0.25f, 0.2f, 0.15f, 0.15f};
	// define the region in which we are supposed to have hits:
	
	REAL xhodo, yhodo, xmin, xmax, ymin, ymax;
	
	// loop on the hits and select hodoscope hits corresponding to the station
	for(int i = 0; i<ic[index].nAH; i++){
		//we only consider hits in the hodoscopes planes corresponding to the station where the tracklet is reconstructed 
		if(hodoplanerange[stID][0]>ic[index].AllHits[i].detectorID || hodoplanerange[stID][1]<ic[index].AllHits[i].detectorID)continue;
		
		//calculate the track position at the hodoscope plane z
		xhodo = planes[ic[index].AllHits[i].detectorID-1].z*tkl.tx+tkl.x0;
		yhodo = planes[ic[index].AllHits[i].detectorID-1].z*tkl.ty+tkl.y0;
		
		//calculate "xmin, xmax, ymin, ymax" in which the track is supposed to pass through; 
		//these are basicially defined as the spatial coverage of the hit hodoscope element (plus some fudge factor for x)
		if(planes[ic[index].AllHits[i].detectorID-1].costheta>0.99){
			xmin = ic[index].AllHits[i].pos-planes[ic[index].AllHits[i].detectorID-1].spacing*(0.5f+fudgefac[stID])-planes[ic[index].AllHits[i].detectorID-1].z*tkl.err_tx+tkl.err_x0;
			xmax = ic[index].AllHits[i].pos+planes[ic[index].AllHits[i].detectorID-1].spacing*(0.5f+fudgefac[stID])+planes[ic[index].AllHits[i].detectorID-1].z*tkl.err_tx+tkl.err_x0;
			ymax = 0.5f*planes[ic[index].AllHits[i].detectorID-1].scaley+planes[ic[index].AllHits[i].detectorID-1].z*tkl.err_ty+tkl.err_y0;
			ymin = -ymax;
		}else{
			xmax = 0.5f*planes[ic[index].AllHits[i].detectorID-1].scalex+planes[ic[index].AllHits[i].detectorID-1].z*tkl.err_tx+tkl.err_x0;
			xmin = -xmax;
			ymin = ic[index].AllHits[i].pos-planes[ic[index].AllHits[i].detectorID-1].spacing-planes[ic[index].AllHits[i].detectorID-1].z*tkl.err_ty+tkl.err_y0;
			ymax = ic[index].AllHits[i].pos+planes[ic[index].AllHits[i].detectorID-1].spacing+planes[ic[index].AllHits[i].detectorID-1].z*tkl.err_ty+tkl.err_y0;
		}

		//fudge = planes[ic[index].AllHits[i].detectorID-1].spacing*fudgefac[stID];
		if(xmin <= xhodo && xhodo <= xmax && ymin <= yhodo && yhodo <= ymax ){
			masked++;
		}
		
	}
	// 
	return masked>0;
}


// function to make the hit pairs in station;
// I assume it will only be called by the tracklet builder
// (not by the main function), so I can make it a "device" function. 
__device__ int make_hitpairs_in_station(gEvent* ic, thrust::pair<int, int>* hitpairs, int* hitidx1, int* hitidx2, short* hitflag1, short* hitflag2, int stID, int projID, gPlane* planes){
	   //TODO: fix memory management: there's a lot of stuff declared here (though it's probably destroyed after so it *might* be fine)
	   int npairs = 0;
	   
	   // I think we assume that by default we want to know where we are
	   int index = threadIdx.x + blockIdx.x * blockDim.x;

	   //declaring arrays for the hit lists
	   for(int i = 0; i<100; i++){
	   	   hitidx1[i] = hitidx2[i] = 0;
		   hitflag1[i] = hitflag2[i] = 0;
	   }

	   //building the lists of hits for each detector plane
	   int detid1 = planes[0].detsuperid[stID][projID]*2;
	   int detid2 = planes[0].detsuperid[stID][projID]*2-1;
	   int superdetid = planes[0].detsuperid[stID][projID];
	   int hitctr1 = 0, hitctr2 = 0;
	   for(int i = 0; i<ic[index].nAH; i++){
	   	  if(ic[index].AllHits[i].detectorID==detid1){
			hitidx1[hitctr1] = i;
			hitctr1++;
		  }
	   	  if(ic[index].AllHits[i].detectorID==detid2){
			hitidx2[hitctr2] = i;
			hitctr2++;
		  }
	   }

	   // pair the hits by position:
	   // if one hit on e.g. x and one hit on x' are closer than
	   // the "spacing" defined for the planes, then the hits can be paired together.
	   int idx1 = -1;
	   int idx2 = -1;
	   for(int i = 0; i<hitctr1; i++){
	   	   idx1++;
		   idx2 = -1;
	   	   for(int j = 0; j<hitctr2; j++){
		   	   idx2++;
			   if( abs(ic[index].AllHits[ hitidx1[idx1] ].pos - ic[index].AllHits[ hitidx2[idx2] ].pos) > planes[0].spacingplane[superdetid] ){
			       continue;
			   }
			   	   
			   hitpairs[npairs] = thrust::make_pair(hitidx1[idx1], hitidx2[idx2]);
			   npairs++;
			   hitflag1[idx1] = 1;
			   hitflag2[idx2] = 1;
	   	   }
	   }
	   // here the hits that cannot be paired to another hit are paired to "nothing"
	   // (but they still have to be paired to be used in the trackletteing)
	   for(int i = 0; i<hitctr1; i++){
	   	   if(hitflag1[i]<1){
			hitpairs[npairs] = thrust::make_pair(hitidx1[i], -1);
			npairs++;
		   }
	   }
	   for(int i = 0; i<hitctr2; i++){
	   	   if(hitflag2[i]<1){
			hitpairs[npairs] = thrust::make_pair(hitidx2[i], -1);
			npairs++;
		   }
	   }
	   	   
	   return npairs;
}



// tracklet in station builder: 
__global__ void gkernel_TrackletinStation(gEvent* ic, gSW* oc, gFitArrays* fitarrays, int stID, gPlane* planes, gFitParams* fitparams) {
	//int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
	//int index = threadIdx.x + blockId * blockDim.x;
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	stID--;
	oc[index].EventID = ic[index].EventID;
	oc[index].nAH = ic[index].nAH;
	
	//if(10000<ic[index].EventID && ic[index].EventID<10050){
	//	for(int m = 0; m<30; m++){
	//		if(planes[m].u_win!=0)printf("index= %d, m = %d, u_win = %1.6f, costheta = %1.6f\n", index, m, planes[m].u_win, planes[m].costheta);
	//	}
	//}
	// loop on hits
	//if( (ic[index].EventID)>10000 && (ic[index].EventID)<10100 ){//just look at a subset with something in it
	//	printf("core idx %d, evt %d: reduced AllHits value : %d\n", (index), ic[index].EventID, (Nhits));
	//}
	// answer is yes, we still have the info from the previous function i.e. running after eR we still benefit from hit reduction;
	// was worth checking, just in case...

	//we don't need pairs of *HITS* necessarily, we just need pairs of indices...
	//thrust::pair<int, int> hitpairs_x[100];
	//thrust::pair<int, int> hitpairs_u[100];
	//thrust::pair<int, int> hitpairs_v[100];

	int nx = make_hitpairs_in_station(ic, oc[index].hitpairs_x, oc[index].hitidx1, oc[index].hitidx2, oc[index].hitflag1, oc[index].hitflag2, stID, 0, planes);
	int nu = make_hitpairs_in_station(ic, oc[index].hitpairs_u, oc[index].hitidx1, oc[index].hitidx2, oc[index].hitflag1, oc[index].hitflag2, stID, 1, planes);
	int nv = make_hitpairs_in_station(ic, oc[index].hitpairs_v, oc[index].hitidx1, oc[index].hitidx2, oc[index].hitflag1, oc[index].hitflag2, stID, 2, planes);
	
	short uidx = stID==0? stID*6 : stID*6+4;
	short vidx = stID==0? stID*6+4 : stID*6;
	
	bool print = false;
	//if(0<=ic[index].EventID && ic[index].EventID<20){
	//	print = true;
	//printf("evt %d, nx = %d, nu = %d, nv = %d, ucostheta(plane %d) = %1.6f, uwin(plane %d) = %1.6f\n", ic[index].EventID, nx, nu, nv, uidx, planes[uidx].costheta, uidx, planes[uidx].u_win);
	//	printf("evt %d, nx = %d, nu = %d, nv = %d\n", ic[index].EventID, nx, nu, nv);
	//}
	//one has to have at least one hit in x, u, v
	if(nx==0 || nu==0 || nv==0)return;
	int n_tkl = oc[index].nTracklets;

	int nhits_tkl = 0;
	int npts = 0;

	REAL fixedpoint[3] = {0, 0, 0};
		
	//X-U combinations first
	for(int i = 0; i< nx; i++){
		double xpos = oc[index].hitpairs_x[i].second>=0 ? 0.5*(ic[index].AllHits[ oc[index].hitpairs_x[i].first ].pos+ic[index].AllHits[ oc[index].hitpairs_x[i].second ].pos): ic[index].AllHits[ oc[index].hitpairs_x[i].first ].pos;
		double umin = xpos*planes[uidx].costheta-planes[uidx].u_win;
		double umax = umin+2*planes[uidx].u_win;
		//if(print){
		//	printf("evt %d, xpos = %1.6f, umin = %1.6f, umax = %1.6f\n", ic[index].EventID, xpos, umin, umax);
		//	printf("evt %d, x1 pos = %1.6f, x2 pos =%1.6f\n", ic[index].EventID, ic[index].AllHits[ oc[index].hitpairs_x[i].first ].pos, 
		//		oc[index].hitpairs_x[i].second >=0 ? ic[index].AllHits[ oc[index].hitpairs_x[i].second ].pos : -1000000);
		//}
		for(int j = 0; j< nu; j++){
			double upos = oc[index].hitpairs_u[j].second>=0 ? 0.5*(ic[index].AllHits[ oc[index].hitpairs_u[j].first ].pos+ic[index].AllHits[ oc[index].hitpairs_u[j].second ].pos): ic[index].AllHits[ oc[index].hitpairs_u[j].first ].pos;

			if(upos<umin || upos>umax)continue;
			//if(print)printf("evt %d, %1.6f <? upos = %1.6f <? %1.6f \n", ic[index].EventID, umin, upos, umax);
			//we chose by convention to start the numbering of the "planes" object arrays to 0 instead of 1, this is why we have to subtract 1 to detectorID to get the correct information
			double z_x = oc[index].hitpairs_x[i].second>=0 ? planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].z_mean : planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].z;
			double z_u = oc[index].hitpairs_u[j].second>=0 ? planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].z_mean : planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].z;
			double z_v = planes[vidx].z_mean;
			//if(ic[index].EventID==0)printf("detid x = %d, detid u = %d, z_x = %1.6f, z_u = %1.6f, z_v = %1.6f\n", ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID, ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID, z_x, z_u, z_v);
			double v_win1 = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].v_win_fac1;
			double v_win2 = fabs(z_u+z_v-2*z_x)*planes[ vidx ].v_win_fac2;
			double v_win3 = fabs(z_v-z_u)*planes[ vidx ].v_win_fac3;
			double v_win = v_win1+v_win2+v_win3+2*planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].spacing;

			double vmin = 2*xpos*planes[uidx].costheta-upos-v_win;
			double vmax = vmin+2*v_win;
			//if(ic[index].EventID==0)printf("vmin = %1.6f, vmax = %1.6f, vwin = %1.6f, vwin1 = %1.6f, vwin2 = %1.6f, vwin3 = %1.6f\n", vmin, vmax, v_win, v_win1, v_win2, v_win3);
			for(int k = 0; k< nv; k++){
				double vpos = oc[index].hitpairs_v[k].second>=0 ? 0.5*(ic[index].AllHits[ oc[index].hitpairs_v[k].first ].pos+ic[index].AllHits[ oc[index].hitpairs_v[k].second ].pos): ic[index].AllHits[ oc[index].hitpairs_v[k].first ].pos;
				//if(ic[index].EventID<20)printf("evt %d: vmin = %1.6f <? vpos = %1.6f <? vmax = %1.6f\n", ic[index].EventID, vmin, vpos, vmax);
				if(vpos<vmin || vpos>vmax)continue;

				oc[index].AllTracklets[n_tkl].stationID = stID;
				nhits_tkl = 0;
				npts = 0;

				if(oc[index].hitpairs_x[i].first>=0){
					if(ic[index].EventID==0)printf("x first hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID, ic[index].AllHits[ oc[index].hitpairs_x[i].first ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_x[i].first ];
					oc[index].AllTracklets[n_tkl].nXHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;
					
					//fitarrays[index].x_array[npts] = ic[index].AllHits[ oc[index].hitpairs_x[i].first ].pos; // costheta = 1.
					//fitarrays[index].dx_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].spacing*0.5f; // costheta = 1.
					//fitarrays[index].y_array[npts] = 0.; // sintheta = 0.
					//fitarrays[index].dy_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].scaley*0.5f; // sintheta = 0.
					//fitarrays[index].z_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].z;
					
					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_x[i].first ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_x[i].first ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_x[i].first ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_x[i].first ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID-1 ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_x[i].second>=0){
					if(ic[index].EventID==0)printf("x second hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID, ic[index].AllHits[ oc[index].hitpairs_x[i].second ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_x[i].second ];
					oc[index].AllTracklets[n_tkl].nXHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;

					//fitarrays[index].x_array[npts] = ic[index].AllHits[ oc[index].hitpairs_x[i].second ].pos; // costheta = 1.
					//fitarrays[index].dx_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].spacing*0.5f; // costheta = 1.
					//fitarrays[index].y_array[npts] = 0.0; // sintheta = 0.
					//fitarrays[index].dy_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].scaley*0.5f; // sintheta = 0.
					//fitarrays[index].z_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].z;

					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_x[i].second ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_x[i].second ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_x[i].second ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_x[i].second ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID-1 ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_u[j].first>=0){
					if(ic[index].EventID==0)printf("u first hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID, ic[index].AllHits[ oc[index].hitpairs_u[j].first ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_u[j].first ];
					oc[index].AllTracklets[n_tkl].nUHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;
					
					//fitarrays[index].x_array[npts] = ic[index].AllHits[ oc[index].hitpairs_u[j].first ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].costheta;
					//fitarrays[index].dx_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].costheta);
					//fitarrays[index].y_array[npts] = ic[index].AllHits[ oc[index].hitpairs_u[j].first ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].sintheta;
					//fitarrays[index].dy_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].sintheta);
					//fitarrays[index].z_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].z;
					
					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_u[j].first ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_u[j].first ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[i].first ].detectorID-1 ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[i].first ].detectorID-1 ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_u[j].first ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_u[j].first ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID-1 ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_u[j].second>=0){
					if(ic[index].EventID==0)printf("u second hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID, ic[index].AllHits[ oc[index].hitpairs_u[j].second ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_u[j].second ];
					oc[index].AllTracklets[n_tkl].nUHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;

					//fitarrays[index].x_array[npts] = ic[index].AllHits[ oc[index].hitpairs_u[j].second ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].costheta;
					//fitarrays[index].dx_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].costheta);
					//fitarrays[index].y_array[npts] = ic[index].AllHits[ oc[index].hitpairs_u[j].second ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].sintheta;
					//fitarrays[index].dy_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].sintheta);
					//fitarrays[index].z_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].z;
					
					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_u[j].second ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_u[j].second ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[i].second ].detectorID-1 ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[i].second ].detectorID-1 ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_u[j].second ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_u[j].second ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID-1 ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_v[k].first>=0){
					if(ic[index].EventID==0)printf("v first hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID, ic[index].AllHits[ oc[index].hitpairs_v[k].first ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_v[k].first ];
					oc[index].AllTracklets[n_tkl].nVHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;
					
					//fitarrays[index].x_array[npts] = ic[index].AllHits[ oc[index].hitpairs_v[k].first ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].costheta;
					//fitarrays[index].dx_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].costheta);
					//fitarrays[index].y_array[npts] = ic[index].AllHits[ oc[index].hitpairs_v[k].first ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].sintheta;
					//fitarrays[index].dy_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].sintheta);
					//fitarrays[index].z_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].z;
					
					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_v[k].first ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_v[k].first ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[i].first ].detectorID-1 ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[i].first ].detectorID-1 ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_v[k].first ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_v[k].first ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID-1 ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_v[k].second>=0){
					if(ic[index].EventID==0)printf("v second hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID, ic[index].AllHits[ oc[index].hitpairs_v[k].second ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_v[k].second ];
					oc[index].AllTracklets[n_tkl].nVHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;
					
					//fitarrays[index].x_array[npts] = ic[index].AllHits[ oc[index].hitpairs_v[k].second ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].costheta;
					//fitarrays[index].dx_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].costheta);
					//fitarrays[index].y_array[npts] = ic[index].AllHits[ oc[index].hitpairs_v[k].second ].pos*planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].sintheta;
					//fitarrays[index].dy_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].spacing*0.5f/fabs(planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].sintheta);
					//fitarrays[index].z_array[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].z;
					
					fitarrays[index].drift_dist[npts] = 0;//ic[index].AllHits[ oc[index].hitpairs_v[k].second ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_v[k].second ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[i].second ].detectorID-1 ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[i].second ].detectorID-1 ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_v[k].second ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_v[k].second ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID-1 ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(npts<4)continue;
				if(oc[index].AllTracklets[n_tkl].nXHits<1 || oc[index].AllTracklets[n_tkl].nUHits<1 || oc[index].AllTracklets[n_tkl].nVHits<1)continue;
				

				fitarrays[index].output_parameters[0] = 0.;
				fitarrays[index].output_parameters[1] = 0.;
				fitarrays[index].output_parameters[2] = 0.;
				fitarrays[index].output_parameters[3] = 0.;

				get_straighttrack_fixedpoint(npts, fitarrays[index].drift_dist, fitarrays[index].resolution, 
								fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z, 
								fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz,
								fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B,
								fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, 
								fixedpoint, fitarrays[index].chi2);

				fitarrays[index].lambda = 0.01f;
				fitarrays[index].chi2prev = 1.e20;
				
				gpufit_algorithm_fitter(fitparams[0].max_iterations,
										fitarrays[index].n_iter, fitarrays[index].iter_failed, fitarrays[index].finished,
										fitarrays[index].state, fitarrays[index].skip, fitarrays[index].singular,
										fitparams[0].tolerance, fitparams[0].nparam,
										fitparams[0].parameter_limits_min, fitparams[0].parameter_limits_max,  
										fitarrays[index].output_parameters, fitarrays[index].prev_parameters, 
										fitarrays[index].chi2, fitarrays[index].chi2prev,
										fitarrays[index].values, fitarrays[index].derivatives, fitarrays[index].gradients,
										fitarrays[index].hessians, fitarrays[index].scaling_vector, 
										fitarrays[index].deltas, fitarrays[index].lambda,
										fitarrays[index].Ainv,
										//fitarrays[index].calc_matrix, fitarrays[index].abs_row, fitarrays[index].abs_row_index,
										npts, fitarrays[index].drift_dist, fitarrays[index].resolution, 
										fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z, 
										fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz);					
				
				oc[index].AllTracklets[n_tkl].x0 = fitarrays[index].output_parameters[0];
				oc[index].AllTracklets[n_tkl].y0 = fitarrays[index].output_parameters[1];
				oc[index].AllTracklets[n_tkl].tx = fitarrays[index].output_parameters[2];
				oc[index].AllTracklets[n_tkl].ty = fitarrays[index].output_parameters[3];
				oc[index].AllTracklets[n_tkl].err_x0 = fitarrays[index].output_parameters_errors[0];
				oc[index].AllTracklets[n_tkl].err_y0 = fitarrays[index].output_parameters_errors[1];
				oc[index].AllTracklets[n_tkl].err_tx = fitarrays[index].output_parameters_errors[2];
				oc[index].AllTracklets[n_tkl].err_ty = fitarrays[index].output_parameters_errors[3];
				oc[index].AllTracklets[n_tkl].chisq = fitarrays[index].chi2;
				
				/*
				if(ic[index].EventID==0)printf("track: x0 = %1.6f +- %1.6f, y0 = %1.6f +- %1.6f, tx = %1.6f +- %1.6f, ty = %1.6f +- %1.6f; chi2 = %1.6f\n", 
					oc[index].AllTracklets[n_tkl].x0, oc[index].AllTracklets[n_tkl].err_x0, 
					oc[index].AllTracklets[n_tkl].y0, oc[index].AllTracklets[n_tkl].err_y0, 
					oc[index].AllTracklets[n_tkl].tx, oc[index].AllTracklets[n_tkl].err_tx, 
					oc[index].AllTracklets[n_tkl].ty, oc[index].AllTracklets[n_tkl].err_ty,
					oc[index].AllTracklets[n_tkl].chisq);
				*/

				if(!match_tracklet_to_hodo(oc[index].AllTracklets[n_tkl], stID, ic, planes))continue;

				
				if(n_tkl>TrackletSizeMax){
					continue;
					//printf("evt %d: n_tkl = %d > %d\n", oc[index].EventID, n_tkl, TrackletSizeMax);
				}
				n_tkl++;
				oc[index].nTKL_stID[stID]++;

			}
		}
	}
	
	//if(print)printf("evt %d number of tracklets %d\n", oc[index].EventID, n_tkl);
	//printf("%d ", n_tkl);
	oc[index].nTracklets = n_tkl;
}

//
//device int match_backtrack_to_proptube(gTracklet tkl, gEvent* ic, gPlane* planes)
//{
//}

//device int resolve_leftright(gTracklet tkl)
//{
//}


// kernel to combine tracklets into back partial tracks
__global__ void gkernel_BackPartialTracks(gEvent* ic, gSW* oc, gFitArrays* fitarrays, gPlane* planes, gFitParams* fitparams){
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	short stID = 5;
	//for the time being, declare in function;
	REAL x_pos[4];
	REAL z_pos[4];
	REAL a, b, sum, det, sx, sy, sxx, syy, sxy;
	short nprop, iprop;
	REAL xExp;
	short nhits_X2, nhits_X3;
	int n_tkl = oc[index].nTracklets;
	int nhits_2 = 0, nhits_3 = 0;
	
	REAL fixedpoint[3] = {0, 0, 0};

	//tracklets in station 3  (D3p, stID 4-1, D3m, stID 5-1)
	for(int i = oc[index].nTKL_stID[2]; i<oc[index].nTKL_stID[2]+oc[index].nTKL_stID[3]+oc[index].nTKL_stID[4]; i++){
		nhits_2 = oc[index].AllTracklets[i].nXHits+oc[index].AllTracklets[i].nUHits+oc[index].AllTracklets[i].nVHits;
		nhits_X2 = 0;
		for(int k = 0; k<nhits_2; k++){
			if(oc[index].AllTracklets[i].hits[k].detectorID==planes[0].dets_x[0] || oc[index].AllTracklets[i].hits[k].detectorID==planes[0].dets_x[1])
			{
				x_pos[nhits_X2] = oc[index].AllTracklets[i].hits[k].pos;
				z_pos[nhits_X2] = planes[oc[index].AllTracklets[i].hits[k].detectorID].z;
				nhits_X2++;
			}
		}
		
		//tracklets in station 2 (D2, stID 3-1)
		for(int j = 0; j<oc[index].nTKL_stID[2]; j++){
			if(fabs(oc[index].AllTracklets[i].tx - oc[index].AllTracklets[j].tx) > TX_MAX || fabs(oc[index].AllTracklets[i].ty - oc[index].AllTracklets[i].ty) > TY_MAX)continue;
			
			nhits_3 = oc[index].AllTracklets[j].nXHits+oc[index].AllTracklets[j].nUHits+oc[index].AllTracklets[j].nVHits;
			nhits_X3 = nhits_X2;
			for(int k = 0; k<nhits_3; k++){
				if(oc[index].AllTracklets[j].hits[k].detectorID==planes[0].dets_x[2] || oc[index].AllTracklets[j].hits[k].detectorID==planes[0].dets_x[3] ||
			   	oc[index].AllTracklets[j].hits[k].detectorID==planes[0].dets_x[4] || oc[index].AllTracklets[j].hits[k].detectorID==planes[0].dets_x[5])
				{
					x_pos[nhits_X3] = oc[index].AllTracklets[i].hits[k].pos;
					z_pos[nhits_X3] = planes[oc[index].AllTracklets[i].hits[k].detectorID].z;
					nhits_X3++;
				}
			}
			
			chi2_simplefit(nhits_X2+nhits_X3, z_pos, x_pos, a, b, sum, det, sx, sy, sxx, syy, sxy);
			if(fabs(a)>X0_MAX || fabs(b)>TX_MAX)continue;
			//prop matching
			nprop = 0;
			for(short ip = 0; ip<4; ip++){
				iprop = 48+ip;
				xExp = a*planes[iprop-1].z+b;
				//loop on hits to find prop
				for(int k = 0; k<ic[index].nAH; k++){
					if(ic[index].AllHits[k].detectorID==iprop){
						if(fabs(ic[index].AllHits[k].pos-xExp)<5.08f){
							nprop++;
							break;
						}
					}
				}
				if(nprop>0)break;
			}
			if(nprop==0)continue;
				
			//Add tracklets: first, get number of hits in each
			oc[index].AllTracklets[n_tkl].stationID = stID;
			oc[index].AllTracklets[n_tkl].nXHits = oc[index].AllTracklets[i].nXHits+oc[index].AllTracklets[j].nXHits;
			oc[index].AllTracklets[n_tkl].nUHits = oc[index].AllTracklets[i].nUHits+oc[index].AllTracklets[j].nUHits;
			oc[index].AllTracklets[n_tkl].nVHits = oc[index].AllTracklets[i].nVHits+oc[index].AllTracklets[j].nVHits;
			
			//add in the new tracklet the hits from the two existing tracklets;
			//also fill the fit arrays
			for(int k = 0; k<nhits_2; k++){
				oc[index].AllTracklets[n_tkl].hits[k]=oc[index].AllTracklets[i].hits[k];
				oc[index].AllTracklets[n_tkl].hitsign[k]=oc[index].AllTracklets[i].hitsign[k];
				
				fitarrays[index].drift_dist[k] = 0.;//oc[index].AllTracklets[n_tkl].hits[k].driftDistance;
				fitarrays[index].resolution[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].spacing/3.4641f;//.resolution;
				fitarrays[index].p1x[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].p1x_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].dp1x * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1);
				fitarrays[index].p1y[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].p1y_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].dp1y * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1);
				fitarrays[index].p1z[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].p1z_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].dp1z * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1);
				fitarrays[index].deltapx[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].deltapx;
				fitarrays[index].deltapy[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].deltapy;
				fitarrays[index].deltapz[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].deltapz;
			}
				
			for(int k = nhits_2; k<nhits_2+nhits_3; k++){
				oc[index].AllTracklets[n_tkl].hits[k]=oc[index].AllTracklets[j].hits[k-nhits_2];
				oc[index].AllTracklets[n_tkl].hitsign[k]=oc[index].AllTracklets[j].hitsign[k-nhits_2];
					fitarrays[index].drift_dist[k] = 0.;//oc[index].AllTracklets[n_tkl].hits[k].driftDistance;
				fitarrays[index].resolution[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].spacing/3.4641f;//.resolution;
				fitarrays[index].p1x[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].p1x_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].dp1x * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1) ;
				fitarrays[index].p1y[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].p1y_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].dp1y * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1) ;
				fitarrays[index].p1z[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].p1z_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].dp1z * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1) ;
				fitarrays[index].deltapx[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].deltapx;
				fitarrays[index].deltapy[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].deltapy;
				fitarrays[index].deltapz[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID-1 ].deltapz;
			}
			
			//then fit:
				
			get_straighttrack_fixedpoint(nhits_2+nhits_3, fitarrays[index].drift_dist, fitarrays[index].resolution, 
							fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z, 
							fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz,
							fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B,
							fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, 
							fixedpoint, fitarrays[index].chi2);
			
			fitarrays[index].output_parameters[0] = b;
			fitarrays[index].output_parameters[2] = a;
			
			// do a rough estimation of the starting fit
			//fitarrays[index].output_parameters[0] = (oc[index].AllTracklets[i].x0+oc[index].AllTracklets[j].x0)*0.5f;
			//fitarrays[index].output_parameters[1] = (oc[index].AllTracklets[i].y0+oc[index].AllTracklets[j].y0)*0.5f;
			//fitarrays[index].output_parameters[2] = (oc[index].AllTracklets[i].tx+oc[index].AllTracklets[j].tx)*0.5f;
			//fitarrays[index].output_parameters[3] = (oc[index].AllTracklets[i].ty+oc[index].AllTracklets[j].ty)*0.5f;
			fitarrays[index].lambda = 0.01f;
			fitarrays[index].chi2prev = 1.e20;
			
			//include fit here:

			gpufit_algorithm_fitter(fitparams[0].max_iterations,
								fitarrays[index].n_iter, fitarrays[index].iter_failed, fitarrays[index].finished,
								fitarrays[index].state, fitarrays[index].skip, fitarrays[index].singular,
								fitparams[0].tolerance, fitparams[0].nparam,
								fitparams[0].parameter_limits_min, fitparams[0].parameter_limits_max,  
								fitarrays[index].output_parameters, fitarrays[index].prev_parameters, 
								fitarrays[index].chi2, fitarrays[index].chi2prev,
								fitarrays[index].values, fitarrays[index].derivatives, fitarrays[index].gradients,
								fitarrays[index].hessians, fitarrays[index].scaling_vector, 
								fitarrays[index].deltas, fitarrays[index].lambda,
								fitarrays[index].Ainv, 
								nhits_2+nhits_3, fitarrays[index].drift_dist, fitarrays[index].resolution, 
								fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z, 
								fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz);					
			
			/**/
			
			if(fitarrays[index].chi2>9000)continue;
			
			oc[index].AllTracklets[n_tkl].x0 = fitarrays[index].output_parameters[0];
			oc[index].AllTracklets[n_tkl].y0 = fitarrays[index].output_parameters[1];
			oc[index].AllTracklets[n_tkl].tx = fitarrays[index].output_parameters[2];
			oc[index].AllTracklets[n_tkl].ty = fitarrays[index].output_parameters[3];
			oc[index].AllTracklets[n_tkl].err_x0 = fitarrays[index].output_parameters_errors[0];
			oc[index].AllTracklets[n_tkl].err_y0 = fitarrays[index].output_parameters_errors[1];
			oc[index].AllTracklets[n_tkl].err_tx = fitarrays[index].output_parameters_errors[2];
			oc[index].AllTracklets[n_tkl].err_ty = fitarrays[index].output_parameters_errors[3];
			oc[index].AllTracklets[n_tkl].chisq = fitarrays[index].chi2;
			
			n_tkl++;
			oc[index].nTKL_stID[stID]++;
		}
	}
	
	oc[index].nTracklets = n_tkl;
	
}


int main(int argn, char * argv[]) {
	
	// initialization: declaration of SRaw event, opening file/tree, affecting rawEvent object to input tree
	// declaring array of gEvent;
	auto start = std::chrono::system_clock::now();
	clock_t cp1 = clock();

	TString inputFile;
	TString inputGeom;
	TString outputFile;
	inputFile = argv[1];
	inputGeom = argv[2];	
	outputFile = argv[3];

	//by default we should use e1039 
	bool e906data = false;
	if(argn>4)e906data = atoi(argv[4]);
	
	cout<<"Running "<<argv[0]<<endl;
	cout<<"Loading "<<argv[1]<<endl;
	cout<<"with geometry: "<<argv[2]<<endl;
	cout<<"Writing "<<argv[3]<<endl;
	
	//Get basic geometry here:
	double u_factor[5] = {5., 5., 5., 15., 15.};
	gPlane plane[nChamberPlanes+nHodoPlanes+nPropPlanes];

	REAL spacingplane_[28] = {0., 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 1.3, 1.3, 1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 4.0, 4.0, 7.0, 7.0, 8.0, 12.0, 12.0, 10.0, 3.0, 3.0, 3.0, 3.0};
	short detsuperid_[7][3] = {{2, 1, 3}, {5, 6, 4}, {8, 9, 7}, {11, 12, 10}, {14, 15, 13}, {25, 26, -1}, {24, 27, -1}}; 
	short dets_x_[6] = {15, 16, 21, 22, 27, 28};
	
	ifstream in_geom(inputGeom.Data());
  	string buffer;
	int ipl, nelem;
	float z, spacing, xoffset, scalex, x0, costheta, scaley, y0, sintheta, resolution, deltaW_;
	float p1x, p1y, p1z, deltapx, deltapy, deltapz, dp1x, dp1y, dp1z;
 	while ( getline(in_geom, buffer) ) {
    	      if (buffer[0] == '#') continue;
	      std::istringstream iss;
	      iss.str(buffer);
	      iss >> ipl >> z >> nelem >> spacing >> xoffset >> scalex >> x0 >> costheta >> scaley >> y0 >> sintheta >> resolution >> p1x >> p1y >> p1z >> deltapx >> deltapy >> deltapz >> dp1x >> dp1y >> dp1z;
	      plane[ipl-1].z = z;
	      plane[ipl-1].nelem = nelem;
	      plane[ipl-1].spacing = spacing;
	      plane[ipl-1].xoffset = xoffset;
	      plane[ipl-1].scalex = scalex;
	      plane[ipl-1].x0 = x0;
	      plane[ipl-1].costheta = costheta;
	      plane[ipl-1].scaley = scaley;
	      plane[ipl-1].y0 = y0;
	      plane[ipl-1].sintheta = sintheta;
	      plane[ipl-1].resolution = resolution;
	      plane[ipl-1].p1x_w1 = p1x;
	      plane[ipl-1].p1y_w1 = p1y;
	      plane[ipl-1].p1z_w1 = p1z;
	      plane[ipl-1].deltapx = deltapx;
	      plane[ipl-1].deltapy = deltapy;
	      plane[ipl-1].deltapz = deltapz;
	      plane[ipl-1].dp1x = dp1x;
	      plane[ipl-1].dp1y = dp1y;
	      plane[ipl-1].dp1z = dp1z;
	      if(ipl>nChamberPlanes+nHodoPlanes){
		for(int k = 0; k<9; k++){
			iss >> deltaW_;
			plane[ipl-1].deltaW_[k] = deltaW_;
		}
	      }else{
		iss >> deltaW_;
		plane[ipl-1].deltaW_[0] = deltaW_;
	      }
	      ipl++;
	}
	
	for(int i = 0; i<5; i++){
		int u_idx = i*6+4;
		if(i==0)u_idx = i*6;
		int x_idx = i*6+2;
		for(int j = 0; j<6; j++){
			int idx = i*6+j;
			plane[idx].z_mean = j%2==0 ? 0.5*(plane[idx].z+plane[idx+1].z):0.5*(plane[idx].z+plane[idx-1].z);
			
			plane[idx].v_win_fac1 = plane[idx].spacing*2*plane[u_idx].costheta;
			plane[idx].v_win_fac2 = plane[u_idx].costheta*TX_MAX;
			plane[idx].v_win_fac3 = fabs(plane[u_idx].sintheta*TY_MAX);
		}
		
		for(int j = 0; j<6; j++){
			int idx = i*6+j;
			plane[idx].u_win = fabs(0.5*plane[u_idx].scaley*plane[u_idx].sintheta) + TX_MAX*fabs((plane[u_idx].z_mean - plane[x_idx].z_mean)*plane[u_idx].costheta) + TY_MAX*fabs((plane[u_idx].z_mean - plane[x_idx].z_mean)*plane[u_idx].sintheta) + 2.*plane[u_idx].spacing + u_factor[i];
		}
		//cout << u_idx << " " << plane[u_idx].u_win << " = " << fabs(0.5*plane[u_idx].scaley*plane[u_idx].sintheta) << " + " << TX_MAX*fabs((plane[u_idx].z_mean - plane[x_idx].z_mean)*plane[u_idx].costheta) << " + " << TY_MAX*fabs((plane[u_idx].z_mean - plane[x_idx].z_mean)*plane[u_idx].sintheta) << " + " << 2.*plane[u_idx].spacing + u_factor[i] << endl;
		//cout << " u costheta " << plane[u_idx].costheta << " u sintheta " << plane[u_idx].sintheta << " x_span " << plane[u_idx].scaley << " spacing " << plane[u_idx].spacing << " z plane_u " << plane[u_idx].z_mean << " z plane_x " << plane[x_idx].z_mean << endl;  
	}
	cout << "Geometry file read out" << endl;
	
	//std::unordered_map<int, double> map_elemPosition[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
	double wire_position[54][400];//Let's keep this: simpler, more robust
	for(int i = 0; i < nChamberPlanes; ++i){
		for(int j = 0; j<28; j++){
			plane[i].spacingplane[j] = spacingplane_[j];
			if(j<6)plane[i].dets_x[j] = dets_x_[j];
			if(j<7){
				for(int k = 0; k<3; k++){
					plane[i].detsuperid[j][k] = detsuperid_[j][k];
				}
			}
		}
		//cout << plane[i].nelem << endl;
      		for(int j = 1; j <= plane[i].nelem; ++j){
          		double pos = (j - (plane[i].nelem+1.)/2.)*plane[i].spacing + plane[i].xoffset + plane[i].x0*plane[i].costheta + plane[i].y0*plane[i].sintheta + plane[i].deltaW_[0];
          		//map_elemPosition[i].insert(posType(j, pos));
			wire_position[i][j] = pos;
		}
	}
	for(int i = nChamberPlanes; i<nChamberPlanes+nHodoPlanes; ++i){
		for(int j = 0; j<28; j++){
			plane[i].spacingplane[j] = spacingplane_[j];
			if(j<6)plane[i].dets_x[j] = dets_x_[j];
			if(j<7){
				for(int k = 0; k<3; k++){
					plane[i].detsuperid[j][k] = detsuperid_[j][k];
				}
			}
		}
		//cout << plane[i].nelem << endl;
	      	for(int j = 1; j <= plane[i].nelem; ++j){
          		double pos = plane[i].x0*plane[i].costheta + plane[i].y0*plane[i].sintheta + plane[i].xoffset + (j - (plane[i].nelem+1)/2.)*plane[i].spacing + plane[i].deltaW_[0];
          		//map_elemPosition[i].insert(posType(j, pos));
			wire_position[i][j] = pos;
		}
	}
	for(int i = nChamberPlanes+nHodoPlanes; i<nChamberPlanes+nHodoPlanes+nPropPlanes; ++i){
		for(int j = 0; j<28; j++){
			plane[i].spacingplane[j] = spacingplane_[j];
			if(j<6)plane[i].dets_x[j] = dets_x_[j];
			if(j<7){
				for(int k = 0; k<3; k++){
					plane[i].detsuperid[j][k] = detsuperid_[j][k];
				}
			}
		}
		//cout << plane[i].nelem << endl;
	      	for(int j = 1; j <= plane[i].nelem; ++j){
          		int moduleID = 8 - int((j - 1)/8);
			//cout << moduleID << endl;
             		double pos = plane[i].x0*plane[i].costheta + plane[i].y0*plane[i].sintheta + plane[i].xoffset + (j - (plane[i].nelem+1)/2.)*plane[i].spacing + plane[i].deltaW_[moduleID];
          		//map_elemPosition[i].insert(posType(j, pos));
			wire_position[i][j] = pos;
		}
		
	}
	
	float par_limits_min[5] = {-X0_MAX, -Y0_MAX, -TX_MAX, -TY_MAX, INVP_MIN};
	float par_limits_max[5] = {+X0_MAX, +Y0_MAX, +TX_MAX, +TY_MAX, INVP_MAX};

	gFitParams fitparams[3];
	fitparams[0].max_iterations = 50;
	fitparams[1].max_iterations = 50;
	fitparams[2].max_iterations = 50;
	
	fitparams[0].nparam = 4;
	fitparams[1].nparam = 4;
	fitparams[2].nparam = 5;
	
	fitparams[0].tolerance = 0.001;
	fitparams[1].tolerance = 0.001;
	fitparams[2].tolerance = 0.001;

	for(int i = 0; i<3; i++){
		for(int j = 0; j<5; j++){
			fitparams[i].parameter_limits_min[j] = par_limits_min[j];
			fitparams[i].parameter_limits_max[j] = par_limits_max[j];
		}
	}
	
	TFile* dataFile = new TFile(inputFile.Data(), "READ");
	TTree* dataTree = 0;// = (TTree *)dataFile->Get("save");
	SRawEvent* rawEvent = new SRawEvent();
	//SQEvent_v1* event = new SQEvent_v1();
	//SQHitVector_v1* hitvec = new SQHitVector_v1();
	
	Int_t     _run_id;
 	Int_t     _spill_id;
	Int_t     _event_id;
   	UShort_t  _trigger;
	Int_t     _qie_presums[4];
	Int_t     _qie_turn_id = 0;
	Int_t     _qie_rf_id = 0;
	Int_t     _qie_rf_inte[33];
	
	std::vector<SQHit*> hit_vec;
	
	
	if(e906data){
		dataTree = (TTree *)dataFile->Get("save");
		dataTree->SetBranchAddress("rawEvent", &rawEvent);
	}else{
		//default option: e1039
		dataTree = (TTree *)dataFile->Get("T");
		dataTree->SetMakeClass(1); //this is necessary to get the tree to read the branch correctly
		dataTree->SetBranchStatus("*", 0);// this speeds up the loop on events

		dataTree->SetBranchStatus("DST.SQEvent._run_id", 1);
		dataTree->SetBranchStatus("DST.SQEvent._spill_id", 1);
		dataTree->SetBranchStatus("DST.SQEvent._event_id", 1);
		dataTree->SetBranchStatus("DST.SQEvent._trigger", 1);
		//dataTree->SetBranchStatus("DST.SQEvent._qie_presums[4]", 1);
		dataTree->SetBranchStatus("DST.SQEvent._qie_turn_id", 1);
		dataTree->SetBranchStatus("DST.SQEvent._qie_rf_id", 1);
		//dataTree->SetBranchStatus("DST.SQEvent._qie_rf_inte[33]", 1);
		
		dataTree->SetBranchAddress("DST.SQEvent._run_id", &_run_id);
		dataTree->SetBranchAddress("DST.SQEvent._spill_id", &_spill_id);
		dataTree->SetBranchAddress("DST.SQEvent._event_id", &_event_id);
		dataTree->SetBranchAddress("DST.SQEvent._trigger", &_trigger);
		//dataTree->SetBranchAddress("DST.SQEvent._qie_presums[4]", &_qie_presums);
		dataTree->SetBranchAddress("DST.SQEvent._qie_turn_id", &_qie_turn_id);
		dataTree->SetBranchAddress("DST.SQEvent._qie_rf_id", &_qie_rf_id);
		//dataTree->SetBranchAddress("DST.SQEvent._qie_rf_inte[33]", &_qie_rf_inte[33]);
		
		dataTree->SetBranchStatus("DST.SQHitVector._vector", 1);
		dataTree->SetBranchAddress("DST.SQHitVector._vector", &hit_vec);
	}
	int nEvtMax = dataTree->GetEntries();
	static gEvent host_gEvent[EstnEvtMax];

	cout << "unfolding " << nEvtMax <<" events" << endl;
	// loop on event: get RawEvent information and load it into gEvent
	for(int i = 0; i < nEvtMax; ++i) {
		if(i%1000==0)cout << i << "/" << nEvtMax <<  endl;
		dataTree->GetEntry(i);
		
		//cout<<"Converting "<<i<<"/"<<nEvtMax<<endl;
		if(e906data){
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
				//host_gEvent[i].AllHits[m].pos=map_elemPosition[(rawEvent->fAllHits[m]).detectorID-1][(rawEvent->fAllHits[m]).elementID];
				host_gEvent[i].AllHits[m].pos=wire_position[(rawEvent->fAllHits[m]).detectorID-1][(rawEvent->fAllHits[m]).elementID];
				host_gEvent[i].AllHits[m].flag=(rawEvent->fAllHits[m]).flag;
			}
			for(int n=0; n<rawEvent->fTriggerHits.size(); n++) {
				host_gEvent[i].TriggerHits[n].index=(rawEvent->fTriggerHits[n]).index;
				host_gEvent[i].TriggerHits[n].detectorID=(rawEvent->fTriggerHits[n]).detectorID;
				host_gEvent[i].TriggerHits[n].elementID=(rawEvent->fTriggerHits[n]).elementID;
				host_gEvent[i].TriggerHits[n].tdcTime=(rawEvent->fTriggerHits[n]).tdcTime;
				host_gEvent[i].TriggerHits[n].driftDistance=(rawEvent->fTriggerHits[n]).driftDistance;
				//host_gEvent[i].TriggerHits[n].pos=map_elemPosition[(rawEvent->fAllHits[n]).detectorID-1][(rawEvent->fAllHits[n]).elementID];
				host_gEvent[i].TriggerHits[n].pos=wire_position[(rawEvent->fAllHits[n]).detectorID-1][(rawEvent->fAllHits[n]).elementID];
				host_gEvent[i].TriggerHits[n].flag=(rawEvent->fTriggerHits[n]).flag;
			}
			// printouts for test
			//if(10000<rawEvent->fEventID&&rawEvent->fEventID<10050){
			//	printf("%d:\n ", rawEvent->fEventID);
			//	for(int l = 1; l<=nChamberPlanes; l++){
			//		printf("%d ", rawEvent->fNHits[l]);
			//	}printf("; %d\n", rawEvent->fAllHits.size());
			//	for(int m = 0; m<=50; m++){
			//		printf("%d, %1.3f;", (rawEvent->fAllHits[m]).detectorID, (rawEvent->fAllHits[m]).pos);
			//	}printf("\n");
			//}
		}else{
			//Default option: e1039
			//if(_event_id<20)cout << " evt: " << _event_id << " nhits = " << hit_vec.size() << endl; 
			host_gEvent[i].RunID = _run_id;
			host_gEvent[i].SpillID = _spill_id;
			host_gEvent[i].EventID = _event_id;
			host_gEvent[i].TriggerBits = _trigger;
			//for(int k = 0; k<4; k++)host_gEvent[i].NRoads[k] = _qie_presums[k];
			host_gEvent[i].TurnID = _qie_turn_id;
			host_gEvent[i].RFID = _qie_rf_id;
			//for(int k = 0; k<33; k++)host_gEvent[i].Intensity[k] = _qie_rf_inte[k];
			
			for(int k = 0; k<nDetectors; k++)host_gEvent[i].NHits[k] = 0;//we will increment those in the vector hit loop
			int ntrighits = 0;
			host_gEvent[i].nAH = hit_vec.size();
			for(int m = 0; m<hit_vec.size(); m++){
				if(hit_vec[m]->get_detector_id()>54){
					//if(_event_id<20)cout << " dark photon plane hit! " << hit_vec[m]->get_detector_id() << endl;  
					continue;
					//dark photon planes; I don't think we care about those for the purpose of online reconstruction... do we?
				}
				host_gEvent[i].NHits[hit_vec[m]->get_detector_id()]++;
				host_gEvent[i].AllHits[m].index=hit_vec[m]->get_hit_id();
				host_gEvent[i].AllHits[m].detectorID=hit_vec[m]->get_detector_id();
				host_gEvent[i].AllHits[m].elementID=hit_vec[m]->get_element_id();
				host_gEvent[i].AllHits[m].tdcTime=hit_vec[m]->get_tdc_time();
				host_gEvent[i].AllHits[m].driftDistance=fabs(hit_vec[m]->get_drift_distance());
				host_gEvent[i].AllHits[m].pos=wire_position[hit_vec[m]->get_detector_id()-1][hit_vec[m]->get_element_id()];
				host_gEvent[i].AllHits[m].flag=(1<<hit_vec[m]->is_in_time());
				//if(host_gEvent[i].EventID<20)cout << " det " << host_gEvent[i].AllHits[m].detectorID << " elem " << host_gEvent[i].AllHits[m].elementID << " time " << host_gEvent[i].AllHits[m].tdcTime << " dd " << host_gEvent[i].AllHits[m].driftDistance << " pos " << host_gEvent[i].AllHits[m].pos << endl;
				if(hit_vec[m]->is_trigger_mask()){
					host_gEvent[i].TriggerHits[ntrighits].index=hit_vec[m]->get_hit_id();
					host_gEvent[i].TriggerHits[ntrighits].detectorID=hit_vec[m]->get_detector_id();
					host_gEvent[i].TriggerHits[ntrighits].elementID=hit_vec[m]->get_element_id();
					host_gEvent[i].TriggerHits[ntrighits].tdcTime=hit_vec[m]->get_tdc_time();
					host_gEvent[i].TriggerHits[ntrighits].driftDistance=fabs(hit_vec[m]->get_drift_distance());
					host_gEvent[i].TriggerHits[ntrighits].pos=wire_position[hit_vec[m]->get_detector_id()-1][hit_vec[m]->get_element_id()];
					host_gEvent[i].TriggerHits[ntrighits].flag=(1<<hit_vec[m]->is_in_time());
					ntrighits++;
				}
			}
			host_gEvent[i].nTH = ntrighits;
		}
	}
	cout << "loaded events" << endl;
	clock_t cp2 = clock();
	
//	for(int i = 0; i < nEvtMax; ++i) {
//		thrust::stable_sort(host_gEvent[i].AllHits, host_gEvent[i].AllHits+host_gEvent[i].nAH, lessthan());
//	}

	clock_t cp3 = clock();
	// evaluate the total size of the gEvent array (and the SW array) for memory allocation 
	// (the memory cannot be dynamically allocated) 
	size_t NBytesAllEvent = EstnEvtMax * sizeof(gEvent);
	size_t NBytesAllSearchWindow = EstnEvtMax * sizeof(gSW);
	size_t NBytesAllPlanes =  nDetectors * sizeof(gPlane);
	size_t NBytesAllFitters = EstnEvtMax * sizeof(gFitArrays);
	size_t NBytesAllFitParam = sizeof(gFitParams)*3;// just to be sure: tracklets, back tracks, global tracks
	
	cout << NBytesAllEvent << " " << NBytesAllSearchWindow << " " << NBytesAllPlanes << " " << NBytesAllFitters << " " << NBytesAllFitParam << endl;
	cout << NBytesAllEvent+NBytesAllSearchWindow+NBytesAllPlanes+NBytesAllFitters+NBytesAllFitParam << endl;
	
	gEvent *host_output_eR = (gEvent*)malloc(NBytesAllEvent);
	gSW *host_output_TKL = (gSW*)malloc(NBytesAllSearchWindow);
	
	// declaring gEvent objects for the device (GPU) to use.
	gEvent *device_gEvent;
	// gEvent *device_output_eR;
	// gEvent *device_input_TKL;
	gSW *device_output_TKL;
	gPlane *device_gPlane;
	gFitParams *device_gFitParams;
	gFitArrays *device_gFitArrays;
	
	//printDeviceStatus();

	// copy of data from host to device: evaluate operation time 
	// Allocating memory for GPU (pointer to allocated device ); check for errors in the process; stops the program if issues encountered
	gpuErrchk( cudaMalloc((void**)&device_gEvent, NBytesAllEvent));
	gpuErrchk( cudaMalloc((void**)&device_output_TKL, NBytesAllSearchWindow));
	//allocating the memory for the planes
	gpuErrchk( cudaMalloc((void**)&device_gPlane, NBytesAllPlanes));
	gpuErrchk( cudaMalloc((void**)&device_gFitArrays, NBytesAllFitters));
	gpuErrchk( cudaMalloc((void**)&device_gFitParams, NBytesAllFitParam));
	//to give ourselves some leeway!
	//gpuErrchk( cudaMalloc((void**)&device_gPlane, NBytesAllPlanes));

	std::size_t free_bytes;
	std::size_t total_bytes;

	CUDA_CHECK_STATUS(cudaMemGetInfo(&free_bytes, &total_bytes));
    	cout << free_bytes << " / " << total_bytes << endl;
	
	// cudaMemcpy(dst, src, count, kind): copies data between host and device:
	// dst: destination memory address; src: source memory address; count: size in bytes; kind: type of transfer
	gpuErrchk( cudaMemcpy(device_gPlane, plane, NBytesAllPlanes, cudaMemcpyHostToDevice));
	gpuErrchk( cudaMemcpy(device_gFitParams, fitparams, NBytesAllFitParam, cudaMemcpyHostToDevice));
	gpuErrchk( cudaMemcpy(device_gEvent, host_gEvent, NBytesAllEvent, cudaMemcpyHostToDevice));
	clock_t cp4 = clock();
		
	// now data is transfered in the device: kernel function for event reconstruction called;
	// note that the function call is made requesting a number of blocks and a number of threads per block
	// in practice we have as many threads total as number of events; 
	auto start_er = std::chrono::system_clock::now();
	gkernel_eR<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent);
	auto end_er = std::chrono::system_clock::now();
	
	// check status of device and synchronize;
	size_t nEvents = EstnEvtMax;
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	// copy result of event reconstruction from device_gEvent to device_input_TKL
	// this input_tkl should be the information that the device uses to reconstruct the tracklets
	// gpuErrchk( cudaMemcpy(device_input_TKL, device_gEvent, NBytesAllEvent, cudaMemcpyDeviceToDevice));

	// shouldn't this function actually be called? should it be the function that puts together tracklets? and then call the fitting???
	// gkernel_TKL<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_input_TKL, device_output_TKL);

	//for(int m = 0; m<30; m++){
	//	if(plane[m].u_win!=0)printf("plane, m = %d, u_win = %1.6f, costheta = %1.6f\n", m, plane[m].u_win, plane[m].costheta);
	//	if(device_gPlane[m].u_win!=0)printf("device_gplane, m = %d, u_win = %1.6f, costheta = %1.6f\n", m, device_gPlane[m].u_win, device_gPlane[m].costheta);
	//}
	
	auto start_tkl2 = std::chrono::system_clock::now();
	// I first want to see if indeed we can reuse the "gEvent" pointer
	int stID = 3;// to make explicit that we are requiring station 3

	//lemme try something...
    	//dim3  threads(1, 1, 1);
    	//dim3  blocks(1, 1, 1);
	
	//threads.y = 5;
    	//threads.z = 4;
    	//threads.x = int(THREADS_PER_BLOCK/threads.y/threads.z);

	//blocks.x = int(EstnEvtMax/threads.x)+1;
	//blocks.y = 1;
	
	//gkernel_TrackletinStation<<<blocks, threads>>>(device_gEvent, device_output_TKL, device_gFitArrays, stID, device_gPlane, device_gFitParams);
	
	gkernel_TrackletinStation<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gFitArrays, stID, device_gPlane, device_gFitParams);
	auto end_tkl2 = std::chrono::system_clock::now();

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	auto start_tkl3 = std::chrono::system_clock::now();
	stID = 4;
	gkernel_TrackletinStation<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gFitArrays, stID, device_gPlane, device_gFitParams);

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	stID = 5;
	gkernel_TrackletinStation<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gFitArrays, stID, device_gPlane, device_gFitParams);
	//cout << endl;
	// check status of device and synchronize again;
	auto end_tkl3 = std::chrono::system_clock::now();
	
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	gkernel_BackPartialTracks<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gFitArrays, device_gPlane, device_gFitParams);
	
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	clock_t cp5 = clock();
	// data transfer from device to host
	gpuErrchk( cudaMemcpy(host_output_eR, device_gEvent, NBytesAllEvent, cudaMemcpyDeviceToHost));

	gpuErrchk( cudaMemcpy(host_output_TKL, device_output_TKL, NBytesAllSearchWindow, cudaMemcpyDeviceToHost));

	// thrust objects: C++ template library based on STL
	// convert raw pointer device_gEvent to device_vector
	// TODO: just don't use raw pointers to begin with
	//thrust::device_ptr<gEvent> d_p_events(device_gEvent);
	//thrust::device_vector<gEvent> d_events(d_p_events, d_p_events + nEvents);
    	
	//std::vector<gEvent> h_events(nEvents);
	//std::copy(d_events.begin(), d_events.end(), h_events.begin());

	//thrust::device_vector<float> d_hit_pos(nEvents);
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
	
	//TODO: need to check that host_gEvent is updated
	cudaMemcpy(host_gEvent, device_gEvent, NBytesAllEvent, cudaMemcpyDeviceToHost);
	cudaFree(device_gEvent);

	cudaMemcpy(host_output_TKL, device_output_TKL, NBytesAllSearchWindow, cudaMemcpyDeviceToHost);
	cudaFree(device_output_TKL);
	
	ofstream out("OutputFile.txt");
	//Write in a file, 
	for(int n = 0; n<nEvtMax; n++){
		if(host_gEvent[n].nAH==0)continue;
		out<<n<<" "<< host_gEvent[n].nAH <<" "<< host_output_TKL[n].nTracklets<<endl;
		for(int k = 1; k<=nDetectors; k++ ){
			out << host_gEvent[n].NHits[k] << " ";
		}out<<endl;
		
		for(int k = 0; k<host_gEvent[n].nAH; k++ ){
			out << host_gEvent[n].AllHits[k].detectorID << " " << host_gEvent[n].AllHits[k].elementID << " " << host_gEvent[n].AllHits[k].driftDistance << endl;
		}
		
		for(int k = 0; k<host_output_TKL[n].nTracklets; k++ ){
			if(isnan(host_output_TKL[n].AllTracklets[k].x0))host_output_TKL[n].AllTracklets[k].x0 = -200;
			if(isnan(host_output_TKL[n].AllTracklets[k].y0))host_output_TKL[n].AllTracklets[k].y0 = -100;
			if(isnan(host_output_TKL[n].AllTracklets[k].tx))host_output_TKL[n].AllTracklets[k].tx = -0.2;
			if(isnan(host_output_TKL[n].AllTracklets[k].ty))host_output_TKL[n].AllTracklets[k].ty = -0.2;
			out << host_output_TKL[n].AllTracklets[k].stationID << " " << host_output_TKL[n].AllTracklets[k].x0 << " " << host_output_TKL[n].AllTracklets[k].y0 << " " << host_output_TKL[n].AllTracklets[k].tx << " " << host_output_TKL[n].AllTracklets[k].ty << " " << endl;
		}
		
	}
	
	//auto end_kernel = std::chrono::system_clock::now();


	delete rawEvent;

	TFile* outFile = new TFile(outputFile.Data(), "RECREATE");
	/*
	//TTree* ORoutput_tree = new TTree("OR_out", "OR_out");
	ORoutput_tree* output = new ORoutput_tree();
	for(int i = 0; i < nEvtMax; ++i) {
		output->Clear();
		for(int k = 1; k<=nDetectors; k++ )output->fNhitsReduced[k] = host_output_eR[i].NHits[k];
		//output->Write();
	}
	output->Write();
	*/

	// printouts for test
	//for(int i = 0; i < nEvtMax; ++i) {
	//	//if(10000<host_gEvent[i].EventID && host_gEvent[i].EventID<10050){
	//	if(host_gEvent[i].EventID==0){
	//		printf("%d:\n ", host_gEvent[i].EventID);
	//		for(int j = 0; j<=nChamberPlanes; j++){
	//			printf("%d ", host_gEvent[i].NHits[j]);
	//		}printf("; %d\n", host_gEvent[i].nAH);
	//		for(int j = 0; j<=host_gEvent[i].nAH; j++){
	//			//if(13<=host_gEvent[i].AllHits[j].detectorID&&host_gEvent[i].AllHits[j].detectorID<=18)
	//			printf("%d, %1.3f; ", host_gEvent[i].AllHits[j].detectorID, host_gEvent[i].AllHits[j].pos);
	//		}
	//		printf("\n");
	//	}
	//}
		
	//for(int i = 0; i < host_gEvent[0].nAH; ++i) {
		//cout<<"D0_1st_wire:" << (host_gEvent[0].NHits[1])<<endl;
		//cout<<"output: "<<(host_gEvent[0].nAH)<<endl;
		//cout<<"output: "<<(device_output_eR)<<endl;
		//cout<<"output: "<<(sizeof(int))<<endl;
		//cout<<"size: "<<i<<endl;
	//}
	// printing the time required for all operations
	clock_t cp6 = clock();
	auto end = std::chrono::system_clock::now();

	auto evt_prep = cp2-cp1;
	auto cp_to_gpu = cp4-cp3;
	auto cp_from_gpu = cp6-cp5;
	
	auto gpu_er = end_er - start_er;
	auto er_tkl = start_tkl2-end_er;
	auto gpu_tkl2 = end_tkl2 - start_tkl2;
	auto tkl2_tkl3 = start_tkl3-end_tkl2;
	auto gpu_tkl3 = end_tkl3 - start_tkl3;
	auto overall = end - start;
	cout<<"Read/prepare events: "<<evt_prep/1000000000.<<endl;
	cout<<"Copy to GPU: "<<cp_to_gpu/1000000000.<<endl;
	cout<<"Copy from GPU: "<<cp_from_gpu/1000000000.<<endl;
	cout<<"event reducing: "<<(gpu_er.count())/1000000000.<<endl;
	cout<<" <-> : "<<(er_tkl.count())/1000000000.<<endl;
	cout<<"st2 trackletting: "<<(gpu_tkl2.count())/1000000000.<<endl;
	cout<<" <-> : "<<(tkl2_tkl3.count())/1000000000.<<endl;
	cout<<"st3 trackletting: "<<(gpu_tkl3.count())/1000000000.<<endl;
	cout<<"Total time: "<<(overall.count())/1000000000.<<endl;
		
	return 0;
}

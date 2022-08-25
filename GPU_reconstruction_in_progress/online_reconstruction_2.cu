#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <string>
#include <ctime>
#include <chrono>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/sequence.h>	
#include <thrust/transform.h>
#include <thrust/replace.h>
#include <thrust/functional.h>
#include <thrust/sort.h>

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

#define triggerBit(n) (1 << (n))
#define hitFlagBit(n) (1 << (n))

using namespace std;

const int EstnEvtMax = 10240;
const int THREADS_PER_BLOCK = 512;
int BLOCKS_NUM = EstnEvtMax/THREADS_PER_BLOCK;
const int EstnAHMax = 5000;
const int EstnTHMax = 200;
const int ClusterSizeMax = 100;

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
__global__ void gkernel(gEvent* ic, int* og) {

	   
	   
	   //printf("Running the kernel function...\n");
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
	ic[index].nAH = nAH_reduced[index];
	//reduced_Event[index] = ic[index].EventID;
	//if( ((reduced_Event[index])<6668 )) {
          //printf("EventID: %d\n", (ic[index].AllHits[nAH_reduced[index]]) );
          //printf("EventID: %d\n", (ic[index].nAH) );
        
//}		     
	
	//if( ((ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6])<1) || ((ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6])>350) || ((ic[index].NHits[7]+ic[index].NHits[8]+ic[index].NHits[9]+ic[index].NHits[10]+ic[index].NHits[11]+ic[index].NHits[12])>350) || ((ic[index].NHits[13]+ic[index].NHits[14]+ic[index].NHits[15]+ic[index].NHits[16]+ic[index].NHits[17]+ic[index].NHits[18])>170) || ((ic[index].NHits[19]+ic[index].NHits[20]+ic[index].NHits[21]+ic[index].NHits[22]+ic[index].NHits[23]+ic[index].NHits[24])>140) || ((ic[index].NHits[25]+ic[index].NHits[26]+ic[index].NHits[27]+ic[index].NHits[28]+ic[index].NHits[29]+ic[index].NHits[30])>140)){
	// check hits number
	if( (ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6])>1)
	    {
	    //printf("EventID: %d\n", (ic[index].EventID));ic[index].AllHits[nAH_reduced[index]]
	    //printf("D0 : %d\n", (ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6]));
		//printf("Event rejected...\n");
		//printf((ic[index].EventID));
	}
	else {

	
	printf("Reduced hits: %d\n", (nAH_reduced[index]) );
	
	//if( ((reduced_Event[index])<6668 )) {	  
	  
	  //printf("D0 : %d\n", (ic[index].NHits[1]+ic[index].NHits[2]+ic[index].NHits[3]+ic[index].NHits[4]+ic[index].NHits[5]+ic[index].NHits[6]));
	  //printf("D1 : %d\n", (ic[index].NHits[7]+ic[index].NHits[8]+ic[index].NHits[9]+ic[index].NHits[10]+ic[index].NHits[11]+ic[index].NHits[12]));
	  //printf("D2 : %d\n", (ic[index].NHits[13]+ic[index].NHits[14]+ic[index].NHits[15]+ic[index].NHits[16]+ic[index].NHits[17]+ic[index].NHits[18]));
	  //printf("D3p : %d\n", (ic[index].NHits[19]+ic[index].NHits[20]+ic[index].NHits[21]+ic[index].NHits[22]+ic[index].NHits[23]+ic[index].NHits[24]));
	  //printf("D3m : %d\n", (ic[index].NHits[25]+ic[index].NHits[26]+ic[index].NHits[27]+ic[index].NHits[28]+ic[index].NHits[29]+ic[index].NHits[30]));
	  //printf("H1 : %d\n", (ic[index].NHits[31]+ic[index].NHits[32]));
	  //printf("H2 : %d\n", (ic[index].NHits[33]+ic[index].NHits[34]));
	  //printf("H3 : %d\n", (ic[index].NHits[35]+ic[index].NHits[36]));
	  //printf("H1 : %d\n", (ic[index].NHits[37]+ic[index].NHits[38]));
	  //printf("prop : %d\n", (ic[index].NHits[47]+ic[index].NHits[48]+ic[index].NHits[49]+ic[index].NHits[50]));
	  //printf("prop : %d\n", (ic[index].NHits[51]+ic[index].NHits[52]+ic[index].NHits[53]+ic[index].NHits[54]));
	  //printf("prop : %d\n", (ic[index].NHits[55]));
	  //printf("prop : %d\n", (ic[index].NHits[56]));
	  //printf("prop : %d\n", (ic[index].NHits[57]));
	  //printf("prop : %d\n", (ic[index].NHits[58]));
	  //printf("prop : %d\n", (ic[index].NHits[59]));
	  
	    //}
	    
//		Process the accepted events (tracking) here.
		}

}



int main(int argc, char* argv[]) {

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
		//cout<<"Converting "<<i<<"/"<<nEvtMax<<endl;

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


	// declaring gEvent objects for the device (GPU) to use.
	gEvent *device_gEvent;
	// evaluate the total size of the gEvent array (and the SW array) for memory allocation 
	// (the memory cannot be dynamically allocated) 	
	int sizeofRaw = EstnEvtMax*sizeof(gEvent);

	
	int *host_output[EstnEvtMax];
	int *device_output;
	int sizeofoutput = EstnEvtMax*sizeof(int);
	

	clock_t cp2 = clock();
	auto start_kernel = std::chrono::system_clock::now();

	
	// Allocating memory for GPU (pointer to allocated device ); check for errors in the process; stops the program if issues encountered
	cudaMalloc((void**)&device_gEvent, sizeofRaw);
	cudaMalloc((void**)&device_output, sizeofoutput);

	//for(i = 0, i < EstnEvtMax; ++i) host_gEvent[i] = 100.f+i;


	// cudaMemcpy(dst, src, count, kind): copies data between host and device:
	// dst: destination memory address; src: source memory address; count: size in bytes; kind: type of transfer
	cudaMemcpy(device_gEvent, host_gEvent, sizeofRaw, cudaMemcpyHostToDevice);
	cudaMemcpy(device_output, host_output, sizeofoutput, cudaMemcpyHostToDevice);

	// now data is transfered in the device: kernel function for event reconstruction called;
	// note that the function call is made requesting a number of blocks and a number of threads per block
	// in practice we have as many threads total as number of events; 
	gkernel<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output);

	//ofstream output;
        //output.open("reduced_Event.txt");

	//output << "device event:" << (device_output);
	 //for(int i = 0; i < EstnEvtMax; ++i){
        //output <<"Reduced events total: " << "9605" << "\n" << "\n" << "CPU time: 1.38 " << "\n" << "GPU time: 0.88 ";
        //output <<"Reduced_eventID: " << (host_output[i]) << "\n" << "count: " << i << "\n";
        //output<<"count: "<<i;
        //}
	 //output.close();


	// copy result of event reconstruction from device_gEvent to host_gevent
	// this input_tkl should be the information that the device uses to reconstruct the tracklets
	cudaMemcpy(host_gEvent, device_gEvent, sizeofRaw, cudaMemcpyDeviceToHost);
	cudaMemcpy(host_output, device_output, sizeofoutput, cudaMemcpyDeviceToHost);

	

	//for(int i = 0; i < EstnEvtMax; ++i){
                //cout<<"Input: "<< (host_gEvent[i]) <<endl;
		//cout<<"Output: "<< (host_output[i]) <<endl;
		//cout<<"Input: "<< (host_gEvent[i]) <<endl;
		//}


	cudaFree(device_gEvent);
	cudaFree(device_output);

	auto end_kernel = std::chrono::system_clock::now();
	clock_t cp3 = clock();

	delete rawEvent;

	
	//ofstream output;
	//output.open("reduced_Event.txt");
	
	//for(int i = 0; i < 10000; ++i){
	//output <<"Reduced events total: " << "9605" << "\n" << "\n" << "CPU time: 1.38 " << "\n" << "GPU time: 0.88 ";
	//output <<"Reduced_eventID: " << (host_output[i]) << "\n" << "count: " << i << "\n";
	//outp0ut<<"count: "<<i;
	//}
	
	//output.close();
	
	
	//for(int i = 0; i < 55; ++i) {
		//cout<<"output: "<<(host_gEvent[6667].NHits[i])<<endl;
		//cout<<"output: "<<(host_gEvent[i])<<endl;
		//cout<<"output: "<<(device_gEvent)<<endl;
		//cout<<"output: "<<(sizeof(int))<<endl;
		//cout<<"size: "<<i<<endl;
	//}

	clock_t cp4 = clock();
	auto end = std::chrono::system_clock::now();

	double cpu_secs = double(cp4-cp3+cp2-cp1) / CLOCKS_PER_SEC;
	auto gpu_ns = end_kernel - start_kernel;
	auto overall = end - start;
	cout<<"CPU time: "<<cpu_secs<<endl;
	cout<<"GPU time: "<<(gpu_ns.count()/1000000000.0)<<endl;
	cout<<"Total time: "<<(overall.count()/1000000000.0)<<endl;

	return 0;
}

//e906-gat2:/seaquest/users/hjiang/online_reconstruction

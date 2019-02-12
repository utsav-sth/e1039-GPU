#include <cstdio>
#include <cstdlib>
#include <iostream>
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

class gHit {
	public:
	int index;
	short detectorID;
	short elementID;
	float tdcTime;
	float driftDistance;
	float pos;
	short flag;
};

class gEvent {
	public:
	int RunID;
	int EventID;
	int SpillID;
	int TriggerBits;
	short TargetPos;
	int TurnID;
	int RFID;
	int Intensity[33];
	short TriggerEmu;
	short NRoads[4];
	int NHits[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
	int nAH;
	int nTH;
	gHit AllHits[EstnAHMax];
	gHit TriggerHits[EstnTHMax];
};

struct lessthan {
	__host__ __device__ bool operator()(const gHit& lhs, const gHit& rhs)
	{
		if(lhs.detectorID < rhs.detectorID)
		{
			return true;
		}
		else if(lhs.detectorID > rhs.detectorID)
		{
			return false;
		}

		if(lhs.elementID < rhs.elementID)
		{
			return true;
		}
		else if(lhs.elementID > rhs.elementID)
		{
			return false;
		}

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

__global__ void gkernel(gEvent* ic, int* og) {
//	printf("Running the kernel function...\n");
	int index = threadIdx.x + blockIdx.x * blockDim.x;

	double w_max[EstnEvtMax];
	double w_min[EstnEvtMax];
	double dt_mean[EstnEvtMax];
	int cluster_iAH_arr_cur[EstnEvtMax];
	int cluster_iAH_arr_size[EstnEvtMax];
	static int cluster_iAH_arr[EstnEvtMax][ClusterSizeMax];
	int uniqueID[EstnEvtMax];
	int uniqueID_curr[EstnEvtMax];
	double tdcTime_curr[EstnEvtMax];
	int iAH[EstnEvtMax];
	int nAH_reduced[EstnEvtMax];

	cluster_iAH_arr_size[index] = 0;
	nAH_reduced[index] = 0;
	for(iAH[index] = 0; iAH[index]<ic[index].nAH; ++iAH[index]) {
		if((ic[index].AllHits[iAH[index]].flag & hitFlagBit(1)) == 0) {
//			printf("Skip out-of-time...\n");
			ic[index].AllHits[iAH[index]].detectorID = 0;
			continue;
		}
		if(ic[index].AllHits[iAH[index]].detectorID < 31 || ic[index].AllHits[iAH[index]].detectorID > 46) {
			uniqueID[index] = ic[index].AllHits[iAH[index]].detectorID*1000 + ic[index].AllHits[iAH[index]].elementID;
			if(uniqueID[index] != uniqueID_curr[index]) {
				uniqueID_curr[index] = uniqueID[index];
				tdcTime_curr[index] = ic[index].AllHits[iAH[index]].tdcTime;
			}
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
		if(ic[index].AllHits[iAH[index]].detectorID <= nChamberPlanes) {
//			printf("%d\n", cluster_iAH_arr_size[index]);
//			printf("Decluster...\n");
			if(cluster_iAH_arr_size[index] == ClusterSizeMax) {
//				printf("Oversized cluster...\n");
			}
			if(cluster_iAH_arr_size[index] == 0) {
				cluster_iAH_arr[index][0] = iAH[index];
				++cluster_iAH_arr_size[index];
			}
			else {
				if((ic[index].AllHits[iAH[index]].detectorID != ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID) || (ic[index].AllHits[iAH[index]].elementID - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].elementID > 1)) {
					if(cluster_iAH_arr_size[index] == 2) {
						w_max[index] = 0.9*0.5*(ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].pos - ic[index].AllHits[cluster_iAH_arr[index][0]].pos);
						w_min[index] = 4.0/9.0*w_max[index];
						if((ic[index].AllHits[cluster_iAH_arr[index][0]].driftDistance > w_max[index] && ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].driftDistance > w_min[index]) || (ic[index].AllHits[cluster_iAH_arr[index][0]].driftDistance > w_min[index] && ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].driftDistance > w_max[index])) {
							if(ic[index].AllHits[cluster_iAH_arr[index][0]].driftDistance > ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].driftDistance) {
//								printf("Skip cluster...\n");
								ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID = 0;
							}
							else {
//								printf("Skip cluster...\n");
								ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID = 0;
							}
						}
						else if((((ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) >= 0.0 && (ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) < 8.0) || ((ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) <= 0.0 && (ic[index].AllHits[cluster_iAH_arr[index][0]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].tdcTime) > -8.0)) && (ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID >= 19 && ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID <= 24)) {
//							printf("Skip cluster...\n");
							ic[index].AllHits[cluster_iAH_arr[index][0]].detectorID = 0;
							ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_size[index]-1]].detectorID = 0;
						}
					}
					if(cluster_iAH_arr_size[index] >= 3) {
						dt_mean[index] = 0.0;
						for(cluster_iAH_arr_cur[index] = 1; cluster_iAH_arr_cur[index] < cluster_iAH_arr_size[index]; ++cluster_iAH_arr_cur[index]) {
							dt_mean[index] += ((ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]-1]].tdcTime) > 0.0 ? (ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]-1]].tdcTime) : (ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]-1]].tdcTime - ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].tdcTime));
						}
						dt_mean[index] = dt_mean[index]/(cluster_iAH_arr_size[index] - 1);
						if(dt_mean[index] < 10.0) {
//							printf("Skip cluster...\n");
							for(cluster_iAH_arr_cur[index] = 0; cluster_iAH_arr_cur[index] < cluster_iAH_arr_size[index]; ++cluster_iAH_arr_cur[index]) {
								ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].detectorID = 0;
							}
						}
						else {
//							printf("Skip cluster...\n");
							for(cluster_iAH_arr_cur[index] = 1; cluster_iAH_arr_cur[index] < cluster_iAH_arr_size[index]; ++cluster_iAH_arr_cur[index]) {
								ic[index].AllHits[cluster_iAH_arr[index][cluster_iAH_arr_cur[index]]].detectorID = 0;
							}
						}
					}
					cluster_iAH_arr_size[index] = 0;
				}
				cluster_iAH_arr[index][cluster_iAH_arr_size[index]] = iAH[index];
				++cluster_iAH_arr_size[index];
			}
		}
	}

	for(iAH[index] = 0; iAH[index]<ic[index].nAH; ++iAH[index]) {
		if(ic[index].AllHits[iAH[index]].detectorID != 0) {
			ic[index].AllHits[nAH_reduced[index]] = ic[index].AllHits[iAH[index]];
			++nAH_reduced[index];
		}
	}
	ic[index].nAH = nAH_reduced[index];

}

int main(int argc, char* argv[]) {

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

	gEvent *device_gEvent;
	int sizeofRaw = EstnEvtMax*sizeof(gEvent);

	int host_output[EstnEvtMax];
	int *device_output;
	int sizeofoutput = EstnEvtMax*sizeof(int);

	clock_t cp2 = clock();
	auto start_kernel = std::chrono::system_clock::now();

	cudaMalloc((void**)&device_gEvent, sizeofRaw);
	cudaMalloc((void**)&device_output, sizeofoutput);
	cudaMemcpy(device_gEvent, host_gEvent, sizeofRaw, cudaMemcpyHostToDevice);
	cudaMemcpy(device_output, host_output, sizeofoutput, cudaMemcpyHostToDevice);

	gkernel<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output);

	cudaMemcpy(host_gEvent, device_gEvent, sizeofRaw, cudaMemcpyDeviceToHost);
	cudaMemcpy(host_output, device_output, sizeofoutput, cudaMemcpyDeviceToHost);
	cudaFree(device_gEvent);
	cudaFree(device_output);

	auto end_kernel = std::chrono::system_clock::now();
	clock_t cp3 = clock();

	delete rawEvent;

	for(int i = 0; i < host_gEvent[0].nAH; ++i) {
		cout<<"output: "<<(host_gEvent[0].AllHits[i].detectorID)<<endl;
	}

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

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
#include "reconstruction_kernels.cuh"

#include "SQEvent_v1.h"
#include "SQHit_v1.h"
#include "SQHitVector_v1.h"


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
	gPlane plane[nChamberPlanes+nHodoPlanes+nPropPlanes+1];

	ifstream in_geom(inputGeom.Data());
  	string buffer;
	int ipl, nelem;
	float z, cellwidth, spacing, xoffset, scalex, x0, x1, x2, costheta, scaley, y0, y1, y2, sintheta, resolution, deltaW_;
	float p1x, p1y, p1z, deltapx, deltapy, deltapz, dp1x, dp1y, dp1z;
 	while ( getline(in_geom, buffer) ) {
    	      if (buffer[0] == '#') continue;
	      std::istringstream iss;
	      iss.str(buffer);
	      iss >> ipl >> z >> nelem >> cellwidth >> spacing >> xoffset >> scalex >> x0 >> x1 >> x2 >> costheta >> scaley >> y0 >> y1 >> y2 >> sintheta >> resolution >> p1x >> p1y >> p1z >> deltapx >> deltapy >> deltapz >> dp1x >> dp1y >> dp1z;
	      plane[ipl].z = z;
	      plane[ipl].nelem = nelem;
	      plane[ipl].cellwidth = cellwidth;
	      plane[ipl].spacing = spacing;
	      plane[ipl].xoffset = xoffset;
	      plane[ipl].scalex = scalex;
	      plane[ipl].x0 = x0;
	      plane[ipl].x1 = x1;
	      plane[ipl].x1 = x2;
	      plane[ipl].costheta = costheta;
	      plane[ipl].scaley = scaley;
	      plane[ipl].y0 = y0;
	      plane[ipl].y1 = y1;
	      plane[ipl].y2 = y2;
	      plane[ipl].sintheta = sintheta;
	      plane[ipl].resolution = resolution;
	      plane[ipl].p1x_w1 = p1x;
	      plane[ipl].p1y_w1 = p1y;
	      plane[ipl].p1z_w1 = p1z;
	      plane[ipl].deltapx = deltapx;
	      plane[ipl].deltapy = deltapy;
	      plane[ipl].deltapz = deltapz;
	      plane[ipl].dp1x = dp1x;
	      plane[ipl].dp1y = dp1y;
	      plane[ipl].dp1z = dp1z;
	      if(ipl>nChamberPlanes+nHodoPlanes){
		for(int k = 0; k<9; k++){
			iss >> deltaW_;
			plane[ipl].deltaW_[k] = deltaW_;
		}
	      }else{
		iss >> deltaW_;
		plane[ipl].deltaW_[0] = deltaW_;
	      }
	      ipl++;
	}
	
	for(int i = 0; i<5; i++){
		int u_idx = i*6+5;
		if(i==0)u_idx = i*6+1;
		int x_idx = i*6+3;
		for(int j = 0; j<6; j++){
			int idx = i*6+j+1;
			plane[idx].z_mean = j%2==0 ? 0.5*(plane[idx].z+plane[idx+1].z):0.5*(plane[idx].z+plane[idx-1].z);
			
			plane[idx].v_win_fac1 = plane[idx].spacing*2*plane[u_idx].costheta;
			plane[idx].v_win_fac2 = plane[u_idx].costheta*TX_MAX;
			plane[idx].v_win_fac3 = plane[u_idx].sintheta*TY_MAX;
		}
		
		for(int j = 0; j<6; j++){
			int idx = i*6+j+1;
			plane[idx].u_win = fabs(0.5*plane[u_idx].scaley*plane[u_idx].sintheta) + TX_MAX*fabs((plane[u_idx].z_mean - plane[x_idx].z_mean)*plane[u_idx].costheta) + TY_MAX*fabs((plane[u_idx].z_mean - plane[x_idx].z_mean)*plane[u_idx].sintheta) + 2.*plane[u_idx].spacing + u_factor[i];
		}
		cout << u_idx << " " << plane[u_idx].u_win << " = " << fabs(0.5*plane[u_idx].scaley*plane[u_idx].sintheta) << " + " << TX_MAX*fabs((plane[u_idx].z_mean - plane[x_idx].z_mean)*plane[u_idx].costheta) << " + " << TY_MAX*fabs((plane[u_idx].z_mean - plane[x_idx].z_mean)*plane[u_idx].sintheta) << " + " << 2.*plane[u_idx].spacing + u_factor[i] << endl;
		cout << " u costheta " << plane[u_idx].costheta << " u sintheta " << plane[u_idx].sintheta << " x_span " << plane[u_idx].scaley << " spacing " << plane[u_idx].spacing << " z plane_u " << plane[u_idx].z_mean << " z plane_x " << plane[x_idx].z_mean << endl;  
	}
	cout << "Geometry file read out" << endl;
	
	double wire_position[55][400];//Let's keep this: simpler, more robust
	for(int i = 1; i <= nChamberPlanes; ++i){
		//cout << plane[i].nelem << endl;
      		for(int j = 1; j <= plane[i].nelem; ++j){
          		double pos = (j - (plane[i].nelem+1.)/2.)*plane[i].spacing + plane[i].xoffset + plane[i].x0*plane[i].costheta + plane[i].y0*plane[i].sintheta + plane[i].deltaW_[0];
			wire_position[i][j] = pos;
			//if(25<=i && i<=30){
			//	 double p1x_ = plane[i].p1x_w1+plane[i].dp1x*(j-1);
			//	 double p2x_ = p1x_+plane[i].deltapx;
			//	 cout << i << " " << j << " " << pos << " " << p1x_ << " " << p2x_ << " " << (p1x_+p2x_)/2. << " " << (p1x_+p2x_)/2.-pos << endl;
			//}
		}
	}
	for(int i = nChamberPlanes+1; i<=nChamberPlanes+nHodoPlanes; ++i){
		//cout << plane[i].nelem << endl;
	      	for(int j = 1; j <= plane[i].nelem; ++j){
          		double pos = plane[i].x0*plane[i].costheta + plane[i].y0*plane[i].sintheta + plane[i].xoffset + (j - (plane[i].nelem+1)/2.)*plane[i].spacing + plane[i].deltaW_[0];
			wire_position[i][j] = pos;
		}
	}
	for(int i = nChamberPlanes+nHodoPlanes+1; i<=nChamberPlanes+nHodoPlanes+nPropPlanes; ++i){
		//cout << plane[i].nelem << endl;
	      	for(int j = 1; j <= plane[i].nelem; ++j){
          		int moduleID = 8 - int((j - 1)/8);
			//cout << moduleID << endl;
             		double pos = plane[i].x0*plane[i].costheta + plane[i].y0*plane[i].sintheta + plane[i].xoffset + (j - (plane[i].nelem+1)/2.)*plane[i].spacing + plane[i].deltaW_[moduleID];
			wire_position[i][j] = pos;
		}
		
	}
	/*
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
	*/
		
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
				host_gEvent[i].AllHits[m].pos=wire_position[(rawEvent->fAllHits[m]).detectorID][(rawEvent->fAllHits[m]).elementID];
				host_gEvent[i].AllHits[m].flag=(rawEvent->fAllHits[m]).flag;
			}
			for(int n=0; n<rawEvent->fTriggerHits.size(); n++) {
				host_gEvent[i].TriggerHits[n].index=(rawEvent->fTriggerHits[n]).index;
				host_gEvent[i].TriggerHits[n].detectorID=(rawEvent->fTriggerHits[n]).detectorID;
				host_gEvent[i].TriggerHits[n].elementID=(rawEvent->fTriggerHits[n]).elementID;
				host_gEvent[i].TriggerHits[n].tdcTime=(rawEvent->fTriggerHits[n]).tdcTime;
				host_gEvent[i].TriggerHits[n].driftDistance=(rawEvent->fTriggerHits[n]).driftDistance;
				host_gEvent[i].TriggerHits[n].pos=wire_position[(rawEvent->fAllHits[n]).detectorID][(rawEvent->fAllHits[n]).elementID];
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
			host_gEvent[i].nAH = 0;//hit_vec.size();
			for(int m = 0; m<hit_vec.size(); m++){
				if(hit_vec[m]->get_detector_id()>54){
					//if(_event_id<20)cout << " dark photon plane hit! " << hit_vec[m]->get_detector_id() << endl;
					continue;
					//dark photon planes; I don't think we care about those for the purpose of online reconstruction... do we?
				}
				host_gEvent[i].nAH++;
				host_gEvent[i].NHits[hit_vec[m]->get_detector_id()]++;
				host_gEvent[i].AllHits[m].index=hit_vec[m]->get_hit_id();
				host_gEvent[i].AllHits[m].detectorID=hit_vec[m]->get_detector_id();
				host_gEvent[i].AllHits[m].elementID=hit_vec[m]->get_element_id();
				host_gEvent[i].AllHits[m].tdcTime=hit_vec[m]->get_tdc_time();
				host_gEvent[i].AllHits[m].driftDistance=fabs(hit_vec[m]->get_drift_distance());
				host_gEvent[i].AllHits[m].pos=wire_position[hit_vec[m]->get_detector_id()][hit_vec[m]->get_element_id()];
				host_gEvent[i].AllHits[m].flag=(1<<hit_vec[m]->is_in_time());
				//if(host_gEvent[i].EventID<20)cout << " det " << host_gEvent[i].AllHits[m].detectorID << " elem " << host_gEvent[i].AllHits[m].elementID << " time " << host_gEvent[i].AllHits[m].tdcTime << " dd " << host_gEvent[i].AllHits[m].driftDistance << " pos " << host_gEvent[i].AllHits[m].pos << endl;
				if(hit_vec[m]->is_trigger_mask()){
					host_gEvent[i].TriggerHits[ntrighits].index=hit_vec[m]->get_hit_id();
					host_gEvent[i].TriggerHits[ntrighits].detectorID=hit_vec[m]->get_detector_id();
					host_gEvent[i].TriggerHits[ntrighits].elementID=hit_vec[m]->get_element_id();
					host_gEvent[i].TriggerHits[ntrighits].tdcTime=hit_vec[m]->get_tdc_time();
					host_gEvent[i].TriggerHits[ntrighits].driftDistance=fabs(hit_vec[m]->get_drift_distance());
					host_gEvent[i].TriggerHits[ntrighits].pos=wire_position[hit_vec[m]->get_detector_id()][hit_vec[m]->get_element_id()];
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
	//size_t NBytesAllSearchWindow = EstnEvtMax * sizeof(gSW);
	size_t NBytesAllOutputEvent = EstnEvtMax * sizeof(gOutputEvent);
	size_t NBytesAllPlanes =  nDetectors * sizeof(gPlane);
	//size_t NBytesAllFitters = EstnEvtMax * sizeof(gFitArrays);
	//size_t NBytesAllFitParam = sizeof(gFitParams)*3;// just to be sure: tracklets, back tracks, global tracks
	//size_t NBytesFitterTools = EstnEvtMax * sizeof(gStraightFitArrays);
	size_t NBytesFitterTools = EstnEvtMax * sizeof(gStraightFitArrays);
	size_t NBytesStraightTrackBuilders = EstnEvtMax * sizeof(gStraightTrackBuilder);
	size_t NBytesFullTrackBuilders = EstnEvtMax * sizeof(gFullTrackBuilder);
	

	//cout << NBytesAllEvent << " " << NBytesAllSearchWindow << " " << NBytesAllPlanes << " " << NBytesAllFitters << " " << NBytesAllFitParam << endl;
	//cout << NBytesAllEvent+NBytesAllSearchWindow+NBytesAllPlanes+NBytesAllFitters+NBytesAllFitParam+NBytesStraightTrackBuilders << endl;
	
	cout << "Total size allocated on GPUs " << NBytesAllEvent+NBytesAllOutputEvent+NBytesAllPlanes+NBytesFitterTools << endl;
	cout << " input events: " << NBytesAllEvent << "; output events: " << NBytesAllOutputEvent << "; straight track builder tools: " << NBytesStraightTrackBuilders
	     << "; fitter tools: " << NBytesFitterTools << "; planes info: " << NBytesAllPlanes << endl;  
	
	gEvent *host_output_eR = (gEvent*)malloc(NBytesAllEvent);
	//gSW *host_output_TKL = (gSW*)malloc(NBytesAllSearchWindow);
	gOutputEvent *host_output_TKL = (gOutputEvent*)malloc(NBytesAllOutputEvent);
	
	// declaring gEvent objects for the device (GPU) to use.
	gEvent *device_gEvent;
	// gEvent *device_output_eR;
	// gEvent *device_input_TKL;
	//gSW *device_output_TKL;
	gOutputEvent *device_output_TKL;
	gPlane *device_gPlane;
	//gFitParams *device_gFitParams;
	//gFitArrays *device_gFitArrays;
	gStraightFitArrays *device_gFitArrays;
	gStraightTrackBuilder *device_gStraightTrackBuilder;
	gFullTrackBuilder *device_gFullTrackBuilder;
	//gFullFitArrays *device_gKalmanFitArrays;

	//printDeviceStatus();

	// copy of data from host to device: evaluate operation time 
	// Allocating memory for GPU (pointer to allocated device ); check for errors in the process; stops the program if issues encountered
	gpuErrchk( cudaMalloc((void**)&device_gEvent, NBytesAllEvent));
	//gpuErrchk( cudaMalloc((void**)&device_output_TKL, NBytesAllSearchWindow));
	gpuErrchk( cudaMalloc((void**)&device_output_TKL, NBytesAllOutputEvent));
	//allocating the memory for the planes
	gpuErrchk( cudaMalloc((void**)&device_gPlane, NBytesAllPlanes));
	//gpuErrchk( cudaMalloc((void**)&device_gFitArrays, NBytesAllFitters));
	//gpuErrchk( cudaMalloc((void**)&device_gFitParams, NBytesAllFitParam));
	gpuErrchk( cudaMalloc((void**)&device_gFitArrays, NBytesFitterTools));
	gpuErrchk( cudaMalloc((void**)&device_gStraightTrackBuilder, NBytesStraightTrackBuilders));

	std::size_t free_bytes;
	std::size_t total_bytes;

	CUDA_CHECK_STATUS(cudaMemGetInfo(&free_bytes, &total_bytes));
    	cout << free_bytes << " / " << total_bytes << endl;
	
	// cudaMemcpy(dst, src, count, kind): copies data between host and device:
	// dst: destination memory address; src: source memory address; count: size in bytes; kind: type of transfer
	gpuErrchk( cudaMemcpy(device_gPlane, plane, NBytesAllPlanes, cudaMemcpyHostToDevice));
	//gpuErrchk( cudaMemcpy(device_gFitParams, fitparams, NBytesAllFitParam, cudaMemcpyHostToDevice));
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
	
	auto start_straight = std::chrono::system_clock::now();
	
	gKernel_XZ_YZ_tracking<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gStraightTrackBuilder, device_gFitArrays, device_gPlane);
	
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	auto end_straight = std::chrono::system_clock::now();
	
	//release here the memory for straight track builders and straight track fitters	
	cudaFree(device_gStraightTrackBuilder);
	
	gpuErrchk( cudaMalloc((void**)&device_gFullTrackBuilder, NBytesFullTrackBuilders));
	
	auto start_global = std::chrono::system_clock::now();

	gKernel_GlobalTracking<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gFullTrackBuilder, device_gFitArrays, device_gPlane, 0);

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	auto end_global = std::chrono::system_clock::now();






#ifdef KTRACKER_REC	
	// copy result of event reconstruction from device_gEvent to device_input_TKL
	// this input_tkl should be the information that the device uses to reconstruct the tracklets
	// gpuErrchk( cudaMemcpy(device_input_TKL, device_gEvent, NBytesAllEvent, cudaMemcpyDeviceToDevice));

	// shouldn't this function actually be called? should it be the function that puts together tracklets? and then call the fitting???
	// gkernel_TKL<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_input_TKL, device_output_TKL);

	//for(int m = 1; m<=30; m++){
	//	if(plane[m].u_win!=0)printf("plane, m = %d, u_win = %1.6f, costheta = %1.6f\n", m, plane[m].u_win, plane[m].costheta);
	//	if(device_gPlane[m].u_win!=0)printf("device_gplane, m = %d, u_win = %1.6f, costheta = %1.6f\n", m, device_gPlane[m].u_win, device_gPlane[m].costheta);
	//}
	
	auto start_tkl2 = std::chrono::system_clock::now();
	// I first want to see if indeed we can reuse the "gEvent" pointer
	int stID = 3;// to make explicit that we are requiring station 3
	
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
#endif
	
	
	
	clock_t cp5 = clock();
	// data transfer from device to host
	gpuErrchk( cudaMemcpy(host_output_eR, device_gEvent, NBytesAllEvent, cudaMemcpyDeviceToHost));
	cudaFree(device_gEvent);
	
	//gpuErrchk( cudaMemcpy(host_output_TKL, device_output_TKL, NBytesAllSearchWindow, cudaMemcpyDeviceToHost));
	gpuErrchk( cudaMemcpy(host_output_TKL, device_output_TKL, NBytesAllOutputEvent, cudaMemcpyDeviceToHost));
	cudaFree(device_output_TKL);

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
	//cudaMemcpy(host_gEvent, device_gEvent, NBytesAllEvent, cudaMemcpyDeviceToHost);

	//cudaMemcpy(host_output_TKL, device_output_TKL, NBytesAllSearchWindow, cudaMemcpyDeviceToHost);
	//cudaMemcpy(host_output_TKL, device_output_TKL, NBytesAllOutputEvent, cudaMemcpyDeviceToHost);
	
	ofstream out("OutputFile.txt");
	//Write in a file, 
	long tklctr = 0;
	for(int n = 0; n<nEvtMax; n++){
		if(host_output_eR[n].nAH==0)continue;
		out<<n<<" "<< host_output_eR[n].nAH <<" "<< host_output_TKL[n].nTracklets<<endl;
		tklctr+= host_output_TKL[n].nTracklets;
		for(int k = 1; k<=nDetectors; k++ ){
			out << host_output_eR[n].NHits[k] << " ";
		}out<<endl;
		
		for(int k = 0; k<host_output_eR[n].nAH; k++ ){
			out << host_output_eR[n].AllHits[k].detectorID << " " << host_output_eR[n].AllHits[k].elementID << " " << host_output_eR[n].AllHits[k].driftDistance << endl;
		}
		
		for(int k = 0; k<host_output_TKL[n].nTracklets; k++ ){
			if(isnan(host_output_TKL[n].AllTracklets[k].x0))host_output_TKL[n].AllTracklets[k].x0 = -200;
			if(isnan(host_output_TKL[n].AllTracklets[k].y0))host_output_TKL[n].AllTracklets[k].y0 = -100;
			if(isnan(host_output_TKL[n].AllTracklets[k].tx))host_output_TKL[n].AllTracklets[k].tx = -0.2;
			if(isnan(host_output_TKL[n].AllTracklets[k].ty))host_output_TKL[n].AllTracklets[k].ty = -0.2;
			out << host_output_TKL[n].AllTracklets[k].stationID << " " << host_output_TKL[n].AllTracklets[k].x0 << " " << host_output_TKL[n].AllTracklets[k].y0 << " " << host_output_TKL[n].AllTracklets[k].tx << " " << host_output_TKL[n].AllTracklets[k].ty << " " << host_output_TKL[n].AllTracklets[k].invP << endl;
		}
		
	}
	
	cout << tklctr << " straight tracks reconstructed" << endl;
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
	//auto er_tkl = start_tkl2-end_er;
	//auto gpu_tkl2 = end_tkl2 - start_tkl2;
	//auto tkl2_tkl3 = start_tkl3-end_tkl2;
	//auto gpu_tkl3 = end_tkl3 - start_tkl3;
	auto overall = end - start;
	cout<<"Read/prepare events: "<<evt_prep/1000000000.<<endl;
	cout<<"Copy to GPU: "<<cp_to_gpu/1000000000.<<endl;
	cout<<"Copy from GPU: "<<cp_from_gpu/1000000000.<<endl;
	cout<<"event reducing: "<<(gpu_er.count())/1000000000.<<endl;
	//cout<<" <-> : "<<(er_tkl.count())/1000000000.<<endl;
	//cout<<"st2 trackletting: "<<(gpu_tkl2.count())/1000000000.<<endl;
	//cout<<" <-> : "<<(tkl2_tkl3.count())/1000000000.<<endl;
	//cout<<"st3 trackletting: "<<(gpu_tkl3.count())/1000000000.<<endl;
	cout<<"Total time: "<<(overall.count())/1000000000.<<endl;
		
	return 0;
}

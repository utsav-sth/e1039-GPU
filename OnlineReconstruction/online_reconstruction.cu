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

#ifdef E1039
#include "SQEvent_v1.h"
#include "SQHit_v1.h"
#include "SQHitVector_v1.h"
#endif

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
	auto cp1 = std::chrono::system_clock::now();

	TString inputFile;
	TString inputGeom;
	TString outputFile;
	inputFile = argv[1];
	inputGeom = argv[2];	
	outputFile = argv[3];

	//by default we should use e1039 
	bool e906data = true;
	//if(argn>4)e906data = atoi(argv[4]);
	
	cout<<"Running "<<argv[0]<<endl;
	cout<<"Loading "<<argv[1]<<endl;
	cout<<"with geometry: "<<argv[2]<<endl;
	cout<<"Writing "<<argv[3]<<endl;
	
	//Get basic geometry here:
	double u_factor[5] = {5., 5., 5., 15., 15.};
	gPlane plane;
	
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
	      plane.z[ipl] = z;
	      plane.nelem[ipl] = nelem;
	      plane.cellwidth[ipl] = cellwidth;
	      plane.spacing[ipl] = spacing;
	      plane.xoffset[ipl] = xoffset;
	      plane.scalex[ipl] = scalex;
	      plane.x0[ipl] = x0;
	      plane.x1[ipl] = x1;
	      plane.x2[ipl] = x2;
	      plane.costheta[ipl] = costheta;
	      plane.scaley[ipl] = scaley;
	      plane.y0[ipl] = y0;
	      plane.y1[ipl] = y1;
	      plane.y2[ipl] = y2;
	      plane.sintheta[ipl] = sintheta;
	      plane.resolution[ipl] = resolution;
	      plane.p1x_w1[ipl] = p1x;
	      plane.p1y_w1[ipl] = p1y;
	      plane.p1z_w1[ipl] = p1z;
	      plane.deltapx[ipl] = deltapx;
	      plane.deltapy[ipl] = deltapy;
	      plane.deltapz[ipl] = deltapz;
	      plane.dp1x[ipl] = dp1x;
	      plane.dp1y[ipl] = dp1y;
	      plane.dp1z[ipl] = dp1z;
	      if(ipl>nChamberPlanes+nHodoPlanes){
		for(int k = 0; k<9; k++){
			iss >> deltaW_;
			plane.deltaW_[ipl*9+k] = deltaW_;
		}
	      }else{
		iss >> deltaW_;
		plane.deltaW_[ipl*9] = deltaW_;
	      }
	      plane.slope_max[ipl] = costheta*TX_MAX+sintheta*TY_MAX;
	      plane.inter_max[ipl] = costheta*X0_MAX+sintheta*Y0_MAX;
	      if(ipl%2==0 && ipl>1){
		double dslope = (plane.resolution[ipl]+plane.resolution[ipl-1])/(plane.z[ipl]-plane.z[ipl-1]);
		double dinter = dslope*plane.z[ipl];
		plane.slope_max[ipl]+= dslope;
		plane.inter_max[ipl]+= dinter;
		plane.slope_max[ipl-1]+= dslope;
		plane.inter_max[ipl-1]+= dinter;
	      }
	      ipl++;
	}
	
	for(int i = 0; i<5; i++){
		int u_idx = i*6+5;
		if(i==0)u_idx = i*6+1;
		int x_idx = i*6+3;
		for(int j = 0; j<6; j++){
			int idx = i*6+j+1;
			plane.z_mean[idx] = j%2==0 ? 0.5*(plane.z[idx]+plane.z[idx+1]):0.5*(plane.z[idx]+plane.z[idx-1]);
			
			plane.v_win_fac1[idx] = plane.spacing[idx]*2*plane.costheta[u_idx];
			plane.v_win_fac2[idx] = plane.costheta[u_idx]*TX_MAX;
			plane.v_win_fac3[idx] = plane.sintheta[u_idx]*TY_MAX;
		}
		
		for(int j = 0; j<6; j++){
			int idx = i*6+j+1;
			plane.u_win[idx] = fabs(0.5*plane.scaley[u_idx]*plane.sintheta[u_idx]) + TX_MAX*fabs((plane.z_mean[u_idx] - plane.z_mean[x_idx])*plane.costheta[u_idx]) + TY_MAX*fabs((plane.z_mean[u_idx] - plane.z_mean[x_idx])*plane.sintheta[u_idx]) + 2.*plane.spacing[u_idx] + u_factor[i];
		}
		cout << u_idx << " " << plane.u_win[u_idx] << " = " << fabs(0.5*plane.scaley[u_idx]*plane.sintheta[u_idx]) << " + " << TX_MAX*fabs((plane.z_mean[u_idx] - plane.z_mean[x_idx])*plane.costheta[u_idx]) << " + " << TY_MAX*fabs((plane.z_mean[u_idx] - plane.z_mean[x_idx])*plane.sintheta[u_idx]) << " + " << 2.*plane.spacing[u_idx] + u_factor[i] << endl;
		cout << " u costheta " << plane.costheta[u_idx] << " u sintheta " << plane.sintheta[u_idx] << " x_span " << plane.scaley[u_idx] << " spacing " << plane.spacing[u_idx] << " z plane_u " << plane.z_mean[u_idx] << " z plane_x " << plane.z_mean[x_idx] << endl;  
	}
	cout << "Geometry file read out" << endl;
	
	double wire_position[55][400];//Let's keep this: simpler, more robust
	for(int i = 1; i <= nChamberPlanes; ++i){
		//cout << plane[i].nelem << endl;
      		for(int j = 1; j <= plane.nelem[i]; ++j){
          		double pos = (j - (plane.nelem[i]+1.)/2.)*plane.spacing[i] + plane.xoffset[i] + plane.x0[i]*plane.costheta[i] + plane.y0[i]*plane.sintheta[i] + plane.deltaW_[i*9];
			wire_position[i][j] = pos;
		}
	}
	for(int i = nChamberPlanes+1; i<=nChamberPlanes+nHodoPlanes; ++i){
		//cout << plane[i].nelem << endl;
	      	for(int j = 1; j <= plane.nelem[i]; ++j){
          		double pos = plane.x0[i]*plane.costheta[i] + plane.y0[i]*plane.sintheta[i] + plane.xoffset[i] + (j - (plane.nelem[i]+1)/2.)*plane.spacing[i] + plane.deltaW_[i*9];
			wire_position[i][j] = pos;
		}
	}
	for(int i = nChamberPlanes+nHodoPlanes+1; i<=nChamberPlanes+nHodoPlanes+nPropPlanes; ++i){
		//cout << plane[i].nelem << endl;
	      	for(int j = 1; j <= plane.nelem[i]; ++j){
          		int moduleID = 8 - int((j - 1)/8);
			//cout << moduleID << endl;
             		double pos = plane.x0[i]*plane.costheta[i] + plane.y0[i]*plane.sintheta[i] + plane.xoffset[i] + (j - (plane.nelem[i]+1)/2.)*plane.spacing[i] + plane.deltaW_[i*9+moduleID];
			wire_position[i][j] = pos;
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

#ifdef E1039
	std::vector<SQHit*> hit_vec;
#endif
	
	if(e906data){
		dataTree = (TTree *)dataFile->Get("save");
		dataTree->SetBranchAddress("rawEvent", &rawEvent);
	}else{
#ifdef E1039
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
#endif
	}
	int nEvtMax = dataTree->GetEntries();
	if(nEvtMax>EstnEvtMax)nEvtMax=EstnEvtMax;

	static gEvent host_gEvent;
	static gEventHitCollections host_gEventHits;
	
	int nAH, nTH;
	
	unsigned int detarrayoffset[nDetectors];
	unsigned int hitarrayoffset[nDetectors];
	//the offset convention for the hit is the following: 
	// 1 <= detID <= 30: EstnEvtMax*nChamberPlanes*5*datasizes::NMaxHitsChambers*(detID-1)
	// 31 <= detID <= 46: EstnEvtMax*nHodoPlanes*4*datasizes::NMaxHitsHodoscopes*(detID-31)
	// 47 <= detID <= 54: EstnEvtMax*nPropPlanes*5*datasizes::NMaxHitsPropTubes*(detID-47)
	for(short k = 1; k<=nChamberPlanes; k++){
		detarrayoffset[k] = EstnEvtMax*(k-1);
		hitarrayoffset[k] = EstnEvtMax*datasizes::NHitsParam*datasizes::NMaxHitsChambers*(k-1);
	}
	for(short k = nChamberPlanes+1; k<=nChamberPlanes+nHodoPlanes; k++){
		detarrayoffset[k] = EstnEvtMax*(k-1);
		hitarrayoffset[k] = EstnEvtMax*datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes*(k-31);
	}
	for(short k = nChamberPlanes+nHodoPlanes+1; k<=nDetectors; k++){
		detarrayoffset[k] = EstnEvtMax*(k-1);
		hitarrayoffset[k] = EstnEvtMax*datasizes::NHitsParam*datasizes::NMaxHitsPropTubes*(k-47);
	}
	
	for(short k = 1; k<=nDetectors; k++){
		cout << detarrayoffset[k] << " " << hitarrayoffset[k] << endl;
	}
	
	short detid;
	int nhits;
	
	cout << "unfolding " << nEvtMax <<" events" << endl;
	// loop on event: get RawEvent information and load it into gEvent
	for(int i = 0; i < nEvtMax; ++i) {
		if(i%1000==0)cout << i << "/" << nEvtMax <<  endl;
		dataTree->GetEntry(i);
		
		//cout<<"Converting "<<i<<"/"<<nEvtMax<<endl;
		if(e906data){
			host_gEvent.RunID[i] = rawEvent->fRunID;
			host_gEvent.EventID[i] = rawEvent->fEventID;
			host_gEvent.SpillID[i] = rawEvent->fSpillID;
			host_gEvent.TriggerBits[i] = rawEvent->fTriggerBits;
			host_gEvent.TargetPos[i] = rawEvent->fTargetPos;
			host_gEvent.TurnID[i] = rawEvent->fTurnID;
			host_gEvent.RFID[i] = rawEvent->fRFID;
			for(int j=0; j<33; j++) {
				host_gEvent.Intensity[i*33+j] = rawEvent->fIntensity[j];
			}
			host_gEvent.TriggerEmu[i] = rawEvent->fTriggerEmu;
			for(int k=0; k<4; k++) {
				host_gEvent.NRoads[i*4+k] = rawEvent->fNRoads[k];
			}
			for(int l=0; l<nDetectors; l++) {
				host_gEvent.NHits[i*nDetectors+l] = rawEvent->fNHits[l];
				if(1 <= l && l <= 30){
					host_gEventHits.NHitsChambers[detarrayoffset[l]+l-1] = rawEvent->fNHits[l];
				}
				if(31 <= l && l <= 46){
					host_gEventHits.NHitsHodo[detarrayoffset[l]+l-31] = rawEvent->fNHits[l];
				}
				if(47 <= l && l <= 54){
					host_gEventHits.NHitsPropTubes[detarrayoffset[l]+l-47] = rawEvent->fNHits[l];
				}
			}
			host_gEvent.nAH[i] = rawEvent->fAllHits.size();
			host_gEvent.nTH[i] = rawEvent->fTriggerHits.size();
			
			for(int m=0; m<rawEvent->fAllHits.size(); m++) {
				detid = (rawEvent->fAllHits[m]).detectorID;
				nhits = rawEvent->fNHits[detid];
				if(1 <= detid && detid <= 30){
					host_gEventHits.HitsChambersRawData[hitarrayoffset[detid]+m] = (rawEvent->fAllHits[m]).elementID;
					host_gEventHits.HitsChambersRawData[hitarrayoffset[detid]+m*nhits] = wire_position[detid][(rawEvent->fAllHits[m]).elementID];
					host_gEventHits.HitsChambersRawData[hitarrayoffset[detid]+m*2*nhits] = (rawEvent->fAllHits[m]).tdcTime;
					host_gEventHits.HitsChambersRawData[hitarrayoffset[detid]+m*3*nhits] = (rawEvent->fAllHits[m]).flag;
					host_gEventHits.HitsChambersRawData[hitarrayoffset[detid]+m*4*nhits] = (rawEvent->fAllHits[m]).driftDistance;
				}

				if(31 <= detid && detid <= 46){
					host_gEventHits.HitsHodoRawData[hitarrayoffset[detid]+m] = (rawEvent->fAllHits[m]).elementID;
					host_gEventHits.HitsHodoRawData[hitarrayoffset[detid]+m*nhits] = wire_position[detid][(rawEvent->fAllHits[m]).elementID];
					host_gEventHits.HitsHodoRawData[hitarrayoffset[detid]+m*2*nhits] = (rawEvent->fAllHits[m]).tdcTime;
					host_gEventHits.HitsHodoRawData[hitarrayoffset[detid]+m*3*nhits] = (rawEvent->fAllHits[m]).flag;
					host_gEventHits.HitsHodoRawData[hitarrayoffset[detid]+m*4*nhits] = (rawEvent->fAllHits[m]).driftDistance;
				}
				
				if(47 <= detid && detid <= 54){
					host_gEventHits.HitsPropTubesRawData[hitarrayoffset[detid]+m] = (rawEvent->fAllHits[m]).elementID;
					host_gEventHits.HitsPropTubesRawData[hitarrayoffset[detid]+m*nhits] = wire_position[detid][(rawEvent->fAllHits[m]).elementID];
					host_gEventHits.HitsPropTubesRawData[hitarrayoffset[detid]+m*2*nhits] = (rawEvent->fAllHits[m]).tdcTime;
					host_gEventHits.HitsPropTubesRawData[hitarrayoffset[detid]+m*3*nhits] = (rawEvent->fAllHits[m]).flag;
					host_gEventHits.HitsPropTubesRawData[hitarrayoffset[detid]+m*4*nhits] = (rawEvent->fAllHits[m]).driftDistance;
				}
				//host_gEvent.AllHits[i*EstnAHMax+m].index=(rawEvent->fAllHits[m]).index;
				//host_gEvent.AllHits[i*EstnAHMax+m].detectorID=(rawEvent->fAllHits[m]).detectorID;
				//host_gEvent.AllHits[i*EstnAHMax+m].elementID=(rawEvent->fAllHits[m]).elementID;
				//host_gEvent.AllHits[i*EstnAHMax+m].tdcTime=(rawEvent->fAllHits[m]).tdcTime;
				//host_gEvent.AllHits[i*EstnAHMax+m].driftDistance=(rawEvent->fAllHits[m]).;
				//host_gEvent.AllHits[i*EstnAHMax+m].pos=wire_position[(rawEvent->fAllHits[m]).detectorID][(rawEvent->fAllHits[m]).elementID];
				//host_gEvent.AllHits[i*EstnAHMax+m].flag=(rawEvent->fAllHits[m]).flag;
			}
			//AFAIK trigger hits are not used anywhere in the online reco...
			//for(int n=0; n<rawEvent->fTriggerHits.size(); n++) {
				//host_gEvent.TriggerHits[i*EstnTHMax+n].index=(rawEvent->fTriggerHits[n]).index;
				//host_gEvent.TriggerHits[i*EstnTHMax+n].detectorID=(rawEvent->fTriggerHits[n]).detectorID;
				//host_gEvent.TriggerHits[i*EstnTHMax+n].elementID=(rawEvent->fTriggerHits[n]).elementID;
				//host_gEvent.TriggerHits[i*EstnTHMax+n].tdcTime=(rawEvent->fTriggerHits[n]).tdcTime;
				//host_gEvent.TriggerHits[i*EstnTHMax+n].driftDistance=(rawEvent->fTriggerHits[n]).driftDistance;
				//host_gEvent.TriggerHits[i*EstnTHMax+n].pos=wire_position[(rawEvent->fAllHits[n]).detectorID][(rawEvent->fAllHits[n]).elementID];
				//host_gEvent.TriggerHits[i*EstnTHMax+n].flag=(rawEvent->fTriggerHits[n]).flag;
			//}
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
#ifdef E1039
			//Default option: e1039
			//if(_event_id<20)cout << " evt: " << _event_id << " nhits = " << hit_vec.size() << endl; 
			host_gEvent.RunID[i] = _run_id;
			host_gEvent.SpillID[i] = _spill_id;
			host_gEvent.EventID[i] = _event_id;
			host_gEvent.TriggerBits[i] = _trigger;
			//for(int k = 0; k<4; k++)host_gEvent[i].NRoads[i*4+k] = _qie_presums[k];
			host_gEvent.TurnID[i] = _qie_turn_id;
			host_gEvent.RFID[i] = _qie_rf_id;
			//for(int k = 0; k<33; k++)host_gEvent.Intensity[i*33+k] = _qie_rf_inte[k];
			
			for(int k = 0; k<nDetectors; k++)host_gEvent[i].NHits[k] = 0;//we will increment those in the vector hit loop
			int ntrighits = 0;
			host_gEvent.nAH[i] = 0;//hit_vec.size();
			for(int m = 0; m<hit_vec.size(); m++){
				if(hit_vec[m]->get_detector_id()>54){
					//if(_event_id<20)cout << " dark photon plane hit! " << hit_vec[m]->get_detector_id() << endl;
					continue;
					//dark photon planes; I don't think we care about those for the purpose of online reconstruction... do we?
				}
				//host_gEvent.nAH[i]++;
				//host_gEvent.NHits[i*nDetectors+hit_vec[m]->get_detector_id()]++;
				//host_gEvent.AllHits[i+EstnAHMax*m].index=hit_vec[m]->get_hit_id();
				//host_gEvent.AllHits[i+EstnAHMax*m].detectorID=hit_vec[m]->get_detector_id();
				//host_gEvent.AllHits[i+EstnAHMax*m].elementID=hit_vec[m]->get_element_id();
				//host_gEvent.AllHits[i+EstnAHMax*m].tdcTime=hit_vec[m]->get_tdc_time();
				//host_gEvent.AllHits[i+EstnAHMax*m].driftDistance=fabs(hit_vec[m]->get_drift_distance());
				//host_gEvent.AllHits[i+EstnAHMax*m].sign_mc=hit_vec[m]->get_drift_distance()/fabs(hit_vec[m]->get_drift_distance());
				//host_gEvent.AllHits[i+EstnAHMax*m].pos=wire_position[hit_vec[m]->get_detector_id()][hit_vec[m]->get_element_id()];
				//host_gEvent.AllHits[i+EstnAHMax*m].flag=(1<<hit_vec[m]->is_in_time());
				
				//if(hit_vec[m]->is_trigger_mask()){
					//host_gEvent.TriggerHits[i*EstnTHMax+ntrighits].index=hit_vec[m]->get_hit_id();
					//host_gEvent.TriggerHits[i*EstnTHMax+ntrighits].detectorID=hit_vec[m]->get_detector_id();
					//host_gEvent.TriggerHits[i*EstnTHMax+ntrighits].elementID=hit_vec[m]->get_element_id();
					//host_gEvent.TriggerHits[i*EstnTHMax+ntrighits].tdcTime=hit_vec[m]->get_tdc_time();
					//host_gEvent.TriggerHits[i*EstnTHMax+ntrighits].driftDistance=fabs(hit_vec[m]->get_drift_distance());
					//host_gEvent.TriggerHits[i*EstnTHMax+ntrighits].pos=wire_position[hit_vec[m]->get_detector_id()][hit_vec[m]->get_element_id()];
					//host_gEvent.TriggerHits[i*EstnTHMax+ntrighits].flag=(1<<hit_vec[m]->is_in_time());
				//	ntrighits++;
				//}
			}
			host_gEvent.nTH[i] = ntrighits;
#endif
		}
	}
	cout << "loaded events" << endl;
	auto cp2 = std::chrono::system_clock::now();

	auto evt_prep = cp2-cp1;
	cout<<"Read/prepare events: "<<evt_prep.count()/1000000000.<<endl;

	// evaluate the total size of the gEvent array (and the SW array) for memory allocation 
	// (the memory cannot be dynamically allocated) 
	size_t NBytesAllEvent = sizeof(gEvent);
	size_t NBytesAllHits = sizeof(gEventHitCollections);
	size_t NBytesAllOutputEvent = sizeof(gOutputEvent);
	size_t NBytesAllPlanes =  sizeof(gPlane);
	//size_t NBytesFitterTools = sizeof(gStraightFitArrays);
	//size_t NBytesStraightTrackBuilders = sizeof(gStraightTrackBuilder);
	//size_t NBytesFullTrackBuilders = sizeof(gFullTrackBuilder);
	//size_t NBytesKalmanFilterTools = sizeof(gStraightFitArrays);
	
	cout << "Total size allocated on GPUs " << NBytesAllEvent+NBytesAllOutputEvent+NBytesAllPlanes+NBytesAllHits << endl;
	cout << " input events: " << NBytesAllEvent << "; output events: " << NBytesAllOutputEvent 
		<< "; raw hits: " << NBytesAllHits 
		//<< "; straight track builder tools: " << NBytesStraightTrackBuilders
	//     << "; fitter tools: " << NBytesFitterTools << "; straight track builders: " << NBytesStraightTrackBuilders 
	//     << "; full track builders: " << NBytesFullTrackBuilders << "; kalman filters: " << NBytesKalmanFilterTools
		<< "; planes info: " << NBytesAllPlanes << endl;  
		
	gEvent* host_output_eR = (gEvent*)malloc(NBytesAllEvent);
	gEvent* host_output_gHits = (gEvent*)malloc(NBytesAllHits);
	gOutputEvent* host_output_TKL = (gOutputEvent*)malloc(NBytesAllOutputEvent);
	
	// declaring gEvent objects for the device (GPU) to use.
	gEvent* device_gEvent;
	gEventHitCollections* device_gHits;
	gOutputEvent* device_output_TKL;
	gPlane* device_gPlane;
	//gStraightFitArrays* device_gFitArrays;
	//gStraightTrackBuilder* device_gStraightTrackBuilder;
	//gFullTrackBuilder* device_gFullTrackBuilder;
	//gKalmanFitArrays* device_gKalmanFitArrays;


	//printDeviceStatus();

	// copy of data from host to device: evaluate operation time 
	// Allocating memory for GPU (pointer to allocated device ); check for errors in the process; stops the program if issues encountered
	gpuErrchk( cudaMalloc((void**)&device_gEvent, NBytesAllEvent));
	gpuErrchk( cudaMalloc((void**)&device_gHits, NBytesAllHits));
	gpuErrchk( cudaMalloc((void**)&device_output_TKL, NBytesAllOutputEvent));
	//allocating the memory for the planes
	gpuErrchk( cudaMalloc((void**)&device_gPlane, NBytesAllPlanes));
	//gpuErrchk( cudaMalloc((void**)&device_gFitArrays, NBytesFitterTools));
	//gpuErrchk( cudaMalloc((void**)&device_gStraightTrackBuilder, NBytesStraightTrackBuilders));
	
	std::size_t free_bytes;
	std::size_t total_bytes;

	CUDA_CHECK_STATUS(cudaMemGetInfo(&free_bytes, &total_bytes));
    	cout << "Current memory foot print: " << free_bytes << " / " << total_bytes << endl;
	
	// cudaMemcpy(dst, src, count, kind): copies data between host and device:
	// dst: destination memory address; src: source memory address; count: size in bytes; kind: type of transfer
	gpuErrchk( cudaMemcpy(device_gPlane, &plane, NBytesAllPlanes, cudaMemcpyHostToDevice));
	gpuErrchk( cudaMemcpy(device_gEvent, &host_gEvent, NBytesAllEvent, cudaMemcpyHostToDevice));
	gpuErrchk( cudaMemcpy(device_gHits, &host_gEventHits, NBytesAllHits, cudaMemcpyHostToDevice));
		
	auto cp3 = std::chrono::system_clock::now();

	auto cp_to_gpu = cp3-cp2;
	cout<<"Copy to GPU: "<<cp_to_gpu.count()/1000000000.<<endl;
		
	// now data is transfered in the device: kernel function for event reconstruction called;
	// note that the function call is made requesting a number of blocks and a number of threads per block
	// in practice we have as many threads total as number of events; 
	//gkernel_eR<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent);
	gkernel_eR<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_gHits);
	
	// check status of device and synchronize;
	size_t nEvents = EstnEvtMax;
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	auto cp4 = std::chrono::system_clock::now();
	auto gpu_er = cp4-cp3;
	cout<<"GPU: event reducing: "<<gpu_er.count()/1000000000.<<endl;
	
	//gKernel_XZ_YZ_tracking_new<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gStraightTrackBuilder, device_gFitArrays, device_gPlane);
	
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	auto cp5 = std::chrono::system_clock::now();
	auto gpu_st = cp5-cp4;
	cout<<"GPU: straight tracking: "<<gpu_st.count()/1000000000.<<endl;

	//release here the memory for straight track builders and straight track fitters	
	//cudaFree( device_gStraightTrackBuilder );
	
	//gpuErrchk( cudaMalloc((void**)&device_gFullTrackBuilder, NBytesFullTrackBuilders));

	cout << "Current memory foot print: " << free_bytes << " / " << total_bytes << endl;
	
	//gKernel_GlobalTrack_building<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gFullTrackBuilder, device_gFitArrays, device_gPlane, 1);

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	//gpuErrchk( cudaMalloc((void**)&device_gKalmanFitArrays, NBytesKalmanFilterTools));
	
	auto cp6 = std::chrono::system_clock::now();
	auto gpu_gt = cp6-cp5;
	cout<<"GPU: global tracking: "<<gpu_gt.count()/1000000000.<<endl;

	//gKernel_GlobalTrack_KalmanFitting<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_output_TKL, device_gKalmanFitArrays, device_gPlane);

	//gpuErrchk( cudaPeekAtLastError() );
	//gpuErrchk( cudaDeviceSynchronize() );
	

	auto cp7 = std::chrono::system_clock::now();
	auto gpu_kf = cp7-cp6;
	cout<<"GPU: kalman filtering: "<<gpu_gt.count()/1000000000.<<endl;

	// data transfer from device to host
	gpuErrchk( cudaMemcpy(host_output_eR, device_gEvent, NBytesAllEvent, cudaMemcpyDeviceToHost));
	cudaFree(device_gEvent);

	gpuErrchk( cudaMemcpy(host_output_gHits, device_gHits, NBytesAllHits, cudaMemcpyDeviceToHost));
	cudaFree(device_gHits);
	
	gpuErrchk( cudaMemcpy(host_output_TKL, device_output_TKL, NBytesAllOutputEvent, cudaMemcpyDeviceToHost));
	cudaFree(device_output_TKL);

	auto cp8 = std::chrono::system_clock::now();
	auto cp_to_cpu = cp8-cp7;
	cout<<"Copy back to CPU: "<<cp_to_cpu.count()/1000000000.<<endl;
//#define TEST 1
#ifdef TEST	
	ofstream out("OutputFile.txt");
	//Write in a file, 
	long tklctr = 0;
	long nEvtsTotal = 0;
	long nEvtsPass = 0;
	for(int n = 0; n<nEvtMax; n++){
		if(host_output_eR->nAH[n]==0)continue;
		nEvtsTotal++;
		if(host_output_eR->HasTooManyHits[n])continue;
		nEvtsPass++;
		out<<n<<" "<< host_output_eR->nAH[n] <<" "<< host_output_TKL->nTracklets[n] <<endl;
		tklctr+= host_output_TKL->nTracklets[n];
		for(int k = 1; k<=nDetectors; k++ ){
			out << host_output_eR->NHits[n*EstnAHMax+k] << " ";
		}out<<endl;
		
		for(int k = 0; k<host_output_eR->nAH[n]; k++ ){
			out << host_output_eR->AllHits[n*EstnAHMax+k].detectorID << " " << host_output_eR->AllHits[n*EstnAHMax+k].elementID << " " << host_output_eR->AllHits[k].driftDistance*host_output_eR->AllHits[n*EstnAHMax+k].sign_mc << endl;
		}
		
		for(int k = 0; k<host_output_TKL->nTracklets[n]; k++ ){
			if(isnan(host_output_TKL->AllTracklets[n*TrackletSizeMax+k].x0))host_output_TKL->AllTracklets[n*TrackletSizeMax+k].x0 = -1000;
			if(isnan(host_output_TKL->AllTracklets[n*TrackletSizeMax+k].y0))host_output_TKL->AllTracklets[n*TrackletSizeMax+k].y0 = -1000;
			if(isnan(host_output_TKL->AllTracklets[n*TrackletSizeMax+k].invP))host_output_TKL->AllTracklets[n*TrackletSizeMax+k].invP = -1.0;
			if(isnan(host_output_TKL->AllTracklets[n*TrackletSizeMax+k].tx))host_output_TKL->AllTracklets[n*TrackletSizeMax+k].tx = -1.0;
			if(isnan(host_output_TKL->AllTracklets[n*TrackletSizeMax+k].ty))host_output_TKL->AllTracklets[n*TrackletSizeMax+k].ty = -1.0;
			out << host_output_TKL->AllTracklets[n*TrackletSizeMax+k].stationID << " " << host_output_TKL->AllTracklets[n*TrackletSizeMax+k].x0 << " " << host_output_TKL->AllTracklets[n*TrackletSizeMax+k].y0 << " " << host_output_TKL->AllTracklets[n*TrackletSizeMax+k].tx << " " << host_output_TKL->AllTracklets[n*TrackletSizeMax+k].ty << " " << host_output_TKL->AllTracklets[n*TrackletSizeMax+k].invP << " " << host_output_TKL->AllTracklets[n*TrackletSizeMax+k].nXHits+host_output_TKL->AllTracklets[n*TrackletSizeMax+k].nUHits+host_output_TKL->AllTracklets[n*TrackletSizeMax+k].nVHits << endl;
			//if(n<100 && host_output_TKL[n].nTracklets>1)cout << n << " " << host_output_TKL[n].AllTracklets[k].stationID << " " << host_output_TKL[n].AllTracklets[k].nXHits<< " " <<host_output_TKL[n].AllTracklets[k].nUHits<< " " <<host_output_TKL[n].AllTracklets[k].nVHits << endl;
			for(int l = 0; l<host_output_TKL->AllTracklets[n*TrackletSizeMax+k].nXHits+host_output_TKL->AllTracklets[n*TrackletSizeMax+k].nUHits+host_output_TKL->AllTracklets[n*TrackletSizeMax+k].nVHits; l++){
				out << host_output_TKL->AllTracklets[n*TrackletSizeMax+k].hits[l].detectorID << " " << host_output_TKL->AllTracklets[n*TrackletSizeMax+k].hits[l].elementID << " " << host_output_TKL->AllTracklets[n*TrackletSizeMax+k].hits[l].driftDistance*host_output_TKL->AllTracklets[n*TrackletSizeMax+k].hitsign[l] << " " << host_output_TKL->AllTracklets[n*TrackletSizeMax+k].hits[l].pos << endl;
			}
		}
		
	}

	cout << tklctr << " total tracks reconstructed" << endl;
	cout << nEvtsPass << " evts with low enough number of hits on " << nEvtsTotal << " events total." << endl; 
	//auto end_kernel = std::chrono::system_clock::now();

#endif		

	delete rawEvent;

	//TFile* outFile = new TFile(outputFile.Data(), "RECREATE");
	//ORoutput_tree* output = new ORoutput_tree();
	//for(int i = 0; i < nEvtMax; ++i) {
	//	output->Clear();
	//	for(int k = 1; k<=nDetectors; k++ )output->fNhitsReduced[k] = host_output_eR[i].NHits[k];
	//	//output->Write();
	//}
	//output->Write();

	auto cp9 = std::chrono::system_clock::now();
	auto write_output = cp9-cp8;
	cout<<"Write Output: "<<write_output.count()/1000000000.<<endl;

	// printing the time required for all operations
	auto end = std::chrono::system_clock::now();
	auto overall = end - start;
	cout<<"Total time: "<<overall.count()/1000000000.<<endl;
	return 0;
}

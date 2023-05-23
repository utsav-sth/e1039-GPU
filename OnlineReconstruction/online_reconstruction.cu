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
#ifdef ROOTSAVE
#include <TRandom.h>
#include <TMatrixD.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TStopwatch.h>
#include <TTimeStamp.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1D.h>
#endif
//#include "LoadInput.h"
#include "OROutput.h"
#include "reconstruction_kernels.cuh"

#ifdef GLDISPLAY
#include "display_utils.cuh"
#endif
//#else

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

#define CUDA_CHECK_STATUS( cuda_function_call ) \
        if (cudaError_t const status = cuda_function_call) \
        { \
            throw std::runtime_error( cudaGetErrorString( status ) ) ; \
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

#ifdef OLDCODE

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

<<<<<<< HEAD
#endif

=======
>>>>>>> e04fbf96b0c9a61c7bd2cf5732b68f3000b4ed2d

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
#ifdef E1039
	if(argn>4)e906data = atoi(argv[4]);
#endif
	
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
	      //TODO: solve the mixing of x1, x2 / y1, y2
	      //iss >> ipl >> z >> nelem >> cellwidth >> spacing >> xoffset >> scalex >> x0 >> x1 >> x2 >> costheta >> scaley >> y0 >> y1 >> y2 >> sintheta >> resolution >> p1x >> p1y >> p1z >> deltapx >> deltapy >> deltapz >> dp1x >> dp1y >> dp1z;
	      iss >> ipl >> z >> nelem >> cellwidth >> spacing >> xoffset >> scalex >> x0 >> y1 >> y2 >> costheta >> scaley >> y0 >> x1 >> x2 >> sintheta >> resolution >> p1x >> p1y >> p1z >> deltapx >> deltapy >> deltapz >> dp1x >> dp1y >> dp1z;
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
	cout << "Geometry file read out" << endl;
	
	double wire_position[55][400];//Let's keep this: simpler, more robust
	for(int i = 1; i <= nChamberPlanes; ++i){
		//cout << plane[i].nelem << endl;
      		for(int j = 1; j <= plane.nelem[i]; ++j){
          		double pos = (j - (plane.nelem[i]+1.)/2.)*plane.spacing[i] + plane.xoffset[i] + plane.x0[i]*plane.costheta[i] + plane.y0[i]*plane.sintheta[i] + plane.deltaW_[i*9];
			wire_position[i][j] = pos;
		}
	}
<<<<<<< HEAD
	for(int i = nChamberPlanes+1; i<=nChamberPlanes+nHodoPlanes; ++i){
		//cout << plane[i].nelem << endl;
	      	for(int j = 1; j <= plane.nelem[i]; ++j){
          		double pos = plane.x0[i]*plane.costheta[i] + plane.y0[i]*plane.sintheta[i] + plane.xoffset[i] + (j - (plane.nelem[i]+1)/2.)*plane.spacing[i] + plane.deltaW_[i*9];
=======
	for(int i = nChamberPlanes; i<nChamberPlanes+nHodoPlanes; ++i){
		cout << plane[i].nelem << endl;
	      	for(int j = 0; j < plane[i].nelem; ++j){
          		double pos = plane[i].x0*plane[i].costheta + plane[i].y0*plane[i].sintheta + plane[i].xoffset + (j - (plane[i].nelem+1)/2.)*plane[i].spacing + plane[i].deltaW_[0];
>>>>>>> e04fbf96b0c9a61c7bd2cf5732b68f3000b4ed2d
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
	
	unsigned int evhitarrayoffset[nDetectors];
	//the offset convention for the hit is the following: for event i: 
	// 1 <= detID <= 30: i*eventhitsize[0]+datasizes::NHitsParam*datasizes::NMaxHitsChambers*(detID-1)
	// 31 <= detID <= 46: +nHodoPlanes*datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes*(detID-31) 
	// 47 <= detID <= 54: +nPropPlanes*datasizes::NHitsParam*datasizes::NMaxHitsPropTubes*(detID-47)
	for(short k = 1; k<=nChamberPlanes; k++){
		evhitarrayoffset[k] = datasizes::NHitsParam*datasizes::NMaxHitsChambers*(k-1);
	}
	for(short k = nChamberPlanes+1; k<=nChamberPlanes+nHodoPlanes; k++){
		evhitarrayoffset[k] = datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes*(k-31);
	}
	for(short k = nChamberPlanes+nHodoPlanes+1; k<=nDetectors; k++){
		evhitarrayoffset[k] = datasizes::NHitsParam*datasizes::NMaxHitsPropTubes*(k-47);
	}
#ifdef DEBUG
	for(short k = 1; k<=nDetectors; k++){
		cout << k << " " << evhitarrayoffset[k] << endl;
	}
#endif
	
	
	short detid;
	
	int nhits;
	//the hit offset is to give an individual offset for each hit detector collection per event
	int hit_ctr[nDetectors];
	int firstevent;
	bool isFPGAtriggered;
	cout << "unfolding " << nEvtMax <<" events" << endl;
	// loop on event: get RawEvent information and load it into gEvent
	for(int i = 0; i < nEvtMax; ++i) {
		if(i%1000==0)cout << i << "/" << nEvtMax <<  endl;
		dataTree->GetEntry(i);
		
		//cout<<"Converting "<<i<<"/"<<nEvtMax<<endl;
		if(e906data){
			if(i==0)firstevent = rawEvent->fEventID;

			//host_gEvent.RunID[i] = rawEvent->fRunID;
			host_gEvent.EventID[i] = rawEvent->fEventID;
			//host_gEvent.SpillID[i] = rawEvent->fSpillID;
			host_gEvent.TriggerBits[i] = rawEvent->fTriggerBits;
			isFPGAtriggered = false;
			for(int k = 0; k<5; k++){
				if( (rawEvent->fTriggerBits & k) != 0 )isFPGAtriggered = true;
			}
			//host_gEvent.TargetPos[i] = rawEvent->fTargetPos;
			//host_gEvent.TurnID[i] = rawEvent->fTurnID;
			//host_gEvent.RFID[i] = rawEvent->fRFID;
			//for(int j=0; j<33; j++) {
			//	host_gEvent.Intensity[i*33+j] = rawEvent->fIntensity[j];
			//}
			//host_gEvent.TriggerEmu[i] = rawEvent->fTriggerEmu;
			for(int k=0; k<4; k++) {
				host_gEvent.NRoads[i*4+k] = rawEvent->fNRoads[k];
			}
			for(int l=1; l<nDetectors; l++) {
				host_gEvent.NHits[i*nDetectors+l-1] = rawEvent->fNHits[l];
				hit_ctr[l] = 0;
				if(1 <= l && l <= 30){
					host_gEventHits.NHitsChambers[(rawEvent->fEventID-firstevent)*nChamberPlanes+l-1] = rawEvent->fNHits[l];
				}
				if(31 <= l && l <= 46){
					host_gEventHits.NHitsHodo[(rawEvent->fEventID-firstevent)*nHodoPlanes+l-31] = rawEvent->fNHits[l];
#ifdef DEBUG
					if(rawEvent->fEventID==debug::EvRef+firstevent)cout << l << " " << rawEvent->fNHits[l] << " " << (rawEvent->fEventID-firstevent)*nHodoPlanes+l-31 << " " 
						<< host_gEventHits.NHitsHodo[(rawEvent->fEventID-firstevent)*nHodoPlanes+l-31] << endl;
#endif
				}
				if(47 <= l && l <= 54){
					host_gEventHits.NHitsPropTubes[(rawEvent->fEventID-firstevent)*nPropPlanes+l-47] = rawEvent->fNHits[l];
#ifdef DEBUG
					if(rawEvent->fEventID==debug::EvRef+firstevent)cout << l << " " << rawEvent->fNHits[l] << " " << (rawEvent->fEventID-firstevent)*nPropPlanes+l-47 << " " 
						<< host_gEventHits.NHitsPropTubes[rawEvent->fEventID*nPropPlanes+l-47] << endl;
#endif
				}
			}
			host_gEvent.nAH[i] = rawEvent->fAllHits.size();
			host_gEvent.nTH[i] = rawEvent->fTriggerHits.size();
						
			for(int m=0; m<rawEvent->fAllHits.size(); m++) {
				detid = (rawEvent->fAllHits[m]).detectorID;
				nhits = rawEvent->fNHits[detid];
#ifdef DEBUG
				if(rawEvent->fEventID==debug::EvRef+firstevent){
					cout << detid << " " << (rawEvent->fAllHits[m]).elementID << " " 
						<< wire_position[detid][(rawEvent->fAllHits[m]).elementID] << " " 
						<< (rawEvent->fAllHits[m]).tdcTime << " " << (rawEvent->fAllHits[m]).flag << " " 
						<< (rawEvent->fAllHits[m]).driftDistance << endl;
				}
#endif				
				if(1 <= detid && detid <= 30){
					if((rawEvent->fEventID-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+(hit_ctr[detid])*4*nhits > EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsChambers)
					cout << rawEvent->fEventID << " " << detid <<  " " << (rawEvent->fEventID-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+(hit_ctr[detid])+4*nhits << " " << EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsChambers << " " << hit_ctr[detid] << endl;
					
#ifdef DEBUG
					if(rawEvent->fEventID==debug::EvRef+firstevent)cout << "hit offsets " << detid << " " << (rawEvent->fEventID-firstevent)*datasizes::eventhitsize[0] << " " << evhitarrayoffset[detid] << " " << (rawEvent->fEventID-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid] << ": " << (rawEvent->fEventID-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+hit_ctr[detid] << " " << (rawEvent->fAllHits[m]).elementID << " " << (rawEvent->fEventID-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+(hit_ctr[detid])+nhits << " " << wire_position[detid][(rawEvent->fAllHits[m]).elementID] << endl;
#endif										
					host_gEventHits.HitsChambersRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+(hit_ctr[detid])] = (float)(rawEvent->fAllHits[m]).elementID;
					host_gEventHits.HitsChambersRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+(hit_ctr[detid])+nhits] = (float)wire_position[detid][(rawEvent->fAllHits[m]).elementID];
					host_gEventHits.HitsChambersRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+(hit_ctr[detid])+2*nhits] = (float)(rawEvent->fAllHits[m]).tdcTime;
					host_gEventHits.HitsChambersRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+(hit_ctr[detid])+3*nhits] = (float)(rawEvent->fAllHits[m]).flag;
					host_gEventHits.HitsChambersRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+(hit_ctr[detid])+4*nhits] = (float)(rawEvent->fAllHits[m]).driftDistance;
					hit_ctr[detid]++;
				}

				if(31 <= detid && detid <= 46){
					if((rawEvent->fEventID-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+(hit_ctr[detid])*4*nhits > EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes)
					cout << rawEvent->fEventID << " " << detid <<  " " << (rawEvent->fEventID-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+(hit_ctr[detid])+4*nhits << " " << EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes << " " << hit_ctr[detid] << endl;
					//if(isFPGAtriggered){
					//	if(33<= detid && detid <= 36)continue;
					//	if(41<= detid && detid <= 44)continue;
					//}else{
					//	if(37<= detid && detid <= 40)continue;
					//	if(32>= detid || detid >= 45)continue;
					//}
#ifdef DEBUG
					if(rawEvent->fEventID==debug::EvRef+firstevent)cout << hit_ctr[detid] << " " << detid << " " << (rawEvent->fAllHits[m]).elementID << " " << wire_position[detid][(rawEvent->fAllHits[m]).elementID] << " " << (rawEvent->fAllHits[m]).tdcTime << " " << (rawEvent->fAllHits[m]).flag << " " << (rawEvent->fAllHits[m]).driftDistance << endl;
#endif					
					host_gEventHits.HitsHodoRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+(hit_ctr[detid])] = (float)(rawEvent->fAllHits[m]).elementID;
					host_gEventHits.HitsHodoRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+(hit_ctr[detid])+nhits] = (float) wire_position[detid][(rawEvent->fAllHits[m]).elementID];
					host_gEventHits.HitsHodoRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+(hit_ctr[detid])+2*nhits] = (float)(rawEvent->fAllHits[m]).tdcTime;
					host_gEventHits.HitsHodoRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+(hit_ctr[detid])+3*nhits] = (float)(rawEvent->fAllHits[m]).flag;
					host_gEventHits.HitsHodoRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+(hit_ctr[detid])*4*nhits] = (float)(rawEvent->fAllHits[m]).driftDistance;
					hit_ctr[detid]++;
				}
				
				if(47 <= detid && detid <= 54){
					if( (rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+(hit_ctr[detid])+4*nhits > EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsPropTubes)
					cout << rawEvent->fEventID << " " << detid <<  " " << (rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+(hit_ctr[detid])*4*nhits << " " << EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsPropTubes << " " << hit_ctr[detid] << endl;

#ifdef DEBUG
					if(rawEvent->fEventID==debug::EvRef+firstevent)cout << "hit offsets " << detid << " " << (rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2] << " " << evhitarrayoffset[detid] << " " << (rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid] << ": " << (rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+hit_ctr[detid] << " " << (rawEvent->fAllHits[m]).elementID << " " << (rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+(hit_ctr[detid])+nhits << " " << wire_position[detid][(rawEvent->fAllHits[m]).elementID] << endl;
#endif
					host_gEventHits.HitsPropTubesRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+(hit_ctr[detid])] = (float)(rawEvent->fAllHits[m]).elementID;
					host_gEventHits.HitsPropTubesRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+(hit_ctr[detid])+nhits] = (float) wire_position[detid][(rawEvent->fAllHits[m]).elementID];
					host_gEventHits.HitsPropTubesRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+(hit_ctr[detid])+2*nhits] = (float) (rawEvent->fAllHits[m]).tdcTime;
					host_gEventHits.HitsPropTubesRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+(hit_ctr[detid])+3*nhits] = (float) (rawEvent->fAllHits[m]).flag;
					host_gEventHits.HitsPropTubesRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+(hit_ctr[detid])+4*nhits] = (float) (rawEvent->fAllHits[m]).driftDistance;
#ifdef DEBUG					
					if(rawEvent->fEventID==debug::EvRef+firstevent)cout << host_gEventHits.HitsPropTubesRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+(hit_ctr[detid])] << " " << host_gEventHits.HitsPropTubesRawData[(rawEvent->fEventID-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+(hit_ctr[detid])+nhits] << endl;
#endif
					hit_ctr[detid]++;
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
			if(i==0)firstevent = _event_id;
			//if(_event_id<20)cout << " evt: " << _event_id << " + " << firstevent << " nhits = " << hit_vec.size() << endl; 
			//host_gEvent.RunID[i] = _run_id;
			//host_gEvent.SpillID[i] = _spill_id;
			host_gEvent.EventID[i] = _event_id;
			host_gEvent.TriggerBits[i] = _trigger;
			//for(int k = 0; k<4; k++)host_gEvent[i].NRoads[i*4+k] = _qie_presums[k];
			//host_gEvent.TurnID[i] = _qie_turn_id;
			//host_gEvent.RFID[i] = _qie_rf_id;
			//for(int k = 0; k<33; k++)host_gEvent.Intensity[i*33+k] = _qie_rf_inte[k];
			
			for(int k = 0; k<nDetectors; k++)host_gEvent.NHits[i*nDetectors+k-1] = 0;//we will increment those in the vector hit loop
			int ntrighits = 0;
			host_gEvent.nAH[i] = 0;//hit_vec.size();
			for(int m = 0; m<hit_vec.size(); m++){
				if(hit_vec[m]->get_detector_id()>54){
					//if(_event_id<20)cout << " dark photon plane hit! " << hit_vec[m]->get_detector_id() << endl;
					continue;
					//dark photon planes; I don't think we care about those for the purpose of online reconstruction... do we?
				}
				host_gEvent.nAH[i]++;
				host_gEvent.NHits[i*nDetectors+hit_vec[m]->get_detector_id()-1]++;
			}
			host_gEvent.nTH[i] = ntrighits;
			
			for(int l=1; l<nDetectors; l++) {
				hit_ctr[l] = 0;//counter needs to be initialized properly...
				if(1 <= l && l <= 30){
#ifdef DEBUG
					if(_event_id==debug::EvRef+firstevent)cout << " det offset ("<< l << "): " << (_event_id-firstevent)*nChamberPlanes+l-1 << endl;
#endif
					host_gEventHits.NHitsChambers[(_event_id-firstevent)*nChamberPlanes+l-1] = host_gEvent.NHits[i*nDetectors+l-1];
				}
				if(31 <= l && l <= 46){
					host_gEventHits.NHitsHodo[(_event_id-firstevent)*nHodoPlanes+l-31] = host_gEvent.NHits[i*nDetectors+l-1];
				}
				if(47 <= l && l <= 54){
					host_gEventHits.NHitsPropTubes[(_event_id-firstevent)*nPropPlanes+l-47] = host_gEvent.NHits[i*nDetectors+l-1];
#ifdef DEBUG
					if(_event_id==debug::EvRef+firstevent)cout << l << " " << rawEvent->fNHits[l] << " " << _event_id*nPropPlanes+l-47 << " " 
						<< host_gEventHits.NHitsPropTubes[_event_id*nPropPlanes+l-47] << endl;
#endif
				}
			}
			//reloop :(
			for(int m = 0; m<hit_vec.size(); m++){
				if(hit_vec[m]->get_detector_id()>54){
					//if(_event_id<20)cout << " dark photon plane hit! " << hit_vec[m]->get_detector_id() << endl;
					continue;
					//dark photon planes; I don't think we care about those for the purpose of online reconstruction... do we?
				}
				detid = hit_vec[m]->get_detector_id();
				nhits = host_gEvent.NHits[i*nDetectors+detid-1];
#ifdef DEBUG
				if(_event_id==debug::EvRef+firstevent){
					cout << detid << " " << hit_vec[m]->get_element_id() << " " 
						<< wire_position[detid][ hit_vec[m]->get_element_id() ] << " " 
						<< hit_vec[m]->get_tdc_time() << " " << (1<<hit_vec[m]->is_in_time()) << " " 
						<< fabs(hit_vec[m]->get_drift_distance()) << endl;
				}
#endif				

				/*
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
					ntrighits++;
				//}
				*/

				//host_gEvent.AllHits[i+EstnAHMax*m].driftDistance=fabs(hit_vec[m]->get_drift_distance());
				
				if(1 <= detid && detid <= 30){
					if((_event_id-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+(hit_ctr[detid])*4*nhits > EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsChambers)
					cout << _event_id << " " << detid <<  " " << (_event_id-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+(hit_ctr[detid])+4*nhits << " " << EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsChambers << " " << hit_ctr[detid] << endl;
					
#ifdef DEBUG
					if(_event_id==debug::EvRef+firstevent)cout << "hit offsets " << detid << " " << (_event_id-firstevent)*datasizes::eventhitsize[0] << " " << evhitarrayoffset[detid] << " " << (_event_id-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid] << ": " << (_event_id-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+hit_ctr[detid] << " " << hit_vec[m]->get_element_id() << " " << (_event_id-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+hit_ctr[detid]+nhits << " " << wire_position[detid][hit_vec[m]->get_element_id()] << endl;
#endif					
					host_gEventHits.HitsChambersRawData[(_event_id-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+hit_ctr[detid]] = (float)hit_vec[m]->get_element_id();
					host_gEventHits.HitsChambersRawData[(_event_id-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+hit_ctr[detid]+nhits] = (float)wire_position[detid][hit_vec[m]->get_element_id()];
					host_gEventHits.HitsChambersRawData[(_event_id-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+hit_ctr[detid]+2*nhits] = (float)hit_vec[m]->get_tdc_time();
					host_gEventHits.HitsChambersRawData[(_event_id-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+hit_ctr[detid]+3*nhits] = (float)(1<<hit_vec[m]->is_in_time());
					host_gEventHits.HitsChambersRawData[(_event_id-firstevent)*datasizes::eventhitsize[0]+evhitarrayoffset[detid]+hit_ctr[detid]+4*nhits] = (float)fabs(hit_vec[m]->get_drift_distance());
					hit_ctr[detid]++;
				}

				if(31 <= detid && detid <= 46){
					if((_event_id-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+hit_ctr[detid]*4*nhits > EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes)
					cout << _event_id << " " << detid <<  " " << (_event_id-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+hit_ctr[detid]+4*nhits << " " << EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes << " " << hit_ctr[detid] << endl;

#ifdef DEBUG
					if(_event_id==debug::EvRef+firstevent)cout << "hit offsets " << detid << " " << (_event_id-firstevent)*datasizes::eventhitsize[1] << " " << evhitarrayoffset[detid] << " " << (_event_id-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid] << ": " << (_event_id-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+hit_ctr[detid] << " " << hit_vec[m]->get_element_id() << " " << (_event_id-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+hit_ctr[detid]+nhits << " el " << hit_vec[m]->get_element_id() << " " << wire_position[detid][hit_vec[m]->get_element_id()] << endl;
#endif
					
					host_gEventHits.HitsHodoRawData[(_event_id-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+hit_ctr[detid]] = (float)hit_vec[m]->get_element_id();
					host_gEventHits.HitsHodoRawData[(_event_id-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+hit_ctr[detid]+nhits] = (float) wire_position[detid][hit_vec[m]->get_element_id()];
					host_gEventHits.HitsHodoRawData[(_event_id-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+hit_ctr[detid]+2*nhits] = (float)hit_vec[m]->get_tdc_time();
					host_gEventHits.HitsHodoRawData[(_event_id-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+hit_ctr[detid]+3*nhits] = (float)(1<<hit_vec[m]->is_in_time());
					host_gEventHits.HitsHodoRawData[(_event_id-firstevent)*datasizes::eventhitsize[1]+evhitarrayoffset[detid]+hit_ctr[detid]+4*nhits] = (float)fabs(hit_vec[m]->get_drift_distance());
					hit_ctr[detid]++;
				}

				if(47 <= detid && detid <= 54){
					if( (_event_id-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+hit_ctr[detid]+4*nhits > EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsPropTubes)
					cout << _event_id << " " << detid <<  " " << (_event_id-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+hit_ctr[detid]*4*nhits << " " << EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsPropTubes << " " << hit_ctr[detid] << endl;

#ifdef DEBUG
					if(_event_id==debug::EvRef+firstevent)cout << "hit offsets " << detid << " " << (_event_id-firstevent)*datasizes::eventhitsize[2] << " " << evhitarrayoffset[detid] << " " << (_event_id-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid] << ": " << (_event_id-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+hit_ctr[detid] << " " << hit_vec[m]->get_element_id() << " " << (_event_id-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+hit_ctr[detid]+nhits << " " << wire_position[detid][hit_vec[m]->get_element_id()] << endl;
#endif
					host_gEventHits.HitsPropTubesRawData[(_event_id-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+hit_ctr[detid]] = (float)hit_vec[m]->get_element_id();
					host_gEventHits.HitsPropTubesRawData[(_event_id-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+hit_ctr[detid]+nhits] = (float) wire_position[detid][hit_vec[m]->get_element_id()];
					host_gEventHits.HitsPropTubesRawData[(_event_id-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+hit_ctr[detid]+2*nhits] = (float) hit_vec[m]->get_tdc_time();
					host_gEventHits.HitsPropTubesRawData[(_event_id-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+hit_ctr[detid]+3*nhits] = (float) (1<<hit_vec[m]->is_in_time());
					host_gEventHits.HitsPropTubesRawData[(_event_id-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+hit_ctr[detid]+4*nhits] = (float) fabs(hit_vec[m]->get_drift_distance());
					hit_ctr[detid]++;
#ifdef DEBUG					
					if(_event_id==debug::EvRef+firstevent)cout << host_gEventHits.HitsPropTubesRawData[(_event_id-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+hit_ctr[detid]] << " " << host_gEventHits.HitsPropTubesRawData[(_event_id-firstevent)*datasizes::eventhitsize[2]+evhitarrayoffset[detid]+hit_ctr[detid]+nhits] << endl;
#endif
				}

			}
			
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
	size_t NBytesAllTracks = sizeof(gEventTrackCollection);
	size_t NBytesHistsArrays = sizeof(gHistsArrays);
	size_t NBytesPlanes = sizeof(gPlane);
	
	cout << "Total size allocated on GPUs " << NBytesAllEvent+NBytesPlanes+NBytesAllHits+NBytesAllTracks+NBytesHistsArrays << endl;
	cout << " input events: " << NBytesAllEvent  
		<< "; raw hits: " << NBytesAllHits << "; tracks " << NBytesAllTracks
		<< "; histograms arrays " << NBytesHistsArrays << "; planes info: " << NBytesPlanes << endl;  
	
	gEvent* host_output_eR = (gEvent*)malloc(NBytesAllEvent);
	gEventHitCollections* host_output_gHits = (gEventHitCollections*)malloc(NBytesAllHits);
	gEventTrackCollection* host_output_gTracks = (gEventTrackCollection*)malloc(NBytesAllTracks);
	gHistsArrays* host_hists = (gHistsArrays*)malloc(NBytesHistsArrays);
	
	float varmin[NVars] = {-X0_MAX, -Y0_MAX, INVP_MIN, -TX_MAX, -TY_MAX, -VXY_MAX, -VXY_MAX, VZ_MIN, -PXY_MAX, -PXY_MAX, PZ_MIN};
	float varmax[NVars] = {+X0_MAX, +Y0_MAX, INVP_MAX, +TX_MAX, +TY_MAX, +VXY_MAX, +VXY_MAX, VZ_MAX, +PXY_MAX, +PXY_MAX, PZ_MAX};
	
#ifdef ROOTSAVE
	string wintitle[NVars] = {"x0", "y0", "invp", "tx", "ty", "vx", "vy", "vz", "px", "py", "pz"};
	TH1D* h1[NVars];
#endif	
	for(int k = 0; k<NVars; k++){	
#ifdef ROOTSAVE
		h1[k] = new TH1D(Form("h%d", k), wintitle[k].c_str(), Nbins_Hists, varmin[k], varmax[k]);
#endif
		host_hists->pts_hw[k] = (varmax[k]-varmin[k])*0.5f/Nbins_Hists;
		for(int l = 0; l<Nbins_Hists; l++){
			host_hists->xpts[k*Nbins_Hists+l] = varmin[k]+host_hists->pts_hw[k]*(2.f*l+1.f);
			if(k==2){
#ifdef DEBUG
				if(l==0)cout << host_hists->pts_hw[k] << " ";
				cout << k << " " << l << " " << k*Nbins_Hists+l << " " << host_hists->xpts[k*Nbins_Hists+l] << endl;
#endif
			}
			host_hists->values[k*Nbins_Hists+l] = 0;
		}
	}
	
	// declaring gEvent objects for the device (GPU) to use.
	gEvent* device_gEvent;
	gEventHitCollections* device_gHits;
	gEventTrackCollection* device_gTracks;
	gHistsArrays* device_gHistsArrays;
	gPlane* device_gPlane;
	
	// copy of data from host to device: evaluate operation time 
	// Allocating memory for GPU (pointer to allocated device ); check for errors in the process; stops the program if issues encountered
	gpuErrchk( cudaMalloc((void**)&device_gEvent, NBytesAllEvent));
	gpuErrchk( cudaMalloc((void**)&device_gHits, NBytesAllHits));
	gpuErrchk( cudaMalloc((void**)&device_gTracks, NBytesAllTracks));
	gpuErrchk( cudaMalloc((void**)&device_gHistsArrays, NBytesHistsArrays));
	//allocating the memory for the planes
	gpuErrchk( cudaMalloc((void**)&device_gPlane, NBytesPlanes));
	
	std::size_t free_bytes;
	std::size_t total_bytes;

	CUDA_CHECK_STATUS(cudaMemGetInfo(&free_bytes, &total_bytes));
    	cout << "Current memory foot print: " << free_bytes << " / " << total_bytes << endl;
	
	// cudaMemcpy(dst, src, count, kind): copies data between host and device:
	// dst: destination memory address; src: source memory address; count: size in bytes; kind: type of transfer
	gpuErrchk( cudaMemcpy(device_gPlane, &plane, NBytesPlanes, cudaMemcpyHostToDevice));
	gpuErrchk( cudaMemcpy(device_gHistsArrays, host_hists, NBytesHistsArrays, cudaMemcpyHostToDevice));
	gpuErrchk( cudaMemcpy(device_gEvent, &host_gEvent, NBytesAllEvent, cudaMemcpyHostToDevice));
	gpuErrchk( cudaMemcpy(device_gHits, &host_gEventHits, NBytesAllHits, cudaMemcpyHostToDevice));
	
#ifdef DEBUG
	gkernel_checkHistos<<<16, Nbins_Hists>>>(device_gHistsArrays, 2);
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	gpuErrchk( cudaMemcpy(host_hists, device_gHistsArrays, NBytesHistsArrays, cudaMemcpyDeviceToHost));
	cout << " x_hw " << host_hists->pts_hw[2] << " ";
	for(int k = 0; k<Nbins_Hists; k++)cout << "  k " << k << " idx " << 2*Nbins_Hists+k << " x " << host_hists->xpts[2*Nbins_Hists+k] << " data " << host_hists->values[2*Nbins_Hists+k]  << endl;
#endif
		
	auto cp3 = std::chrono::system_clock::now();

	auto cp_to_gpu = cp3-cp2;
	cout<<"Copy to GPU: "<<cp_to_gpu.count()/1000000000.<<endl;
		
	// now data is transfered in the device: kernel function for event reconstruction called;
	// note that the function call is made requesting a number of blocks and a number of threads per block
	// in practice we have as many threads total as number of events; 
	gkernel_eR<<<BLOCKS_NUM,8>>>(device_gHits, device_gEvent->HasTooManyHits);
	
	// check status of device and synchronize;
	size_t nEvents = EstnEvtMax;
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	auto cp4 = std::chrono::system_clock::now();
	auto gpu_er = cp4-cp3;
	cout<<"GPU: event reducing: "<<gpu_er.count()/1000000000.<<endl;
<<<<<<< HEAD

	//for test only
#ifdef DEBUG
#ifdef GLDISPLAY
	//since the transfer takes so little time, it can be considered to do it all on CPU, if needed...
	//runDisplay(argn, argv, device_gHistsArrays->xpts, device_gHistsArrays->values);

	gpuErrchk( cudaMemcpy(host_hists, device_gHistsArrays, NBytesHistsArrays, cudaMemcpyDeviceToHost));
	cudaFree(device_gHistsArrays);
	//TEST ONLY
	for(int k = 0; k<128; k++)host_hists->values[k] = 1000.*sin(k*3.141592653/180.);
	//
	//for(int k = 0; k<128; k++)cout << "  k " << k << " data " << host_hists->values[k] 
	//				//<< " " << host_hists->xpts[k] << " " << (int)1000*sin(host_hists->xpts[k]) 
	//				<< "  ";
	//cout << endl;
	runDisplay(argn, argv, host_hists->values);
#endif
#endif

	gKernel_XZ_tracking<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(
		device_gHits,
		device_gTracks,
		device_gPlane->z,
		device_gPlane->spacing,
#ifdef USE_DET_RESOL
		device_gPlane->resolution,
#endif
#ifdef DEBUG
		device_gEvent->EventID,
#endif
		device_gEvent->nTracklets,
		device_gEvent->HasTooManyHits);
=======
	
	gKernel_XZ_YZ_tracking<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gStraightTrackBuilder, device_gFitArrays, device_gPlane);
>>>>>>> e04fbf96b0c9a61c7bd2cf5732b68f3000b4ed2d
	
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

#ifdef DEBUG
	gKernel_check_tracks<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gTracks, device_gEvent->HasTooManyHits, debug::EvRef);
#endif

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	auto cp5 = std::chrono::system_clock::now();
	auto gpu_stx = cp5-cp4;
	cout<<"GPU: XZ straight tracking: "<<gpu_stx.count()/1000000000.<<endl;

	gKernel_YZ_tracking<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(
		device_gHits,
		device_gTracks,
		device_gPlane,
#ifdef DEBUG
		device_gEvent->EventID,
#endif
		device_gEvent->HasTooManyHits);

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

#ifdef DEBUG
	gKernel_check_tracks<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gTracks, device_gEvent->HasTooManyHits, debug::EvRef);
	
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
#endif
	
	gKernel_BackTrackCleaning<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gTracks, device_gEvent->HasTooManyHits);
		
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	auto cp6 = std::chrono::system_clock::now();
	auto gpu_sty = cp6-cp5;
	cout<<"GPU: YZ straight tracking: "<<gpu_sty.count()/1000000000.<<endl;
	
	gKernel_Global_tracking<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(
		device_gHits,
		device_gTracks,
		device_gPlane,
#ifdef DEBUG
		device_gEvent->EventID,
#endif
		device_gEvent->HasTooManyHits);

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

#ifdef DEBUG
	gKernel_check_tracks<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gTracks, device_gEvent->HasTooManyHits, debug::EvRef);

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
#endif

	gKernel_GlobalTrackCleaning<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(
		device_gHits,
		device_gTracks,
		device_gPlane->z,
#ifdef DEBUG
		device_gEvent->EventID,
#endif
		device_gEvent->HasTooManyHits);

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	auto cp7 = std::chrono::system_clock::now();
	auto gpu_gt = cp7-cp6;
	cout<<"GPU: global tracking: "<<gpu_gt.count()/1000000000.<<endl;


	gKernel_Vertexing<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(
		device_gTracks,
		device_gPlane->z,
#ifdef DEBUG
		device_gEvent->EventID,
#endif
		device_gEvent->HasTooManyHits);

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );


#ifdef DEBUG
	gKernel_check_tracks<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gTracks, device_gEvent->HasTooManyHits, debug::EvRef);
#endif

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	auto cp8 = std::chrono::system_clock::now();
	auto gpu_vtx = cp8-cp7;
	cout<<"GPU: Vertexing: "<<gpu_vtx.count()/1000000000.<<endl;

	
	gKernel_fill_display_histograms<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gTracks, device_gHistsArrays, device_gEvent->HasTooManyHits);
	
#ifdef ROOTSAVE
	gpuErrchk( cudaMemcpy(host_hists, device_gHistsArrays, NBytesHistsArrays, cudaMemcpyDeviceToHost));
	cudaFree(device_gHistsArrays);
	
	for(int k = 0; k<NVars; k++){
		for(int l = 0; l<Nbins_Hists; l++){
			h1[k]->SetBinContent(l+1, host_hists->values[k*Nbins_Hists+l]);
		}
	}
	
	TCanvas* C1 = new TCanvas("C1", "display OR", 1350, 800);
	C1->Divide(3, 2);
	for(int k = 0; k<5; k++){
		C1->cd(k+1);
		h1[k]->Draw("hist");
	}
	C1->SaveAs("DISPLAY.pdf");
#endif		
	//The display should come after the vertexing - but before the copy of the output back to the CPU.
#ifdef GLDISPLAY
	//since the transfer takes so little time, it can be considered to do it all on CPU, if needed...
	//runDisplay(argn, argv, device_gHistsArrays->xpts, device_gHistsArrays->values);

	gpuErrchk( cudaMemcpy(host_hists, device_gHistsArrays, NBytesHistsArrays, cudaMemcpyDeviceToHost));
	//TEST ONLY
	//for(int k = 0; k<128; k++)host_hists->values[k] = 1000.*sin(k*3.141592653/180.);
	//
#ifdef DEBUG
	for(int k = 0; k<Nbins_Hists; k++)cout << "  k " << k << " " << host_hists->xpts[2*Nbins_Hists+k] << " data " << host_hists->values[2*Nbins_Hists+k] << endl;
#endif
	cout << endl;
	runDisplay(argn, argv, host_hists->values);
	cudaFree(device_gHistsArrays);
#endif		
	//gKernel_GlobalTrack_KalmanFitting<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_output_TKL, device_gKalmanFitArrays, device_gPlane);

	auto cp9 = std::chrono::system_clock::now();
	auto gpu_dis = cp9-cp8;
	cout<<"GPU: Histograming and displaying: "<<gpu_dis.count()/1000000000.<<endl;

	// data transfer from device to host
	gpuErrchk( cudaMemcpy(host_output_eR, device_gEvent, NBytesAllEvent, cudaMemcpyDeviceToHost));
	cudaFree(device_gEvent);

	gpuErrchk( cudaMemcpy(host_output_gHits, device_gHits, NBytesAllHits, cudaMemcpyDeviceToHost));
	cudaFree(device_gHits);

	gpuErrchk( cudaMemcpy(host_output_gTracks, device_gTracks, NBytesAllTracks, cudaMemcpyDeviceToHost));
	cudaFree(device_gHits);

#ifndef ROOTSAVE
	gpuErrchk( cudaMemcpy(host_hists, device_gHistsArrays, NBytesHistsArrays, cudaMemcpyDeviceToHost));
	cudaFree(device_gHistsArrays);
#endif
	
	auto cp10 = std::chrono::system_clock::now();
	auto cp_to_cpu = cp10-cp9;
	cout<<"Copy back to CPU: "<<cp_to_cpu.count()/1000000000.<<endl;
	
	int nGood = EstnEvtMax;
	for(int i = 0; i<EstnEvtMax; i++)if(host_output_eR->HasTooManyHits[i])nGood--;
	cout << nGood << " events over " << EstnEvtMax << endl;
	
	unsigned nhits_total;
	unsigned int tkl_coll_offset;
	unsigned int array_thread_offset;
	unsigned int tklmult_idx;
	int ntkl;
	int nTracklets;
	int nhits_tkl;
	//ofstream out("OutputFile.txt");
	//Write in a file, 
	long tklctr = 0;
	long nEvtsTotal = 0;
	long nEvtsPass = 0;
	/*
	for(int n = 0; n<nEvtMax; n++){
		if(host_output_eR->nAH[n]==0)continue;
		nEvtsTotal++;
		if(host_output_eR->HasTooManyHits[n])continue;
		nEvtsPass++;

		nhits_total = 0;
		for(int k = 1; k<=nChamberPlanes; k++ )nhits_total+= host_output_gHits->NHitsChambers[n*nChamberPlanes+k-1];
		for(int k = nChamberPlanes+1; k<=nChamberPlanes+nHodoPlanes; k++ )nhits_total+= host_output_gHits->NHitsHodo[n*nHodoPlanes+k-nChamberPlanes-1];
		for(int k = nChamberPlanes+nHodoPlanes+1; k<nDetectors; k++ )nhits_total+= host_output_gHits->NHitsPropTubes[n*nPropPlanes+k-nChamberPlanes-nHodoPlanes-1];

		tkl_coll_offset = n*datasizes::TrackSizeMax*datasizes::NTracksParam;
		nTracklets = 0;
		for(int m = 0; m<THREADS_PER_BLOCK; m++){
			tklmult_idx = n*THREADS_PER_BLOCK+m;
			nTracklets+= host_output_gTracks->NTracks[tklmult_idx];
		}
		out<<n<<" "<< nhits_total <<" "<<nTracklets <<endl;
#ifdef DEBUG
		if(n==debug::EvRef)cout << n<<" "<< host_output_eR->nAH[n] <<" "<< nTracklets <<endl;
#endif
		tklctr+= nTracklets;
		
		for(int k = 1; k<=nChamberPlanes; k++ )out << host_output_gHits->NHitsChambers[n*nChamberPlanes+k-1] << " ";
		for(int k = nChamberPlanes+1; k<=nChamberPlanes+nHodoPlanes; k++ )out << host_output_gHits->NHitsHodo[n*nHodoPlanes+k-nChamberPlanes-1] << " ";
		for(int k = nChamberPlanes+nHodoPlanes+1; k<nDetectors; k++ )out << host_output_gHits->NHitsPropTubes[n*nPropPlanes+k-nChamberPlanes-nHodoPlanes-1] << " ";
		out << endl;
		for(int k = 1; k<=nChamberPlanes; k++ ){
			nhits = host_output_gHits->NHitsChambers[n*nChamberPlanes+k-1];
			for(int l = 0; l<nhits; l++){
				out << k << " " << host_output_gHits->HitsChambersRawData[n*datasizes::eventhitsize[0]+evhitarrayoffset[k]+l] << " " 
				    << host_output_gHits->HitsChambersRawData[n*datasizes::eventhitsize[0]+evhitarrayoffset[k]+l+4*nhits] << endl;
			}
		}
		for(int k = nChamberPlanes+1; k<=nChamberPlanes+nHodoPlanes; k++ ){
			nhits = host_output_gHits->NHitsHodo[n*nHodoPlanes+k-nChamberPlanes-1];
			for(int l = 0; l<nhits; l++){
				out << k << " " << host_output_gHits->HitsHodoRawData[n*datasizes::eventhitsize[1]+evhitarrayoffset[k]+l] << " " 
				    << host_output_gHits->HitsHodoRawData[n*datasizes::eventhitsize[1]+evhitarrayoffset[k]+l+4*nhits] << endl;
			}
		}
		for(int k = nChamberPlanes+nHodoPlanes+1; k<nDetectors; k++ ){
			nhits = host_output_gHits->NHitsPropTubes[n*nPropPlanes+k-nChamberPlanes-nHodoPlanes-1];
			for(int l = 0; l<nhits; l++){
				out << k << " " << host_output_gHits->HitsPropTubesRawData[n*datasizes::eventhitsize[2]+evhitarrayoffset[k]+l] << " " 
				    << host_output_gHits->HitsPropTubesRawData[n*datasizes::eventhitsize[2]+evhitarrayoffset[k]+l+4*nhits] << endl;
			}
		}
		
		for(int m = 0; m<THREADS_PER_BLOCK; m++){
			tklmult_idx = n*THREADS_PER_BLOCK+m;
			array_thread_offset = m*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
			nTracklets = host_output_gTracks->NTracks[tklmult_idx];
			for(int k = 0; k<nTracklets; k++ ){
				out << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam] << " " //stid
				    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+7] << " " //x0
				    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+8] << " " //y0
				    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+5] << " " //tx
				    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+6] << " " //ty
				    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+9] << " " //invP
				    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+2] << endl; //Nhits
				nhits_tkl = host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+2];
				for(int l = 0; l<nhits_tkl; l++){
					out << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+16+l] << " " //detid
					    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+34+l] << " " //elid (chan)
					    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+70+l] << " " //drift
					    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+52+l] << endl; //pos
				}
				out << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+106] << " " 
				    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+107] << " " 
				    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+108] << " " 
				    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+109] << " " 
				    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+110] << " " 
				    << host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+111] << endl;
			}
		}
	}
	*/
	//TFile* outFile = new TFile(outputFile.Data(), "RECREATE");
	ORoutput_tree* output = new ORoutput_tree(outputFile.Data());
	for(int n = 0; n < nEvtMax; ++n) {
		output->Clear();
		nEvtsTotal++;
		if(host_output_eR->nAH[n]==0)continue;
		if(host_output_eR->HasTooManyHits[n]){
#ifdef DEBUG
			cout << " event " << host_output_eR->EventID[n] << endl;
#endif
			continue;
		}
		nEvtsPass++;

		nhits_total = 0;
		//for(int k = 1; k<=nChamberPlanes; k++ )nhits_total+= host_output_gHits->NHitsChambers[n*nChamberPlanes+k-1];
		//for(int k = nChamberPlanes+1; k<=nChamberPlanes+nHodoPlanes; k++ )nhits_total+= host_output_gHits->NHitsHodo[n*nHodoPlanes+k-nChamberPlanes-1];
		//for(int k = nChamberPlanes+nHodoPlanes+1; k<nDetectors; k++ )nhits_total+= host_output_gHits->NHitsPropTubes[n*nPropPlanes+k-nChamberPlanes-nHodoPlanes-1];
		
		output->fEventID = host_output_eR->EventID[n];
				
		for(int k = 1; k<=nChamberPlanes; k++ ){
			nhits = host_output_gHits->NHitsChambers[n*nChamberPlanes+k-1];
			nhits_total+= nhits;
			for(int l = 0; l<nhits; l++){
				output->fHitDetID.push_back(k);
				output->fHitChan.push_back(host_output_gHits->HitsChambersRawData[n*datasizes::eventhitsize[0]+evhitarrayoffset[k]+l]);
				output->fHitPos.push_back(host_output_gHits->HitsChambersRawData[n*datasizes::eventhitsize[0]+evhitarrayoffset[k]+l+1*nhits]);
				output->fHitTDC.push_back(host_output_gHits->HitsChambersRawData[n*datasizes::eventhitsize[0]+evhitarrayoffset[k]+l+2*nhits]);
				output->fHitDrift.push_back(host_output_gHits->HitsChambersRawData[n*datasizes::eventhitsize[0]+evhitarrayoffset[k]+l+4*nhits]);
			}
		}
		for(int k = nChamberPlanes+1; k<=nChamberPlanes+nHodoPlanes; k++ ){
			nhits = host_output_gHits->NHitsHodo[n*nHodoPlanes+k-nChamberPlanes-1];
			nhits_total+= nhits;

			for(int l = 0; l<nhits; l++){
				output->fHitDetID.push_back(k);
				output->fHitChan.push_back(host_output_gHits->HitsHodoRawData[n*datasizes::eventhitsize[1]+evhitarrayoffset[k]+l]);
				output->fHitPos.push_back(host_output_gHits->HitsHodoRawData[n*datasizes::eventhitsize[1]+evhitarrayoffset[k]+l+1*nhits]);
				output->fHitTDC.push_back(host_output_gHits->HitsHodoRawData[n*datasizes::eventhitsize[1]+evhitarrayoffset[k]+l+2*nhits]);
				output->fHitDrift.push_back(host_output_gHits->HitsHodoRawData[n*datasizes::eventhitsize[1]+evhitarrayoffset[k]+l+4*nhits]);
			}
		}
		for(int k = nChamberPlanes+nHodoPlanes+1; k<nDetectors; k++ ){
			nhits = host_output_gHits->NHitsPropTubes[n*nPropPlanes+k-nChamberPlanes-nHodoPlanes-1];
			nhits_total+= nhits;
			for(int l = 0; l<nhits; l++){
				output->fHitDetID.push_back(k);
				output->fHitChan.push_back((int)host_output_gHits->HitsPropTubesRawData[n*datasizes::eventhitsize[2]+evhitarrayoffset[k]+l]);
				output->fHitPos.push_back(host_output_gHits->HitsPropTubesRawData[n*datasizes::eventhitsize[2]+evhitarrayoffset[k]+l+1*nhits]);
				output->fHitTDC.push_back(host_output_gHits->HitsPropTubesRawData[n*datasizes::eventhitsize[2]+evhitarrayoffset[k]+l+2*nhits]);
				output->fHitDrift.push_back(host_output_gHits->HitsPropTubesRawData[n*datasizes::eventhitsize[2]+evhitarrayoffset[k]+l+4*nhits]);
			}
		}
		output->fNHits = nhits_total;
		
		tkl_coll_offset = n*datasizes::TrackSizeMax*datasizes::NTracksParam;
		//nTracklets = 0;
		//for(int m = 0; m<THREADS_PER_BLOCK; m++){
		//	tklmult_idx = n*THREADS_PER_BLOCK+m;
		//	nTracklets+= host_output_gTracks->NTracks[tklmult_idx];
		//}
		
		nTracklets = 0;
		for(int m = 0; m<THREADS_PER_BLOCK; m++){
			tklmult_idx = n*THREADS_PER_BLOCK+m;
			array_thread_offset = m*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;

			ntkl = host_output_gTracks->NTracks[tklmult_idx];
			nTracklets+= ntkl;
			for(int k = 0; k<ntkl; k++ ){
				nhits_tkl = (int)host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+2];
				output->fTrackStID.push_back((int)host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam]);
				output->fTrackNHits.push_back(nhits_tkl);
				output->fTrackChi2.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+3]);
				output->fTrackX0.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+7]);
				output->fTrackY0.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+8]);
				output->fTrackTX.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+5]);
				output->fTrackTY.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+6]);
				output->fTrackInvP.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+9]);
				output->fTrackErrX0.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+12]);
				output->fTrackErrY0.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+13]);
				output->fTrackErrTX.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+10]);
				output->fTrackErrTY.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+11]);
				output->fTrackErrInvP.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+14]);
				output->fTrackCharge.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+15]);
				for(int l = 0; l<nhits_tkl; l++){
					output->fTrackHitsDetID.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+16+l]);
					output->fTrackHitsChan.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+34+l]);
					output->fTrackHitsPos.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+52+l]);
					output->fTrackHitsDrift.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+70+l]);
					output->fTrackHitsSign.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+88+l]);
				}
				output->fTrackVx.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+106]);
				output->fTrackVy.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+107]);
				output->fTrackVz.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+108]);
				output->fTrackPx.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+109]);
				output->fTrackPy.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+110]);
				output->fTrackPz.push_back(host_output_gTracks->TracksRawData[tkl_coll_offset+array_thread_offset+k*datasizes::NTracksParam+111]);
			}
		}
		output->fNTracks = nTracklets;
		tklctr+= nTracklets;
		output->FillTree();
	}
	//output->Write();
	output->Close();

	cout << tklctr << " total tracks reconstructed" << endl;
	cout << nEvtsPass << " evts with low enough number of hits on " << nEvtsTotal << " events total." << endl; 

	delete rawEvent;

	auto cp11 = std::chrono::system_clock::now();
	auto write_output = cp11-cp10;
	cout<<"Write Output: "<<write_output.count()/1000000000.<<endl;

	// printing the time required for all operations
	auto end = std::chrono::system_clock::now();
	auto overall = end - start;
	cout<<"Total time: "<<overall.count()/1000000000.<<endl;
	return 0;
}

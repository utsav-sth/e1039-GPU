#include "reconstruction_helper.cuh"
#include "trackanalyticalminimizer.cuh"

// kernel functions: 
// CUDA C++ extends C++ by allowing the programmer to define C++ functions, called kernels, that, when called, 
// are executed N times in parallel by N different CUDA threads, as opposed to only once like regular C++ functions. 

//called over 8192 blocks 8 threads each
__global__ void gkernel_eR(gEventHitCollections* hitcolls, bool* hastoomanyhits) {
	//const int index = threadIdx.x + blockIdx.x * blockDim.x;
	assert(blockIdx.x>=0 && blockIdx.x<EstnEvtMax);
	
	__shared__ unsigned short station_multiplicity[THREADS_PER_BLOCK][5];//0: D0; 1: D2; 2: D3p; 3: D3m; 4: proptubes
	
	for(int i = 0; i<5; i++)station_multiplicity[threadIdx.x][i] = 0;
	
	float hitarraycopy[600];
	int nvalues;
	
	short hitflag[200];
	int nhits, hitidx;
	//const int //thread: threadIdx.x;
	//Load the hit collections: 3 wire chamber planes, 2 hodoscope planes, and 1 prop planes.
	// calculate the collection offsets:
	//
	const unsigned int detid_chambers[3] = {
		geometry::eff_detid_chambers[threadIdx.x],
		geometry::eff_detid_chambers[threadIdx.x+8],
		geometry::eff_detid_chambers[threadIdx.x+16]
	};
	const unsigned int nhits_chambers[3] = {
		hitcolls->NHitsChambers[blockIdx.x*nChamberPlanes+detid_chambers[0]-1], 
		hitcolls->NHitsChambers[blockIdx.x*nChamberPlanes+detid_chambers[1]-1], 
		hitcolls->NHitsChambers[blockIdx.x*nChamberPlanes+detid_chambers[2]-1] 
	};
	const unsigned int offsets_hitcoll_chambers[3] = {
		blockIdx.x*datasizes::eventhitsize[0]+datasizes::NHitsParam*datasizes::NMaxHitsChambers*(detid_chambers[0]-1), 
		blockIdx.x*datasizes::eventhitsize[0]+datasizes::NHitsParam*datasizes::NMaxHitsChambers*(detid_chambers[1]-1),
		blockIdx.x*datasizes::eventhitsize[0]+datasizes::NHitsParam*datasizes::NMaxHitsChambers*(detid_chambers[2]-1)
	};
	
	const unsigned int detid_hodo[2] = {
		31+threadIdx.x,
		31+threadIdx.x+8
	};
	
	const unsigned int nhits_hodo[2] = {
		hitcolls->NHitsHodo[blockIdx.x*nHodoPlanes+threadIdx.x], 
		hitcolls->NHitsHodo[blockIdx.x*nHodoPlanes+(threadIdx.x+8)]
	};
	const unsigned int offsets_hitcoll_hodo[2] = {
		blockIdx.x*datasizes::eventhitsize[1]+datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes*threadIdx.x, 
		blockIdx.x*datasizes::eventhitsize[1]+datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes*(threadIdx.x+8)
	};
	
	const unsigned int detid_prop = 47+threadIdx.x;
	const unsigned int nhits_prop = hitcolls->NHitsPropTubes[blockIdx.x*nPropPlanes+threadIdx.x];
	const unsigned int offset_hitcoll_prop = blockIdx.x*datasizes::eventhitsize[2]+datasizes::NHitsParam*datasizes::NMaxHitsPropTubes*threadIdx.x;
	
	const gHits hitcoll_prop(hitcolls->HitsPropTubesRawData, nhits_prop, offset_hitcoll_prop);
	nhits = event_reduction(hitcoll_prop, hitflag, detid_prop, nhits_prop);
	
	hitcolls->NHitsPropTubes[blockIdx.x*nPropPlanes+threadIdx.x] = nhits;
	hitidx = 0;
	for(int k = 0; k<nhits_prop; k++){
		if(hitflag[k]>0){
			hitarraycopy[hitidx] = hitcoll_prop.chan(k);
			hitarraycopy[hitidx+nhits] = hitcoll_prop.pos(k);
			hitarraycopy[hitidx+nhits*2] = hitcoll_prop.tdc(k);
			hitarraycopy[hitidx+nhits*3] = hitcoll_prop.flag(k);
			hitarraycopy[hitidx+nhits*4] = hitcoll_prop.drift(k);
			hitidx++;
		}
	}
	nvalues = datasizes::NHitsParam*nhits;
	for(int k = 0; k<nvalues; k++)hitcolls->HitsPropTubesRawData[offset_hitcoll_prop+k] = hitarraycopy[k];
	
	station_multiplicity[threadIdx.x][4]+= nhits;
			
	int stid;
	for(int i = 0; i<3; i++){
		if(detid_chambers[i]<=12){
			stid = 0;
		}else{
			stid = (detid_chambers[i]-7)/6;
		}
		gHits hitcoll_chambers = gHits(hitcolls->HitsChambersRawData, nhits_chambers[i], offsets_hitcoll_chambers[i]);
		nhits = event_reduction(hitcoll_chambers, hitflag, detid_chambers[i], nhits_chambers[i]);
		
		hitcolls->NHitsChambers[blockIdx.x*nChamberPlanes+detid_chambers[i]-1] = nhits;
		hitidx = 0;
		for(int k = 0; k<nhits_chambers[i]; k++){
			if(hitflag[k]>0){
				hitarraycopy[hitidx] = hitcoll_chambers.chan(k);
				hitarraycopy[hitidx+nhits] = hitcoll_chambers.pos(k);
				hitarraycopy[hitidx+nhits*2] = hitcoll_chambers.tdc(k);
				hitarraycopy[hitidx+nhits*3] = hitcoll_chambers.flag(k);
				hitarraycopy[hitidx+nhits*4] = hitcoll_chambers.drift(k);
				hitidx++;
			}
			//if(blockIdx.x==debug::EvRef)printf("thread %d detid %d chan %1.1f, drift %1.4f, flag %d, hitidx %d hitarray_chan %1.1f \n", threadIdx.x, detid_chambers[i], hitcoll_chambers.chan(k), hitcoll_chambers.drift(k), hitflag[k], hitidx-1, hitarraycopy[hitidx-1]);
		}
		nvalues = datasizes::NHitsParam*nhits;
		for(int k = 0; k<nvalues; k++)hitcolls->HitsChambersRawData[offsets_hitcoll_chambers[i]+k] = hitarraycopy[k];

		station_multiplicity[threadIdx.x][stid]+=nhits;
	}
	for(int i = 0; i<2; i++){
		gHits hitcoll_hodo = gHits(hitcolls->HitsHodoRawData, nhits_hodo[i], offsets_hitcoll_hodo[i]);
		nhits = event_reduction(hitcoll_hodo, hitflag, detid_hodo[i], nhits_hodo[i]);
		
		hitcolls->NHitsHodo[blockIdx.x*nHodoPlanes+detid_hodo[i]-31] = nhits;
		hitidx = 0;
		for(int k = 0; k<nhits_hodo[i]; k++){
			if(hitflag[k]>0){
				hitarraycopy[hitidx] = hitcoll_hodo.chan(k);
				hitarraycopy[hitidx+nhits] = hitcoll_hodo.pos(k);
				hitarraycopy[hitidx+nhits*2] = hitcoll_hodo.tdc(k);
				hitarraycopy[hitidx+nhits*3] = hitcoll_hodo.flag(k);
				hitarraycopy[hitidx+nhits*4] = hitcoll_hodo.drift(k);
				hitidx++;
			}
		}
		nvalues = datasizes::NHitsParam*nhits;
		for(int k = 0; k<nvalues; k++)hitcolls->HitsHodoRawData[offsets_hitcoll_hodo[i]+k] = hitarraycopy[k];
	}
	__syncthreads();
	
	unsigned short station_mult[5] = {0, 0, 0, 0, 0};
	for(short j = 0; j<5; j++){
		for(short i = 0; i<THREADS_PER_BLOCK; i++){
			station_mult[j]+= station_multiplicity[i][j];
		}
	}
	
	hastoomanyhits[blockIdx.x] = false;
	if(station_mult[0]>selection::MaxD0Multiplicity)hastoomanyhits[blockIdx.x] = true;
	if(station_mult[1]>selection::MaxD2Multiplicity)hastoomanyhits[blockIdx.x] = true;
	if(station_mult[2]>selection::MaxD3Multiplicity)hastoomanyhits[blockIdx.x] = true;
	if(station_mult[3]>selection::MaxD3Multiplicity)hastoomanyhits[blockIdx.x] = true;
	//if(station_mult[4]>selection::MaxPropMultiplicity)hastoomanyhits[blockIdx.x] = true;
	
#ifdef DEBUG
	if(blockIdx.x!=debug::EvRef)hastoomanyhits[blockIdx.x] = true;

	if(blockIdx.x==debug::EvRef){
		for(int k = 0; k<3; k++){
		int nhits_chambers_ = hitcolls->NHitsChambers[blockIdx.x*nChamberPlanes+detid_chambers[k]-1];
		gHits hitcoll_chambers = gHits(hitcolls->HitsChambersRawData, nhits_chambers_, offsets_hitcoll_chambers[k]);
		printf(" det offset %d array offset %d nhits_chambers(%d) = %d \n", blockIdx.x*nChamberPlanes+detid_chambers[k]-1, offsets_hitcoll_chambers[k], detid_chambers[k], nhits_chambers_);
		for(int i = 0; i<nhits_chambers_; i++){
			printf("%d %d %1.0f %1.4f %1.4f %1.0f %1.4f\n", offsets_hitcoll_chambers[k], detid_chambers[k], hitcoll_chambers.chan(i), hitcoll_chambers.pos(i), hitcoll_chambers.tdc(i), hitcoll_chambers.flag(i), hitcoll_chambers.drift(i));
			//printf("%d %d %d %d %d %d %d %d %d %1.0f %1.4f\n", ev, offsets_hitcoll_chambers[k], detid_chambers[k], i, offsets_hitcoll_chambers[k]+i, offsets_hitcoll_chambers[k]+i+nhits_chambers[k], offsets_hitcoll_chambers[k]+i+nhits_chambers[k]*2, offsets_hitcoll_chambers[k]+i+nhits_chambers*3, offsets_hitcoll_chambers[k]+i+nhits_chambers[k]*4, hitcolls->HitsChambersRawData[offsets_hitcoll_chambers[k]+i], hitcolls->HitsChambersRawData[offsets_hitcoll_chambers[k]+i+nhits_chambers[2]]);
		}
		}
		int nhits_hodo_ = hitcolls->NHitsHodo[blockIdx.x*nHodoPlanes+threadIdx.x];
		gHits hitcoll_hodo = gHits(hitcolls->HitsHodoRawData, nhits_hodo_, offsets_hitcoll_hodo[0]);
		printf(" det offset %d array offset %d nhits_hodo(%d) = %d \n", blockIdx.x*nHodoPlanes+threadIdx.x, offsets_hitcoll_hodo[0], detid_hodo[0], nhits_hodo_);
		for(int i = 0; i<nhits_hodo_; i++){
			printf("%d %d %1.0f %1.4f %1.4f %1.0f %1.4f\n", offsets_hitcoll_hodo[0], 31+threadIdx.x, hitcoll_hodo.chan(i), hitcoll_hodo.pos(i), hitcoll_hodo.tdc(i), hitcoll_hodo.flag(i), hitcoll_hodo.drift(i));
			//printf("%d %d %d %d %d %d %d %d %d %1.0f %1.4f\n", ev, offsets_hitcoll_hodo[0], 31+threadIdx.x, i, offsets_hitcoll_hodo[0]+i, offsets_hitcoll_hodo[0]+i+nhits_hodo[0], offsets_hitcoll_hodo[0]+i+nhits_hodo[0]*2, offsets_hitcoll_hodo[0]+i+nhits_hodo[0]*3, offsets_hitcoll_hodo[0]+i+nhits_hodo[0]*4, hitcolls->HitsHodoRawData[offsets_hitcoll_hodo[0]+i], hitcolls->HitsHodoRawData[offsets_hitcoll_hodo[0]+i+nhits_hodo[0]]);

		}
		int nhits_prop_ = hitcolls->NHitsPropTubes[blockIdx.x*nPropPlanes+threadIdx.x];
		gHits hitcoll_prop2(hitcolls->HitsPropTubesRawData, nhits_prop_, offset_hitcoll_prop);
		printf(" det offset %d array offset %d nhits_prop(%d) = %d \n", blockIdx.x*nPropPlanes+threadIdx.x, offset_hitcoll_prop, 47+threadIdx.x, nhits_prop_);
		for(int i = 0; i<nhits_prop_; i++){
			printf("%d %d %1.0f %1.4f %1.4f %1.0f %1.4f\n", offset_hitcoll_prop, 47+threadIdx.x, hitcoll_prop2.chan(i), hitcoll_prop2.pos(i), hitcoll_prop2.tdc(i), hitcoll_prop2.flag(i), hitcoll_prop2.drift(i));
			//printf("%d %d %d %d %d %d %d %d %d %1.0f %1.4f\n", ev, offset_hitcoll_prop, 47+threadIdx.x, i, offset_hitcoll_prop+i, offset_hitcoll_prop+i+nhits_prop, offset_hitcoll_prop+i+nhits_prop*2, offset_hitcoll_prop+i+nhits_prop*3, offset_hitcoll_prop+i+nhits_prop*4, hitcolls->HitsPropTubesRawData[offset_hitcoll_prop+i], hitcolls->HitsPropTubesRawData[offset_hitcoll_prop+i+nhits_prop]);

		}
		if(threadIdx.x==0)for(int i = 0; i<5; i++)printf("thread %d station %d mult %d \n ", threadIdx.x, i, station_mult[i]);
	}
#endif
}



////////////////////////////////////
//
//          XZ TRACKING
//
////////////////////////////////////

__global__ void gKernel_XZ_tracking(
	gEventHitCollections* hitcolls,
	gEventTrackCollection* tklcoll,
	const float* z_array,
	const float* spacing_array,
#ifdef USE_DET_RESOL
	const float* res_array,
#endif
	int* nTracklets,
#ifdef DEBUG
	   const int* eventID,
#endif
	bool* hastoomanyhits)
//(gTrackingXZparams* parameters)			
{
	if(hastoomanyhits[blockIdx.x]){
#ifdef DEBUG
		if(threadIdx.x==0)printf("Evt %d discarded, too many hits\n", eventID[blockIdx.x]);
#endif
		return;
	}
	
	const int nbins_st2 = 7;//28/4
	//const int nbins_st3 = 28;//58/2
	//const int nbins_st3 = 14;
	const int nbins_st3 = 7;
	
	// 8 threads...
	// thread 0: bin0_st2 = 0/2*7, st3 = 3; thread 1: bin0_st2 = (1/2 = 0)*7, st3 = 1; thread 2: bin0_st2 = 2/2*7, st3 = 3; thread 3: bin0_st2 = (3/2 = 1)*7, st3 = 4; ...
	// const int bin0_st2 = (threadIdx.x/2)*nbins_st2;
	// const int bin0_st3 = 0;
	// 16 threads
	// const int bin0_st2 = (threadIdx.x/4)*nbins_st2;
	// const int bin0_st3 = (threadIdx.x%4/2)*nbins_st3;
	// 32 threads
	// thread 0: bin0_st2 = 0/8*7 = 0, st3 = 3; bin0_st3 = (0%8/2 = 0)*7 = 0; thread 1: bin0_st2 = (1/8 = 0)*7 = 0, st3 = 4, bin0_st3 = (1%8/2 = 0)*7 = 0; 
	// thread 2: bin0_st2 = 2/8*7 = 0, st3 = 3; bin0_st3 = (2%8/2 = 1)*7 = 7; thread 3: bin0_st2 = (3/8 = 0)*7 = 0, st3 = 4; bin0_st3 = (3%8/2 = 1)*7 = 7; 
	const int bin0_st2 = (threadIdx.x/8)*nbins_st2;
	const int bin0_st3 = (threadIdx.x%8/2)*nbins_st3;
	const int st3 = 3+threadIdx.x%2;//check d3p for even threads, d3m for odd threads...
	
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf(" thread %d bin0 st2 %d bin0 st3 %d \n", threadIdx.x, bin0_st2, bin0_st3);
#endif
	
	short hitflag1[100];
	short hitflag2[100];
	
	// As a starting point we will process 8 bins per thread: x acceptance divided by 4 * 2 bins for d3p and d3m 
	int nhitpairs_x2[nbins_st2];
	int nhitpairs_x3[nbins_st3];
        //pairs in station 2
        thrust::pair<int, int> hitpairs_x2[nbins_st2*geometry::MaxHitsProj[0]];
        //pairs in station 3
        thrust::pair<int, int> hitpairs_x3[nbins_st3*geometry::MaxHitsProj[0]];
	
	unsigned int offset_hitcoll;
	
	short stid, projid, detid, detoff;
	
	projid = 0;
	stid = 2;
	detoff = 1;
	short detid_list[4];
	
	//Get the required hit collections here!
	
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[0] = detid;
	int nhits_st2x;
	const gHits hits_st2x = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2x);
	const float z_st2x = z_array[detid];
#ifdef USE_DET_RESOL
	const float res_st2x = res_array[detid];
#else
	const float res_st2x = spacing_array[detid]*InvSqrt12;
#endif
	
	detid-= 1;
	detid_list[1] = detid;
	int nhits_st2xp;
	const gHits hits_st2xp = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2xp);
	const float z_st2xp = z_array[detid];
#ifdef USE_DET_RESOL
	const float res_st2xp = res_array[detid];
#else
	const float res_st2xp = spacing_array[detid]*InvSqrt12;
#endif
	
	make_hitpairs_in_station_bins(hits_st2x, nhits_st2x, hits_st2xp, nhits_st2xp, hitpairs_x2, nhitpairs_x2, bin0_st2, nbins_st2, hitflag1, hitflag2, stid, projid);

	stid = st3;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[2] = detid;
	int nhits_st3x;
	const gHits hits_st3x = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3x);
#ifdef DEBUG	
	if(blockIdx.x==debug::EvRef && st3==4)printf("block %d detid %d nhits %d, vs %d \n", blockIdx.x, detid, nhits_st3x, hitcolls->NHitsChambers[blockIdx.x*nChamberPlanes+detid-1]);
#endif
	const float z_st3x = z_array[detid];
#ifdef USE_DET_RESOL
	const float res_st3x = res_array[detid];
#else
	const float res_st3x = spacing_array[detid]*InvSqrt12;
#endif

	detid-= 1;
	int nhits_st3xp;
	detid_list[3] = detid;
	const gHits hits_st3xp = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3xp);
#ifdef DEBUG	
	if(blockIdx.x==debug::EvRef && st3==4)printf("block %d detid %d nhits %d, vs %d \n", blockIdx.x, detid, nhits_st3xp, hitcolls->NHitsChambers[blockIdx.x*nChamberPlanes+detid-1]);
#endif
	const float z_st3xp = z_array[detid];
#ifdef USE_DET_RESOL
	const float res_st3xp = res_array[detid];
#else
	const float res_st3xp = spacing_array[detid]*InvSqrt12;
#endif
	
	make_hitpairs_in_station_bins(hits_st3x, nhits_st3x, hits_st3xp, nhits_st3xp, hitpairs_x3, nhitpairs_x3, bin0_st3, nbins_st3, hitflag1, hitflag2, stid, projid);

#ifdef DEBUG	
	if(blockIdx.x==debug::EvRef && st3==4){
		printf("nhits %d %d \n", nhits_st3x, nhits_st3xp);
		for(int i = 0; i<nbins_st3; i++){
			printf("thread %d npairs %d \n", threadIdx.x, nhitpairs_x3[i]);
			for(int j = 0; j<nhitpairs_x3[i]; j++){
				printf("thread %d < %d , %d >\n", threadIdx.x, hitpairs_x3[i+nbins_st3*j].first, hitpairs_x3[i+nbins_st3*j].second);
				if(hitpairs_x3[i+nbins_st3*j].first>=0)printf("thread %d bin %d st3x, hit %d chan %1.0f pos %1.4f \n", threadIdx.x, i, hitpairs_x3[i+nbins_st3*j].first, hits_st3x.chan(hitpairs_x3[i+nbins_st3*j].first), hits_st3x.pos(hitpairs_x3[i+nbins_st3*j].first));
				if(hitpairs_x3[i+nbins_st3*j].second>=0)printf("thread %d bin %d st3xp, hit %d chan %1.0f pos %1.4f \n", threadIdx.x, i, hitpairs_x3[i+nbins_st3*j].second, hits_st3xp.chan(hitpairs_x3[i+nbins_st3*j].second), hits_st3xp.pos(hitpairs_x3[i+nbins_st3*j].second) );
			}
		}
	}
	
	if(debug::EvRef==blockIdx.x)printf("thread %d dets %d %d %d %d \n", threadIdx.x, detid_list[0], detid_list[1], detid_list[2], detid_list[3]);
#endif
	
	stid = 6-1;
	detid = geometry::detsuperid[stid][projid]*2;
	int nhits_p1x1;
	const gHits hits_p1x1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1x1);
	const float z_p1x1 = z_array[detid];

	detid-= 1;
	int nhits_p1x2;
	const gHits hits_p1x2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1x2);
	const float z_p1x2 = z_array[detid];
	
	stid = 7-1;
	detid = geometry::detsuperid[stid][projid]*2;
	int nhits_p2x1;
	const gHits hits_p2x1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2x1);
	const float z_p2x1 = z_array[detid];

	detid-= 1;
	int nhits_p2x2;
	const gHits hits_p2x2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2x2);
	const float z_p2x2 = z_array[detid];

	// hodoscope hits
	stid = 1;//2-1
	detid = geometry::hodoplanesx[stid][0];
	int nhits_h2x1;
	const gHits hits_h2x1 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h2x1);

	detid = geometry::hodoplanesx[stid][1];
	int nhits_h2x2;
	const gHits hits_h2x2 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h2x2);
	
	stid = 2;//3-1
	detid = geometry::hodoplanesx[stid][0];
	int nhits_h3x1;
	const gHits hits_h3x1 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h3x1);

	detid = geometry::hodoplanesx[stid][1];
	int nhits_h3x2;
	const gHits hits_h3x2 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h3x2);

	stid = 3;//4-1
	detid = geometry::hodoplanesx[stid][0];
	int nhits_h4x1;
	const gHits hits_h4x1 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h4x1);

	detid = geometry::hodoplanesx[stid][1];
	int nhits_h4x2;
	const gHits hits_h4x2 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h4x2);

	bool maskhodo[4];


	bool bin_overflows = false;
	for(int bin = 0; bin<nbins_st2; bin++){
		if(nhitpairs_x2[bin]>geometry::MaxHitsProj[0])bin_overflows = true;
#ifdef DEBUG
		if(nhitpairs_x2[bin])printf("evt %d bin %d nhits x2 = %d\n", eventID[blockIdx.x], bin, nhitpairs_x2[bin]);
#endif
	}
	for(int bin = 0; bin<nbins_st3; bin++){
		if(nhitpairs_x3[bin]>geometry::MaxHitsProj[0])bin_overflows = true;
#ifdef DEBUG
		if(nhitpairs_x3[bin])printf("evt %d bin %d nhits x3 = %d\n", eventID[blockIdx.x], bin, nhitpairs_x3[bin]);
#endif
	}
	if(bin_overflows){
		hastoomanyhits[blockIdx.x] = true;
		return;
	}
	
	//variables for 2D track fit
	short detID[4];
	float X[4];
	float errX[4];
	float Z[4];
	
	float A_[4];
	float Ainv_[4];
	float B_[2];
	float Par[2];
	float ParErr[2];
	float chi2;

	float x0, tx, errx0, errtx;
	
	//Arrays for other basic hit info
	short elID[4];
	float drift[4];
#ifdef FULLCODE
	float tdc[4];
#endif
	
	short sign[4];
	
	//small varaibles for prop tube matching
	short nprop;//, iprop, idet;
	float xExp, ipos;
	bool checknext;
	
	short nhits_x;
	short nhits_x2, nhits_x3;
	short i_x, i_x2, i_x3;
	int i_hit;
	
	// if parallel:
	// const short nbins_st3 = geometry::N_WCHitsBins[2];
	const short nbins_total = nbins_st2*nbins_st3;
	short bin2, bin3;
	int binId;
	
	int nx2, nx3;
	int ncomb_x;
		
	int n_goodxz;

	//tracklet container://limit the access to the large tracklet arrays.
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	const unsigned int tklmult_idx = blockIdx.x*THREADS_PER_BLOCK+threadIdx.x;

	tklcoll->NTracks[tklmult_idx] = 0;

	__shared__ unsigned int list_of_threads[THREADS_PER_BLOCK];
	__shared__ unsigned int ntkl_per_thread[THREADS_PER_BLOCK];
	for(int k = 0; k<THREADS_PER_BLOCK; k++){
		ntkl_per_thread[k] = 0;
		list_of_threads[k] = 0;
	}
	__shared__ short thread_min[THREADS_PER_BLOCK];//thread with lowest number of tracks: each thread
	__shared__ bool addtrack[THREADS_PER_BLOCK];//flag to mark threads requesting to add a new track in another thread.
	__shared__ bool threadbusy[THREADS_PER_BLOCK];//flag to indicate a busy thread - which means we need to find a new one
	unsigned int ntkl_min;
	unsigned int array_offset;
	unsigned int nthreads_busy;
	unsigned int nslots_available;
	
	for(short i = 0; i<nbins_total; i++){
		bin2 = i%nbins_st2;
		bin3 = (i-bin2)/nbins_st2;
		binId = threadIdx.x+THREADS_PER_BLOCK*i;
		//printf("bin %d %d\n", bin2, bin3);
		// if parallel:
		// bin3+=nbins_st3*i_thread;

		nx2 = nhitpairs_x2[bin2];
		nx3 = nhitpairs_x3[bin3];

		addtrack[threadIdx.x] = false;
		threadbusy[threadIdx.x] = false;

		if(nx2 == 0 || nx3==0) continue;
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef && st3==4)printf("bin %d %d, nx %d %d \n", bin2, bin3, nx2, nx3);
#endif		
		//for(int k = 0; k<THREADS_PER_BLOCK; k++)ntkl_per_bin[k][i] = 0;

		// evaluating the number of combinations;
		ncomb_x = nx2*nx3;
		
		n_goodxz = 0;
		
		for(i_x = 0; i_x<ncomb_x; i_x++){
			i_x2 = i_x%nx2;
			i_x3 = (i_x-i_x2)/nx2;
			
			for(int k = 0; k<THREADS_PER_BLOCK; k++){
				list_of_threads[k] = 0;
				addtrack[k] = false;
				threadbusy[k] = false;
			}
			//addtrack[threadIdx.x] = false;
			//threadbusy[threadIdx.x] = false;
			
			nhits_x = 0;
			
			if(hitpairs_x2[bin2+nbins_st2*i_x2].first>=0){
				detID[nhits_x] = detid_list[0];
				i_hit = hitpairs_x2[bin2+nbins_st2*i_x2].first;
				X[nhits_x] = hits_st2x.pos(i_hit);
				errX[nhits_x] = res_st2x;
				Z[nhits_x] = z_st2x;
				elID[nhits_x] = (short)hits_st2x.chan(i_hit);
#ifdef FULLCODE
				tdc[nhits_x] = hits_st2x.tdc(i_hit);
#endif
				drift[nhits_x] = hits_st2x.drift(i_hit);
				sign[nhits_x] = 0;
				nhits_x++;
			}
			if(hitpairs_x2[bin2+nbins_st2*i_x2].second>=0){
				detID[nhits_x] = detid_list[1];
				i_hit = hitpairs_x2[bin2+nbins_st2*i_x2].second;
				X[nhits_x] = hits_st2xp.pos(i_hit);
				errX[nhits_x] = res_st2xp;
				Z[nhits_x] = z_st2xp;
				elID[nhits_x] = (short)hits_st2xp.chan(i_hit);
#ifdef FULLCODE
				tdc[nhits_x] = hits_st2xp.tdc(i_hit);
#endif
				drift[nhits_x] = hits_st2xp.drift(i_hit);
				sign[nhits_x] = 0;
				nhits_x++;
			}
			
			nhits_x2 = nhits_x;
			if(nhits_x2==0) continue;

#ifdef DEBUG
			if(blockIdx.x==debug::EvRef && st3==4) printf("thread %d in loop (i_x2 %d i_x3 %d): < %d , %d >\n", threadIdx.x, i_x2, i_x3, hitpairs_x3[bin3+nbins_st3*i_x3].first, hitpairs_x3[bin3+nbins_st3*i_x3].second);
#endif
			
			if(hitpairs_x3[bin3+nbins_st3*i_x3].first>=0){
				detID[nhits_x] = detid_list[2];
				i_hit = hitpairs_x3[bin3+nbins_st3*i_x3].first;
				X[nhits_x] = hits_st3x.pos(i_hit);
				errX[nhits_x] = res_st3x;
				Z[nhits_x] = z_st3x;
				elID[nhits_x] = (short)hits_st3x.chan(i_hit);
#ifdef FULLCODE
				tdc[nhits_x] = hits_st3x.tdc(i_hit);
#endif
				drift[nhits_x] = hits_st3x.drift(i_hit);
				sign[nhits_x] = 0;
				nhits_x++;
			}
			if(hitpairs_x3[bin3+nbins_st3*i_x3].second>=0){
				detID[nhits_x] = detid_list[3];
				i_hit = hitpairs_x3[bin3+nbins_st3*i_x3].second;
				X[nhits_x] = hits_st3xp.pos(i_hit);
				errX[nhits_x] = res_st3xp;
				Z[nhits_x] = z_st3xp;
				elID[nhits_x] = (short)hits_st3xp.chan(i_hit);
#ifdef FULLCODE
				tdc[nhits_x] = hits_st3xp.tdc(i_hit);
#endif
				drift[nhits_x] = hits_st3xp.drift(i_hit);
				sign[nhits_x] = 0;
				nhits_x++;
			}
			
			nhits_x3 = nhits_x-nhits_x2;
			if(nhits_x3==0) continue;

			
			fit_2D_track(nhits_x, X, Z, errX, A_, Ainv_, B_, Par, ParErr, chi2);
			x0 = Par[0];
			tx = Par[1];
			errx0 = ParErr[0];
			errtx = ParErr[1];

			if(fabs(x0)>X0_MAX+2*errx0 || fabs(tx)>TX_MAX+2*errtx)continue;
			
			resolve_leftright_xhits(x0, tx, errx0, errtx, nhits_x, detID, X, drift, sign, z_array, 150.);
			resolve_single_leftright_xhits(x0, tx, nhits_x, detID, X, sign, z_array);
			
			fit_2D_track_drift(nhits_x, X, drift, sign, Z, errX, A_, Ainv_, B_, Par, ParErr, chi2);
			x0 = Par[0];
			tx = Par[1];
			errx0 = ParErr[0];
			errtx = ParErr[1];

			if(fabs(x0)>X0_MAX || fabs(tx)>TX_MAX)continue;
			//if(fabs(x0+tx*geometry::Z_KMAG_BEND)>geometry::X_KMAG_BEND)continue;
			//the following does not figure in the ktracker selection... yet it is fairly effective
			//if(fabs(x0)<10.0)continue;
			//if(x0<0 && tx<-x0*0.0008)continue;
			//if(x0>0 && tx>-x0*0.0008)continue;

			
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("evt %d x0 = %1.4f, tx = %1.4f \n", blockIdx.x, x0, tx);
#endif
			n_goodxz++;
			
			//prop matching
			nprop = 0;
			checknext = true;

			//we want at least one hit in p1x1 or p1x2...
			for(int n = 0; n<nhits_p1x1; n++){
				ipos = hits_p1x1.pos(n);
				xExp = tx*z_p1x1+x0;
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("evt %d p1x1, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
				if(fabs(ipos-xExp)<5.08f){
					nprop++;
					checknext = false;
					break;
				}
			}
			if(checknext){
				for(int n = 0; n<nhits_p1x2; n++){
					ipos = hits_p1x2.pos(n);
					xExp = tx*z_p1x2+x0;
#ifdef DEBUG
					if(blockIdx.x<=debug::EvRef)printf("evt %d p1x2, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
					if(fabs(ipos-xExp)<5.08f){
						nprop++;
//						checknext = false;
						break;
					}
				}
			}
			//and at least one hit in p2x1 or p2x2
			checknext = true;
			
			if(checknext){
				for(int n = 0; n<nhits_p2x1; n++){
					ipos = hits_p2x1.pos(n);
					xExp = tx*z_p2x1+x0;
#ifdef DEBUG
					if(blockIdx.x<=debug::EvRef)printf("evt %d p2x1, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
					if(fabs(ipos-xExp)<7.62f){
						nprop++;
						checknext = false;
						break;
					}
				}
			}
			if(checknext){
				for(int n = 0; n<nhits_p2x2; n++){
					ipos = hits_p2x2.pos(n);
					xExp = tx*z_p2x2+x0;
#ifdef DEBUG
					if(blockIdx.x<=debug::EvRef)printf("evt %d p2x2, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
					if(fabs(ipos-xExp)<7.62f){
						nprop++;
//						checknext = false;
						break;
					}
				}
			}
			//
			if(nprop<selection::NpropXYhitsMin)continue;
			
			//hodoscope
			stid = 1;//2-1

			maskhodo[stid] = 0;
			detid = geometry::hodoplanesx[stid][0];
			maskhodo[stid] = match_track_XZ_to_hodo(stid, detid, nhits_h2x1, hits_h2x1, x0, tx, errx0, errtx, z_array);
			if(!maskhodo[stid]){
				detid = geometry::hodoplanesx[stid][1];
				maskhodo[stid] = maskhodo[stid] || match_track_XZ_to_hodo(stid, detid, nhits_h2x2, hits_h2x2, x0, tx, errx0, errtx, z_array);
			}
			if(!maskhodo[stid])continue;
			
			stid = 2;//3-1
			maskhodo[stid] = 0;
			detid = geometry::hodoplanesx[stid][0];
			maskhodo[stid] = match_track_XZ_to_hodo(stid, detid, nhits_h3x1, hits_h3x1, x0, tx, errx0, errtx, z_array);
			if(!maskhodo[stid]){
				detid = geometry::hodoplanesx[stid][1];
				maskhodo[stid] = maskhodo[stid] || match_track_XZ_to_hodo(stid, detid, nhits_h3x2, hits_h3x2, x0, tx, errx0, errtx, z_array);
			}
			if(!maskhodo[stid])continue;

			stid = 3;//4-1
			maskhodo[stid] = 0;
			detid = geometry::hodoplanesx[stid][0];
			maskhodo[stid] = match_track_XZ_to_hodo(stid, detid, nhits_h4x1, hits_h4x1, x0, tx, errx0, errtx, z_array);
			if(!maskhodo[stid]){
				detid = geometry::hodoplanesx[stid][1];
				maskhodo[stid] = maskhodo[stid] || match_track_XZ_to_hodo(stid, detid, nhits_h4x2, hits_h4x2, x0, tx, errx0, errtx, z_array);
			}
			if(!maskhodo[stid])continue;
			
			//we can probably afford to spare time for synchronization here since XZ is extremely fast!
			addtrack[threadIdx.x] = true;

#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("thread %d wants to add a track!\n", threadIdx.x);
#endif
			ntkl_min = 100000;
			thread_min[threadIdx.x] = -1;
			nthreads_busy = 0;
			nslots_available = datasizes::TrackSizeMax;
			__syncthreads();
			for(int k = 0; k<THREADS_PER_BLOCK; k++){
				//list_of_threads[k] = 0;
				//looking for thread with least number of tracks (thread_min)
				//look only for threads with same "parity"
				if(ntkl_min>ntkl_per_thread[k] && st3==3+k%2){
					ntkl_min = ntkl_per_thread[k];
					thread_min[threadIdx.x] = k;
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("k %d st3 %d 3+kmod2 %d ntkl_min %d ntkl_per_thread %d thread_min(%d) = %d \n", k, st3, 3+k%2, ntkl_min, ntkl_per_thread[k], threadIdx.x, thread_min[threadIdx.x]);
#endif
					//threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
				}
				nslots_available-= ntkl_per_thread[k];
				// here we have a first "pick" for alternate thread.
				// listing the threads which are "busy" with a "addtrack" request here.
				if(addtrack[k]){
					list_of_threads[nthreads_busy] = k;
					nthreads_busy++;
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("thread? %d, nthreads_busy %d, list_of_threads = %d \n", k, nthreads_busy, list_of_threads[nthreads_busy-1]);
#endif
				}

			}
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("thread_min(%d) = %d \n", threadIdx.x, thread_min[threadIdx.x]);
			
			if(blockIdx.x==debug::EvRef)if( thread_min[threadIdx.x]%2 != threadIdx.x%2 )printf("!!! thread_min(%d) = %d \n", threadIdx.x, thread_min[threadIdx.x]);
			if(blockIdx.x==debug::EvRef){
				printf("number of slots available %d\n", nslots_available);
				for(int m = 0; m<nthreads_busy; m++)printf(" m %d, list_of_threads[m] = %d \n", m, list_of_threads[m]);
			}
#endif
			//if the number of tracks to add is greater than the number of slots available, we stop - the event is too full. 
			if(nslots_available<nthreads_busy){
				printf("WARNING: Block %d cannot store anymore tracks! ending the event with a high multiplicity flag! \n", blockIdx.x);
				hastoomanyhits[blockIdx.x] = true;
				break;
			}
			__syncthreads();
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("number of threads requesting to add a track (thread %d): %d  \n", threadIdx.x, nthreads_busy);
#endif
			//first, assign a unique "thread min" for each full thread so that they don't step on each other...
			// we have the following info:
			//   the list of threads that want to add a track in another thread,
			//   the number of tracks in each thread
			//   the "thread_min" i.e. the first alternate thread candidate for all threads. 
			threadbusy[thread_min[list_of_threads[0]]] = true;
			for(int l = 1; l<nthreads_busy; l++){
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("(thread %d) list_of_threads %d thread_min %d thread_min[l] %d =? thread_min[l-1] %d \n", threadIdx.x, l, list_of_threads[l], thread_min[list_of_threads[l]], thread_min[list_of_threads[l-1]]);
#endif
				if(thread_min[list_of_threads[l]]==thread_min[list_of_threads[l-1]] || // if the current "alternate" thread is the same as the previous one 
					thread_min[list_of_threads[l]]==thread_min[list_of_threads[0]] || // or if it's the same as the first one
					st3!=3+thread_min[list_of_threads[l]]%2 ){ // or if the thread_min is not compatible with the actual thread 
					threadbusy[thread_min[list_of_threads[l]]] = true; // the thread is marked as busy and cannot be used
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("(thread %d) thread %d busy\n", threadIdx.x, thread_min[list_of_threads[l]]);
#endif
				}
				
				//if the thread is busy, find another thread!
				if(threadbusy[thread_min[threadIdx.x]]){
					ntkl_min = 100000;
					for(int k = 0; k<THREADS_PER_BLOCK; k++){
						//we want another thread that is the least populated possible, that is not already busy, and that is compatible with the actual thread...
						if(ntkl_min>ntkl_per_thread[k] && !threadbusy[k] && thread_min[list_of_threads[l]]%2==k%2){ 
							ntkl_min = ntkl_per_thread[k];
							thread_min[list_of_threads[l]] = k;
							//thread_min[threadIdx.x] = k;
							threadbusy[k] = true;
#ifdef DEBUG
							if(blockIdx.x==debug::EvRef)printf("k %d ntkl_min %d list_of_threads %d thread_min(listofthreads[l]) = %d, thread_min(thread)\n", k, ntkl_min, list_of_threads[l], thread_min[list_of_threads[l]], thread_min[threadIdx.x]);
#endif
						}
					}
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("thread %d ntkl_min %d thread_min(thread) = %d\n", threadIdx.x, ntkl_min, thread_min[threadIdx.x]);
#endif
				}
			}
			__syncthreads();
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("alt thread chosen for thread %d: %d \n", threadIdx.x, thread_min[threadIdx.x]);
#endif
			array_offset = thread_min[threadIdx.x]*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;

#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("actual thread %d store thread %d offset %d stid %d local bin %d, bin0 st2 %d, bin0 st3 %d, st3 %d, ntkl_per_thread %d \n",  threadIdx.x, thread_min[threadIdx.x], tkl_coll_offset+array_offset, binId, i, bin0_st2, bin0_st3, st3, ntkl_per_thread[thread_min[threadIdx.x]]);

			if( threadIdx.x%2 != thread_min[threadIdx.x]%2 )printf(" !!! actual thread %d  store thread %d, st3 %d st_thread %d\n", threadIdx.x, thread_min[threadIdx.x], st3, 3+thread_min[threadIdx.x]%2);
			
			for(int kk = 0; kk<nthreads_busy; kk++){
				for(int ll = 0; ll<nthreads_busy; ll++){
					if(kk==ll)continue;
					if(thread_min[list_of_threads[kk]]==thread_min[list_of_threads[ll]] && ntkl_per_thread[list_of_threads[kk]]==ntkl_per_thread[list_of_threads[ll]]){
						printf("!!! alt thread chosen for thread %d: %d (n_tkl = %d) == alt thread chosen for thread %d: %d (n_tkl = %d) \n", list_of_threads[kk], thread_min[list_of_threads[kk]], ntkl_per_thread[list_of_threads[kk]], list_of_threads[ll], thread_min[list_of_threads[ll]], ntkl_per_thread[list_of_threads[ll]]);
					}
				}
			}
#endif
			
			tklcoll->setStationID(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], (float)binId);
			tklcoll->setThreadID(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], (float)threadIdx.x);
			tklcoll->setnHits(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], (float)nhits_x);
			tklcoll->setTx(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], tx);
			tklcoll->setX0(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], x0);
			tklcoll->setErrTx(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], errtx);
			tklcoll->setErrX0(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], errx0);
			
			for(int n = 0; n<nhits_x; n++){
				tklcoll->setHitDetID(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, detID[n]);
				tklcoll->setHitChan(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, elID[n]);
				tklcoll->setHitPos(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, X[n]);
				tklcoll->setHitDrift(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, drift[n]);
				tklcoll->setHitSign(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, sign[n]);
#ifdef FULLCODE					
				tklcoll->setHitTDC(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, tdc[n]);
				tklcoll->setHitResidual(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, 0.0);
#endif
			}
			ntkl_per_thread[thread_min[threadIdx.x]]++;
			tklcoll->NTracks[blockIdx.x*THREADS_PER_BLOCK+thread_min[threadIdx.x]]++;
			__syncthreads();
		}// end loop on hits
	}//end loop on bins
	//evaluate number of tracklets

	__syncthreads();
	
	int N_tracklets = 0;

#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf("block %d thread %d tracklets per thread: %d \n", blockIdx.x, threadIdx.x, ntkl_per_thread[threadIdx.x]);
#endif
	//__shared__ unsigned int array_thread_offset[THREADS_PER_BLOCK];
	for(int k = 0; k<THREADS_PER_BLOCK; k++){
		if(tklcoll->NTracks[tklmult_idx]!=ntkl_per_thread[threadIdx.x])printf("block %d, thread %d n tracks for thread ? %d %d \n", blockIdx.x, threadIdx.x, tklcoll->NTracks[tklmult_idx], ntkl_per_thread[threadIdx.x]);
		N_tracklets+= ntkl_per_thread[k];
	}
	if(ntkl_per_thread[threadIdx.x]>datasizes::TrackSizeMax/THREADS_PER_BLOCK){
#ifdef DEBUG
		printf("block %d thread %d tracklets per thread: %d \n", blockIdx.x, threadIdx.x, ntkl_per_thread[threadIdx.x]);
#endif
		hastoomanyhits[blockIdx.x] = true;
	}
	nTracklets[blockIdx.x] = N_tracklets;
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf(" Ntracklets %d \n", N_tracklets);
#endif
	//at the end like that it's probably fine...
	if(N_tracklets>=datasizes::TrackSizeMax){
		printf("block %d thread %d tracklets total %d \n", blockIdx.x, threadIdx.x, N_tracklets);
		hastoomanyhits[blockIdx.x] = true;
		tklcoll->NTracks[blockIdx.x] = blockIdx.x*THREADS_PER_BLOCK+threadIdx.x;
	}
#ifdef DEBUG
	else{
		for(int m = 0; m<ntkl_per_thread[threadIdx.x]; m++){
			if(blockIdx.x==debug::EvRef)printf(" thread: %d, m %d, offset %d stid/bin?  %1.0f thread %1.0f nhits %1.0f x0 %1.4f tx %1.4f \n", 
				threadIdx.x, m, tkl_coll_offset+array_thread_offset, 
				tklcoll->TracksRawData[tkl_coll_offset+array_thread_offset+m*datasizes::NTracksParam],  
				tklcoll->TracksRawData[tkl_coll_offset+array_thread_offset+m*datasizes::NTracksParam+1],  
				tklcoll->TracksRawData[tkl_coll_offset+array_thread_offset+m*datasizes::NTracksParam+2], 
				tklcoll->TracksRawData[tkl_coll_offset+array_thread_offset+m*datasizes::NTracksParam+7], 
				tklcoll->TracksRawData[tkl_coll_offset+array_thread_offset+m*datasizes::NTracksParam+5]
				);
		}
	}
#endif

}


////////////////////////////////////
//
//          YZ TRACKING
//
////////////////////////////////////

__global__ void gKernel_YZ_tracking(
	gEventHitCollections* hitcolls,
	gEventTrackCollection* tklcoll,
	const gPlane* planes,
#ifdef DEBUG
	int* eventID,
#endif
	bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x]){
#ifdef DEBUG
		if(threadIdx.x==0)printf("Evt %d discarded, too many hits\n", eventID[blockIdx.x]);
#endif
		return;
	}
	const int nbins_st2 = 7;
	//const int nbins_st3 = 28;
	//const int nbins_st3 = 14;
	const int nbins_st3 = 7;
	
	const int st3 = 3+threadIdx.x%2;//check d3p for even threads, d3m for odd threads...
	
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf(" thread %d bin0 st2 %d bin0 st3 %d \n", threadIdx.x, bin0_st2, bin0_st3);
#endif

	short hitidx1[100];
	short hitidx2[100];
	short hitflag1[100];
	short hitflag2[100];
	
        //pairs in station 2
        thrust::pair<int, int> hitpairs_u2[100];
        //pairs in station 3
        thrust::pair<int, int> hitpairs_u3[100];

        //pairs in station 2
        thrust::pair<int, int> hitpairs_v2[100];
        //pairs in station 3
        thrust::pair<int, int> hitpairs_v3[100];
	
	
	unsigned int offset_hitcoll;
	
	short stid, projid, detid, detoff, chan;
	
	projid = 1;
	stid = 2;
	detoff = 1;
	short detid_list[8];
	float det_spacing[8];
	
	//Get the required hit collections here!
	
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[0] = detid;
	det_spacing[0] = planes->spacing[detid];
	int nhits_st2u;
	const gHits hits_st2u = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2u);
	const float z_st2u = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st2u = planes->resolution[detid];
#else
	const float res_st2u = planes->spacing[detid]*InvSqrt12;
#endif
	
	detid-= 1;
	detid_list[1] = detid;
	det_spacing[1] = planes->spacing[detid];
	int nhits_st2up;
	const gHits hits_st2up = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2up);
	const float z_st2up = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st2up = planes->resolution[detid];
#else
	const float res_st2up = planes->spacing[detid]*InvSqrt12;
#endif
	
	projid = 2;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[2] = detid;
	det_spacing[2] = planes->spacing[detid];
	int nhits_st2v;
	const gHits hits_st2v = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2v);
	const float z_st2v = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st2v = planes->resolution[detid];
#else
	const float res_st2v = planes->spacing[detid]*InvSqrt12;
#endif
	
	detid-= 1;
	detid_list[3] = detid;
	det_spacing[3] = planes->spacing[detid];
	int nhits_st2vp;
	const gHits hits_st2vp = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2vp);
	const float z_st2vp = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st2vp = planes->resolution[detid];
#else
	const float res_st2vp = planes->spacing[detid]*InvSqrt12;
#endif
		
	projid = 1;
	stid = st3;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[4] = detid;
	det_spacing[4] = planes->spacing[detid];
	int nhits_st3u;
	const gHits hits_st3u = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3u);
	const float z_st3u = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st3u = planes->resolution[detid];
#else
	const float res_st3u = planes->spacing[detid]*InvSqrt12;
#endif

	detid-= 1;
	int nhits_st3up;
	detid_list[5] = detid;
	det_spacing[5] = planes->spacing[detid];
	const gHits hits_st3up = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3up);
	const float z_st3up = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st3up = planes->resolution[detid];
#else
	const float res_st3up = planes->spacing[detid]*InvSqrt12;
#endif
	
	projid = 2;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[6] = detid;
	det_spacing[6] = planes->spacing[detid];
	int nhits_st3v;
	const gHits hits_st3v = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3v);
	const float z_st3v = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st3v = planes->resolution[detid];
#else
	const float res_st3v = planes->spacing[detid]*InvSqrt12;
#endif
	
	detid-= 1;
	int nhits_st3vp;
	detid_list[7] = detid;
	det_spacing[7] = planes->spacing[detid];
	const gHits hits_st3vp = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3vp);
	const float z_st3vp = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st3vp = planes->resolution[detid];
#else
	const float res_st3vp = planes->spacing[detid]*InvSqrt12;
#endif

#ifdef PROP_Y_MATCH	
	int nprop;
	float yExp, ipos;
	bool checknext;

	// proportional tube hits
	projid = 1;
	stid = 6-1;
	detid = geometry::detsuperid[stid][projid]*2;
	int nhits_p1y1;
	const gHits hits_p1y1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1y1);
	const float z_p1y1 = planes->z[detid];

	detid-= 1;
	int nhits_p1y2;
	const gHits hits_p1y2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1y2);
	const float z_p1y2 = planes->z[detid];
	
	stid = 7-1;
	detid = geometry::detsuperid[stid][projid]*2;
	int nhits_p2y1;
	const gHits hits_p2y1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2y1);
	const float z_p2y1 = planes->z[detid];

	detid-= 1;
	int nhits_p2y2;
	const gHits hits_p2y2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2y2);
	const float z_p2y2 = planes->z[detid];
#endif
	
	// hodoscope hits
	stid = 1;//2-1
	detid = geometry::hodoplanesx[stid][0];
	int nhits_h2x1;
	const gHits hits_h2x1 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h2x1);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h2x1 %d \n", nhits_h2x1);
	for(int l = 0; l<nhits_h2x1; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h2x1.chan(l), hits_h2x1.pos(l));
	}
#endif


	detid = geometry::hodoplanesx[stid][1];
	int nhits_h2x2;
	const gHits hits_h2x2 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h2x2);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h2x2 %d \n", nhits_h2x2);
	for(int l = 0; l<nhits_h2x2; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h2x2.chan(l), hits_h2x2.pos(l));
	}
#endif

#ifdef HODO_Y_MATCH
	detid = geometry::hodoplanesy[stid][0];
	int nhits_h2y1;
	const gHits hits_h2y1 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h2y1);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h2y1 %d \n", nhits_h2y1);
	for(int l = 0; l<nhits_h2y1; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h2y1.chan(l), hits_h2y1.pos(l));
	}
#endif

	detid = geometry::hodoplanesy[stid][1];
	int nhits_h2y2;
	const gHits hits_h2y2 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h2y2);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h2y2 %d \n", nhits_h2y2);
	for(int l = 0; l<nhits_h2y2; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h2y2.chan(l), hits_h2y2.pos(l));
	}
#endif
#endif
	
	stid = 2;//3-1
	detid = geometry::hodoplanesx[stid][0];
	int nhits_h3x1;
	const gHits hits_h3x1 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h3x1);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h3x1 %d \n", nhits_h3x1);
	for(int l = 0; l<nhits_h3x1; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h3x1.chan(l), hits_h3x1.pos(l));
	}
#endif

	detid = geometry::hodoplanesx[stid][1];
	int nhits_h3x2;
	const gHits hits_h3x2 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h3x2);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h3x2 %d \n", nhits_h3x2);
	for(int l = 0; l<nhits_h3x2; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h3x2.chan(l), hits_h3x2.pos(l));
	}
#endif

#ifdef HODO_Y_MATCH
	detid = geometry::hodoplanesy[stid][0];
	int nhits_h3y1;
	const gHits hits_h3y1 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h3y1);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h3y1 %d \n", nhits_h3y1);
	for(int l = 0; l<nhits_h3y1; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h3y1.chan(l), hits_h3y1.pos(l));
	}
#endif

	detid = geometry::hodoplanesy[stid][1];
	int nhits_h3y2;
	const gHits hits_h3y2 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h3y2);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h3y2 %d \n", nhits_h3y2);
	for(int l = 0; l<nhits_h3y2; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h3y2.chan(l), hits_h3y2.pos(l));
	}
#endif
#endif
	
	stid = 3;//4-1
	detid = geometry::hodoplanesx[stid][0];
	int nhits_h4x1;
	const gHits hits_h4x1 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h4x1);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h4x1 %d \n", nhits_h4x1);
	for(int l = 0; l<nhits_h4x1; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h4x1.chan(l), hits_h4x1.pos(l));
	}
#endif

	detid = geometry::hodoplanesx[stid][1];
	int nhits_h4x2;
	const gHits hits_h4x2 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h4x2);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h4x2 %d \n", nhits_h4x2);
	for(int l = 0; l<nhits_h4x2; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h4x2.chan(l), hits_h4x2.pos(l));
	}
#endif
	
#ifdef HODO_Y_MATCH
	detid = geometry::hodoplanesy[stid][0];
	int nhits_h4y1;
	const gHits hits_h4y1 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h4y1);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h4y1 %d \n", nhits_h4y1);
	for(int l = 0; l<nhits_h4y1; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h4y1.chan(l), hits_h4y1.pos(l));
	}
#endif

	detid = geometry::hodoplanesy[stid][1];
	int nhits_h4y2;
	const gHits hits_h4y2 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h4y2);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h4y2 %d \n", nhits_h4y2);
	for(int l = 0; l<nhits_h4y2; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h4y2.chan(l), hits_h4y2.pos(l));
	}
#endif
#endif

#ifdef DEBUG
	if(debug::EvRef==blockIdx.x)printf("thread %d dets %d %d %d %d %d %d %d %d \n", threadIdx.x, detid_list[0], detid_list[1], detid_list[2], detid_list[3], detid_list[4], detid_list[5], detid_list[6], detid_list[7]);
#endif
	
	bool maskhodo[4];
	
	//variables for 2D track fit
	short detID[8];
	float Y[8];
	float errY[8];
	float Z[8];
	
	float A_[4];
	float Ainv_[4];
	float B_[2];
	float Par[2];
	float ParErr[2];
	float chi2;

	float y0, ty;
	float err_y0, err_ty;
	
	//Arrays for other basic hit info
	short elID[8];
	float pos[8];
	float drift[8];
#ifdef FULLCODE
	float tdc[8];
#endif
	short sign[12];

	short nhits_uv;
	
	short nhits_u2, nhits_u3;
	short nhits_v2, nhits_v3;
	
	const short nbins_total = nbins_st2*nbins_st3;
	short bin2, bin3;
	
	int nu2, nv2, nu3, nv3;
	
	int ncomb_uv3, ncomb_uv2;
	
	short i_u2, i_u3, i_v2, i_v3;
	
	int i_uv3, i_uv2, i_hit;

	float y, err_y;
	
	float chi2min = 10000.1f;

	int localbin;

	//arrays for chi2
	float dd[12];
	float res[12];
	float p1x[12];
	float p1y[12];
	float p1z[12];
	float dpx[12];
	float dpy[12];
	float dpz[12];
	
	//TODO get Hodo hits stations 2, 3, 4...
	
	//get the tracks...
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	__shared__ float besttrackYZdata[THREADS_PER_BLOCK][datasizes::NTracksParam];
	
	//gTracklet tkl;
	unsigned int Ntracks;
	const gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);

#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf("actual thread %d, offset %d, ntracks %d\n", threadIdx.x, tkl_coll_offset+array_thread_offset, Ntracks ); 
#endif
		
	int trackthread;	
	float x0, tx;
	float err_x0, err_tx;
	float xmin, xmax;
	int nhits_x;
	bool update_track;
	//
	for(int i = 0; i<Ntracks; i++){
		trackthread = (int)Tracks.threadID(i);
#ifdef DEBUG
		if(blockIdx.x!=debug::EvRef){
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 3);
			continue;
		}
		//if(blockIdx.x==debug::EvRef && Tracks.nHits(i)==3 && Tracks.hits_chan(i, 0)==76 && Tracks.hits_chan(i, 2)==105)
		if( (threadIdx.x%2==0 && (int)Tracks.threadID(i)%2==1) || (threadIdx.x%2==1 && (int)Tracks.threadID(i)%2==0) ){
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 3);
			continue;
		}
		//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));

		if(blockIdx.x==debug::EvRef)printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
#endif
		x0 = Tracks.x0(i);
		tx = Tracks.tx(i);
		err_x0 = Tracks.err_x0(i);
		err_tx = Tracks.err_tx(i);
		nhits_x =  Tracks.nHits(i);
		update_track = false;
		chi2min = 10000.1f;
		
		//get the corresponding local bin!
		//binId = threadIdx.x+THREADS_PER_BLOCK*localbin;
		localbin = (Tracks.stationID(i)-trackthread)/THREADS_PER_BLOCK;
		
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("actual thread %d x0 %1.4f tx %1.4f nhits %d track thread %d stid %1.0f local bin %d \n", threadIdx.x, x0, tx, nhits_x, trackthread, Tracks.stationID(i), localbin);
#endif		
		stid = 2;
		projid = 1;
		//nx1 = make_hitpairs_in_station(hits_st1x, nhits_st1x, hits_st1xp, nhits_st1xp, hitpairs_x1, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);
		xmin = min((tx-err_tx)*z_st2u, (tx-err_tx)*z_st2up)+x0-err_x0-det_spacing[0];
		xmax = max((tx+err_tx)*z_st2u, (tx+err_tx)*z_st2up)+x0+err_x0+det_spacing[0];
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("x0 %1.4f tx %1.4f, xmin %1.4f xmax %1.4f z_st2u %1.4f  z_st2up %1.4f \n", x0, tx, xmin, xmax, z_st2u, z_st2up);
#endif		
		nu2 = make_hitpairs_in_station(hits_st2u, nhits_st2u, hits_st2up, nhits_st2up, hitpairs_u2, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid, planes, xmin, xmax);

		projid = 2;
		xmin = min((tx-err_tx)*z_st2v, (tx-err_tx)*z_st2vp)+x0-err_x0-det_spacing[2];
		xmax = max((tx+err_tx)*z_st2v, (tx+err_tx)*z_st2vp)+x0+err_x0+det_spacing[2];
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("x0 %1.4f tx %1.4f, xmin %1.4f xmax %1.4f z_st2v %1.4f  z_st2vp %1.4f \n", x0, tx, xmin, xmax, z_st2v, z_st2vp);
#endif		
		nv2 = make_hitpairs_in_station(hits_st2v, nhits_st2v, hits_st2vp, nhits_st2vp, hitpairs_v2, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid, planes, xmin, xmax);

		stid = st3;
		projid = 1;
		xmin = min((tx-err_tx)*z_st3u, (tx-err_tx)*z_st3up)+x0-err_x0-det_spacing[4];
		xmax = max((tx+err_tx)*z_st3u, (tx+err_tx)*z_st3up)+x0+err_x0+det_spacing[4];
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("x0 %1.4f tx %1.4f, xmin %1.4f xmax %1.4f z_st3u %1.4f  z_st3up %1.4f \n", x0, tx, xmin, xmax, z_st3u, z_st3up);
#endif		
		nu3 = make_hitpairs_in_station(hits_st3u, nhits_st3u, hits_st3up, nhits_st3up, hitpairs_u3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid, planes, xmin, xmax);

		projid = 2;
		projid = 2;
		xmin = min((tx-err_tx)*z_st3v, (tx-err_tx)*z_st3vp)+x0-err_x0-det_spacing[6];
		xmax = max((tx+err_tx)*z_st3v, (tx+err_tx)*z_st3vp)+x0+err_x0+det_spacing[6];
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("x0 %1.4f tx %1.4f, xmin %1.4f xmax %1.4f z_st3v %1.4f  z_st3vp %1.4f \n", x0, tx, xmin, xmax, z_st3v, z_st3vp);
#endif		
		nv3 = make_hitpairs_in_station(hits_st3v, nhits_st3v, hits_st3vp, nhits_st3vp, hitpairs_v3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid, planes, xmin, xmax);
		
		tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 3);
		
		//retrieve the bins in st2, st3
		bin2 = localbin%nbins_st2;
		bin3 = (localbin-bin2)/nbins_st2;
		
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("nu2 %d nu3 %d nv2 %d nv3 %d\n", nu2, nu3, nv2, nv3);
#endif		
		if(nu2 == 0 || nu3==0) continue;
		if(nv2 == 0 || nv3==0) continue;
		
		ncomb_uv3 = nu3*nv3;
		ncomb_uv2 = nu2*nv2;
		
		for(int m = 0; m<nhits_x; m++){
			detid = (Tracks.hits_detid(i, m));
			chan = (Tracks.hits_chan(i, m));
			dd[m] = (Tracks.hits_drift(i, m));
			sign[m] = (Tracks.hits_sign(i, m));
			res[m] = planes->resolution[detid];
			p1x[m] = x_bep(detid, chan, planes);
			p1y[m] = y_bep(detid, chan, planes);
			p1z[m] = z_bep(detid, chan, planes);
			dpx[m] = planes->deltapx[detid];
			dpy[m] = planes->deltapy[detid];
			dpz[m] = planes->deltapz[detid];
		}

#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf(" evt %d x0 %1.4f +- %1.4f tx %1.4f +- %1.4f\n", blockIdx.x, x0, err_x0, tx, err_tx);
#endif
		
				
		ty = 0;
		
		for(i_uv3 = 0; i_uv3<ncomb_uv3; i_uv3++){
			i_u3 = i_uv3%nu3;
			i_v3 = (i_uv3-i_u3)/nu3;
			
			nhits_uv = 0;
			//bin3+nbins_st3*i_u3
			if(hitpairs_u3[i_u3].first>=0){
				detID[nhits_uv] = detid_list[4];
				i_hit = hitpairs_u3[i_u3].first;
				Z[nhits_uv] = z_st3u;
				elID[nhits_uv] = (short)hits_st3u.chan(i_hit);
#ifdef FULLCODE
				tdc[nhits_uv] = hits_st3u.tdc(i_hit);
#endif
				sign[nhits_x+nhits_uv] = 0;
				drift[nhits_uv] = hits_st3u.drift(i_hit);
				pos[nhits_uv] = hits_st3u.pos(i_hit);
				if(calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y)){
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
					res[nhits_x+nhits_uv] = planes->resolution[detID[nhits_uv]];
					dd[nhits_x+nhits_uv] = drift[nhits_uv];
					p1x[nhits_x+nhits_uv] = x_bep(detID[nhits_uv], elID[nhits_uv], planes);
					p1y[nhits_x+nhits_uv] = y_bep(detID[nhits_uv], elID[nhits_uv], planes);
					p1z[nhits_x+nhits_uv] = z_bep(detID[nhits_uv], elID[nhits_uv], planes);
					dpx[nhits_x+nhits_uv] = planes->deltapx[detID[nhits_uv]];
					dpy[nhits_x+nhits_uv] = planes->deltapy[detID[nhits_uv]];
					dpz[nhits_x+nhits_uv] = planes->deltapz[detID[nhits_uv]];
					
					Y[nhits_uv] = y;
					errY[nhits_uv] = err_y;
					nhits_uv++;
					ty+= y/z_st3u;
					
				}
			}
			if(hitpairs_u3[i_u3].second>=0){
				detID[nhits_uv] = detid_list[5];
				i_hit = hitpairs_u3[i_u3].second;
				Z[nhits_uv] = z_st3up;
				elID[nhits_uv] = (short)hits_st3up.chan(i_hit);
#ifdef FULLCODE
				tdc[nhits_uv] = hits_st3up.tdc(i_hit);
#endif
				sign[nhits_x+nhits_uv] = 0;
				drift[nhits_uv] = hits_st3up.drift(i_hit);
				pos[nhits_uv] = hits_st3up.pos(i_hit);
				if(calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y)){
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
					res[nhits_x+nhits_uv] = planes->resolution[detID[nhits_uv]];
					dd[nhits_x+nhits_uv] = drift[nhits_uv];
					p1x[nhits_x+nhits_uv] = x_bep(detID[nhits_uv], elID[nhits_uv], planes);
					p1y[nhits_x+nhits_uv] = y_bep(detID[nhits_uv], elID[nhits_uv], planes);
					p1z[nhits_x+nhits_uv] = z_bep(detID[nhits_uv], elID[nhits_uv], planes);
					dpx[nhits_x+nhits_uv] = planes->deltapx[detID[nhits_uv]];
					dpy[nhits_x+nhits_uv] = planes->deltapy[detID[nhits_uv]];
					dpz[nhits_x+nhits_uv] = planes->deltapz[detID[nhits_uv]];

					Y[nhits_uv] = y;
					errY[nhits_uv] = err_y;
					nhits_uv++;
					ty+= y/z_st3up;
				}
			}
			
			nhits_u3 = nhits_uv;
			if(nhits_u3==0) continue;
			
			if(hitpairs_v3[i_v3].first>=0){
				detID[nhits_uv] = detid_list[6];
				i_hit = hitpairs_v3[i_v3].first;
				Z[nhits_uv] = z_st3v;
				elID[nhits_uv] = (short)hits_st3v.chan(i_hit);
#ifdef FULLCODE
				tdc[nhits_uv] = hits_st3v.tdc(i_hit);
#endif
				sign[nhits_x+nhits_uv] = 0;
				drift[nhits_uv] = hits_st3v.drift(i_hit);
				pos[nhits_uv] = hits_st3v.pos(i_hit);
				if(calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y)){
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
					res[nhits_x+nhits_uv] = planes->resolution[detID[nhits_uv]];
					dd[nhits_x+nhits_uv] = drift[nhits_uv];
					p1x[nhits_x+nhits_uv] = x_bep(detID[nhits_uv], elID[nhits_uv], planes);
					p1y[nhits_x+nhits_uv] = y_bep(detID[nhits_uv], elID[nhits_uv], planes);
					p1z[nhits_x+nhits_uv] = z_bep(detID[nhits_uv], elID[nhits_uv], planes);
					dpx[nhits_x+nhits_uv] = planes->deltapx[detID[nhits_uv]];
					dpy[nhits_x+nhits_uv] = planes->deltapy[detID[nhits_uv]];
					dpz[nhits_x+nhits_uv] = planes->deltapz[detID[nhits_uv]];

					Y[nhits_uv] = y;
					errY[nhits_uv] = err_y;
					nhits_uv++;
					ty+= y/z_st3v;
				}
			}
			if(hitpairs_v3[i_v3].second>=0){
				detID[nhits_uv] = detid_list[7];
				i_hit = hitpairs_v3[i_v3].second;
				Z[nhits_uv] = z_st3vp;
				elID[nhits_uv] = (short)hits_st3vp.chan(i_hit);
#ifdef FULLCODE
				tdc[nhits_uv] = hits_st3vp.tdc(i_hit);
#endif
				sign[nhits_x+nhits_uv] = 0;
				drift[nhits_uv] = hits_st3vp.drift(i_hit);
				pos[nhits_uv] = hits_st3vp.pos(i_hit);
				if(calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y)){
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
					res[nhits_x+nhits_uv] = planes->resolution[detID[nhits_uv]];
					dd[nhits_x+nhits_uv] = drift[nhits_uv];
					p1x[nhits_x+nhits_uv] = x_bep(detID[nhits_uv], elID[nhits_uv], planes);
					p1y[nhits_x+nhits_uv] = y_bep(detID[nhits_uv], elID[nhits_uv], planes);
					p1z[nhits_x+nhits_uv] = z_bep(detID[nhits_uv], elID[nhits_uv], planes);
					dpx[nhits_x+nhits_uv] = planes->deltapx[detID[nhits_uv]];
					dpy[nhits_x+nhits_uv] = planes->deltapy[detID[nhits_uv]];
					dpz[nhits_x+nhits_uv] = planes->deltapz[detID[nhits_uv]];

					Y[nhits_uv] = y;
					errY[nhits_uv] = err_y;
					nhits_uv++;
					ty+= y/z_st3vp;
				}
			}
			
			nhits_v3 = nhits_uv-nhits_u3;
			if(nhits_v3==0) continue;

			ty = ty/nhits_uv;
			
			for(i_uv2 = 0; i_uv2<ncomb_uv2; i_uv2++){
				i_u2 = i_uv2%nu2;
				i_v2 = (i_uv2-i_u2)/nu2;

				nhits_uv = nhits_u3+nhits_v3;
				
				if(hitpairs_u2[i_u2].first>=0){
					detID[nhits_uv] = detid_list[0];
					i_hit = hitpairs_u2[i_u2].first;
					Z[nhits_uv] = z_st2u;
					elID[nhits_uv] = (short)hits_st2u.chan(i_hit);
#ifdef FULLCODE
					tdc[nhits_uv] = hits_st2u.tdc(i_hit);
#endif
					sign[nhits_x+nhits_uv] = 0;
					drift[nhits_uv] = hits_st2u.drift(i_hit);
					pos[nhits_uv] = hits_st2u.pos(i_hit);
					if(calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y)){
#ifdef DEBUG
						if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
						//if( fabs(y-ty*z_st2u)/err_y<=20.0f ){//since ty is *very* rough let's be generous...				
							res[nhits_x+nhits_uv] = planes->resolution[detID[nhits_uv]];
							dd[nhits_x+nhits_uv] = drift[nhits_uv];
							p1x[nhits_x+nhits_uv] = x_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1y[nhits_x+nhits_uv] = y_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1z[nhits_x+nhits_uv] = z_bep(detID[nhits_uv], elID[nhits_uv], planes);
							dpx[nhits_x+nhits_uv] = planes->deltapx[detID[nhits_uv]];
							dpy[nhits_x+nhits_uv] = planes->deltapy[detID[nhits_uv]];
							dpz[nhits_x+nhits_uv] = planes->deltapz[detID[nhits_uv]];

							Y[nhits_uv] = y;
							errY[nhits_uv] = err_y;
							nhits_uv++;
						//}
					}
				}
				if(hitpairs_u2[i_u2].second>=0){
					detID[nhits_uv] = detid_list[1];
					i_hit = hitpairs_u2[i_u2].second;
					Z[nhits_uv] = z_st2up;
					elID[nhits_uv] = (short)hits_st2up.chan(i_hit);
#ifdef FULLCODE
					tdc[nhits_uv] = hits_st2up.tdc(i_hit);
#endif
					sign[nhits_x+nhits_uv] = 0;
					drift[nhits_uv] = hits_st2up.drift(i_hit);
					pos[nhits_uv] = hits_st2up.pos(i_hit);
					if(calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y)){
#ifdef DEBUG
						if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
						//if( fabs(y-ty*z_st2up)/err_y<=20.0f ){//since ty is *very* rough let's be generous...				
							res[nhits_x+nhits_uv] = planes->resolution[detID[nhits_uv]];
							dd[nhits_x+nhits_uv] = drift[nhits_uv];
							p1x[nhits_x+nhits_uv] = x_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1y[nhits_x+nhits_uv] = y_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1z[nhits_x+nhits_uv] = z_bep(detID[nhits_uv], elID[nhits_uv], planes);
							dpx[nhits_x+nhits_uv] = planes->deltapx[detID[nhits_uv]];
							dpy[nhits_x+nhits_uv] = planes->deltapy[detID[nhits_uv]];
							dpz[nhits_x+nhits_uv] = planes->deltapz[detID[nhits_uv]];

							Y[nhits_uv] = y;
							errY[nhits_uv] = err_y;
							nhits_uv++;
						//}
					}
				}
				
				nhits_u2 = nhits_uv-nhits_u3-nhits_v3;
				if(nhits_u2==0) continue;
				
				if(hitpairs_v2[i_v2].first>=0){
					detID[nhits_uv] = detid_list[2];
					i_hit = hitpairs_v2[i_v2].first;
					Z[nhits_uv] = z_st2v;
					elID[nhits_uv] = (short)hits_st2v.chan(i_hit);
#ifdef FULLCODE
					tdc[nhits_uv] = hits_st2v.tdc(i_hit);
#endif
					sign[nhits_x+nhits_uv] = 0;
					drift[nhits_uv] = hits_st2v.drift(i_hit);
					pos[nhits_uv] = hits_st2v.pos(i_hit);
					if(calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y)){
#ifdef DEBUG
						if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
						//if( fabs(y-ty*z_st2v)/err_y<=20.0f ){//since ty is *very* rough let's be generous...				
							res[nhits_x+nhits_uv] = planes->resolution[detID[nhits_uv]];
							dd[nhits_x+nhits_uv] = drift[nhits_uv];
							p1x[nhits_x+nhits_uv] = x_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1y[nhits_x+nhits_uv] = y_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1z[nhits_x+nhits_uv] = z_bep(detID[nhits_uv], elID[nhits_uv], planes);
							dpx[nhits_x+nhits_uv] = planes->deltapx[detID[nhits_uv]];
							dpy[nhits_x+nhits_uv] = planes->deltapy[detID[nhits_uv]];
							dpz[nhits_x+nhits_uv] = planes->deltapz[detID[nhits_uv]];
							
							Y[nhits_uv] = y;
							errY[nhits_uv] = err_y;
							nhits_uv++;
						//}
					}
				}
				if(hitpairs_v2[i_v2].second>=0){
					detID[nhits_uv] = detid_list[3];
					i_hit = hitpairs_v2[i_v2].second;
					Z[nhits_uv] = z_st2vp;
					elID[nhits_uv] = (short)hits_st2vp.chan(i_hit);
#ifdef FULLCODE
					tdc[nhits_uv] = hits_st2vp.tdc(i_hit);
#endif
					sign[nhits_x+nhits_uv] = 0;
					drift[nhits_uv] = hits_st2vp.drift(i_hit);
					pos[nhits_uv] = hits_st2vp.pos(i_hit);
					if(calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y)){
#ifdef DEBUG
						if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
						//if( fabs(y-ty*z_st2vp)/err_y<=20.0f ){//since ty is *very* rough let's be generous...				
							res[nhits_x+nhits_uv] = planes->resolution[detID[nhits_uv]];
							dd[nhits_x+nhits_uv] = drift[nhits_uv];
							p1x[nhits_x+nhits_uv] = x_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1y[nhits_x+nhits_uv] = y_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1z[nhits_x+nhits_uv] = z_bep(detID[nhits_uv], elID[nhits_uv], planes);
							dpx[nhits_x+nhits_uv] = planes->deltapx[detID[nhits_uv]];
							dpy[nhits_x+nhits_uv] = planes->deltapy[detID[nhits_uv]];
							dpz[nhits_x+nhits_uv] = planes->deltapz[detID[nhits_uv]];
							
							Y[nhits_uv] = y;
							errY[nhits_uv] = err_y;
							nhits_uv++;
						//}
					}
				}
				
				nhits_v2 = nhits_uv-nhits_u3-nhits_v3-nhits_u2;
				if(nhits_v2==0) continue;
				
				fit_2D_track(nhits_uv, Y, Z, errY, A_, Ainv_, B_, Par, ParErr, chi2);
				
				y0 = Par[0];
				ty = Par[1];

				err_y0 = ParErr[0];
				err_ty = ParErr[1];
				
				if(fabs(y0)>Y0_MAX+2*err_y0 || fabs(ty)>TY_MAX+2*err_ty)continue;
				//if(fabs(y0+ty*geometry::Z_KMAG_BEND)>geometry::Y_KMAG_BEND)continue;

				//LR ambiguity resolution
				resolve_leftright_newhits(x0, tx, y0, ty, err_x0, err_tx, err_y0, err_ty, nhits_uv, detID, pos, drift, sign+nhits_x, planes, 150.);
				resolve_single_leftright_newhits(x0, tx, y0, ty, nhits_uv, detID, pos, sign+nhits_x, planes);
				
				for(int ll = 0; ll<nhits_uv; ll++){
					Y[ll]+= sign[ll]*drift[ll]*planes->sintheta[detID[ll]];
				}
				fit_2D_track(nhits_uv, Y, Z, errY, A_, Ainv_, B_, Par, ParErr, chi2);
				
				y0 = Par[0];
				ty = Par[1];

				err_y0 = ParErr[0];
				err_ty = ParErr[1];
				
				if(fabs(y0)>Y0_MAX || fabs(ty)>TY_MAX)continue;
				
#ifdef PROP_Y_MATCH
				//prop matching
				nprop = 0;
				checknext = true;
				
				for(int n = 0; n<nhits_p1y1; n++){
					ipos = hits_p1y1.pos(n);
					yExp = ty*z_p1y1+y0;
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("evt %d p1y1, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
					if(fabs(ipos-yExp)<5.08f){
						nprop++;
						checknext = false;
						break;
					}
				}
				if(checknext){
					for(int n = 0; n<nhits_p1y2; n++){
						ipos = hits_p1y2.pos(n);
						yExp = ty*z_p1y2+y0;
#ifdef DEBUG
						if(blockIdx.x<=debug::EvRef)printf("evt %d p1x2, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
						if(fabs(ipos-yExp)<5.08f){
							nprop++;
//							checknext = false;
							break;
						}
					}
				}
				checknext = true;
				
//				if(checknext){
				for(int n = 0; n<nhits_p2y1; n++){
					ipos = hits_p2y1.pos(n);
					yExp = ty*z_p2y1+y0;
#ifdef DEBUG
					if(blockIdx.x<=debug::EvRef)printf("evt %d p2x1, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
					if(fabs(ipos-yExp)<10.16f){
						nprop++;
						checknext = false;
						break;
					}
				}
//				}
				if(checknext){
					for(int n = 0; n<nhits_p2y2; n++){
						ipos = hits_p2y2.pos(n);
						yExp = ty*z_p2y2+y0;
#ifdef DEBUG
						if(blockIdx.x<=debug::EvRef)printf("evt %d p2y2, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
						if(fabs(ipos-yExp)<10.16f){
							nprop++;
//							checknext = false;
							break;
						}
					}
				}
				if(nprop<selection::NpropXYhitsMin)continue;
#endif	

#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf(" evt %d x0 %1.4f +- %1.4f tx %1.4f +- %1.4f y0 %1.4f +- %1.4f ty %1.4f +- %1.4f \n", blockIdx.x, x0, err_x0, tx, err_tx, y0, err_y0, ty, err_tx);
#endif
				//hodoscope matching
				stid = 1;//2-1

				maskhodo[stid] = 0;
				detid = geometry::hodoplanesx[stid][0];
				maskhodo[stid] = match_tracklet_to_hodo(stid, detid, nhits_h2x1, hits_h2x1, x0, y0, tx, ty, err_x0, err_y0, err_tx, err_ty, planes);
				if(!maskhodo[stid]){
					detid = geometry::hodoplanesx[stid][1];
					maskhodo[stid] = maskhodo[stid] || match_tracklet_to_hodo(stid, detid, nhits_h2x2, hits_h2x2, x0, y0, tx, ty, err_x0, err_y0, err_tx, err_ty, planes);
#ifdef HODO_Y_MATCH
					if(!maskhodo[stid]){
						detid = geometry::hodoplanesy[stid][0];
						maskhodo[stid] = maskhodo[stid] || match_tracklet_to_hodo(stid, detid, nhits_h2y1, hits_h2y1, x0, y0, tx, ty, err_x0, err_y0, err_tx, err_ty, planes);
						if(!maskhodo[stid]){
							detid = geometry::hodoplanesy[stid][1];
							maskhodo[stid] = maskhodo[stid] || match_tracklet_to_hodo(stid, detid, nhits_h2y2, hits_h2y2, x0, y0, tx, ty, err_x0, err_y0, err_tx, err_ty, planes);
						}
					}
#endif
				}

				if(!maskhodo[stid])continue;
				
				stid = 2;//3-1
				maskhodo[stid] = 0;
				detid = geometry::hodoplanesx[stid][0];
				maskhodo[stid] = match_tracklet_to_hodo(stid, detid, nhits_h3x1, hits_h3x1, x0, y0, tx, ty, err_x0, err_y0, err_tx, err_ty, planes);
				if(!maskhodo[stid]){
					detid = geometry::hodoplanesx[stid][1];
					maskhodo[stid] = maskhodo[stid] || match_tracklet_to_hodo(stid, detid, nhits_h3x2, hits_h3x2, x0, y0, tx, ty, err_x0, err_y0, err_tx, err_ty, planes);
#ifdef HODO_Y_MATCH
					if(!maskhodo[stid]){
						detid = geometry::hodoplanesy[stid][0];
						maskhodo[stid] = maskhodo[stid] || match_tracklet_to_hodo(stid, detid, nhits_h3y1, hits_h3y1, x0, y0, tx, ty, err_x0, err_y0, err_tx, err_ty, planes);
						if(!maskhodo[stid]){
							detid = geometry::hodoplanesy[stid][1];
							maskhodo[stid] = maskhodo[stid] || match_tracklet_to_hodo(stid, detid, nhits_h3y2, hits_h3y2, x0, y0, tx, ty, err_x0, err_y0, err_tx, err_ty, planes);
						}
					}
#endif
				}

				if(!maskhodo[stid])continue;
				
				stid = 3;//4-1
				maskhodo[stid] = 0;
				detid = geometry::hodoplanesx[stid][0];
				maskhodo[stid] = match_tracklet_to_hodo(stid, detid, nhits_h4x1, hits_h4x1, x0, y0, tx, ty, err_x0, err_y0, err_tx, err_ty, planes);
				if(!maskhodo[stid]){
					detid = geometry::hodoplanesx[stid][1];
					maskhodo[stid] = maskhodo[stid] || match_tracklet_to_hodo(stid, detid, nhits_h4x2, hits_h4x2, x0, y0, tx, ty, err_x0, err_y0, err_tx, err_ty, planes);
#ifdef HODO_Y_MATCH
					if(!maskhodo[stid]){
						detid = geometry::hodoplanesy[stid][0];
						maskhodo[stid] = maskhodo[stid] || match_tracklet_to_hodo(stid, detid, nhits_h4y1, hits_h4y1, x0, y0, tx, ty, err_x0, err_y0, err_tx, err_ty, planes);
						if(!maskhodo[stid]){
							detid = geometry::hodoplanesy[stid][1];
							maskhodo[stid] = maskhodo[stid] || match_tracklet_to_hodo(stid, detid, nhits_h4y2, hits_h4y2, x0, y0, tx, ty, err_x0, err_y0, err_tx, err_ty, planes);
						}
					}
#endif
				}
				if(!maskhodo[stid])continue;
				
				
				//chi2 evaluation of track candidate
				chi2 = chi2_track(nhits_x+nhits_uv, dd, sign, res, p1x, p1y, p1z, dpx, dpy, dpz, x0, y0, tx, ty);
				//if(chi2>1000000)continue;
#ifdef DEBUG
//				if(blockIdx.x==debug::EvRef && Tracks.hits_chan(i, 0)==29 && Tracks.hits_chan(i, 1)==30 && Tracks.hits_chan(i, 2)==11 && Tracks.hits_chan(i, 3)==11 && elID[0]==13 && elID[1]==13 && elID[2]==30 && elID[3]==30 && elID[4]==36 && elID[5]==36 && elID[6]==40 && elID[7]==41 ){
//				chi2 = chi2_track(nhits_x+nhits_uv, dd, sign, res, p1x, p1y, p1z, dpx, dpy, dpz, x0, y0, tx, ty, true);
//					printf("thread %d i %d x0 %1.4f y0 %1.4f tx %1.4f ty %1.4f chi2 %1.4f \n", threadIdx.x, i, x0, y0, tx, ty, chi2);
//				}
#endif
				if(chi2<chi2min){//establish a criterion to keep trakcs with the highest number of hitss
					chi2min = chi2;
					update_track = true;
					besttrackYZdata[threadIdx.x][2] = nhits_uv;
					besttrackYZdata[threadIdx.x][3] = chi2;
					besttrackYZdata[threadIdx.x][6] = ty;
					besttrackYZdata[threadIdx.x][8] = y0;
					besttrackYZdata[threadIdx.x][11] = err_ty;
					besttrackYZdata[threadIdx.x][13] = err_y0;
				
					for(int n = 0; n<nhits_uv;n++){
					besttrackYZdata[threadIdx.x][16+n] = detID[n];
					besttrackYZdata[threadIdx.x][34+n] = elID[n];
					besttrackYZdata[threadIdx.x][52+n] = pos[n];
					besttrackYZdata[threadIdx.x][70+n] = drift[n];
					besttrackYZdata[threadIdx.x][88+n] = sign[nhits_x+n];
#ifdef FULLCODE
					besttrackYZdata[threadIdx.x][106+n] = tdc[n];
					besttrackYZdata[threadIdx.x][124+n] = 0;
#endif
					}
				}
				
			}
		}
		
		if(update_track){
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 5);
			nhits_x = Tracks.nHits(i);
			nhits_uv = besttrackYZdata[threadIdx.x][2];
			tklcoll->setnHits(tkl_coll_offset+array_thread_offset, i, nhits_x+besttrackYZdata[threadIdx.x][2]);
			tklcoll->setChisq(tkl_coll_offset+array_thread_offset, i, besttrackYZdata[threadIdx.x][3]);
			tklcoll->setTy(tkl_coll_offset+array_thread_offset, i, besttrackYZdata[threadIdx.x][6]);
			tklcoll->setY0(tkl_coll_offset+array_thread_offset, i, besttrackYZdata[threadIdx.x][8]);
			tklcoll->setErrTy(tkl_coll_offset+array_thread_offset, i, besttrackYZdata[threadIdx.x][11]);
			tklcoll->setErrY0(tkl_coll_offset+array_thread_offset, i, besttrackYZdata[threadIdx.x][13]);
			
			for(int n = 0; n<nhits_uv; n++){
				tklcoll->setHitDetID(tkl_coll_offset+array_thread_offset, i, nhits_x+n, besttrackYZdata[threadIdx.x][16+n]);
				tklcoll->setHitChan(tkl_coll_offset+array_thread_offset, i, nhits_x+n, besttrackYZdata[threadIdx.x][34+n]);
				tklcoll->setHitPos(tkl_coll_offset+array_thread_offset, i, nhits_x+n, besttrackYZdata[threadIdx.x][52+n]);
				tklcoll->setHitDrift(tkl_coll_offset+array_thread_offset, i, nhits_x+n, besttrackYZdata[threadIdx.x][70+n]);
				tklcoll->setHitSign(tkl_coll_offset+array_thread_offset, i, nhits_x+n, besttrackYZdata[threadIdx.x][88+n]);
#ifdef FULLCODE
				tklcoll->setHitTDC(tkl_coll_offset+array_thread_offset, i, nhits_x+n, besttrackYZdata[threadIdx.x][106+n]);
				tklcoll->setHitResidual(tkl_coll_offset+array_thread_offset, i, nhits_x+n, besttrackYZdata[threadIdx.x][124+n]);
#endif
			}
		}
		
		
	}//end loop on existing tracklets
	
	
	
}


////////////////////////////////////
//
//          GLOBAL TRACKING
//
////////////////////////////////////

__global__ void gKernel_Global_tracking(
	gEventHitCollections* hitcolls,
	gEventTrackCollection* tklcoll,
	const gPlane* planes,
#ifdef DEBUG
	int* eventID,
#endif
	bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x]){
#ifdef DEBUG
		if(threadIdx.x==0)printf("Evt %d discarded, too many hits\n", eventID[blockIdx.x]);
#endif
		return;
	}
	
	short hitidx1[100];
	short hitidx2[100];
	short hitflag1[100];
	short hitflag2[100];
	
	int nx1, nu1, nv1;
	
        //pairs in station 1
        thrust::pair<int, int> hitpairs_x1[100];
        thrust::pair<int, int> hitpairs_u1[100];
        thrust::pair<int, int> hitpairs_v1[100];
	
	//get the hits
	unsigned int offset_hitcoll;
	
	short stid, projid, detid, detoff;
	
	projid = 0;
	stid = 0;
	detoff = 1;
	short detid_list[8];
	
	//Get the required hit collections here!
	projid = 0;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[0] = geometry::eff_detid_chambers[detid-1];
	int nhits_st1x;
	const gHits hits_st1x = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1x);
	const float z_st1x = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st1x = planes->resolution[detid];
#else
	const float res_st1x = planes->spacing[detid]*InvSqrt12;
#endif	
	
	detid-= 1;
	detid_list[1] = geometry::eff_detid_chambers[detid-1];
	int nhits_st1xp;
	const gHits hits_st1xp = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1xp);
	const float z_st1xp = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st1xp = planes->resolution[detid];
#else
	const float res_st1xp = planes->spacing[detid]*InvSqrt12;
#endif	

	projid = 1;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[2] = geometry::eff_detid_chambers[detid-1];
	int nhits_st1u;
	const gHits hits_st1u = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1u);
	const float z_st1u = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st1u = planes->resolution[detid];
#else
	const float res_st1u = planes->spacing[detid]*InvSqrt12;
#endif	
	
	detid-= 1;
	detid_list[3] = geometry::eff_detid_chambers[detid-1];
	int nhits_st1up;
	const gHits hits_st1up = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1up);
	const float z_st1up = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st1up = planes->resolution[detid];
#else
	const float res_st1up = planes->spacing[detid]*InvSqrt12;
#endif	
	
	projid = 2;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[4] = geometry::eff_detid_chambers[detid-1];
	int nhits_st1v;
	const gHits hits_st1v = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1v);
	const float z_st1v = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st1v = planes->resolution[detid];
#else
	const float res_st1v = planes->spacing[detid]*InvSqrt12;
#endif	
	
	detid-= 1;
	detid_list[5] = geometry::eff_detid_chambers[detid-1];
	int nhits_st1vp;
	const gHits hits_st1vp = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1vp);
	const float z_st1vp = planes->z[detid];
#ifdef USE_DET_RESOL
	const float res_st1vp = planes->resolution[detid];
#else
	const float res_st1vp = planes->spacing[detid]*InvSqrt12;
#endif	

	// TODO: load the hodoscope hits
	stid = 0;//1-1
	detid = geometry::hodoplanesx[stid][0];
	int nhits_h1x1;
	const gHits hits_h1x1 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h1x1);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h1x1 %d \n", nhits_h1x1);
	for(int l = 0; l<nhits_h1x1; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h1x1.chan(l), hits_h1x1.pos(l));
	}
#endif
	detid = geometry::hodoplanesx[stid][1];
	int nhits_h1x2;
	const gHits hits_h1x2 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h1x2);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h1x2 %d \n", nhits_h1x2);
	for(int l = 0; l<nhits_h1x2; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h1x2.chan(l), hits_h1x2.pos(l));
	}
#endif
	
#ifdef HODO_Y_MATCH
	detid = geometry::hodoplanesy[stid][0];
	int nhits_h1y1;
	const gHits hits_h1y1 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h1y1);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h1y1 %d \n", nhits_h1y1);
	for(int l = 0; l<nhits_h1y1; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h1y1.chan(l), hits_h1y1.pos(l));
	}
#endif

	detid = geometry::hodoplanesy[stid][1];
	int nhits_h1y2;
	const gHits hits_h1y2 = hitcolls->hitshodos(blockIdx.x, detid, nhits_h1y2);
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef && threadIdx.x==0){printf("nhits h1y2 %d \n", nhits_h1y2);
	for(int l = 0; l<nhits_h1y2; l++)printf("det %d chan %1.0f pos %1.4f \n", detid, hits_h1y2.chan(l), hits_h1y2.pos(l));
	}
#endif
#endif	
	bool maskhodo;
	
	//Sagitta ratio
	float pos_exp[3];
	float window[3];
	
	//full track calculation
	float x0_st1, tx_st1;
	float errx0_st1, errtx_st1;

	float x0, tx;
	float errx0, errtx;
	float y0, ty;
	float erry0, errty;
	
	float invP, errinvP;
	short charge;
	
	bool update_track;
	short nhits_st23, nhits_st1;
	
	//variables for the loop
	int i_x, i_u, i_v, i_uv;
	int ncomb_uv;
	int nhits_x, nhits_u, nhits_v, nhits_uv;
	int i_hit;
	
	//variables for 2D track fit
	short detID[6];
	short elID[6];
	float pos[6];
	float drift[6];
#ifdef FULLCODE
	float tdc[6];
#endif
	short sign[6];
	
	float Y[12];
	float errY[12];
	float Z_[12];
	short nyhits = 0;
	
	float X[3];
	float errX[3];
	float Z[3];
	
	float y, err_y;
			
	float A_[4];
	float Ainv_[4];
	float B_[2];
	float Par[2];
	float ParErr[2];
	float chi2;	
	
	float res[6];
	float p1x[6];
	float p1y[6];
	float p1z[6];
	float dpx[6];
	float dpy[6];
	float dpz[6];

	float chi2min = 10000.1f;

	//get the tracks...
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	__shared__ float besttrackdata[THREADS_PER_BLOCK][datasizes::NTracksParam];
	
	unsigned int Ntracks;
	const gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);
	
	for(int i = 0; i<Ntracks; i++){
		//tkl = Tracks.getTracklet(i);
		update_track = false;
		nhits_st23 = Tracks.nHits(i);
		
		if(Tracks.stationID(i)<5)continue;
		chi2min = 10000.1f;
		
		x0 = Tracks.x0(i);
		tx = Tracks.tx(i);
		errx0 = Tracks.err_x0(i);
		errtx = Tracks.err_tx(i);

		y0 = Tracks.y0(i);
		ty = Tracks.ty(i);
		erry0 = Tracks.err_y0(i);
		errty = Tracks.err_ty(i);
		
		SagittaRatioInStation1(x0, tx, y0, ty, Tracks.get_lasthitdetid(i), pos_exp, window, planes->z, planes->costheta, planes->sintheta);
		
		projid = 0;
		nx1 = make_hitpairs_in_station(hits_st1x, nhits_st1x, hits_st1xp, nhits_st1xp, hitpairs_x1, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);
		
		projid = 1;
		nu1 = make_hitpairs_in_station(hits_st1u, nhits_st1u, hits_st1up, nhits_st1up, hitpairs_u1, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);
		
		projid = 2;
		nv1 = make_hitpairs_in_station(hits_st1v, nhits_st1v, hits_st1vp, nhits_st1vp, hitpairs_v1, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);

#ifdef DEBUG
nx1, nu1, nv1);
		if(blockIdx.x==debug::EvRef){
			printf("x0 %1.4f y0 %1.4f tx %1.4f ty %1.4f nhits %1.0f detid %d, nx1 %d nu1 %d nv1 %d \n", x0, y0, tx, ty, Tracks.nHits(i), Tracks.get_lasthitdetid(i), nx1, nu1, nv1);
			for(int ll=0; ll<3; ll++)printf(" %d pos %1.4f window %1.4f  =>  %1.4f < %1.4f \n ", ll, pos_exp[ll], window[ll], pos_exp[ll]-window[ll], pos_exp[ll]+window[ll]);
		}
#endif
		
		if(nx1==0 || nu1==0 || nv1==0)continue;
		
		nyhits = 0;
		for(int k = 0; k<nhits_st23; k++){
			if(planes->costheta[(short)Tracks.hits_detid(i, k)]>0.99)continue;
			
			if(calculate_y_uvhit(Tracks.hits_detid(i, k), Tracks.hits_chan(i, k), Tracks.hits_drift(i, k), Tracks.hits_sign(i, k), x0, tx, planes, y, err_y)){
				Y[nyhits] = y;
				errY[nyhits] = err_y;
				Z_[nyhits] = planes->z[(short)Tracks.hits_detid(i, k)];
				nyhits++;
			}
		}
		
		ncomb_uv = nu1*nv1;
		
		//triple loop on hits
		for(i_x = 0; i_x<nx1; i_x++){
			nhits_x = 0;
			if(hitpairs_x1[i_x].first>=0){
				detID[nhits_x] = detid_list[0];
				i_hit = hitpairs_x1[i_x].first;
				X[nhits_x] = hits_st1x.pos(i_hit);
				pos[nhits_x] = hits_st1x.pos(i_hit);
				errX[nhits_x] = res_st1x;
				Z[nhits_x] = z_st1x;
				elID[nhits_x] = (short)hits_st1x.chan(i_hit);
#ifdef FULLCODE
				tdc[nhits_x] = hits_st1x.tdc(i_hit);
#endif
				sign[nhits_x] = 0;
				drift[nhits_x] = hits_st1x.drift(i_hit);
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_x], elID[nhits_x], pos[nhits_x], drift[nhits_x]);
#endif
				res[nhits_x] = planes->resolution[detID[nhits_x]];
				p1x[nhits_x] = x_bep(detID[nhits_x], elID[nhits_x], planes);
				p1y[nhits_x] = y_bep(detID[nhits_x], elID[nhits_x], planes);
				p1z[nhits_x] = z_bep(detID[nhits_x], elID[nhits_x], planes);
				dpx[nhits_x] = planes->deltapx[detID[nhits_x]];
				dpy[nhits_x] = planes->deltapy[detID[nhits_x]];
				dpz[nhits_x] = planes->deltapz[detID[nhits_x]];
				
				nhits_x++;
			}
			if(hitpairs_x1[i_x].second>=0){
				detID[nhits_x] = detid_list[1];
				i_hit = hitpairs_x1[i_x].second;
				X[nhits_x] = hits_st1xp.pos(i_hit);
				pos[nhits_x] = hits_st1xp.pos(i_hit);
				errX[nhits_x] = res_st1xp;
				Z[nhits_x] = z_st1xp;
				elID[nhits_x] = (short)hits_st1xp.chan(i_hit);
#ifdef FULLCODE
				tdc[nhits_x] = hits_st1xp.tdc(i_hit);
#endif
				sign[nhits_x] = 0;
				drift[nhits_x] = hits_st1xp.drift(i_hit);
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_x], elID[nhits_x], pos[nhits_x], drift[nhits_x]);
#endif
				res[nhits_x] = planes->resolution[detID[nhits_x]];
				p1x[nhits_x] = x_bep(detID[nhits_x], elID[nhits_x], planes);
				p1y[nhits_x] = y_bep(detID[nhits_x], elID[nhits_x], planes);
				p1z[nhits_x] = z_bep(detID[nhits_x], elID[nhits_x], planes);
				dpx[nhits_x] = planes->deltapx[detID[nhits_x]];
				dpy[nhits_x] = planes->deltapy[detID[nhits_x]];
				dpz[nhits_x] = planes->deltapz[detID[nhits_x]];
				
				nhits_x++;
			}
						
			if(nhits_x==0) continue;
			
			X[nhits_x] = x0+tx*geometry::Z_KMAG_BEND;
			errX[nhits_x] = errx0+errtx*geometry::Z_KMAG_BEND;
			Z[nhits_x] = geometry::Z_KMAG_BEND;

#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("x0 %1.4f y0 %1.4f tx %1.4f ty %1.4f nhits x %d first hit det %d chan %d pos %1.4f last hit det %d chan %d pos %1.4f, zmag x pos %1.4f \n", x0, y0, tx, ty, nhits_x, detID[0], elID[0], X[0], detID[nhits_x-1], elID[nhits_x-1], X[nhits_x-1], X[nhits_x]);
#endif
			
			fit_2D_track(nhits_x+1, X, Z, errX, A_, Ainv_, B_, Par, ParErr, chi2);
			x0_st1 = Par[0];
			tx_st1 = Par[1];

			errx0_st1 = Par[0];
			errtx_st1 = Par[1];

			resolve_leftright_xhits(x0_st1, tx_st1, errx0_st1, errtx_st1, nhits_x, detID, X, drift, sign, planes->z, 150.);
			resolve_single_leftright_xhits(x0_st1, tx_st1, nhits_x, detID, X, sign, planes->z);
			
			fit_2D_track_drift(nhits_x, X, drift, sign, Z, errX, A_, Ainv_, B_, Par, ParErr, chi2);
			x0_st1 = Par[0];
			tx_st1 = Par[1];

			errx0_st1 = Par[0];
			errtx_st1 = Par[1];
			
			for(short l = 0; l<nhits_x; l++){
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("l %d hit sign %d drift %1.4f \n", l, sign[l], drift[l]);
#endif
				X[l]+= sign[l]*drift[l];
			}
			fit_2D_track(nhits_x+1, X, Z, errX, A_, Ainv_, B_, Par, ParErr, chi2);
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("before: x0 %1.4f, tx %1.4f, after: x0 %1.4f tx %1.4f \n", x0_st1, tx_st1, Par[0], Par[1]);
#endif
			x0_st1 = Par[0];
			tx_st1 = Par[1];

			//charge = calculate_charge(tx, x0);
			invP = calculate_invP_charge(tx, tx_st1, charge);
			errinvP = calculate_invP_error(errtx, errtx_st1);
						
			if(invP>INVP_MAX+errinvP || invP<INVP_MIN-errinvP)continue;
			
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef){
				printf("thread %d tx_st1 %1.4f tx %1.4f invP %1.4f  %1.4f charge %d \n", threadIdx.x,  tx_st1, tx, (tx_st1 - tx) / geometry::PT_KICK_KMAG, invP, charge);
			}
			if(blockIdx.x==debug::EvRef)printf("x0 %1.4f y0 %1.4f tx %1.4f ty %1.4f invP %1.4f x0_st1 %1.4f tx_st1 %1.4f charge %d nhits x %d first hit det %d chan %d last hit det %d chan %d \n", x0, y0, tx, ty, invP, x0_st1, tx_st1, charge, nhits_x, detID[0], elID[0], detID[nhits_x-1], elID[nhits_x-1]);
#endif

							
			//add the UV hits
			for(i_uv = 0; i_uv<ncomb_uv; i_uv++){
				i_u = i_uv%nu1;
				i_v = (i_uv-i_u)/nu1;
				
				nhits_uv = 0;
				if(hitpairs_u1[i_u].first>=0){
					detID[nhits_x+nhits_uv] = detid_list[2];
					i_hit = hitpairs_u1[i_u].first;
					elID[nhits_x+nhits_uv] = (short)hits_st1u.chan(i_hit);
#ifdef FULLCODE
					tdc[nhits_x+nhits_uv] = hits_st1u.tdc(i_hit);
#endif
					sign[nhits_x+nhits_uv] = 0;
					drift[nhits_x+nhits_uv] = hits_st1u.drift(i_hit);
					pos[nhits_x+nhits_uv] = hits_st1u.pos(i_hit);
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %d pos %1.4f drift %1.4f\n", detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], pos[nhits_x+nhits_uv], drift[nhits_x+nhits_uv]);
#endif

					if(calculate_y_uvhit(detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], drift[nhits_x+nhits_uv], 0, x0_st1, tx_st1, planes, y, err_y)){
#ifdef DEBUG
						if(blockIdx.x==debug::EvRef){
							printf("det %d chan %d pos %1.4f drift %1.4f\n", detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], pos[nhits_x+nhits_uv], drift[nhits_x+nhits_uv]);
							printf(" ( %1.4f - %1.4f )^2 / (%1.4f^2 + %1.4f^2) = %1.4f \n", y, y_trk(y0, ty, z_st1u), err_y, err_y_trk(erry0, errty, z_st1u), (y-y_trk(y0, ty, z_st1u))*(y-y_trk(y0, ty, z_st1u)) / (err_y*err_y+err_y_trk(erry0, errty, z_st1u)*err_y_trk(erry0, errty, z_st1u)) );
						}
#endif
						//if( (y-y_trk(y0, ty, z_st1u))*(y-y_trk(y0, ty, z_st1u)) / (err_y*err_y+err_y_trk(erry0, errty, z_st1u)*err_y_trk(erry0, errty, z_st1u))<=12.0f ){
							res[nhits_x+nhits_uv] = planes->resolution[detID[nhits_uv]];
							p1x[nhits_x+nhits_uv] = x_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1y[nhits_x+nhits_uv] = y_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1z[nhits_x+nhits_uv] = z_bep(detID[nhits_uv], elID[nhits_uv], planes);
							dpx[nhits_x+nhits_uv] = planes->deltapx[detID[nhits_uv]];
							dpy[nhits_x+nhits_uv] = planes->deltapy[detID[nhits_uv]];
							dpz[nhits_x+nhits_uv] = planes->deltapz[detID[nhits_uv]];

							Y[nyhits+nhits_uv] = y;
							errY[nyhits+nhits_uv] = err_y;
							Z_[nyhits+nhits_uv] = z_st1u;
							
							nhits_uv++;
						//}
					}
				}
				if(hitpairs_u1[i_u].second>=0){
					detID[nhits_x+nhits_uv] = detid_list[3];
					i_hit = hitpairs_u1[i_u].second;
					elID[nhits_x+nhits_uv] = (short)hits_st1up.chan(i_hit);
#ifdef FULLCODE
					tdc[nhits_x+nhits_uv] = hits_st1up.tdc(i_hit);
#endif
					sign[nhits_x+nhits_uv] = 0;
					drift[nhits_x+nhits_uv] = hits_st1up.drift(i_hit);
					pos[nhits_x+nhits_uv] = hits_st1up.pos(i_hit);
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %d pos %1.4f drift %1.4f\n", detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], pos[nhits_x+nhits_uv], drift[nhits_x+nhits_uv]);
#endif
					if(calculate_y_uvhit(detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], drift[nhits_x+nhits_uv], 0, x0_st1, tx_st1, planes, y, err_y)){
#ifdef DEBUG
						if(blockIdx.x==debug::EvRef){
							printf("det %d chan %d pos %1.4f drift %1.4f\n", detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], pos[nhits_x+nhits_uv], drift[nhits_x+nhits_uv]);
							printf(" ( %1.4f - %1.4f )^2 / (%1.4f^2 + %1.4f^2) = %1.f \n", y, y_trk(y0, ty, z_st1up), err_y, err_y_trk(erry0, errty, z_st1up), (y-y_trk(y0, ty, z_st1up))*(y-y_trk(y0, ty, z_st1up))/(err_y*err_y+err_y_trk(erry0, errty, z_st1up)*err_y_trk(erry0, errty, z_st1up)) );
						}
#endif
						//if( (y-y_trk(y0, ty, z_st1up))*(y-y_trk(y0, ty, z_st1up))/(err_y*err_y+err_y_trk(erry0, errty, z_st1up)*err_y_trk(erry0, errty, z_st1up))<=12.0f ){
							res[nhits_x+nhits_uv] = planes->resolution[detID[nhits_uv]];
							p1x[nhits_x+nhits_uv] = x_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1y[nhits_x+nhits_uv] = y_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1z[nhits_x+nhits_uv] = z_bep(detID[nhits_uv], elID[nhits_uv], planes);
							dpx[nhits_x+nhits_uv] = planes->deltapx[detID[nhits_uv]];
							dpy[nhits_x+nhits_uv] = planes->deltapy[detID[nhits_uv]];
							dpz[nhits_x+nhits_uv] = planes->deltapz[detID[nhits_uv]];

							Y[nyhits+nhits_uv] = y;
							errY[nyhits+nhits_uv] = err_y;
							Z_[nyhits+nhits_uv] = z_st1up;
							
							nhits_uv++;
						//}
					}
				}
				nhits_u = nhits_uv;
				
				if(nhits_u==0)continue;
				
				if(hitpairs_v1[i_v].first>=0){
					detID[nhits_x+nhits_uv] = detid_list[4];
					i_hit = hitpairs_v1[i_v].first;
					elID[nhits_x+nhits_uv] = (short)hits_st1v.chan(i_hit);
#ifdef FULLCODE
					tdc[nhits_x+nhits_uv] = hits_st1v.tdc(i_hit);
#endif
					sign[nhits_x+nhits_uv] = 0;
					drift[nhits_x+nhits_uv] = hits_st1v.drift(i_hit);
					pos[nhits_x+nhits_uv] = hits_st1v.pos(i_hit);
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %d pos %1.4f drift %1.4f\n", detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], pos[nhits_x+nhits_uv], drift[nhits_x+nhits_uv]);
#endif
					if(calculate_y_uvhit(detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], drift[nhits_x+nhits_uv], 0, x0_st1, tx_st1, planes, y, err_y)){
#ifdef DEBUG
						if(blockIdx.x==debug::EvRef){
							printf("det %d chan %d pos %1.4f drift %1.4f\n", detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], pos[nhits_x+nhits_uv], drift[nhits_x+nhits_uv]);
							printf(" ( %1.4f - %1.4f )^2 / (%1.4f^2 + %1.4f^2) = %1.f \n", y, y_trk(y0, ty, z_st1v), err_y, err_y_trk(erry0, errty, z_st1v), (y-y_trk(y0, ty, z_st1v))*(y-y_trk(y0, ty, z_st1v))/(err_y*err_y+err_y_trk(erry0, errty, z_st1v)*err_y_trk(erry0, errty, z_st1v)) );
						}
#endif
						//if( (y-y_trk(y0, ty, z_st1v))*(y-y_trk(y0, ty, z_st1v))/(err_y*err_y+err_y_trk(erry0, errty, z_st1v)*err_y_trk(erry0, errty, z_st1v))<=12.0f ){
							res[nhits_x+nhits_uv] = planes->resolution[detID[nhits_uv]];
							p1x[nhits_x+nhits_uv] = x_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1y[nhits_x+nhits_uv] = y_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1z[nhits_x+nhits_uv] = z_bep(detID[nhits_uv], elID[nhits_uv], planes);
							dpx[nhits_x+nhits_uv] = planes->deltapx[detID[nhits_uv]];
							dpy[nhits_x+nhits_uv] = planes->deltapy[detID[nhits_uv]];
							dpz[nhits_x+nhits_uv] = planes->deltapz[detID[nhits_uv]];

							Y[nyhits+nhits_uv] = y;
							errY[nyhits+nhits_uv] = err_y;
							Z_[nyhits+nhits_uv] = z_st1v;
							
							nhits_uv++;
						//}
					}
				}
				if(hitpairs_v1[i_v].second>=0){
					detID[nhits_x+nhits_uv] = detid_list[5];
					i_hit = hitpairs_v1[i_v].second;
					elID[nhits_x+nhits_uv] = (short)hits_st1vp.chan(i_hit);
#ifdef FULLCODE
					tdc[nhits_x+nhits_uv] = hits_st1vp.tdc(i_hit);
#endif
					sign[nhits_x+nhits_uv] = 0;
					drift[nhits_x+nhits_uv] = hits_st1vp.drift(i_hit);
					pos[nhits_x+nhits_uv] = hits_st1vp.pos(i_hit);
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %d pos %1.4f drift %1.4f\n", detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], pos[nhits_x+nhits_uv], drift[nhits_x+nhits_uv]);
#endif
					if(calculate_y_uvhit(detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], drift[nhits_x+nhits_uv], 0, x0_st1, tx_st1, planes, y, err_y)){
#ifdef DEBUG
						if(blockIdx.x==debug::EvRef){
							printf(" ( %1.4f - %1.4f )^2 / (%1.4f^2 + %1.4f^2) = %1.f \n", y, y_trk(y0, ty, z_st1vp), err_y, err_y_trk(erry0, errty, z_st1vp), (y-y_trk(y0, ty, z_st1vp))*(y-y_trk(y0, ty, z_st1vp))/(err_y*err_y+err_y_trk(erry0, errty, z_st1vp)*err_y_trk(erry0, errty, z_st1vp)) );
						}
#endif
						//if( (y-y_trk(y0, ty, z_st1vp))*(y-y_trk(y0, ty, z_st1vp))/(err_y*err_y+err_y_trk(erry0, errty, z_st1vp)*err_y_trk(erry0, errty, z_st1vp))<=12.0f ){
							res[nhits_x+nhits_uv] = planes->resolution[detID[nhits_uv]];
							p1x[nhits_x+nhits_uv] = x_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1y[nhits_x+nhits_uv] = y_bep(detID[nhits_uv], elID[nhits_uv], planes);
							p1z[nhits_x+nhits_uv] = z_bep(detID[nhits_uv], elID[nhits_uv], planes);
							dpx[nhits_x+nhits_uv] = planes->deltapx[detID[nhits_uv]];
							dpy[nhits_x+nhits_uv] = planes->deltapy[detID[nhits_uv]];
							dpz[nhits_x+nhits_uv] = planes->deltapz[detID[nhits_uv]];

							Y[nyhits+nhits_uv] = y;
							errY[nyhits+nhits_uv] = err_y;
							Z_[nyhits+nhits_uv] = z_st1vp;
							
							nhits_uv++;
						//}
					}
				}
				nhits_v = nhits_uv-nhits_u;
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("x0 %1.4f y0 %1.4f tx %1.4f ty %1.4f invP %1.4f nhits uv %d first hit det %d chan %d last hit det %d chan %d \n", x0, y0, tx, ty, invP, nhits_uv, detID[nhits_x+nhits_u], elID[nhits_x+nhits_u], detID[nhits_x+nhits_uv-1], elID[nhits_x+nhits_uv-1]);
#endif
				if(nhits_v==0)continue;
				
				fit_2D_track(nyhits+nhits_uv, Y, Z_, errY, A_, Ainv_, B_, Par, ParErr, chi2);
				
				y0 = Par[0];
				ty = Par[1];

				erry0 = ParErr[0];
				errty = ParErr[1];
				
				if(fabs(y0)>Y0_MAX+2*erry0 || fabs(ty)>TY_MAX+2*errty)continue;

				//resolve left right...
				resolve_leftright_newhits(x0_st1, tx_st1, y0, ty, errx0_st1, errtx_st1, erry0, errty, nhits_x+nhits_uv, detID, pos, drift, sign, planes, 150.);
				resolve_single_leftright_newhits(x0_st1, tx_st1, y0, ty, nhits_x+nhits_uv, detID, pos, sign, planes);
				
				for(int ll = 0; ll<nhits_uv; ll++){
					Y[nyhits+ll]+= sign[ll]*drift[ll]*planes->sintheta[detID[nhits_x+nhits_uv]];
				}
				
				fit_2D_track(nyhits+nhits_uv, Y, Z_, errY, A_, Ainv_, B_, Par, ParErr, chi2);
				
				y0 = Par[0];
				ty = Par[1];

				erry0 = ParErr[0];
				errty = ParErr[1];
				
				if(fabs(y0)>Y0_MAX || fabs(ty)>TY_MAX)continue;
	
				// matching hodoscope here!
				stid = 0;//1-1
				maskhodo = 0;
				detid = geometry::hodoplanesx[stid][0];

#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf(" evt %d x0_st1 %1.4f +- %1.4f tx_st1 %1.4f +- %1.4f y0 %1.4f +- %1.4f ty %1.4f +- %1.4f invP %1.4f +- %1.4f \n", blockIdx.x, x0_st1, errx0_st1, tx_st1, errtx_st1, y0, erry0, ty, errty, invP, errinvP);
#endif

				maskhodo = match_tracklet_to_hodo(stid, detid, nhits_h1x1, hits_h1x1, x0_st1, y0, tx_st1, ty, errx0_st1, erry0, errtx_st1, errty, planes);
				if(!maskhodo){
					detid = geometry::hodoplanesx[stid][1];
					maskhodo = maskhodo || match_tracklet_to_hodo(stid, detid, nhits_h1x2, hits_h1x2, x0_st1, y0, tx_st1, ty, errx0_st1, erry0, errtx_st1, errty, planes);
				}
	
#ifdef HODO_Y_MATCH
				if(!maskhodo){
					detid = geometry::hodoplanesy[stid][0];
					maskhodo = match_tracklet_to_hodo(stid, detid, nhits_h1y1, hits_h1y1, x0_st1, y0, tx_st1, ty, errx0_st1, erry0, errtx_st1, errty, planes);
					if(!maskhodo){
						detid = geometry::hodoplanesy[stid][1];
						maskhodo = maskhodo || match_tracklet_to_hodo(stid, detid, nhits_h1y2, hits_h1y2, x0_st1, y0, tx_st1, ty, errx0_st1, erry0, errtx_st1, errty, planes);
					}
				}
				if(!maskhodo)continue;
#endif
				
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf(" x0_st1 %1.4f tx_st1 %1.4f y0 %1.4f ty %1.4f \n", x0_st1, tx_st1, y0, ty);
				
#endif
				
				//chi2 fit...
				chi2 = Tracks.chisq(i)+chi2_track(nhits_x+nhits_uv, drift, sign, res, p1x, p1y, p1z, dpx, dpy, dpz, x0_st1, y0, tx_st1, ty);
				if(chi2>10000.f)continue;
#ifdef DEBUG
//				if(blockIdx.x==debug::EvRef && Tracks.hits_chan(i, 0)==29 && Tracks.hits_chan(i, 1)==30 && Tracks.hits_chan(i, 2)==11 && Tracks.hits_chan(i, 3)==11 && Tracks.hits_chan(i, 4)==13 && Tracks.hits_chan(i, 5)==13 && Tracks.hits_chan(i, 6)==30 && Tracks.hits_chan(i, 7)==30 && Tracks.hits_chan(i, 8)==36 && Tracks.hits_chan(i, 9)==36 && Tracks.hits_chan(i, 10)==40 && Tracks.hits_chan(i, 11)==41 ){
//					chi2 = chi2_track(nhits_x+nhits_uv, drift, sign, res, p1x, p1y, p1z, dpx, dpy, dpz, x0_st1, y0, tx_st1, ty, true);
//					printf("x0 %1.4f y0 %1.4f tx %1.4f ty %1.4f invP %1.4f chi2 %1.4f \n", x0, y0, tx, ty, invP, chi2);
//				}
#endif				
				if(chi2<chi2min){
					chi2min = chi2;
					update_track = true;
					besttrackdata[threadIdx.x][2] = nhits_x+nhits_uv;
					besttrackdata[threadIdx.x][3] = chi2;
					besttrackdata[threadIdx.x][9] = invP;
					besttrackdata[threadIdx.x][14] = errinvP;
					besttrackdata[threadIdx.x][15] = charge;
					
					for(int m = 0; m<nhits_x+nhits_uv;m++){
						besttrackdata[threadIdx.x][16+m] = detID[m];
						besttrackdata[threadIdx.x][34+m] = elID[m];
						besttrackdata[threadIdx.x][52+m] = pos[m];
						besttrackdata[threadIdx.x][70+m] = drift[m];
						besttrackdata[threadIdx.x][88+m] = sign[m];
#ifdef FULLCODE
						besttrackdata[threadIdx.x][105+m] = tdc[m];
						besttrackdata[threadIdx.x][124+m] = 0;
#endif
					}
				}

				
			}
			
		}//end
		
		if(update_track){
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 6);
			nhits_st1 = besttrackdata[threadIdx.x][2];
			tklcoll->setnHits(tkl_coll_offset+array_thread_offset, i, nhits_st23+besttrackdata[threadIdx.x][2]);
			tklcoll->setChisq(tkl_coll_offset+array_thread_offset, i, besttrackdata[threadIdx.x][3]);
			tklcoll->setinvP(tkl_coll_offset+array_thread_offset, i, besttrackdata[threadIdx.x][9]);
			tklcoll->setErrinvP(tkl_coll_offset+array_thread_offset, i, besttrackdata[threadIdx.x][14]);
			tklcoll->setCharge(tkl_coll_offset+array_thread_offset, i, besttrackdata[threadIdx.x][15]);
				
			for(int n = 0; n<nhits_st1; n++){
				tklcoll->setHitDetID(tkl_coll_offset+array_thread_offset, i, nhits_st23+n, besttrackdata[threadIdx.x][16+n]);
				tklcoll->setHitChan(tkl_coll_offset+array_thread_offset, i, nhits_st23+n, besttrackdata[threadIdx.x][34+n]);
				tklcoll->setHitPos(tkl_coll_offset+array_thread_offset, i, nhits_st23+n, besttrackdata[threadIdx.x][52+n]);
				tklcoll->setHitDrift(tkl_coll_offset+array_thread_offset, i, nhits_st23+n, besttrackdata[threadIdx.x][70+n]);
				tklcoll->setHitSign(tkl_coll_offset+array_thread_offset, i, nhits_st23+n, besttrackdata[threadIdx.x][88+n]);
#ifdef FULLCODE
				tklcoll->setHitTDC(tkl_coll_offset+array_thread_offset, i, nhits_st23+n, besttrackdata[threadIdx.x][106+n]);
				tklcoll->setHitResidual(tkl_coll_offset+array_thread_offset, i, nhits_st23+n, besttrackdata[threadIdx.x][124+n]);
#endif
			}
		}
	}//end tracklets loop
}

// --------------- //
//
// Track vertexing 
//
// --------------- //

__global__ void gKernel_Vertexing(
	gEventTrackCollection* tklcoll,
	const float* z_array,
#ifdef DEBUG
	int* eventID,
#endif
	bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x]){
#ifdef DEBUG
		if(threadIdx.x==0)printf("Evt %d discarded, too many hits\n", eventID[blockIdx.x]);
#endif
		return;
	}
	
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	unsigned int Ntracks;
	const gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);
	
	short detid_first;
	
	float z_0;
	float state_0[5];
	//float cov_mat_0[25];
	
	float tx_i, tx_f;
	float ptot_i, ptot_b, ptot_f;
	float pz_0, pz_f;
	float traj1[3];
	float traj2[3];
	float pos_b[3];
	
	float x0, tx, x0_st1, tx_st1, y0, ty, invP, charge;
	
	int ix, iy, iz;
	int ix_m1, iy_m1, iz_m1;

	float pos_array[3*(geometry::NSTEPS_FMAG+geometry::NSTEPS_TARGET)];
	float mom_array[3*(geometry::NSTEPS_FMAG+geometry::NSTEPS_TARGET)];
	
	float dz;
	int step, step_x, step_y;
	
	float dca2_min, dca_xmin, dca_ymin;
	float dca2, dca_x, dca_y;

	//"output" observables	
	float dump_pos[3];
	float dump_mom[3];
	
	float vertex_pos[3];
	float vertex_mom[3];
	
	for(int i = 0; i<Ntracks; i++){
		if(Tracks.stationID(i)<6)continue;
		
		x0 = Tracks.x0(i);
		y0 = Tracks.y0(i);
		tx = Tracks.tx(i);
		ty = Tracks.ty(i);
		invP = Tracks.invP(i); 
		charge = Tracks.charge(i); 
		detid_first = Tracks.get_firsthitdetid(i);
		
		calculate_x0_tx_st1(x0, tx, invP, charge, x0_st1, tx_st1);
		
		z_0 = z_array[detid_first];
		
		state_0[0] = charge*invP*sqrtf(1.f+tx_st1*tx_st1+ty*ty);
		state_0[1] = tx_st1;
		state_0[2] = ty;
		state_0[3] = x_trk(x0_st1, tx_st1, z_0);
		state_0[4] = y_trk(y0, ty, z_0);

#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("front state %1.4f %1.4f %1.4f %1.4f %1.4f, %d z_0 =  %1.4f \n", state_0[0], state_0[1], state_0[2], state_0[3], state_0[4], detid_first, z_0);
#endif
		
		pos_array[0] = state_0[3]+state_0[1]*(geometry::FMAG_LENGTH-z_0);
		pos_array[1] = state_0[4]+state_0[2]*(geometry::FMAG_LENGTH-z_0);
		pos_array[2] = geometry::FMAG_LENGTH;
		
		pz_0 = 1.f/(state_0[0]*sqrtf(1.f+state_0[1]*state_0[1]+state_0[2]*state_0[2]));
		
		mom_array[0] = pz_0*state_0[1];
		mom_array[1] = pz_0*state_0[2];
		mom_array[2] = pz_0;
		
#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf(" FMAG_LENGTH-z_0 %1.4f pos(0) %1.4f %1.4f %1.4f  mom(0) %1.4f %1.4f %1.4f \n", geometry::FMAG_LENGTH-z_0, pos_array[0], pos_array[1], pos_array[2], mom_array[0], mom_array[1], mom_array[2] );
		if(blockIdx.x==debug::EvRef)printf("charge %1.0f\n", charge);
#endif

		step = 1;
		for(; step <= geometry::NSTEPS_FMAG; ++step){
			ix = step*3;
			iy = step*3+1;
			iz = step*3+2;
			ix_m1 = (step-1)*3;
			iy_m1 = (step-1)*3+1;
			iz_m1 = (step-1)*3+2;
			
			tx_i = mom_array[ix_m1]/mom_array[iz_m1];
			tx_f = tx_i + 2.f*charge*geometry::PTKICK_UNIT*geometry::STEP_FMAG/sqrt(mom_array[ix_m1]*mom_array[ix_m1]+mom_array[iz_m1]*mom_array[iz_m1]);
			
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("tx_i = %1.4f, tx_f = %1.4f, PTKICK_UNIT %1.4f, charge %1.0f, sqrt %1.4f \n", tx_i, tx_f, geometry::PTKICK_UNIT, charge, sqrt(mom_array[ix_m1]*mom_array[ix_m1]+mom_array[iz_m1]*mom_array[iz_m1]));
#endif
			
			traj1[0] = tx_i*geometry::STEP_FMAG;
			traj1[1] = ty*geometry::STEP_FMAG;
			traj1[2] = geometry::STEP_FMAG;
			
			pos_b[0] = pos_array[ix_m1]-traj1[0];
			pos_b[1] = pos_array[iy_m1]-traj1[1];
			pos_b[2] = pos_array[iz_m1]-traj1[2];
			
			ptot_i = sqrtf(mom_array[ix_m1]*mom_array[ix_m1]+mom_array[iy_m1]*mom_array[iy_m1]+mom_array[iz_m1]*mom_array[iz_m1]);
			ptot_b = ptot_i;
			if(pos_b[2] > geometry::FMAG_HOLE_LENGTH || pos_b[0]*pos_b[0]+pos_b[1]*pos_b[1]>geometry::FMAG_HOLE_RADIUS){
				ptot_b+= (geometry::DEDX_UNIT_0 + geometry::DEDX_UNIT_1*ptot_i + geometry::DEDX_UNIT_2*ptot_i*ptot_i + geometry::DEDX_UNIT_3*ptot_i*ptot_i*ptot_i + geometry::DEDX_UNIT_4*ptot_i*ptot_i*ptot_i*ptot_i)*sqrtf( traj1[0]*traj1[0] + traj1[1]*traj1[1] + traj1[2]*traj1[2]);
			}
			
			traj2[0] = tx_f*geometry::STEP_FMAG;
			traj2[1] = ty*geometry::STEP_FMAG;
			traj2[2] = geometry::STEP_FMAG;
			
			pos_array[ix] = pos_b[0]-traj2[0];
			pos_array[iy] = pos_b[1]-traj2[1];
			pos_array[iz] = pos_b[2]-traj2[2];
			
			ptot_f = ptot_b;
			if(pos_array[iz] > geometry::FMAG_HOLE_LENGTH || pos_b[ix]*pos_b[ix]+pos_b[iy]*pos_b[iy]>geometry::FMAG_HOLE_RADIUS){
				ptot_f+= (geometry::DEDX_UNIT_0 + geometry::DEDX_UNIT_1*ptot_b + geometry::DEDX_UNIT_2*ptot_b*ptot_b + geometry::DEDX_UNIT_3*ptot_b*ptot_b*ptot_b + geometry::DEDX_UNIT_4*ptot_b*ptot_b*ptot_b*ptot_b)*sqrtf( traj2[0]*traj2[0] + traj2[1]*traj2[1] + traj2[2]*traj2[2]);
			}
			
			pz_f = ptot_f/sqrtf(1.f+tx_f*tx_f+ty*ty);

			mom_array[ix] = pz_f*tx_f;
			mom_array[iy] = pz_f*ty;
			mom_array[iz] = pz_f;
			
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("%d %d %d %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f \n", ix, iy, iz, pos_array[ix], pos_array[iy], pos_array[iz], mom_array[ix], mom_array[iy], mom_array[iz] );
#endif

			dz = pos_b[2] - geometry::Z_DUMP;
			if(fabs(dz)<geometry::STEP_FMAG){
				if(dz<0){
					dump_pos[0] = pos_b[0]+dz*tx_f;
					dump_pos[1] = pos_b[1]+dz*ty;
					dump_pos[2] = pos_b[2]+dz;

					dump_mom[0] = mom_array[ix];
					dump_mom[1] = mom_array[iy];
					dump_mom[2] = mom_array[iz];
				}else{
					dump_pos[0] = pos_b[0]+dz*tx_i;
					dump_pos[1] = pos_b[1]+dz*ty;
					dump_pos[2] = pos_b[2]+dz;
					
					dump_mom[0] = mom_array[ix_m1];
					dump_mom[1] = mom_array[iy_m1];
					dump_mom[2] = mom_array[iz_m1];
				}
			}
		}//end loop on FMAG steps
		
		for(; step<=geometry::NSTEPS_FMAG+geometry::NSTEPS_TARGET; ++step){
			ix = step*3;
			iy = step*3+1;
			iz = step*3+2;
			ix_m1 = (step-1)*3;
			iy_m1 = (step-1)*3+1;
			iz_m1 = (step-1)*3+2;

			tx_i = mom_array[ix_m1]/mom_array[iz_m1];

			traj1[0] = tx_i*geometry::STEP_TARGET;
			traj1[1] = ty*geometry::STEP_TARGET;
			traj1[2] = geometry::STEP_TARGET;

			mom_array[ix] = mom_array[ix_m1];
			mom_array[iy] = mom_array[iy_m1];
			mom_array[iz] = mom_array[iz_m1];
			
			pos_array[ix] = pos_array[ix_m1]-traj1[0];
			pos_array[iy] = pos_array[iy_m1]-traj1[1];
			pos_array[iz] = pos_array[iz_m1]-traj1[2];
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("%d %d %d %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f \n", ix, iy, iz, pos_array[ix], pos_array[iy], pos_array[iz], mom_array[ix], mom_array[iy], mom_array[iz] );
#endif
		}//end loop on TARGET steps
		
		//now?
		dca2_min = 1.e18;
		dca_xmin = 1.e9;
		dca_ymin = 1.e9;
		
		step = geometry::NSTEPS_FMAG+geometry::NSTEPS_TARGET;
		step_x = step;
		step_y = step;
		
		for(int j = 0; j<=geometry::NSTEPS_FMAG+geometry::NSTEPS_TARGET; j++){
			ix = j*3;
			iy = j*3+1;
			iz = j*3+2;
			if(geometry::FMAGSTR*charge*mom_array[ix] < 0.) continue;
			
			dca2 = (pos_array[ix]-geometry::X_BEAM)*(pos_array[ix]-geometry::X_BEAM) + (pos_array[iy]-geometry::Y_BEAM)*(pos_array[iy]-geometry::Y_BEAM);

			if(dca2<dca2_min){
				dca2_min  =  dca2;
				step = j;
			}
			
			dca_x = fabs(pos_array[ix]-geometry::X_BEAM);
			if(dca_xmin>dca_x){
				dca_xmin  =  dca_x;
				step_x = j;
			}
			
			dca_y = fabs(pos_array[iy]-geometry::Y_BEAM);
			if(dca_ymin>dca_y){
				dca_ymin  =  dca_y;
				step_y = j;
			}
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("%d %d %d %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f %1.4f \n", ix, iy, iz, pos_array[ix], pos_array[iy], pos_array[iz], mom_array[ix], mom_array[iy], mom_array[iz], dca2, dca_x, dca_y );
#endif
		}
		
		ix = step*3;
		iy = step*3+1;
		iz = step*3+2;
		
		vertex_pos[0] = pos_array[ix];
		vertex_pos[1] = pos_array[iy];
		vertex_pos[2] = pos_array[iz];
		
		vertex_mom[0] = mom_array[ix];
		vertex_mom[1] = mom_array[iy];
		vertex_mom[2] = mom_array[iz];

		tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 7);//vertexing has been done...

		tklcoll->setVtxPos(tkl_coll_offset+array_thread_offset, i, vertex_pos);
		tklcoll->setVtxMom(tkl_coll_offset+array_thread_offset, i, vertex_mom);
		
		
	}//end loop on tracks
	
	
}


// --------------------------------------------- //
//
// !!! track cleaning kernels !!!
// 
// --------------------------------------------- //


__global__ void gKernel_TrackBadHitRemoval(
	gEventTrackCollection* tklcoll,
	const short track_stid_ref,
	const gPlane* planes,
	const bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x])return;
	
	//tracks
	unsigned int Ntracks;
	const gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);
	
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	int i;
	short ihit, nhitsi;
	//float res_curr, res_max, res_max2;
	//short hit_rm_idx; 
	short hit_rm_neighbor;
	float cut;
	
	float x0, errx0, y0, erry0, tx, errtx, ty, errty, x0_st1, tx_st1;
	float y, err_y;
	short charge;
	float invP, errinvP;
	
	const short ndof = track_stid_ref>5? 5 : 4;
	
	bool isupdated;
	//short signflip[18];
	//for(short k = 0; k<18; k++){
	//	signflip[k] = 0;
	//}
	float resid_[18];	
	
	short stid;
	short detid[18];
	short chan[18];
	short sign[18];
	float pos[18];
	float drift[18];
	
	float X_[4];
	float errX_[4];
	float Z_23[4];
	float X_st1_[3];
	float errX_st1_[3];
	float Z_1[3];
	float Y_[12];
	float errY_[12];
	float Z_[12];
	
	float A_[4];
	float Ainv_[4];
	float B_[2]; 
	float Par[2];
	float ParErr[2];
	float chi2;
	
#ifdef FULLCODE
	float tdc[18];
#endif
	short nhits_st23, nhits_st1, nhits_x, nhits_uv, nhits, nhits_x_st23, nhits_x_st1;
	
	for(i = 0; i<Ntracks; i++){
//#ifdef DEBUG
		if(blockIdx.x==debug::EvRef)printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
//#endif
		if(Tracks.stationID(i)<track_stid_ref)continue;
		nhitsi = Tracks.nHits(i);
		nhits_st23 = 0;
		nhits_st1 = 0;
		
		x0 = Tracks.x0(i);
		y0 = Tracks.y0(i);
		tx = Tracks.tx(i);
		ty = Tracks.ty(i);

		//res_max = -1.;
		//res_max2 = -1.;
		//hit_rm_idx = -1;
		hit_rm_neighbor = -1;
		
		if(track_stid_ref>5)calculate_x0_tx_st1(x0, tx, Tracks.invP(i), Tracks.charge(i), x0_st1, tx_st1);
		
		isupdated = false;
		//first fill in the arrays... that makes sense...
		for(ihit = 0; ihit<nhitsi; ihit++){
			detid[ihit] = Tracks.hits_detid(i, ihit);
			chan[ihit] = Tracks.hits_chan(i, ihit);
			sign[ihit] = Tracks.hits_sign(i, ihit);
			pos[ihit] = Tracks.hits_pos(i, ihit);
			drift[ihit] = Tracks.hits_drift(i, ihit);
#ifdef FULLCODE
			tdc[ihit] = Tracks.hits_tdc(i, ihit);
#endif

			if(detid[ihit]>12){
				nhits_st23++;
				resid_[ihit] = residual(detid[ihit], chan[ihit], drift[ihit], sign[ihit], planes, x0, y0, tx, ty);
			}else{
				resid_[ihit] = residual(detid[ihit], chan[ihit], drift[ihit], sign[ihit], planes, x0_st1, y0, tx_st1, ty);
				nhits_st1++;
			}
//#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("hit %d detid %d drift %1.4f resid %1.4f sign %d \n", ihit, detid[ihit], drift[ihit], resid_[ihit], sign[ihit]);
//#endif
		}
		//then reloop *once* to apply the logic....
		for(ihit = 0; ihit<nhitsi; ihit++){
			
			hit_rm_neighbor = detid[ihit]%2==0? ihit+1: ihit-1;
			if( abs( detid[hit_rm_neighbor] - detid[ihit] )>1 ||  detid[hit_rm_neighbor]<0 ||  hit_rm_neighbor>=nhitsi )hit_rm_neighbor = -1;
			
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf(" hit %d det %d neighbor %d det %d \n", ihit, detid[ihit], hit_rm_neighbor, detid[hit_rm_neighbor]);
#endif			
			stid = (detid[ihit]-1)/6-1;
			if(stid<0 && detid[ihit]>=0)stid = 0;
			
			cut = sign[ihit]==0 ? selection::rejectwin[stid]+fabs(drift[ihit]) : selection::rejectwin[stid];

			if( fabs(resid_[ihit]) > cut ){
				//if the neighbor hit id good, is larger than current hit and has a larger residual, skip because we will do it next 
				if(hit_rm_neighbor>ihit && resid_[hit_rm_neighbor]>resid_[ihit]){
//#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("Check next hit instead\n");
//#endif
					continue;
				}

//#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("cut %1.4f hit %d detid %d resid %1.4f drift %1.4f sign %d \n", cut, ihit, detid[ihit], resid_[ihit], drift[ihit], sign[ihit]);
//#endif
						
				isupdated = true;
				
				if(fabs(resid_[ihit]-2*drift[ihit]*sign[ihit]) < cut){
//#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("hit %d det %d sign changed \n", ihit, detid[ihit]);
//#endif
					sign[ihit] = -sign[ihit];
					if(hit_rm_neighbor>=0)sign[hit_rm_neighbor] = 0;
				}else{
					detid[ihit] = -1;
					if(hit_rm_neighbor<0){
						// if we have two hits in the same chamber that are not good, then the track cannot be kept 
						// we don't need to do the follow up either
//#ifdef DEBUG
						if(blockIdx.x==debug::EvRef)printf("evt %d thread %d track removed \n", blockIdx.x, threadIdx.x);
//#endif
						tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, Tracks.stationID(i)-2);
						isupdated = false;
						break;
					}
					if(detid[ihit]>12){
						nhits_st23--;
					}else{
						nhits_st1--;
					}
				}
				
				
			}
		}
		if(isupdated){
//#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("evt %d thread %d x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f stid %1.0f \n", blockIdx.x, threadIdx.x, x0, tx, y0, ty, Tracks.stationID(i));
//#endif			
			nhits_x = 0;
			nhits_uv = 0;
			nhits_x_st23 = 0;
			nhits_x_st1 = 0;
			
			resolve_single_leftright_newhits(x0, tx, y0, ty, nhits_st23, detid, pos, sign, planes);
			resolve_single_leftright_newhits(x0_st1, tx_st1, y0, ty, nhits_st1, detid+nhits_st23, pos+nhits_st23, sign+nhits_st23, planes);

/*
			for(ihit = 0; ihit<nhitsi; ihit++){
				if(detid[ihit]>=0 && (detid[ihit]-1)%6==2 || (detid[ihit]-1)%6==3 ){
					if( detid[ihit]>12 ){
						Z_23[nhits_x_st23] = planes->z[detid[ihit]];
						X_[nhits_x_st23] = pos[ihit];
						errX_[nhits_x_st23] = planes->spacing[detid[ihit]]*InvSqrt12;
						nhits_x_st23++;
					}else{
						Z_1[nhits_x_st1] = planes->z[detid[ihit]];
						X_st1_[nhits_x_st1] = pos[ihit];
						errX_st1_[nhits_x_st1] = planes->spacing[detid[ihit]]*InvSqrt12;
						nhits_x_st1++;
					}
				}
			}
			
			fit_2D_track(nhits_x_st23, X_, Z_23, errX_, A_, Ainv_, B_, Par, ParErr, chi2);
			x0 = Par[0];
			errx0 = ParErr[0];
			tx = Par[1];
			errtx = ParErr[1];
			
			Z_1[nhits_x_st1] = geometry::Z_KMAG_BEND;
			X_st1_[nhits_x_st1] = (x0, tx, Z_1[nhits_x_st1]);
			errX_st1_[nhits_x_st1] = (errx0, errtx, Z_1[nhits_x_st1]);
			
			fit_2D_track(nhits_x_st1+1, X_st1_, Z_1, errX_st1_, A_, Ainv_, B_, Par, ParErr, chi2);
			x0_st1 = Par[0];
			tx_st1 = Par[1];
			
			invP = calculate_invP_charge(tx, tx_st1, charge);
			errinvP = calculate_invP_error(errtx, ParErr[1]);
			
			for(ihit = 0; ihit<nhitsi; ihit++){
				if(detid[ihit]>=0 && (detid[ihit]-1)%6!=2 && (detid[ihit]-1)%6!=3 ){
					if( detid[ihit]>12 ){
						calculate_y_uvhit(detid[ihit], chan[ihit], drift[ihit], sign[ihit], x0, tx, planes, y, err_y);
					}else{
						calculate_y_uvhit(detid[ihit], chan[ihit], drift[ihit], sign[ihit], x0_st1, tx_st1, planes, y, err_y);
					}
					Y_[nhits_uv] = y;
					errY_[nhits_uv] = err_y;
					Z_[nhits_uv] = planes->z[detid[ihit]];
					nhits_uv++;
				}
			}
			
			fit_2D_track(nhits_uv, Y_, Z_1, errX_st1_, A_, Ainv_, B_, Par, ParErr, chi2);
			y0 = Par[0];
			erry0 = ParErr[0];
			ty = Par[1];
			errty = ParErr[1];
			//update the track parameters:
			
			tklcoll->setTx(tkl_coll_offset+array_thread_offset, i, tx);
			tklcoll->setX0(tkl_coll_offset+array_thread_offset, i, x0);
			tklcoll->setErrTx(tkl_coll_offset+array_thread_offset, i, errtx);
			tklcoll->setErrX0(tkl_coll_offset+array_thread_offset, i, errx0);
			tklcoll->setTy(tkl_coll_offset+array_thread_offset, i, ty);
			tklcoll->setY0(tkl_coll_offset+array_thread_offset, i, y0);
			tklcoll->setErrTy(tkl_coll_offset+array_thread_offset, i, errty);
			tklcoll->setErrY0(tkl_coll_offset+array_thread_offset, i, erry0);
			tklcoll->setinvP(tkl_coll_offset+array_thread_offset, i, invP);
			tklcoll->setErrinvP(tkl_coll_offset+array_thread_offset, i, errinvP);
*/
			
			chi2 = 0;
			nhits = 0;
			for(ihit = 0; ihit<nhitsi; ihit++){
				if(detid[ihit]>=0){
					resid_[nhits] = residual(detid[ihit], chan[ihit], drift[ihit], sign[ihit], planes, x0, y0, tx, ty)/planes->spacing[detid[ihit]]/InvSqrt12;
					chi2+= resid_[nhits]*resid_[nhits]/planes->spacing[detid[ihit]]/planes->spacing[detid[ihit]]*InvSqrt12*InvSqrt12;
					tklcoll->setHitDetID(tkl_coll_offset+array_thread_offset, i, nhits, detid[ihit]);
					tklcoll->setHitChan(tkl_coll_offset+array_thread_offset, i, nhits, chan[ihit]);
					tklcoll->setHitPos(tkl_coll_offset+array_thread_offset, i, nhits, pos[ihit]);
					tklcoll->setHitDrift(tkl_coll_offset+array_thread_offset, i, nhits, drift[ihit]);
					tklcoll->setHitSign(tkl_coll_offset+array_thread_offset, i, nhits, sign[ihit]);
#ifdef FULLCODE
					tklcoll->setHitTDC(tkl_coll_offset+array_thread_offset, i, nhits, tdc[ihit]);
					tklcoll->setHitResidual(tkl_coll_offset+array_thread_offset, i, nhits, resid_[ihit]);
#endif
					nhits++;
				}
			}
			tklcoll->setnHits(tkl_coll_offset+array_thread_offset, i, nhits);
			tklcoll->setChisq(tkl_coll_offset+array_thread_offset, i, chi2);
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("evt %d thread %d x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f stid %1.0f \n", blockIdx.x, threadIdx.x, x0, tx, y0, ty, Tracks.stationID(i));
#endif

		}

		/*
		isupdated = true;
		
		while(isupdated){
			isupdated = false; 
			//removal of "bad hits":
			// first: evaluation of residual for each hit, to find the largest residual
			for(ihit = 0; ihit<nhitsi; ihit++){
				res_curr = resid_[ihit];
				
				if(res_max<res_curr){
					res_max = resid_[ihit];
					res_max2 = fabs(resid_[ihit]-2*drift[ihit]*sign[ihit]);
					
					hit_rm_idx = ihit;
					hit_rm_neighbor = detid[ihit]%2==0? ihit+1: ihit-1;
					if( abs( detid[hit_rm_neighbor] - detid[ihit] )>1 ||  detid[hit_rm_neighbor]<0 )hit_rm_neighbor = -1;
				}
			}
		
			if(hit_rm_idx<0)break;
				
			if(sign==0 && Tracks.chisq(i)/(nhitsi-ndof) < selection::chi2dofmax )break;
						
			cut = abs(sign[hit_rm_idx])*drift[hit_rm_idx] + planes->spacing[detid[hit_rm_idx]]*InvSqrt12;
		
			if(res_max > cut){
				if(res_max2 < cut && signflip[hit_rm_idx]<2){
					//flip the sign
					sign[hit_rm_idx] = -sign[hit_rm_idx];
					signflip[hit_rm_idx]++;
					if(hit_rm_neighbor>=0)sign[hit_rm_neighbor] = 0;
				}else{
					detid[hit_rm_idx] = -1;
					if(hit_rm_neighbor<0){
						tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, Tracks.stationID(i)-2);
						break;
					}
				}
			}
			isupdated = true;
		}
		*/
		/*
		if(update_track){
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 6);
			nhits_st1 = besttrackdata[threadIdx.x][2];
			tklcoll->setnHits(tkl_coll_offset+array_thread_offset, i, nhits_st23+nhits_st1);
			tklcoll->setChisq(tkl_coll_offset+array_thread_offset, i, besttrackdata[threadIdx.x][3]);
			tklcoll->setinvP(tkl_coll_offset+array_thread_offset, i, besttrackdata[threadIdx.x][9]);
			tklcoll->setErrinvP(tkl_coll_offset+array_thread_offset, i, besttrackdata[threadIdx.x][14]);
			tklcoll->setCharge(tkl_coll_offset+array_thread_offset, i, besttrackdata[threadIdx.x][15]);
				
		}
		*/
	}
}


// --------------------------------------------- //
//
// function to clean the tracks after XZ-YZ processing
// 
// --------------------------------------------- //

__global__ void gKernel_BackTrackCleaning(
	gEventTrackCollection* tklcoll,
	const bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x])return;
	
	//tracks
	unsigned int Ntracks;
	const gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);
	
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	int i, j; 
	short ihit, jhit;
	short nhitsi, nhitsj;
	
	
	int nhits_common = 0;	
		
	for(i = 0; i<Ntracks; i++){
		//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
		if(Tracks.stationID(i)<5)continue;
		nhitsi = Tracks.nHits(i);
		nhits_common = 0;
		
		for(j = i+1; j<Ntracks; j++){
		//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
			if(Tracks.stationID(j)<5)continue;
			nhitsj = Tracks.nHits(j);
			
			for(ihit = 0; ihit<nhitsi; ihit++){
				for(jhit = 0; jhit<nhitsj; jhit++){
					if( Tracks.hits_detid(j, jhit)==Tracks.hits_detid(i, ihit) && Tracks.hits_chan(j, jhit)==Tracks.hits_chan(i, ihit) ){
						nhits_common++;
						break;
					}
				}
			}
			
			if(Tracks.chisq(i)/(nhitsi-4) < Tracks.chisq(j)/(nhitsj-4)){
				if(nhits_common>nhitsi*0.3333f ){
					tklcoll->setStationID(tkl_coll_offset+array_thread_offset, j, 3);
				}
			
			}else{
				if(nhits_common>nhitsj*0.3333f ){
					tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 3);
				}
			
			} 

		}	

	}
}




// --------------------------------------------- //
//
// function to clean the tracks after full processing
// 
// --------------------------------------------- //

__global__ void gKernel_GlobalTrackCleaning(
	gEventHitCollections* hitcolls,
	gEventTrackCollection* tklcoll,
	const float* z_array,
	const bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x])return;
	
	int i, j;

	float x0, tx;
	float y0, ty;
	
	float invP;
	
	//proportional tubes
	short stid, projid, detid;
	short nprop;
	float xExp, yExp, ipos;
	bool checknext;
	float prop_pos[2];
	float prop_z[2];
	float a, b;
	bool goodsegmentx;
	bool goodsegmenty;

	int n, m;
	
	projid = 0;
	stid = 6-1;
	detid = geometry::detsuperid[stid][projid]*2;
	int nhits_p1x1;
	const gHits hits_p1x1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1x1);
	const float z_p1x1 = z_array[detid];

	detid-= 1;
	int nhits_p1x2;
	const gHits hits_p1x2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1x2);
	const float z_p1x2 = z_array[detid];
	
	stid = 7-1;
	detid = geometry::detsuperid[stid][projid]*2;
	int nhits_p2x1;
	const gHits hits_p2x1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2x1);
	const float z_p2x1 = z_array[detid];

	detid-= 1;
	int nhits_p2x2;
	const gHits hits_p2x2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2x2);
	const float z_p2x2 = z_array[detid];
	
	// proportional tube hits
	projid = 1;
	stid = 6-1;
	detid = geometry::detsuperid[stid][projid]*2;
	int nhits_p1y1;
	const gHits hits_p1y1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1y1);
	const float z_p1y1 = z_array[detid];

	detid-= 1;
	int nhits_p1y2;
	const gHits hits_p1y2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1y2);
	const float z_p1y2 = z_array[detid];
	
	stid = 7-1;
	detid = geometry::detsuperid[stid][projid]*2;
	int nhits_p2y1;
	const gHits hits_p2y1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2y1);
	const float z_p2y1 = z_array[detid];

	detid-= 1;
	int nhits_p2y2;
	const gHits hits_p2y2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2y2);
	const float z_p2y2 = z_array[detid];
	
	//tracks
	unsigned int Ntracks;
	const gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);
	
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	float cut = 0.03;
	
	for(i = 0; i<Ntracks; i++){
		//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
		if(Tracks.stationID(i)<6)continue;

		if(Tracks.chisq(i)/(Tracks.nHits(i)-5)>selection::chi2dofmax){//  || Tracks.chisq(i)>100.f
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 5);
			continue;
		}
		
		for(j = i+1; j<Ntracks; j++){
		//printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
			if(Tracks.stationID(j)<6)continue;
			
				
			if(Tracks.chisq(i)/(Tracks.nHits(i)-5) < Tracks.chisq(j)/(Tracks.nHits(j)-5)){
				if( (Tracks.invP(i)-Tracks.invP(j))/Tracks.invP(i) < selection::merge_thres)  tklcoll->setStationID(tkl_coll_offset+array_thread_offset, j, 5);
			}else{
				if( (Tracks.invP(i)-Tracks.invP(i))/Tracks.invP(j) < selection::merge_thres)  tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 5);
			} 

		}

		x0 = Tracks.x0(i);
		y0 = Tracks.y0(i);
		tx = Tracks.tx(i);
		ty = Tracks.ty(i);
		invP = Tracks.invP(i);
						
		//cut = max(invP*0.11825f, 0.00643f-0.00009f/invP+0.00000046f/invP/invP);
		cut = 0.03f;
		
		goodsegmentx = false;
		// *very rough* tracklet segment
		// loop on first plane first
		for(n = 0; n<nhits_p1x1; n++){
			ipos = hits_p1x1.pos(n);
			xExp = tx*z_p1x1+x0;
			prop_pos[0] = ipos;
			prop_z[0] = z_p1x1;
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("evt %d p1x1, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
			if(fabs(ipos-xExp)<5.08f){
				nprop++;
				for(m = 0; m<nhits_p2x1; m++){
					ipos = hits_p2x1.pos(m);
					xExp = tx*z_p2x1+x0;
					prop_pos[1] = ipos;
					prop_z[1] = z_p2x1;
#ifdef DEBUG
					if(blockIdx.x<=debug::EvRef)printf("evt %d p2x1, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
					if(fabs(ipos-xExp)>7.62f)continue;

					a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
					//if(blockIdx.x<=debug::EvRef && st3==3 && hits_p1x1.chan(n)==21 && hits_p2x1.chan(m)==19)printf("a %1.4f tx %1.4f, chan p1x1 %1.0f pos %1.4f, chan p2x1 %1.0f  pos %1.4f \n", a, tx, hits_p1x1.chan(n), prop_pos[0], hits_p2x1.chan(m), prop_pos[1] );
					if(fabs(a)>TX_MAX)continue;
					if(fabs(a-tx)>cut)continue;
					b = prop_pos[0]-a*prop_z[0];
					//if(blockIdx.x<=debug::EvRef && st3==3 && hits_p1x1.chan(n)==21 && hits_p2x1.chan(m)==19)printf("a %1.4f tx %1.4f, b %1.4f, x0 %1.4f, chan p1x1 %1.0f, chan p2x1 %1.0f \n", a, tx, b, x0, hits_p1x1.chan(n), hits_p2x1.chan(m) );
					if(fabs(b)>X0_MAX)continue;
					if(fabs( (tx*2028.19f+x0)-(a*2028.19f+b) )>3.0f )continue;
					//if(blockIdx.x<=debug::EvRef && threadIdx.x%2==0 && hits_p1x1.chan(n)==21 && hits_p2x1.chan(m)==19)printf("thread %d a %1.4f tx %1.4f, b %1.4f, x0 %1.4f, pos_tkl %1.4f pos_seg %1.4f, detID[0] %d, detID[2] %d elID[0-3] %d %d %d %d \n", threadIdx.x, a, tx, b, x0, tx*2028.19f+x0, a*2028.19f+b, detID[0], detID[2], elID[0], elID[1], elID[2], elID[3] );
					
					nprop++;
					checknext = false;
					goodsegmentx = true;
					break;
				}
				if(!goodsegmentx){
					for(m = 0; m<nhits_p2x2; m++){
					ipos = hits_p2x2.pos(m);
					xExp = tx*z_p2x2+x0;
					prop_pos[1] = ipos;
					prop_z[1] = z_p2x2;
#ifdef DEBUG
					if(blockIdx.x<=debug::EvRef)printf("evt %d p2x2, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
					if(fabs(ipos-xExp)>7.62f)continue;

					a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
					if(fabs(a)>TX_MAX)continue;
					if(fabs(a-tx)>cut)continue;
					b = prop_pos[0]-a*prop_z[0];
					if(fabs(b)>X0_MAX)continue;
					if(fabs( (tx*2028.19f+x0)-(a*2028.19f+b) )>3.0f )continue;
						
					nprop++;
					checknext = false;
					goodsegmentx = true;
					break;
					}
				}
			}
		}
		
		if(!goodsegmentx){
			for(n = 0; n<nhits_p1x2; n++){
				ipos = hits_p1x2.pos(n);
				xExp = tx*z_p1x2+x0;
				prop_pos[0] = ipos;
				prop_z[0] = z_p1x2;
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("evt %d p1x1, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
				if(fabs(ipos-xExp)<5.08f){
					nprop++;
					for(m = 0; m<nhits_p2x1; m++){
						ipos = hits_p2x1.pos(m);
						xExp = tx*z_p2x1+x0;
						prop_pos[1] = ipos;
						prop_z[1] = z_p2x1;
#ifdef DEBUG
						if(blockIdx.x<=debug::EvRef)printf("evt %d p2x1, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
						if(fabs(ipos-xExp)>7.62f)continue;
	
						a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
						if(fabs(a)>TX_MAX)continue;
						if(fabs(a-tx)>cut)continue;
						b = prop_pos[0]-a*prop_z[0];
						if(fabs(b)>X0_MAX)continue;
						if(fabs( (tx*2028.19f+x0)-(a*2028.19f+b) )>3.0f )continue;
							
						nprop++;
						checknext = false;
						goodsegmentx = true;
						break;
					}
					if(!goodsegmentx){
						for(m = 0; m<nhits_p2x2; m++){
						ipos = hits_p2x2.pos(m);
						xExp = tx*z_p2x2+x0;
						prop_pos[1] = ipos;
						prop_z[1] = z_p2x2;
#ifdef DEBUG
						if(blockIdx.x<=debug::EvRef)printf("evt %d p2x2, ipos = %1.4f, xExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, xExp, ipos-xExp);
#endif
						if(fabs(ipos-xExp)>7.62f)continue;
	
						a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
						if(fabs(a)>TX_MAX)continue;
						if(fabs(a-tx)>cut)continue;
						b = prop_pos[0]-a*prop_z[0];
						if(fabs(b)>X0_MAX)continue;
						if(fabs( (tx*2028.19f+x0)-(a*2028.19f+b) )>3.0f )continue;
							
						nprop++;
						checknext = false;
						goodsegmentx = true;
						break;
						}
					}
				}
			}
		}
		if(!goodsegmentx){
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 5);
			continue;
		}
		
		goodsegmenty = false;
		// *very rough* tracklet segment
		// loop on first plane first
		for(n = 0; n<nhits_p1y1; n++){
			ipos = hits_p1y1.pos(n);
			yExp = ty*z_p1y1+y0;
			prop_pos[0] = ipos;
			prop_z[0] = z_p1y1;
#ifdef DEBUG
			if(blockIdx.x==debug::EvRef)printf("evt %d p1y1, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
			if(fabs(ipos-yExp)<5.08f){
				nprop++;
				for(m = 0; m<nhits_p2y1; m++){
					ipos = hits_p2y1.pos(m);
					yExp = ty*z_p2y1+y0;
					prop_pos[1] = ipos;
					prop_z[1] = z_p2y1;
#ifdef DEBUG
					if(blockIdx.x<=debug::EvRef)printf("evt %d p2y1, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
					if(fabs(ipos-yExp)>10.16f)continue;

					a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
					if(fabs(a)>TY_MAX)continue;
					if(fabs(a-ty)>cut)continue;
					b = prop_pos[0]-a*prop_z[0];
					if(fabs(b)>Y0_MAX)continue;
					if(fabs( (ty*2028.19f+y0)-(a*2028.19f+b) )>3.0f )continue;
				
					nprop++;
					checknext = false;
					goodsegmenty = true;
					break;
				}
				if(!goodsegmenty){
					for(m = 0; m<nhits_p2y2; m++){
						ipos = hits_p2y2.pos(m);
						yExp = ty*z_p2y2+y0;
						prop_pos[1] = ipos;
						prop_z[1] = z_p2y2;
#ifdef DEBUG
						if(blockIdx.x<=debug::EvRef)printf("evt %d p2y2, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
						if(fabs(ipos-yExp)>10.16f)continue;
		
						a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
						if(fabs(a)>TY_MAX)continue;
						if(fabs(a-ty)>cut)continue;
						b = prop_pos[0]-a*prop_z[0];
						if(fabs(b)>Y0_MAX)continue;
						if(fabs( (ty*2028.19f+y0)-(a*2028.19f+b) )>3.0f )continue;
							
						nprop++;
						checknext = false;
						goodsegmenty = true;
						break;
					}
				}
			}
		}
			
		if(!goodsegmenty){
			for(n = 0; n<nhits_p1y2; n++){
				ipos = hits_p1y2.pos(n);
				yExp = ty*z_p1y2+y0;
				prop_pos[0] = ipos;
				prop_z[0] = z_p1y2;
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("evt %d p1y1, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
				if(fabs(ipos-yExp)<5.08f){
					nprop++;
					for(m = 0; m<nhits_p2y1; m++){
						ipos = hits_p2y1.pos(m);
						yExp = ty*z_p2y1+y0;
						prop_pos[1] = ipos;
						prop_z[1] = z_p2y1;
#ifdef DEBUG
						if(blockIdx.x<=debug::EvRef)printf("evt %d p2y1, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
						if(fabs(ipos-yExp)>10.16f)continue;
		
						a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
						if(fabs(a)>TY_MAX)continue;
						if(fabs(a-ty)>cut)continue;
						b = prop_pos[0]-a*prop_z[0];
						if(fabs(b)>Y0_MAX)continue;
						if(fabs( (ty*2028.19f+y0)-(a*2028.19f+b) )>3.0f )continue;
							
						nprop++;
						checknext = false;
						goodsegmenty = true;
						break;
					}
					if(!goodsegmenty){
						for(m = 0; m<nhits_p2y2; m++){
							ipos = hits_p2y2.pos(m);
							yExp = ty*z_p2y2+y0;
							prop_pos[1] = ipos;
							prop_z[1] = z_p2y2;
#ifdef DEBUG
							if(blockIdx.x<=debug::EvRef)printf("evt %d p2y2, ipos = %1.4f, yExp = %1.4f, diff = %1.4f \n", blockIdx.x, ipos, yExp, ipos-yExp);
#endif
							if(fabs(ipos-yExp)>10.16f)continue;
	
							a = (prop_pos[1]-prop_pos[0])/(prop_z[1]-prop_z[0]);
							if(fabs(a)>TY_MAX)continue;
							if(fabs(a-ty)>cut)continue;
							b = prop_pos[0]-a*prop_z[0];
							if(fabs(b)>Y0_MAX)continue;
							if(fabs( (ty*2028.19f+y0)-(a*2028.19f+b) )>3.0f )continue;
								
							nprop++;
							checknext = false;
							goodsegmenty = true;
							break;
						}
					}
				}
			}
		}
		if(!goodsegmenty){
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 5);
			continue;
		}
	
	}
	
	
	
}



// --------------------------------------------- //
//
// simple track printing function for debugging. 
// 
// --------------------------------------------- //
__global__ void gKernel_check_tracks(gEventTrackCollection* tklcoll, const bool* hastoomanyhits, const int blockID)
{
	if(hastoomanyhits[blockIdx.x])return;
	
	unsigned int nTracks;
	
	if(blockIdx.x==blockID){
		const gTracks Tracks = tklcoll->tracks(blockID, threadIdx.x, nTracks);
		
		for(int i = 0; i<nTracks; i++){
			printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
		}
	}
	
}



__global__ void gkernel_checkHistos(gHistsArrays *HistsArrays, const short hist_num){
	const int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if(blockIdx.x==hist_num || hist_num==-1){
		if(threadIdx.x==0)printf(" hw %1.4f ", HistsArrays->pts_hw[blockIdx.x]);
		printf(" thread %d idx %d xpts %1.4f values %1.4f\n", threadIdx.x, idx, HistsArrays->xpts[idx], HistsArrays->values[idx]);
	}
}


// --------------------------------------------- //
//
// function to fill the "histograms"/arrays to be displayed
// 
// --------------------------------------------- //
__global__ void gKernel_fill_display_histograms(gEventTrackCollection* tklcoll, gHistsArrays *HistsArrays, const bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x])return;
	
	unsigned int nTracks;
	float vars[NVars];
	float x_hw, x;
	
	int bin;
	int i, j, k;
	
	//__shared__ float values_[NVars*Nbins_Hists];
	//for(j = 0; j<Nbins_Hists*NVars; j++)values_[j] = 0;
	
	const gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, nTracks);
	
	for(i = 0; i<nTracks; i++){
		if(Tracks.stationID(i)<7 || Tracks.chisq(i)/(Tracks.nHits(i)-5)>selection::chi2dofmax)continue;// || Tracks.chisq(i)>100.f
		
		vars[0] = Tracks.x0(i);
		vars[1] = Tracks.y0(i);
		vars[2] = Tracks.invP(i);
		vars[3] = Tracks.tx(i);
		vars[4] = Tracks.ty(i);
		
		vars[5] = Tracks.vx(i);
		vars[6] = Tracks.vy(i);
		vars[7] = Tracks.vz(i);
		vars[8] = Tracks.px(i);
		vars[9] = Tracks.py(i);
		vars[10] = Tracks.pz(i);
		
		for(j = 0; j<Nbins_Hists; j++){
			//if(blockIdx.x==debug::EvRef && threadIdx)printf("%d %d %1.4f %1.4f \n", j, threadIdx.x, HistsArrays->xpts[j], HistsArrays->pts_hw[k]);
			for(k = 0; k<NVars; k++){
				bin = j+Nbins_Hists*k;
				x_hw = HistsArrays->pts_hw[k];
				x = HistsArrays->xpts[bin];
				if(x-x_hw<=vars[k] && vars[k]<x+x_hw){
					HistsArrays->values[bin]+=1.f;
					//values_[bin]+=1.f;
#ifdef DEBUG
					if(k==2)printf(" bin %d  %1.4f < %1.4f < %1.4f : %1.4f; ", bin, x-x_hw, vars[k], x+x_hw, HistsArrays->values[bin]);
#endif
				}
			}
		}
	}
	__syncthreads();
	//cudaDeviceSynchronize();
	//for(j = 0; j<Nbins_Hists*NVars; j++){
	//	HistsArrays->values[j]+= values_[j];
	//}
}


#ifdef OLDCODE


// ------------------------------
// Track Kalman fitting 
// ------------------------------

__global__ void gKernel_GlobalTrack_KalmanFitting(gOutputEvent* oc, gKalmanFitArrays* fitarrays, const gPlane* planes)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	if(oc[index].HasTooManyHits)return;
	
	int nhits_tkl;
	
	for(int k = 0; k<oc[index].nTracklets; k++){
		if(oc[index].AllTracklets[k].stationID<6)continue;
		
		nhits_tkl = oc[index].AllTracklets[k].nXHits+oc[index].AllTracklets[k].nUHits+oc[index].AllTracklets[k].nVHits;
		
		//initialize state
		//initialize_arrays(oc[index].AllTracklets[k], fitarrays[index]);
		
		
			if(oc[index].EventID<20)
		for(int i = 0; i<nhits_tkl; i++){
printf("hit %d, z = %1.3f \n", i, planes[oc[index].AllTracklets[k].hits[i].detectorID].z);
			// first predict the state at the hit z
			predict_state(oc[index].AllTracklets[k], planes[oc[index].AllTracklets[k].hits[i].detectorID].z, fitarrays[index]);
			update_state(oc[index].AllTracklets[k], oc[index].AllTracklets[k].hits[i], fitarrays[index], planes[oc[index].AllTracklets[k].hits[i].detectorID]);
		}
	}
	
	//
}


#endif





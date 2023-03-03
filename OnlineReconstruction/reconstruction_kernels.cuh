#include "reconstruction_helper.cuh"
//#include "tracknumericalminimizer.cuh"
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
	if(station_mult[0]>250)hastoomanyhits[blockIdx.x] = true;
	if(station_mult[1]>200)hastoomanyhits[blockIdx.x] = true;
	if(station_mult[2]>150)hastoomanyhits[blockIdx.x] = true;
	if(station_mult[3]>150)hastoomanyhits[blockIdx.x] = true;
	if(station_mult[4]>250)hastoomanyhits[blockIdx.x] = true;

#ifdef DEBUG
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
	const float* res_array,
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
	const float res_st2x = res_array[detid];
	
	detid-= 1;
	detid_list[1] = detid;
	int nhits_st2xp;
	const gHits hits_st2xp = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2xp);
	const float z_st2xp = z_array[detid];
	const float res_st2xp = res_array[detid];
	
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
	const float res_st3x = res_array[detid];

	detid-= 1;
	int nhits_st3xp;
	detid_list[3] = detid;
	const gHits hits_st3xp = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3xp);
#ifdef DEBUG	
	if(blockIdx.x==debug::EvRef && st3==4)printf("block %d detid %d nhits %d, vs %d \n", blockIdx.x, detid, nhits_st3xp, hitcolls->NHitsChambers[blockIdx.x*nChamberPlanes+detid-1]);
#endif
	const float z_st3xp = z_array[detid];
	const float res_st3xp = res_array[detid];

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

	float x0, tx;
	
	//Arrays for other basic hit info
	short elID[4];
	float drift[4];
	float tdc[4];
	
	//small varaibles for prop tube matching
	short nprop, iprop, idet;
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
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackletSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackletSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	const unsigned int tklmult_idx = blockIdx.x*THREADS_PER_BLOCK+threadIdx.x;

	tklcoll->NTracks[tklmult_idx] = 0;

	__shared__ unsigned int ntkl_per_thread[THREADS_PER_BLOCK];
	for(int k = 0; k<THREADS_PER_BLOCK; k++)ntkl_per_thread[k] = 0;
	__shared__ short thread_min[THREADS_PER_BLOCK];//thread with lowest number of tracks: each thread
	__shared__ bool addtrack[THREADS_PER_BLOCK];//flag to mark threads requesting to add a new track in another thread.
	__shared__ bool threadbusy[THREADS_PER_BLOCK];
	unsigned int ntkl_min;
	unsigned int array_offset;
	unsigned int nthreads_full;
	unsigned int nslots_available;
	//__shared__
	unsigned int list_of_threads[THREADS_PER_BLOCK];
	
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
			
			addtrack[threadIdx.x] = false;
			threadbusy[threadIdx.x] = false;
			
			nhits_x = 0;
			
			if(hitpairs_x2[bin2+nbins_st2*i_x2].first>=0){
				detID[nhits_x] = detid_list[0];
				i_hit = hitpairs_x2[bin2+nbins_st2*i_x2].first;
				X[nhits_x] = hits_st2x.pos(i_hit);
				errX[nhits_x] = res_st2x;
				Z[nhits_x] = z_st2x;
				elID[nhits_x] = (short)hits_st2x.chan(i_hit);
				tdc[nhits_x] = hits_st2x.tdc(i_hit);
				drift[nhits_x] = hits_st2x.drift(i_hit);
				nhits_x++;
			}
			if(hitpairs_x2[bin2+nbins_st2*i_x2].second>=0){
				detID[nhits_x] = detid_list[1];
				i_hit = hitpairs_x2[bin2+nbins_st2*i_x2].second;
				X[nhits_x] = hits_st2xp.pos(i_hit);
				errX[nhits_x] = res_st2xp;
				Z[nhits_x] = z_st2xp;
				elID[nhits_x] = (short)hits_st2xp.chan(i_hit);
				tdc[nhits_x] = hits_st2xp.tdc(i_hit);
				drift[nhits_x] = hits_st2xp.drift(i_hit);
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
				tdc[nhits_x] = hits_st3x.tdc(i_hit);
				drift[nhits_x] = hits_st3x.drift(i_hit);
				nhits_x++;
			}
			if(hitpairs_x3[bin3+nbins_st3*i_x3].second>=0){
				detID[nhits_x] = detid_list[3];
				i_hit = hitpairs_x3[bin3+nbins_st3*i_x3].second;
				X[nhits_x] = hits_st3xp.pos(i_hit);
				errX[nhits_x] = res_st3xp;
				Z[nhits_x] = z_st3xp;
				elID[nhits_x] = (short)hits_st3xp.chan(i_hit);
				tdc[nhits_x] = hits_st3xp.tdc(i_hit);
				drift[nhits_x] = hits_st3xp.drift(i_hit);
				nhits_x++;
			}
			
			nhits_x3 = nhits_x-nhits_x2;
			if(nhits_x3==0) continue;

			
			fit_2D_track(nhits_x, X, Z, errX, A_, Ainv_, B_, Par, ParErr, chi2);
			x0 = Par[0];
			tx = Par[1];
			
			if(fabs(x0)>1.05*X0_MAX || fabs(tx)>1.05*TX_MAX)continue;

			n_goodxz++;
			
			//prop matching
			nprop = 0;
			checknext = true;

			for(int n = 0; n<nhits_p1x1; n++){
				//ipos = hits_p1x1.pos(n);
				xExp = tx*z_p1x1+x0;
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
					if(fabs(ipos-xExp)<5.08f){
						nprop++;
						checknext = false;
						break;
					}
				}
			}
			if(checknext){
				for(int n = 0; n<nhits_p1x2; n++){
					ipos = hits_p1x2.pos(n);
					xExp = tx*z_p1x2+x0;
					if(fabs(ipos-xExp)<5.08f){
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
					if(fabs(ipos-xExp)<5.08f){
						nprop++;
						checknext = false;
						break;
					}
				}
			}
			if(nprop==0)continue;

#ifdef TEST
			if(ntkl_per_thread[threadIdx.x]<datasizes::TrackletSizeMax/THREADS_PER_BLOCK){
				//what if we try to fill the large arrays straight from here?...
				tklcoll->setStationID(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], binId);
				tklcoll->setThreadID(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], threadIdx.x);
				tklcoll->setnHits(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], nhits_x);
				tklcoll->setTx(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], tx);
				tklcoll->setX0(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], x0);
				tklcoll->setErrTx(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], ParErr[1]);
				tklcoll->setErrX0(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], ParErr[0]);
				
				for(int n = 0; n<nhits_x; n++){
					tklcoll->setHitDetID(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], n, detID[n]);
					tklcoll->setHitChan(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], n, elID[n]);
					tklcoll->setHitPos(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], n, X[n]);
					tklcoll->setHitDrift(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], n, drift[n]);
					tklcoll->setHitSign(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], n, 0.0);
#ifdef FULLCODE
					tklcoll->setHitTDC(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], n, tdc[n]);
					tklcoll->setHitResidual(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], n, 0.0);
#endif
				}
				
				tklcoll->NTracks[tklmult_idx]++;

#ifdef DEBUG
				if(blockIdx.x==debug::EvRef){
					printf("thread %d n_ %d  i %d bin2 %d bin3 %d\n", threadIdx.x, ntkl_per_thread[threadIdx.x], i, bin2, bin3);
					printf("thread %d n_ %d offset %d bin %d x0 %1.4f tx %1.4f, nhits %d \n", threadIdx.x, ntkl_per_thread[threadIdx.x], tkl_coll_offset+array_thread_offset, binId, x0, tx, nhits_x);
					//for(int m = 0; m<nhits_x; m++){
					//	printf("thread %d bin %d hit %d det %d chan %d pos %1.4f tdc %1.4f\n", threadIdx.x, binId, m, detID[m], elID[m], X[m], tdc[m]); 
					//}
				}
				if(ntkl_per_thread[threadIdx.x]>=datasizes::TrackletSizeMax)printf("block %d thread %d bin %d ntkl_per_thread %d \n", blockIdx.x, threadIdx.x, binId, ntkl_per_thread[threadIdx.x]);
				if(nhits_x >= datasizes::MaxHitsPerTrack)printf("block %d thread %d bin %d nhits x %d \n", blockIdx.x, threadIdx.x, binId, nhits_x);
#endif
				ntkl_per_thread[threadIdx.x]++;
			}else{
#endif			
				//we can probably afford to spare time for synchronization here since XZ is extremely fast!
				addtrack[threadIdx.x] = true;
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("thread %d wants to add a track!\n", threadIdx.x);
#endif
				ntkl_min = 100000;
				thread_min[threadIdx.x] = -1;
				nthreads_full = 0;
				nslots_available = datasizes::TrackletSizeMax;
				__syncthreads();
				for(int k = 0; k<THREADS_PER_BLOCK; k++){
					if(ntkl_min>ntkl_per_thread[k] && st3==3+k%2){
						ntkl_min = ntkl_per_thread[k];
						thread_min[threadIdx.x] = k;
						//threadIdx.x*datasizes::TrackletSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
					}
					nslots_available-= ntkl_per_thread[k];
					
					if(addtrack[k]){
#ifdef DEBUG
						if(blockIdx.x==debug::EvRef)printf("thread? %d\n", k);
#endif
						list_of_threads[nthreads_full] = k;
						nthreads_full++;
					}
				}
				if(nslots_available<nthreads_full){
					printf("WARNING: Block %d cannot store anymore tracks! ending the event with a high multiplicity flag! \n", blockIdx.x);
					hastoomanyhits[blockIdx.x] = true;
					break;
				}
				__syncthreads();
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("number of threads requesting to add a track (thread %d): %d  \n", threadIdx.x, nthreads_full);
#endif
				//first, assign a unique "thread min" for each full thread so that they don't step on each other...
				// we have the following info:
				//   the list of threads that want to add a track in another thread,
				//   the number of tracks in each thread
				//   the min_thread
				threadbusy[thread_min[list_of_threads[0]]] = true;
				for(int l = 1; l<nthreads_full; l++){
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("(%d) %d %d %d =? %d \n", threadIdx.x, l, list_of_threads[l], thread_min[list_of_threads[l]], thread_min[list_of_threads[l-1]]);
#endif
					if(thread_min[list_of_threads[l]]==thread_min[list_of_threads[l-1]] ||
						thread_min[list_of_threads[l]]==thread_min[list_of_threads[0]] ||
						st3!=3+thread_min[list_of_threads[l]]%2 ){
						threadbusy[thread_min[list_of_threads[l]]] = true;
#ifdef DEBUG
						if(blockIdx.x==debug::EvRef)printf("(thread %d) thread %d busy\n", threadIdx.x, thread_min[list_of_threads[l]]);
#endif
					}
					
					if(threadbusy[thread_min[threadIdx.x]]){
						ntkl_min = 100000;
						for(int k = 0; k<THREADS_PER_BLOCK; k++){
							if(ntkl_min>ntkl_per_thread[k] && !threadbusy[k]){
								ntkl_min = ntkl_per_thread[k];
								thread_min[list_of_threads[l]] = k;
								threadbusy[k] = true;
							}
						}
					}
				}
				__syncthreads();
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("alt thread chosen for thread %d: %d \n", threadIdx.x, thread_min[threadIdx.x]);
#endif
				array_offset = thread_min[threadIdx.x]*datasizes::TrackletSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;

#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("actual thread %d store thread %d offset %d stid %d local bin %d, bin0 st2 %d, bin0 st3 %d, st3 %d \n", threadIdx.x, thread_min[threadIdx.x], tkl_coll_offset+array_offset, binId, i, bin0_st2, bin0_st3, st3);
#endif
				tklcoll->setStationID(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], (float)binId);
				tklcoll->setThreadID(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], (float)threadIdx.x);
				tklcoll->setnHits(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], (float)nhits_x);
				tklcoll->setTx(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], tx);
				tklcoll->setX0(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], x0);
				tklcoll->setErrTx(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], ParErr[1]);
				tklcoll->setErrX0(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], ParErr[0]);
				
				for(int n = 0; n<nhits_x; n++){
					tklcoll->setHitDetID(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, detID[n]);
					tklcoll->setHitChan(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, elID[n]);
					tklcoll->setHitPos(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, X[n]);
					tklcoll->setHitDrift(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, drift[n]);
					tklcoll->setHitSign(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, 0.0);
#ifdef FULLCODE					
					tklcoll->setHitTDC(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, tdc[n]);
					tklcoll->setHitResidual(tkl_coll_offset+array_offset, ntkl_per_thread[thread_min[threadIdx.x]], n, 0.0);
#endif
				}
				ntkl_per_thread[thread_min[threadIdx.x]]++;
				tklcoll->NTracks[blockIdx.x*THREADS_PER_BLOCK+thread_min[threadIdx.x]]++;
				__syncthreads();
#ifdef TEST
			}
#endif
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
	if(ntkl_per_thread[threadIdx.x]>datasizes::TrackletSizeMax/THREADS_PER_BLOCK){
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
	if(N_tracklets>=datasizes::TrackletSizeMax){
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

//is it even necessary???
///////////////////////////////////////////////////////
//
// tracklets reordering/spreading over threads:
//
///////////////////////////////////////////////////////


__global__ void gKernel_SpreadTracksOverThreads(
	gEventTrackCollection* tklcoll,
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
	// const int bin0_st2 = (threadIdx.x/8)*nbins_st2;
	// const int bin0_st3 = (threadIdx.x%8/2)*nbins_st3;
	const int st3 = 3+threadIdx.x%2;//check d3p for even threads, d3m for odd threads...
	
#ifdef DEBUG
	if(blockIdx.x==debug::EvRef)printf(" thread %d bin0 st2 %d bin0 st3 %d \n", threadIdx.x, bin0_st2, bin0_st3);
#endif

	short hitidx1[100];
	short hitidx2[100];
	short hitflag1[100];
	short hitflag2[100];
	
	// As a starting point we will process 8 bins per thread: x acceptance divided by 4 * 2 bins for d3p and d3m 
	//int nhitpairs_u2[nbins_st2];
	//int nhitpairs_u3[nbins_st3];
        //pairs in station 2
        //thrust::pair<int, int> hitpairs_u2[nbins_st2*geometry::MaxHitsProj[1]];
        //pairs in station 3
        //thrust::pair<int, int> hitpairs_u3[nbins_st3*geometry::MaxHitsProj[1]];

	//int nhitpairs_v2[nbins_st2];
	//int nhitpairs_v3[nbins_st3];
        //pairs in station 2
        //thrust::pair<int, int> hitpairs_v2[nbins_st2*geometry::MaxHitsProj[2]];
        //pairs in station 3
        //thrust::pair<int, int> hitpairs_v3[nbins_st3*geometry::MaxHitsProj[2]];

        //pairs in station 2
        thrust::pair<int, int> hitpairs_u2[100];
        //pairs in station 3
        thrust::pair<int, int> hitpairs_u3[100];

        //pairs in station 2
        thrust::pair<int, int> hitpairs_v2[100];
        //pairs in station 3
        thrust::pair<int, int> hitpairs_v3[100];
	
	
	unsigned int offset_hitcoll;
	
	short stid, projid, detid, detoff;
	
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
	const float res_st2u = planes->spacing[detid];
	
	detid-= 1;
	detid_list[1] = detid;
	det_spacing[1] = planes->spacing[detid];
	int nhits_st2up;
	const gHits hits_st2up = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2up);
	const float z_st2up = planes->z[detid];
	const float res_st2up = planes->spacing[detid];
	
	//make_hitpairs_in_station_bins(hits_st2u, nhits_st2u, hits_st2up, nhits_st2up, hitpairs_u2, nhitpairs_u2, bin0_st2, nbins_st2, hitflag1, hitflag2, stid, projid);
	
	projid = 2;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[2] = detid;
	det_spacing[2] = planes->spacing[detid];
	int nhits_st2v;
	const gHits hits_st2v = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2v);
	const float z_st2v = planes->z[detid];
	const float res_st2v = planes->spacing[detid];
	
	detid-= 1;
	detid_list[3] = detid;
	det_spacing[3] = planes->spacing[detid];
	int nhits_st2vp;
	const gHits hits_st2vp = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2vp);
	const float z_st2vp = planes->z[detid];
	const float res_st2vp = planes->spacing[detid];
	
	//make_hitpairs_in_station_bins(hits_st2v, nhits_st2v, hits_st2vp, nhits_st2vp, hitpairs_v2, nhitpairs_v2, bin0_st2, nbins_st2, hitflag1, hitflag2, stid, projid);
		
	projid = 1;
	stid = st3;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[4] = detid;
	det_spacing[4] = planes->spacing[detid];
	int nhits_st3u;
	const gHits hits_st3u = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3u);
	const float z_st3u = planes->z[detid];
	const float res_st3u = planes->spacing[detid];

	detid-= 1;
	int nhits_st3up;
	detid_list[5] = detid;
	det_spacing[5] = planes->spacing[detid];
	const gHits hits_st3up = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3up);
	const float z_st3up = planes->z[detid];
	const float res_st3up = planes->spacing[detid];

	//make_hitpairs_in_station_bins(hits_st3u, nhits_st3u, hits_st3up, nhits_st3up, hitpairs_u3, nhitpairs_u3, bin0_st3, nbins_st3, hitflag1, hitflag2, stid, projid);
	
	projid = 2;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[6] = detid;
	det_spacing[6] = planes->spacing[detid];
	int nhits_st3v;
	const gHits hits_st3v = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3v);
	const float z_st3v = planes->z[detid];
	const float res_st3v = planes->spacing[detid];

	detid-= 1;
	int nhits_st3vp;
	detid_list[7] = detid;
	det_spacing[7] = planes->spacing[detid];
	const gHits hits_st3vp = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3vp);
	const float z_st3vp = planes->z[detid];
	const float res_st3vp = planes->spacing[detid];

	//make_hitpairs_in_station_bins(hits_st3v, nhits_st3v, hits_st3vp, nhits_st3vp, hitpairs_v3, nhitpairs_v3, bin0_st3, nbins_st3, hitflag1, hitflag2, stid, projid);
	
	/*
	bool bin_overflows = false;
	for(int bin = 0; bin<nbins_st2; bin++){
		if(nhitpairs_u2[bin]>geometry::MaxHitsProj[1])bin_overflows = true;
		if(nhitpairs_v2[bin]>geometry::MaxHitsProj[2])bin_overflows = true;
#ifdef DEBUG
		if(nhitpairs_u2[bin])printf("evt %d bin %d nhits u2 = %d\n", eventID[blockIdx.x], bin, nhitpairs_u2[bin]);
		if(nhitpairs_v2[bin])printf("evt %d bin %d nhits v2 = %d\n", eventID[blockIdx.x], bin, nhitpairs_v2[bin]);
#endif
	}
	for(int bin = 0; bin<nbins_st3; bin++){
		if(nhitpairs_u3[bin]>geometry::MaxHitsProj[1])bin_overflows = true;
		if(nhitpairs_v3[bin]>geometry::MaxHitsProj[2])bin_overflows = true;
#ifdef DEBUG
		if(nhitpairs_u3[bin])printf("evt %d bin %d nhits u3 = %d\n", eventID[blockIdx.x], bin, nhitpairs_u3[bin]);
		if(nhitpairs_v3[bin])printf("evt %d bin %d nhits v3 = %d\n", eventID[blockIdx.x], bin, nhitpairs_v3[bin]);
#endif
	}
	if(bin_overflows){
		hastoomanyhits[blockIdx.x] = true;
		return;
	}
	*/
	
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
	
	//Arrays for other basic hit info
	float elID[8];
	float pos[8];
	float drift[8];
	//float tdc[8];

	short nhits_uv;
	
	short nhits_u2, nhits_u3;
	short nhits_v2, nhits_v3;
	
	// if parallel:
	// const short nbins_st3 = geometry::N_WCHitsBins[2];
	const short nbins_total = nbins_st2*nbins_st3;
	short bin2, bin3;
	
	int nu2, nv2, nu3, nv3;
	
	int ncomb_uv3, ncomb_uv2;
	
	short i_u2, i_u3, i_v2, i_v3;
	
	int i_uv3, i_uv2, i_hit;

	float y, err_y;
	
	float chi2min = 10000.1f;

	int localbin;
	
	//TODO get Hodo hits stations 2, 3, 4...
	
	//get the tracks...
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackletSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackletSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
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
#ifdef DEBUG
		if(threadIdx.x!=(int)Tracks.threadID(i)){
			printf("does that actually ever happen???? ");
		}
#endif
		trackthread = (int)Tracks.threadID(i);
		if(blockIdx.x==debug::EvRef)printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f nhits %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f invP %1.4f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.nHits(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.invP(i));
		//tkl = Tracks.getTracklet(i);
		x0 = Tracks.x0(i);
		tx = Tracks.tx(i);
		err_x0 = Tracks.err_x0(i);
		err_tx = Tracks.err_tx(i);
		nhits_x =  Tracks.nHits(i);
		update_track = false;
		
		//get the corresponding local bin!
		//binId = threadIdx.x+THREADS_PER_BLOCK*localbin;
		//localbin = (tkl.stationID-tkl.threadID)/THREADS_PER_BLOCK;
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
		
		//x0 = tkl.x0;
		//tx = tkl.tx;
		ncomb_uv3 = nu3*nv3;
		ncomb_uv2 = nu2*nv2;

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
				//tdc[nhits_uv] = hits_st3u.tdc(i_hit);
				drift[nhits_uv] = hits_st3u.drift(i_hit);
				pos[nhits_uv] = hits_st3u.pos(i_hit);
				calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
				Y[nhits_uv] = y;
				errY[nhits_uv] = err_y;
				nhits_uv++;
				ty+= y/z_st3u;
			}
			if(hitpairs_u3[i_u3].second>=0){
				detID[nhits_uv] = detid_list[5];
				i_hit = hitpairs_u3[i_u3].second;
				Z[nhits_uv] = z_st3up;
				elID[nhits_uv] = (short)hits_st3up.chan(i_hit);
				//tdc[nhits_uv] = hits_st3up.tdc(i_hit);
				drift[nhits_uv] = hits_st3up.drift(i_hit);
				pos[nhits_uv] = hits_st3up.pos(i_hit);
				calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
				Y[nhits_uv] = y;
				errY[nhits_uv] = err_y;
				nhits_uv++;
				ty+= y/z_st3up;
			}
			
			nhits_u3 = nhits_uv;
			if(nhits_u3==0) continue;
			
			if(hitpairs_v3[i_v3].first>=0){
				detID[nhits_uv] = detid_list[6];
				i_hit = hitpairs_v3[i_v3].first;
				Z[nhits_uv] = z_st3v;
				elID[nhits_uv] = (short)hits_st3v.chan(i_hit);
				//tdc[nhits_uv] = hits_st3v.tdc(i_hit);
				drift[nhits_uv] = hits_st3v.drift(i_hit);
				pos[nhits_uv] = hits_st3v.pos(i_hit);
				calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
				Y[nhits_uv] = y;
				errY[nhits_uv] = err_y;
				nhits_uv++;
				ty+= y/z_st3v;
			}
			if(hitpairs_v3[i_v3].second>=0){
				detID[nhits_uv] = detid_list[7];
				i_hit = hitpairs_v3[i_v3].second;
				Z[nhits_uv] = z_st3vp;
				elID[nhits_uv] = (short)hits_st3vp.chan(i_hit);
				//tdc[nhits_uv] = hits_st3vp.tdc(i_hit);
				drift[nhits_uv] = hits_st3vp.drift(i_hit);
				pos[nhits_uv] = hits_st3vp.pos(i_hit);
				calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
				Y[nhits_uv] = y;
				errY[nhits_uv] = err_y;
				nhits_uv++;
				ty+= y/z_st3vp;
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
					//tdc[nhits_uv] = hits_st2u.tdc(i_hit);
					drift[nhits_uv] = hits_st2u.drift(i_hit);
					pos[nhits_uv] = hits_st2u.pos(i_hit);
					calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
					//if( (y-ty*z_st2u)/err_y<3.f ){//since ty is *very* rough let's be generous...				
						Y[nhits_uv] = y;
						errY[nhits_uv] = err_y;
						nhits_uv++;
					//}
				}
				if(hitpairs_u2[i_u2].second>=0){
					detID[nhits_uv] = detid_list[1];
					i_hit = hitpairs_u2[i_u2].second;
					Z[nhits_uv] = z_st2up;
					elID[nhits_uv] = (short)hits_st2up.chan(i_hit);
					//tdc[nhits_uv] = hits_st2up.tdc(i_hit);
					drift[nhits_uv] = hits_st2up.drift(i_hit);
					pos[nhits_uv] = hits_st2up.pos(i_hit);
					calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
					//if( fabs(y-ty*z_st2up)/err_y<3.f ){//since ty is *very* rough let's be generous...				
						Y[nhits_uv] = y;
						errY[nhits_uv] = err_y;
						nhits_uv++;
					//}
				}
				
				nhits_u2 = nhits_uv-nhits_u3-nhits_v3;
				if(nhits_u2==0) continue;
				
				if(hitpairs_v2[i_v2].first>=0){
					detID[nhits_uv] = detid_list[2];
					i_hit = hitpairs_v2[i_v2].first;
					Z[nhits_uv] = z_st2v;
					elID[nhits_uv] = (short)hits_st2v.chan(i_hit);
					//tdc[nhits_uv] = hits_st2v.tdc(i_hit);
					drift[nhits_uv] = hits_st2v.drift(i_hit);
					pos[nhits_uv] = hits_st2v.pos(i_hit);
					calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
					//if( fabs(y-ty*z_st2v)/err_y<3.f ){//since ty is *very* rough let's be generous...				
						Y[nhits_uv] = y;
						errY[nhits_uv] = err_y;
						nhits_uv++;
					//}
				}
				if(hitpairs_v2[i_v2].second>=0){
					detID[nhits_uv] = detid_list[3];
					i_hit = hitpairs_v2[i_v2].second;
					Z[nhits_uv] = z_st2vp;
					elID[nhits_uv] = (short)hits_st2vp.chan(i_hit);
					//tdc[nhits_uv] = hits_st2vp.tdc(i_hit);
					drift[nhits_uv] = hits_st2vp.drift(i_hit);
					pos[nhits_uv] = hits_st2vp.pos(i_hit);
					calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_uv], elID[nhits_uv], pos[nhits_uv], drift[nhits_uv]);
#endif
					//if( fabs(y-ty*z_st2vp)/err_y<3.f ){//since ty is *very* rough let's be generous...				
						Y[nhits_uv] = y;
						errY[nhits_uv] = err_y;
						nhits_uv++;
					//}
				}
				
				nhits_v2 = nhits_uv-nhits_u3-nhits_v3-nhits_u2;
				if(nhits_v2==0) continue;
				
				fit_2D_track(nhits_uv, Y, Z, errY, A_, Ainv_, B_, Par, ParErr, chi2);
				
				y0 = Par[0];
				ty = Par[1];
				
				//TODO: hodoscope matching
				
				//TODO: chi2 evaluation of track candidate
				
				//TODO: LR ambiguity resolution
				
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef){
					printf("thread %d bin %d y0 %1.4f ty %1.4f, nhits %d, u3 %d, v3 %d, u2 %d, v2 %d\n", threadIdx.x, localbin, y0, ty, nhits_x, nhits_u3, nhits_v3, nhits_u2, nhits_v2);
				}
#endif
				//if(chi2>chi2min && ){
				//chi2min = chi2
				update_track = true;
				besttrackYZdata[threadIdx.x][2] = nhits_uv;
				besttrackYZdata[threadIdx.x][3] = chi2;
				besttrackYZdata[threadIdx.x][6] = ty;
				besttrackYZdata[threadIdx.x][8] = y0;
				besttrackYZdata[threadIdx.x][11] = ParErr[1];
				besttrackYZdata[threadIdx.x][13] = ParErr[0];
				
				for(int n = 0; n<nhits_uv;n++){
				besttrackYZdata[threadIdx.x][16+n] = detID[n];
				besttrackYZdata[threadIdx.x][34+n] = elID[n];
				besttrackYZdata[threadIdx.x][52+n] = pos[n];
				besttrackYZdata[threadIdx.x][70+n] = drift[n];
				besttrackYZdata[threadIdx.x][88+n] = 0;
#ifdef FULLCODE
				besttrackYZdata[threadIdx.x][106+n] = tdc[n];
				besttrackYZdata[threadIdx.x][124+n] = 0;
#endif
				}
				//}
				
			}
		}
		
		if(update_track){
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 5);
			nhits_x = Tracks.nHits(i);
			nhits_uv = besttrackYZdata[threadIdx.x][2];
			tklcoll->setnHits(tkl_coll_offset+array_thread_offset, i, nhits_x+nhits_uv);
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
	detid_list[0] = detid;
	int nhits_st1x;
	const gHits hits_st1x = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1x);
	const float z_st1x = planes->z[detid];
	const float res_st1x = planes->spacing[detid];
	
	detid-= 1;
	detid_list[1] = detid;
	int nhits_st1xp;
	const gHits hits_st1xp = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1xp);
	const float z_st1xp = planes->z[detid];
	const float res_st1xp = planes->spacing[detid];
	
	projid = 1;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[2] = detid;
	int nhits_st1u;
	const gHits hits_st1u = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1u);
	const float z_st1u = planes->z[detid];
	const float res_st1u = planes->spacing[detid];
	
	detid-= 1;
	detid_list[3] = detid;
	int nhits_st1up;
	const gHits hits_st1up = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1up);
	const float z_st1up = planes->z[detid];
	const float res_st1up = planes->spacing[detid];
	
	projid = 2;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[4] = detid;
	int nhits_st1v;
	const gHits hits_st1v = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1v);
	const float z_st1v = planes->z[detid];
	const float res_st1v = planes->spacing[detid];
	
	detid-= 1;
	detid_list[5] = detid;
	int nhits_st1vp;
	const gHits hits_st1vp = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1vp);
	const float z_st1vp = planes->z[detid];
	const float res_st1vp = planes->spacing[detid];

	// TODO: load the hodoscope hits
	//detid = nChamberPlanes+1;
	//int nhits_h1x1;
	//const gHits hits_h1x1 = hitcolls->hitshodo(blockIdx.x, detid, nhits_h1x1);
	
	
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
	float elID[6];
	float pos[6];
	float drift[6];
	//float tdc[6];

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
	
	//get the tracks...
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackletSizeMax*datasizes::NTracksParam;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackletSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	
	__shared__ float besttrackdata[THREADS_PER_BLOCK][datasizes::NTracksParam];
	
	unsigned int Ntracks;
	const gTracks Tracks = tklcoll->tracks(blockIdx.x, threadIdx.x, Ntracks);
	
	for(int i = 0; i<Ntracks; i++){
		//tkl = Tracks.getTracklet(i);
		update_track = false;
		nhits_st23 = Tracks.nHits(i);
		
		if(Tracks.stationID(i)<5)continue;
		
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
		if(blockIdx.x==debug::EvRef)printf("x0 %1.4f y0 %1.4f tx %1.4f ty %1.4f nhits %1.0f detid %d, nx1 %d nu1 %d nv1 %d \n", x0, y0, tx, ty, Tracks.nHits(i), Tracks.get_lasthitdetid(i), nx1, nu1, nv1);
#endif
				
		if(nx1==0 || nu1==0 || nv1==0)continue;
		
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
				//tdc[nhits_x] = hits_st1x.tdc(i_hit);
				drift[nhits_x] = hits_st1x.drift(i_hit);
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_x], elID[nhits_x], pos[nhits_x], drift[nhits_x]);
#endif
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
				//tdc[nhits_x] = hits_st1xp.tdc(i_hit);
				drift[nhits_x] = hits_st1xp.drift(i_hit);
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_x], elID[nhits_x], pos[nhits_x], drift[nhits_x]);
#endif
				nhits_x++;
			}
			
			if(nhits_x==0) continue;
			
			X[nhits_x] = x0+tx*geometry::Z_KMAG_BEND;
			errX[nhits_x] = errx0+errtx*geometry::Z_KMAG_BEND;
			Z[nhits_x] = geometry::Z_KMAG_BEND;
			
			fit_2D_track(nhits_x+1, X, Z, errX, A_, Ainv_, B_, Par, ParErr, chi2);
			x0_st1 = Par[0];
			tx_st1 = Par[1];

			errx0_st1 = Par[0];
			errtx_st1 = Par[1];
			
			invP = calculate_invP_charge(tx, tx_st1, charge);
			errinvP = calculate_invP_error(errtx, errtx_st1);

#ifdef DEBUG
			if(blockIdx.x==debug::EvRef){
				printf("thread %d tx_st1 %1.4f tx %1.4f invP %1.4f  %1.4f charge %d \n", threadIdx.x,  tx_st1, tx, (tx_st1 - tx) / geometry::PT_KICK_KMAG, invP, charge);
			}
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
					//tdc[nhits_x+nhits_uv] = hits_st1u.tdc(i_hit);
					drift[nhits_x+nhits_uv] = hits_st1u.drift(i_hit);
					pos[nhits_x+nhits_uv] = hits_st1u.pos(i_hit);
					calculate_y_uvhit(detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], drift[nhits_x+nhits_uv], 0, x0, tx, planes, y, err_y);
					//if( (y-y_trk(y0, ty, z_st1u))*(y-y_trk(y0, ty, z_st1u))/(err_y*err_y+err_y_trk(erry0, errty, z_st1u)*err_y_trk(erry0, errty, z_st1u))<2.0 )
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], pos[nhits_x+nhits_uv], drift[nhits_x+nhits_uv]);
#endif
					nhits_uv++;
				}
				if(hitpairs_u1[i_u].second>=0){
					detID[nhits_x+nhits_uv] = detid_list[3];
					i_hit = hitpairs_u1[i_u].second;
					Z[nhits_x+nhits_uv] = z_st1up;
					elID[nhits_x+nhits_uv] = (short)hits_st1up.chan(i_hit);
					//tdc[nhits_x+nhits_uv] = hits_st1up.tdc(i_hit);
					drift[nhits_x+nhits_uv] = hits_st1up.drift(i_hit);
					pos[nhits_x+nhits_uv] = hits_st1up.pos(i_hit);
					calculate_y_uvhit(detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], drift[nhits_x+nhits_uv], 0, x0, tx, planes, y, err_y);
					//if( (y-y_trk(y0, ty, z_st1up))*(y-y_trk(y0, ty, z_st1up))/(err_y*err_y+err_y_trk(erry0, errty, z_st1up)*err_y_trk(erry0, errty, z_st1up))<2.0 )
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], pos[nhits_x+nhits_uv], drift[nhits_x+nhits_uv]);
#endif
					nhits_uv++;
				}
				nhits_u = nhits_uv;
				if(nhits_u==0)continue;
				
				if(hitpairs_v1[i_v].first>=0){
					detID[nhits_x+nhits_uv] = detid_list[4];
					i_hit = hitpairs_v1[i_v].first;
					Z[nhits_x+nhits_uv] = z_st1v;
					elID[nhits_x+nhits_uv] = (short)hits_st1v.chan(i_hit);
					//tdc[nhits_x+nhits_uv] = hits_st1v.tdc(i_hit);
					drift[nhits_x+nhits_uv] = hits_st1v.drift(i_hit);
					pos[nhits_x+nhits_uv] = hits_st1v.pos(i_hit);
					calculate_y_uvhit(detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], drift[nhits_x+nhits_uv], 0, x0, tx, planes, y, err_y);
					//if( (y-y_trk(y0, ty, z_st1v))*(y-y_trk(y0, ty, z_st1v))/(err_y*err_y+err_y_trk(erry0, errty, z_st1v)*err_y_trk(erry0, errty, z_st1v))<2.0 )
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], pos[nhits_x+nhits_uv], drift[nhits_x+nhits_uv]);
#endif
					nhits_uv++;
				}
				if(hitpairs_v1[i_v].second>=0){
					detID[nhits_x+nhits_uv] = detid_list[5];
					i_hit = hitpairs_v1[i_v].second;
					Z[nhits_x+nhits_uv] = z_st1vp;
					elID[nhits_x+nhits_uv] = (short)hits_st1vp.chan(i_hit);
					//tdc[nhits_x+nhits_uv] = hits_st1vp.tdc(i_hit);
					drift[nhits_x+nhits_uv] = hits_st1vp.drift(i_hit);
					pos[nhits_x+nhits_uv] = hits_st1vp.pos(i_hit);
					calculate_y_uvhit(detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], drift[nhits_x+nhits_uv], 0, x0, tx, planes, y, err_y);
					//if( (y-y_trk(y0, ty, z_st1vp))*(y-y_trk(y0, ty, z_st1vp))/(err_y*err_y+err_y_trk(erry0, errty, z_st1vp)*err_y_trk(erry0, errty, z_st1vp))<2.0 )
#ifdef DEBUG
					if(blockIdx.x==debug::EvRef)printf("det %d chan %1.0f pos %1.4f drift %1.4f\n", detID[nhits_x+nhits_uv], elID[nhits_x+nhits_uv], pos[nhits_x+nhits_uv], drift[nhits_x+nhits_uv]);
#endif
					nhits_uv++;
				}
				nhits_v = nhits_uv-nhits_u;
				if(nhits_v==0)continue;
				
				//TODO: matching hodoscope, 
				
				
				//TODO: resolve left right...
				
				
				//TODO: chi2 fit...

#ifdef DEBUG
				if(blockIdx.x==debug::EvRef){
					printf("thread %d invP %1.4f charge %d nhits %d \n", threadIdx.x, invP, charge, nhits_st23);
				}
#endif				
				//if(chi2>chi2min){
				//chi2min = chi2
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
				besttrackdata[threadIdx.x][88+m] = 0;
#ifdef FULLCODE
				besttrackdata[threadIdx.x][105+m] = tdc[m];
				besttrackdata[threadIdx.x][124+m] = 0;
#endif
				}
				//}

				
			}
			
		}//end
		
		if(update_track){
			tklcoll->setStationID(tkl_coll_offset+array_thread_offset, i, 6);
			nhits_st1 = besttrackdata[threadIdx.x][2];
			tklcoll->setnHits(tkl_coll_offset+array_thread_offset, i, nhits_st23+nhits_st1);
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



// simple track printing function for debugging. 

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



#ifdef OLDCODE

// function to match a tracklet to a hodoscope hit
__device__ int match_tracklet_to_hodo(const gTracklet tkl, const int stID, gEvent* ic, const gPlane* planes)
{
	const int index = threadIdx.x + blockIdx.x * blockDim.x;
	const int idxoff_global = index*EstnAHMax;
	int masked = 0;//false
	// first, define the search region, and foremost, the planes on which we define this search region, which depends on the station ID we're looking at
	// define the region in which we are supposed to have hits:

	//printf(" stID %d hodo plane[0] %d [1] %d \n", stID, geometry::hodoplanerange[stID][0], geometry::hodoplanerange[stID][1]);
	//printf(" x0 %1.4f +- %1.4f, y0 %1.4f +- %1.4f, tx %1.4f +- %1.4f, ty %1.4f +- %1.4f \n", tkl.x0, tkl.err_x0, tkl.y0, tkl.err_y0, tkl.tx, tkl.err_tx, tkl.ty, tkl.err_ty);		
	
	float xhodo, yhodo, err_x, err_y, xmin, xmax, ymin, ymax;
		
	// loop on the hits and select hodoscope hits corresponding to the station
	for(int i = 0; i<ic->nAH[index]; i++){
		//we only consider hits in the hodoscopes planes corresponding to the station where the tracklet is reconstructed 
		if(geometry::hodoplanerange[stID][0]>ic->AllHits[idxoff_global+i].detectorID || geometry::hodoplanerange[stID][1]<ic->AllHits[idxoff_global+i].detectorID)continue;
		
		//calculate the track position at the hodoscope plane z
		xhodo = planes->z[ic->AllHits[idxoff_global+i].detectorID]*tkl.tx+tkl.x0;
		yhodo = planes->z[ic->AllHits[idxoff_global+i].detectorID]*tkl.ty+tkl.y0;
		
		err_x = 3.f*(fabs(planes->z[ic->AllHits[idxoff_global+i].detectorID]*tkl.err_tx)+tkl.err_x0);
		err_y = 3.f*(fabs(planes->z[ic->AllHits[idxoff_global+i].detectorID]*tkl.err_ty)+tkl.err_y0);

		//printf(" det %d elem %d z_hodo %1.4f x_hodo %1.4f y_hodo %1.4f err x %1.4f err y %1.4f \n", ic->AllHits[idxoff_global+i].detectorID, ic->AllHits[idxoff_global+i].elementID, planes[ic->AllHits[idxoff_global+i].detectorID].z, xhodo, yhodo, err_x, err_y);

		//calculate "xmin, xmax, ymin, ymax" in which the track is supposed to pass through; 
		//these are basicially defined as the spatial coverage of the hit hodoscope element (plus some fudge factor for x)
		if(planes->costheta[ic->AllHits[idxoff_global+i].detectorID]>0.99){
			xmin = ic->AllHits[idxoff_global+i].pos-planes->cellwidth[ic->AllHits[idxoff_global+i].detectorID]*0.5f;
			xmax = ic->AllHits[idxoff_global+i].pos+planes->cellwidth[ic->AllHits[idxoff_global+i].detectorID]*0.5f;
			
			ymin = planes->y1[ic->AllHits[idxoff_global+i].detectorID];
			ymax = planes->y2[ic->AllHits[idxoff_global+i].detectorID];
			
			xmin-=(xmax-xmin)*geometry::hodofudgefac[stID];
			xmax+=(xmax-xmin)*geometry::hodofudgefac[stID];
			
			ymin-=(ymax-ymin)*geometry::hodofudgefac[stID];
			ymax+=(ymax-ymin)*geometry::hodofudgefac[stID];
			
			ymin+=planes->y0[ic->AllHits[idxoff_global+i].detectorID];
			ymax+=planes->y0[ic->AllHits[idxoff_global+i].detectorID];
		}else{
			xmin = planes->x1[ic->AllHits[idxoff_global+i].detectorID];
			xmax = planes->x2[ic->AllHits[idxoff_global+i].detectorID];
			
			ymin = ic->AllHits[idxoff_global+i].pos-planes->cellwidth[ic->AllHits[idxoff_global+i].detectorID]*0.5f;
			ymax = ic->AllHits[idxoff_global+i].pos+planes->cellwidth[ic->AllHits[idxoff_global+i].detectorID]*0.5f;
			
			xmin-=(xmax-xmin)*geometry::hodofudgefac[stID];
			xmax+=(xmax-xmin)*geometry::hodofudgefac[stID];
			
			ymin-=(ymax-ymin)*geometry::hodofudgefac[stID];
			ymax+=(ymax-ymin)*geometry::hodofudgefac[stID];
		}

		//printf(" xmin %1.4f xmax %1.4f, ymin %1.4f ymax %1.4f xfudge %1.4f\n", xmin, xmax, ymin, ymax, (xmax-xmin)*0.15 );
		err_x+= (xmax-xmin)*0.15;

		xmin-= err_x;
		xmax+= err_x;
		ymin-= err_y;
		ymax+= err_y;
		
		//printf(" xmin %1.4f xmax %1.4f, ymin %1.4f ymax %1.4f \n", xmin, xmax, ymin, ymax);
		if(xmin <= xhodo && xhodo <= xmax && ymin <= yhodo && yhodo <= ymax ){
			masked++;
			break;
		}
		
	}
	// 
	return masked>0;
}




__device__ float refit_backpartialtrack_with_drift(gTracklet& tkl, gStraightFitArrays* fitarray, const gPlane* planes){
	const int index = threadIdx.x + blockIdx.x * blockDim.x;
	const int idxoff_array = index*MaxHitsPerTrack;
	
	float A_[4];
	float Ainv_[4];
	float B_[2];
	float Par[4];
	float ParErr[4];
	float chi2;
	
	//X hits are stored first, so we fit them first;
	for(int i = 0; i<tkl.nXHits; i++){
		FillFitArrays_X(i, tkl.hits[i], tkl.hitsign[i], fitarray, planes);
		FillChi2Arrays(i, tkl.hits[i], tkl.hitsign[i], fitarray, planes);
	}
	fit_2D_track(tkl.nXHits, fitarray->x_array+idxoff_array, fitarray->z_array+idxoff_array, fitarray->dx_array+idxoff_array, A_, Ainv_, B_, Par, ParErr, chi2);

	tkl.x0 = Par[0];
	tkl.err_x0 = ParErr[0];
	
	tkl.tx = Par[1];
	tkl.err_tx = ParErr[1];
	
	float y, err_y;
	
	for(int i = tkl.nXHits; i<tkl.nXHits+tkl.nUHits+tkl.nVHits; i++){
		if( calculate_y_uvhit(y, err_y, tkl.hits[i], tkl.hitsign[i], tkl, planes) ){
			FillFitArrays_UV(i-tkl.nXHits, tkl.hits[i], fitarray, planes, y, err_y);
			FillChi2Arrays(i, tkl.hits[i], tkl.hitsign[i], fitarray, planes);
		}
	}
	fit_2D_track(tkl.nXHits, fitarray->x_array+idxoff_array, fitarray->z_array+idxoff_array, fitarray->dx_array+idxoff_array, A_, Ainv_, B_, Par, ParErr, chi2);
	
	tkl.y0 = Par[0];
	tkl.err_y0 = ParErr[0];
	
	tkl.ty = Par[1];
	tkl.err_ty = ParErr[1];
	
	Par[0] = tkl.x0;
	Par[1] = tkl.y0;
	Par[2] = tkl.tx;
	Par[3] = tkl.ty;

	chi2_straight(tkl.nXHits+tkl.nUHits+tkl.nVHits, fitarray->drift_dist+idxoff_array, fitarray->resolution+idxoff_array,
			fitarray->p1x+idxoff_array, fitarray->p1y+idxoff_array, fitarray->p1z+idxoff_array,
			fitarray->deltapx+idxoff_array, fitarray->deltapy+idxoff_array, fitarray->deltapz+idxoff_array,
			Par, chi2);
	
	tkl.chisq = chi2;
	return chi2;
}


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

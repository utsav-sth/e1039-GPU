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
	if(station_mult[3]>120)hastoomanyhits[blockIdx.x] = true;
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

__global__ void gKernel_XZ_tracking(gEventHitCollections* hitcolls, gEventTrackCollection* tklcoll, const float* z_array, const float* res_array, const int* eventID, bool* hastoomanyhits)
//(gTrackingXZparams* parameters)			
{
	if(hastoomanyhits[blockIdx.x]){
#ifdef DEBUG
		if(threadIdx.x==0)printf("Evt %d discarded, too many hits\n", eventID[blockIdx.x]);
#endif
		return;
	}
	
	const int nbins_st2 = 7;//28/4
	const int nbins_st3 = 29;//58/2
	
	// thread 0: bin0_st2 = 0/2*7, st3 = 3; thread 1: bin0_st2 = (1/2 = 0)*7, st3 = 1; thread 2: bin0_st2 = 2/2*7, st3 = 3; thread 3: bin0_st2 = (3/2 = 1)*7, st3 = 4; 
	const int bin0_st2 = (threadIdx.x/2)*nbins_st2;
	const int bin0_st3 = 0;
	int st3 = 3+threadIdx.x%2;//check d3p for even threads, d3m for odd threads...

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
		//if(nhitpairs_x2[bin])printf("evt %d bin %d nhits x2 = %d\n", ic[index].EventID, bin, nhitpairs_x2[bin]);
	}
	for(int bin = 0; bin<nbins_st3; bin++){
		if(nhitpairs_x3[bin]>geometry::MaxHitsProj[0])bin_overflows = true;
		//if(nhitpairs_x3[bin])printf("evt %d bin %d nhits x3 = %d\n", ic[index].EventID, bin, nhitpairs_x3[bin]);
	}
	if(bin_overflows){
		hastoomanyhits[blockIdx.x] = true;
		return;
	}
	
	//tracklet container://limit the access to the large tracklet arrays.
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackletSizeMax*datasizes::NTracksParam;
	
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

	//static float tkl_data_local[datasizes::TrackletSizeMax][datasizes::NTracksParam];
	//const unsigned int array_thread_offset_2 = threadIdx.x*datasizes::TrackletSizeMax/THREADS_PER_BLOCK;
		
	__shared__ unsigned int ntkl_per_thread[THREADS_PER_BLOCK];
	for(int k = 0; k<THREADS_PER_BLOCK; k++)ntkl_per_thread[k] = 0;
	const unsigned int array_thread_offset = threadIdx.x*datasizes::TrackletSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	//const unsigned int nhits_thread_max = datasizes::TrackletSizeMax/THREADS_PER_BLOCK;
	const unsigned int tklmult_idx = blockIdx.x*THREADS_PER_BLOCK+threadIdx.x;

	tklcoll->NTracks[tklmult_idx] = 0;
	
	for(short i = 0; i<nbins_total; i++){
		bin2 = i%nbins_st2;
		bin3 = (i-bin2)/nbins_st2;
		binId = threadIdx.x+THREADS_PER_BLOCK*i;
		//printf("bin %d %d\n", bin2, bin3);
		// if parallel:
		// bin3+=nbins_st3*i_thread;

		nx2 = nhitpairs_x2[bin2];
		nx3 = nhitpairs_x3[bin3];
		
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
			
			//for(int k = 0; k<THREADS_PER_BLOCK; k++)ntkl[k] = 0;
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

			if(ntkl_per_thread[threadIdx.x]<datasizes::TrackletSizeMax/THREADS_PER_BLOCK){
				//what if we try to fill the large arrays straight from here?...
				tklcoll->setStationID(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], binId);
				tklcoll->setThreadID(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], threadIdx.x);
				tklcoll->setnHits(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], nhits_x);
				tklcoll->setTx(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], tx);
				tklcoll->setX0(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], x0);
				tklcoll->setErrTx(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], ParErr[1]);
				tklcoll->setErrX0(tkl_coll_offset+array_thread_offset, ntkl_per_thread[threadIdx.x], ParErr[0]);
				tklcoll->NTracks[tklmult_idx]++;
#ifdef DEBUG
				if(blockIdx.x==debug::EvRef){
					printf("thread %d n_ %d offset %d bin %d x0 %1.4f tx %1.4f, nhits %d \n", threadIdx.x, ntkl_per_thread[threadIdx.x], tkl_coll_offset+array_thread_offset, binId, x0, tx, nhits_x);
					//for(int m = 0; m<nhits_x; m++){
					//	printf("thread %d bin %d hit %d det %d chan %d pos %1.4f tdc %1.4f\n", threadIdx.x, binId, m, detID[m], elID[m], X[m], tdc[m]); 
					//}
				}
				if(ntkl_per_thread[threadIdx.x]>=datasizes::TrackletSizeMax)printf("block %d thread %d bin %d ntkl_per_thread %d \n", blockIdx.x, threadIdx.x, binId, ntkl_per_thread[threadIdx.x]);
				if(nhits_x >= datasizes::MaxHitsPerTrack)printf("block %d thread %d bin %d nhits x %d \n", blockIdx.x, threadIdx.x, binId, nhits_x);
#endif
			}

			ntkl_per_thread[threadIdx.x]++;
		}// end loop on hits
	}//end loop on bins
	//evaluate number of tracklets

	__syncthreads();
	
	int N_tracklets = 0;

	//__shared__ unsigned int array_thread_offset[THREADS_PER_BLOCK];
	for(int k = 0; k<THREADS_PER_BLOCK; k++){
		N_tracklets+= ntkl_per_thread[k];
		//array_thread_offset[k] = k*datasizes::TrackletSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK;
	}
	if(ntkl_per_thread[threadIdx.x]>datasizes::TrackletSizeMax/THREADS_PER_BLOCK){
		printf("block %d thread %d tracklets per thread: %d \n", blockIdx.x, threadIdx.x, ntkl_per_thread[threadIdx.x]);
		hastoomanyhits[blockIdx.x] = true;
	}

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


////////////////////////////////////
//
//          YZ TRACKING
//
////////////////////////////////////

__global__ void gKernel_YZ_tracking(gEventHitCollections* hitcolls, gEventTrackCollection* tklcoll, const gPlane* planes, int* eventID, bool* hastoomanyhits)
{
	if(hastoomanyhits[blockIdx.x]){
#ifdef DEBUG
		if(threadIdx.x==0)printf("Evt %d discarded, too many hits\n", eventID[blockIdx.x]);
#endif
		return;
	}
	
	const int nbins_st2 = 7;//28/4
	const int nbins_st3 = 29;//58/2
	
	// thread 0: bin0_st2 = 0/2*7, st3 = 3; thread 1: bin0_st2 = (1/2 = 0)*7, st3 = 1; thread 2: bin0_st2 = 2/2*7, st3 = 3; thread 3: bin0_st2 = (3/2 = 1)*7, st3 = 4; 
	const int bin0_st2 = (threadIdx.x/2)*nbins_st2;
	const int bin0_st3 = 0;
	int st3 = 3+threadIdx.x%2;//check d3p for even threads, d3m for odd threads...

	short hitflag1[100];
	short hitflag2[100];
	
	// As a starting point we will process 8 bins per thread: x acceptance divided by 4 * 2 bins for d3p and d3m 
	int nhitpairs_u2[nbins_st2];
	int nhitpairs_u3[nbins_st3];
        //pairs in station 2
        thrust::pair<int, int> hitpairs_u2[nbins_st2*geometry::MaxHitsProj[1]];
        //pairs in station 3
        thrust::pair<int, int> hitpairs_u3[nbins_st3*geometry::MaxHitsProj[1]];

	int nhitpairs_v2[nbins_st2];
	int nhitpairs_v3[nbins_st3];
        //pairs in station 2
        thrust::pair<int, int> hitpairs_v2[nbins_st2*geometry::MaxHitsProj[2]];
        //pairs in station 3
        thrust::pair<int, int> hitpairs_v3[nbins_st3*geometry::MaxHitsProj[2]];

	
	unsigned int offset_hitcoll;
	
	short stid, projid, detid, detoff;
	
	projid = 1;
	stid = 2;
	detoff = 1;
	short detid_list[8];
	
	//Get the required hit collections here!
	
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[0] = detid;
	int nhits_st2u;
	const gHits hits_st2u = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2u);
	const float z_st2u = planes->z[detid];
	const float res_st2u = planes->spacing[detid];
	
	detid-= 1;
	detid_list[1] = detid;
	int nhits_st2up;
	const gHits hits_st2up = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2up);
	const float z_st2up = planes->z[detid];
	const float res_st2up = planes->spacing[detid];
	
	make_hitpairs_in_station_bins(hits_st2u, nhits_st2u, hits_st2up, nhits_st2up, hitpairs_u2, nhitpairs_u2, bin0_st2, nbins_st2, hitflag1, hitflag2, stid, projid);
	
	projid = 2;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[2] = detid;
	int nhits_st2v;
	const gHits hits_st2v = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2v);
	const float z_st2v = planes->z[detid];
	const float res_st2v = planes->spacing[detid];
	
	detid-= 1;
	detid_list[3] = detid;
	int nhits_st2vp;
	const gHits hits_st2vp = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2vp);
	const float z_st2vp = planes->z[detid];
	const float res_st2vp = planes->spacing[detid];
	
	make_hitpairs_in_station_bins(hits_st2v, nhits_st2v, hits_st2vp, nhits_st2vp, hitpairs_v2, nhitpairs_v2, bin0_st2, nbins_st2, hitflag1, hitflag2, stid, projid);
		
	projid = 1;
	stid = st3;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[4] = detid;
	int nhits_st3u;
	const gHits hits_st3u = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3u);
	const float z_st3u = planes->z[detid];
	const float res_st3u = planes->spacing[detid];

	detid-= 1;
	int nhits_st3up;
	detid_list[5] = detid;
	const gHits hits_st3up = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3up);
	const float z_st3up = planes->z[detid];
	const float res_st3up = planes->spacing[detid];

	make_hitpairs_in_station_bins(hits_st3u, nhits_st3u, hits_st3up, nhits_st3up, hitpairs_v3, nhitpairs_v3, bin0_st3, nbins_st3, hitflag1, hitflag2, stid, projid);

	projid = 2;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[6] = detid;
	int nhits_st3v;
	const gHits hits_st3v = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3v);
	const float z_st3v = planes->z[detid];
	const float res_st3v = planes->spacing[detid];

	detid-= 1;
	int nhits_st3vp;
	detid_list[7] = detid;
	const gHits hits_st3vp = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3vp);
	const float z_st3vp = planes->z[detid];
	const float res_st3vp = planes->spacing[detid];

	make_hitpairs_in_station_bins(hits_st3v, nhits_st3v, hits_st3vp, nhits_st3vp, hitpairs_v3, nhitpairs_v3, bin0_st3, nbins_st3, hitflag1, hitflag2, stid, projid);
	
	
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
	float elID[4];
	float pos[4];
	float drift[4];
	float tdc[4];

	short nhits_uv;
	
	short nhits_u2, nhits_u3;
	short nhits_v2, nhits_v3;
	
	// if parallel:
	// const short nbins_st3 = geometry::N_WCHitsBins[2];
	const short nbins_total = nbins_st2*nbins_st3;
	short bin2, bin3;
	
	int nu2, nv2, nu3, nv3;
	
	int ncomb_uv;
	
	short i_u2, i_u3, i_v2, i_v3;
	
	int i_uv, i_hit;

	float y, err_y;
	
	float chi2min = 10000.1f;

	int localbin;
	
	static float tkl_data_local[datasizes::TrackletSizeMax][142];
	
	//TODO get Hodo hits stations 2, 3, 4...
	
	//get the tracks...
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackletSizeMax;
	const int N_tracklets = tklcoll->NTracks[blockIdx.x];
	
	gTracks tracks(tklcoll->TracksRawData, N_tracklets, tkl_coll_offset);
	
	//gTracklet tkl;
	float x0, tx;
	int nhits_x;
	bool update_track;
	//
	for(int i = 0; i<N_tracklets; i++){
		if(threadIdx.x!=tracks.threadID(i))continue;
		//tkl = tracks.getTracklet(i);
		x0 = tracks.x0(i);
		tx = tracks.tx(i);
		update_track = false;
		//get the corresponding local bin!
		//binId = threadIdx.x+THREADS_PER_BLOCK*localbin;
		//localbin = (tkl.stationID-tkl.threadID)/THREADS_PER_BLOCK;
		localbin = (tracks.stationID(i)-tracks.threadID(i))/THREADS_PER_BLOCK;
		
		//retrieve the bins in st2, st3
		bin2 = localbin%nbins_st2;
		bin3 = (localbin-bin2)/nbins_st2;

		nu2 = nhitpairs_u2[bin2];
		nu3 = nhitpairs_u3[bin3];
		
		if(nu2 == 0 || nu3==0) continue;

		nv2 = nhitpairs_v2[bin2];
		nv3 = nhitpairs_v3[bin3];
		
		if(nv2 == 0 || nv3==0) continue;
		
		//x0 = tkl.x0;
		//tx = tkl.tx;
		
		for(i_uv = 0; i_uv<ncomb_uv; i_uv++){
			nhits_uv = 0;
			i_u2 = i_uv%nu2;
			i_v2 = ((i_uv-i_u2)/nu2)%nv2;
			i_u3 = (((i_uv-i_u2)/nu2-i_v2)/nv2)%nu3;
			i_v3 = ((((i_uv-i_u2)/nu2)-i_v2)/nv2-i_u3)/nu3;
			
			nhits_uv = 0;
			
			if(hitpairs_u2[bin2+nbins_st2*i_u2].first>=0){
				detID[nhits_uv] = detid_list[0];
				i_hit = hitpairs_u2[bin2+nbins_st2*i_u2].first;
				Z[nhits_uv] = z_st2u;
				elID[nhits_uv] = (short)hits_st2u.chan(i_hit);
				tdc[nhits_uv] = hits_st2u.tdc(i_hit);
				drift[nhits_uv] = hits_st2u.drift(i_hit);
				pos[nhits_uv] = hits_st2u.pos(i_hit);
				calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
				Y[nhits_uv] = y;
				errY[nhits_uv] = err_y;
				nhits_uv++;
			}
			if(hitpairs_u2[bin2+nbins_st2*i_u2].second>=0){
				detID[nhits_uv] = detid_list[1];
				i_hit = hitpairs_u2[bin2+nbins_st2*i_u2].second;
				Z[nhits_uv] = z_st2up;
				elID[nhits_uv] = (short)hits_st2up.chan(i_hit);
				tdc[nhits_uv] = hits_st2up.tdc(i_hit);
				drift[nhits_uv] = hits_st2up.drift(i_hit);
				pos[nhits_uv] = hits_st2up.pos(i_hit);
				calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
				Y[nhits_uv] = y;
				errY[nhits_uv] = err_y;
				nhits_uv++;
			}
			
			nhits_u2 = nhits_uv;
			if(nhits_u2==0) continue;

			if(hitpairs_v2[bin2+nbins_st2*i_v2].first>=0){
				detID[nhits_uv] = detid_list[2];
				i_hit = hitpairs_v2[bin2+nbins_st2*i_v2].first;
				Z[nhits_uv] = z_st2v;
				elID[nhits_uv] = (short)hits_st2v.chan(i_hit);
				tdc[nhits_uv] = hits_st2v.tdc(i_hit);
				drift[nhits_uv] = hits_st2v.drift(i_hit);
				pos[nhits_uv] = hits_st2v.pos(i_hit);
				calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
				Y[nhits_uv] = y;
				errY[nhits_uv] = err_y;
				nhits_uv++;
			}
			if(hitpairs_v2[bin2+nbins_st2*i_v2].second>=0){
				detID[nhits_uv] = detid_list[3];
				i_hit = hitpairs_v2[bin2+nbins_st2*i_v2].second;
				Z[nhits_uv] = z_st2vp;
				elID[nhits_uv] = (short)hits_st2vp.chan(i_hit);
				tdc[nhits_uv] = hits_st2vp.tdc(i_hit);
				drift[nhits_uv] = hits_st2vp.drift(i_hit);
				pos[nhits_uv] = hits_st2vp.pos(i_hit);
				calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
				Y[nhits_uv] = y;
				errY[nhits_uv] = err_y;
				nhits_uv++;
			}
			
			nhits_v2 = nhits_uv-nhits_u2;
			if(nhits_v2==0) continue;

			if(hitpairs_u3[bin3+nbins_st3*i_u3].first>=0){
				detID[nhits_uv] = detid_list[0];
				i_hit = hitpairs_u3[bin3+nbins_st3*i_u3].first;
				Z[nhits_uv] = z_st3u;
				elID[nhits_uv] = (short)hits_st3u.chan(i_hit);
				tdc[nhits_uv] = hits_st3u.tdc(i_hit);
				drift[nhits_uv] = hits_st3u.drift(i_hit);
				pos[nhits_uv] = hits_st3u.pos(i_hit);
				calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
				Y[nhits_uv] = y;
				errY[nhits_uv] = err_y;
				nhits_uv++;
			}
			if(hitpairs_u3[bin3+nbins_st3*i_u3].second>=0){
				detID[nhits_uv] = detid_list[1];
				i_hit = hitpairs_u3[bin3+nbins_st3*i_u3].second;
				Z[nhits_uv] = z_st3up;
				elID[nhits_uv] = (short)hits_st3up.chan(i_hit);
				tdc[nhits_uv] = hits_st3up.tdc(i_hit);
				drift[nhits_uv] = hits_st3up.drift(i_hit);
				pos[nhits_uv] = hits_st3up.pos(i_hit);
				calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
				Y[nhits_uv] = y;
				errY[nhits_uv] = err_y;
				nhits_uv++;
			}
			
			nhits_u3 = nhits_uv-nhits_u2-nhits_v2;
			if(nhits_u3==0) continue;

			if(hitpairs_v3[bin3+nbins_st3*i_v3].first>=0){
				detID[nhits_uv] = detid_list[3];
				i_hit = hitpairs_v3[bin3+nbins_st3*i_v3].first;
				Z[nhits_uv] = z_st3v;
				elID[nhits_uv] = (short)hits_st3v.chan(i_hit);
				tdc[nhits_uv] = hits_st3v.tdc(i_hit);
				drift[nhits_uv] = hits_st3v.drift(i_hit);
				pos[nhits_uv] = hits_st3v.pos(i_hit);
				calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
				Y[nhits_uv] = y;
				errY[nhits_uv] = err_y;
				nhits_uv++;
			}
			if(hitpairs_v3[bin3+nbins_st3*i_v3].second>=0){
				detID[nhits_uv] = detid_list[3];
				i_hit = hitpairs_v3[bin3+nbins_st3*i_v3].second;
				Z[nhits_uv] = z_st3vp;
				elID[nhits_uv] = (short)hits_st3vp.chan(i_hit);
				tdc[nhits_uv] = hits_st3vp.tdc(i_hit);
				drift[nhits_uv] = hits_st3vp.drift(i_hit);
				pos[nhits_uv] = hits_st3vp.pos(i_hit);
				calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
				Y[nhits_uv] = y;
				errY[nhits_uv] = err_y;
				nhits_uv++;
			}
			
			nhits_v3 = nhits_uv-nhits_u2-nhits_v2-nhits_u3;
			if(nhits_v3==0) continue;
			
			fit_2D_track(nhits_uv, Y, Z, errY, A_, Ainv_, B_, Par, ParErr, chi2);
			
			y0 = Par[0];
			ty = Par[1];

			//TODO: hodoscope matching

			//TODO: chi2 evaluation of track candidate

			//TODO: LR ambiguity resolution
			
			if(blockIdx.x==debug::EvRef){
				printf("thread %d bin %d x0 %1.4f tx %1.4f, nhits %d \n", threadIdx.x, localbin, y0, ty, nhits_x);
			}
			
			//if(chi2>chi2min){
			//chi2min = chi2
			update_track = true;
			tkl_data_local[i][2] = nhits_uv;
			tkl_data_local[i][6] = ty;
			tkl_data_local[i][8] = y0;
			tkl_data_local[i][11] = ParErr[1];
			tkl_data_local[i][13] = ParErr[0];
			
			for(int m = 0; m<nhits_uv;m++){
				tkl_data_local[i][16+m] = detID[m];
				tkl_data_local[i][16+m+nhits_x] = elID[m];
				tkl_data_local[i][16+m+nhits_x*2] = pos[m];
				tkl_data_local[i][16+m+nhits_x*3] = tdc[m];
				tkl_data_local[i][16+m+nhits_x*4] = drift[m];
				tkl_data_local[i][105+m] = 0;
				tkl_data_local[i][124+m] = 0;
			}
			//}
			
		}
		
		if(update_track){
			nhits_x = tracks.nHits(i);
			nhits_uv = tkl_data_local[i][2];
			/*
			tklcoll->Tracklets[tkl_coll_offset+i].stationID=5;
			tklcoll->Tracklets[tkl_coll_offset+i].nHits=nhits_x+nhits_uv;
			tklcoll->Tracklets[tkl_coll_offset+i].ty=tkl_data_local[i][6];
			tklcoll->Tracklets[tkl_coll_offset+i].y0=tkl_data_local[i][8];
			tklcoll->Tracklets[tkl_coll_offset+i].err_ty=tkl_data_local[i][11];
			tklcoll->Tracklets[tkl_coll_offset+i].err_y0=tkl_data_local[i][13];
			for(int n = 0; n<nhits_uv;n++){
				tklcoll->Tracklets[tkl_coll_offset+i].hits[nhits_x+n].detectorID = tkl_data_local[i][16+n];
				tklcoll->Tracklets[tkl_coll_offset+i].hits[nhits_x+n].elementID = tkl_data_local[i][16+n+nhits_uv];
				tklcoll->Tracklets[tkl_coll_offset+i].hits[nhits_x+n].pos = tkl_data_local[i][16+n+nhits_uv*2];
				tklcoll->Tracklets[tkl_coll_offset+i].hits[nhits_x+n].tdcTime = tkl_data_local[i][16+n+nhits_uv*3];
				tklcoll->Tracklets[tkl_coll_offset+i].hits[nhits_x+n].driftDistance = tkl_data_local[i][16+n+nhits_uv*4];
				tklcoll->Tracklets[tkl_coll_offset+i].hitsign[nhits_x+n] = tkl_data_local[i][105+n];
				tklcoll->Tracklets[tkl_coll_offset+i].residual[nhits_x+n] = tkl_data_local[i][124+n];
			}
			*/
		}else{
			/*
			tklcoll->Tracklets[tkl_coll_offset+i].stationID=3;
			*/
		}
		
		
	}//end loop on existing tracklets
	
	
	
}


////////////////////////////////////
//
//          GLOBAL TRACKING
//
////////////////////////////////////

__global__ void gKernel_Global_tracking(gEventHitCollections* hitcolls, gEventTrackCollection* tklcoll, const gPlane* planes, int* eventID, bool* hastoomanyhits)
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
	detid_list[0] = detid;
	int nhits_st1u;
	const gHits hits_st1u = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1u);
	const float z_st1u = planes->z[detid];
	const float res_st1u = planes->spacing[detid];
	
	detid-= 1;
	detid_list[1] = detid;
	int nhits_st1up;
	const gHits hits_st1up = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1up);
	const float z_st1up = planes->z[detid];
	const float res_st1up = planes->spacing[detid];
	
	projid = 2;
	detid = geometry::detsuperid[stid][projid]*2;
	detid_list[2] = detid;
	int nhits_st1v;
	const gHits hits_st1v = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1v);
	const float z_st1v = planes->z[detid];
	const float res_st1v = planes->spacing[detid];
	
	detid-= 1;
	detid_list[3] = detid;
	int nhits_st1vp;
	const gHits hits_st1vp = hitcolls->hitschambers(blockIdx.x, geometry::eff_detid_chambers[detid-1], nhits_st1vp);
	const float z_st1vp = planes->z[detid];
	const float res_st1vp = planes->spacing[detid];

	//get the tracks...
	const unsigned int tkl_coll_offset = blockIdx.x*datasizes::TrackletSizeMax;
	const int N_tracklets = tklcoll->NTracks[blockIdx.x];
	
	gTracks tracks(tklcoll->TracksRawData, N_tracklets, tkl_coll_offset);
	
	//gTracklet tkl;
	
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
	short nhits, nhits_new;
	
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
	float tdc[6];

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

	float tkl_data_local[datasizes::TrackletSizeMax][142];
		
	for(int i = 0; i<N_tracklets; i++){
		//tkl = tracks.getTracklet(i);
		update_track = false;
		/*
		if(threadIdx.x!=tkl.threadID)continue;
		x0 = tkl.x0;
		tx = tkl.tx;
		errx0 = tkl.err_x0;
		errtx = tkl.err_tx;

		y0 = tkl.err_y0;
		ty = tkl.err_ty;
		erry0 = tkl.err_y0;
		errty = tkl.err_ty;
		*/
		
		if(threadIdx.x!=tracks.threadID(i))continue;
		
		x0 = tracks.x0(i);
		tx = tracks.tx(i);
		errx0 = tracks.err_x0(i);
		errtx = tracks.err_tx(i);

		y0 = tracks.y0(i);
		ty = tracks.ty(i);
		erry0 = tracks.err_y0(i);
		errty = tracks.err_ty(i);

		SagittaRatioInStation1(x0, tx, y0, ty, tracks.hits_detid(i, tracks.nHits(i)-1), pos_exp, window, planes->z, planes->costheta, planes->sintheta);
		
		projid = 0;
		nx1 = make_hitpairs_in_station(hits_st1x, nhits_st1x, hits_st1xp, nhits_st1xp, hitpairs_x1, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);
		
		projid = 1;
		nu1 = make_hitpairs_in_station(hits_st1u, nhits_st1u, hits_st1up, nhits_st1up, hitpairs_u1, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);
		
		projid = 2;
		nv1 = make_hitpairs_in_station(hits_st1v, nhits_st1v, hits_st1vp, nhits_st1vp, hitpairs_v1, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);
		
		if(nx1==0 || nu1==0 || nv1==0)continue;
		
		ncomb_uv = nu1*nv1;
		
		//triple loop on hits
		for(i_x = 0; i_x<nx1; i_x++){
			nhits_x = 0;
			if(hitpairs_x1[i_x].first>=0){
				detID[nhits_x] = detid_list[0];
				i_hit = hitpairs_x1[i_x].first;
				X[nhits_x] = hits_st1x.pos(i_hit);
				errX[nhits_x] = res_st1x;
				Z[nhits_x] = z_st1x;
				elID[nhits_x] = (short)hits_st1x.chan(i_hit);
				tdc[nhits_x] = hits_st1x.tdc(i_hit);
				drift[nhits_x] = hits_st1x.drift(i_hit);
				nhits_x++;
			}
			if(hitpairs_x1[i_x].second>=0){
				detID[nhits_x] = detid_list[1];
				i_hit = hitpairs_x1[i_x].second;
				X[nhits_x] = hits_st1xp.pos(i_hit);
				errX[nhits_x] = res_st1xp;
				Z[nhits_x] = z_st1xp;
				elID[nhits_x] = (short)hits_st1xp.chan(i_hit);
				tdc[nhits_x] = hits_st1xp.tdc(i_hit);
				drift[nhits_x] = hits_st1xp.drift(i_hit);
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
			
			//add the UV hits
			for(i_uv = 0; i_uv<ncomb_uv; i_uv++){
				i_u = i_uv%nu1;
				i_v = (i_uv-i_u)/nu1;
				
				nhits_uv = 0;
				if(hitpairs_u1[i_u].first>=0){
					detID[nhits_x+nhits_uv] = detid_list[0];
					i_hit = hitpairs_u1[i_u].first;
					elID[nhits_x+nhits_uv] = (short)hits_st1u.chan(i_hit);
					tdc[nhits_x+nhits_uv] = hits_st1u.tdc(i_hit);
					drift[nhits_x+nhits_uv] = hits_st1u.drift(i_hit);
					pos[nhits_x+nhits_uv] = hits_st1u.pos(i_hit);
					calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
					if( (y-y_trk(y0, ty, z_st1u))*(y-y_trk(y0, ty, z_st1u))/(err_y*err_y+err_y_trk(erry0, errty, z_st1u)*err_y_trk(erry0, errty, z_st1u))<2.0 )
						nhits_uv++;
				}
				if(hitpairs_u1[i_u].second>=0){
					detID[nhits_uv] = detid_list[1];
					i_hit = hitpairs_u1[i_u].second;
					Z[nhits_uv] = z_st1up;
					elID[nhits_uv] = (short)hits_st1up.chan(i_hit);
					tdc[nhits_uv] = hits_st1up.tdc(i_hit);
					drift[nhits_uv] = hits_st1up.drift(i_hit);
					pos[nhits_uv] = hits_st1up.pos(i_hit);
					calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
					if( (y-y_trk(y0, ty, z_st1up))*(y-y_trk(y0, ty, z_st1up))/(err_y*err_y+err_y_trk(erry0, errty, z_st1up)*err_y_trk(erry0, errty, z_st1up))<2.0 )
						nhits_uv++;
				}
				nhits_u = nhits_uv;
				if(nhits_u==0)continue;
				
				if(hitpairs_v1[i_v].first>=0){
					detID[nhits_x+nhits_uv] = detid_list[0];
					i_hit = hitpairs_v1[i_v].first;
					elID[nhits_x+nhits_uv] = (short)hits_st1v.chan(i_hit);
					tdc[nhits_x+nhits_uv] = hits_st1v.tdc(i_hit);
					drift[nhits_x+nhits_uv] = hits_st1v.drift(i_hit);
					pos[nhits_x+nhits_uv] = hits_st1v.pos(i_hit);
					calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
					if( (y-y_trk(y0, ty, z_st1v))*(y-y_trk(y0, ty, z_st1v))/(err_y*err_y+err_y_trk(erry0, errty, z_st1v)*err_y_trk(erry0, errty, z_st1v))<2.0 )
						nhits_uv++;
				}
				if(hitpairs_v1[i_v].second>=0){
					detID[nhits_uv] = detid_list[1];
					i_hit = hitpairs_v1[i_v].second;
					Z[nhits_uv] = z_st1vp;
					elID[nhits_uv] = (short)hits_st1vp.chan(i_hit);
					tdc[nhits_uv] = hits_st1vp.tdc(i_hit);
					drift[nhits_uv] = hits_st1vp.drift(i_hit);
					pos[nhits_uv] = hits_st1vp.pos(i_hit);
					calculate_y_uvhit(detID[nhits_uv], elID[nhits_uv], drift[nhits_uv], 0, x0, tx, planes, y, err_y);
					if( (y-y_trk(y0, ty, z_st1vp))*(y-y_trk(y0, ty, z_st1vp))/(err_y*err_y+err_y_trk(erry0, errty, z_st1vp)*err_y_trk(erry0, errty, z_st1vp))<2.0 )
						nhits_uv++;
				}
				nhits_v = nhits_uv-nhits_u;
				if(nhits_v==0)continue;
				
				//TODO: matching hodoscope, 
				
				
				//TODO: resolve left right...
				
				
				//if(chi2>chi2min){
				//chi2min = chi2
				update_track = true;
				tkl_data_local[i][2] = nhits_x+nhits_uv;
				tkl_data_local[i][9] = invP;
				tkl_data_local[i][14] = errinvP;
				tkl_data_local[i][15] = charge;
				
				for(int m = 0; m<nhits_x+nhits_uv;m++){
				tkl_data_local[i][16+m] = detID[m];
				tkl_data_local[i][16+m+nhits_x] = elID[m];
				tkl_data_local[i][16+m+nhits_x*2] = pos[m];
				tkl_data_local[i][16+m+nhits_x*3] = tdc[m];
				tkl_data_local[i][16+m+nhits_x*4] = drift[m];
				tkl_data_local[i][105+m] = 0;
				tkl_data_local[i][124+m] = 0;
				}
				//}

				
			}
			
		}//end
		
		if(update_track){
			nhits = tracks.nHits(i);
			nhits_new = tkl_data_local[i][2];
			/*
			tklcoll->Tracklets[tkl_coll_offset+i].stationID=6;
			tklcoll->Tracklets[tkl_coll_offset+i].nHits=nhits+nhits_new;
			tklcoll->Tracklets[tkl_coll_offset+i].invP=tkl_data_local[i][9];
			tklcoll->Tracklets[tkl_coll_offset+i].err_invP=tkl_data_local[i][14];
			tklcoll->Tracklets[tkl_coll_offset+i].err_y0=tkl_data_local[i][15];
			for(int n = 0; n<nhits_new;n++){
				tklcoll->Tracklets[tkl_coll_offset+i].hits[nhits+n].detectorID = tkl_data_local[i][16+n];
				tklcoll->Tracklets[tkl_coll_offset+i].hits[nhits+n].elementID = tkl_data_local[i][16+n+nhits_new];
				tklcoll->Tracklets[tkl_coll_offset+i].hits[nhits+n].pos = tkl_data_local[i][16+n+nhits_new*2];
				tklcoll->Tracklets[tkl_coll_offset+i].hits[nhits+n].tdcTime = tkl_data_local[i][16+n+nhits_new*3];
				tklcoll->Tracklets[tkl_coll_offset+i].hits[nhits+n].driftDistance = tkl_data_local[i][16+n+nhits_new*4];
				tklcoll->Tracklets[tkl_coll_offset+i].hitsign[nhits+n] = tkl_data_local[i][105+n];
				tklcoll->Tracklets[tkl_coll_offset+i].residual[nhits+n] = tkl_data_local[i][124+n];
			}
			*/
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
			printf("thread %d tracklet %d thread %1.0f bin/stid %1.0f x0 %1.4f tx %1.4f y0 %1.4f ty %1.4f nhits %1.0f: \n", threadIdx.x, i, Tracks.threadID(i), Tracks.stationID(i), Tracks.x0(i), Tracks.tx(i), Tracks.y0(i), Tracks.ty(i), Tracks.nHits(i));
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


__global__ void gKernel_XZ_YZ_tracking_new(gEvent* ic, gOutputEvent* oc, gStraightTrackBuilder* straighttrackbuilder, gStraightFitArrays* fitarrays, const gPlane* planes, bool best_candyz_only = true)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	//printf("index %d, thread %d, block %d \n", index, threadIdx.x, blockIdx.x);

	oc->EventID[index] = ic->EventID[index];
	oc->nAH[index] = ic->nAH[index];
	oc->HasTooManyHits[index] = hastoomanyhits[index];
	oc->nTracklets[index] = 0;
	
	if(oc->HasTooManyHits[index])return;
	
	const int idxoff_global = index*EstnAHMax;
	const int idxoff_tkl = index*TrackletSizeMax;
	
	const short nbins_st2 = geometry::N_WCHitsBins[1];
	const short nbins_st3 = geometry::N_WCHitsBins[2]+geometry::N_WCHitsBins[3];
	
	const int st2x_offset = nbins_st2*geometry::MaxHitsProj[0]*index;
	const int st2uv_offset = nbins_st2*geometry::MaxHitsProj[1]*index;
	
	const int st3x_offset = nbins_st3*geometry::MaxHitsProj[0]*index;
	const int st3uv_offset = nbins_st3*geometry::MaxHitsProj[1]*index;
	
	short stid, projid;
	
	unsigned int hitidx1[100];
	unsigned int hitidx2[100];
	short hitflag1[100];
	short hitflag2[100];

	int nhitpairs_x2[28];
	int nhitpairs_u2[28];
	int nhitpairs_v2[28];
	
	int nhitpairs_x3[58];
	int nhitpairs_u3[58];
	int nhitpairs_v3[58];
	
	float A_[4];
	float Ainv_[4];
	float B_[2];
	float Par[4];
	float ParErr[4];
	float chi2;
	
        // 1- get hit pairs in bins in st2, st3:
	// the pairs for st3p and st3m are now in the same array. to allow for parallelization of d3p/d3m if needed;
	projid = 0;
	//D2: stid = 3-1
	stid = 2;
	//make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_x2, straighttrackbuilder->nhitpairs_x2, straighttrackbuilder->hitidx1, straighttrackbuilder->hitidx2, straighttrackbuilder->hitflag1, straighttrackbuilder->hitflag2, stid, projid);
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_x2, nhitpairs_x2, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);

	//D3p: stid = 4-1
	stid = 3;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_x3, nhitpairs_x3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);

	//D3p: stid = 5-1
	stid = 4;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_x3, nhitpairs_x3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
		
#ifdef DEBUG
	for(int bin = 0; bin<58; bin++){
	if(nhitpairs_x3[bin])printf("evt %d bin %d nhits x3 = %d\n", ic[index].EventID, bin, nhitpairs_x3[bin]);
	for(int i = 0; i<nhitpairs_x3[bin]; i++){
		if(hitpairs_x3[bin+28*i].first>-1)printf(" evt %d bin %d  i %d detid1 %d elid1 %d\n", ic[index].EventID, bin, i, 
				ic->AllHits[idxoff_global+hitpairs_x3[bin+29*i].first].detectorID, ic->AllHits[idxoff_global+hitpairs_x3[bin+29*i].first].elementID);
		if(hitpairs_x3[bin+28*i].second>-1)printf(" evt %d bin %d i %d detid2 %d elid2 %d\n", ic[index].EventID, bin, i, 
				ic->AllHits[idxoff_global+hitpairs_x3[bin+29*i].second].detectorID, ic->AllHits[idxoff_global+hitpairs_x3[bin+29*i].second].elementID);
	}
	}
#endif
	
	projid = 1;
	stid = 2;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_u2, nhitpairs_u2, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
	
	stid = 3;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_u3, nhitpairs_u3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
	
	stid = 4;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_u3, nhitpairs_u3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
		
	projid = 2;
	stid = 2;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_v2, nhitpairs_v2, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
		
	stid = 3;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_v3, nhitpairs_v3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
		
	stid = 4;
	make_hitpairs_in_station_bins(ic, straighttrackbuilder->hitpairs_v3, nhitpairs_v3, hitidx1, hitidx2, hitflag1, hitflag2, stid, projid);
	
	bool bin_overflows = false;
	for(int bin = 0; bin<28; bin++){
		if(nhitpairs_x2[bin]>geometry::MaxHitsProj[0])bin_overflows = true;
		//if(nhitpairs_x2[bin])printf("evt %d bin %d nhits x2 = %d\n", ic[index].EventID, bin, nhitpairs_x2[bin]);
		if(nhitpairs_u2[bin]>geometry::MaxHitsProj[1])bin_overflows = true;
		//if(nhitpairs_u2[bin])printf("evt %d bin %d nhits u2 = %d\n", ic[index].EventID, bin, nhitpairs_u2[bin]);
		if(nhitpairs_v2[bin]>geometry::MaxHitsProj[2])bin_overflows = true;
		//if(nhitpairs_v2[bin])printf("evt %d bin %d nhits v2 = %d\n", ic[index].EventID, bin, nhitpairs_v2[bin]);
	}
	for(int bin = 0; bin<58; bin++){
		if(nhitpairs_x3[bin]>geometry::MaxHitsProj[0])bin_overflows = true;
		//if(nhitpairs_x3[bin])printf("evt %d bin %d nhits x3 = %d\n", ic[index].EventID, bin, nhitpairs_x3[bin]);
		if(nhitpairs_u3[bin]>geometry::MaxHitsProj[1])bin_overflows = true;
		//if(nhitpairs_u3[bin])printf("evt %d bin %d nhits u3 = %d\n", ic[index].EventID, bin, nhitpairs_u3[bin]);
		if(nhitpairs_v3[bin]>geometry::MaxHitsProj[2])bin_overflows = true;
		//if(nhitpairs_v3[bin])printf("evt %d bin %d nhits v3 = %d\n", ic[index].EventID, bin, nhitpairs_v3[bin]);
	}
	if(bin_overflows){
		printf("Evt %d discarded, too many hits\n", ic[index].EventID);
		return;
	}
	
	// declare variables useful for the loop 
	int ntkl = 0;
	
	//for the time being, declare in function;
	short nprop, iprop, idet;
	float xExp, ipos;
	
	short nhits_x;
	short nhits_uv;
	short nhits_v;
	
	short nhits_x2, nhits_x3;
	short nhits_u2, nhits_u3;
	short nhits_v2, nhits_v3;
	
	// if parallel:
	// const short nbins_st3 = geometry::N_WCHitsBins[2];
	const short nbins_total = nbins_st2*nbins_st3;
	short bin2, bin3;
	
	int nx2, nu2, nv2;
	int nx3, nu3, nv3;
	
	int ncomb_x, ncomb_uv;
	
	short i_x2, i_x3;
	short i_u2, i_u3, i_v2, i_v3;
	
	int i_x, i_uv, i_hit;

	float y, err_y;
	float ty;
	
	float chi2min = 10000.1f;
	
	int n_goodxz;
	//int hitidx;
	//gHit hit;
	gTrack2D trackXZ, trackYZ, besttrackYZ;
	
	//loop on bins FIRST
	for(short i = 0; i<nbins_total; i++){
		bin2 = i%nbins_st2;
		bin3 = (i-bin2)/nbins_st2;
		//printf("bin %d %d\n", bin2, bin3);
		// if parallel:
		// bin3+=nbins_st3*i_thread;

		nx2 = nhitpairs_x2[bin2];
		nx3 = nhitpairs_x3[bin3];
		
		//printf("nx %d %d \n", nx2, nx3);
		if(nx2 == 0 || nx3==0) continue;
		
		nu2 = nhitpairs_u2[bin2];
		nu3 = nhitpairs_u3[bin3];

		//printf("nu %d %d \n", nu2, nu3);
		if(nu2 == 0 || nu3==0) continue;
		
		nv2 = nhitpairs_v2[bin2];
		nv3 = nhitpairs_v3[bin3];

		//printf("nv %d %d \n", nv2, nv3);
		if(nv2 == 0 || nv3==0) continue;
		
		// evaluating the number of combinations;
		ncomb_x = nx2*nx3;
		ncomb_uv = nu2*nu3*nv2*nv3;
		
		n_goodxz = 0;
				
		for(i_x = 0; i_x<ncomb_x; i_x++){
			i_x2 = i_x%nx3;
			i_x3 = (i_x-i_x2)/nx3;
			
			nhits_v = -1;
			nhits_x = 0;

			if(straighttrackbuilder->hitpairs_x2[st2x_offset+bin2+nbins_st2*i_x2].first>=0){
				trackXZ.hitlist[nhits_x] = straighttrackbuilder->hitpairs_x2[st2x_offset+bin2+nbins_st2*i_x2].first;
				FillFitArrays_X(nhits_x, ic->AllHits[trackXZ.hitlist[nhits_x]], 0, fitarrays, planes);
				nhits_x++;
			}
			if(straighttrackbuilder->hitpairs_x2[st2x_offset+bin2+nbins_st2*i_x2].second>=0){
				trackXZ.hitlist[nhits_x] = straighttrackbuilder->hitpairs_x2[st2x_offset+bin2+nbins_st2*i_x2].second;
				FillFitArrays_X(nhits_x, ic->AllHits[trackXZ.hitlist[nhits_x]], 0, fitarrays, planes);
				nhits_x++;
			}
			
			nhits_x2 = nhits_x;
			if(nhits_x2==0) continue;
			
			if(straighttrackbuilder->hitpairs_x3[st3x_offset+bin3+nbins_st3*i_x3].first>=0){
				trackXZ.hitlist[nhits_x] = straighttrackbuilder->hitpairs_x3[st3x_offset+bin3+nbins_st3*i_x3].first;
				FillFitArrays_X(nhits_x, ic->AllHits[trackXZ.hitlist[nhits_x]], 0, fitarrays, planes);
				nhits_x++;
			}
			if(straighttrackbuilder->hitpairs_x3[st3x_offset+bin3+nbins_st3*i_x3].second>=0){
				trackXZ.hitlist[nhits_x] = straighttrackbuilder->hitpairs_x3[st3x_offset+bin3+nbins_st3*i_x3].second;
				FillFitArrays_X(nhits_x, ic->AllHits[trackXZ.hitlist[nhits_x]], 0, fitarrays, planes);
				nhits_x++;
			}
			
			nhits_x3 = nhits_x-nhits_x2;
			if(nhits_x3==0) continue;
			
			fit_2D_track(nhits_x, fitarrays->x_array, fitarrays->z_array, fitarrays->dx_array, A_, Ainv_, B_, Par, ParErr, chi2);
			if(fabs(Par[0])>1.05*X0_MAX || fabs(Par[1])>1.05*TX_MAX)continue;
			
			trackXZ.x_0 = Par[0];
			trackXZ.err_x_0 = ParErr[0];
			trackXZ.tx_ = Par[1];
			trackXZ.err_tx_ = ParErr[1];
			trackXZ.nhits = nhits_x;
			
			n_goodxz++;
			
			//prop matching
			nprop = 0;
			for(int n = 0; n<ic->nAH[index]; n++){
				idet = ic->AllHits[idxoff_global+n].detectorID;
				ipos = ic->AllHits[idxoff_global+n].pos;
				for(short ip = 0; ip<4; ip++){
					iprop = 48+ip;
					xExp = trackXZ.tx_*planes->z[iprop]+trackXZ.x_0;
					//loop on hits to find prop
					if(idet==iprop){
						if(fabs(ipos-xExp)<5.08f){
							nprop++;
							break;
						}
					}
				}
				if(nprop>0)break;
			}
			if(nprop==0)continue;

			//partially fill the tracklet xz info:
			oc->AllTracklets[idxoff_tkl+ntkl].stationID = 5;
			oc->AllTracklets[idxoff_tkl+ntkl].x0 = trackXZ.x_0;
			oc->AllTracklets[idxoff_tkl+ntkl].err_x0 = trackXZ.err_x_0;
			oc->AllTracklets[idxoff_tkl+ntkl].tx = trackXZ.tx_;
			oc->AllTracklets[idxoff_tkl+ntkl].err_tx = trackXZ.err_tx_;

			for(i_hit = 0; i_hit<nhits_x; i_hit++){
				FillChi2Arrays(i_hit, ic->AllHits[idxoff_global+trackXZ.hitlist[i_hit]], 0, fitarrays, planes);
				oc->AllTracklets[idxoff_tkl+ntkl].hits[i_hit] = ic->AllHits[idxoff_global+trackXZ.hitlist[i_hit]];
			}
			oc->AllTracklets[idxoff_tkl+ntkl].nXHits = nhits_x;
			
			for(i_uv = 0; i_uv<ncomb_uv; i_uv++){
				nhits_uv = 0;
				i_u2 = i_uv%nu2;
				i_v2 = ((i_uv-i_u2)/nu2)%nv2;
				i_u3 = (((i_uv-i_u2)/nu2-i_v2)/nv2)%nu3;
				i_v3 = ((((i_uv-i_u2)/nu2)-i_v2)/nv2-i_u3)/nu3;
				
				if(straighttrackbuilder->hitpairs_u2[st2uv_offset+bin2+nbins_st2*i_u2].first>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_u2[st2uv_offset+bin2+nbins_st2*i_u2].first], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_u2[st2uv_offset+bin2+nbins_st2*i_u2].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder->hitpairs_u2[st2uv_offset+bin2+nbins_st2*i_u2].second>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_u2[st2uv_offset+bin2+nbins_st2*i_u2].second], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_u2[st2uv_offset+bin2+nbins_st2*i_u2].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}
				
				nhits_u2 = nhits_uv;
				if(nhits_u2==0) continue;
				
				if(straighttrackbuilder->hitpairs_v2[st2uv_offset+bin2+nbins_st2*i_v2].first>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_v2[st2uv_offset+bin2+nbins_st2*i_v2].first], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_v2[st2uv_offset+bin2+nbins_st2*i_v2].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder->hitpairs_v2[st2uv_offset+bin2+nbins_st2*i_v2].second>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_v2[st2uv_offset+bin2+nbins_st2*i_v2].second], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_v2[st2uv_offset+bin2+nbins_st2*i_v2].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}

				nhits_v2 = nhits_uv-nhits_u2;
				if(nhits_v2==0) continue;
				
				if(straighttrackbuilder->hitpairs_u3[st3uv_offset+bin3+nbins_st3*i_u3].first>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_u3[st3uv_offset+bin3+nbins_st3*i_u3].first], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_u3[st3uv_offset+bin3+nbins_st3*i_u3].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder->hitpairs_u3[st3uv_offset+bin3+nbins_st3*i_u3].second>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_u3[st3uv_offset+bin3+nbins_st3*i_u3].second], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_u3[st3uv_offset+bin3+nbins_st3*i_u3].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}

				nhits_u3 = nhits_uv-nhits_u2-nhits_v2;
				if(nhits_u3==0) continue;

				if(straighttrackbuilder->hitpairs_v3[st3uv_offset+bin3+nbins_st3*i_v3].first>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_v3[st3uv_offset+bin3+nbins_st3*i_v3].first], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_v3[st3uv_offset+bin3+nbins_st3*i_v3].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder->hitpairs_v3[st3uv_offset+bin3+nbins_st3*i_v3].second>=0){
					if(calculate_y_uvhit(y, err_y, ic->AllHits[straighttrackbuilder->hitpairs_v3[st3uv_offset+bin3+nbins_st3*i_v3].second], 0, trackXZ, planes)){
					trackYZ.hitlist[nhits_uv] = straighttrackbuilder->hitpairs_v3[st3uv_offset+bin3+nbins_st3*i_v3].first;
					FillFitArrays_UV(nhits_uv, ic->AllHits[trackYZ.hitlist[nhits_uv]], fitarrays, planes, y, err_y);
					nhits_uv++;
					}
				}

				nhits_v3 = nhits_uv-nhits_u2-nhits_v2-nhits_u3;
				if(nhits_v3==0) continue;

				fit_2D_track(nhits_uv, fitarrays->y_array, fitarrays->z_array, fitarrays->dy_array, A_, Ainv_, B_, Par, ParErr, chi2);
				
				trackYZ.x_0 = Par[0];
				trackYZ.err_x_0 = ParErr[0];
				trackYZ.tx_ = Par[1];
				trackYZ.err_tx_ = ParErr[1];
				trackYZ.nhits = nhits_uv;
				
				// now evaluate the back track candidate.
				oc->AllTracklets[idxoff_tkl+ntkl].y0 = trackYZ.x_0;
				oc->AllTracklets[idxoff_tkl+ntkl].err_y0 = trackYZ.err_x_0;
				oc->AllTracklets[idxoff_tkl+ntkl].ty = trackYZ.tx_;
				oc->AllTracklets[idxoff_tkl+ntkl].err_ty = trackYZ.err_tx_;

				
				// filter with hodoscope matching before evaluating chi2.
				if(!match_tracklet_to_hodo(oc[index].AllTracklets[ntkl], 2, ic, planes))continue;
				if(!match_tracklet_to_hodo(oc[index].AllTracklets[ntkl], 3, ic, planes) &&
				 !match_tracklet_to_hodo(oc[index].AllTracklets[ntkl], 4, ic, planes))continue;
				
				for(i_hit = nhits_x; i_hit<nhits_x+nhits_uv; i_hit++){
					oc->AllTracklets[idxoff_tkl+ntkl].hits[i_hit] = ic->AllHits[idxoff_global+trackYZ.hitlist[i_hit-nhits_x]];
				}
				oc->AllTracklets[idxoff_tkl+ntkl].nUHits = nhits_u2+nhits_u3;
				oc->AllTracklets[idxoff_tkl+ntkl].nVHits = nhits_v2+nhits_v3;
				
				resolve_leftright(oc[index].AllTracklets[ntkl], planes, 40.);
				resolve_leftright(oc[index].AllTracklets[ntkl], planes, 150.);
				resolve_single_leftright(oc[index].AllTracklets[ntkl], planes);
			
				chi2 = refit_backpartialtrack_with_drift(oc[index].AllTracklets[ntkl], fitarrays, planes);

				resolve_leftright(oc[index].AllTracklets[ntkl], planes, 40.);
				resolve_leftright(oc[index].AllTracklets[ntkl], planes, 150.);
				resolve_single_leftright(oc[index].AllTracklets[ntkl], planes);
			
				for(i_hit = 0; i_hit<nhits_x; i_hit++){
					trackXZ.hitsign[i_hit] = oc->AllTracklets[idxoff_tkl+ntkl].hitsign[i_hit];
				}
				for(i_hit = nhits_x; i_hit<nhits_x+nhits_uv; i_hit++){
					trackYZ.hitsign[i_hit-nhits_x] = oc->AllTracklets[idxoff_tkl+ntkl].hitsign[i_hit];
				}
				
				if(best_candyz_only){
					if(chi2<chi2min){
						chi2min = chi2;
						besttrackYZ.x_0 = oc->AllTracklets[idxoff_tkl+ntkl].y0;
						besttrackYZ.err_x_0 = oc->AllTracklets[idxoff_tkl+ntkl].err_y0;
						besttrackYZ.tx_ = oc->AllTracklets[idxoff_tkl+ntkl].ty;
						besttrackYZ.err_tx_ = oc->AllTracklets[idxoff_tkl+ntkl].err_ty;
						besttrackYZ.nhits = nhits_uv;
						nhits_v = nhits_v2+nhits_v3;
						for(i_hit = 0; i_hit<nhits_uv; i_hit++){
							besttrackYZ.hitlist[i_hit] = trackYZ.hitlist[i_hit];
							besttrackYZ.hitsign[i_hit] = trackYZ.hitsign[i_hit];
						}
					}
				}else{
					if(ntkl<TrackletSizeMax){
						ntkl++;
						oc->AllTracklets[idxoff_tkl+ntkl].chisq = chi2;
						oc->AllTracklets[idxoff_tkl+ntkl].stationID = 5;
						oc->AllTracklets[idxoff_tkl+ntkl].x0 = trackXZ.x_0;
						oc->AllTracklets[idxoff_tkl+ntkl].err_x0 = trackXZ.err_x_0;
						oc->AllTracklets[idxoff_tkl+ntkl].tx = trackXZ.tx_;
						oc->AllTracklets[idxoff_tkl+ntkl].err_tx = trackXZ.err_tx_;
						oc->AllTracklets[idxoff_tkl+ntkl].y0 = besttrackYZ.x_0;
						oc->AllTracklets[idxoff_tkl+ntkl].err_y0= besttrackYZ.err_x_0;
						oc->AllTracklets[idxoff_tkl+ntkl].ty = besttrackYZ.tx_;
						oc->AllTracklets[idxoff_tkl+ntkl].err_ty = besttrackYZ.err_tx_;
					}
				}
			}// end loop on uv hits

			if(best_candyz_only){
				if(nhits_v<0 || chi2min>=10000.f)continue;
			
				oc->AllTracklets[idxoff_tkl+ntkl].chisq = chi2min;
				oc->AllTracklets[idxoff_tkl+ntkl].y0 = besttrackYZ.x_0;
				oc->AllTracklets[idxoff_tkl+ntkl].err_y0= besttrackYZ.err_x_0;
				oc->AllTracklets[idxoff_tkl+ntkl].ty = besttrackYZ.tx_;
				oc->AllTracklets[idxoff_tkl+ntkl].err_ty = besttrackYZ.err_tx_;
			
				nhits_uv = trackYZ.nhits;
			
				for(i_hit = 0; i_hit<nhits_uv; i_hit++){
					oc->AllTracklets[idxoff_tkl+ntkl].hits[i_hit+nhits_x] = ic->AllHits[idxoff_global+trackYZ.hitlist[i_hit]];
				}
				if(ntkl<TrackletSizeMax)ntkl++;	
			}
		}// end loop on x hits
		//if(n_goodxz==0)printf("bin2 %d bin3 %d\n", bin2, bin3);
	}//end loop on bins
		
	oc->nTracklets[index] = ntkl;
}


// ------------------------------
// Global Track candidates
// ------------------------------


__global__ void gKernel_GlobalTrack_building(gEvent* ic, gOutputEvent* oc, gFullTrackBuilder* fulltrackbuilder, gStraightFitArrays* fitarrays, const gPlane* planes, bool ktracker)
{
	int index = threadIdx.x + blockIdx.x * blockDim.x;

	//if(oc[index].HasTooManyHits)return;
	
	int N_tkl = oc[index].nTracklets;
	
	int nx1, nu1, nv1;
	float pos_exp[3];
	float window[3];
	short projid;
	short stid = 0;
	short nhits_X1, nhits_U1, nhits_V1, nhits_1;
	//float y, err_y;
	
	float umin1, umax1;
	float vmin1, vmax1;
	float xpos1, upos1, vpos1;
	
	float z_x1, z_u1, z_v1;
	float v_win1, v_win2, v_win3, v_win;
	
	if(ktracker){
	
	for(int i = 0; i<oc[index].nTracklets; i++){
		SagittaRatioInStation1(oc[index].AllTracklets[i], pos_exp, window, planes);
		
		projid = 0;
		nx1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_x1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);
		projid = 1;
		nu1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_u1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);
		projid = 2;
		nv1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_v1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, pos_exp[projid]-window[projid], pos_exp[projid]+window[projid]);
		
		if(nx1==0 || nu1==0 || nv1==0)continue;
		
		//triple loop on hits
		for(int j = 0; j<nx1; j++){
			nhits_X1 = 0;
			//in x, we can already make a first evaluation of tx_st1, and a first evaluation of Pinv
			xpos1 = fulltrackbuilder[index].hitpairs_x1[j].second>=0 ? 0.5f*(ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].pos+ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].pos) : ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].pos;
			
			umin1 = xpos1*planes[geometry::detsuperid[0][1]*2].costheta-planes[geometry::detsuperid[0][1]*2].u_win;
			umax1 = umin1+2*planes[geometry::detsuperid[0][1]*2].u_win;
			
			if(fulltrackbuilder[index].hitpairs_x1[j].first>=0){
				fitarrays->x_array[nhits_X1] = ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].pos;
				fitarrays->dx_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].detectorID].spacing/3.4641f;
				fitarrays->z_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].detectorID].z;
				fulltrackbuilder[index].hitlist[nhits_X1] = fulltrackbuilder[index].hitpairs_x1[j].first;
				nhits_X1++;
			}
			if(fulltrackbuilder[index].hitpairs_x1[j].second>=0){
				fitarrays->x_array[nhits_X1] = ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].pos;
				fitarrays->dx_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].detectorID].spacing/3.4641f;
				fitarrays->z_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].detectorID].z;
				fulltrackbuilder[index].hitlist[nhits_X1] = fulltrackbuilder[index].hitpairs_x1[j].second;
				nhits_X1++;
			}
			
			//fit to obtain the st1 track segment: this time fitting from KMag, not from the origin:
			fitarrays->x_array[nhits_X1] = oc[index].AllTracklets[i].x0+oc[index].AllTracklets[i].tx*geometry::Z_KMAG_BEND;
			fitarrays->dx_array[nhits_X1] = oc[index].AllTracklets[i].err_x0+oc[index].AllTracklets[i].err_tx*geometry::Z_KMAG_BEND;
			fitarrays->z_array[nhits_X1] = geometry::Z_KMAG_BEND;
			
			fit_2D_track(nhits_X1+1, fitarrays->x_array, fitarrays->z_array, fitarrays->dx_array, fitarrays->A, fitarrays->Ainv, fitarrays->B, fitarrays->output_parameters, fitarrays->output_parameters_errors, fitarrays->chi2_2d);
			
			fulltrackbuilder[index].TrackXZ_st1.x_0 = fitarrays->output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.err_x_0 = fitarrays->output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.tx_ = fitarrays->output_parameters[1];
			fulltrackbuilder[index].TrackXZ_st1.err_tx_ = fitarrays->output_parameters[1];
			
			for(int k = 0; k<nu1; k++){
				nhits_U1 = nhits_X1;
				upos1 = fulltrackbuilder[index].hitpairs_u1[k].second>=0 ? 0.5f*(ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_u1[k].first].pos+ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_u1[k].second].pos) : ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_u1[k].first].pos;
				
				if(upos1<umin1 || upos1>umax1)continue;
				
				z_x1 = fulltrackbuilder[index].hitpairs_x1[j].second>=0 ? planes[geometry::detsuperid[0][0]*2].z_mean : planes[geometry::detsuperid[0][0]*2].z; 
				z_u1 = fulltrackbuilder[index].hitpairs_u1[k].second>=0 ? planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_u1[k].first].detectorID].z_mean : planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_u1[k].first].detectorID].z; 
				z_v1 = planes[geometry::detsuperid[0][2]*2].z_mean;
				
				v_win1 = planes[ ic->AllHits[idxoff_global+ fulltrackbuilder[index].hitpairs_u1[k].first ].detectorID].v_win_fac1;
				v_win2 = fabs(z_u1+z_v1-2*z_x1)*planes[ geometry::detsuperid[0][2]*2 ].v_win_fac2;
				v_win3 = fabs((z_v1-z_u1)*planes[ geometry::detsuperid[0][2]*2 ].v_win_fac3);
				v_win = v_win1+v_win2+v_win3+2*planes[ ic->AllHits[idxoff_global+ fulltrackbuilder[index].hitpairs_u1[k].first ].detectorID].spacing;
			
				vmin1 = 2*xpos1*planes[geometry::detsuperid[0][1]*2].costheta-upos1-v_win;
				vmax1 = vmin1+2*v_win;
				
				if(fulltrackbuilder[index].hitpairs_u1[k].first>=0){
					fulltrackbuilder[index].hitlist[nhits_U1] = fulltrackbuilder[index].hitpairs_u1[k].first;
					nhits_U1++;
				}
				if(fulltrackbuilder[index].hitpairs_u1[k].second>=0){
					fulltrackbuilder[index].hitlist[nhits_U1] = fulltrackbuilder[index].hitpairs_u1[k].second;
					nhits_U1++;
				}
				
				for(int l = 0; l<nv1; l++){
					nhits_V1 = nhits_U1;
					vpos1 = fulltrackbuilder[index].hitpairs_v1[l].second>=0 ? 0.5f*(ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_v1[l].first].pos+ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_v1[l].second].pos) : ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_v1[l].first].pos;
					
					if(vpos1<vmin1 || vpos1>vmax1)continue;
					
					if(fulltrackbuilder[index].hitpairs_v1[l].first>=0){
						fulltrackbuilder[index].hitlist[nhits_V1] = fulltrackbuilder[index].hitpairs_v1[l].first;
						nhits_V1++;
					}
					if(fulltrackbuilder[index].hitpairs_v1[l].second>=0){
						fulltrackbuilder[index].hitlist[nhits_V1] = fulltrackbuilder[index].hitpairs_v1[l].second;
						nhits_V1++;
					}
					
					//building here the global track candidate
					
					oc[index].AllTracklets[N_tkl].stationID = 6;
					oc[index].AllTracklets[N_tkl].nXHits = oc[index].AllTracklets[i].nXHits + nhits_X1;
					oc[index].AllTracklets[N_tkl].nUHits = oc[index].AllTracklets[i].nUHits + nhits_U1-nhits_X1;
					oc[index].AllTracklets[N_tkl].nVHits = oc[index].AllTracklets[i].nVHits + nhits_V1-nhits_U1;
					
					oc[index].AllTracklets[N_tkl].x0 = oc[index].AllTracklets[i].x0;
					oc[index].AllTracklets[N_tkl].err_x0 = oc[index].AllTracklets[i].err_x0;
					oc[index].AllTracklets[N_tkl].y0 = oc[index].AllTracklets[i].y0;
					oc[index].AllTracklets[N_tkl].err_y0 = oc[index].AllTracklets[i].err_y0;
					oc[index].AllTracklets[N_tkl].tx = oc[index].AllTracklets[i].tx;
					oc[index].AllTracklets[N_tkl].err_tx = oc[index].AllTracklets[i].err_tx;
					oc[index].AllTracklets[N_tkl].ty = oc[index].AllTracklets[i].ty;
					oc[index].AllTracklets[N_tkl].err_ty = oc[index].AllTracklets[i].err_ty;
					
					calculate_invP(oc[index].AllTracklets[N_tkl], fulltrackbuilder[index].TrackXZ_st1);
					
					nhits_1 = 0;
					for(int n = 0; n<oc[index].AllTracklets[i].nXHits+oc[index].AllTracklets[i].nUHits+oc[index].AllTracklets[i].nVHits; n++){
						oc[index].AllTracklets[N_tkl].hits[nhits_1] = oc[index].AllTracklets[i].hits[n];
						nhits_1++;
					}
					for(int n = 0; n<nhits_V1; n++){
						oc[index].AllTracklets[N_tkl].hits[nhits_1] = ic->AllHits[idxoff_global+fulltrackbuilder[index].hitlist[n]];
						nhits_1++;
					}
					
					resolve_leftright(oc[index].AllTracklets[N_tkl], planes, 75.);
					resolve_leftright(oc[index].AllTracklets[N_tkl], planes, 150.);
					resolve_single_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes);
					
					//refit after???
					
					N_tkl++;
				}
			}//end loop on u hits
		}//end loop on x hits
	}//end loop on straight tracks
	
	}else{
	
	float x_st1_mean[2];
	float x_st1_width[2];
	
	float xmin, xmax;
	
	short hyp;
	
	projid = 1;
	projid = 2;

	for(int i = 0; i<oc[index].nTracklets; i++){
		hyp = 0;
		extrapolate_track_position_st1(oc[index].AllTracklets[i], x_st1_mean, x_st1_width, planes, hyp);
		xmin = min(x_st1_mean[0]-x_st1_width[0], x_st1_mean[1]-x_st1_width[1]);
		xmax = max(x_st1_mean[0]+x_st1_width[0], x_st1_mean[1]+x_st1_width[1]);
		
		projid = 0;
		nx1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_x1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);
		
		if(nx1==0){
			hyp = 1;
			extrapolate_track_position_st1(oc[index].AllTracklets[i], x_st1_mean, x_st1_width, planes, hyp);
			xmin = min(x_st1_mean[0]-x_st1_width[0], x_st1_mean[1]-x_st1_width[1]);
			xmax = max(x_st1_mean[0]+x_st1_width[0], x_st1_mean[1]+x_st1_width[1]);
			projid = 0;
			nx1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_x1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);
		}
		
		if(nx1==0)continue;
		//triple loop on hits
		for(int j = 0; j<nx1; j++){
			nhits_X1 = 0;
			
			if(fulltrackbuilder[index].hitpairs_x1[j].first>=0){
				fitarrays->x_array[nhits_X1] = ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].pos;
				fitarrays->dx_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].detectorID].spacing/3.4641f;
				fitarrays->z_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].first].detectorID].z;
				fulltrackbuilder[index].hitlist[nhits_X1] = fulltrackbuilder[index].hitpairs_x1[j].first;
				nhits_X1++;
			}
			if(fulltrackbuilder[index].hitpairs_x1[j].second>=0){
				fitarrays->x_array[nhits_X1] = ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].pos;
				fitarrays->dx_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].detectorID].spacing/3.4641f;
				fitarrays->z_array[nhits_X1] = planes[ic->AllHits[idxoff_global+fulltrackbuilder[index].hitpairs_x1[j].second].detectorID].z;
				fulltrackbuilder[index].hitlist[nhits_X1] = fulltrackbuilder[index].hitpairs_x1[j].second;
				nhits_X1++;
			}
			
			//fit to obtain the st1 track segment: this time fitting from KMag, not from the origin:
			fitarrays->x_array[nhits_X1] = oc[index].AllTracklets[i].x0+oc[index].AllTracklets[i].tx*geometry::Z_KMAG_BEND;
			fitarrays->dx_array[nhits_X1] = oc[index].AllTracklets[i].err_x0+oc[index].AllTracklets[i].err_tx*geometry::Z_KMAG_BEND;
			fitarrays->z_array[nhits_X1] = geometry::Z_KMAG_BEND;
			
			fit_2D_track(nhits_X1+1, fitarrays->x_array, fitarrays->z_array, fitarrays->dx_array, fitarrays->A, fitarrays->Ainv, fitarrays->B, fitarrays->output_parameters, fitarrays->output_parameters_errors, fitarrays->chi2_2d);
			
			fulltrackbuilder[index].TrackXZ_st1.x_0 = fitarrays->output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.err_x_0 = fitarrays->output_parameters[0];
			fulltrackbuilder[index].TrackXZ_st1.tx_ = fitarrays->output_parameters[1];
			fulltrackbuilder[index].TrackXZ_st1.err_tx_ = fitarrays->output_parameters[1];

			projid = 1;
			find_xmin_xmax_in_chamber(xmin, xmax, fulltrackbuilder[index].TrackXZ_st1, stid, projid, planes);
			nu1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_u1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);
			
			projid = 2;
			find_xmin_xmax_in_chamber(xmin, xmax, fulltrackbuilder[index].TrackXZ_st1, stid, projid, planes);
			nv1 = make_hitpairs_in_station(ic, fulltrackbuilder[index].hitpairs_v1, fulltrackbuilder[index].hitidx1, fulltrackbuilder[index].hitidx2, fulltrackbuilder[index].hitflag1, fulltrackbuilder[index].hitflag2, stid, projid, planes, xmin, xmax);
			
			
			//
			for(int k = 0; k<nu1; k++){
				nhits_U1 = nhits_X1;
				if(fulltrackbuilder[index].hitpairs_u1[k].first>=0){
					fulltrackbuilder[index].hitlist[nhits_U1] = fulltrackbuilder[index].hitpairs_u1[k].first;
					nhits_U1++;
				}
				if(fulltrackbuilder[index].hitpairs_u1[k].second>=0){
					fulltrackbuilder[index].hitlist[nhits_U1] = fulltrackbuilder[index].hitpairs_u1[k].second;
					nhits_U1++;
				}
				
				for(int l = 0; l<nv1; l++){
					nhits_V1 = nhits_U1;

					if(fulltrackbuilder[index].hitpairs_v1[l].first>=0){
						fulltrackbuilder[index].hitlist[nhits_V1] = fulltrackbuilder[index].hitpairs_v1[l].first;
						nhits_V1++;
					}
					if(fulltrackbuilder[index].hitpairs_v1[l].second>=0){
						fulltrackbuilder[index].hitlist[nhits_V1] = fulltrackbuilder[index].hitpairs_v1[l].second;
						nhits_V1++;
					}
					//printf("nhits V1 = %d\n", nhits_V1-nhits_U1);

					//building here the global track candidate
					oc[index].AllTracklets[N_tkl].stationID = 6;
					oc[index].AllTracklets[N_tkl].nXHits = oc[index].AllTracklets[i].nXHits + nhits_X1;
					oc[index].AllTracklets[N_tkl].nUHits = oc[index].AllTracklets[i].nUHits + nhits_U1-nhits_X1;
					oc[index].AllTracklets[N_tkl].nVHits = oc[index].AllTracklets[i].nVHits + nhits_V1-nhits_U1;
					
					oc[index].AllTracklets[N_tkl].x0 = oc[index].AllTracklets[i].x0;
					oc[index].AllTracklets[N_tkl].err_x0 = oc[index].AllTracklets[i].err_x0;
					oc[index].AllTracklets[N_tkl].y0 = oc[index].AllTracklets[i].y0;
					oc[index].AllTracklets[N_tkl].err_y0 = oc[index].AllTracklets[i].err_y0;
					oc[index].AllTracklets[N_tkl].tx = oc[index].AllTracklets[i].tx;
					oc[index].AllTracklets[N_tkl].err_tx = oc[index].AllTracklets[i].err_tx;
					oc[index].AllTracklets[N_tkl].ty = oc[index].AllTracklets[i].ty;
					oc[index].AllTracklets[N_tkl].err_ty = oc[index].AllTracklets[i].err_ty;
					
					calculate_invP(oc[index].AllTracklets[N_tkl], fulltrackbuilder[index].TrackXZ_st1);
					
					nhits_1 = 0;
					for(int n = 0; n<oc[index].AllTracklets[i].nXHits+oc[index].AllTracklets[i].nUHits+oc[index].AllTracklets[i].nVHits; n++){
						oc[index].AllTracklets[N_tkl].hits[nhits_1] = oc[index].AllTracklets[i].hits[n];
						nhits_1++;
					}
					for(int n = 0; n<nhits_V1; n++){
						oc[index].AllTracklets[N_tkl].hits[nhits_1] = ic->AllHits[idxoff_global+fulltrackbuilder[index].hitlist[n]];
						nhits_1++;
					}

					resolve_leftright(oc[index].AllTracklets[N_tkl], planes, 75.);
					resolve_leftright(oc[index].AllTracklets[N_tkl], planes, 150.);
					
					//refit after???
					
					N_tkl++;
				}//end loop on v hits
			}//end loop on u hits
		}//end loop on x_hits
		//printf("nTracklets = %d \n", N_tkl);
	}//end loop on tracklets
	}	

	oc[index].nTracklets = N_tkl;
	
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

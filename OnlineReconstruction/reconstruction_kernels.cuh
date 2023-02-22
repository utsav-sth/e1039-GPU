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
		for(int k = 0; k<nvalues; k++)hitcolls->HitsChambersRawData[offset_hitcoll_prop+k] = hitarraycopy[k];

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
		for(int k = 0; k<nvalues; k++)hitcolls->HitsHodoRawData[offset_hitcoll_prop+k] = hitarraycopy[k];
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
	if(blockIdx.x+1==2044){
		int nhits_chambers_ = hitcolls->NHitsChambers[blockIdx.x*nChamberPlanes+detid_chambers[2]-1];
		gHits hitcoll_chambers = gHits(hitcolls->HitsChambersRawData, nhits_chambers_, offsets_hitcoll_chambers[2]);
		printf(" det offset %d array offset %d nhits_chambers(%d) = %d \n", blockIdx.x*nChamberPlanes+detid_chambers[2]-1, offsets_hitcoll_chambers[2], detid_chambers[2], nhits_chambers_);
		for(int i = 0; i<nhits_chambers_; i++){
			printf("%d %d %1.0f %1.4f %1.4f %1.0f %1.4f\n", offsets_hitcoll_chambers[2], detid_chambers[2], hitcoll_chambers.chan(i), hitcoll_chambers.pos(i), hitcoll_chambers.tdc(i), hitcoll_chambers.flag(i), hitcoll_chambers.drift(i));
			//printf("%d %d %d %d %d %d %d %d %d %1.0f %1.4f\n", ev, offsets_hitcoll_chambers[2], detid_chambers[2], i, offsets_hitcoll_chambers[2]+i, offsets_hitcoll_chambers[2]+i+nhits_chambers[2], offsets_hitcoll_chambers[2]+i+nhits_chambers[2]*2, offsets_hitcoll_chambers[2]+i+nhits_chambers*3, offsets_hitcoll_chambers[2]+i+nhits_chambers[2]*4, hitcolls->HitsChambersRawData[offsets_hitcoll_chambers[2]+i], hitcolls->HitsChambersRawData[offsets_hitcoll_chambers[2]+i+nhits_chambers[2]]);

		}}

		int nhits_hodo_ = hitcolls->NHitsHodo[blockIdx.x*nHodoPlanes+threadIdx.x];
		gHits hitcoll_hodo = gHits(hitcolls->HitsHodoRawData, nhits_hodo_, offsets_hitcoll_hodo[0]);
		printf(" det offset %d array offset %d nhits_hodo(%d) = %d \n", blockIdx.x*nChamberPlanes+threadIdx.x, offsets_hitcoll_hodo[0], detid_hodo[0], nhits_hodo_);
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
	
//#ifdef DEBUG
	
//#endif
}


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
	int bin0_st2 = (threadIdx.x/2)*7;
	int st3 = 3+threadIdx.x%2;//check d3p for even threads, d3m for odd threads...

	short hitflag1[100];
	short hitflag2[100];
	
	// As a starting point we will process 8 bins per thread: x acceptance divided by 4 * 2 bins for d3p and d3m 
	int nhitpairs_x2[nbins_st2];
	int nhitpairs_x3[29];
        //pairs in station 2
        thrust::pair<int, int> hitpairs_x2[nbins_st2*geometry::MaxHitsProj[0]];//7*10
        //pairs in station 3
        thrust::pair<int, int> hitpairs_x3[nbins_st3*geometry::MaxHitsProj[0]];//29*10

	
	unsigned int offset_hitcoll;
	
	short stid, projid, detid, detoff;
	
	projid = 0;
	stid = 2-1;
	detoff = 1;
	
	detid = geometry::detsuperid[stid][projid];
	int nhits_st2x;
	const gHits hits_st2x = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2x);
	const float z_st2x = z_array[detid];
	const float res_st2x = res_array[detid];
	
	detid-= 1;
	int nhits_st2xp;
	const gHits hits_st2xp = hitcolls->hitschambers(blockIdx.x, detid, nhits_st2xp);
	const float z_st2xp = z_array[detid];
	const float res_st2xp = res_array[detid];

	make_hitpairs_in_station_bins(hits_st2x, nhits_st2x, hits_st2xp, nhits_st2xp, hitpairs_x2, nhitpairs_x2, bin0_st2, nbins_st2, hitflag1, hitflag2, stid, projid);
	
	stid = st3-1;
	detid = geometry::detsuperid[stid][projid];
	int nhits_st3x;
	const gHits hits_st3x = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3x);
	const float z_st3x = z_array[detid];
	const float res_st3x = res_array[detid];

	detid-= 1;
	int nhits_st3xp;
	const gHits hits_st3xp = hitcolls->hitschambers(blockIdx.x, detid, nhits_st3xp);
	const float z_st3xp = z_array[detid];
	const float res_st3xp = res_array[detid];


	make_hitpairs_in_station_bins(hits_st3x, nhits_st3x, hits_st3xp, nhits_st3xp, hitpairs_x3, nhitpairs_x3, 0, nbins_st3, hitflag1, hitflag2, stid, projid);
	
	stid = 6-1;
	detid = geometry::detsuperid[stid][projid];
	int nhits_p1x1;
	const gHits hits_p1x1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1x1);
	const float z_p1x1 = z_array[detid];
	const float res_p1x1 = res_array[detid];

	detid-= 1;
	int nhits_p1x2;
	const gHits hits_p1x2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p1x2);
	const float z_p1x2 = z_array[detid];
	const float res_p1x2 = res_array[detid];
	
	stid = 7-1;
	detid = geometry::detsuperid[stid][projid];
	int nhits_p2x1;
	const gHits hits_p2x1 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2x1);
	const float z_p2x1 = z_array[detid];
	const float res_p2x1 = res_array[detid];

	detid-= 1;
	int nhits_p2x2;
	const gHits hits_p2x2 = hitcolls->hitsprop(blockIdx.x, detid, nhits_p2x2);
	const float z_p2x2 = z_array[detid];
	const float res_p2x2 = res_array[detid];
	
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
		//hastoomanyhits[blockIdx.x] = true;
	}
	
	float X[4];
	float errX[4];
	float Z[4];
	float A_[4];
	
	float Ainv_[4];
	float B_[2];
	float Par[2];
	float ParErr[2];
	float chi2;
	
	short nprop, iprop, idet;
	float xExp, ipos;
	
	short nhits_x;
	short nhits_x2, nhits_x3;
	
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
	
}

/*
__global__ void gKernel_XZ_tracking_prep(const gEventHitCollections* hitcolls, bool* hastoomanyhits, const gPlane* planes)
{
	if(hastoomanyhits[blockIdx.x]){
		return;
	}
	
	gTrackingXZparams param;
	

	const unsigned int detid_prop = 47+threadIdx.x;
	const unsigned int nhits_prop = hitcolls->NHitsPropTubes[blockIdx.x*nPropPlanes+threadIdx.x];
	const unsigned int offset_hitcoll_prop = blockIdx.x*datasizes::eventhitsize[2]+datasizes::NHitsParam*datasizes::NMaxHitsPropTubes*threadIdx.x;
		
	unsigned int detid, nhits, offset_hitcoll;
	
	short stid;
	short projid = 0;
	
	short detoff = 1;
	stid = 1;
	for(int i = 0; i<2; i++){
		detid = detsuperid[stid][projid]-i;
		nhits = hitcolls->NHitsPropTubes[blockIdx.x*nPropPlanes+detid-detoff];
		offset_hitcoll = blockIdx.x*datasizes::eventhitsize[0]+datasizes::NHitsParam*datasizes::NMaxHitsChambers*(detid-detoff);
		param.hits_st2x[i] = gHits(hitcolls->HitsChambersRawData, nhits, offset_hitcoll);

		param.z_array[i] = planes.z[detid];
		param.res_array[i] = planes.spacing[detid];
	}
	
	stid = 2;
	for(int i = 0; i<4; i++){
		if(i==2)stid+=1;
		detid = detsuperid[stid][projid]-i%2;
		nhits = hitcolls->NHitsPropTubes[blockIdx.x*nPropPlanes+detid-detoff];
		offset_hitcoll = blockIdx.x*datasizes::eventhitsize[0]+datasizes::NHitsParam*datasizes::NMaxHitsChambers*(detid-detoff);
		param.hits_st3x[i] = gHits(hitcolls->HitsChambersRawData, nhits, offset_hitcoll);

		param.z_array[i+2] = planes.z[detid];
		param.res_array[i+2] = planes.spacing[detid];
	}
	
	stid = 5;
	detoff = 47;
	for(int i = 0; i<4; i++){
		if(i==2)stid+=1;
		detid = detsuperid[stid][projid]-i%2;
		nhits = hitcolls->NHitsPropTubes[blockIdx.x*nPropPlanes+detoff];
		offset_hitcoll = blockIdx.x*datasizes::eventhitsize[2]+datasizes::NHitsParam*datasizes::NMaxHitsPropTubes*detoff;
		param.hits_px[i] = gHits(hitcolls->HitsPropTubesRawData, nhits, offset_hitcoll);

		param.z_array[i+6] = planes.z[detid];
		param.res_array[i+6] = planes.spacing[detid];
	}
	
	
	gKernel_XZ_tracking<<<1,THREADS_PER_BLOCK>>>(param);
}
*/

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

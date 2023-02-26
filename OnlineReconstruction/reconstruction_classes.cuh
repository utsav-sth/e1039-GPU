#include "reconstruction_constants.h"

//clone of LoadEvent::Hit:
struct gHit {
	public:
	//int index; // global hit index in the hit array
	short detectorID; // ID of the detector: one ID for each DC wire plane (30 total), hodoscope plane (16 total), proportional tube plane (8 total).
	short elementID; // ID of the element in the detector: wire/slat/tube number
	float tdcTime; // raw TDC time from the DAQ 
	float driftDistance; // calculated drift distance from RT profile (supplied in database) IF tdcTime between tmin and tmax defined for detector; 
	//short sign_mc;//temp
	float pos; // position in the projection of the detector (e.g. X in a X plane, etc)
	//short flag; // 1: in time; 2: hodo mask; 3: trigger mask
};

//it may be beneficial to have several classes of hits...
struct gHits {
	public:
	const unsigned int NHitsTotal;
	float* m_hitdata;
	
	//convention: offset: chan (element ID) 0; pos 1; tdc 2; flag 3; drift 4;
		
	__host__ __device__ gHits(float* basedata, const unsigned total_number_of_hits, const unsigned offset = 0) :
    		m_hitdata(reinterpret_cast<float*>(basedata) + offset), NHitsTotal(total_number_of_hits)
		{
			static_assert(sizeof(float) == sizeof(unsigned));
			assert((((size_t) basedata) & sizeof(float)) == 0);
		}
		
	__host__ __device__ inline float chan(const unsigned index) const
		{
			assert(index < NHitsTotal);
			return m_hitdata[index];
		}
	
	__host__ __device__ inline float pos(const unsigned index) const
		{
			assert(index < NHitsTotal);
			return m_hitdata[NHitsTotal + index];
		}
		
	__host__ __device__ inline float tdc(const unsigned index) const
		{
			assert(index < NHitsTotal);
			return m_hitdata[NHitsTotal*2 + index];
		}
		
	__host__ __device__ inline float flag(const unsigned index) const
		{
			assert(index < NHitsTotal);
			return m_hitdata[NHitsTotal*3 + index];
		}
	
	__host__ __device__ inline float drift(const unsigned index) const
		{
			assert(index < NHitsTotal);
			return m_hitdata[NHitsTotal*4 + index];
		}
};


struct gTracklet {
      public:
      __device__ gTracklet(){
	nHits = 0;
      }
	            
      short stationID;//0
      short threadID;//1
      short nHits;//2
      float chisq;//3
      float chisq_vtx;//4

      
      float tx;//5
      float ty;//6
      float x0;//7
      float y0;//8
      float invP;//9
      
      float err_tx;//10
      float err_ty;//11
      float err_x0;//12
      float err_y0;//13
      float err_invP;//14
      
      short charge;//15

      //maybe we can be a bit more sober in memory, and just require hit "event" index?
      gHit hits[datasizes::MaxHitsPerTrack];// array of all hits:16-105
      short hitsign[datasizes::MaxHitsPerTrack];//106-123
      float residual[datasizes::MaxHitsPerTrack];//124-141
};

struct gEvent {
	public:
	//int RunID[EstnEvtMax]; // Run Number
	int EventID[EstnEvtMax]; // Event number
	//int SpillID[EstnEvtMax]; // Spill number
	int TriggerBits[EstnEvtMax]; // hash of the trigger bits: 0-4: MATRIX1-5; 5-9: NIM1-5;
	short TargetPos[EstnEvtMax]; // target position: proxy for target ID?
	//int TurnID[EstnEvtMax]; // => related to beam intensity
	//int RFID[EstnEvtMax]; // => related to beam intensity
	//int Intensity[EstnEvtMax*33]; //  16 before, one onset, and 16 after
	//short TriggerEmu[EstnEvtMax]; // 1 if MC event
	short NRoads[EstnEvtMax*4]; // 0, positive top; 1, positive bottom; 2, negative top; 3, negative bottom
	int NHits[EstnEvtMax*nDetectors]; // number of hits in each detector plane
	int nAH[EstnEvtMax]; // size of AllHits
	int nTH[EstnEvtMax]; // size of TriggerHits
	//limit the max hit multiplicity for unreduced events to 2 times of what it is for reduced events
	//float HitsChambersReducedData[EstnEvtMax*nChamberPlanes*5*datasizes::NMaxHitsChambers];
	//float HitsPropTubesReducedData[EstnEvtMax*nChamberPlanes*5*datasizes::NMaxHitsPropTubes*2];
	//float HitsHodoReducedData[EstnEvtMax*nChamberPlanes*4*datasizes::NMaxHitsHodoscopes];
	//gHit AllHits[EstnEvtMax*EstnAHMax]; // array of all hits
	//gHit TriggerHits[EstnEvtMax*EstnTHMax]; // array of trigger hits
	int nTracklets[EstnEvtMax];
	bool HasTooManyHits[EstnEvtMax];//bool to flag an event with too many hits
};

struct gEventHitCollections {
	public:
	//TODO: add offset calculation functions!!!
	
	unsigned int NHitsChambers[EstnEvtMax*nChamberPlanes];
	float HitsChambersRawData[EstnEvtMax*nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsChambers];
	
	unsigned int NHitsHodo[EstnEvtMax*nHodoPlanes]; 
	float HitsHodoRawData[EstnEvtMax*nHodoPlanes*datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes];
	
	unsigned int NHitsPropTubes[EstnEvtMax*nPropPlanes]; 
	float HitsPropTubesRawData[EstnEvtMax*nPropPlanes*datasizes::NHitsParam*datasizes::NMaxHitsPropTubes];
	
	__device__ const gHits hitschambers(const unsigned int event, const short detid, int& nhits) {
		nhits = NHitsChambers[event*nChamberPlanes+detid-1];
		return gHits(HitsChambersRawData, nhits, event*datasizes::eventhitsize[0]+datasizes::NHitsParam*datasizes::NMaxHitsChambers*(detid-1) );
	}

	__device__ const gHits hitshodos(const unsigned int event, const short detid, int& nhits){
		nhits = NHitsHodo[event*nHodoPlanes+detid-nChamberPlanes-1];
		return gHits(HitsHodoRawData, nhits, event*datasizes::eventhitsize[1]+datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes*(detid-nChamberPlanes-1) );
	}

	__device__ const gHits hitsprop(const unsigned int event, const short detid, int& nhits){
		nhits = NHitsPropTubes[event*nPropPlanes+detid-nChamberPlanes-nHodoPlanes-1];
		return gHits(HitsPropTubesRawData, nhits, event*datasizes::eventhitsize[2]+datasizes::NHitsParam*datasizes::NMaxHitsPropTubes*(detid-nChamberPlanes-nHodoPlanes-1) );
	}

};

struct gTracks {
	public:
	const unsigned int NTracksTotal;
	const unsigned int TrackSize;
	float* m_trackdata;
	
	__host__ __device__ gTracks(float* basedata, const unsigned total_number_of_tracks, const unsigned offset = 0) :
		m_trackdata(reinterpret_cast<float*>(basedata) + offset), NTracksTotal(total_number_of_tracks), TrackSize(datasizes::NTracksParam)
		{
			static_assert(sizeof(float) == sizeof(unsigned));
			assert((((size_t) basedata) & sizeof(float)) == 0);
		}

	__host__ __device__ inline float stationID(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index];
		}

	__host__ __device__ inline float threadID(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+1];
		}

	__host__ __device__ inline float nHits(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+2];
		}
	
	__host__ __device__ inline float chisq(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+3];
		}
	
	__host__ __device__ inline float chisq_vtx(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+4];
		}

	//track parameters
	__host__ __device__ inline float tx(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+5];
		}

	__host__ __device__ inline float ty(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+6];
		}

	__host__ __device__ inline float x0(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+7];
		}

	__host__ __device__ inline float y0(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+8];
		}

	__host__ __device__ inline float invP(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+9];
		}

	__host__ __device__ inline float err_tx(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+10];
		}

	__host__ __device__ inline float err_ty(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+11];
		}

	__host__ __device__ inline float err_x0(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+12];
		}

	__host__ __device__ inline float err_y0(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+13];
		}

	__host__ __device__ inline float err_invP(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+14];
		}

	__host__ __device__ inline float charge(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+15];
		}

	__host__ __device__ inline float hits_detid(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+16 + ihit ];
		}

	__host__ __device__ inline float hits_elid(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+34 + ihit ];
		}

	__host__ __device__ inline float hits_pos(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+52 + ihit ];
		}

	__host__ __device__ inline float hits_tdc(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+70 + ihit ];
		}

	__host__ __device__ inline float hits_drift(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+88 + ihit ];
		}
	
	__host__ __device__ inline float hits_sign(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+106 + ihit ];
		}
	
	__host__ __device__ inline float hits_residual(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+124 + ihit ];
		}
	
	
};

struct gTracklets {
	public:
	const unsigned int ntkl;
	gTracklet* tkl_list;
	
	
	__host__ __device__ gTracklets(gTracklet* baselist, const unsigned total_number_of_tracks, const unsigned offset = 0) :
		tkl_list(baselist + offset), ntkl(total_number_of_tracks)
		{
			static_assert(sizeof(float) == sizeof(unsigned));
			assert((((size_t) baselist) & sizeof(gTracklet)) == 0);
		}
	
	__host__ __device__ inline gTracklet getTracklet(const unsigned index) const
		{
			assert(index < ntkl);
			return tkl_list[index];
		}
		
	__host__ __device__ inline void setTracklet(const unsigned index, gTracklet tkl) const
		{
			assert(index < ntkl);
			tkl_list[index] = tkl;
		}
};

struct gEventTrackCollection{
	unsigned short NTracks[EstnEvtMax*THREADS_PER_BLOCK];
	float TracksRawData[EstnEvtMax*datasizes::TrackletSizeMax*datasizes::NTracksParam];
	//gTracklet Tracklets[EstnEvtMax*datasizes::TrackletSizeMax];
	//__device__ const gTracklets tracks(const unsigned int event, int& ntracks) {
	//	ntracks = NTracks[event];
	//	return gTracklets(Tracklets, ntracks, event*datasizes::TrackletSizeMax);
	//}

	__device__ const gTracks tracks(const unsigned int event, const unsigned threadID, unsigned int& ntracks) {
		ntracks = NTracks[event*THREADS_PER_BLOCK+threadID];
		return gTracks(TracksRawData, ntracks, event*datasizes::TrackletSizeMax*datasizes::NTracksParam+threadID*datasizes::TrackletSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK);
	}
	
	__device__ void setStationID(const unsigned int evt_offset, const unsigned int itrack, const float stid) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam] = stid;
	}
	__device__ void setThreadID(const unsigned int evt_offset, const unsigned int itrack, const float threadid) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+1] = threadid;
	}
	__device__ void setnHits(const unsigned int evt_offset, const unsigned int itrack, const float nhits) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+2] = nhits;
	}
	__device__ void setChisq(const unsigned int evt_offset, const unsigned int itrack, const float chisq) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+3] = chisq;
	}
	__device__ void setChisqVtx(const unsigned int evt_offset, const unsigned int itrack, const float chisq_vtx) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+4] = chisq_vtx;
	}
	//track parameters
	__device__ void setTx(const unsigned int evt_offset, const unsigned int itrack, const float tx) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+5] = tx;
	}
	__device__ void setTy(const unsigned int evt_offset, const unsigned int itrack, const float ty) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+6] = ty;
	}
	__device__ void setX0(const unsigned int evt_offset, const unsigned int itrack, const float x0) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+7] = x0;
	}
	__device__ void setY0(const unsigned int evt_offset, const unsigned int itrack, const float y0) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+8] = y0;
	}
	__device__ void setinvP(const unsigned int evt_offset, const unsigned int itrack, const float invp) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+9] = invp;
	}
	__device__ void setErrTx(const unsigned int evt_offset, const unsigned int itrack, const float errtx) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+10] = errtx;
	}
	__device__ void setErrTy(const unsigned int evt_offset, const unsigned int itrack, const float errty) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+11] = errty;
	}
	__device__ void setErrX0(const unsigned int evt_offset, const unsigned int itrack, const float errx0) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+12] = errx0;
	}
	__device__ void setErrY0(const unsigned int evt_offset, const unsigned int itrack, const float erry0) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+13] = erry0;
	}
	__device__ void setErrinvP(const unsigned int evt_offset, const unsigned int itrack, const float errinvp) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+14] = errinvp;
	}
	
	/*
	__device__ void setStationID(const unsigned int event, const unsigned int itrack, const float stid) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack] = stid;
	}
	__device__ void setThreadID(const unsigned int event, const unsigned int itrack, const float threadid) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]] = threadid;
	}
	__device__ void setnHits(const unsigned int event, const unsigned int itrack, const float nhits) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*2] = nhits;
	}
	__device__ void setChisq(const unsigned int event, const unsigned int itrack, const float chisq) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*3] = chisq;
	}
	__device__ void setChisqVtx(const unsigned int event, const unsigned int itrack, const float chisq_vtx) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*4] = chisq_vtx;
	}
	//track parameters
	__device__ void setTx(const unsigned int event, const unsigned int itrack, const float tx) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*5] = tx;
	}
	__device__ void setTy(const unsigned int event, const unsigned int itrack, const float ty) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*6] = ty;
	}
	__device__ void setX0(const unsigned int event, const unsigned int itrack, const float x0) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*7] = x0;
	}
	__device__ void setY0(const unsigned int event, const unsigned int itrack, const float y0) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*8] = y0;
	}
	__device__ void setinvP(const unsigned int event, const unsigned int itrack, const float invp) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*9] = invp;
	}
	__device__ void setErrTx(const unsigned int event, const unsigned int itrack, const float errtx) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*10] = errtx;
	}
	__device__ void setErrTy(const unsigned int event, const unsigned int itrack, const float errty) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*11] = errty;
	}
	__device__ void setErrX0(const unsigned int event, const unsigned int itrack, const float errx0) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*12] = errx0;
	}
	__device__ void setErrY0(const unsigned int event, const unsigned int itrack, const float erry0) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*13] = erry0;
	}
	__device__ void setErrinvP(const unsigned int event, const unsigned int itrack, const float errinvp) {
		TracksRawData[event*datasizes::TrackletSizeMax*datasizes::NTracksParam+itrack*NTracks[event]*14] = errinvp;
	}
	*/
};


//geometry carrier
struct gPlane {
      public:
      float z[nDetectors];
      int nelem[nDetectors];
      float cellwidth[nDetectors];
      float spacing[nDetectors];
      float xoffset[nDetectors];
      float scalex[nDetectors];
      float x0[nDetectors];
      float x1[nDetectors];
      float x2[nDetectors];
      float costheta[nDetectors];
      float scaley[nDetectors];
      float y0[nDetectors];
      float y1[nDetectors];
      float y2[nDetectors];
      float sintheta[nDetectors];
      float resolution[nDetectors];
      float deltaW_[nDetectors*9];
      float z_mean[nDetectors];
      float u_win[nDetectors];
      float v_win_fac1[nDetectors];
      float v_win_fac2[nDetectors];
      float v_win_fac3[nDetectors];
      float p1x_w1[nDetectors];
      float p1y_w1[nDetectors];
      float p1z_w1[nDetectors];
      float deltapx[nDetectors];
      float deltapy[nDetectors];
      float deltapz[nDetectors];
      float dp1x[nDetectors];
      float dp1y[nDetectors];
      float dp1z[nDetectors];
      float slope_max[nDetectors];
      float inter_max[nDetectors];
};

#ifdef OLDCODE

struct gTrack2D {
      public:
      // note: x_ is to be understood as either x or y...
      float tx_;
      float x_0;
      
      float err_tx_;
      float err_x_0;
      
      float chisq;
              
      short nhits;
      int hitlist[8];
      short hitsign[8];
      
};

struct gTrackingXZparams{
	public:
	gTrackingXZparams();
	//utils
	float z_array[10];
	float res_array[10];
	
	//input
	gHits hits_st2x;
	gHits hits_st2xp;
	gHits hits_st3px;
	gHits hits_st3pxp;
	gHits hits_st3mx;
	gHits hits_st3mxp;

	gHits hits_p1x1;
	gHits hits_p1x2;
	gHits hits_p2x1;
	gHits hits_p2x2;
	gHits hits_st2x[2];
	gHits hits_st3x[4];
	gHits hits_px[4];
	
	//output
	gTracklets Tracks;
};

struct gOutputEvent {
public:
	int EventID[EstnEvtMax];
	int nAH[EstnEvtMax];
	bool HasTooManyHits[EstnEvtMax];//bool to flag an event with too many hits
	int nTracklets[EstnEvtMax];
//	gTracklet AllTracklets[EstnEvtMax*TrackletSizeMax];
};
#endif


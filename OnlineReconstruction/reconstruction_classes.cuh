#include "reconstruction_constants.h"

//clone of LoadEvent::Hit:
struct gHit {
	public:
	int index; // global hit index in the hit array
	short detectorID; // ID of the detector: one ID for each DC wire plane (30 total), hodoscope plane (16 total), proportional tube plane (8 total).
	short elementID; // ID of the element in the detector: wire/slat/tube number
	float tdcTime; // raw TDC time from the DAQ 
	float driftDistance; // calculated drift distance from RT profile (supplied in database) IF tdcTime between tmin and tmax defined for detector; 
	short sign_mc;//temp
	float pos; // position in the projection of the detector (e.g. X in a X plane, etc)
	short flag; // 1: in time; 2: hodo mask; 3: trigger mask
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
	            
      short stationID;
      short nHits;
      float chisq;
      float chisq_vtx;

      
      float tx;
      float ty;
      float x0;
      float y0;
      float invP;
      
      float err_tx;
      float err_ty;
      float err_x0;
      float err_y0;
      float err_invP;
      
      short charge;

      //maybe we can be a bit more sober in memory, and just require hit "event" index?
      gHit hits[datasizes::MaxHitsPerTrack];// array of all hits
      short hitsign[datasizes::MaxHitsPerTrack];
      float residual[datasizes::MaxHitsPerTrack];
};

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

struct gEvent {
	public:
	int RunID[EstnEvtMax]; // Run Number
	int EventID[EstnEvtMax]; // Event number
	int SpillID[EstnEvtMax]; // Spill number
	int TriggerBits[EstnEvtMax]; // hash of the trigger bits: 0-4: MATRIX1-5; 5-9: NIM1-5;
	short TargetPos[EstnEvtMax]; // target position: proxy for target ID?
	int TurnID[EstnEvtMax]; // => related to beam intensity
	int RFID[EstnEvtMax]; // => related to beam intensity
	int Intensity[EstnEvtMax*33]; //  16 before, one onset, and 16 after
	short TriggerEmu[EstnEvtMax]; // 1 if MC event
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
	
	__device__ const gHits hitschambers(const unsigned int event, const short detid, int &nhits) {
		nhits = NHitsChambers[event*nChamberPlanes+detid-1];
		return gHits(HitsChambersRawData, nhits, event*datasizes::eventhitsize[0]+datasizes::NHitsParam*datasizes::NMaxHitsChambers*(detid-1) );
	}

	__device__ const gHits hitshodos(const unsigned int event, const short detid, int &nhits){
		nhits = NHitsHodo[event*nHodoPlanes+detid-31];
		return gHits(HitsHodoRawData, nhits, event*datasizes::eventhitsize[1]+datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes*(detid-31) );
	}

	__device__ const gHits hitsprop(const unsigned int event, const short detid, int &nhits){
		nhits = NHitsPropTubes[event*nPropPlanes+detid-47];
		return gHits(HitsPropTubesRawData, nhits, event*datasizes::eventhitsize[2]+datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes*(detid-47) );
	}

};

/*
struct gTracks {
	public:
	const unsigned int NTracksTotal;
	float* m_trackdata;
	float* m_hitdata;
	
	__host__ __device__ gTracks(float* basedata, float* hitdata, const unsigned total_number_of_tracks, const unsigned offset = 0) :
		m_trackdata(basedata + offset), m_hitdata(hitdata + offset), NTracksTotal(total_number_of_tracks)
		{
			static_assert(sizeof(float) == sizeof(unsigned));
			assert((((size_t) basedata) & sizeof(float)) == 0);
		}
};
*/

struct gEventTrackCollection{
	unsigned int NTracks[EstnEvtMax];
	gTracklet Tracklets[EstnEvtMax*datasizes::TrackletSizeMax];
};


struct gTracklets {
	public:
	const unsigned int ntkl;
	gTracklet* tkl_list;
	
	__host__ __device__ gTracklets(gTracklet* baselist, gTracklet* tkl_list, const unsigned total_number_of_tracks, const unsigned offset = 0) :
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

/*
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
*/


struct gOutputEvent {
public:
	int EventID[EstnEvtMax];
	int nAH[EstnEvtMax];
	bool HasTooManyHits[EstnEvtMax];//bool to flag an event with too many hits
	int nTracklets[EstnEvtMax];
//	gTracklet AllTracklets[EstnEvtMax*TrackletSizeMax];
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

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
struct gWCHits {
	public:
	const unsigned int NHitsTotal;
	float* m_hitdata;
	
	//convention: offset: chan (element ID) 0; pos 1; tdc 2; flag 3; drift 4;
		
	__host__ __device__ gWCHits(float* basedata, const unsigned total_number_of_hits, const unsigned offset = 0) :
    		m_hitdata(basedata + offset), NHitsTotal(total_number_of_hits)
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

//it may be beneficial to have several classes of hits...
struct gHodoHits {
	public:
	const unsigned int NHitsTotal;
	float* m_hitdata;
	
	//convention: offset: chan (element ID) 0; pos 1; tdc 2; flag 3;
	//maybe we could revise that?
	
	__host__ __device__ gHodoHits(float* basedata, const unsigned total_number_of_hits, const unsigned offset = 0) :
    		m_hitdata(basedata + offset), NHitsTotal(total_number_of_hits)
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
};

struct gHitPairs{
	
};



struct gTracklet {
      public:
      __device__ gTracklet(){
	nXHits = nUHits = nVHits = 0;
      }
      __device__ short nHits(){
        return nXHits + nUHits + nVHits;
      } 
            
      short stationID;
      short nXHits;
      short nUHits;
      short nVHits;

      float chisq;
      float chisq_vtx;

      //maybe we can be a bit more sober in memory, and just require hit "event" index?
      gHit hits[MaxHitsPerTrack];// array of all hits
      short hitsign[MaxHitsPerTrack];
      
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
      float residual[MaxHitsPerTrack];
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
	
	//float HitsChambersRawData[EstnEvtMax*nChamberPlanes*5*datasizes::NMaxHitsChambers*2];
	//float HitsPropTubesRawData[EstnEvtMax*nChamberPlanes*5*datasizes::NMaxHitsPropTubes*2];
	//float HitsHodoRawData[EstnEvtMax*nChamberPlanes*4*datasizes::NMaxHitsHodoscopes*2];
	gHit AllHits[EstnEvtMax*EstnAHMax]; // array of all hits
	gHit TriggerHits[EstnEvtMax*EstnTHMax]; // array of trigger hits
	bool HasTooManyHits[EstnEvtMax];//bool to flag an event with too many hits
};

// one per chamber?
struct gEventHitCollections {
	public:
	unsigned int NHitsChambers[EstnEvtMax*nChamberPlanes];
	float HitsChambersReducedData[EstnEvtMax*nChamberPlanes*5*datasizes::NMaxHitsChambers];
	
	unsigned int NHitsPropTubes[EstnEvtMax*nChamberPlanes]; 
	float HitsPropTubesReducedData[EstnEvtMax*nChamberPlanes*5*datasizes::NMaxHitsPropTubes];

	unsigned int NHitsHodo[EstnEvtMax*nChamberPlanes]; 
	float HitsHodoReducedData[EstnEvtMax*nChamberPlanes*4*datasizes::NMaxHitsHodoscopes];
};

struct gEventReducerInput{
	
	
};


struct gTracks2D {
	public:
	const unsigned int NTracksTotal;
	float* m_trackdata;
	float* m_hitdata;
	
	__host__ __device__ gTracks2D(float* basedata, float* hitdata, const unsigned total_number_of_tracks, const unsigned offset = 0) :
		m_trackdata(basedata + offset), m_hitdata(hitdata + offset), NTracksTotal(total_number_of_tracks)
		{
			static_assert(sizeof(float) == sizeof(unsigned));
			assert((((size_t) basedata) & sizeof(float)) == 0);
		}
	
};



/*
struct gChamberHitColl{
	public:
	int stID;
	int projID;//X: 0, U: 1, V: 2
	thrust::pair<gHit, gHit> hitpairs[datasizes::NMaxHitPairsChambers];
};

struct gHodoHitColl{
	public:
	int stID;
	int projID;//X, 0; Y: 1
	gHit Hits[datasizes::NMaxHitsHodoscopes];
};

struct gPropHitColl{
	public:
	int stID;
	int projID;//X, 0; Y: 1
	gHit Hits[datasizes::NMaxHitsPropTubes];
};
*/

struct gFullTrackBuilder{
public:
	gTrack2D TrackXZ_st1[EstnEvtMax];
	
      	int hitlist[EstnEvtMax*MaxHitsPerTrack];
      	short hitsign[EstnEvtMax*MaxHitsPerTrack];
	
	thrust::pair<int, int> hitpairs_x1[EstnEvtMax*100];
	thrust::pair<int, int> hitpairs_u1[EstnEvtMax*100];
	thrust::pair<int, int> hitpairs_v1[EstnEvtMax*100];
	
	//util arrays for pair making
//	int hitidx1[EstnEvtMax*100];
//	int hitidx2[EstnEvtMax*100];
//	short hitflag1[EstnEvtMax*100];
//	short hitflag2[EstnEvtMax*100];
};


struct gStraightTrackBuilder{
public:
//	gTrack2D trackXZ[EstnEvtMax];
//	gTrack2D trackYZ[EstnEvtMax];
//	gTrack2D besttrackYZ[EstnEvtMax];
	
        //pairs in station 2
        thrust::pair<int, int> hitpairs_x2[EstnEvtMax*280];//28*10
        thrust::pair<int, int> hitpairs_u2[EstnEvtMax*1120];//28*40
        thrust::pair<int, int> hitpairs_v2[EstnEvtMax*1120];
        //pairs in station 3
        thrust::pair<int, int> hitpairs_x3[EstnEvtMax*580];//29*10*2
        thrust::pair<int, int> hitpairs_u3[EstnEvtMax*2320];//29*40*2
        thrust::pair<int, int> hitpairs_v3[EstnEvtMax*2320];
	
//	int nhitpairs_x2[EstnEvtMax*28];
//	int nhitpairs_u2[EstnEvtMax*28];
//	int nhitpairs_v2[EstnEvtMax*28];
	
//	int nhitpairs_x3[EstnEvtMax*58];
//	int nhitpairs_u3[EstnEvtMax*58];
//	int nhitpairs_v3[EstnEvtMax*58];
	
//	//util arrays for pair making
//	int hitidx1[EstnEvtMax*100];
//	int hitidx2[EstnEvtMax*100];
//	short hitflag1[EstnEvtMax*100];
//	short hitflag2[EstnEvtMax*100];
};

struct gStraightFitArrays {
public:
      int npoints[EstnEvtMax];
      float drift_dist[EstnEvtMax*MaxHitsPerTrack]; // hit drift distance
      float resolution[EstnEvtMax*MaxHitsPerTrack]; // detector resolution
      
      float p1x[EstnEvtMax*MaxHitsPerTrack];// x bottom end point of the wire hit 
      float p1y[EstnEvtMax*MaxHitsPerTrack];// y bottom end point of the wire hit 
      float p1z[EstnEvtMax*MaxHitsPerTrack];// z bottom end point of the wire hit 
      
      float deltapx[EstnEvtMax*MaxHitsPerTrack];// x distance between bottom and top end points of the wire hit 
      float deltapy[EstnEvtMax*MaxHitsPerTrack];// y distance between bottom and top end points of the wire hit 
      float deltapz[EstnEvtMax*MaxHitsPerTrack];// z distance between bottom and top end points of the wire hit 
      
//      float output_parameters[EstnEvtMax*4];
//      float output_parameters_errors[EstnEvtMax*4];
//      float chi2_2d[EstnEvtMax];
//      float chi2[EstnEvtMax];

      float x_array[EstnEvtMax*MaxHitsPerTrack];// x position arrays
      float y_array[EstnEvtMax*MaxHitsPerTrack];// y position arrays
      float z_array[EstnEvtMax*MaxHitsPerTrack];// z position arrays
      float dx_array[EstnEvtMax*MaxHitsPerTrack];// x position uncertainty
      float dy_array[EstnEvtMax*MaxHitsPerTrack];// x position uncertainty
      
//      float A[EstnEvtMax*4];// matrix: max size 2x2
//      float Ainv[EstnEvtMax*4];// inverted matrix
//      float B[EstnEvtMax*2];// input vector
};

struct gKalmanFitArrays{
public:
	float state[EstnEvtMax*5];// 5-vector: x0, y0, tx, ty, invP
	float Cov[EstnEvtMax*25];// symmetric 5x5 matrix: C00 = err_x0, C11 = err_y0, C22 = err_tx, C33 = err_ty, C44 = err_invP
	float H[EstnEvtMax*2];
	float K[EstnEvtMax*5];// matrix: max size 5x5, but we can use this unique array for all possible sizes
	float KCResKt[EstnEvtMax*25];// matrix 5x5, result of tensor product of K*K
	float chi2[EstnEvtMax];// chi2
};

struct gOutputEvent {
public:
	int EventID[EstnEvtMax];
	int nAH[EstnEvtMax];
	bool HasTooManyHits[EstnEvtMax];//bool to flag an event with too many hits
	int nTracklets[EstnEvtMax];
	gTracklet AllTracklets[EstnEvtMax*TrackletSizeMax];
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

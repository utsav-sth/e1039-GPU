#define nChamberPlanes 30
#define nHodoPlanes 16
#define nPropPlanes 8
#define nDetectors (nChamberPlanes+nHodoPlanes+nPropPlanes)
#define Epsilon 0.00001f

#define triggerBit(n) (1 << (n))
#define hitFlagBit(n) (1 << (n))

using namespace std;

const int EstnEvtMax = 10240;
const int THREADS_PER_BLOCK = 512;
int BLOCKS_NUM = EstnEvtMax/THREADS_PER_BLOCK;
const int EstnAHMax = 5000;
const int EstnTHMax = 200;
const int ClusterSizeMax = 100;
const int Track2DSizeMax = 100;
const int TrackletSizeMax = 500;
const int MaxHitsPerTrack = 18;

const double TX_MAX = 0.15;
const double TY_MAX = 0.1;
const double X0_MAX = 150;
const double Y0_MAX = 50;
const double INVP_MAX = 0.2;
const double INVP_MIN = 0.01;


namespace geometry{
	__device__ constexpr float spacingplane[28] = {0., 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 1.3, 1.3, 1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 4.0, 4.0, 7.0, 7.0, 8.0, 12.0, 12.0, 10.0, 3.0, 3.0, 3.0, 3.0};
	__device__ constexpr short detsuperid[7][3] = {{2, 1, 3}, {5, 6, 4}, {8, 9, 7}, {11, 12, 10}, {14, 15, 13}, {25, 26, -1}, {24, 27, -1}}; 
	__device__ constexpr short planetype[15] = {1, 0, 2, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1};
	__device__ constexpr short detmap_zincr[31] = {-1,  0,  1,  3,  2,  4,  5,  0,  1,  3,  2,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 12, 13, 14, 15, 16, 17};
	__device__ constexpr short dets_x[6] = {15, 16, 21, 22, 27, 28};
	__device__ constexpr short hodoplanerange[5][2] = {{31, 34}, {31, 34}, {35, 38}, {39, 40}, {41, 46}};// range of planes to look for hits
	__device__ constexpr float hodofudgefac[5] = {0.25, 0.25, 0.2, 0.15, 0.0};
	__device__ constexpr float Z_TARGET = -300;
	__device__ constexpr float Z_DUMP = 42;
	__device__ constexpr float SAGITTA_TARGET_CENTER = 1.85;
	__device__ constexpr float SAGITTA_TARGET_WIDTH = 0.25;
	__device__ constexpr float SAGITTA_DUMP_CENTER = 1.5;
	__device__ constexpr float SAGITTA_DUMP_WIDTH = 0.3;
	__device__ constexpr float PT_KICK_KMAG = 0.4016;
	__device__ constexpr float Z_KMAG_BEND = 1064.26;
	__device__ constexpr short lrpossibility[4][2] = {{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
	__device__ constexpr short N_WCHitsBins[4] = {32, 28, 29, 29};
	__device__ constexpr short MaxHitsProj[3] = {10, 40, 40};
	__device__ constexpr short WCHitsBins[4][3][2][32] = {
        	{{{1, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86, 91, 96, 101, 106, 111, 116, 121, 126, 131, 136, 141, 146, 151, 156}, 
		  {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160}}, // d0x
        	 {{1, 1, 1, 1, 4, 9, 13, 18, 23, 28, 33, 38, 43, 47, 52, 57, 62, 67, 72, 76, 81, 86, 91, 96, 101, 106, 110, 115, 120, 125, 130, 135}, 
		  {66, 71, 76, 81, 86, 91, 95, 100, 105, 110, 115, 120, 124, 129, 134, 139, 144, 149, 154, 158, 163, 168, 173, 178, 183, 188, 192, 197, 201, 201, 201, 201}}, // d0u
        	 {{1, 1, 1, 3, 8, 13, 18, 22, 27, 32, 37, 42, 47, 52, 56, 61, 66, 71, 76, 81, 86, 90, 95, 100, 105, 110, 115, 120, 124, 129, 134, 139, }, 
		  {63, 68, 73, 78, 83, 87, 92, 97, 102, 107, 112, 117, 121, 126, 131, 136, 141, 146, 151, 155, 160, 165, 170, 175, 180, 184, 189, 194, 199, 201, 201, 201}}}, // d0v
        	{{{1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 109, -1, -1, -1, -1}, 
		  {4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, -1, -1, -1, -1}}, // d2x
        	 {{1, 1, 1, 1, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, -1, -1, -1, -1}, 
		  {31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99, 103, 107, 111, 115, 119, 123, 127, 128, 128, 128, -1, -1, -1, -1}}, // d2u
        	 {{1, 1, 1, 3, 7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99, -1, -1, -1, -1}, 
		  {31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99, 103, 107, 111, 115, 119, 123, 127, 128, 128, 128, -1, -1, -1, -1, }}}, // d2v
        	{{{1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 109, 113, -1, -1, -1}, 
		  {4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, -1, -1, -1}}, // d3px
        	 {{1, 4, 7, 11, 15, 19, 23, 27, 31, 35, 38, 42, 46, 50, 54, 58, 62, 66, 70, 73, 77, 81, 85, 89, 93, 97, 101, 104, 108, -1, -1, -1}, 
		  {26, 30, 34, 38, 42, 46, 50, 53, 57, 61, 65, 69, 73, 77, 81, 84, 88, 92, 96, 100, 104, 108, 112, 115, 119, 123, 127, 131, 134, -1, -1, -1}}, // d3pu
        	 {{1, 4, 7, 11, 15, 19, 23, 27, 31, 35, 38, 42, 46, 50, 54, 58, 62, 66, 70, 73, 77, 81, 85, 89, 93, 97, 101, 104, 108, -1, -1, -1}, 
		  {26, 30, 34, 38, 42, 46, 49, 53, 57, 61, 65, 69, 73, 77, 80, 84, 88, 92, 96, 100, 104, 108, 112, 115, 119, 123, 127, 131, 134, -1, -1, -1}}}, // d3pv
        	{{{1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 109, 113, -1, -1, -1}, 
		  {4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 112, 116, -1, -1, -1}}, // d3px
        	 {{1, 4, 8, 11, 15, 19, 23, 27, 31, 35, 39, 43, 46, 50, 54, 58, 62, 66, 70, 74, 77, 81, 85, 89, 93, 97, 101, 105, 109, -1, -1, -1}, 
		  {27, 31, 34, 38, 42, 46, 50, 54, 58, 62, 65, 69, 73, 77, 81, 85, 89, 93, 96, 100, 104, 108, 112, 116, 120, 124, 128, 131, 134, -1, -1, -1}}, // d3pu
        	 {{1, 4, 8, 12, 16, 20, 24, 28, 31, 35, 39, 43, 47, 51, 55, 59, 62, 66, 70, 74, 78, 82, 86, 90, 94, 97, 101, 105, 109, -1, -1, -1}, 
		  {27, 31, 34, 38, 42, 46, 50, 54, 58, 62, 65, 69, 73, 77, 81, 85, 89, 93, 97, 100, 104, 108, 112, 116, 120, 124, 128, 131, 134, -1, -1, -1}}} // d3pv
        	};
}

//clone of LoadEvent::Hit:
class gHit {
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

namespace extrapolation_tools{
	__device__ constexpr float straight_st1_det_extrap[2][2][2] = { { {-0.004794, 170.823}, {-0.00207589,  -170.346} }, 
									{ {0.244445, 166.351},  {0.031227, -171.986} } };
	
	__device__ constexpr float invP_x0_[2][2] = {{-0.00422085, 0.00107737}, {0.00157655, 0.000549662}};
	__device__ constexpr float err_invP_x0[2] = {0.0041388368, 0.0061518968};
}


class gTracklet {
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

class gTrack2D {
      public:
      // note: x_ is to be understood as either x or y...
      float tx_;
      float x_0;
      
      float err_tx_;
      float err_x_0;
      
      float chisq;
              
      short nhits;
      int hitlist[MaxHitsPerTrack];
      short hitsign[MaxHitsPerTrack];
      
};

class gEvent {
	public:
	gEvent(){
		RunID = EventID = SpillID = -1;
	}
	int RunID; // Run Number
	int EventID; // Event number
	int SpillID; // Spill number
	int TriggerBits; // hash of the trigger bits: 0-4: MATRIX1-5; 5-9: NIM1-5;
	short TargetPos; // target position: proxy for target ID?
	int TurnID; // => related to beam intensity
	int RFID; // => related to beam intensity
	int Intensity[33]; //  16 before, one onset, and 16 after
	short TriggerEmu; // 1 if MC event
	short NRoads[4]; // 0, positive top; 1, positive bottom; 2, negative top; 3, negative bottom
	int NHits[nDetectors+1]; // number of hits in each detector plane
	int nAH; // size of AllHits
	int nTH; // size of TriggerHits
	gHit AllHits[EstnAHMax]; // array of all hits
	gHit TriggerHits[EstnTHMax]; // array of trigger hits
	bool HasTooManyHits;//bool to flag an event with too many hits
};

class gFullTrackBuilder{
public:
	gTrack2D TrackXZ_st1;
	
      	int hitlist[MaxHitsPerTrack];
      	short hitsign[MaxHitsPerTrack];
	
	thrust::pair<int, int> hitpairs_x1[100];
	thrust::pair<int, int> hitpairs_u1[100];
	thrust::pair<int, int> hitpairs_v1[100];
	
	//util arrays for pair making
	int hitidx1[100];
	int hitidx2[100];
	short hitflag1[100];
	short hitflag2[100];
};


class gStraightTrackBuilder{
public:
	gTrack2D trackXZ;
	gTrack2D trackYZ;
	gTrack2D besttrackYZ;
	
        //pairs in station 2
        thrust::pair<int, int> hitpairs_x2[280];//28*10
        thrust::pair<int, int> hitpairs_u2[1120];//28*40
        thrust::pair<int, int> hitpairs_v2[1120];
        //pairs in station 3
        thrust::pair<int, int> hitpairs_x3[580];//29*10*2
        thrust::pair<int, int> hitpairs_u3[2320];//29*40*2
        thrust::pair<int, int> hitpairs_v3[2320];
	
	int nhitpairs_x2[28];
	int nhitpairs_u2[28];
	int nhitpairs_v2[28];
	
	int nhitpairs_x3[58];
	int nhitpairs_u3[58];
	int nhitpairs_v3[58];
	
	//util arrays for pair making
	int hitidx1[100];
	int hitidx2[100];
	short hitflag1[100];
	short hitflag2[100];
};

class gStraightFitArrays {
public:
      int npoints;
      float drift_dist[MaxHitsPerTrack]; // hit drift distance
      float resolution[MaxHitsPerTrack]; // detector resolution
      
      float p1x[MaxHitsPerTrack];// x bottom end point of the wire hit 
      float p1y[MaxHitsPerTrack];// y bottom end point of the wire hit 
      float p1z[MaxHitsPerTrack];// z bottom end point of the wire hit 
      
      float deltapx[MaxHitsPerTrack];// x distance between bottom and top end points of the wire hit 
      float deltapy[MaxHitsPerTrack];// y distance between bottom and top end points of the wire hit 
      float deltapz[MaxHitsPerTrack];// z distance between bottom and top end points of the wire hit 
      
      float output_parameters[4];
      float output_parameters_errors[4];
      float chi2_2d;
      float chi2;

      float x_array[MaxHitsPerTrack];// x position arrays
      float y_array[MaxHitsPerTrack];// y position arrays
      float z_array[MaxHitsPerTrack];// z position arrays
      float dx_array[MaxHitsPerTrack];// x position uncertainty
      float dy_array[MaxHitsPerTrack];// x position uncertainty
      
      float A[4];// matrix: max size 2x2
      float Ainv[4];// inverted matrix
      float B[2];// input vector
};

class gKalmanFitArrays{
public:
	float state[5];// 5-vector: x0, y0, tx, ty, invP
	float Cov[25];// symmetric 5x5 matrix: C00 = err_x0, C11 = err_y0, C22 = err_tx, C33 = err_ty, C44 = err_invP
	float H[2];
	float K[5];// matrix: max size 5x5, but we can use this unique array for all possible sizes
	float KCResKt[25];// matrix 5x5, result of tensor product of K*K
	float chi2;// chi2
};

class gOutputEvent {
public:
	int EventID;
	int nAH;
	bool HasTooManyHits;//bool to flag an event with too many hits
	int nTracklets;
	gTracklet AllTracklets[TrackletSizeMax+1];//save one slot for a new candidate
};

//geometry carrier
class gPlane {
      public:
      float z;
      int nelem;
      float cellwidth;
      float spacing;
      float xoffset;
      float scalex;
      float x0;
      float x1;
      float x2;
      float costheta;
      float scaley;
      float y0;
      float y1;
      float y2;
      float sintheta;
      float resolution;
      float deltaW_[9];
      float z_mean;
      float u_win;
      float v_win_fac1;
      float v_win_fac2;
      float v_win_fac3;
      float p1x_w1;
      float p1y_w1;
      float p1z_w1;
      float deltapx;
      float deltapy;
      float deltapz;
      float dp1x;
      float dp1y;
      float dp1z;
      float slope_max;
      float inter_max;
};

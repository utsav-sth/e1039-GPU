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
const int TrackletSizeMax = 200;

const double TX_MAX = 0.15;
const double TY_MAX = 0.1;
const double X0_MAX = 150;
const double Y0_MAX = 50;
const double INVP_MAX = 0.2;
const double INVP_MIN = 0.01;

namespace geometry{
	__device__ constexpr float spacingplane[28] = {0., 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 1.3, 1.3, 1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 4.0, 4.0, 7.0, 7.0, 8.0, 12.0, 12.0, 10.0, 3.0, 3.0, 3.0, 3.0};
	__device__ constexpr short detsuperid[7][3] = {{2, 1, 3}, {5, 6, 4}, {8, 9, 7}, {11, 12, 10}, {14, 15, 13}, {25, 26, -1}, {24, 27, -1}}; 
	__device__ constexpr short dets_x[6] = {15, 16, 21, 22, 27, 28};
}

//clone of LoadEvent::Hit:
class gHit {
	public:
	int index; // global hit index in the hit array
	short detectorID; // ID of the detector: one ID for each DC wire plane (30 total), hodoscope plane (16 total), proportional tube plane (8 total).
	short elementID; // ID of the element in the detector: wire/slat/tube number
	float tdcTime; // raw TDC time from the DAQ 
	float driftDistance; // calculated drift distance from RT profile (supplied in database) IF tdcTime between tmin and tmax defined for detector; 
	float pos; // position in the projection of the detector (e.g. X in a X plane, etc)
	short flag; // 1: in time; 2: hodo mask; 3: trigger mask
};

class gTracklet {
      public:
      gTracklet(){
	nXHits = nUHits = nVHits = 0;
      }
      
      short stationID;
      short nXHits;
      short nUHits;
      short nVHits;

      float chisq;
      float chisq_vtx;

      //maybe we can be a bit more sober in memory, and just require hit "event" index?
      gHit hits[nDetectors];// array of all hits
      short hitsign[nDetectors];
      
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
      
      float residual[nChamberPlanes];
};

class gFitParams {
public:
      int max_iterations;
      short nparam;
      
      float parameter_limits_min[5];
      float parameter_limits_max[5];
      
      float tolerance;
};




class gFitArrays {
public:
      int npoints;
      float drift_dist[nChamberPlanes]; // hit drift distance
      float resolution[nChamberPlanes]; // detector resolution
      
      float p1x[nChamberPlanes];// x bottom end point of the wire hit 
      float p1y[nChamberPlanes];// y bottom end point of the wire hit 
      float p1z[nChamberPlanes];// z bottom end point of the wire hit 
      
      float deltapx[nChamberPlanes];// x distance between bottom and top end points of the wire hit 
      float deltapy[nChamberPlanes];// y distance between bottom and top end points of the wire hit 
      float deltapz[nChamberPlanes];// z distance between bottom and top end points of the wire hit 
      
      float prev_parameters[5];
      float output_parameters[5];
      float output_parameters_errors[5];
      float chi2prev;
      float chi2;

      float values[nChamberPlanes];
      float derivatives[5*nChamberPlanes];
      float gradients[5];
      float hessians[25];
      
      float scaling_vector[5];
      float deltas[5];
      
      float lambda;//scale factor
      
      //float calc_matrix[25];
      //float abs_row[25];
      //int abs_row_index[25];

      int n_iter;
      short iter_failed;
      short finished;
      short state;
      short skip;
      short singular;
      
      //float x_array[nChamberPlanes];// x position arrays
      //float y_array[nChamberPlanes];// y position arrays
      //float z_array[nChamberPlanes];// z position arrays
      //float dx_array[nChamberPlanes];// x position uncertainty
      //float dy_array[nChamberPlanes];// x position uncertainty

      float A[25];// matrix: max size 5x5, but we can use this unique array for all possible sizes
      float Ainv[25];// matrix
      float B[5];// input vector

      //float output_parameters_steps[5];
      //float doublederivatives[5];
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
};

//Output class
class gSW {
public:
	gSW(){
		EventID = -1;
		nAH = 0;
		nTracklets = 0;
	}
	thrust::pair<int, int> hitpairs_x[100];
	thrust::pair<int, int> hitpairs_u[100];
	thrust::pair<int, int> hitpairs_v[100];
	int hitidx1[100];
	int hitidx2[100];
	short hitflag1[100];
	short hitflag2[100];
	int EventID;
	int nAH;
	int nTracklets;
	//really tempted to replace tracklet array with array of IDs
	gTracklet AllTracklets[TrackletSizeMax];
	short nTKL_stID[7];//0: D0; 1: D1; 2: D2; 3: D3p; 4: D3m; 5: back partial; 6: global
};

//geometry carrier
class gPlane {
      public:
      float z;
      int nelem;
      float spacing;
      float xoffset;
      float scalex;
      float x0;
      float costheta;
      float scaley;
      float y0;
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
};

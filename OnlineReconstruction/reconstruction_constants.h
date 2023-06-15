//#define GLDISPLAY 1
#define E1039 1
//#define DEBUG 1
#define USE_DET_RESOL 1
#define KMAG_ON 1
//#define SAVE_ALL_TRACKS 1
#define REFINED_ER
#define PROP_Y_MATCH 1
//#define HODO_Y_MATCH 1
#define nChamberPlanes 30
#define nHodoPlanes 16
#define nPropPlanes 8
#define nDetectors (nChamberPlanes+nHodoPlanes+nPropPlanes+1)

#define triggerBit(n) (1 << (n))
#define hitFlagBit(n) (1 << (n))

using namespace std;

const int EstnEvtMax = 10000;
const int THREADS_PER_BLOCK = 32;//16;//8;
// eight threads per block for ER: do 3 chambers, 2 hodoscopes, 1 prop tube per thread! 
// eight threads per block for XZ tracking, YZ tracking: do 7 bins in st2 * 29 bins in st3!
int BLOCKS_NUM = EstnEvtMax;///THREADS_PER_BLOCK;
const int EstnAHMax = 4096;
const int EstnTHMax = 256;
const int ClusterSizeMax = 192;

const int NVars = 11;
const int Nbins_Hists = 128;
const float TX_MAX = 0.15;
const float TY_MAX = 0.1;
const float X0_MAX = 150;//cm
const float Y0_MAX = 50;//cm
const float INVP_MAX = 0.2;//GeV-1
const float INVP_MIN = 0.01;//GeV-1
const float VXY_MAX = 20.0;//cm
const float VZ_MIN = -500.0;//cm
const float VZ_MAX = +300.0;//cm
const float PXY_MAX = 20.0;//GeV
const float PZ_MIN = 0.0;//GeV
const float PZ_MAX = 150.0;//GeV

const float InvSqrt12 = 0.288675135;

#ifdef EXTRASTUFF
__device__ constexpr float gauss_quantiles[101] = {-3.003f, -2.32567f, -2.05594f, -1.88212f, -1.75025f, -1.64236f, -1.55245f, -1.47453f, -1.4026f, -1.34266f, -1.28272f, -1.22877f, -1.17483f, -1.12687f, -1.07892f, -1.03696f, -0.995005f, -0.953047f, -0.917083f, -0.881119f, -0.839161f, -0.809191f, -0.773227f, -0.737263f, -0.707293f, -0.677323f, -0.641359f, -0.611389f, -0.581419f, -0.551449f, -0.527473f, -0.497502f, -0.467532f, -0.437562f, -0.413586f, -0.383616f, -0.35964f, -0.32967f, -0.305694f, -0.281718f, -0.251748f, -0.227772f, -0.203796f, -0.173826f, -0.14985f, -0.125874f, -0.101898f, -0.0719281f, -0.047952f, -0.023976f, 0.0f, 0.023976f, 0.047952f, 0.0779221f, 0.101898f, 0.125874f, 0.14985f, 0.173826f, 0.203796f, 0.227772f, 0.251748f, 0.281718f, 0.305694f, 0.32967f, 0.35964f, 0.383616f, 0.413586f, 0.437562f, 0.467532f, 0.497502f, 0.527473f, 0.551449f, 0.581419f, 0.611389f, 0.641359f, 0.677323f, 0.707293f, 0.737263f, 0.773227f, 0.809191f, 0.839161f, 0.881119f, 0.917083f, 0.953047f, 0.995005f, 1.03696f, 1.07892f, 1.12687f, 1.17483f, 1.22877f, 1.28272f, 1.34266f, 1.4026f, 1.47453f, 1.55245f, 1.64236f, 1.75025f, 1.88212f, 2.04995f, 2.32567f, 3.003};

const float deltax0_ = 0.221572;
const float deltax0_sigma = 0.472572;
const float deltay0_ = 0.0387412;
const float deltay0_sigma = 0.770639;
#endif

namespace geometry{
#ifdef E1039
	__device__ constexpr short eff_detid_chambers[24] = { 1,  2,  3,  4,  5,  6, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
#else
	__device__ constexpr short eff_detid_chambers[24] = { 7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
#endif
	__device__ constexpr float spacingplane[28] = {0., 0.40, 0.40, 0.40, 0.40, 0.40, 0.40, 1.3, 1.3, 1.3, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 4.0, 4.0, 7.0, 7.0, 8.0, 12.0, 12.0, 10.0, 3.0, 3.0, 3.0, 3.0};
	__device__ constexpr short detsuperid[7][3] = {{2, 1, 3}, {5, 6, 4}, {8, 9, 7}, {11, 12, 10}, {14, 15, 13}, {25, 24, -1}, {26, 27, -1}}; 
	__device__ constexpr short planetype[15] = {1, 0, 2, 1, 0, 2, 2, 0, 1, 2, 0, 1, 2, 0, 1};
	//__device__ constexpr short sagitta_detid[3] = {1};
	//__device__ constexpr short detmap_zincr[31] = {-1,  0,  1,  3,  2,  4,  5,  0,  1,  3,  2,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 12, 13, 14, 15, 16, 17};
	//__device__ constexpr short dets_x[6] = {15, 16, 21, 22, 27, 28};
	//__device__ constexpr short hodoplanerange[5][2] = {{31, 34}, {31, 34}, {35, 38}, {39, 40}, {41, 46}};// range of planes to look for hits
	__device__ constexpr short hodoplanesx[4][2] = {{31, 32}, {37, 38}, {39, 40}, {45, 46}};// range of planes to look for hits
	__device__ constexpr short hodoplanesy[4][2] = {{33, 34}, {35, 36}, {41, 42}, {43, 44}};// range of planes to look for hits
	__device__ constexpr float hodofudgefac[4] = {0.25, 0.2, 0.15, 0.0};
	__device__ constexpr float hodoXcellwidth[4] = {7.32, 13.0f, 14.5875, 19.6475};
	__device__ constexpr float X_BEAM = 0.;
	__device__ constexpr float Y_BEAM = 0.;
	__device__ constexpr float Z_DUMP = 42.;
	__device__ constexpr float SAGITTA_TARGET_CENTER = 1.85;
	__device__ constexpr float SAGITTA_TARGET_WIDTH = 0.25;
	__device__ constexpr float SAGITTA_DUMP_CENTER = 1.5;
	__device__ constexpr float SAGITTA_DUMP_WIDTH = 0.3;
//	__device__ constexpr float X_KMAG_BEND = 144.78;
//	__device__ constexpr float Y_KMAG_BEND = 101.6;
#ifdef E1039
	__device__ constexpr float PT_KICK_KMAG = -0.3819216;//PT_KICK_MAG*KMAGSTR = 0.4016* -0.951;
#ifdef KMAG_ON
	__device__ constexpr float KMAGSTR = -0.951;
#else
	__device__ constexpr float KMAGSTR = 0.0;
#endif
	__device__ constexpr float FMAGSTR = -1.054;
	//__device__ constexpr float PT_KICK_FMAG = -3.066086;// PT_KICK_FMAG*FMAGSTR/FMAG_LENGTH = 2.909*-1.054 = -3.066086
	__device__ constexpr float PTKICK_UNIT = -0.006096568;// PT_KICK_FMAG*FMAGSTR/FMAG_LENGTH = 2.909*-1.054/502.92 = -0.006096568
	__device__ constexpr float Z_KMAG_BEND = 1064.26;
	__device__ constexpr float Z_FMAG_BEND = 251.4;
	__device__ constexpr float Z_TARGET = -300.;
#else
	__device__ constexpr float PT_KICK_KMAG = -0.4141;//PT_KICK_MAG*KMAGSTR = 0.404* -1.025;
#ifdef KMAG_ON
	__device__ constexpr float KMAGSTR = -1.025;
#else
	__device__ constexpr float KMAGSTR = 0.0;
#endif
	__device__ constexpr float FMAGSTR = -1.044;
	__device__ constexpr float PTKICK_UNIT = -0.006038726;//-0.PT_KICK_FMAG*FMAGSTR/FMAG_LENGTH = 2.909*-1.044/502.92 = -0.006038726
	__device__ constexpr float Z_KMAG_BEND = 1041.8;
	__device__ constexpr float Z_FMAG_BEND = 251.4;
	__device__ constexpr float Z_TARGET = -129.54;
#endif
	__device__ constexpr float FMAG_LENGTH = 502.92;
	__device__ constexpr float FMAG_HOLE_LENGTH = 27.94;
	__device__ constexpr float FMAG_HOLE_RADIUS = 1.27;
	__device__ constexpr int NSTEPS_FMAG = 100;
	__device__ constexpr float STEP_FMAG = 2.5146;// FMAG_LENGTH/NSTEPS_FMAG/2 = 502.92/100/2 = 2.5146 
	__device__ constexpr float DEDX_UNIT_0 = 0.014282073;// 7.18274/502.92;
	__device__ constexpr float DEDX_UNIT_1 = 7.1869681e-05;// 0.0361447/502.92;
	__device__ constexpr float DEDX_UNIT_2 =-1.4279150e-06;//-0.000718127/502.92;
	__device__ constexpr float DEDX_UNIT_3 = 1.5853655e-08;// 7.97312e-06/502.92;
	__device__ constexpr float DEDX_UNIT_4 =-6.0741470e-11;//-3.05481e-08/502.92;
	__device__ constexpr int NSTEPS_TARGET = 200;
	__device__ constexpr float STEP_TARGET = 2.5;// |Z_UPSTREAM|/NSTEPS_TARGET = 500/100 = 5. 
	__device__ constexpr short lrpossibility[4][2] = {{1, 1}, {1, -1}, {-1, 1}, {-1, -1}};
	__device__ constexpr short N_WCHitsBins[4] = {32, 28, 28, 28};
	__device__ constexpr short MaxHitsProj[3] = {20, 80, 80};
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
        	{{{1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 109, -1, -1, -1, -1}, 
		  {4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 116, -1, -1, -1, -1}}, // d3px
        	 {{1, 4, 7, 11, 15, 19, 23, 27, 31, 35, 38, 42, 46, 50, 54, 58, 62, 66, 70, 73, 77, 81, 85, 89, 93, 97, 101, 104, -1, -1, -1, -1}, 
		  {26, 30, 34, 38, 42, 46, 50, 53, 57, 61, 65, 69, 73, 77, 81, 84, 88, 92, 96, 100, 104, 108, 112, 115, 119, 123, 127, 134, -1, -1, -1, -1}}, // d3pu
        	 {{1, 4, 7, 11, 15, 19, 23, 27, 31, 35, 38, 42, 46, 50, 54, 58, 62, 66, 70, 73, 77, 81, 85, 89, 93, 97, 101, 104, -1, -1, -1, -1}, 
		  {26, 30, 34, 38, 42, 46, 49, 53, 57, 61, 65, 69, 73, 77, 80, 84, 88, 92, 96, 100, 104, 108, 112, 115, 119, 123, 127, 134, -1, -1, -1, -1}}}, // d3pv
        	{{{1, 5, 9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61, 65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 109, -1, -1, -1, -1}, 
		  {4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 76, 80, 84, 88, 92, 96, 100, 104, 108, 116, -1, -1, -1, -1}}, // d3px
        	 {{1, 4, 8, 11, 15, 19, 23, 27, 31, 35, 39, 43, 46, 50, 54, 58, 62, 66, 70, 74, 77, 81, 85, 89, 93, 97, 101, 105, -1, -1, -1, -1}, 
		  {27, 31, 34, 38, 42, 46, 50, 54, 58, 62, 65, 69, 73, 77, 81, 85, 89, 93, 96, 100, 104, 108, 112, 116, 120, 124, 128, 134, -1, -1, -1, -1}}, // d3pu
        	 {{1, 4, 8, 12, 16, 20, 24, 28, 31, 35, 39, 43, 47, 51, 55, 59, 62, 66, 70, 74, 78, 82, 86, 90, 94, 97, 101, 105, -1, -1, -1, -1}, 
		  {27, 31, 34, 38, 42, 46, 50, 54, 58, 62, 65, 69, 73, 77, 81, 85, 89, 93, 97, 100, 104, 108, 112, 116, 120, 124, 128, 134, -1, -1, -1, -1}}} // d3pv
        	};//TODO: reduce this monster to just what is needed/used
}

namespace datasizes{
	__host__ __device__ constexpr int NHitsParam = 5;
#ifdef TESTMOMENTUM
	__host__ __device__ constexpr int NTracksParam = 150;
#else
	__host__ __device__ constexpr int NTracksParam = 148;
#endif
	__host__ __device__ constexpr int NMaxHitsChambers = 200;
	__host__ __device__ constexpr int NMaxHitsHodoscopes = 60;
	__host__ __device__ constexpr int NMaxHitsPropTubes = 120;
	__host__ __device__ constexpr unsigned int eventhitsize[3] = {
		nChamberPlanes*datasizes::NHitsParam*datasizes::NMaxHitsChambers, 
		nHodoPlanes*datasizes::NHitsParam*datasizes::NMaxHitsHodoscopes, 
		nPropPlanes*datasizes::NHitsParam*datasizes::NMaxHitsPropTubes
	};
	__host__ __device__ constexpr int TrackSizeMax = 960;
	__host__ __device__ constexpr int MaxHitsPerTrack = 18;
}

namespace selection{
	__host__ __device__ constexpr short MaxD0Multiplicity = 140;
	__host__ __device__ constexpr short MaxD2Multiplicity = 80;
	__host__ __device__ constexpr short MaxD3Multiplicity = 60;
	__host__ __device__ constexpr short MaxPropMultiplicity = 1250;
	
	__host__ __device__ constexpr float chi2dofmax = 250;
	__host__ __device__ constexpr short NpropXYhitsMin = 2;
	__host__ __device__ constexpr float merge_thres = 0.015;
	__host__ __device__ constexpr float rejectwin[4] = {1.5, 2.0, 1.8, 1.8};
	__host__ __device__ constexpr float x_vtx_cut = 15.0;
	__host__ __device__ constexpr float y_vtx_cut = 25.0;
	__host__ __device__ constexpr float tx_vtx_cut = 0.15;
	__host__ __device__ constexpr float ty_vtx_cut = 0.1;
#ifdef E1039
	//__device__ constexpr float rejectwin[4] = {0.12, 0.15, 0.16, 0.14};
#else
	//__device__ constexpr float rejectwin[4] = {0.13, 0.18, 0.16, 0.16};
#endif
}

namespace extrapolation_tools{
	__device__ constexpr float straight_st1_det_extrap[2][2][2] = { { {-0.004794, 170.823}, {-0.00207589,  -170.346} }, 
									{ {0.244445, 166.351},  {0.031227, -171.986} } };
	
	__device__ constexpr float invP_x0_[2][2] = {{-0.00422085, 0.00107737}, {0.00157655, 0.000549662}};
	__device__ constexpr float err_invP_x0[2] = {0.0041388368, 0.0061518968};
}

namespace debug{
  __host__ __device__ constexpr unsigned int EvRef = 13;
}

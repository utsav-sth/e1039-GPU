#include "reconstruction_constants.h"

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


struct gTracklet {
	public:
	float* m_trackletdata;
	
	__host__ __device__ gTracklet(float* basedata, const unsigned offset = 0) :
		m_trackletdata(reinterpret_cast<float*>(basedata) + offset)
		{
			static_assert(sizeof(float) == sizeof(unsigned));
			assert((((size_t) basedata) & sizeof(float)) == 0);
		}
	
	__host__ __device__ inline float stationID() const
		{
			return m_trackletdata[0];
		}

	__host__ __device__ inline float threadID() const
		{
			return m_trackletdata[1];
		}

	__host__ __device__ inline float nHits() const
		{
			return m_trackletdata[2];
		}
	
	__host__ __device__ inline float chisq() const
		{
			return m_trackletdata[3];
		}
	
	__host__ __device__ inline float chisq_vtx() const
		{
			return m_trackletdata[4];
		}

	//track parameters
	__host__ __device__ inline float tx() const
		{
			return m_trackletdata[5];
		}

	__host__ __device__ inline float ty() const
		{
			return m_trackletdata[6];
		}

	__host__ __device__ inline float x0() const
		{
			return m_trackletdata[7];
		}

	__host__ __device__ inline float y0() const
		{
			return m_trackletdata[8];
		}

	__host__ __device__ inline float invP() const
		{
			return m_trackletdata[9];
		}

	__host__ __device__ inline float err_tx() const
		{
			return m_trackletdata[10];
		}

	__host__ __device__ inline float err_ty() const
		{
			return m_trackletdata[11];
		}

	__host__ __device__ inline float err_x0() const
		{
			return m_trackletdata[12];
		}

	__host__ __device__ inline float err_y0() const
		{
			return m_trackletdata[13];
		}

	__host__ __device__ inline float err_invP() const
		{
			return m_trackletdata[14];
		}

	__host__ __device__ inline float charge() const
		{
			return m_trackletdata[15];
		}

	__host__ __device__ inline float hits_detid(const unsigned ihit) const
		{
			return m_trackletdata[16 + ihit ];
		}

	__host__ __device__ inline float hits_chan(const unsigned ihit) const
		{
			return m_trackletdata[34 + ihit ];
		}

	__host__ __device__ inline float hits_pos(const unsigned ihit) const
		{
			return m_trackletdata[52 + ihit ];
		}

	__host__ __device__ inline float hits_drift(const unsigned ihit) const
		{
			return m_trackletdata[70 + ihit ];
		}
	
	__host__ __device__ inline float hits_sign(const unsigned ihit) const
		{
			return m_trackletdata[88 + ihit ];
		}
	__host__ __device__ inline float hits_tdc(const unsigned ihit) const
		{
			return m_trackletdata[124 + ihit ];
		}
	__host__ __device__ inline float hits_residual(const unsigned ihit) const
		{
			return m_trackletdata[124 + ihit ];
		}
	__host__ __device__ inline unsigned int get_lasthitdetid() const
		{
			const int nhits = (int)nHits();
			int detid;
			int detid_max = -1;
			for(int i = 0; i<nhits; i++){
				detid = (int)hits_detid(i);
				if(detid>detid_max)detid_max = detid;
			}
			return detid_max;		 
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

	//vertex parameters...
	__host__ __device__ inline float vx(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+16];
		}
	__host__ __device__ inline float vy(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+17];
		}
	__host__ __device__ inline float vz(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+18];
		}
	__host__ __device__ inline float px(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+19];
		}
	__host__ __device__ inline float py(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+20];
		}
	__host__ __device__ inline float pz(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+21];
		}
	
	//then hits...
	__host__ __device__ inline float hits_detid(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+22 + ihit ];
		}

	__host__ __device__ inline float hits_chan(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+40 + ihit ];
		}

	__host__ __device__ inline float hits_pos(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+58 + ihit ];
		}

	__host__ __device__ inline float hits_drift(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+76 + ihit ];
		}
	
	__host__ __device__ inline float hits_sign(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+94 + ihit ];
		}
//additional hit parameters
	__host__ __device__ inline float hits_tdc(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+112 + ihit ];
		}
	__host__ __device__ inline float hits_residual(const unsigned index, const unsigned ihit) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+130 + ihit ];
		}
#ifdef ALLVTXVARS
#ifdef TESTMOMENTUM
	__host__ __device__ inline float invp_tgt(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+148];
		}
#endif
//additional vertex parameters
	__host__ __device__ inline float tgt_vx(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+148];
		}
	__host__ __device__ inline float tgt_vy(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+149];
		}
	__host__ __device__ inline float tgt_vz(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+150];
		}
	__host__ __device__ inline float tgt_px(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+151];
		}
	__host__ __device__ inline float tgt_py(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+152];
		}
	__host__ __device__ inline float tgt_pz(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+153];
		}	
	__host__ __device__ inline float dump_vx(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+154];
		}
	__host__ __device__ inline float dump_vy(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+155];
		}
	__host__ __device__ inline float dump_vz(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+156];
		}
	__host__ __device__ inline float dump_px(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+157];
		}
	__host__ __device__ inline float dump_py(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+158];
		}
	__host__ __device__ inline float dump_pz(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return m_trackdata[TrackSize*index+159];
		}	
#endif
	

	__host__ __device__ inline unsigned int get_lasthitdetid(const unsigned index, const float* z_array) const
		{
			assert(index < NTracksTotal);
			const int nhits = (int)nHits(index);
			short detid;
			short detid_max = 1;
			float z_max = -10000.f;
			for(short i = 0; i<nhits; i++){
				detid = (short)hits_detid(index, i);
				if(detid>0 && z_array[detid]>z_max){
					detid_max = detid;
					z_max = z_array[detid];
				}
			}
			return detid_max;		 
		}

	__host__ __device__ inline unsigned int get_firsthitdetid(const unsigned index, const float* z_array) const
		{
			assert(index < NTracksTotal);
			const int nhits = (int)nHits(index);
			short detid;
			short detid_min = 30;
			float z_min = +10000.f;
			for(short i = 0; i<nhits; i++){
				detid = (short)hits_detid(index, i);
				if(detid>0 && z_array[detid]<z_min){
					detid_min = detid;
					z_min = z_array[detid];
				}
			}
			return detid_min;
		}
	
	__host__ __device__ inline gTracklet Track(const unsigned index) const
		{
			assert(index < NTracksTotal);
			return gTracklet(m_trackdata, TrackSize*index);
		}
};


struct gEventTrackCollection{
	unsigned short NTracks[EstnEvtMax*THREADS_PER_BLOCK];
	float TracksRawData[EstnEvtMax*datasizes::TrackSizeMax*datasizes::NTracksParam];
	//unsigned short NRecTracks[EstnEvtMax*THREADS_PER_BLOCK];
	//float RecTracksRawData[EstnEvtMax*datasizes::TrackSizeMax*datasizes::NRecTracksParam];
	//gTracklet Tracklets[EstnEvtMax*datasizes::TrackSizeMax];
	//__device__ const gTracklets tracks(const unsigned int event, int& ntracks) {
	//	ntracks = NTracks[event];
	//	return gTracklets(Tracklets, ntracks, event*datasizes::TrackSizeMax);
	//}

	__device__ const gTracks tracks(const unsigned int event, const unsigned threadID, unsigned int& ntracks) {
		ntracks = NTracks[event*THREADS_PER_BLOCK+threadID];
		return gTracks(TracksRawData, ntracks, event*datasizes::TrackSizeMax*datasizes::NTracksParam+threadID*datasizes::TrackSizeMax*datasizes::NTracksParam/THREADS_PER_BLOCK);
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
	__device__ void setCharge(const unsigned int evt_offset, const unsigned int itrack, const float charge) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+15] = charge;
	}
	//vertex info
	__device__ void setVtxPos(const unsigned int evt_offset, const unsigned int itrack, const float* pos) {
		for(short i = 0; i<3; i++)TracksRawData[evt_offset+itrack*datasizes::NTracksParam+16+i] = pos[i];
	}
	__device__ void setVtxMom(const unsigned int evt_offset, const unsigned int itrack, const float* mom) {
		for(short i = 0; i<3; i++)TracksRawData[evt_offset+itrack*datasizes::NTracksParam+19+i] = mom[i];
	}
	//hit info
	__device__ void setHitDetID(const unsigned int evt_offset, const unsigned int itrack, const unsigned int ihit, const float detid) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+22+ihit] = detid;
	}
	__device__ void setHitChan(const unsigned int evt_offset, const unsigned int itrack, const unsigned int ihit, const float chan) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+40+ihit] = chan;
	}
	__device__ void setHitPos(const unsigned int evt_offset, const unsigned int itrack, const unsigned int ihit, const float pos) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+58+ihit] = pos;
	}
	__device__ void setHitDrift(const unsigned int evt_offset, const unsigned int itrack, const unsigned int ihit, const float drift) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+76+ihit] = drift;
	}
	__device__ void setHitSign(const unsigned int evt_offset, const unsigned int itrack, const unsigned int ihit, const float sign) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+94+ihit] = sign;
	}
	__device__ void setHitTDC(const unsigned int evt_offset, const unsigned int itrack, const unsigned int ihit, const float tdc) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+112+ihit] = tdc;
	}
	__device__ void setHitResidual(const unsigned int evt_offset, const unsigned int itrack, const unsigned int ihit, const float resid) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+130+ihit] = resid;
	}
#ifdef ALLVTXVARS
#ifdef TEST_MOMENTUM
	__device__ void setInvPTarget(const unsigned int evt_offset, const unsigned int itrack, const float invp_tgt) {
		TracksRawData[evt_offset+itrack*datasizes::NTracksParam+148] = invp_tgt;
	}
#endif
	__device__ void setTargetPos(const unsigned int evt_offset, const unsigned int itrack, const float* pos) {
		for(short i = 0; i<3; i++)TracksRawData[evt_offset+itrack*datasizes::NTracksParam+148+i] = pos[i];
	}
	__device__ void setTargetMom(const unsigned int evt_offset, const unsigned int itrack, const float* mom) {
		for(short i = 0; i<3; i++)TracksRawData[evt_offset+itrack*datasizes::NTracksParam+151+i] = mom[i];
	}
	__device__ void setDumpPos(const unsigned int evt_offset, const unsigned int itrack, const float* pos) {
		for(short i = 0; i<3; i++)TracksRawData[evt_offset+itrack*datasizes::NTracksParam+154+i] = pos[i];
	}
	__device__ void setDumpMom(const unsigned int evt_offset, const unsigned int itrack, const float* mom) {
		for(short i = 0; i<3; i++)TracksRawData[evt_offset+itrack*datasizes::NTracksParam+157+i] = mom[i];
	}
#endif
};


struct gDimuons{
	public:
	const unsigned int NDimuonsTotal;
	const unsigned int DimuonSize;
	float* m_dimuondata;
	
	__host__ __device__ gDimuons(float* basedata, const unsigned total_number_of_dimuons, const unsigned offset = 0) :
		m_dimuondata(reinterpret_cast<float*>(basedata) + offset), NDimuonsTotal(total_number_of_dimuons), DimuonSize(datasizes::NDimuonParam)
		{
			static_assert(sizeof(float) == sizeof(unsigned));
			assert((((size_t) basedata) & sizeof(float)) == 0);
		}
		
	__host__ __device__ inline float mass(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index ];
		}
	__host__ __device__ inline float pT(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+1 ];
		}
	__host__ __device__ inline float xF(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+2 ];
		}
	__host__ __device__ inline float x1(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+3 ];
		}
	__host__ __device__ inline float x2(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+4 ];
		}
	__host__ __device__ inline float costheta(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+5 ];
		}
	__host__ __device__ inline float phi(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+6 ];
		}
	__host__ __device__ inline float mass_single(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+7 ];
		}
	__host__ __device__ inline float chi2_single(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+8 ];
		}	
	
	__host__ __device__ inline float dim_vx(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+9 ];
		}
	__host__ __device__ inline float dim_vy(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+10 ];
		}
	__host__ __device__ inline float dim_vz(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+11 ];
		}

	__host__ __device__ inline float pos_vx(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+12 ];
		}
	__host__ __device__ inline float pos_vy(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+13 ];
		}
	__host__ __device__ inline float pos_vz(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+14 ];
		}
	__host__ __device__ inline float neg_vx(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+15 ];
		}
	__host__ __device__ inline float neg_vy(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+16 ];
		}
	__host__ __device__ inline float neg_vz(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+17 ];
		}
	
	__host__ __device__ inline float p_pos_x(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+18 ];
		}
	__host__ __device__ inline float p_pos_y(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+19 ];
		}
	__host__ __device__ inline float p_pos_z(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+20 ];
		}
	__host__ __device__ inline float p_pos_E(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+21 ];
		}
	
	__host__ __device__ inline float p_neg_x(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+22 ];
		}
	__host__ __device__ inline float p_neg_y(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+23 ];
		}
	__host__ __device__ inline float p_neg_z(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+24 ];
		}
	__host__ __device__ inline float p_neg_E(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+25 ];
		}

	__host__ __device__ inline float p_pos_single_x(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+26 ];
		}
	__host__ __device__ inline float p_pos_single_y(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+27 ];
		}
	__host__ __device__ inline float p_pos_single_z(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+28 ];
		}
	__host__ __device__ inline float p_pos_single_E(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+29 ];
		}
	
	__host__ __device__ inline float p_neg_single_x(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+30 ];
		}
	__host__ __device__ inline float p_neg_single_y(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+31 ];
		}
	__host__ __device__ inline float p_neg_single_z(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+32 ];
		}
	__host__ __device__ inline float p_neg_single_E(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+33 ];
		}
	
	__host__ __device__ inline float chi2_vtx(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+34 ];
		}
	__host__ __device__ inline float chi2_kf(const unsigned index) const
		{
			assert(index < NDimuonsTotal);
			return m_dimuondata[DimuonSize*index+35 ];
		}
//    //Track projections at different location
//    TVector3 proj_target_pos;
//    TVector3 proj_dump_pos;
//    TVector3 proj_target_neg;
//    TVector3 proj_dump_neg;

//    //Chisq of three test position
//    Double_t chisq_target;
//    Double_t chisq_dump;
//    Double_t chisq_upstream;
};

struct gEventDimuonCollection{
	unsigned short NDimuons[EstnEvtMax];
	float DimuonsRawData[EstnEvtMax*datasizes::DimuonSizeMax*datasizes::NDimuonParam];

	__device__ const gDimuons Dimuons(const unsigned int event, unsigned int& ndimuons) {
		ndimuons = NDimuons[event];
		return gDimuons(DimuonsRawData, ndimuons, event*datasizes::DimuonSizeMax*datasizes::NDimuonParam);
	}
	
	__device__ void setMass(const unsigned int evt_offset, const unsigned int idim, const float mass) {
		DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam] = mass;
	}
	__device__ void setpT(const unsigned int evt_offset, const unsigned int idim, const float pT) {
		DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+1] = pT;
	}
	__device__ void setxF(const unsigned int evt_offset, const unsigned int idim, const float xF) {
		DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+2] = xF;
	}
	__device__ void setx1(const unsigned int evt_offset, const unsigned int idim, const float x1) {
		DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+3] = x1;
	}
	__device__ void setx2(const unsigned int evt_offset, const unsigned int idim, const float x2) {
		DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+4] = x2;
	}
	__device__ void setcostheta(const unsigned int evt_offset, const unsigned int idim, const float costheta) {
		DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+5] = costheta;
	}
	__device__ void setphi(const unsigned int evt_offset, const unsigned int idim, const float phi) {
		DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+6] = phi;
	}
	__device__ void setMassSingle(const unsigned int evt_offset, const unsigned int idim, const float mass_single) {
		DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+7] = mass_single;
	}
	__device__ void setChi2Single(const unsigned int evt_offset, const unsigned int idim, const float chi2_single) {
		DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+8] = chi2_single;
	}
	
	__device__ void setDimVtx(const unsigned int evt_offset, const unsigned int idim, const float* pos) {
		for(short i = 0; i<3; i++)DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+9+i] = pos[i];
	}
	__device__ void setPosVtx(const unsigned int evt_offset, const unsigned int idim, const float* pos) {
		for(short i = 0; i<3; i++)DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+12+i] = pos[i];
	}
	__device__ void setNegVtx(const unsigned int evt_offset, const unsigned int idim, const float* pos) {
		for(short i = 0; i<3; i++)DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+15+i] = pos[i];
	}
	
	__device__ void setPos4Vector(const unsigned int evt_offset, const unsigned int idim, const float* fourvec) {
		for(short i = 0; i<4; i++)DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+18+i] = fourvec[i];
	}
	__device__ void setNeg4Vector(const unsigned int evt_offset, const unsigned int idim, const float* fourvec) {
		for(short i = 0; i<4; i++)DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+22+i] = fourvec[i];
	}
	__device__ void setPosSingle4Vector(const unsigned int evt_offset, const unsigned int idim, const float* fourvec) {
		for(short i = 0; i<4; i++)DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+26+i] = fourvec[i];
	}
	__device__ void setNegSingle4Vector(const unsigned int evt_offset, const unsigned int idim, const float* fourvec) {
		for(short i = 0; i<4; i++)DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+30+i] = fourvec[i];
	}
	
	__device__ void setChi2Vtx(const unsigned int evt_offset, const unsigned int idim, const float chi2_vtx) {
		DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+34] = chi2_vtx;
	}
	__device__ void setChi2KF(const unsigned int evt_offset, const unsigned int idim, const float chi2_KF) {
		DimuonsRawData[evt_offset+idim*datasizes::NDimuonParam+35] = chi2_KF;
	}
};

struct gHistsArrays{
	public:
	short nbins[histotools::ntriggers*histotools::nvars_total];
	float pts_hw[histotools::ntriggers*histotools::nvars_total];
	float xpts[histotools::ntriggers*histotools::nvars_total*histotools::nbins_max];
	float values[histotools::ntriggers*histotools::nvars_total*histotools::nbins_max];
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
      //float deltaW_[nDetectors*9];
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


/*
//Histograms
struct gHist1D {
	public:
	const unsigned int m_nbins;
	const float m_binhw;
	float* m_bincenter;
	float* m_bincontent;
	
	__host__ __device__ gHist1D(const int nbins, const float xmin, const float xmax) :
    		m_nbins(nbins),  m_binhw( (xmax-xmin)/nbins )
    		{
			for(int i = 0; i<m_nbins; i++){
				m_bincenter[i] = xmin+m_binhw*(i+0.5f);
				m_bincontent[i] = 0;
			}
		}

	__host__ __device__ gHist1D(const int nbins, const float binhw, float* xpts, float* values) :
		m_nbins(nbins),  m_binhw(binhw), m_bincenter(reinterpret_cast<float*>(xpts)), m_bincontent(reinterpret_cast<float*>(values)){}
	
	__device__ void Fill(const float x)
		{
			for(int i = 0; i<m_nbins; i++){
				if(m_bincenter[i]-m_binhw <= x && x < m_bincenter[i]+m_binhw);
				m_bincontent[i]+= 1.f;
			}
		}
};
*/


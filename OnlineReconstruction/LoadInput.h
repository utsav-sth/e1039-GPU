#include <iostream>
#include <vector>
#include <TObject.h>
#include <TROOT.h>

#define nChamberPlanes 30
#define nHodoPlanes 16
#define nPropPlanes 8

class Hit: public TObject
{
public:
   Hit();

    Int_t index;
    Short_t detectorID;
    Short_t elementID;
    Float_t tdcTime;
    Float_t driftDistance;
    Float_t pos;
    UShort_t flag;
    ClassDef(Hit, 1)
};

class SRawEvent: public TObject
{
public:
    SRawEvent();
    ~SRawEvent();

    Int_t fRunID;
    Int_t fEventID;
    Int_t fSpillID;

    Int_t fTriggerBits;

    Short_t fTargetPos;

    Int_t fTurnID;
    Int_t fRFID;
    Int_t fIntensity[33];

    Short_t fTriggerEmu;
    Short_t fNRoads[4];

    Int_t fNHits[nChamberPlanes+nHodoPlanes+nPropPlanes+1];
    std::vector<Hit> fAllHits;
    std::vector<Hit> fTriggerHits;
    ClassDef(SRawEvent, 1)
}; 

class SQHit: public TObject
{
public:
  SQHit();
  ~SQHit() {}
  //SQHit*        Clone() const {return (new SQHit_v1(*this));}
  
  int _hit_id;                   ///< hitID
  short _detector_id;            ///< mapping from detector name to ID
  short _element_id;             ///< elementID
  short _level;                  ///< level of v1495 (meaningful only for trigger hit)

  float _tdc_time;               ///< tdcTime
  float _drift_distance;         ///< driftDistance
  float _pos;                    ///< pos?

  unsigned short _flag;          ///< bits collection  
  ClassDef(SQHit, 1)
};

class SQHitVector: public TObject
{
public:
  SQHitVector();
  ~SQHitVector() {}

  //I will add the rest if I feel like it is needed
  //SQHitVector(const SQHitVector& hitmap);


  //SQHitVector& operator=(const SQHitVector& hitmap);
  
  //const SQHit* at(const size_t idkey) const;
  //SQHit*       at(const size_t idkey);
  
  //std::vector<SQHit*>::const_iterator begin()     const {return _vector.begin();}
  //std::vector<SQHit*>::const_iterator   end()     const {return _vector.end();}
  
  std::vector<SQHit*> _vector;
};

class SQEvent: public TObject
{
public:
  SQEvent();
  ~SQEvent() {}
  
  	int _run_id;
	int _spill_id;
	int _event_id;
	int _coda_event_id;

	unsigned short _trigger; //< NIM[1-5], MATRIX[1-5]

	int _raw_matrix[5];

	int _after_inh_matrix[5];

	int _data_quality;

	int _vme_time;

        int _qie_presums[4];
        int _qie_trig_cnt;
        int _qie_turn_id;
        int _qie_rf_id;
        int _qie_rf_inte[33];

        short _flag_v1495;
        short _n_board_qie;
        short _n_board_v1495;
        short _n_board_taiwan;
        short _n_board_trig_b;
        short _n_board_trig_c;

    ClassDef(SQEvent, 1)
};

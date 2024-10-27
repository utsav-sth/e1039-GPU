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


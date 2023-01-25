#include <iostream>
#include <vector>
#include <TObject.h>
#include <TROOT.h>
#include "LoadInput.h"

ClassImp(Hit)
ClassImp(SRawEvent)

Hit::Hit() : index(-1), detectorID(-1), flag(0)
{
}

SRawEvent::SRawEvent() : fRunID(-1), fEventID(-1), fSpillID(-1), fTriggerBits(-1), fTriggerEmu(-1)
{
    for(Int_t i = 0; i < nChamberPlanes+nHodoPlanes+nPropPlanes+1; i++)
    {
        fNHits[i] = 0;
    }
}

SRawEvent::~SRawEvent()
{
}

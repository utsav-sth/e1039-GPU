#include <iostream>
#include <vector>
#include <TObject.h>
#include <TROOT.h>
#include "LoadInput.h"

ClassImp(Hit)
ClassImp(SRawEvent)
ClassImp(SQHit)
ClassImp(SQEvent)

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

SQHit::SQHit() : _hit_id(-1), _detector_id(-1), _flag(0)
{
  
}

SQEvent::SQEvent() : _run_id(-1), _event_id(-1), _spill_id(-1), _coda_event_id(-1), _trigger(-1)
{
    // for(Int_t i = 0; i < nChamberPlanes+nHodoPlanes+nPropPlanes+1; i++)
    // {
    //     fNHits[i] = 0;
    // }
}

SQHitVector::SQHitVector()
: _vector() {
}

//I will add the rest if I feel like it is needed
// SQHitVector::SQHitVector(const SQHitVector& hitvector)
//   : _vector() {
//   for (std::vector<SQHit*>::const_iterator iter = hitvector.begin();
//       iter != hitvector.end();
//       ++iter) {
//     const SQHit *hit = *iter;
//     _vector.push_back(hit->Clone());
//   }
// }

// SQHitVector& SQHitVector::operator=(const SQHitVector& hitvector) {
//   Reset();
//   for (std::vector<SQHit*>::const_iterator iter = hitvector.begin();
//       iter != hitvector.end();
//       ++iter) {
//     const SQHit *hit = *iter;
//     _vector.push_back(hit->Clone());
//   }
//   return *this;
// }

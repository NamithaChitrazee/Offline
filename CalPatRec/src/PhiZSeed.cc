#include "Offline/CalPatRec/inc/PhiZSeed.hh"
#include <numeric>

namespace mu2e {

// Add hit to the seed
void PhiZSeed::AddHit(HitData_t* Hd) {
    //if(!Hd) return;
    fHitList.push_back(Hd);
    fNHits++;

    fNStrawHits += 1;           // assuming one straw per HitData_t
}

// Clear the seed
void PhiZSeed::Clear() {
    fHitList.clear();
    fNHits = 0;
    fNStrawHits = 0;
    fChi2Par = 0;
    fGood = 1;
    fPhiZIndex = -1;
}


} // namespace mu2e

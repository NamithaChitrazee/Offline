#ifndef CalPatRec_PhiZCandidate_hh
#define CalPatRec_PhiZCandidate_hh

#include "Offline/CalPatRec/inc/PhiZSeed.hh"
#include "Offline/CalPatRec/inc/PhiZFinder_structures.hh"

namespace mu2e {

struct PhiZCandidate {
    int fIndex;
    int fMask;
    int fFirstStation;
    int fLastStation;
    PhiZSeed* fSeed[kNStations];
    int fNSeeds;
    int fNHits;
    int fNStrawHits;

    PhiZCandidate();
    PhiZCandidate(int Index, PhiZSeed* Seed);

    void AddSeed(PhiZSeed* Seed);
};

} // namespace mu2e
#endif

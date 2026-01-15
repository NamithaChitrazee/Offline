#ifndef CalPatRec_PhiZSeed_hh
#define CalPatRec_PhiZSeed_hh

#include <vector>
#include "Offline/CalPatRec/inc/HitData_t.hh"

namespace mu2e {

using CalPatRec::HitData_t;

class PhiZSeed {
public:
    PhiZSeed(int Index)
      : fIndex(Index),
        fStation(0),
        fGood(1),
        fNHits(0),
        fNStrawHits(0),
        fChi2Par(0),
        fPhiZIndex(-1) {}

    ~PhiZSeed() {}

    // Add hit to the seed
    void AddHit(HitData_t* Hd);

    // Clear all hits
    void Clear();

    // Accessors
    HitData_t* Hit(int i) const { return fHitList[i]; }
    int nHits() const { return fNHits; }
    int nStrawHits() const { return fNStrawHits; }

    int Index() const { return fIndex; }
    int Station() const { return fStation; }
    int Good() const { return fGood; }

    void SetStation(int Station) { fStation = Station; }
    void SetPhiZIndex(int Index) { fPhiZIndex = Index; }

    void Init(HitData_t* Hd0) { AddHit(Hd0); }


private:
    int fIndex;          // seed index
    int fStation;        // station number
    int fGood;           // flag

    int fNHits;          // total number of hits
    int fNStrawHits;     // total number of straw hits

    float fChi2Par;      // chi2 of the hits

    int fPhiZIndex;      // candidate index

    std::vector<HitData_t*> fHitList; // store all hits in the seed
};

} // namespace mu2e
#endif

#include "Offline/CalPatRec/inc/PhiZCandidate.hh"

namespace mu2e {

PhiZCandidate::PhiZCandidate() :
    fIndex(-1), fMask(0), fFirstStation(999), fLastStation(-1),
    fNSeeds(0), fNHits(0), fNStrawHits(0) {
    for(int s=0; s<kNStations; ++s) fSeed[s] = nullptr;
}

PhiZCandidate::PhiZCandidate(int Index, PhiZSeed* Seed) : PhiZCandidate() {
    fIndex = Index;
    if (Seed) AddSeed(Seed);
}

void PhiZCandidate::AddSeed(PhiZSeed* Seed) {
    int station = Seed->Station();
    fSeed[station] = Seed;
    fNSeeds++;

    if(fFirstStation > station) fFirstStation = station;
    if(fLastStation < station) fLastStation = station;

    fNHits += Seed->nHits();
    fNStrawHits += Seed->nStrawHits();
}
} // namespace mu2e

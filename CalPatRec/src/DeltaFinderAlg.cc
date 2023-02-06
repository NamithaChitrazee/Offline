///////////////////////////////////////////////////////////////////////////////
// DeltaFinderAlg uses a face-based (instead of a panel-based) data organization
///////////////////////////////////////////////////////////////////////////////
#include "Offline/CalPatRec/inc/DeltaFinderAlg.hh"

#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

namespace mu2e {

  using namespace DeltaFinderTypes;

  using DeltaFinderTypes::PanelZ_t;

  DeltaFinderAlg::DeltaFinderAlg(const fhicl::Table<DeltaFinderAlg::Config>& config, Data_t* Data) :
    _timeBin               (config().timeBin()          ),
    _minHitTime            (config().minHitTime()       ),
    _maxDeltaEDep          (config().maxDeltaEDep()     ),
    _maxSeedEDep           (config().maxSeedEDep()      ),
    _minProtonSeedEDep     (config().minProtonSeedEDep()),
    _minNSeeds             (config().minNSeeds()        ),
    _minDeltaNHits         (config().minDeltaNHits()    ),
    _maxEleHitEnergy       (config().maxEleHitEnergy()  ),
    _minT                  (config().minimumTime()      ), // nsec
    _maxT                  (config().maximumTime()      ), // nsec
    _maxHitSeedDt          (config().maxHitSeedDt()     ), // nsec
    _maxChi2Seed           (config().maxChi2Seed()      ),
    _scaleTwo              (config().scaleTwo()         ),
    _maxChi2Radial         (config().maxChi2Radial()    ),
    _maxChi2All            (config().maxChi2All()       ),
    _maxChi2SeedDelta      (config().maxChi2SeedDelta() ),
    _rCore                 (config().rCore()            ),
    _maxGap                (config().maxGap()           ),
    _sigmaR                (config().sigmaR()           ),
    _sigmaR2               (_sigmaR*_sigmaR             ),
    _maxDriftTime          (config().maxDriftTime()     ),
    _maxSeedDt             (config().maxSeedDt()        ),
    _maxHitDt              (config().maxHitDt()         ),
    _maxStrawDt            (config().maxStrawDt()       ),
    _maxDtDs               (config().maxDtDs()          ),
    _maxDtDc               (config().maxDtDc()          ),
    _debugLevel            (config().debugLevel()        ),
    _diagLevel             (config().diagLevel()         ),
    _printErrors           (config().printErrors()       ),
    _testOrder             (config().testOrder()         ),
    _testHitMask           (config().testHitMask()       ),
    _goodHitMask           (config().goodHitMask()       ),
    _bkgHitMask            (config().bkgHitMask()        ),
    _updateSeedCOG         (config().updateSeedCOG()     )
  {

    _data    = Data;

    printf("DeltaFinderAlg created\n");
  }
//-----------------------------------------------------------------------------
// make sure the two hits used to make a new seed are not a part of an already found seed
//-----------------------------------------------------------------------------
  int DeltaFinderAlg::checkDuplicates(int Station, int Face1, const HitData_t* Hit1, int Face2, const HitData_t* Hit2) {

    int rc(0);

    int nseeds = _data->NSeeds(Station);
    for (int i=0; i<nseeds; i++) {
      DeltaSeed* seed = _data->deltaSeed(Station,i);
      if ((seed->HitData(Face1) == Hit1) and (seed->HitData(Face2) == Hit2)) {
//-----------------------------------------------------------------------------
// 'seed' contains both Hit1 and Hit2 in faces Face1 and Face2 correspondingly, done
//-----------------------------------------------------------------------------
        rc = 1;
                                                                          break;
      }
    }
    return rc;
  }

//-----------------------------------------------------------------------------
// input 'Seed' has two hits , non parallel wires
// try to add more close hits to it (one hit per face)
//-----------------------------------------------------------------------------
  void DeltaFinderAlg::completeSeed(DeltaSeed* Seed) {

    assert ((Seed->SFace(1) >= 0) and (Seed->NHits() == 2));
//-----------------------------------------------------------------------------
// loop over remaining faces, 'f2' - face in question
// bins are 50ns wide, so only need to loop over hits in 3 of them
//-----------------------------------------------------------------------------
    int first_tbin(0), last_tbin(_maxT/_timeBin), max_bin(_maxT/_timeBin);

    int station = Seed->Station();
    float tseed = Seed->TMean();

    int time_bin = (int) (tseed/_timeBin);
    if (time_bin >       0) first_tbin = time_bin-1;
    if (time_bin < max_bin) last_tbin  = time_bin+1;

    float xseed   = Seed->Xc ();
    float yseed   = Seed->Yc ();
    float rho     = sqrt(xseed*xseed+yseed*yseed);
    float seed_nx = xseed/rho;
    float seed_ny = yseed/rho;

    for (int face=0; face<kNFaces; face++) {
      if (Seed->fFaceProcessed[face] == 1)                            continue;
//-----------------------------------------------------------------------------
// face is different from the two first faces used
//-----------------------------------------------------------------------------
      FaceZ_t* fz    = &_data->fFaceData[station][face];
      int      ftbin = first_tbin;
      int      ltbin = last_tbin;

      while ((ftbin<ltbin) and (fz->fFirst[ftbin] < 0)) ftbin++;
      while ((ltbin>ftbin) and (fz->fFirst[ltbin] < 0)) ltbin--;

      int first = fz->fFirst[ftbin];
      if (first < 0) continue;
      int last  = fz->fLast [ltbin];

      HitData_t* closest_hit(nullptr);
      float      best_chi2  (_maxChi2Radial);

      for (int ih=first; ih<=last; ih++) {
//-----------------------------------------------------------------------------
// don't try to re-use hits which already made it into a 3- and 4-hit seeds
//-----------------------------------------------------------------------------
        HitData_t* hd      = &fz->fHitData[ih];
        if (hd->Used() >= 3)                                       continue;
        float corr_time    = hd->fCorrTime;
        const ComboHit* ch = hd->fHit;

        int ip = ch->strawId().panel()/2;
        Pzz_t* pz  = fz->Panel(ip);
//-----------------------------------------------------------------------------
// check seed-panel overlap in phi
//-----------------------------------------------------------------------------
        float ns_dot_np = seed_nx*pz->nx+seed_ny*pz->ny;
        if (ns_dot_np < 0.5)                                        continue;
//-----------------------------------------------------------------------------
// 2017-10-05 PM: consider all hits
// hit time should be consistent with the already existing times - the difference
// between any two measured hit times should not exceed _maxDriftTime
// (_maxDriftTime represents the maximal drift time in the straw, should there be some tolerance?)
// FIXME need some timing checks here !
//-----------------------------------------------------------------------------
        if (corr_time-Seed->T0Max() > _maxHitSeedDt)                break;
        if (Seed->T0Min()-corr_time > _maxHitSeedDt)                continue;

        float xc(xseed), yc(yseed), chi2_par(0), chi2_perp(0), chi2(0);

        if (_updateSeedCOG != 0) {
//-----------------------------------------------------------------------------
// try updating the seed candidate coordinates with the hit added
// to see if that could speed the code up by improvind the efficiency
// of picking up the right hits
// from now on, the default
//-----------------------------------------------------------------------------
          double nr   = ch->pos().x()*pz->wy-ch->pos().y()*pz->wx;

          double nx2m = (Seed->fSnx2+pz->wx*pz->wx);
          double nxym = (Seed->fSnxy+pz->wx*pz->wy);
          double ny2m = (Seed->fSny2+pz->wy*pz->wy);
          double nxrm = (Seed->fSnxr+pz->wx*nr    );
          double nyrm = (Seed->fSnyr+pz->wy*nr    );

          double d    = nx2m*ny2m-nxym*nxym;

          xc          = (nyrm*nx2m-nxrm*nxym)/d;
          yc          = (nyrm*nxym-nxrm*ny2m)/d;

          float dx    = ch->pos().x()-xc;
          float dy    = ch->pos().y()-yc;

          float dxy_dot_w = dx*pz->wx+dy*pz->wy;
          float dxy_dot_n = dx*pz->nx+dy*pz->ny;
          float drho      = fmax(fabs(dxy_dot_n)-_rCore,0);
          chi2_par        = (dxy_dot_w*dxy_dot_w)/(hd->fSigW2+_sigmaR2);
          chi2_perp       = (drho*drho)/_sigmaR2;
          chi2            = chi2_par+chi2_perp;
          if (chi2 < best_chi2) {

            float seed_chi2_par(0), seed_chi2_perp(0);

            seedChi2(Seed,xc,yc,seed_chi2_par,seed_chi2_perp);

            chi2_par  += seed_chi2_par;
            chi2_perp += seed_chi2_perp;
            chi2       = (chi2_par+chi2_perp)/(Seed->NHits()+1);
          }
        }
        else {
//-----------------------------------------------------------------------------
// approximate the chi2 based only on the new hit residual from the seed
// by itself, this calculation is faster, however the approach has a disadvantage
// of not being symmetric over the hits, which could introduce ghost seeds and
// a slowdown due to that.... the outcome is difficult to predict,
// need to check experimentally
// split into wire parallel and perpendicular components
// assume wire direction vector is normalized to 1
//-----------------------------------------------------------------------------
          float dx        = ch->pos().x()-xc;
          float dy        = ch->pos().y()-yc;

          float dxy_dot_w = dx*pz->wx+dy*pz->wy;
          float dxy_dot_n = dx*pz->nx+dy*pz->ny;
          float drho      = fmax(fabs(dxy_dot_n)-_rCore,0);

          chi2_par        = (dxy_dot_w*dxy_dot_w)/(hd->fSigW2+_sigmaR2);
          chi2_perp       = (drho*drho)/_sigmaR2;
//-----------------------------------------------------------------------------
// in this case, only one hit contributes to the chi^2 , no renormalization
//-----------------------------------------------------------------------------
          chi2            = chi2_par+chi2_perp;
        }

        if (chi2 < best_chi2) {
                                        // new best hit
          closest_hit = hd;
          best_chi2   = chi2;
        }
      }

      if (closest_hit) {
//-----------------------------------------------------------------------------
// add hit
//-----------------------------------------------------------------------------
        closest_hit->fChi2Min = best_chi2;
        Seed->AddHit(closest_hit,face);
      }
    }
//-----------------------------------------------------------------------------
// update seed time and X and Y coordinates, accurate knowledge of Z is not very relevant
// also define the number of hits on a found seed
//-----------------------------------------------------------------------------
    Seed->CalculateCogAndChi2(_rCore,_sigmaR2);
    for (int face=0; face<kNFaces; face++) {
      HitData_t* hd = Seed->HitData(face);
      if (hd) hd->fUsed = Seed->NHits();
    }
  }

//-----------------------------------------------------------------------------
// the process moves upstream
//-----------------------------------------------------------------------------
  void DeltaFinderAlg::connectSeeds() {

    for (int is=kNStations-1; is>=0; is--) {
//-----------------------------------------------------------------------------
// 1. loop over existing compton seeds and match them to existing delta candidates
//    number of deltas may increase as the execution moves from one station to another
//-----------------------------------------------------------------------------
      int ndelta = _data->nDeltaCandidates();
      int nseeds = _data->NComptonSeeds(is);
      for (int ids=0; ids<nseeds; ids++) {
        DeltaSeed* seed = _data->ComptonSeed(is,ids);

        if (! seed->Good() )                                          continue;
        if (seed->Used()   )                                          continue;
//-----------------------------------------------------------------------------
// first, loop over existing delta candidates and try to associate the seed
// with one of them
//-----------------------------------------------------------------------------
        DeltaCandidate* closest(nullptr);
        float           chi2min (_maxChi2SeedDelta);

        for (int idc=0; idc<ndelta; idc++) { // this loop creates new deltas, do not overloop
          DeltaCandidate* dc = _data->deltaCandidate(idc);
//-----------------------------------------------------------------------------
// skip candidates already merged with others
//-----------------------------------------------------------------------------
          if (dc->Active() == 0      )                                continue;
          if (dc->seed[is] != nullptr)                                continue;
          int first = dc->FirstStation();
//-----------------------------------------------------------------------------
// a delta candidate starting from a seed in a previous station may already have
// a seed in this station found, skip such a candidate
//-----------------------------------------------------------------------------
          if (first == is)                                            continue;
          int gap  = first-is;
          if (gap > _maxGap)                                          continue;
          float t0 = dc->T0(is);
          assert(t0 > 0);
          float dt = t0-seed->TMean();
          if (fabs(dt) > _maxSeedDt+_maxDtDs*gap)                     continue;
//-----------------------------------------------------------------------------
// the time is OK'ish - checks should be implemented more accurately (FIXME)
// look at the coordinates
//-----------------------------------------------------------------------------
          float chi2 = seedDeltaChi2(seed,dc);
          if (chi2 < chi2min) {
//-----------------------------------------------------------------------------
// everything matches, new closest delta candidate
//-----------------------------------------------------------------------------
            closest = dc;
            chi2min = chi2;
          }
        }
//-----------------------------------------------------------------------------
// if a DeltaSeed has been "attached" to a DeltaCandidate, this is it.
//-----------------------------------------------------------------------------
        if (closest) {
          closest->AddSeed(seed,is);
          seed->fChi2Delta = chi2min;
                                                                      continue;
        }
//-----------------------------------------------------------------------------
// DeltaSeed has not been linked to any existing delta candidate, create
// a new delta candidate and see if it is good enough
//-----------------------------------------------------------------------------
        DeltaCandidate delta(ndelta,seed,is);
//-----------------------------------------------------------------------------
// first, try to extend it backwards, in a compact way, to pick missing
// seeds or hits
//-----------------------------------------------------------------------------
        for (int is2=is+1; is2<kNStations; is2++) {
          recoverStation(&delta,delta.fLastStation,is2,1,1);
//-----------------------------------------------------------------------------
// continue only if found something, allow one empty station
//-----------------------------------------------------------------------------
          if (is2-delta.fLastStation > 2) break;
        }
//-----------------------------------------------------------------------------
// next, try to extend it upstream, by one step, use of seeds is allowed
//-----------------------------------------------------------------------------
        if (is > 0) {
          recoverStation(&delta,is,is-1,1,1);
        }
//-----------------------------------------------------------------------------
// store only delta candidates with hits in 2 stations or more - as the station
// lever arm is about 8 cm, a normal track segment may look like a delta
// for each station define expected T0min and T0max
// to keep _minNSeeds=2 need to look forward...
// 'addDeltaCandidate' makes a deep copy
//-----------------------------------------------------------------------------
        if (delta.fNSeeds >= _minNSeeds) {
          _data->addDeltaCandidate(&delta);
        }
        else {
//-----------------------------------------------------------------------------
// mark all seeds as unassigned, this should be happening only in 1-seed case
//-----------------------------------------------------------------------------
          for (int is=delta.fFirstStation; is<=delta.fLastStation; is++) {
            DeltaSeed* ds = delta.seed[is];
            if (ds) ds->fDeltaIndex = -1;
          }
        }
      }
//-----------------------------------------------------------------------------
// all seeds in a given station processed
// loop over existing delta candidates which do not have seeds in a given station
// and see if can pick up single hits
// do not extrapolate forward by more than two (?) empty stations - leave it
// as a constant for the moment
// 'ndelta' doesn't count deltas starting in this station, but there is no need
// to loop over them either
//-----------------------------------------------------------------------------
      for (int idc=0; idc<ndelta; idc++) {
        DeltaCandidate* dc = _data->deltaCandidate(idc);
        int first = dc->FirstStation();
        if (first == is )                                             continue;
        if (first-is > 3)                                             continue;
//-----------------------------------------------------------------------------
// if a delta candidate has been created in this routine, time limits
// may not be defined. Make sure they are
//-----------------------------------------------------------------------------
        recoverStation(dc,first,is,1,0);
      }
    }
//-----------------------------------------------------------------------------
// at this point we have a set of delta candidates, which might need to be merged
//-----------------------------------------------------------------------------
    mergeDeltaCandidates();
  }

//-----------------------------------------------------------------------------
// find delta electron seeds in 'Station' with hits in faces 'f' and 'f+1'
// do not consider proton hits with eDep > _minHtEnergy
//-----------------------------------------------------------------------------
  void DeltaFinderAlg::findSeeds(int Station, int Face) {

    FaceZ_t* fz1 = &_data->fFaceData[Station][Face];
    int      nh1 = fz1->fHitData.size();
//-----------------------------------------------------------------------------
// modulo misalignments, panels in stations 2 and 3 are oriented exactly the same
// way as in stations 0 and 1, etc
//-----------------------------------------------------------------------------
    for (int h1=0; h1<nh1; ++h1) {
//-----------------------------------------------------------------------------
// hit has not been used yet to start a seed, however it could've been used as a second seed
//-----------------------------------------------------------------------------
      HitData_t*      hd1 = &fz1->fHitData[h1];
      if (hd1->Used() >= 3)                                           continue;
      const ComboHit* ch1 = hd1->fHit;
      int seed_found = 0;
//-----------------------------------------------------------------------------
// panels 0,2,4 are panels 0,1,2 in the first  (#0) face of a plane
// panels 1,3,5 are panels 0,1,2 in the second (#1) face
//-----------------------------------------------------------------------------
      int    ip1 = ch1->strawId().panel() / 2;
      Pzz_t* pz1 = fz1->Panel(ip1);
      float  wx1 = pz1->wx;
      float  wy1 = pz1->wy;

      float  x1  = ch1->pos().x();
      float  y1  = ch1->pos().y();
      float  t1  = ch1->time();
//-----------------------------------------------------------------------------
// figure out the first and the last timing bins to loop over
// loop over 3 bins (out of > 20) - the rest cant contain hits of interest
//-----------------------------------------------------------------------------
      int time_bin = (int) t1/_timeBin;
      int first_tbin(0), last_tbin(_maxT/_timeBin), max_bin(_maxT/_timeBin);

      if (time_bin >       0) first_tbin = time_bin-1;
      if (time_bin < max_bin) last_tbin  = time_bin+1;
//-----------------------------------------------------------------------------
// loop over 'next' faces
// timing bins may be empty...
//-----------------------------------------------------------------------------
      for (int f2=Face+1; f2<kNFaces; f2++) {
        FaceZ_t* fz2   = &_data->fFaceData[Station][f2];
        float    zc    = (fz1->z+fz2->z)/2;

        int      ftbin = first_tbin;
        int      ltbin = last_tbin;

        while ((ftbin<ltbin) and (fz2->fFirst[ftbin] < 0)) ftbin++;
        while ((ltbin>ftbin) and (fz2->fFirst[ltbin] < 0)) ltbin--;
        int first = fz2->fFirst[ftbin];
        if (first < 0)                                             continue;

        int last  = fz2->fLast [ltbin];
        for (int h2=first; h2<=last; h2++) {
          HitData_t*      hd2 = &fz2->fHitData[h2];
          if (hd2->Used() >= 3)                                    continue;
          const ComboHit* ch2 = hd2->fHit;
          float t2 = ch2->time();
          float dt = t2-t1;

          if (dt >  _maxDriftTime) break;
          if (dt < -_maxDriftTime) continue;
//-----------------------------------------------------------------------------
// 'ip2' - panel index within its face
// check overlap in phi between the panels coresponding to the wires - 120 deg
//-----------------------------------------------------------------------------
          int    ip2  = ch2->strawId().panel() / 2;
          Pzz_t* pz2  = fz2->Panel(ip2);
          float  n1n2 = pz1->nx*pz2->nx+pz1->ny*pz2->ny;
          if (n1n2 < -0.5)                                          continue;
//-----------------------------------------------------------------------------
// hits are consistent in time,
//-----------------------------------------------------------------------------
          float x2    = ch2->pos().x();
          float y2    = ch2->pos().y();

          double wx2   = pz2->wx;
          double wy2   = pz2->wy;
          double w1w2  = wx1*wx2+wy1*wy2;
          double q12   = 1-w1w2*w1w2;
//-----------------------------------------------------------------------------
// hits are ordered in time, so if ct2-ct > _maxDriftTime, can proceed with the next panel
//-----------------------------------------------------------------------------
// intersect the two straws, we need coordinates of the intersection point and
// two distances from hits to the intersection point, 4 numbers in total
//-----------------------------------------------------------------------------
          double r12n1 = (x1-x2)*wx1+(y1-y2)*wy1;
          double r12n2 = (x1-x2)*wx2+(y1-y2)*wy2;

          double wd1   = -(r12n2*w1w2-r12n1)/q12;

          float  xc    = x1-wx1*wd1;
          float  yc    = y1-wy1*wd1;

          double wd2   = -(r12n2-w1w2*r12n1)/q12;
//-----------------------------------------------------------------------------
// require both hits to be close enough to the intersection point
//-----------------------------------------------------------------------------
          float chi2_hd1 = wd1*wd1/hd1->fSigW2;
          float chi2_hd2 = wd2*wd2/hd2->fSigW2;

          if (chi2_hd1 > _maxChi2Seed)                                continue;
          if (chi2_hd2 > _maxChi2Seed)                                continue;
//-----------------------------------------------------------------------------
// this may be used with some scale factor sf < 2
//-----------------------------------------------------------------------------
          if ((chi2_hd1+chi2_hd2) > _scaleTwo*_maxChi2Seed)           continue;
//-----------------------------------------------------------------------------
// check whether there already is a seed containing both hits
//-----------------------------------------------------------------------------
          int is_duplicate  = checkDuplicates(Station,Face,hd1,f2,hd2);
          if (is_duplicate)                               continue;
//-----------------------------------------------------------------------------
// new seed : an intersection of two wires coresponsing to close in time combo hits
//-----------------------------------------------------------------------------
          hd1->fChi2Min     = chi2_hd1;
          hd2->fChi2Min     = chi2_hd2;

          DeltaSeed* seed   = _data->NewDeltaSeed(Station,Face,hd1,f2,hd2);

          seed->CofM.SetXYZ(xc,yc,zc);
//-----------------------------------------------------------------------------
// mark both hits as a part of a seed, so they would not be used individually
// - see HitData_t::Used()
//-----------------------------------------------------------------------------
          hd1->fSeed  = seed;
          hd2->fSeed  = seed;
//-----------------------------------------------------------------------------
// complete search for hits of this seed, mark it BAD (or 'not-LEE') if a proton
// in principle, could place "high-charge" seeds into a separate list
// that should improve the performance
// if a seed EDep > _maxSeedEDep       (5 keV), can't be a low energy electron (LEE)
// if a seed EDep > _minProtonSeedEDep (3 keV), could be a proton
//-----------------------------------------------------------------------------
          completeSeed(seed);

          if (seed->Chi2TotN() > _maxChi2Seed) {
//-----------------------------------------------------------------------------
// discard found seed
//-----------------------------------------------------------------------------
            seed->fGood = -3000-seed->fIndex;
          }
          else {
            if (seed->EDep() > _maxSeedEDep) {
              seed->fGood = -2000-seed->fIndex;
            }
            else {
              _data->AddComptonSeed(seed,Station);
            }

            if (seed->EDep() > _minProtonSeedEDep) {
              _data->AddProtonSeed(seed,Station);
            }
            seed_found = seed->NHits();
          }
//-----------------------------------------------------------------------------
// if found seed has hits in 3 or 4 faces, use next first hit
//-----------------------------------------------------------------------------
          if (seed_found >= 3) break;
        }
        if (seed_found >= 3) break;
      }
    }
  }

//-----------------------------------------------------------------------------
// TODO: update the time as more hits are added
//-----------------------------------------------------------------------------
  void DeltaFinderAlg::findSeeds() {

    for (int s=0; s<kNStations; ++s) {
      for (int face=0; face<kNFaces-1; face++) {
//-----------------------------------------------------------------------------
// find seeds starting from 'face' in a given station 's'
//-----------------------------------------------------------------------------
        findSeeds(s,face);
      }
      pruneSeeds(s);
    }
  }

//-----------------------------------------------------------------------------
// merge Delta Candidates : check for duplicates !
//-----------------------------------------------------------------------------
  int DeltaFinderAlg::mergeDeltaCandidates() {
    int rc(0);
    float max_d2(20*20);  // mm^2, to be adjusted FIXME

    int ndelta = _data->nDeltaCandidates();

    for (int i1=0; i1<ndelta-1; i1++) {
      DeltaCandidate* dc1 = _data->deltaCandidate(i1);
      if (dc1->Active() == 0)                                          continue;
      float x1 = dc1->CofM.x();
      float y1 = dc1->CofM.y();
      for (int i2=i1+1; i2<ndelta; i2++) {
        DeltaCandidate* dc2 = _data->deltaCandidate(i2);
        if (dc2->Active() == 0)                                        continue;
        float x2 = dc2->CofM.x();
        float y2 = dc2->CofM.y();
        float d2 = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
        if (d2 > max_d2)                                               continue;
//-----------------------------------------------------------------------------
// too lazy to extrapolate the time to the same Z...  ** FIXME
//-----------------------------------------------------------------------------
        float t1 = dc1->T0(dc1->LastStation());
        float t2 = dc2->T0(dc2->FirstStation());
//-----------------------------------------------------------------------------
// time check could be done more intelligently - compare time at the same Z
//-----------------------------------------------------------------------------
        if (fabs(t1-t2) > _maxDtDc)                                    continue;

        int dds = dc2->FirstStation()-dc1->LastStation();
        if (dds < 0) {
          if (_printErrors) {
            printf("ERROR in DeltaFinderAlg::%s:",__func__);
            printf("i1, i2, dc1->LastStation, dc2->FirstStation: %2i %2i %2i %2i \n",
                   i1,i2,dc1->LastStation(),dc2->FirstStation());
          }
                                                                       continue;
        }
        else if (dds > _maxGap)                                        continue;
//-----------------------------------------------------------------------------
// merge two delta candidates, not too far from each other in Z,
// leave dc1 active, mark dc2 as not active
//-----------------------------------------------------------------------------
        dc1->MergeDeltaCandidate(dc2,_printErrors);
        dc2->SetIndex(-1000-dc1->Index());
      }
    }
    return rc;
  }

//------------------------------------------------------------------------------
// I'd love to use the hit flags, however that is confusing:
// - hits with very large deltaT get placed to the middle of the wire and not flagged,
// - however, some hits within the fiducial get flagged with the ::radsel flag...
// use only "good" hits
//-----------------------------------------------------------------------------
  int DeltaFinderAlg::orderHits() {
    ChannelID cx, co;
//-----------------------------------------------------------------------------
// vector of pointers to CH, ordered in time. Initial list is not touched
//-----------------------------------------------------------------------------
    _data->_v.resize(_data->_nComboHits);

    for (int i=0; i<_data->_nComboHits; i++) {
      _data->_v[i] = &(*_data->chcol)[i];
    }

    std::sort(_data->_v.begin(), _data->_v.end(),
              [](const ComboHit*& a, const ComboHit*& b) { return a->time() < b->time(); });
//-----------------------------------------------------------------------------
// at this point hits in '_v' are already ordered in time
//-----------------------------------------------------------------------------
    for (int ih=0; ih<_data->_nComboHits; ih++) {
      const ComboHit* ch = _data->_v[ih];

      const StrawHitFlag* flag   = &ch->flag();
      if (_testHitMask && (! flag->hasAllProperties(_goodHitMask) || flag->hasAnyProperty(_bkgHitMask)) ) continue;

      // float corr_time    = ch->correctedTime();

      cx.Station                 = ch->strawId().station();
      cx.Plane                   = ch->strawId().plane() % 2;
      cx.Face                    = -1;
      cx.Panel                   = ch->strawId().panel();
//-----------------------------------------------------------------------------
// get Z-ordered location
//-----------------------------------------------------------------------------
      Data_t::orderID(&cx, &co);

      int os       = co.Station;
      int of       = co.Face;
      int op       = co.Panel;

      if (_printErrors) {
        if ((os < 0) || (os >= kNStations     )) printf(" >>> ERROR: wrong station number: %i\n",os);
        if ((of < 0) || (of >= kNFaces        )) printf(" >>> ERROR: wrong face    number: %i\n",of);
        if ((op < 0) || (op >= kNPanelsPerFace)) printf(" >>> ERROR: wrong panel   number: %i\n",op);
      }
//-----------------------------------------------------------------------------
// prototype face-based hit storage
// hits are already time-ordered - that makes it easy to define fFirst
// for each face, define multiple time bins and indices of the first and the last
// hits in each bin
//-----------------------------------------------------------------------------
      FaceZ_t* fz  = &_data->fFaceData[os][of];
      int      loc = fz->fHitData.size();

      fz->fHitData.push_back(HitData_t(ch,of));
      int time_bin = int (ch->time()/_timeBin) ;

      if (fz->fFirst[time_bin] < 0) fz->fFirst[time_bin] = loc;
      fz->fLast[time_bin] = loc;
    }

    return 0;
  }

//-----------------------------------------------------------------------------
// some of found seeds could be duplicates or ghosts
// in case two DeltaSeeds share the first seed hit, leave only the best one
// the seeds we're loooping over have been reconstructed within the same station
// also reject seeds with Chi2Tot > _maxChi2Tot=10
//-----------------------------------------------------------------------------
  void DeltaFinderAlg::pruneSeeds(int Station) {

    int nseeds =  _data->NSeeds(Station);

    for (int i1=0; i1<nseeds-1; i1++) {
      DeltaSeed* ds1 = _data->deltaSeed(Station,i1);
      if (ds1->fGood < 0)                                             continue;

      if (ds1->Chi2AllN() > _maxChi2All) {
        ds1->fGood = -1000-i1;
                                                                      continue;
      }

      float tmean1 = ds1->TMean();

      for (int i2=i1+1; i2<nseeds; i2++) {
        DeltaSeed* ds2 = _data->deltaSeed(Station,i2);
        if (ds2->fGood < 0)                                           continue;

        if (ds2->Chi2AllN() > _maxChi2All) {
          ds2->fGood = -1000-i2;
                                                                      continue;
        }

        float tmean2 = ds2->TMean();

        if (fabs(tmean1-tmean2) > _maxSeedDt)                         continue;
//-----------------------------------------------------------------------------
// the two segments are close in time , both have acceptable chi2's
// *FIXME* didn't check distance !!!!!
// so far, allow duplicates during the search
// the two DeltaSeeds share could have significantly overlapping hit content
//-----------------------------------------------------------------------------
        int noverlap            = 0;
        int nfaces_with_overlap = 0;
        for (int face=0; face<kNFaces; face++) {
          const HitData_t* hh1 = ds1->HitData(face);
          const HitData_t* hh2 = ds2->HitData(face);
          if (hh1 and (hh1 == hh2)) {
            noverlap            += 1;
            nfaces_with_overlap += 1;
          }
        }

        if (nfaces_with_overlap > 1) {
//-----------------------------------------------------------------------------
// overlap significant, leave in only one DeltaSeed - which one?
//-----------------------------------------------------------------------------
          if (ds1->fNFacesWithHits > ds2->fNFacesWithHits) {
            ds2->fGood = -1000-i1;
          }
          else if (ds2->fNFacesWithHits > ds1->fNFacesWithHits) {
            ds1->fGood = -1000-i2;
            break;
          }
          else {
//-----------------------------------------------------------------------------
//both seeds have the same number of hits - compare chi2's
//-----------------------------------------------------------------------------
            if (ds1->Chi2AllN() <  ds2->Chi2AllN()) {
              ds2->fGood = -1000-i1;
            }
            else {
              ds1->fGood = -1000-i2;
              break;
            }
          }
        }
        else if (nfaces_with_overlap > 0) {
//-----------------------------------------------------------------------------
// only one overlapping hit
// special treatment of 2-hit seeds to reduce the number of ghosts
//-----------------------------------------------------------------------------
          if (ds1->fNFacesWithHits == 2) {
            if (ds2->fNFacesWithHits > 2) {
              ds1->fGood = -1000-i2;
              break;
            }
            else {
//-----------------------------------------------------------------------------
// the second seed also has 2 faces with hits
//-----------------------------------------------------------------------------
              if (ds1->Chi2AllN() <  ds2->Chi2AllN()) ds2->fGood = -1000-i1;
              else {
                ds1->fGood = -1000-i2;
                break;
              }
            }
          }
          else {
//-----------------------------------------------------------------------------
// the first seed has N>2 hits
//-----------------------------------------------------------------------------
            if (ds2->fNFacesWithHits == 2) {
//-----------------------------------------------------------------------------
// the 2nd seed has only 2 hits and there is an overlap
//-----------------------------------------------------------------------------
              ds2->fGood = -1000-i1;
            }
            else {
//-----------------------------------------------------------------------------
// the second seed also has N>2 faces with hits, but there is only one overlap
// leave both seeds in
//-----------------------------------------------------------------------------
            }
          }
        }
      }
    }
  }

//------------------------------------------------------------------------------
// start from looking at the "holes" in the seed pattern
// delta candidates in the list are already required to have at least 2 segments
// extend them outwards by one station
//-----------------------------------------------------------------------------
  int DeltaFinderAlg::recoverMissingHits() {

    int ndelta = _data->nDeltaCandidates();
    for (int idelta=0; idelta<ndelta; idelta++) {
      DeltaCandidate* dc = _data->deltaCandidate(idelta);
//-----------------------------------------------------------------------------
// don't extend candidates made out of one segment - but there is no such
// start from the first station to define limits
//-----------------------------------------------------------------------------
      int s1 = dc->fFirstStation;
      int s2 = dc->fLastStation-1;
      int last(-1);
//-----------------------------------------------------------------------------
// first check inside "holes"
//-----------------------------------------------------------------------------
      for (int i=s1; i<=s2; i++) {
        if (dc->seed[i] != nullptr) {
          last  = i;
          continue;
        }
//-----------------------------------------------------------------------------
// define expected T0 limits
//-----------------------------------------------------------------------------
        recoverStation(dc,last,i,1,0);
      }

      last  = dc->fFirstStation;
      for (int i=last-1; i>=0; i--) {
//-----------------------------------------------------------------------------
// skip empty stations
//-----------------------------------------------------------------------------
        if (dc->fFirstStation -i > _maxGap) break;
        recoverStation(dc,dc->fFirstStation,i,1,0);
      }

      last = dc->fLastStation;
      for (int i=last+1; i<kNStations; i++) {
//-----------------------------------------------------------------------------
// skip empty stations
//-----------------------------------------------------------------------------
        if (i-dc->fLastStation > _maxGap) break;
        recoverStation(dc,dc->fLastStation,i,1,0);
      }
    }

    return 0;
  }

//-----------------------------------------------------------------------------
// return 1 if a seed has been found , 0 otherwise
//-----------------------------------------------------------------------------
  int DeltaFinderAlg::recoverSeed(DeltaCandidate* Delta, int LastStation, int Station) {
    // int rc(0);
                                        // predicted time range for this station
    float tdelta = Delta->T0(Station);

    float      chi2min(_maxChi2SeedDelta);
    DeltaSeed* closest_seed(nullptr);

    float dt   = _maxSeedDt + _maxDtDs*fabs(Station-LastStation);

    int nseeds = _data->NComptonSeeds(Station);
    for (int i=0; i<nseeds; i++) {
      DeltaSeed* seed =  _data->ComptonSeed(Station,i);
      if (seed->Good() == 0)                                          continue;
      if (seed->Used()     )                                          continue;
//-----------------------------------------------------------------------------
// one might need some safety here, but not the _maxDriftTime
//-----------------------------------------------------------------------------
      if (fabs(tdelta-seed->TMean()) > dt)                            continue;

      float chi2 = seedDeltaChi2(seed,Delta);

      if (chi2 < chi2min) {
                                        // new best seed
        closest_seed = seed;
        chi2min      = chi2;
      }
    }

    if (closest_seed) {
//-----------------------------------------------------------------------------
// the closest seed found, add it to the delta candidate and exit
// it is marked as associated with the delta candidate in DeltaCandidate::AddSeed
//-----------------------------------------------------------------------------
      Delta->AddSeed(closest_seed,Station);
      closest_seed->fChi2Delta = chi2min;
    }

    return (closest_seed != nullptr);
  }

//------------------------------------------------------------------------------
// try to recover hits of a 'Delta' candidate in a given 'Station'
// 'Delta' doesn't have hits in this station, check all hits here
// when predicting time, use the same value of Z for both layers of a given face
// return 1 if something has been found
//-----------------------------------------------------------------------------
  int DeltaFinderAlg::recoverStation(DeltaCandidate* Delta, int LastStation, int Station, int UseUsedHits, int RecoverSeeds) {

                                        // predicted time range for this station
    float tdelta   = Delta->T0(Station);
    float xdelta   = Delta->Xc();
    float ydelta   = Delta->Yc();
    float rdelta   = sqrt(xdelta*xdelta+ydelta*ydelta);
    float delta_nx = xdelta/rdelta;
    float delta_ny = ydelta/rdelta;

    int first_tbin(0), last_tbin(_maxT/_timeBin), max_bin(_maxT/_timeBin);

    int time_bin   = (int) (tdelta/_timeBin);

    if (time_bin >       0) first_tbin = time_bin-1;
    if (time_bin < max_bin) last_tbin  = time_bin+1;
//-----------------------------------------------------------------------------
// first, loop over the existing seeds - need for the forward step
//-----------------------------------------------------------------------------
    if (RecoverSeeds) {
      int closest_seed_found = recoverSeed(Delta,LastStation,Station);
      if (closest_seed_found == 1) return 1;
    }
//-----------------------------------------------------------------------------
// no seeds found, look for single hits in the 'Station'
//-----------------------------------------------------------------------------
    DeltaSeed*  new_seed (nullptr);

    float dt_hit = _maxHitDt+_maxDtDs*fabs(Station-LastStation);

    for (int face=0; face<kNFaces; face++) {
      FaceZ_t* fz = &_data->fFaceData[Station][face];
      int ftbin = first_tbin;
      int ltbin = last_tbin;

      while ((ftbin<ltbin) and (fz->fFirst[ftbin] < 0)) ftbin++;
      while ((ltbin>ftbin) and (fz->fFirst[ltbin] < 0)) ltbin--;
      int first = fz->fFirst[ftbin];
      if (first < 0) continue;
      int last  = fz->fLast [ltbin];
//-----------------------------------------------------------------------------
// loop over hits in this face
//-----------------------------------------------------------------------------
      for (int h=first; h<=last; h++) {
        HitData_t* hd = &fz->fHitData[h];
        if (hd->Used() >= 3)                                          continue;
        if ((UseUsedHits == 0) and hd->Used())                        continue;
//-----------------------------------------------------------------------------
// don't skip hits already included into seeds - a two-hit stereo seed
// could be random
//-----------------------------------------------------------------------------
        float corr_time = hd->fCorrTime;
//-----------------------------------------------------------------------------
// can gain a little bit by checking the time, leave the check in
// predicted time is the particle time, the drift time should be larger
//-----------------------------------------------------------------------------
        if (corr_time > tdelta+dt_hit)                              break;
        if (corr_time < tdelta-dt_hit)                              continue;

        const ComboHit*  ch = hd->fHit;
        int    ip   = ch->strawId().panel()/2;
        Pzz_t* pz   = fz->Panel(ip);
        float  n1n2 = pz->nx*delta_nx+pz->ny*delta_ny;
//-----------------------------------------------------------------------------
// 0.5 corresponds to delta(phi) = 60 deg
//-----------------------------------------------------------------------------
        if (n1n2 < 0.5)                                               continue;
//-----------------------------------------------------------------------------
// the hit is consistent with the Delta in phi and time, check more accurately
//-----------------------------------------------------------------------------
        double dx       = ch->pos().x()-xdelta;
        double dy       = ch->pos().y()-ydelta;

        double dw       = dx*pz->wx+dy*pz->wy;           // distance along the wire
        double dr       = dx*pz->nx+dy*pz->ny;

        float chi2_par  = (dw*dw)/(hd->fSigW2+_sigmaR2);
        float ddr       = fmax(fabs(dr)-_rCore,0);
        float chi2_perp = (ddr*ddr)/_sigmaR2;
        float chi2      = chi2_par + chi2_perp;

        if (chi2 >= _maxChi2Radial)                                 continue;
        if (hd->Used()) {
//-----------------------------------------------------------------------------
// hit is a part of a seed. if the seed has 2 or less hits, don't check the chi2
// - that could be a random overlap
// if the seed has 3 or more hits, check the chi2
//-----------------------------------------------------------------------------
          int nh = hd->fSeed->NHits();
          if ((nh >= 3) and (chi2 > hd->fChi2Min))                   continue;
        }
//-----------------------------------------------------------------------------
// new hit needs to be added, create a special 1-hit seed for that
// in most cases, expect this seed not to have the second hit, but it may
// such a seed has its own CofM undefined
// ** FIXME ..in principle, at this point may want to check if the hit was used...
//-----------------------------------------------------------------------------
        if (new_seed == nullptr) {
          hd->fChi2Min = chi2;
          new_seed = _data->NewDeltaSeed(Station,face,hd,-1,nullptr);
        }
        else {
          if (face == new_seed->SFace(0)) {
//-----------------------------------------------------------------------------
// another close hit in the same panel, choose the best
//-----------------------------------------------------------------------------
            if (chi2 >= new_seed->HitData(face)->fChi2Min)          continue;
//-----------------------------------------------------------------------------
// new best hit in the same face
//-----------------------------------------------------------------------------
            hd->fChi2Min = chi2;
            new_seed->ReplaceFirstHit(hd);
          }
          else {
//-----------------------------------------------------------------------------
// more than one hit added in the hit pickup mode. The question is why a seed,
// constructed out of those two hits has not been added (or created)
//-----------------------------------------------------------------------------
            if (_printErrors) {
              printf("ERROR in DeltaFinderAlg::recoverStation: ");
              printf("station=%2i - shouldn\'t be getting here, printout of new_seed and hd follows\n",Station);
              printf("chi2_par, chi2_perp, chi2: %8.2f %8.2f %8.2f\n",chi2_par, chi2_perp, chi2);

              printf("DELTA:\n");
              _data->printDeltaCandidate(Delta,"");
              printf("SEED:\n");
              _data->printDeltaSeed(new_seed,"");
              printf("HIT:\n");
              _data->printHitData  (hd      ,"");
            }

            new_seed->AddHit(hd,face);
          }

          if (corr_time < new_seed->fMinHitTime) new_seed->fMinHitTime = corr_time;
          if (corr_time > new_seed->fMaxHitTime) new_seed->fMaxHitTime = corr_time;
        }
      }
      if (new_seed) new_seed->fFaceProcessed[face] = 1;
    }
//-----------------------------------------------------------------------------
// station is processed, see if anything has been found
// some parameters of seeds found in a recovery mode are not defined because
// there was no pre-seeding, for example
//-----------------------------------------------------------------------------
    int rc(0);
    if (new_seed) {
      int face0                       = new_seed->SFace(0);
      new_seed->HitData(face0)->fSeed = new_seed;

      Delta->AddSeed(new_seed,Station);
      rc = 1;
    }
                                        // return 1 if hits were added
    return rc;
  }

//-----------------------------------------------------------------------------
  void  DeltaFinderAlg::run() {
    orderHits();
//-----------------------------------------------------------------------------
// loop over all stations and find delta seeds - 2-3-4 combo hit stubs
// a seed is always a stereo object
//-----------------------------------------------------------------------------
    findSeeds();
//-----------------------------------------------------------------------------
// connect seeds and create delta candidates
// at this stage, extend seeds to pick up single its in neighbor stations
// for single hits do not allo gaps
//-----------------------------------------------------------------------------
    connectSeeds();
//-----------------------------------------------------------------------------
// for existing delta candidates, pick up single gap-separated single hits
// no new candidates is created at this step
//-----------------------------------------------------------------------------
    recoverMissingHits();
//-----------------------------------------------------------------------------
// after recovery of missing hits, it is possible that some of delta candidates
// may need to be merged - try again
//-----------------------------------------------------------------------------
    mergeDeltaCandidates();
//-----------------------------------------------------------------------------
// finally, mark deltas
//-----------------------------------------------------------------------------
    int ndeltas = _data->nDeltaCandidates();

    for (int i=0; i<ndeltas; i++) {
      DeltaCandidate* dc = _data->deltaCandidate(i);
//-----------------------------------------------------------------------------
// skip merged in delta candidates
// also require a delta candidate to have at least 5 hits
// do not consider proton stub candidates (those with <EDep> > 0.004)
//-----------------------------------------------------------------------------
      if (dc->Active() == 0)                                          continue;

      if (dc->NHits () < _minDeltaNHits) dc->fMask |= DeltaCandidate::kNHitsBit;
      if (dc->EDep  () > _maxDeltaEDep ) dc->fMask |= DeltaCandidate::kEDepBit;
    }
  }

//-----------------------------------------------------------------------------
  void DeltaFinderAlg::seedChi2(DeltaSeed* Seed, float Xc, float Yc, float& Chi2Par, float& Chi2Perp) {
    Chi2Par  = 0;
    Chi2Perp = 0;

    for (int face=0; face<kNFaces; face++) {
      const HitData_t* hd = Seed->HitData(face);
      if (hd == nullptr)                                              continue;
      float dx = hd->fHit->pos().x()-Xc;
      float dy = hd->fHit->pos().y()-Yc;
//-----------------------------------------------------------------------------
// split into wire parallel and perpendicular components
//-----------------------------------------------------------------------------
      const XYZVectorF& wdir = hd->fHit->wdir();

      float dxy_dot_w = dx*wdir.x()+dy*wdir.y();
      float dxy_dot_n = dx*wdir.y()-dy*wdir.x();

      float  chi2_par = (dxy_dot_w*dxy_dot_w)/(_sigmaR2+hd->fSigW2);
      float drr       = fmax(fabs(dxy_dot_n)-_rCore,0);
      float chi2_perp = (drr*drr)/_sigmaR2;

      Chi2Par        += chi2_par;
      Chi2Perp       += chi2_perp;
    }
  }

//-----------------------------------------------------------------------------
// in principle, can just use Delta Xc,Yc
//-----------------------------------------------------------------------------
  float DeltaFinderAlg::seedDeltaChi2(DeltaSeed* Seed, DeltaCandidate* Delta) {

    // double nxym = (Delta->fSnxy+Seed->fSnxy);
    // double nx2m = (Delta->fSnx2+Seed->fSnx2);
    // double ny2m = (Delta->fSny2+Seed->fSny2);
    // double nxrm = (Delta->fSnxr+Seed->fSnxr);
    // double nyrm = (Delta->fSnyr+Seed->fSnyr);

    // double d    = nx2m*ny2m-nxym*nxym;
    // float  xc   = (nyrm*nx2m-nxrm*nxym)/d;
    // float  yc   = (nyrm*nxym-nxrm*ny2m)/d;

    float  xc   = Delta->Xc();
    float  yc   = Delta->Yc();

    float chi2, chi2_par, chi2_perp;
    seedChi2(Seed, xc,yc, chi2_par, chi2_perp);
    chi2 = (chi2_par+chi2_perp)/Seed->NHits();

    return chi2;
  }

}
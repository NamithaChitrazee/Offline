///////////////////////////////////////////////////////////////////////////////
// PhiZSeedFinderAlg uses a face-based (instead of a panel-based) data organization
///////////////////////////////////////////////////////////////////////////////
#include "Offline/CalPatRec/inc/ChannelID.hh"
#include "Offline/CalPatRec/inc/PhiZSeedFinderAlg.hh"

using CalPatRec::ChannelID;

namespace mu2e {

  using namespace PhiZSeedFinderTypes;

  PhiZSeedFinderAlg::PhiZSeedFinderAlg(const fhicl::Table<PhiZSeedFinderAlg::Config>& config, Data_t* Data) :
    _debugLevel            (config().debugLevel()       ),
    _diagLevel             (config().diagLevel()        ),
    _testOrder             (config().testOrder()        )
  {

    _data    = Data;
//-----------------------------------------------------------------------------
// as the algorithm is supposed to work only on time clusters, make sure that
// only one time bin is defined.
// could think of further simplification down the road
//-----------------------------------------------------------------------------
    _timeBin = 2000;

    printf("PhiZSeedFinderAlg created\n");
  }

//------------------------------------------------------------------------------
// all hits belong to the same time cluster, order them in Z and in time
//-----------------------------------------------------------------------------
  int PhiZSeedFinderAlg::orderHits(const TimeCluster* Tc) {
/*
   std::cout << "PhiZSeedFinderAlg_orderHits " << nseeds << "\n"

    ChannelID cx, co;
//-----------------------------------------------------------------------------
// vector of pointers to CH, ordered in time. Initial list is not touched
//-----------------------------------------------------------------------------
    _data->_nComboHits = Tc->nhits();

    _data->_v.resize(_data->_nComboHits);

    for (int i=0; i<_data->_nComboHits; i++) {
      const StrawHitIndex& ind = Tc->hits().at(i);
      _data->_v[i] = &(*_data->chcol)[ind];
    }
    std::cout<<"nComboHits = "<<_data->_nComboHits<<std::endl;
    //std::sort(_data->_v.begin(), _data->_v.end(),
    //          [](const ComboHit*& a, const ComboHit*& b) { return a->time() < b->time(); });
    std::sort(_data->_v.begin(), _data->_v.end(),
          [](const ComboHit*& a, const ComboHit*& b) {
              return a->pos().z() < b->pos().z();
          });

//-----------------------------------------------------------------------------
// at this point hits in '_v' are already ordered in time
//-----------------------------------------------------------------------------
    for (int ih=0; ih<_data->_nComboHits; ih++) {
      const ComboHit* ch = _data->_v[ih];

      cx.Station                 = ch->strawId().station();
      cx.Plane                   = ch->strawId().plane() % 2;
      cx.Face                    = -1;
      cx.Panel                   = ch->strawId().panel();
//-----------------------------------------------------------------------------
// get Z-ordered location
//-----------------------------------------------------------------------------
      ChannelID::orderID(&cx, &co);

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
*/
    return 0;
  }

//-----------------------------------------------------------------------------
// find delta electron seeds in 'Station' with hits in faces 'f' and 'f+1'
// do not consider proton hits with eDep > _minHtEnergy
//-----------------------------------------------------------------------------
void PhiZSeedFinderAlg::findSeeds(int Station, int Face) {

    std::cout << "=====================================\n";
    std::cout << "findSeeds\nstation/face " << Station << ", face " << Face << "\n";
    std::cout << "=====================================\n";

    auto* fz1 = _data->faceData(Station, Face);
    int nh1 = fz1->fHitData.size();

    if (nh1 == 0) {
        std::cout << "No hits in station " << Station << ", face " << Face << std::endl;
        return;
    }

    // Create a new seed
    PhiZSeed* seed = _data->newPhiZSeed();
    seed->SetStation(Station); // if you still want to record station

    for (int h1 = 0; h1 < nh1; ++h1) {
        HitData_t* hd = &fz1->fHitData[h1];
        const auto& pos = hd->fHit->pos();

        //Preselect hits
        int hitIndex = hd->fHit->index(0);  // <-- corrected
        //double phi   = hd->fHit->phi();     // already exists
        double x = pos.x();
        double y = pos.x();
        double z = pos.x();
        if(hitIndex > 10000) continue;
        if(sqrt(x*x + y*y) < 300) continue;
        if(abs(z) > 1520) continue;

        std::cout << "Hit #" << h1
                  << " | x=" << pos.x()
                  << " y=" << pos.y()
                  << " z=" << pos.z() << std::endl;

        if (h1 == 0) {
            seed->Init(hd);
        } else {
            seed->AddHit(hd);
        }
    }

    std::cout << "PhiZSeed created at station " << Station
              << " with " << seed->nHits() << " hits.\n";



//-----------------------------------------------------------------------------
// hit has not been used yet to start a seed, however it could've been used as a second seed
//-----------------------------------------------------------------------------
/*      HitData_t*      hd1 = &fz1->fHitData[h1];
      if (hd1->Used() >= 3)                                           continue;

      float  wx1 = hd1->fWx;
      float  wy1 = hd1->fWy;
      float  x1  = hd1->fX;
      float  y1  = hd1->fY;

      const ComboHit* ch1 = hd1->fHit;
      int   seed_found    = 0;
//-----------------------------------------------------------------------------
// panels 0,2,4 are panels 0,1,2 in the first  (#0) face of a plane
// panels 1,3,5 are panels 0,1,2 in the second (#1) face
//-----------------------------------------------------------------------------
      int    ip1 = ch1->strawId().panel() / 2;
      Pzz_t* pz1 = fz1->Panel(ip1);
//-----------------------------------------------------------------------------
// figure out the first and the last timing bins to loop over
// loop over 3 bins (out of > 20) - the rest cant contain hits of interest
//-----------------------------------------------------------------------------
      float  t1       = ch1->time();
      int    time_bin = (int) t1/_timeBin;

      int    first_tbin(0), last_tbin(_maxT/_timeBin), max_bin(_maxT/_timeBin);

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
        if (first < 0)                                                continue;

        int last  = fz2->fLast [ltbin];
        for (int h2=first; h2<=last; h2++) {
          HitData_t*      hd2 = &fz2->fHitData[h2];
          if (hd2->Used() >= 3)                                       continue;
          const ComboHit* ch2 = hd2->fHit;
          float t2 = ch2->time();
          float dt = t2-t1;

          if (dt < -_maxDriftTime)                                    continue;
          if (dt >  _maxDriftTime)                                    break;
//-----------------------------------------------------------------------------
// the following check relies on the TOT... not quite sure yet
// however, it also makes sense to require that both pulses have a reasonable width,
// so leave it in for the moment
//-----------------------------------------------------------------------------
          float dtcorr = hd1->fCorrTime-hd2->fCorrTime;
          if (fabs(dtcorr) > _maxDriftTime)                           continue;
//-----------------------------------------------------------------------------
// 'ip2' - panel index within its face
// check overlap in phi between the panels coresponding to the wires - 120 deg
//-----------------------------------------------------------------------------
          int    ip2  = ch2->strawId().panel() / 2;
          Pzz_t* pz2  = fz2->Panel(ip2);
          float  n1n2 = pz1->nx*pz2->nx+pz1->ny*pz2->ny;
          if (n1n2 < -0.5)                                            continue;
//-----------------------------------------------------------------------------
// hits are consistent in time,
//-----------------------------------------------------------------------------
          float x2     = hd2->fX;
          float y2     = hd2->fY;

          double wx2   = hd2->fWx;
          double wy2   = hd2->fWy;
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

          if (chi2_hd1 > _maxChi2Par)                                 continue;
          if (chi2_hd2 > _maxChi2Par)                                 continue;
//-----------------------------------------------------------------------------
// this may be used with some scale factor sf < 2
// the following line is a provision for future...
//-----------------------------------------------------------------------------
          float chi2_time = (dtcorr*dtcorr)/sigma_dt_2;
          float chi2_tot  = chi2_time+(chi2_hd1+chi2_hd2)/2;

          if (chi2_tot > _maxChi2Seed)                                continue;
//-----------------------------------------------------------------------------
// check whether there already is a seed containing both hits
//-----------------------------------------------------------------------------
          int is_duplicate  = checkDuplicates(Station,Face,hd1,f2,hd2);
          if (is_duplicate)                                           continue;
//-----------------------------------------------------------------------------
// new seed : an intersection of two wires coresponsing to close in time combo hits
//-----------------------------------------------------------------------------
          hd1->fChi2Min     = chi2_hd1;
          hd2->fChi2Min     = chi2_hd2;

          // DeltaSeed* seed   = _data->NewDeltaSeed(Station,hd1,hd2,xc,yc,zc);
          DeltaSeed* seed   = _data->newDeltaSeed(Station);
          seed->Init(hd1,hd2,xc,yc,zc);
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
// if the seed EDep > _maxSeedEDep       (5 keV), can't be a low energy electron (LEE)
// if  seed EDep > _minProtonSeedEDep (3 keV), could be a proton
//-----------------------------------------------------------------------------
          completeSeed(seed);

          if (seed->Chi2TotN() > _maxChi2Seed) {
//-----------------------------------------------------------------------------
// discard found seed
//-----------------------------------------------------------------------------
            seed->fGood = -3000-seed->fIndex;
          }
          else {
//-----------------------------------------------------------------------------
// lists of proton and compton seeds are not mutually exclusive -
// some (3 keV < EDep < 5 keV) could be either
//-----------------------------------------------------------------------------
            int n_high_edep_hits(0);
            for (int face=0; face<kNFaces; face++) {
              HitData_t* hit = seed->HitData(face);
              if (hit and (hit->fHit->energyDep() > _minProtonHitEDep)) n_high_edep_hits++;
            }
            seed->SetNHighEDepHits(n_high_edep_hits);

            if (seed->EDep() > _maxSeedEDep)        seed->fGood = -2000-seed->fIndex;
            else {
              _data->AddComptonSeed(seed,Station);
            }

            if ((seed->EDep() > _minProtonSeedEDep) and (n_high_edep_hits >= 2)) {
//-----------------------------------------------------------------------------
// make sure the number of proton-like hits, hits above proton E(min) is at least two
//-----------------------------------------------------------------------------
              _data->AddProtonSeed (seed,Station);
            }

            seed_found = seed->nHits();
          }
//-----------------------------------------------------------------------------
// if found seed has hits in 3 or 4 faces, use next first hit
//-----------------------------------------------------------------------------
          if (seed_found >= 3) break;
        }
        if (seed_found >= 3) break;
      }
*/
  }


//-----------------------------------------------------------------------------
// TODO: update the time as more hits are added
//-----------------------------------------------------------------------------
  void PhiZSeedFinderAlg::findSeeds() {

    std::cout<<"findSeeds "<<std::endl;
    std::cout<<"kNStations = "<<kNStations<<std::endl;
    for (int s= kNStations -1; s>=0; s--) {
      std::cout<<"kNFaces = "<<kNFaces<<std::endl;
      for (int face=0; face<kNFaces; face++) {

      std::cout<<"station/face = "<<s<<"/"<<face<<std::endl;
//-----------------------------------------------------------------------------
// find seeds starting from 'face' in a given station 's'
//-----------------------------------------------------------------------------
        findSeeds(s,face);
      }
      //pruneSeeds(s);
    }
    std::cout<<"FindSeeds END"<<std::endl;
  }



//-----------------------------------------------------------------------------
  void  PhiZSeedFinderAlg::run(const TimeCluster* Tc) {
    orderHits(Tc);
    //-----------------------------------------------------------------------------
    // loop over all stations and find delta seeds - 2-3-4 combo hit stubs
    // a seed is always a stereo object
    //-----------------------------------------------------------------------------
    //findSeeds();

    /*int nseeds = _data->nSeeds();
    int count = 0;
    std::cout << "Total seeds found: " << nseeds << "\n";
for (int i = 0; i < _data->nSeeds(); ++i) {
    PhiZSeed* seed = _data->seed(i);
    std::cout << "Seed #" << i
              << " (station " << seed->Station() << ")"
              << " with " << seed->nHits() << " hits\n";

    for (int j = 0; j < seed->nHits(); ++j) {
          count++;
        const auto* hd = seed->Hit(j);
        const auto& pos = hd->fHit->pos();

        int hitIndex = hd->fHit->index(0);  // <-- corrected
        double phi   = hd->fHit->phi();     // already exists

        std::cout << "   hit " << j
                  << " index=" << hitIndex
                  << " station=" << seed->Station()
                  << " x=" << pos.x()
                  << " y=" << pos.y()
                  << " z=" << pos.z()
                  << " phi=" << phi
                  << "\n";
    }
}

std::cout<<"count = "<<count<<std::endl;
*/












  }

}

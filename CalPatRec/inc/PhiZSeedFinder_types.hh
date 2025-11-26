#ifndef CalPatRec_PhiZSeedFinder_types_hh
#define CalPatRec_PhiZSeedFinder_types_hh

namespace art {
  class Event;
}

namespace fhicl {
  class ParameterSet;
}

#include "TObject.h"
#include "TClonesArray.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"

#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/RecoDataProducts/inc/HelixHit.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/CalPatRec/inc/ChannelID.hh"


#include "Offline/Mu2eUtilities/inc/ManagedList.hh"
#include "Offline/CalPatRec/inc/CalPatRec_enums.hh"
#include "Offline/CalPatRec/inc/Pzz_t.hh"
#include "Offline/CalPatRec/inc/HitData_t.hh"
#include "Offline/CalPatRec/inc/PhiZFinder_structures.hh"
#include "Offline/CalPatRec/inc/PhiZSeed.hh"
#include "Offline/CalPatRec/inc/PhiZCandidate.hh"

using CalPatRec::ChannelID;
using CalPatRec::Pzz_t;
using mu2e::CalPatRec::HitData_t;

namespace mu2e {
  class Panel;
  class SimParticle;
  class PhiZSeedFinderAlg;
//-----------------------------------------------------------------------------
// delta-electron seed: structure within the station
// doesn't own anything, no need to delete any pinters
//-----------------------------------------------------------------------------
  namespace PhiZSeedFinderTypes {

    struct Config {
      fhicl::Atom<std::string> tool_type             {fhicl::Name("tool_type"             ), fhicl::Comment("tool type: PhiZSeedFinderDiag")  };
      fhicl::Atom<int>         mcTruth               {fhicl::Name("mcTruth"               ), fhicl::Comment("MC truth")                       };
      fhicl::Atom<int>         diagLevel             {fhicl::Name("diagLevel"             ), fhicl::Comment("diagnostic level")               };
      fhicl::Atom<bool>        mcDiag                {fhicl::Name("mcDiag"                ), fhicl::Comment("MC diag")                        };
      fhicl::Atom<int>         printOTracker         {fhicl::Name("printOTracker"         ), fhicl::Comment("print ordered Tracker")          };
      fhicl::Atom<int>         printComboHits        {fhicl::Name("printComboHits"        ), fhicl::Comment("print combo hits")               };
      fhicl::Atom<int>         printGoodComboHits    {fhicl::Name("printGoodComboHits"    ), fhicl::Comment("print good combo hits")          };
      fhicl::Atom<int>         printShcol            {fhicl::Name("printShcol"            ), fhicl::Comment("if 1, print shColl"             )};

      fhicl::Table<McUtilsToolBase::Config> mcUtils  {fhicl::Name("mcUtils"               ), fhicl::Comment("MC Diag plugin"                 )};
    };

    struct FaceZ_t {
      int                     fID;         // 3*face+panel, for pre-calculating overlaps

      std::vector<HitData_t>  fHitData;
      int                     fFirst[100];   // ** FIXME - need larger dimension for off-spill cosmics...
      int                     fLast [100];

      Pzz_t                   fPanel[3];
      double                  z;           //

      Pzz_t*                  Panel(int I) { return &fPanel[I]; }

      int                     nHits        () { return fHitData.size(); }
      HitData_t*              hitData      (int I) { return &fHitData     [I]; }
    };

    struct Data_t {
      const art::Event*             event;
      const Tracker*                tracker;
      const DiskCalorimeter*        calorimeter;

      art::InputTag                 chCollTag;
      art::InputTag                 tcCollTag;

      const ComboHitCollection*     chcol;
      const TimeClusterCollection*  tccol;
      const StrawDigiMCCollection*  sdmcColl;

      PhiZSeedFinderAlg*            _finder;

      int                           debugLevel;              // printout level

      int                           _nTimeClusters;
      int                           _nComboHits;
      int                           _nStrawHits;
      std::vector<const ComboHit*>  _v;                      // sorted

      //ManagedList<PhiZSeed>        fListOfSeeds       [kNStations];
      ManagedList<PhiZSeed>        fListOfSeeds;
      std::vector<PhiZCandidate>   fListOfPhiZCandidates;
//-----------------------------------------------------------------------------
// try to avoid looping over panels
//-----------------------------------------------------------------------------
      FaceZ_t                       fFaceData   [kNStations][kNFaces];
      int                           stationUsed [kNStations];
//-----------------------------------------------------------------------------
// station #2 is the same as station #0 etc...
//-----------------------------------------------------------------------------
      int                           panelOverlap[2][12][12];
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
      Data_t();
      ~Data_t();

      PhiZCandidate*  phizCandidate (int I)             { return &fListOfPhiZCandidates [I]; }
      int  nPhiZCandidates ()        { return fListOfPhiZCandidates.size (); }

      FaceZ_t* faceData   (int Station, int Face) { return &fFaceData[Station][Face]; }

      //class in PhiZSeedFinder_types
      void     InitEvent(const art::Event* Evt, int DebugLevel);
      void     InitGeometry();

      void addPhiZCandidate(PhiZCandidate* PhiZ) { fListOfPhiZCandidates.push_back(*PhiZ); }

      /*PhiZSeed*  newPhiZSeed(int Station) {
        PhiZSeed* dc = fListOfSeeds[Station].New();
        dc->SetStation(Station);
        return dc;
      }*/
      PhiZSeed* newPhiZSeed() {
        PhiZSeed* dc = fListOfSeeds.New();
        return dc;
      }
      int nSeeds() { return fListOfSeeds.N(); }
      PhiZSeed* seed(int i) { return fListOfSeeds.at(i); }

//-----------------------------------------------------------------------------
// diagInfo for PhiZSeed
//-----------------------------------------------------------------------------
      TimeClusterCollection*  _tccolnew;// new TimeCluster of segment candidates
//-----------------------------------------------------------------------------
// diagInfo for MC Info
//-----------------------------------------------------------------------------
      std::vector<double> h_MCnStrawHitsPerParticle;
      std::vector<double> h_MCnComboHitsPerParticle;
      std::vector<double> h_MCMom;
      std::vector<double> h_MCTransverseMom;
      std::vector<double> h_MCTanLambda;

//-----------------------------------------------------------------------------
// diagInfo for HelixFind
//-----------------------------------------------------------------------------
      HelixSeedCollection* _hscolnew;// new HelixSeed of segment candidates
      std::vector<double> h_lineFitter_chi2Dof[3];// chi2Dof of linear fit for each small segment candidates before merging [/per small segment candidate]
      std::vector<double> h_circleFitter_chi2Dof[3];// chi2Dof of circle fit with combohits [/per segment candidate]
      std::vector<int> h_circleFitter_nhits[3];// chi2Dof of circle fit with combohits [/per segment candidate]
      std::vector<double> h_diffradius;//diff between = (helix fit radius - MC radius);
      std::vector<double> h_fraction;//fraction = number of missed hits in a heilx - number of expected hits in a helix;
    };

  }// namespace PhiZSeedFinderTypes
}// namespace mu2e
#endif

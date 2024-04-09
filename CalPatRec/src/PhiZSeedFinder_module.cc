/////////////////////////////////////////////////////////////////////////////
// P.Murat
//
// flag combo hits as 'delta'  - hits of identified low energy electrons
//                and 'proton' - hits of identified protons/deuterons
//
// always writes out a ComboHitCollection with correct flags,
// to be used by downstream modules
//
// WriteFilteredComboHits = 0: write out all hits
//                        = 1: write out only hits not flagged as 'delta' or 'proton'
//                             (to be used in trigger)
//
// parameter defaults: CalPatRec/fcl/prolog.fcl
//////////////////////////////////////////////////////////////////////////////
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"

#include "art/Utilities/make_tool.h"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

// conditions
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"

// data
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"

// diagnostics
#include "Offline/CalPatRec/inc/ChannelID.hh"
#include "Offline/CalPatRec/inc/PhiZSeedFinder_types.hh"
#include "Offline/CalPatRec/inc/PhiZSeedFinderAlg.hh"

//ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include <TROOT.h>
#include "TLine.h"
#include "TEllipse.h"
#include <TPaveText.h>
#include <TSystem.h>
#include "TLatex.h"
#include <TF1.h>
#include <TStyle.h>

// C++
#include <algorithm>
#include <cmath>
#include <vector>

using namespace std;

using CalPatRec::ChannelID;

namespace mu2e {

  using namespace PhiZSeedFinderTypes;

  class PhiZSeedFinder: public art::EDProducer {
  public:

    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>   sschCollTag            {Name("sschCollTag"        )     , Comment("SS ComboHit collection tag" ) };
      fhicl::Atom<art::InputTag>   chCollTag              {Name("chCollTag"          )     , Comment("ComboHit collection tag"    ) };
      fhicl::Atom<art::InputTag>   tcCollTag              {Name("tcCollTag"          )     , Comment("time cluster collection tag") };
      fhicl::Atom<int>             debugLevel              {Name("debugLevel"        )     , Comment("debug level"                ) };
      fhicl::Atom<int>             diagLevel              {Name("diagLevel"          )     , Comment("diag level"                  ) };
      fhicl::Atom<int>             printErrors            {Name("printErrors"        )     , Comment("print errors"                ) };
      fhicl::Atom<int>             writeFilteredComboHits {Name("writeFilteredComboHits"), Comment("0: write all CH, 1: write filtered CH") };
      fhicl::Atom<int>             writeStrawHits          {Name("writeStrawHits"    )     , Comment("1: write all SH, new flags" ) };
      fhicl::Atom<int>             testOrder              {Name("testOrder"          )     , Comment("1: test order"              ) };
      fhicl::Atom<bool>             testHitMask            {Name("testHitMask"        )     , Comment("true: test hit mask"        ) };
      fhicl::Sequence<std::string> goodHitMask            {Name("goodHitMask"        )     , Comment("good hit mask"              ) };
      fhicl::Sequence<std::string> bkgHitMask              {Name("bkgHitMask"        )     , Comment("background hit mask"        ) };

      fhicl::Table<PhiZSeedFinderTypes::Config> diagPlugin      {Name("diagPlugin"      ), Comment("Diag plugin"           ) };
      fhicl::Table<PhiZSeedFinderAlg::Config>    finderParameters{Name("finderParameters"), Comment("finder alg parameters" ) };
    };

//-----------------------------------------------------------------------------
// PhiZSeedFinder Constructors
//-----------------------------------------------------------------------------
     struct ev5_HitsInNthStation {
          float phi;
          int strawhits;
          float x;
          float y;
          float z;
          int station;
          int plane;
          int face;
          int panel;
          int hitID;
      };

       struct ev5_Segment {
          double deltaphi;
          double z;
          double alpha;
          double beta;
          double chiNDF;
          int station;
          int reference_point;
          bool usedforfit;
      };

  protected:
//-----------------------------------------------------------------------------
// talk-to parameters: input collections and algorithm parameters
//-----------------------------------------------------------------------------
    art::InputTag    _sschCollTag;
    art::InputTag    _chCollTag;
    art::InputTag    _tcCollTag;                 // time cluster coll tag

    int              _writeFilteredComboHits;   // write filtered combo hits
    // int             _writeStrawHitFlags;        // obsolete
    int              _writeStrawHits;           // write out filtered (?) straw hits

    int              _debugLevel;
    int              _diagLevel;
    int              _printErrors;
    int              _testOrder;

    StrawHitFlag    _bkgHitMask;

    std::unique_ptr<ModuleHistToolBase> _hmanager;
//-----------------------------------------------------------------------------
// cache event/geometry objects
//-----------------------------------------------------------------------------
    const ComboHitCollection*     _sschColl ;

    const Tracker*               _tracker;
    const DiskCalorimeter*       _calorimeter;

    PhiZSeedFinderTypes::Data_t  _data;               // all data used
    int                           _testOrderPrinted;

    PhiZSeedFinderAlg*           _finder;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    explicit PhiZSeedFinder(const art::EDProducer::Table<Config>& config);

  private:

    bool         findData             (const art::Event&  Evt);
//-----------------------------------------------------------------------------
// overloaded methods of the module class
//-----------------------------------------------------------------------------
    void         beginJob() override;
    void         beginRun(art::Run& ARun) override;
    void         endJob  () override;
    void         produce (art::Event& E ) override;

//-----------------------------------------------------------------------------
// PhiZSeedFinder unique functions (by h.kitagawa)
//-----------------------------------------------------------------------------
    void ev5_FillHitsInTimeCluster(const TimeCluster* Tc, std::vector<ev5_HitsInNthStation>& ComboHitsInCluster);
    void ev5_SegmentSearchInTriplet(const std::vector<ev5_HitsInNthStation> HitsInTimeCluster, std::vector<std::vector<ev5_HitsInNthStation>>& best_triplet_segments, std::vector<std::vector<ev5_Segment>>& diag_best_triplet_segments, double thre_residual);
    double ev5_DeltaPhi(double x1, double y1, double x2, double y2);
    double ev5_ParticleDirection(double x1, double y1, double x2, double y2);
    double ev5_ResidualDeltaPhi(double alpha, double beta, double z1, double phi1, double z2, double phi2);
    bool ev5_TripletQuality(const std::vector<ev5_Segment> diag_hit);
    void ev5_select_best_segments_step_01(const std::vector<ev5_HitsInNthStation>& HitsInTimeCluster, const std::vector<std::vector<ev5_HitsInNthStation>>& segment_candidates, const std::vector<std::vector<ev5_Segment>>& diag_segment_candidates, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi);
    void ev5_select_best_segments_step_02(int station, const std::vector<ev5_HitsInNthStation>& HitsInTimeCluster, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi);
    void ev5_select_best_segments_step_03(const std::vector<ev5_HitsInNthStation>& HitsInTimeCluster, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi);
    void ev5_fit_slope(const std::vector<ev5_Segment>& hit_diag, double& alpha, double& beta, double& chindf);
    void ev5_select_best_segments_step_04(const std::vector<ev5_HitsInNthStation>& HitsInTimeCluster, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi);
    void ev5_select_best_segments_step_05(const std::vector<ev5_HitsInNthStation>& HitsInTimeCluster, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi);
    void ev5_select_best_segments_step_06(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi);
    void ev5_select_best_segments_step_07(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi, int& NumberOfSegments);

//-----------------------------------------------------------------------------
// function
//-----------------------------------------------------------------------------
    void initTimeCluster(TimeCluster& tc);

  };

//-----------------------------------------------------------------------------
  PhiZSeedFinder::PhiZSeedFinder(const art::EDProducer::Table<Config>& config):
    art::EDProducer{config},
    _sschCollTag           (config().sschCollTag()       ),
    _chCollTag             (config().chCollTag()         ),
    _tcCollTag             (config().tcCollTag()         ),
    _debugLevel             (config().debugLevel()         ),
    _diagLevel             (config().diagLevel()         ),
    _printErrors           (config().printErrors()       ),
    _testOrder             (config().testOrder()         ),
    _bkgHitMask             (config().bkgHitMask()         )
  {

    consumes<TimeClusterCollection>(_tcCollTag);
    consumes<ComboHitCollection>   (_chCollTag);

    produces<TimeClusterCollection>();

    _finder = new PhiZSeedFinderAlg(config().finderParameters,&_data);

    _testOrderPrinted = 0;

    if (_diagLevel != 0) _hmanager = art::make_tool  <ModuleHistToolBase>(config().diagPlugin,"diagPlugin");
    else                 _hmanager = std::make_unique<ModuleHistToolBase>();

    _data.chCollTag       = _chCollTag;
    _data.tcCollTag       = _tcCollTag;
    _data._finder         = _finder;          // for diagnostics
  }

  //-----------------------------------------------------------------------------
  void PhiZSeedFinder::beginJob() {
    if (_diagLevel > 0) {
      art::ServiceHandle<art::TFileService> tfs;
      _hmanager->bookHistograms(tfs);
    }
  }

  //-----------------------------------------------------------------------------
  void PhiZSeedFinder::endJob() {
  }

//-----------------------------------------------------------------------------
// create a Z-ordered representation of the tracker
//-----------------------------------------------------------------------------
  void PhiZSeedFinder::beginRun(art::Run& aRun) {

    _data.InitGeometry();
//-----------------------------------------------------------------------------
// it is enough to print that once
//-----------------------------------------------------------------------------
    if (_testOrder && (_testOrderPrinted == 0)) {
      ChannelID::testOrderID  ();
      ChannelID::testdeOrderID();
      _testOrderPrinted = 1;
    }

    if (_diagLevel != 0) _hmanager->debug(&_data,1);
  }

//-----------------------------------------------------------------------------
  bool PhiZSeedFinder::findData(const art::Event& Evt) {

    auto tccH    = Evt.getValidHandle<mu2e::TimeClusterCollection>(_tcCollTag);
    _data.tccol = tccH.product();

    auto chcH = Evt.getValidHandle<mu2e::ComboHitCollection>(_chCollTag);
    _data.chcol = chcH.product();

    return (_data.tccol != nullptr) and (_data.chcol != nullptr);
  }

//-----------------------------------------------------------------------------
    void PhiZSeedFinder::ev5_FillHitsInTimeCluster(const TimeCluster* Tc, std::vector<ev5_HitsInNthStation>& ComboHitsInCluster){
      /*if (tc.size() > 0) {
        for(size_t ipeak=0; ipeak<tc.size(); ipeak++) {
          const std::vector<StrawHitIndex>& TC = tc[ipeak];
          int nh = TC.size();
          for(int i=0; i<nh; i++){
          int ind = TC[i];
          const mu2e::ComboHit* ch = &_chcol->at(ind);
          }
        }
      }*/
    std::cout<<"Tc.size() = "<<Tc->nhits()<<std::endl;
    const std::vector<StrawHitIndex>& ordchcol = Tc->hits();
    int nh = ordchcol.size();
    for (int ih=0; ih<nh; ih++) {
      int ind = ordchcol[ih];
      const ComboHit* ch = &_data.chcol->at(ind);
      ev5_HitsInNthStation hitsincluster;
      hitsincluster.hitID = ch->_sid.asUint16();
      hitsincluster.phi = ch->phi();
      hitsincluster.x = ch->pos().x();
      hitsincluster.y = ch->pos().y();
      hitsincluster.z = ch->pos().z();
      hitsincluster.strawhits = ch->_nsh;
      //hitsincluster.station = ch.strawId().station();
      hitsincluster.station = ch->strawId().station();

      ComboHitsInCluster.push_back(hitsincluster);
    }

  /*  for (int ih=0; ih<_data._nComboHits; ih++) {
      const ComboHit* ch = _data._v[ih];
      ev5_HitsInNthStation hitsincluster;
      hitsincluster.hitID = ch->_sid.asUint16();
      hitsincluster.phi = ch->phi();
      hitsincluster.x = ch->pos().x();
      hitsincluster.y = ch->pos().y();
      hitsincluster.z = ch->pos().z();
      hitsincluster.strawhits = ch->_nsh;
      //hitsincluster.station = ch.strawId().station();
      hitsincluster.station = ch->strawId().station();

      ComboHitsInCluster.push_back(hitsincluster);
    }*/

    //sort the vector in ascending order of the z-coordinate
    std::sort(ComboHitsInCluster.begin(), ComboHitsInCluster.end(), [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b) { return a.z < b.z; } );

}

//-----------------------------------------------------------------------------
      void PhiZSeedFinder::ev5_SegmentSearchInTriplet(const std::vector<ev5_HitsInNthStation> HitsInTimeCluster, std::vector<std::vector<ev5_HitsInNthStation>>& best_triplet_segments, std::vector<std::vector<ev5_Segment>>& diag_best_triplet_segments, double thre_residual) {


      std::vector<std::vector<ev5_HitsInNthStation>> all_ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>> all_ThisIsBestSegment_Diag;
      all_ThisIsBestSegment.clear();
      all_ThisIsBestSegment_Diag.clear();

      // number of stations
      // search slope in 3 consecutive stations
      // Loop from 0 to 16 station
      int nstation = 18;
      for(int n=0; n<nstation-2; n++){

          // Take combo hits in n-th, (n+1)-th, (n+2)-th stations
          std::vector<ev5_HitsInNthStation> Hits_In_Station[3];
          for(int k=0; k<3;k++) Hits_In_Station[k].clear();
          for(int j=0; j<(int)HitsInTimeCluster.size(); j++){
              ev5_HitsInNthStation hitsin_nthstation;
              hitsin_nthstation.phi          = HitsInTimeCluster.at(j).phi;
              hitsin_nthstation.strawhits    = HitsInTimeCluster.at(j).station;
              hitsin_nthstation.x            = HitsInTimeCluster.at(j).x;
              hitsin_nthstation.y            = HitsInTimeCluster.at(j).y;
              hitsin_nthstation.z            = HitsInTimeCluster.at(j).z;
              hitsin_nthstation.station      = HitsInTimeCluster.at(j).station;
              hitsin_nthstation.plane        = HitsInTimeCluster.at(j).plane;
              hitsin_nthstation.face        = HitsInTimeCluster.at(j).face;
              hitsin_nthstation.panel        = HitsInTimeCluster.at(j).panel;
              hitsin_nthstation.hitID        = HitsInTimeCluster.at(j).hitID;
              if(n == HitsInTimeCluster.at(j).station) Hits_In_Station[0].push_back(hitsin_nthstation);
              if(n+1 == HitsInTimeCluster.at(j).station) Hits_In_Station[1].push_back(hitsin_nthstation);
              if(n+2 == HitsInTimeCluster.at(j).station) Hits_In_Station[2].push_back(hitsin_nthstation);
          }

        //2 >= ComboHits [/station]
        size_t nCHsInStn_1 = Hits_In_Station[0].size();
        size_t nCHsInStn_2 = Hits_In_Station[1].size();
        size_t nCHsInStn_3 = Hits_In_Station[2].size();
        int nHitInStation = (int)nCHsInStn_1 + (int)nCHsInStn_2 + (int)nCHsInStn_3;
        if(!(nHitInStation >= 5)) continue;
        // (1st, 2nd, 3rd) = (CH>=1, CH>=1, CH>=1)
        if(!((int)nCHsInStn_1 >= 1)) continue;
        if(!((int)nCHsInStn_2 >= 1)) continue;
        if(!((int)nCHsInStn_3 >= 1)) continue;


/*      std::cout << "-----------------------------------" << std::endl;
      std::cout << "              Station               " << n        << std::endl;
      std::cout << "-----------------------------------" << std::endl;
        // Print informations for n-th and (n+1)-th stations
        std::cout<<"Station = "<<n<<std::endl;
       for(size_t p=0; p<Hits_In_Station[0].size(); p++){
          std::cout<<"phi/nStrawHits/station/plane/face/panel/x/y/z = "<<Hits_In_Station[0].at(p).phi<<"/"<<Hits_In_Station[0].at(p).strawhits<<"/"<<Hits_In_Station[0].at(p).station<<"/"<<Hits_In_Station[0].at(p).plane<<"/"<<Hits_In_Station[0].at(p).face<<"/"<<Hits_In_Station[0].at(p).panel<<"/"<<Hits_In_Station[0].at(p).x<<"/"<<Hits_In_Station[0].at(p).y<<"/"<<Hits_In_Station[0].at(p).z<<std::endl;
        }
        std::cout<<"Station = "<<n+1<<std::endl;
       for(size_t p=0; p<Hits_In_Station[1].size(); p++){
          std::cout<<"phi/nStrawHits/station/plane/face/panel/x/y/z = "<<Hits_In_Station[1].at(p).phi<<"/"<<Hits_In_Station[1].at(p).strawhits<<"/"<<Hits_In_Station[1].at(p).station<<"/"<<Hits_In_Station[1].at(p).plane<<"/"<<Hits_In_Station[1].at(p).face<<"/"<<Hits_In_Station[1].at(p).panel<<"/"<<Hits_In_Station[1].at(p).x<<"/"<<Hits_In_Station[1].at(p).y<<"/"<<Hits_In_Station[1].at(p).z<<std::endl;
        }
        std::cout<<"Station = "<<n+2<<std::endl;
       for(size_t p=0; p<Hits_In_Station[2].size(); p++){
          std::cout<<"phi/nStrawHits/station/plane/face/panel/x/y/z = "<<Hits_In_Station[2].at(p).phi<<"/"<<Hits_In_Station[2].at(p).strawhits<<"/"<<Hits_In_Station[2].at(p).station<<"/"<<Hits_In_Station[2].at(p).plane<<"/"<<Hits_In_Station[2].at(p).face<<"/"<<Hits_In_Station[2].at(p).panel<<"/"<<Hits_In_Station[2].at(p).x<<"/"<<Hits_In_Station[2].at(p).y<<"/"<<Hits_In_Station[2].at(p).z<<std::endl;
        }
*/
        //-------------------------------------------------------------
        //        Find segments in 3 consecutive stations
        //-------------------------------------------------------------
        //Segment candidates with MC true info
        std::vector<std::vector<ev5_HitsInNthStation>> segment_candidates;
        segment_candidates.clear();
        //Segment candidates with diagnostic MC true info
        std::vector<std::vector<ev5_Segment>>  diag_segment_candidates;
        diag_segment_candidates.clear();


        //Find segment candidate
        //Loop 1st station
        for(int j=0; j<(int)nCHsInStn_1; j++){
          //Loop 3rd station
          for(int k=0; k<(int)nCHsInStn_3; k++){
            int flag_hit[3] = {0};
            // Scalar product between 1st and 3rd hit
            // first_Hit = {x, y, z, phi}
            double first_Hit[4] = {Hits_In_Station[0].at(j).x, Hits_In_Station[0].at(j).y, Hits_In_Station[0].at(j).z, Hits_In_Station[0].at(j).phi};
            double last_Hit[4] = {Hits_In_Station[2].at(k).x, Hits_In_Station[2].at(k).y, Hits_In_Station[2].at(k).z, Hits_In_Station[2].at(k).phi};
            double DeltaPhi = ev5_DeltaPhi(first_Hit[0], first_Hit[1], last_Hit[0], last_Hit[1]);
            double ref_sign = ev5_ParticleDirection(first_Hit[0], first_Hit[1], last_Hit[0], last_Hit[1]);//+ is eletron, - is positve particle
            //std::cout<<"DeltaPhi between 1st ST/3rd ST= "<<DeltaPhi<<std::endl;
            //std::cout<<"particle has positive/negative track = "<<ref_sign<<std::endl;

            //deltaPhi cut on 1st and 3rd hit
            if(DeltaPhi > 2.0) continue; //1st hit and 3rd hit in the triplet should be within DeltaPhi < 2.0[rad]

            // Calculate the slope between 1st and 3rd hit
            double phi[3] = {0.0, 0.0, ref_sign*DeltaPhi};
            double z[3] = {first_Hit[2], 0.0, last_Hit[2]};
            double mean_phi = (phi[0] + phi[2])/2.0;
            double mean_z = (z[0] + z[2])/2.0;
            double m_n[3] = {0.0};
            double m_d[3] = {0.0};
            m_n[0] = (z[0] - mean_z)*(phi[0] - mean_phi);
            m_n[2] = (z[2] - mean_z)*(phi[2] - mean_phi);
            m_d[0] = pow(z[0] - mean_z, 2);
            m_d[2] = pow(z[2] - mean_z, 2);
            double slope_alpha = (m_n[0] + m_n[2])/(m_d[0] + m_d[2]);
            double slope_beta = phi[0] - slope_alpha*z[0];
            //std::cout<<"phi1/phi3 = "<<first_Hit[3]<<"/"<<last_Hit[3]<<std::endl;
            //std::cout<<"b/m = "<<b<<"/"<<m<<std::endl;
            //std::cout<<"slope_alpha/beta = "<<slope_alpha<<"/"<<slope_beta<<std::endl;
            //std::cout<<"a = "<<(phi[2]-phi[0])/(fabs(z[2]-z[0]))<<std::endl;

            //Fill 1st hit info
            std::vector<ev5_HitsInNthStation> hit_candidates;
            hit_candidates.clear();
            hit_candidates.push_back(Hits_In_Station[0].at(j));
            std::vector<ev5_Segment> SegmentInTripletStation;
            SegmentInTripletStation.clear();
            ev5_Segment hit_first;
            hit_first.deltaphi = 0.0;
            hit_first.z = Hits_In_Station[0].at(j).z;
            hit_first.alpha = slope_alpha;
            hit_first.beta = slope_beta;
            hit_first.station = Hits_In_Station[0].at(j).station;
            hit_first.usedforfit = true;
            hit_first.reference_point = 1;
            SegmentInTripletStation.push_back(hit_first);
            flag_hit[0] = 1;
            //Fill 3rd hit info
            ev5_Segment hit_last;
            hit_last.deltaphi = phi[2];
            hit_last.z = Hits_In_Station[2].at(k).z;
            hit_last.alpha = slope_alpha;
            hit_last.beta = slope_beta;
            hit_last.station = Hits_In_Station[2].at(k).station;
            hit_last.usedforfit = true;
            SegmentInTripletStation.push_back(hit_last);
            hit_candidates.push_back(Hits_In_Station[2].at(k));
            flag_hit[2] = 1;

            //Check hits in 1st station, if they configure the segment
            //For 1st station
            //std::cout<<" 1st station "<<std::endl;
            for(int l=0; l<(int)nCHsInStn_1; l++){
              if(l==j) continue;
              double middle_Hit[3] = {Hits_In_Station[0].at(l).x, Hits_In_Station[0].at(l).y, Hits_In_Station[0].at(l).z};
              double DeltaPhi = ev5_DeltaPhi(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1]);
              double sign = ev5_ParticleDirection(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1] );
              //std::cout<<"Phi 1 = "<<DeltaPhi<<std::endl;
              phi[1] = sign*DeltaPhi;
              z[1] = middle_Hit[2];
              ev5_Segment hit1;
              hit1.deltaphi = phi[1];
              hit1.z = Hits_In_Station[0].at(l).z;
              hit1.alpha = slope_alpha;
              hit1.beta = slope_beta;
              hit1.station = Hits_In_Station[0].at(l).station;
              hit1.usedforfit = false;
              //std::cout<<"Phi/z = "<<phi[1]<<"/"<<z[1]<<std::endl;

              double residual_phi = ev5_ResidualDeltaPhi(slope_alpha, slope_beta, first_Hit[2], phi[0], middle_Hit[2], phi[1]);
              //std::cout<<"1 residual_phi = "<<residual_phi<<std::endl;
              if(residual_phi < thre_residual){
                hit_candidates.push_back(Hits_In_Station[0].at(l));
                SegmentInTripletStation.push_back(hit1);
                flag_hit[0] = 1;
              }
            }

            //Check hits in 2nd station, if they configure the segment
            //For 2nd station
            //std::cout<<" 2nd station "<<std::endl;
            for(int l=0; l<(int)nCHsInStn_2; l++){
              double middle_Hit[3] = {Hits_In_Station[1].at(l).x, Hits_In_Station[1].at(l).y, Hits_In_Station[1].at(l).z};
              double DeltaPhi = ev5_DeltaPhi(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1]);
              double sign = ev5_ParticleDirection(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1] );
              //std::cout<<"Phi 2 = "<<DeltaPhi<<std::endl;
              phi[1] = sign*DeltaPhi;
              z[1] = middle_Hit[2];
              ev5_Segment hit2;
              hit2.deltaphi = phi[1];
              hit2.z = Hits_In_Station[1].at(l).z;
              hit2.alpha = slope_alpha;
              hit2.beta = slope_beta;
              hit2.station = Hits_In_Station[1].at(l).station;
              hit2.usedforfit = false;
              //std::cout<<"Phi/z = "<<phi[1]<<"/"<<z[1]<<std::endl;

              double residual_phi = ev5_ResidualDeltaPhi(slope_alpha, slope_beta, first_Hit[2], phi[0], middle_Hit[2], phi[1]);
              //std::cout<<"2 residual_phi = "<<residual_phi<<std::endl;
              if(residual_phi < thre_residual){
                hit_candidates.push_back(Hits_In_Station[1].at(l));
                SegmentInTripletStation.push_back(hit2);
                flag_hit[1] = 1;
              }
            }

            //Check hits in 3rd station, if they configure the segment
            //For 3rd station
            //std::cout<<" 3rd station "<<std::endl;
            for(int l=0; l<(int)nCHsInStn_3; l++){
              if(l==k) continue;
              double middle_Hit[3] = {Hits_In_Station[2].at(l).x, Hits_In_Station[2].at(l).y, Hits_In_Station[2].at(l).z};
              double DeltaPhi = ev5_DeltaPhi(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1]);
              double sign = ev5_ParticleDirection(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1] );
              //std::cout<<"Phi 3 = "<<DeltaPhi<<std::endl;
              phi[1] = sign*DeltaPhi;
              z[1] = middle_Hit[2];
              ev5_Segment hit3;
              hit3.deltaphi = phi[1];
              hit3.z = Hits_In_Station[2].at(l).z;
              hit3.alpha = slope_alpha;
              hit3.beta = slope_beta;
              hit3.station = Hits_In_Station[2].at(l).station;
              hit3.usedforfit = false;
              //std::cout<<"Phi/z = "<<phi[1]<<"/"<<z[1]<<std::endl;

              double residual_phi = ev5_ResidualDeltaPhi(slope_alpha, slope_beta, first_Hit[2], phi[0], middle_Hit[2], phi[1]);
              //std::cout<<"3 residual_phi = "<<residual_phi<<std::endl;
              if(residual_phi < thre_residual){
                hit_candidates.push_back(Hits_In_Station[2].at(l));
                SegmentInTripletStation.push_back(hit3);
                flag_hit[2] = 1;
              }
            }

            //If segments is founded and have enough hits, save the segment info
            int hitInSegment = (int)hit_candidates.size();
            //std::cout<<"Total hit in slope = "<<hitInSegment<<std::endl;
            if(hitInSegment < 5) continue;//at least ComboHits >= 5 in the segment
            if(flag_hit[0] == 0 or flag_hit[1] == 0 or flag_hit[2] == 0) continue;//segmnet should have at least ComboHits >= 1[/station] and at least 3 consecutive stations
            bool segment_quality = ev5_TripletQuality(SegmentInTripletStation);//at least 3 consecutive stations is required and good segemnet quality
            if(segment_quality == 0) continue;//Good quality = 1, Bad quality = 0
            segment_candidates.push_back(hit_candidates);
            diag_segment_candidates.push_back(SegmentInTripletStation);
          }
        }//end segment candidate search

      if(segment_candidates.size() == 0) continue;
      if(diag_segment_candidates.size() == 0) continue;

      //---------------------------------------------------------------------
      // Select several best candidates in 3 consecutive stations
      //---------------------------------------------------------------------
      std::vector<std::vector<ev5_HitsInNthStation>> ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>> ThisIsBestSegment_Diag;
      ThisIsBestSegment.clear();
      ThisIsBestSegment_Diag.clear();
      // remove duplicate segments in 3 station based on hitID and fit each slope
      ev5_select_best_segments_step_01(HitsInTimeCluster, segment_candidates, diag_segment_candidates, ThisIsBestSegment, ThisIsBestSegment_Diag, nHitInStation, thre_residual);
      // collect remaining hits in 3 station and add those hits to the segment
      ev5_select_best_segments_step_02(n, HitsInTimeCluster, ThisIsBestSegment, ThisIsBestSegment_Diag, nHitInStation, thre_residual);
      // extend the slope to neighboring and then, add hits to the segment and fit the segment
      ev5_select_best_segments_step_03(HitsInTimeCluster, ThisIsBestSegment, ThisIsBestSegment_Diag, nHitInStation, thre_residual);
      //if Segments >= 2, remove duplicate segments based on hitID
      ev5_select_best_segments_step_04(HitsInTimeCluster, ThisIsBestSegment, ThisIsBestSegment_Diag, nHitInStation, thre_residual);
      //if Segments >= 2, remove duplicate segments based on slope value and fraction of overlapped hits
      ev5_select_best_segments_step_05(HitsInTimeCluster, ThisIsBestSegment, ThisIsBestSegment_Diag, nHitInStation, thre_residual);


      //-------------------------------
      // Fill only one best candidate
      //-------------------------------
      for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
        best_triplet_segments.push_back(ThisIsBestSegment.at(j));
        diag_best_triplet_segments.push_back(ThisIsBestSegment_Diag.at(j));
      }
      for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
        all_ThisIsBestSegment.push_back(ThisIsBestSegment.at(j));
        all_ThisIsBestSegment_Diag.push_back(ThisIsBestSegment_Diag.at(j));
      }
      } // loop end on stations

      // remove duplicate segments based on hitID. remove segments based on based on slope value, fraction of overlapped hits, and Chi2/NDF
      ev5_select_best_segments_step_06(all_ThisIsBestSegment, all_ThisIsBestSegment_Diag, thre_residual);
      // merge segments based on slope value, phi range. remove overlapped hits in the segments based on hitID
      int number_of_merged_segments = 0;
      ev5_select_best_segments_step_07(all_ThisIsBestSegment, all_ThisIsBestSegment_Diag, thre_residual, number_of_merged_segments);
      //std::cout<<"number_of_merged_segments = "<<number_of_merged_segments<<std::endl;


      //fill
      best_triplet_segments = all_ThisIsBestSegment;
      diag_best_triplet_segments = all_ThisIsBestSegment_Diag;

}//end ev5_SegmentSearchInTriplet


//-----------------------------------------------------------------------------
    double PhiZSeedFinder::ev5_DeltaPhi(double x1, double y1, double x2, double y2){

      double a[2] = {x1, y1};
      double b[2] = {x2, y2};
      double c = (a[0]*b[0]+a[1]*b[1])/(sqrt(a[0]*a[0]+a[1]*a[1])*sqrt(b[0]*b[0]+b[1]*b[1]));
      double deltaphi = acos(c);//[rad]

      return deltaphi;

    }

//--------------------------------------------------------------------------------//
    double PhiZSeedFinder::ev5_ParticleDirection(double x1, double y1, double x2, double y2){

      //cross product: check if particle is left-hand or right-hand
      double a[2] = {x1, y1};
      double b[2] = {x2, y2};
      double cross_product = a[0]*b[1] - a[1]*b[0];
      int sign = 0;
      if(cross_product < 0) sign = -1;//negative: non CE particle, 2Pi[rad] -> O[rad]
      else sign = 1;//positive: CE like particle, O[rad] -> 2Pi[rad]
      //std::cout<<"sign = "<<sign<<std::endl;

       return sign;
    }

//--------------------------------------------------------------------------------//
  double PhiZSeedFinder::ev5_ResidualDeltaPhi(double alpha, double beta, double z1, double phi1, double z2, double phi2){

    double y_guess = 0.0;
    //y_guess = alpha*fabs(z2 - z1) + phi1;//Y = aX + b
    y_guess = alpha*z2 + beta;//Y = aX + b
    //if(z1 > z2) y_guess = -y_guess;
    double residual = fabs(phi2 - y_guess);
    //std::cout<<"Y_guess = "<<y_guess<<std::endl;
    //std::cout<<"residual = "<<residual<<std::endl;
    //if(residual > 1.0 or residual < -1.0) std::cout<<"ev5_ResidualDeltaPhi_kanben "<<" "<< event <<std::endl;
    //ev5_residual.push_back(phi2 - y_guess);

    return fabs(residual);
  }


//--------------------------------------------------------------------------------//
  bool PhiZSeedFinder::ev5_TripletQuality(const std::vector<ev5_Segment> diag_hit){

    std::vector<ev5_Segment> hit = diag_hit;


    std::vector<int> nstations;
    for(int i=0; i<(int)hit.size(); i++){
      int station = hit[i].station;
      nstations.push_back(station);
    }
    //std::cout << "before removed nstations.size() = " << (int)nstations.size()<<std::endl;
    //Sort the vector in increasing order
    std::sort(nstations.begin(), nstations.end());
    //Remove duplicates
    nstations.erase(std::unique(nstations.begin(), nstations.end()), nstations.end());
    //std::cout << "after removed nstations.size() = " << (int)nstations.size()<<std::endl;

    // re-sort vector in ascending order of z-coordinate
    std::sort(hit.begin(), hit.end(), [](const ev5_Segment& a, const ev5_Segment& b){
     return a.z < b.z;
    });

    /*for(int i=0; i<(int)hit.size(); i++){
      int station = hit[i].station;
      double z = hit[i].z;
      double deltaphi = hit[i].deltaphi;
      //std::cout<<"station/z/deltaphi = "<<station<<"/"<<z<<"/"<<deltaphi<<std::endl;
    }*/


    std::vector<double> DeltaPhi;
    for(int i=0; i<(int)nstations.size()-1; i++){
      double phi[2] = {-9999.9, -9999.9};
      for(int j=0; j<(int)hit.size(); j++){
        //std::cout<<"nstations[i]/hit[j].station/hit[j].deltaphi = "<<nstations[i]<<"/"<<hit[j].station<<"/"<<hit[j].deltaphi<<std::endl;
        if(nstations[i] == hit[j].station) phi[0] = hit[j].deltaphi;
      }
      //std::cout<<"phi[0]/phi[1] = "<<phi[0]<<"/"<<phi[1]<<std::endl;
      for(int j=0; j<(int)hit.size(); j++){
        //std::cout<<"nstations[i+1]/hit[j].station = "<<nstations[i+1]<<"/"<<hit[j].station<<"/"<<hit[j].deltaphi<<std::endl;
        if(nstations[i+1] == hit[j].station){
          phi[1] = hit[j].deltaphi;
          break;
        }
      }
      if(phi[0] < -900. or phi[1] < -900.) std::cout<<"gomi_desuyo"<<std::endl;
      //std::cout<<"phi[0]/phi[1] = "<<phi[0]<<"/"<<phi[1]<<std::endl;
      double diff_phi = 999.9;
      if(phi[0]*phi[1] < 0){//either phi is positive or negative
        phi[0] = fabs(phi[0]);
        phi[1] = fabs(phi[1]);
        diff_phi = phi[1] + phi[0];
      }
      else{//phi1 = 0 and phi2 = +(-), or phi1 = +(-) and phi2 = +(-)
        diff_phi = fabs(phi[0] - phi[1]);
      }

      DeltaPhi.push_back(diff_phi);
      //std::cout<<"diff_phi = "<<diff_phi<<std::endl;
      if(diff_phi > 10) std::cout<<"gomi_desuyo"<<std::endl;
  }


    bool QualityIsGood = 1;//Good = 1, Bad = 0
    for(int i=0; i<(int)DeltaPhi.size(); i++){
      for(int j=0; j<(int)DeltaPhi.size(); j++){
        if(i==j) continue;

        //std::cout<<"DeltaPhi[i]/DeltaPhi[i]*2/DeltaPhi[j] = "<<DeltaPhi[i]<<"/"<<DeltaPhi[i]*2<<"/"<<DeltaPhi[j]<<std::endl;
        if(DeltaPhi[i]*2 < DeltaPhi[j] and DeltaPhi[j] > 0.4) QualityIsGood = 0;
      }
    }

    //std::cout<<"QualityIsGood = "<<QualityIsGood<<std::endl;
    return QualityIsGood;
  }

//-----------------------------------------------------------------------------
  void PhiZSeedFinder::ev5_select_best_segments_step_01(const std::vector<ev5_HitsInNthStation>& HitsInTimeCluster, const std::vector<std::vector<ev5_HitsInNthStation>>& segment_candidates, const std::vector<std::vector<ev5_Segment>>& diag_segment_candidates, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi){

    //std::cout << "-----------------------------------" << std::endl;
    //std::cout << "-----------------------------------" << std::endl;
    //std::cout << " ev5_select_best_segments_step_01     " << std::endl;
    //std::cout << "-----------------------------------" << std::endl;
    //std::cout << "-----------------------------------" << std::endl;
    //-----------------------------------------------------------------
    //    Remove same combination in the same 3 consecutive stations
    //-----------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
    //  (1). Check all hitIDs in each segments and remove segments if there is duplicate
    //---------------------------------------------------------------------------------------------------------
    std::vector<std::vector<int>> hitID_list;
    hitID_list.clear();
    for(int i=0; i<(int)segment_candidates.size(); i++){
      std::vector<int> temp_hitID_list;
      temp_hitID_list.clear();
      for(int j=0; j<(int)segment_candidates.at(i).size(); j++){
        temp_hitID_list.push_back(segment_candidates.at(i).at(j).hitID);
      }
      //Sort temp_hitID_list in increasing order
      sort(temp_hitID_list.begin(), temp_hitID_list.end());
      hitID_list.push_back(temp_hitID_list);
    }
    // Remove duplicates
    std::sort(hitID_list.begin(), hitID_list.end());
    hitID_list.erase(std::unique(hitID_list.begin(), hitID_list.end()), hitID_list.end());

    //---------------------------------------------------------------------------------------------------------
    // (2). Delete if there is overlap: example, 0-1-2 and 0-1-2-4, 1-2-4 and 0-1-2-4, in this case 0-1-2-4 will remain
    //---------------------------------------------------------------------------------------------------------
    std::vector<int> delete_index;
    delete_index.clear();
    for (int i = 0; i < (int)hitID_list.size(); i++) {
      //std::cout << "i = " << i << std::endl;
      for (int j = 0; j <(int)hitID_list.size(); j++) {
        if(i == j) continue;
        if(hitID_list.at(j) == hitID_list.at(i)){
          //std::cout << "found j = " << j << std::endl;
          continue;
        }
        int size = 0;
        for (int k = 0; k <(int)hitID_list.at(j).size(); k++) {
          for (int l = 0; l <(int)hitID_list.at(i).size(); l++) {
            if(hitID_list.at(i).at(l) == hitID_list.at(j).at(k)) size++;
          }
        }
        if(size == (int)hitID_list.at(i).size()) delete_index.push_back(i);
      }
    }

    //  make a new hitID after removing the overlap event
    std::vector<std::vector<int>> new_hitID_list;
    for(int i = 0; i <(int)hitID_list.size(); i++){
      int go = 1;
      for(int j=0; j<(int)delete_index.size(); j++){
        if(delete_index.at(j) == i) go = 0;
      }
      if(go == 1) new_hitID_list.push_back(hitID_list.at(i));
    }

    //---------------------------------------------------------------------------------------------------------
    // (3)  find an endex in "new hitID" correspond to the segment_candidates
    //---------------------------------------------------------------------------------------------------------
    std::vector<int> find_index;
    find_index.clear();
    for (int i = 0; i < (int)new_hitID_list.size(); i++) {
      //std::cout << "i = " << i << std::endl;
        bool flag = 0;
        int index_for_segment = 0;
        for(int j=0; j<(int)segment_candidates.size(); j++){
          int count = 0;
          if(new_hitID_list.at(i).size() != segment_candidates.at(j).size()) continue;
          for (int k = 0; k<(int)new_hitID_list.at(i).size(); k++) {
              for(int l=0; l<(int)segment_candidates.at(j).size(); l++){
                if(new_hitID_list.at(i).at(k) == segment_candidates.at(j).at(l).hitID) count++;
              }
          }
          if(count == (int)new_hitID_list.at(i).size()){
            flag = 1;
            index_for_segment = j;
            break;
          }
        }
        if(flag == 1) find_index.push_back(index_for_segment);
    }

    //---------------------------------------------------------------------------------------------------------
    // Fit segments
    //---------------------------------------------------------------------------------------------------------
    int nSegmentInTriplet = (int)diag_segment_candidates.size();
    //select the best candidate
    for(int i=0; i<nSegmentInTriplet; i++){

      bool go = 0;
      for(int j=0; j<(int)find_index.size(); j++){
        if(find_index[j] == i) go = 1;
      }
      if(go != 1) continue;

      //Plot graph: Phi vs. Z of each segment
      TGraph *graph = new TGraph();
      int p = 0;
      int nthStation = 0;
      double z_max = diag_segment_candidates.at(i).at(0).z;
      double z_min = diag_segment_candidates.at(i).at(0).z;
      double m = diag_segment_candidates.at(i).at(0).alpha;
      int nComboHitsSegment = (int)diag_segment_candidates.at(i).size();
      for(int j=0; j<nComboHitsSegment; j++){
        // z as x-axis and phi as y-axis
        double z = diag_segment_candidates.at(i).at(j).z;
        double phi = diag_segment_candidates.at(i).at(j).deltaphi;
        graph->SetPoint(p++, z, phi);

        if(z_min > z) z_min = z;
        if(z_max < z) z_max = z;
        if(nthStation < diag_segment_candidates.at(i).at(j).station) nthStation = diag_segment_candidates.at(i).at(j).station;
      }

      //Fitting function: f(x) = a*x + b
      TF1* fitFunc = new TF1("linearFit", "[0]*x + [1]", graph->GetXaxis()->GetXmin(), graph->GetXaxis()->GetXmax());
      fitFunc->SetParameter(0, m);
      fitFunc->SetParameter(1, 0.0);
      fitFunc->SetLineColor(kBlue);
      fitFunc->SetLineStyle(2);
      graph->Fit(fitFunc, "R");

      //Print fitting result
      //std::cout << "   Slope (a): " << fitFunc->GetParameter(0) << " +/- " << fitFunc->GetParError(0) << std::endl;
      //std::cout << "   Intercept (b): " << fitFunc->GetParameter(1) << " +/- " << fitFunc->GetParError(1) << std::endl;
      //std::cout << "   Chi-square: " << fitFunc->GetChisquare() << std::endl;
      //std::cout << "   NDF: " << fitFunc->GetNDF() << std::endl;



      TCanvas *canvas = new TCanvas("canvas", "My TGraph", 800, 600);
      canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
      graph->SetMarkerColor(kRed);
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(0.5);
      graph->Draw("AP");
      graph->GetXaxis()->SetTitle("Z [mm]");
      graph->GetYaxis()->SetTitle("#Delta #Phi [rad]");
      // Add title at the top
      TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
      title->SetFillColor(0);
      title->SetTextAlign(22);
      title->Draw("same");
      // Set the range
      graph->GetXaxis()->SetLimits(z_min-200, z_max+200);
      graph->GetYaxis()->SetRangeUser(-2.0, 2.0);

      // Draw dphi/dX in text
      TLatex dphidx;
      dphidx.SetTextSize(0.03);
      dphidx.SetTextAlign(13);
      int ntot = nCH;
      dphidx.DrawLatexNDC(0.15, 0.88, Form("Station (%d - %d - %d), CHs = %d", nthStation-2, nthStation-1, nthStation, ntot));
      dphidx.DrawLatexNDC(0.15, 0.85, Form("#frac{d#phi}{dZ} = %f, fitted with 2 CHs #candidate %d", m, i));
      dphidx.DrawLatexNDC(0.15, 0.78, Form("CHs in Segment = %d", (int)diag_segment_candidates.at(i).size()));
      dphidx.DrawLatexNDC(0.15, 0.73, Form("#Delta #phi = %.1f", threshold_deltaphi));
      dphidx.DrawLatexNDC(0.15, 0.70, "Fit Results: ");
      dphidx.DrawLatexNDC(0.15, 0.67, Form("Slope (a): %f +/- %f", fitFunc->GetParameter(0), fitFunc->GetParError(0)));
      dphidx.DrawLatexNDC(0.15, 0.64, Form("Intercept (b): %f +/- %f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
      dphidx.DrawLatexNDC(0.15, 0.61, Form("Chi-square: %f", fitFunc->GetChisquare()));
      dphidx.DrawLatexNDC(0.15, 0.58, Form("NDF: %d", fitFunc->GetNDF()));

      //select best candidate
      //std::cout<<"(fitFunc->GetChisquare()/fitFunc->GetNDF())*100.0 = "<<(fitFunc->GetChisquare()/fitFunc->GetNDF())*100.0<<std::endl;
      double chiNDF = fitFunc->GetChisquare()/fitFunc->GetNDF();

      //push back segment
      ThisIsBestSegment.push_back(segment_candidates.at(i));
      ThisIsBestSegment_Diag.push_back(diag_segment_candidates.at(i));
      int index = (int)ThisIsBestSegment_Diag.size()-1;
      for(int j=0; j<(int)ThisIsBestSegment_Diag.at(index).size(); j++){
        ThisIsBestSegment_Diag.at(index).at(j).alpha = fitFunc->GetParameter(0);
        ThisIsBestSegment_Diag.at(index).at(j).beta = fitFunc->GetParameter(1);
        ThisIsBestSegment_Diag.at(index).at(j).chiNDF = chiNDF*100.0;
      }


      //delete
      delete canvas;
      delete graph;
      delete title;
      //Plot end

    }//end nSegmentInTriplet



  }//end ev5_select_best_segments_step_01

//-----------------------------------------------------------------------------
  void PhiZSeedFinder::ev5_select_best_segments_step_02(int station, const std::vector<ev5_HitsInNthStation>& HitsInTimeCluster, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi){

      std::vector<std::vector<ev5_HitsInNthStation>> segments = ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>>& diag_segments = ThisIsBestSegment_Diag;

      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << " ev5_select_best_segments_step_02  " << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;

      //---------------------------------------------------------------------------------------------
      // Search remaining ComboHit candidates in 3 consecutive stations
      //---------------------------------------------------------------------------------------------
      //obtain the minimum station in the segment
      int min_station = station;
      //std::cout << "min_station = " <<min_station<< std::endl;

      //Segment candidates
      std::vector<std::vector<ev5_HitsInNthStation>> new_segments;
      new_segments.clear();
      std::vector<std::vector<ev5_Segment>>  new_diag_segments;
      new_diag_segments.clear();

      //Take segment
      for(int i=0; i<(int)segments.size(); i++){

        //Fill all hit info of "i-th" segment
        std::vector<ev5_HitsInNthStation> hit_candidates;
        std::vector<ev5_Segment> hit_diag_candidates;
        hit_candidates.clear();
        hit_diag_candidates.clear();
        for(int j=0; j<(int)segments.at(i).size(); j++){
          hit_candidates.push_back(segments.at(i).at(j));
          hit_diag_candidates.push_back(diag_segments.at(i).at(j));
        }

        double slope_alpha = hit_diag_candidates.at(0).alpha;//Slope (a)
        double slope_beta = hit_diag_candidates.at(0).beta;//Intercept (b)
        //std::cout<<"slope_alpha/beta = "<<slope_alpha<<"/"<<slope_beta<<std::endl;

        // Take remaining ComboHits in 3 consecutive stations
        for(int k=min_station; k<=min_station+2; k++){

        std::vector<ev5_HitsInNthStation> Hits_In_Station;
        Hits_In_Station.clear();
        for(int j=0; j<(int)HitsInTimeCluster.size(); j++){
            ev5_HitsInNthStation hitsin_nthstation;
            hitsin_nthstation.phi          = HitsInTimeCluster.at(j).phi;
            hitsin_nthstation.strawhits    = HitsInTimeCluster.at(j).station;
            hitsin_nthstation.x            = HitsInTimeCluster.at(j).x;
            hitsin_nthstation.y            = HitsInTimeCluster.at(j).y;
            hitsin_nthstation.z            = HitsInTimeCluster.at(j).z;
            hitsin_nthstation.station      = HitsInTimeCluster.at(j).station;
            hitsin_nthstation.plane        = HitsInTimeCluster.at(j).plane;
            hitsin_nthstation.face        = HitsInTimeCluster.at(j).face;
            hitsin_nthstation.panel        = HitsInTimeCluster.at(j).panel;
            hitsin_nthstation.hitID        = HitsInTimeCluster.at(j).hitID;
            if(k != HitsInTimeCluster.at(j).station) continue;
            int flag_alreadyUsed = 0;
            for (const auto &element : hit_candidates) {
              if(element.hitID == HitsInTimeCluster.at(j).hitID) flag_alreadyUsed = 1;
            }
            if(flag_alreadyUsed == 1) continue;
            Hits_In_Station.push_back(hitsin_nthstation);
        }
        size_t nCHsInStn_1 = Hits_In_Station.size();
        if(!((int)nCHsInStn_1 >= 1)) break;//at leat 1 ComboHits

          /*for(size_t p=0; p<Hits_In_Station.size(); p++){
            std::cout<<"phi/nStrawHits/station/plane/face/panel/x/y/z = "<<Hits_In_Station.at(p).phi<<"/"<<Hits_In_Station.at(p).strawhits<<"/"<<Hits_In_Station.at(p).station<<"/"<<Hits_In_Station.at(p).plane<<"/"<<Hits_In_Station.at(p).face<<"/"<<Hits_In_Station.at(p).panel<<"/"<<Hits_In_Station.at(p).x<<"/"<<Hits_In_Station.at(p).y<<"/"<<Hits_In_Station.at(p).z<<std::endl;
          }*/


            // reference point: deltaphi = 0, first_Hit[4] = {x, y, z, phi}
            double first_Hit[4] = {9999.9, 9999.9, 9999.9, 9999.9};
            double reference_phi = -999.9;
            for(int p=0;p<(int)hit_diag_candidates.size(); p++){
              if(hit_diag_candidates.at(p).reference_point == 1){
                first_Hit[0] = hit_candidates[p].x;
                first_Hit[1] = hit_candidates[p].y;
                first_Hit[2] = hit_candidates[p].z;
                first_Hit[3] = hit_candidates[p].phi;
              }
              reference_phi = hit_diag_candidates[p].deltaphi;
            }
            //if(first_Hit[0] > 9999. or first_Hit[2] > 9999. or first_Hit[2] > 9999. or first_Hit[3] > 9999.) std::cout<<"first_Hit_wrong_check_please"<<std::endl;


            //Check hits in "n-th" station, if they configure the segment
            //For "n-th" station
            //std::cout<<" station = "<<k<<std::endl;
            double phi[2] = {reference_phi, 0.0};
            for(int l=0; l<(int)nCHsInStn_1; l++){
              double middle_Hit[3] = {Hits_In_Station.at(l).x, Hits_In_Station.at(l).y, Hits_In_Station.at(l).z};
              double DeltaPhi = ev5_DeltaPhi(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1]);
              double sign = ev5_ParticleDirection(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1] );
              //std::cout<<"DeltaPhi between 1st ST/3rd ST= "<<DeltaPhi<<std::endl;
              //std::cout<<"particle has positive/negative track = "<<sign<<std::endl;
              //std::cout<<"Phi 1 = "<<DeltaPhi<<std::endl;
              phi[1] = sign*DeltaPhi;
              ev5_Segment hit1;
              hit1.deltaphi = phi[1];
              hit1.z = Hits_In_Station.at(l).z;
              hit1.alpha = slope_alpha;
              hit1.station = Hits_In_Station.at(l).station;
              hit1.usedforfit = false;
              //std::cout<<"Phi/z = "<<phi[1]<<"/"<<middle_Hit[2]<<std::endl;

              double residual_phi = ev5_ResidualDeltaPhi(slope_alpha, slope_beta, first_Hit[2], phi[0], middle_Hit[2], phi[1]);
              //std::cout<<"1 residual_phi = "<<residual_phi<<std::endl;
              if(residual_phi < threshold_deltaphi){
                hit_candidates.push_back(Hits_In_Station.at(l));
                hit_diag_candidates.push_back(hit1);
              }
            }

        }//end loop on station

        new_segments.push_back(hit_candidates);
        new_diag_segments.push_back(hit_diag_candidates);

      }//end loop on segment


    //Re-Fill
    ThisIsBestSegment.clear();
    ThisIsBestSegment_Diag.clear();
    ThisIsBestSegment = new_segments;
    ThisIsBestSegment_Diag = new_diag_segments;


  }//end ev5_select_best_segments_step_02

//--------------------------------------------------------------------------------//

  void PhiZSeedFinder::ev5_select_best_segments_step_03(const std::vector<ev5_HitsInNthStation>& HitsInTimeCluster, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi){

      std::vector<std::vector<ev5_HitsInNthStation>> segments = ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>>& diag_segments = ThisIsBestSegment_Diag;

      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << " ev5_select_best_segments_step_03     " << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;

      //---------------------------------------------------------------------------------------------
      // Search surrounding ComboHit candidates in neighboring consecutive stations
      // For example: 3 consecutive station (2, 3, 4)
      // station : 0, 5, 6, 7...
      //---------------------------------------------------------------------------------------------
      //If hits found in neighboring consecutive stations, set up the flag
      int size = (int)segments.size();
      int *hit_found;
      hit_found = new int[size];
      for(int i=0; i<size; i++){
        hit_found[i] = 0;
      }

      //obtain the minimum station in the segment
      int min_station = 999;
      for(int i=0; i<(int)segments.size(); i++){
        for(int j=0; j<(int)segments.at(i).size(); j++){
          int station = segments.at(i).at(j).station;
          if(min_station > station) min_station = station;
        }
      }
      //if(min_station == 999) std::cout<<"nthstation_is_999_something_went_wrong!!"<<std::endl;


      //Segment candidates
      std::vector<std::vector<ev5_HitsInNthStation>> new_segment;
      new_segment.clear();
      std::vector<std::vector<ev5_Segment>>  new_diag_segment;
      new_diag_segment.clear();



      //Take segment
      for(int i=0; i<(int)segments.size(); i++){

        //Fill all hit info
        std::vector<ev5_HitsInNthStation> hit_candidates;
        std::vector<ev5_Segment> hit_diag_candidates;
        hit_candidates.clear();
        hit_diag_candidates.clear();
        for(int j=0; j<(int)segments.at(i).size(); j++){
            hit_candidates.push_back(segments.at(i).at(j));
            hit_diag_candidates.push_back(diag_segments.at(i).at(j));
        }

        // Take combo hits in (n-i)-th stations
        int n = min_station;
        int loop = n;
        int count = 1;
        for(int k=1; 0<=loop-k; k++){
        if(count != k) break;
        //std::cout << "-----------------------------------" << std::endl;
        //std::cout << "            (n-i)-th               " << std::endl;
        //std::cout << "-----------------------------------" << std::endl;
        //std::cout << "station = " << loop-k<< std::endl;

        double slope_alpha = 0.0;//Slope (a)
        double slope_beta = 0.0;//Intercept (b)
        double ChiNDF = 0.0;
        ev5_fit_slope(hit_diag_candidates, slope_alpha, slope_beta, ChiNDF);
        //std::cout<<"slope_alpha/beta = "<<slope_alpha<<"/"<<slope_beta<<std::endl;


        std::vector<ev5_HitsInNthStation> Hits_In_Station;
        Hits_In_Station.clear();
        for(int j=0; j<(int)HitsInTimeCluster.size(); j++){
            ev5_HitsInNthStation hitsin_nthstation;
            hitsin_nthstation.phi          = HitsInTimeCluster.at(j).phi;
            hitsin_nthstation.strawhits    = HitsInTimeCluster.at(j).station;
            hitsin_nthstation.x            = HitsInTimeCluster.at(j).x;
            hitsin_nthstation.y            = HitsInTimeCluster.at(j).y;
            hitsin_nthstation.z            = HitsInTimeCluster.at(j).z;
            hitsin_nthstation.station      = HitsInTimeCluster.at(j).station;
            hitsin_nthstation.plane        = HitsInTimeCluster.at(j).plane;
            hitsin_nthstation.face        = HitsInTimeCluster.at(j).face;
            hitsin_nthstation.panel        = HitsInTimeCluster.at(j).panel;
            hitsin_nthstation.hitID        = HitsInTimeCluster.at(j).hitID;
            if(loop-k == HitsInTimeCluster.at(j).station) Hits_In_Station.push_back(hitsin_nthstation);
        }
          size_t nCHsInStn_1 = Hits_In_Station.size();
          if(!((int)nCHsInStn_1 >= 1)) break;

          /*for(size_t p=0; p<Hits_In_Station.size(); p++){
            std::cout<<"phi/nStrawHits/station/plane/face/panel/x/y/z = "<<Hits_In_Station.at(p).phi<<"/"<<Hits_In_Station.at(p).strawhits<<"/"<<Hits_In_Station.at(p).station<<"/"<<Hits_In_Station.at(p).plane<<"/"<<Hits_In_Station.at(p).face<<"/"<<Hits_In_Station.at(p).panel<<"/"<<Hits_In_Station.at(p).x<<"/"<<Hits_In_Station.at(p).y<<"/"<<Hits_In_Station.at(p).z<<std::endl;
          }*/

          std::vector<ev5_HitsInNthStation> LeftHit = segments.at(i);
          // re-sort vector in ascending order of z-coordinate
          std::sort(LeftHit.begin(), LeftHit.end(), [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b){
            return a.z < b.z;
          });
          //Check the left most hits in the segment

          double first_Hit[4] = {LeftHit.at(0).x, LeftHit.at(0).y, LeftHit.at(0).z, LeftHit.at(0).phi};
          //Check hits in (n-i)th station, if they configure the segment
          //For (n-i)th station
          int hit_yes = 0;
           for(int l=0; l<(int)nCHsInStn_1; l++){
              double middle_Hit[3] = {Hits_In_Station.at(l).x, Hits_In_Station.at(l).y, Hits_In_Station.at(l).z};
              double DeltaPhi = ev5_DeltaPhi(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1]);
              double sign = ev5_ParticleDirection(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1] );
              //Calculate the slope between 1st and 3rd hit
              double phi[2] = {0.0, sign*DeltaPhi};
              ev5_Segment hit;
              hit.deltaphi = phi[1];
              hit.z = Hits_In_Station.at(l).z;
              hit.alpha = slope_alpha;
              hit.station = Hits_In_Station.at(l).station;
              hit.usedforfit = false;

              double residual_phi = ev5_ResidualDeltaPhi(slope_alpha, slope_beta, first_Hit[2], phi[0], middle_Hit[2], phi[1]);
              //std::cout<<"residual_phi = "<<residual_phi<<std::endl;
              double  thre_residual = 0.4;
              if(residual_phi < thre_residual){
              //std::cout<<"found!!!!!"<<std::endl;
                hit_found[i]++;
                hit_yes = 1;
                hit_candidates.push_back(Hits_In_Station.at(l));
                hit_diag_candidates.push_back(hit);
              }
           }
          if(hit_yes == 1) count++;

        }//end (n-i)-th stations


        // Take combo hits in (n+i)-th stations
        loop = n+2;
        int nstation = 18;
        count = 1;
        for(int k=1; loop+k < nstation; k++){
        if(count != k) break;
        //std::cout << "-----------------------------------" << std::endl;
        //std::cout << "            (n+i)-th               " << std::endl;
        //std::cout << "-----------------------------------" << std::endl;
        //std::cout << "station = " << loop+k  << std::endl;
        double slope_alpha = 0.0;//Slope (a)
        double slope_beta = 0.0;//Intercept (b)
        double ChiNDF = 0.0;
        ev5_fit_slope(hit_diag_candidates, slope_alpha, slope_beta, ChiNDF);
        //std::cout<<"slope_alpha/beta = "<<slope_alpha<<"/"<<slope_beta<<std::endl;

        std::vector<ev5_HitsInNthStation> Hits_In_Station;
        Hits_In_Station.clear();
        for(int j=0; j<(int)HitsInTimeCluster.size(); j++){
            ev5_HitsInNthStation hitsin_nthstation;
            hitsin_nthstation.phi          = HitsInTimeCluster.at(j).phi;
            hitsin_nthstation.strawhits    = HitsInTimeCluster.at(j).station;
            hitsin_nthstation.x            = HitsInTimeCluster.at(j).x;
            hitsin_nthstation.y            = HitsInTimeCluster.at(j).y;
            hitsin_nthstation.z            = HitsInTimeCluster.at(j).z;
            hitsin_nthstation.station      = HitsInTimeCluster.at(j).station;
            hitsin_nthstation.plane        = HitsInTimeCluster.at(j).plane;
            hitsin_nthstation.face        = HitsInTimeCluster.at(j).face;
            hitsin_nthstation.panel        = HitsInTimeCluster.at(j).panel;
            hitsin_nthstation.hitID        = HitsInTimeCluster.at(j).hitID;
            if(loop+k == HitsInTimeCluster.at(j).station) Hits_In_Station.push_back(hitsin_nthstation);
        }

          size_t nCHsInStn_1 = Hits_In_Station.size();
          if(!((int)nCHsInStn_1 >= 1)) break;
          //std::cout << "loop/k " << loop<<"/"<<k<< std::endl;

          /*for(size_t p=0 ; p<Hits_In_Station.size(); p++){
            std::cout<<"phi/nStrawHits/station/plane/face/panel/x/y/z = "<<Hits_In_Station.at(p).phi<<"/"<<Hits_In_Station.at(p).strawhits<<"/"<<Hits_In_Station.at(p).station<<"/"<<Hits_In_Station.at(p).plane<<"/"<<Hits_In_Station.at(p).face<<"/"<<Hits_In_Station.at(p).panel<<"/"<<Hits_In_Station.at(p).x<<"/"<<Hits_In_Station.at(p).y<<"/"<<Hits_In_Station.at(p).z<<std::endl;
          }*/

          std::vector<ev5_HitsInNthStation> RightHit = segments.at(i);
          // re-sort vector in ascending order of z-coordinate
          std::sort(RightHit.begin(), RightHit.end(), [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b){
            return a.z < b.z;
          });
          //Check the left most hits in the segment

          //int index = (int)RightHit.size()-1;
          //double first_Hit[4] = {RightHit.at(index).x, RightHit.at(index).y, RightHit.at(index).z, RightHit.at(index).phi};
          double first_Hit[4] = {RightHit.at(0).x, RightHit.at(0).y, RightHit.at(0).z, RightHit.at(0).phi};
          //Check hits in (n-i)th station, if they configure the segment
          int hit_yes = 0;
           for(int l=0; l<(int)nCHsInStn_1; l++){
              double middle_Hit[3] = {Hits_In_Station.at(l).x, Hits_In_Station.at(l).y, Hits_In_Station.at(l).z};
              double DeltaPhi = ev5_DeltaPhi(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1]);
              double sign = ev5_ParticleDirection(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1] );
              //Calculate the slope between 1st and 3rd hit
              double phi[2] = {0.0, sign*DeltaPhi};
              ev5_Segment hit;
              hit.deltaphi = phi[1];
              hit.z = Hits_In_Station.at(l).z;
              hit.alpha = slope_alpha;
              hit.station = Hits_In_Station.at(l).station;
              hit.usedforfit = false;

              double residual_phi = ev5_ResidualDeltaPhi(slope_alpha, slope_beta, first_Hit[2], phi[0], middle_Hit[2], phi[1]);
              //std::cout<<"residual_phi = "<<residual_phi<<std::endl;
              double  thre_residual = 0.4;
              if(residual_phi < thre_residual){
                hit_found[i]++;
                hit_yes = 1;
                hit_candidates.push_back(Hits_In_Station.at(l));
                hit_diag_candidates.push_back(hit);
              }
           }
          if(hit_yes == 1) count++;

        }//end (n+i)-th station

          double slope_alpha = 0.0;//Slope (a)
          double slope_beta = 0.0;//Intercept (b)
          double ChiNDF = 0.0;
          ev5_fit_slope(hit_diag_candidates, slope_alpha, slope_beta, ChiNDF);
          for(int j=0; j<(int)hit_diag_candidates.size(); j++){
            hit_diag_candidates[j].alpha = slope_alpha;
            hit_diag_candidates[j].chiNDF = ChiNDF;
          }
          new_segment.push_back(hit_candidates);
          new_diag_segment.push_back(hit_diag_candidates);

      }//end segment loop



    //Re-Fill
    ThisIsBestSegment.clear();
    ThisIsBestSegment_Diag.clear();
    ThisIsBestSegment = new_segment;
    ThisIsBestSegment_Diag = new_diag_segment;


  }//end ev5_select_best_segments_step_03

//--------------------------------------------------------------------------------//
    void PhiZSeedFinder::ev5_fit_slope(const std::vector<ev5_Segment>& hit_diag, double& alpha, double& beta, double& chindf){

      //Plot graph: Phi vs. Z of each segment
      TGraph *graph = new TGraph();
      int p = 0;
      int nthStation = 0;
      double z_max = -9999.;
      double z_min = 9999.;
      int nComboHitsSegment = (int)hit_diag.size();
      for(int j=0; j<nComboHitsSegment; j++){
        // z as x-axis and phi as y-axis
        double z = hit_diag.at(j).z;
        double phi = hit_diag.at(j).deltaphi;
        graph->SetPoint(p++, z, phi);
        if(z_min > z) z_min = z;
        if(z_max < z) z_max = z;
        if(nthStation < hit_diag.at(j).station) nthStation = hit_diag.at(j).station;
      }

      //Fitting function: f(x) = a*x + b
      TF1* fitFunc = new TF1("linearFit", "[0]*x + [1]", graph->GetXaxis()->GetXmin(), graph->GetXaxis()->GetXmax());
      double slope_alpha = hit_diag.at(0).alpha;
      double slope_beta = 0.0;
      fitFunc->SetParameter(0, slope_alpha);
      fitFunc->SetParameter(1, slope_beta);
      fitFunc->SetLineColor(kBlue);
      fitFunc->SetLineStyle(2);
      graph->Fit(fitFunc, "R");

      //Print fitting result
      /*std::cout << "Fit Results:" << std::endl;
      std::cout << "Slope (a): " << fitFunc->GetParameter(0) << " +/- " << fitFunc->GetParError(0) << std::endl;
      std::cout << "Intercept (b): " << fitFunc->GetParameter(1) << " +/- " << fitFunc->GetParError(1) << std::endl;
      std::cout << "Chi-square *100.: " << fitFunc->GetChisquare()*100.0 << std::endl;
      std::cout << "NDF: " << fitFunc->GetNDF() << std::endl;
      */


      TCanvas *canvas = new TCanvas("canvas", "My TGraph", 800, 600);
      canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
      graph->SetMarkerColor(kRed);
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(0.5);
      graph->Draw("AP");
      graph->GetXaxis()->SetTitle("Z [mm]");
      graph->GetYaxis()->SetTitle("#Delta #Phi [rad]");
      // Add title at the top
      TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
      title->SetFillColor(0);
      title->SetTextAlign(22);
      title->Draw("same");
      // Set the range
      graph->GetXaxis()->SetLimits(z_min-200, z_max+200);
      graph->GetYaxis()->SetRangeUser(-2.0, 2.0);

      // Draw dashed line y = 0
      //TLine *zeroLine = new TLine(z[0]-200, 0, z[2]+200, 0);
      //zeroLine->SetLineColor(kBlack);
      //zeroLine->SetLineStyle(2);
      //zeroLine->SetLineWidth(2);
      //zeroLine->Draw("same");

      // Draw dphi/dX in text
      TLatex dphidx;
      dphidx.SetTextSize(0.03);
      dphidx.SetTextAlign(13);
      dphidx.DrawLatexNDC(0.15, 0.88, Form("Station (%d - %d - %d)", nthStation-2, nthStation-1, nthStation));
      dphidx.DrawLatexNDC(0.15, 0.85, Form("#frac{d#phi}{dZ} = %f, fitted with 2 CHs #candidate %d", slope_alpha, 999));
      dphidx.DrawLatexNDC(0.15, 0.78, Form("CHs in Segment = %d", (int)hit_diag.size()));
      dphidx.DrawLatexNDC(0.15, 0.70, "Fit Results: ");
      dphidx.DrawLatexNDC(0.15, 0.67, Form("Slope (a): %f +/- %f", fitFunc->GetParameter(0), fitFunc->GetParError(0)));
      dphidx.DrawLatexNDC(0.15, 0.64, Form("Intercept (b): %f +/- %f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
      dphidx.DrawLatexNDC(0.15, 0.61, Form("Chi-square: %f", fitFunc->GetChisquare()));
      dphidx.DrawLatexNDC(0.15, 0.58, Form("NDF: %d", fitFunc->GetNDF()));

      //select best candidate
      double ChiNDF = fitFunc->GetChisquare()/fitFunc->GetNDF();
      //std::cout<<"chiNDF*100.0 = "<<ChiNDF*100.0<<std::endl;


      //delete
      delete canvas;
      delete graph;
      delete title;
      //delete zeroLine;
      //Plot end

      //return fitting values
      alpha = fitFunc->GetParameter(0);
      beta = fitFunc->GetParameter(1);
      chindf = ChiNDF;
  }

//--------------------------------------------------------------------------------//
  void PhiZSeedFinder::ev5_select_best_segments_step_04(const std::vector<ev5_HitsInNthStation>& HitsInTimeCluster, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi){

      if(!(ThisIsBestSegment.size() >= 2)) return;
      if(!(ThisIsBestSegment_Diag.size() >= 2)) return;

      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << " ev5_select_best_segments_step_04     " << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
    //-----------------------------------------------------------------
    //    Remove same combination in the segments(# of Segments >= 2)
    //-----------------------------------------------------------------
    //---------------------------------------------------------------------------------------------------------
    //  (1). Check all hitIDs in each segments and remove segments if there is duplicate
    //---------------------------------------------------------------------------------------------------------
    std::vector<std::vector<int>> hitID_list;
    hitID_list.clear();
    for(int i=0; i<(int)ThisIsBestSegment.size(); i++){
      std::vector<int> temp_hitID_list;
      temp_hitID_list.clear();
      for(int j=0; j<(int)ThisIsBestSegment.at(i).size(); j++){
        temp_hitID_list.push_back(ThisIsBestSegment.at(i).at(j).hitID);
      }
      //Sort temp_hitID_list in increasing order
      sort(temp_hitID_list.begin(), temp_hitID_list.end());
      hitID_list.push_back(temp_hitID_list);
    }

    // Remove duplicates
    std::sort(hitID_list.begin(), hitID_list.end());
    hitID_list.erase(std::unique(hitID_list.begin(), hitID_list.end()), hitID_list.end());


    //---------------------------------------------------------------------------------------------------------
    // (2). delete if there is overlap: example, 0-1-2 and 0-1-2-4, 1-2-4 and 0-1-2-4, in this case 0-1-2-4 will remain
    //---------------------------------------------------------------------------------------------------------
    std::vector<int> delete_index;
    delete_index.clear();
    for (int i = 0; i < (int)hitID_list.size(); i++) {
      //std::cout << "i = " << i << std::endl;
      for (int j = 0; j <(int)hitID_list.size(); j++) {
        if(i == j) continue;
        if(hitID_list.at(j) == hitID_list.at(i)){
          //std::cout << "found j = " << j << std::endl;
          continue;
        }
        int size = 0;
        for (int k = 0; k <(int)hitID_list.at(j).size(); k++) {
          for (int l = 0; l <(int)hitID_list.at(i).size(); l++) {
            if(hitID_list.at(i).at(l) == hitID_list.at(j).at(k)) size++;
          }
        }
        if(size == (int)hitID_list.at(i).size()) delete_index.push_back(i);
      }
    }

    //  make a new hitID after removing the overlap event
    std::vector<std::vector<int>> new_hitID_list;
    for(int i = 0; i <(int)hitID_list.size(); i++){
      int go = 1;
      for(int j=0; j<(int)delete_index.size(); j++){
        if(delete_index.at(j) == i) go = 0;
      }
      if(go == 1) new_hitID_list.push_back(hitID_list.at(i));
    }


    //---------------------------------------------------------------------------------------------------------
    // (3)  find an endex in "new hitID" correspond to the ThisIsBestSegment
    //---------------------------------------------------------------------------------------------------------
    //std::cout<<"Step (3-0)"<<std::endl;
    std::vector<int> find_index;
    find_index.clear();
    for (int i = 0; i < (int)new_hitID_list.size(); i++) {
        bool flag = 0;
        int index_for_segment = 0;
        for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
          int count = 0;
          if(new_hitID_list.at(i).size() != ThisIsBestSegment.at(j).size()) continue;
          for (int k = 0; k<(int)new_hitID_list.at(i).size(); k++) {
              for(int l=0; l<(int)ThisIsBestSegment.at(j).size(); l++){
                if(new_hitID_list.at(i).at(k) == ThisIsBestSegment.at(j).at(l).hitID) count++;
              }
          }
          if(count == (int)new_hitID_list.at(i).size()){
            flag = 1;
            index_for_segment = j;
            break;
          }
        }
        if(flag == 1) find_index.push_back(index_for_segment);

    }


  }//end ev5_select_best_segments_step_04


//--------------------------------------------------------------------------------//
  void PhiZSeedFinder::ev5_select_best_segments_step_05(const std::vector<ev5_HitsInNthStation>& HitsInTimeCluster, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi){

      if(!(ThisIsBestSegment.size() >= 2)) return;
      if(!(ThisIsBestSegment_Diag.size() >= 2)) return;

    //std::cout << "-----------------------------------" << std::endl;
    //std::cout << "-----------------------------------" << std::endl;
    //std::cout << " ev5_select_best_segments_step_05     " << std::endl;
    //std::cout << "-----------------------------------" << std::endl;
    //std::cout << "-----------------------------------" << std::endl;

    //---------------------------------------------------------------------------------------------------------
    // Select good segments and remove unwanted segments
    //---------------------------------------------------------------------------------------------------------
    std::vector<std::vector<ev5_HitsInNthStation>> temp_ThisIsBestSegment;
    std::vector<std::vector<ev5_Segment>> temp_ThisIsBestSegment_Diag;
    temp_ThisIsBestSegment = ThisIsBestSegment;
    temp_ThisIsBestSegment_Diag = ThisIsBestSegment_Diag;


    std::vector<int> delete_segmentID;
    delete_segmentID.clear();
    for(int i=0; i<(int)temp_ThisIsBestSegment.size(); i++){
      //std::cout << "temp i =  " << i<< std::endl;
      for(int j=0; j<(int)temp_ThisIsBestSegment.size(); j++){
        if(i==j) continue;
        //std::cout << "temp j =  " << j<< std::endl;
        double alpha[2] = {temp_ThisIsBestSegment_Diag.at(i).at(0).alpha, temp_ThisIsBestSegment_Diag.at(j).at(0).alpha};
        //std::cout << "alpha[0]*alpha[1] =    " <<alpha[0]*alpha[1]<< std::endl;
        double alpha_diff = 0.0;
        double alpha_sigma = 0.0003199;//obtained from  "ev5_plot_AlphaDiagBestSegment"
        if(alpha[0]*alpha[1] > 0){//both slopes are negitice or positive
          alpha_diff = fabs(alpha[0]  - alpha[1]);
        }
        if(alpha[0]*alpha[1] < 0){//either slope is positive or negative
          alpha[0] = fabs(alpha[0]);
          alpha[1] = fabs(alpha[1]);
          alpha_diff = alpha[0] + alpha[1];
        }
        if(alpha_diff > alpha_sigma*5) continue;

        double chiNDF[2] = {temp_ThisIsBestSegment_Diag.at(i).at(0).chiNDF, temp_ThisIsBestSegment_Diag.at(j).at(0).chiNDF};
        int TotalHitinSegments[2] = {(int)temp_ThisIsBestSegment_Diag.at(i).size(), (int)temp_ThisIsBestSegment_Diag.at(j).size()};
        double fraction[2] = {0.0};
        int overlappedHits = 0;
        for(int k=0; k<(int)temp_ThisIsBestSegment.at(i).size(); k++){
          for(int l=0; l<(int)temp_ThisIsBestSegment.at(j).size(); l++){
            int hitID[2] = {temp_ThisIsBestSegment.at(i).at(k).hitID, temp_ThisIsBestSegment.at(j).at(l).hitID};
            if(hitID[0] == hitID[1]) overlappedHits++;
          }
        }
        fraction[0] = (double)overlappedHits/TotalHitinSegments[0];
        fraction[1] = (double)overlappedHits/TotalHitinSegments[1];
        //std::cout << "TotalHitinSegments[0]/[1] =  " <<TotalHitinSegments[0]<<"/"<<TotalHitinSegments[1]<<std::endl;
        //std::cout << "overlappedHits =  " <<overlappedHits<<std::endl;
        //std::cout << "fraction[0]/[1] =  " <<fraction[0]<<"/"<<fraction[1]<<std::endl;


        //(1) 2 segment are not identical but ComboHits are overlapped and can remove 1 segemnt
        if(fraction[0] <= 0.3 and fraction[1]  >= 0.7){
          delete_segmentID.push_back(i);
        }
        //(2) 2 segment are identical and ComboHits are overlapped so select only 1 segment
        if(fraction[0] >= 0.7 and fraction[1]  >= 0.7){
          if(chiNDF[0] < chiNDF[1]) delete_segmentID.push_back(j);
          else delete_segmentID.push_back(i);
        }
        //(3) 2 segments are not identical and ComboHits are not overlappedcan so these 2 segments can be candidates
      }
    }


      // Sort the delete_segmentID vector in increasing order
      std::sort(delete_segmentID.begin(), delete_segmentID.end());
      //Remove duplicates
      delete_segmentID.erase(std::unique(delete_segmentID.begin(), delete_segmentID.end()), delete_segmentID.end());

    ThisIsBestSegment.clear();
    ThisIsBestSegment_Diag.clear();
    for(int i=0; i<(int)temp_ThisIsBestSegment.size(); i++){
      int flag = 0;
      for(int j=0; j<(int)delete_segmentID.size(); j++){
        if(delete_segmentID[j] == i) flag = 1;
      }
      if(flag == 1) continue;
      //std::cout << "hit_found = " << hit_found[i] <<std::endl;
      ThisIsBestSegment.push_back(temp_ThisIsBestSegment.at(i));
      ThisIsBestSegment_Diag.push_back(temp_ThisIsBestSegment_Diag.at(i));
    }

    //std::cout<<"ThisIsBestSegment size = "<<ThisIsBestSegment.size()<<std::endl;

  }//end ev5_select_best_segments_step_05

//-----------------------------------------------------------------------------
  void PhiZSeedFinder::ev5_select_best_segments_step_06(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi){


      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << " ev5_select_best_segments_step_06     " << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;

    //-----------------------------------------------------------------
    //    Remove identical segments
    //-----------------------------------------------------------------

    //---------------------------------------------------------------------------------------------------------
    //  (1). Check all hitIDs in each segments and remove segments if there is duplicate
    //---------------------------------------------------------------------------------------------------------
    //std::cout<<"Print the intial hitID_list "<<std::endl;
    std::vector<std::vector<int>> hitID_list;
    hitID_list.clear();
    for(int i=0; i<(int)all_ThisIsBestSegment.size(); i++){
      std::vector<int> temp_hitID_list;
      temp_hitID_list.clear();
      for(int j=0; j<(int)all_ThisIsBestSegment.at(i).size(); j++){
        temp_hitID_list.push_back(all_ThisIsBestSegment.at(i).at(j).hitID);
      }
      //Sort temp_hitID_list in increasing order
      sort(temp_hitID_list.begin(), temp_hitID_list.end());
      hitID_list.push_back(temp_hitID_list);
    }

    // Remove duplicates
    std::sort(hitID_list.begin(), hitID_list.end());
    hitID_list.erase(std::unique(hitID_list.begin(), hitID_list.end()), hitID_list.end());

    //---------------------------------------------------------------------------------------------------------
    // (2). delete if there is overlap: example, 0-1-2 and 0-1-2-4, 1-2-4 and 0-1-2-4, in this case 0-1-2-4 will remain
    //---------------------------------------------------------------------------------------------------------
    std::vector<int> delete_index;
    delete_index.clear();
    for (int i = 0; i < (int)hitID_list.size(); i++) {
      //std::cout << "i = " << i << std::endl;
      for (int j = 0; j <(int)hitID_list.size(); j++) {
        if(i == j) continue;
        if(hitID_list.at(j) == hitID_list.at(i)){
          //std::cout << "found j = " << j << std::endl;
          continue;
        }
        int size = 0;
        for (int k = 0; k <(int)hitID_list.at(j).size(); k++) {
          for (int l = 0; l <(int)hitID_list.at(i).size(); l++) {
            if(hitID_list.at(i).at(l) == hitID_list.at(j).at(k)) size++;
          }
        }
        if(size == (int)hitID_list.at(i).size()) delete_index.push_back(i);
      }
    }

    //  make a new hitID after removing the overlap event
    std::vector<std::vector<int>> new_hitID_list;
    for(int i = 0; i <(int)hitID_list.size(); i++){
      int go = 1;
      for(int j=0; j<(int)delete_index.size(); j++){
        if(delete_index.at(j) == i) go = 0;
      }
      if(go == 1) new_hitID_list.push_back(hitID_list.at(i));
    }

    //---------------------------------------------------------------------------------------------------------
    // (3)  find an endex in "new hitID" correspond to the all_ThisIsBestSegment
    //---------------------------------------------------------------------------------------------------------
    //std::cout<<"Step (3-0)"<<std::endl;
    std::vector<int> find_index;
    find_index.clear();
    for (int i = 0; i < (int)new_hitID_list.size(); i++) {
      //std::cout << "i = " << i << std::endl;
        bool flag = 0;
        int index_for_segment = 0;
        for(int j=0; j<(int)all_ThisIsBestSegment.size(); j++){
          int count = 0;
          if(new_hitID_list.at(i).size() != all_ThisIsBestSegment.at(j).size()) continue;
          for (int k = 0; k<(int)new_hitID_list.at(i).size(); k++) {
              for(int l=0; l<(int)all_ThisIsBestSegment.at(j).size(); l++){
                if(new_hitID_list.at(i).at(k) == all_ThisIsBestSegment.at(j).at(l).hitID) count++;
              }
          }
          if(count == (int)new_hitID_list.at(i).size()){
            flag = 1;
            index_for_segment = j;
            break;
          }
        }
        if(flag == 1) find_index.push_back(index_for_segment);
    }
    //find endex is only the parameter used for next step
    //std::cout<<"(int)find_index.size() = "<<(int)find_index.size()<<std::endl;

    std::vector<std::vector<ev5_HitsInNthStation>> temp_ThisIsBestSegment;
    std::vector<std::vector<ev5_Segment>> temp_ThisIsBestSegment_Diag;
    temp_ThisIsBestSegment.clear();
    temp_ThisIsBestSegment_Diag.clear();
    //select the best candidate
    for(int i=0; i<(int)all_ThisIsBestSegment.size(); i++){
      bool go = 0;
      for(int j=0; j<(int)find_index.size(); j++){
        if(find_index[j] == i) go = 1;
      }
      if(go != 1) continue;
      //push back segment
      temp_ThisIsBestSegment.push_back(all_ThisIsBestSegment.at(i));
      temp_ThisIsBestSegment_Diag.push_back(all_ThisIsBestSegment_Diag.at(i));
    }

    //---------------------------------------------------------------------------------------------------------
    // Select good segments and remove unwanted segments
    //---------------------------------------------------------------------------------------------------------

    std::vector<int> delete_segmentID;
    delete_segmentID.clear();
    for(int i=0; i<(int)temp_ThisIsBestSegment.size(); i++){
      //std::cout << "temp i =  " << i<< std::endl;
      for(int j=0; j<(int)temp_ThisIsBestSegment.size(); j++){
        if(i==j) continue;
        //std::cout << "temp j =  " << j<< std::endl;
        double alpha[2] = {temp_ThisIsBestSegment_Diag.at(i).at(0).alpha, temp_ThisIsBestSegment_Diag.at(j).at(0).alpha};
        //std::cout << "alpha[0]/alpha[1] =    " <<alpha[0]<<"/"<<alpha[1]<< std::endl;
        double alpha_diff = 0.0;
        double alpha_sigma = 0.0003199;//obtained from  "ev5_plot_AlphaDiagBestSegment"
        if(alpha[0]*alpha[1] > 0){//both slopes are negitice or positive
          alpha_diff = fabs(alpha[0]  - alpha[1]);
        }
        if(alpha[0]*alpha[1] < 0){//either slope is positive or negative
          alpha[0] = fabs(alpha[0]);
          alpha[1] = fabs(alpha[1]);
          alpha_diff = alpha[0] + alpha[1];
        }
    //    if(alpha_diff > alpha_sigma*5) continue;
        //std::cout << "alpha_diff =   " <<alpha_diff<< std::endl;

        double chiNDF[2] = {temp_ThisIsBestSegment_Diag.at(i).at(0).chiNDF, temp_ThisIsBestSegment_Diag.at(j).at(0).chiNDF};
        int TotalHitinSegments[2] = {(int)temp_ThisIsBestSegment_Diag.at(i).size(), (int)temp_ThisIsBestSegment_Diag.at(j).size()};
        double fraction[2] = {0.0};
        int overlappedHits = 0;
        for(int k=0; k<(int)temp_ThisIsBestSegment.at(i).size(); k++){
          for(int l=0; l<(int)temp_ThisIsBestSegment.at(j).size(); l++){
            int hitID[2] = {temp_ThisIsBestSegment.at(i).at(k).hitID, temp_ThisIsBestSegment.at(j).at(l).hitID};
            if(hitID[0] == hitID[1]) overlappedHits++;
          }
        }
        fraction[0] = (double)overlappedHits/TotalHitinSegments[0];
        fraction[1] = (double)overlappedHits/TotalHitinSegments[1];
        //std::cout << "TotalHitinSegments[0]/[1] =  " <<TotalHitinSegments[0]<<"/"<<TotalHitinSegments[1]<<std::endl;
        //std::cout << "overlappedHits =  " <<overlappedHits<<std::endl;
        //std::cout << "fraction[0]/[1] =  " <<fraction[0]<<"/"<<fraction[1]<<std::endl;

        if(alpha_diff > alpha_sigma*5) continue;
        //(1) 2 segment are not identical but ComboHits are overlapped and can remove 1 segemnt
        if(fraction[0] <= 0.3 and fraction[1]  >= 0.6){
          delete_segmentID.push_back(i);
        }

        //(2) 2 segment are identical and ComboHits are overlapped so select only 1 segment
        if(fraction[0] >= 0.7 and fraction[1]  >= 0.7){
          if(chiNDF[0] < chiNDF[1]) delete_segmentID.push_back(j);
          else delete_segmentID.push_back(i);
        }
      }
    }


      //std::cout << "before removed delete_segmentID.size() = " << (int)delete_segmentID.size()<<std::endl;
      // Sort the delete_segmentID vector in increasing order
      std::sort(delete_segmentID.begin(), delete_segmentID.end());
      //Remove duplicates
      delete_segmentID.erase(std::unique(delete_segmentID.begin(), delete_segmentID.end()), delete_segmentID.end());
      //std::cout << "after removed delete_segmentID.size() = " << (int)delete_segmentID.size()<<std::endl;

    all_ThisIsBestSegment.clear();
    all_ThisIsBestSegment_Diag.clear();
    for(int i=0; i<(int)temp_ThisIsBestSegment.size(); i++){
      int flag = 0;
      for(int j=0; j<(int)delete_segmentID.size(); j++){
        if(delete_segmentID[j] == i) flag = 1;
      }
      if(flag == 1) continue;
      //std::cout << "hit_found = " << hit_found[i] <<std::endl;
      all_ThisIsBestSegment.push_back(temp_ThisIsBestSegment.at(i));
      all_ThisIsBestSegment_Diag.push_back(temp_ThisIsBestSegment_Diag.at(i));
    }

    //std::cout<<"ThisIsBestSegment size = "<<all_ThisIsBestSegment.size()<<std::endl;



    //---------------------------------------------------------------------------------------------------------
    // Fit segments
    //---------------------------------------------------------------------------------------------------------
    //select the best candidate
    for(int i=0; i<(int)all_ThisIsBestSegment.size(); i++){

      //Plot graph: Phi vs. Z of each segment
      TGraph *graph = new TGraph();
      int p = 0;
      int nthStation = 999;
      double z_max = all_ThisIsBestSegment_Diag.at(i).at(0).z;
      double z_min = all_ThisIsBestSegment_Diag.at(i).at(0).z;
      double m = all_ThisIsBestSegment_Diag.at(i).at(0).alpha;
      int nComboHitsSegment = (int)all_ThisIsBestSegment_Diag.at(i).size();
      for(int j=0; j<nComboHitsSegment; j++){
        // z as x-axis and phi as y-axis
        double z = all_ThisIsBestSegment_Diag.at(i).at(j).z;
        double phi = all_ThisIsBestSegment_Diag.at(i).at(j).deltaphi;
        graph->SetPoint(p++, z, phi);

        if(z_min > z) z_min = z;
        if(z_max < z) z_max = z;
        if(nthStation > all_ThisIsBestSegment_Diag.at(i).at(j).station) nthStation = all_ThisIsBestSegment_Diag.at(i).at(j).station;
      }

      //Fitting function: f(x) = a*x + b
      TF1* fitFunc = new TF1("linearFit", "[0]*x + [1]", graph->GetXaxis()->GetXmin(), graph->GetXaxis()->GetXmax());
      fitFunc->SetParameter(0, m);
      fitFunc->SetParameter(1, 0.0);
      fitFunc->SetLineColor(kBlue);
      fitFunc->SetLineStyle(2);
      graph->Fit(fitFunc, "R");

      //Print fitting result
      //std::cout << "Fit Results:" << std::endl;
      //std::cout << "   Slope (a): " << fitFunc->GetParameter(0) << " +/- " << fitFunc->GetParError(0) << std::endl;
      //std::cout << "   Intercept (b): " << fitFunc->GetParameter(1) << " +/- " << fitFunc->GetParError(1) << std::endl;
      //std::cout << "   Chi-square: " << fitFunc->GetChisquare() << std::endl;
      //std::cout << "   NDF: " << fitFunc->GetNDF() << std::endl;



      TCanvas *canvas = new TCanvas("canvas", "My TGraph", 800, 600);
      canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
      graph->SetMarkerColor(kRed);
      graph->SetMarkerStyle(20);
      graph->SetMarkerSize(0.5);
      graph->Draw("AP");
      graph->GetXaxis()->SetTitle("Z [mm]");
      graph->GetYaxis()->SetTitle("#Delta #Phi [rad]");
      // Add title at the top
      TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
      //int run = eventId.run();
      //int subrun = eventId.subRun();
      //std::stringstream eventStringStream;
      //eventStringStream << "run: " << run << " subRun: " << subrun << " event: " << eventNumber;
      //title->AddText(Form("Z vs. Phi (Event, TC, ST-1st-2nd-3rd, #candidate) = (%d, %d, %d-%d-%d, %d) (CE only)", eventNumber, TC, nthStation-2, nthStation-1, nthStation, i));
      title->SetFillColor(0);
      title->SetTextAlign(22);
      title->Draw("same");
      // Set the range
      graph->GetXaxis()->SetLimits(z_min-200, z_max+200);
      graph->GetYaxis()->SetRangeUser(-2.0, 2.0);

      // Draw dashed line y = 0
      //TLine *zeroLine = new TLine(z[0]-200, 0, z[2]+200, 0);
      //zeroLine->SetLineColor(kBlack);
      //zeroLine->SetLineStyle(2);
      //zeroLine->SetLineWidth(2);
      //zeroLine->Draw("same");

      // Draw dphi/dX in text
      TLatex dphidx;
      dphidx.SetTextSize(0.03);
      dphidx.SetTextAlign(13);
      dphidx.DrawLatexNDC(0.15, 0.88, Form("Station (%d - %d - %d)", nthStation, nthStation+1, nthStation+2));
      //dphidx.DrawLatexNDC(0.15, 0.85, Form("#frac{d#phi}{dZ} = %f, fitted with 2 CHs #candidate %d", m, i));
      dphidx.DrawLatexNDC(0.15, 0.78, Form("CHs in Segment = %d", (int)all_ThisIsBestSegment_Diag.at(i).size()));
      dphidx.DrawLatexNDC(0.15, 0.73, Form("#Delta #phi = %.1f", threshold_deltaphi));
      dphidx.DrawLatexNDC(0.15, 0.70, "Fit Results: ");
      dphidx.DrawLatexNDC(0.15, 0.67, Form("Slope (a): %f +/- %f", fitFunc->GetParameter(0), fitFunc->GetParError(0)));
      dphidx.DrawLatexNDC(0.15, 0.64, Form("Intercept (b): %f +/- %f", fitFunc->GetParameter(1), fitFunc->GetParError(1)));
      dphidx.DrawLatexNDC(0.15, 0.61, Form("Chi-square: %f", fitFunc->GetChisquare()));
      dphidx.DrawLatexNDC(0.15, 0.58, Form("NDF: %d", fitFunc->GetNDF()));

      //select best candidate
      //double chiNDF = fitFunc->GetChisquare()/fitFunc->GetNDF();
      //std::cout<<"chiNDF*100.0 = "<<chiNDF*100.0<<std::endl;

      //delete
      delete canvas;
      delete graph;
      delete title;
      //delete zeroLine;
      //Plot end

    }//end


  }//end ev5_select_best_segments_step_06

//-----------------------------------------------------------------------------
  void PhiZSeedFinder::ev5_select_best_segments_step_07(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi, int& NumberOfSegments){


      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << " ev5_select_best_segments_step_07  " << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;

      std::vector<std::vector<ev5_HitsInNthStation>> best_segments = all_ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>> diag_best_segments = all_ThisIsBestSegment_Diag;
      NumberOfSegments = 0;

      //std::cout << "best_segments = " <<best_segments.size()<< std::endl;
    /*for(int i=0; i<(int)best_segments.size(); i++){
      std::cout << "ComboHits in segments = " <<best_segments.at(i).size()<< std::endl;
      for(int j=0; j<(int)best_segments.at(i).size(); j++){
        std::cout<<"No/hitID/x/y/z/phi = "<<i<<"/"<<best_segments.at(i).at(j).hitID<<"/"<<best_segments.at(i).at(j).x<<"/"<<best_segments.at(i).at(j).y<<"/"<<best_segments.at(i).at(j).z<<"/"<<best_segments.at(i).at(j).phi<<std::endl;
      }
    }  */

    //-----------------------------------------------------------------
    //    Fill slope value "alpha" and stations from each segments
    //-----------------------------------------------------------------
    int nSegmentsInTC = (int)diag_best_segments.size();
    std::vector<double> alpha;
    std::vector<std::vector<int>> nstations;
    std::vector<std::vector<double>> phi_window;
    for(int i=0; i<nSegmentsInTC; i++){
      //std::cout<<"best_segments.at("<<i<<").size() = "<<diag_best_segments.at(i).size()<<std::endl;
      for(int j=0; j<(int)diag_best_segments.at(i).size(); j++){
        alpha.push_back(diag_best_segments.at(i).at(j).alpha);
        break;
      }
      std::vector<int> station;
      station.clear();
      std::vector<double> phi;
      phi.clear();
      double hit_phi[2] = {9999.9, -9999.9};//[0] = min, [1] = max
      for(int j=0; j<(int)diag_best_segments.at(i).size(); j++){
        if(hit_phi[0] > best_segments.at(i).at(j).phi) hit_phi[0] = best_segments.at(i).at(j).phi;
        if(hit_phi[1] < best_segments.at(i).at(j).phi) hit_phi[1] = best_segments.at(i).at(j).phi;
        //std::cout<<"phi = "<<best_segments.at(i).at(j).phi<<std::endl;
        station.push_back(diag_best_segments.at(i).at(j).station);
      }
      phi.push_back(hit_phi[0]);//phi[rad] min
      phi.push_back(hit_phi[1]);//phi[rad] max
      phi_window.push_back(phi);
      nstations.push_back(station);

    }

    /*for(int i=0; i<(int)alpha.size(); i++){
        std::cout<<"alpha["<<i<<"] = "<<alpha.at(i)<<std::endl;
      for(int j=0; j<(int)nstations.at(i).size(); j++){
        std::cout<<"station["<<i<<"] = "<<nstations.at(i).at(j)<<std::endl;
      }
    }
    for(int i=0; i<(int)phi_window.size(); i++){
        std::cout<<"phi_window["<<i<<"] = "<<phi_window.at(i).size()<<std::endl;
      for(int j=0; j<(int)phi_window.at(i).size(); j++){
        std::cout<<"phi_window["<<i<<"] = "<<phi_window.at(i).at(j)<<std::endl;
      }
    }*/



    //---------------------------------------------------------------------------------------------------------
    //  (1). Find Combinations
    //---------------------------------------------------------------------------------------------------------
    std::vector<std::vector<int>> ncandidate;
    for(int i=0; i<(int)alpha.size(); i++){
      std::vector<int> combi;
      combi.clear();
      combi.push_back(i);
      for(int j=0; j<(int)alpha.size(); j++){
        if(i == j) continue;

        //station cut
        //segmnet i and j should not have the same station
        int flag = 0;
        for(int k=0; k<(int)nstations.at(i).size(); k++){
          for(int l=0; l<(int)nstations.at(j).size(); l++){
            if(nstations.at(i).at(k) == nstations.at(j).at(l)) flag = 1;
          }
        }
        if(flag == 1) continue;

        //phi cut
        //segmnet i and j should not have the same station
        double phi_0[4] = {0.0};
        double phi_1[2] = {0.0};
        double phi_2[4] = {0.0};
        double phi_3[2] = {0.0};
        int phi_flag[4] ={0, 0, 0, 0};
        //For segment 1
        //std::cout<<"alpha = "<<alpha.at(i)<<"/"<<alpha.at(j)<<std::endl;
        //std::cout<<"phi_window = "<<phi_window.at(i).at(0)<<"/"<<phi_window.at(i).at(1)<<std::endl;
            if((phi_window.at(i).at(0) < -2.0 and 2.0 < phi_window.at(i).at(1)) or (phi_window.at(i).at(1) < -2.0 and 2.0 < phi_window.at(i).at(0))){
              double min = -9999.9;
              double max = 9999.9;
              for(int k=0; k<(int)best_segments.at(i).size(); k++){
                if(best_segments.at(i).at(k).phi < 0.0){
                  if(min < best_segments.at(i).at(k).phi) min = best_segments.at(i).at(k).phi;
                }
                if(best_segments.at(i).at(k).phi > 0.0){
                  if(max > best_segments.at(i).at(k).phi) max = best_segments.at(i).at(k).phi;
                }
              }
              phi_0[0] = -M_PI;
              phi_0[1] = min;
              phi_0[2] = max;
              phi_0[3] = M_PI;
              phi_flag[0] = 1;
              //std::cout<<"min/max = "<<min<<"/"<<max<<std::endl;
            }else{
              phi_1[0] = phi_window.at(i).at(0);
              phi_1[1] = phi_window.at(i).at(1);
              phi_flag[1] = 1;
            }
        //For segment 2
        //std::cout<<"phi_window = "<<phi_window.at(j).at(0)<<"/"<<phi_window.at(j).at(1)<<std::endl;
            if((phi_window.at(j).at(0) < -2.0 and 2.0 < phi_window.at(j).at(1)) or (phi_window.at(j).at(1) < -2.0 and 2.0 < phi_window.at(j).at(0))){
              double min = -9999.9;
              double max = 9999.9;
              for(int k=0; k<(int)best_segments.at(j).size(); k++){
                if(best_segments.at(j).at(k).phi < 0.0){
                  if(min < best_segments.at(j).at(k).phi) min = best_segments.at(j).at(k).phi;
                }
                if(best_segments.at(j).at(k).phi > 0.0){
                  if(max > best_segments.at(j).at(k).phi) max = best_segments.at(j).at(k).phi;
                }
              }
              phi_2[0] = -M_PI;
              phi_2[1] = min;
              phi_2[2] = max;
              phi_2[3] = M_PI;
              phi_flag[2] = 1;
              //std::cout<<"min/max = "<<min<<"/"<<max<<std::endl;
            }else{
              phi_3[0] = phi_window.at(j).at(0);
              phi_3[1] = phi_window.at(j).at(1);
              phi_flag[3] = 1;
            }



        if(phi_flag[0] == 1 and phi_flag[3] == 1){
          if(phi_0[1] < phi_3[0] and phi_3[1] < phi_0[2]) flag = 1;
        }
        if(phi_flag[1] == 1 and phi_flag[2] == 1){
          if(phi_2[1] < phi_1[0] and phi_1[1] < phi_2[2]) flag = 1;
        }
        if(phi_flag[1] == 1 and phi_flag[3] == 1){
          if(phi_1[1] < phi_3[0] and phi_1[1] < phi_3[0]) flag = 1;
          if(phi_1[0] > phi_3[1] and phi_1[0] > phi_3[1]) flag = 1;
        }

        if(flag == 1) continue;
        //std::cout<<"Pass PhiWindow"<<flag<<std::endl;



        double alpha_diff = 9999.9;
        if(alpha.at(i)*alpha.at(j) > 0) alpha_diff = fabs(alpha.at(i) - alpha.at(j));
        if(alpha.at(i)*alpha.at(j) < 0) alpha_diff = fabs(alpha.at(i)) + fabs(alpha.at(j));

        //double alpha_sigma = 0.0004659;//obtained from  "ev5_plot_AlphaDiagBestSegment"
        //double alpha_sigma = 0.0004659;//obtained from  "ev5_plot_AlphaDiagBestSegment"
        double alpha_sigma = 0.0003199;//obtained from  "ev5_plot_AlphaDiagBestSegment"
        //std::cout<<"alpha["<<i<<"] = "<<alpha.at(i)<<"/ alpha["<<j<<"] = "<<alpha.at(j)<<std::endl;
        //std::cout<<"alpha_diff = "<<alpha_diff<<std::endl;
        if(alpha_diff < alpha_sigma*5){//5 sigma region
        //if(alpha_diff < alpha_sigma){//5 sigma region
          combi.push_back(j);
        }
      }
      ncandidate.push_back(combi);
    }


    //Sort ncandidate in increasing order:
    //Example: before {1, 2, 3}, {1, 2, 3, 3}, {3, 2, 1}, {4, 2, 1}
    //after {1, 2, 3}, {1, 2, 3, 3}, {1, 2, 3}, {1, 2, 4}
    for (auto& row : ncandidate) {
      std::sort(row.begin(), row.end());
    }

    // Sort ncandidate again
    // Example: before {1, 2, 3}, {1, 2, 3, 3}, {1, 2, 3}, {1, 2, 4}
    // after {1, 2, 3}, {1, 2, 3}, {1, 2, 3, 3}, {1, 2, 4}
    std::sort(ncandidate.begin(), ncandidate.end());


    // Delete duplicate
    //Example: before {1, 2, 3}, {1, 2, 3}, {1, 2, 3, 3}, {1, 2, 4}
    //after {1, 2, 3}, {1, 2, 3}, {1, 2, 3}, {1, 2, 4}
    for (int j = 0; j < static_cast<int>(ncandidate.size()); j++) {
      for (int i = 0; i < static_cast<int>(ncandidate.at(j).size()); i++) {
        ncandidate.at(j).erase(std::unique(ncandidate.at(j).begin(), ncandidate.at(j).end()), ncandidate.at(j).end());
       }
    }

    //delete if there is duplicate
    //Example: before {1, 2, 3}, {1, 2, 3, 3}, {1, 2, 3}, {1, 2, 4}
    //after {1, 2, 3}, {1, 2, 3, 3}, {1, 2, 4}
    ncandidate.erase(std::unique(ncandidate.begin(), ncandidate.end()), ncandidate.end());

    // delete if there is overlap: example, 0-1-2 and 0-1-2-4, 1-2-4 and 0-1-2-4, in this case 0-1-2-4 will remain
    std::vector<int> delete_index;
    delete_index.clear();
    for (int i = 0; i < static_cast<int>(ncandidate.size()); i++) {
    //std::cout << "i = " << i << std::endl;
    for (int j = 0; j < static_cast<int>(ncandidate.size()); j++) {
        if(i == j) continue;
        if(ncandidate.at(j) == ncandidate.at(i)) continue;
        int size = 0;
        for (int k = 0; k < static_cast<int>(ncandidate.at(j).size()); k++) {
          for (int l = 0; l < static_cast<int>(ncandidate.at(i).size()); l++) {
            if(ncandidate.at(i).at(l) == ncandidate.at(j).at(k)) size++;
          }
        }
        if(size == (int)ncandidate.at(i).size()) delete_index.push_back(i);
      }
    }
    std::vector<std::vector<int>> new_ncandidate;
    new_ncandidate.clear();
    for(int i = 0; i < static_cast<int>(ncandidate.size()); i++){
      int go = 1;
      for(int j=0; j<(int)delete_index.size(); j++){
        if(delete_index.at(j) == i) go = 0;
      }
      if(go == 1) new_ncandidate.push_back(ncandidate.at(i));
    }
    ncandidate.clear();
    ncandidate = new_ncandidate;

    std::vector<std::vector<ev5_HitsInNthStation>> hitsIn_best_segments;
    hitsIn_best_segments.clear();
    for(int i = 0; i < static_cast<int>(ncandidate.size()); i++){
      std::vector<ev5_HitsInNthStation> hits;
      hits.clear();
      for(int j = 0; j < static_cast<int>(ncandidate.at(i).size()); j++){
        int index = ncandidate.at(i).at(j);
        //hits = best_segments.at(index);
        for(int k = 0; k < static_cast<int>(best_segments.at(index).size()); k++){
          hits.push_back(best_segments.at(index).at(k));
        }
      }
      hitsIn_best_segments.push_back(hits);
    }

    //Remove duplicates from each inner vector in hitsIn_best_segments
    //Sort ncandidate in increasing order:
    //Example: before {11, 12, 13}, {21, 12, 23, 21, 33}, {32, 7, 7, 8, 9}
    //after {11, 12, 13}, {12, 21, 21, 23, 33}, {7, 7, 8, 9, 32}
    // then
    //Example: before  {11, 12, 13}, {12, 21, 21, 23, 33}, {7, 7, 8, 9, 32}
    //after {11, 12, 13}, {12, 21, 23, 33}, {7, 8, 9, 32}
    for(auto& hits : hitsIn_best_segments){
      std::sort(hits.begin(), hits.end(), [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b) {
        return a.hitID < b.hitID;
      });
      // Erase duplicates
      auto last = std::unique(hits.begin(), hits.end(), [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b) {
        return a.hitID == b.hitID;
      });
      hits.erase(last, hits.end());
    }


    all_ThisIsBestSegment = hitsIn_best_segments;
    NumberOfSegments = (int)hitsIn_best_segments.size();

  }//end ev5_select_best_segments_step_07

//-----------------------------------------------------------------------------
// Function to create the time clusters
//-----------------------------------------------------------------------------
  void PhiZSeedFinder::initTimeCluster(TimeCluster& tc){
    int nstrs = tc._strawHitIdxs.size();
    tc._nsh = 0;
    float tacc(0),tacc2(0),xacc(0),yacc(0),zacc(0),weight(0);
    for (int i=0; i<nstrs; i++) {
      int loc = tc._strawHitIdxs[i];
      const ComboHit* ch = &_data.chcol->at(loc);
      float htime = ch->correctedTime();
      float hwt = ch->nStrawHits();
      weight += hwt;
      tacc   += htime*hwt;
      tacc2  += htime*htime*hwt;
      xacc   += ch->pos().x()*hwt;
      yacc   += ch->pos().y()*hwt;
      zacc   += ch->pos().z()*hwt;
      tc._nsh += ch->nStrawHits();
    }
    tacc/=weight;
    tacc2/=weight;
    xacc/=weight;
    yacc/=weight;
    zacc/=weight;

    tc._t0._t0    = tacc;
    tc._t0._t0err = sqrtf(tacc2-tacc*tacc);
    tc._pos        = XYZVectorF(xacc, yacc, zacc);

  }

























//-----------------------------------------------------------------------------
  void PhiZSeedFinder::produce(art::Event& Event) {
    if (_debugLevel) printf("* >>> PhiZSeedFinder::produce  event number: %10i\n",Event.event());
//-----------------------------------------------------------------------------
// clear memory in the beginning of event processing and cache event pointer
//-----------------------------------------------------------------------------
    _data.InitEvent(&Event,_debugLevel);
//-----------------------------------------------------------------------------
// process event
//-----------------------------------------------------------------------------
    if (! findData(Event)) {
      const char* message = "mu2e::PhiZSeedFinder_module::produce: data missing or incomplete";
      throw cet::exception("RECO")<< message << endl;
    }

    _data._nTimeClusters = _data.tccol->size();
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
    std::unique_ptr<TimeClusterCollection> tccol1(new TimeClusterCollection);

//-----------------------------------------------------------------------------
// run delta finder, it also finds proton time clusters
//-----------------------------------------------------------------------------
  std::cout<<"PhiZSeed_kita"<<std::endl;
    for (int i=0; i<_data._nTimeClusters; i++) {
      const TimeCluster* tc = &_data.tccol->at(i);
      _finder->run(tc);
      //--------------------
      // PhiZ finder start
      //--------------------
      std::vector<ev5_HitsInNthStation> HitsInTimeCluster;
      HitsInTimeCluster.clear();
      ev5_FillHitsInTimeCluster(tc, HitsInTimeCluster);
      //for(int j=0; j<(int)HitsInTimeCluster.size(); j++){
      //  std::cout<<"hitID/x/y/z/phi/stations/strawhits = "<<HitsInTimeCluster[j].hitID<<"/"<<HitsInTimeCluster[j].x<<"/"<<HitsInTimeCluster[j].y<<"/"<<HitsInTimeCluster[j].z<<"/"<<HitsInTimeCluster[j].phi<<"/"<<HitsInTimeCluster[j].station<<"/"<<HitsInTimeCluster[j].strawhits<<std::endl;
      //}
      double thre_residual = 0.4;// +/-0.4[rad] delta phi window for the slope
      //-----------------------------------------
      // segment search in 3 consecutive stations
      //-----------------------------------------
      std::vector<std::vector<ev5_HitsInNthStation>> best_triplet_segments;
      best_triplet_segments.clear();
      std::vector<std::vector<ev5_Segment>>  diag_best_triplet_segments;
      diag_best_triplet_segments.clear();
      ev5_SegmentSearchInTriplet(HitsInTimeCluster, best_triplet_segments, diag_best_triplet_segments, thre_residual);

      //-----------------------------------------------------------------------------
      // make a new time cluster
      //-----------------------------------------------------------------------------
      for(int j=0; j<(int)best_triplet_segments.size(); j++){
      const std::vector<StrawHitIndex>& ordchcol = tc->hits();
      int nh = ordchcol.size();
      std::vector<int> hitindex;
      for(int ih=0; ih<nh; ih++){
        int ind = ordchcol[ih];
        const ComboHit* ch = &_data.chcol->at(ind);
        int hitID[2];
        hitID[0] = ch->_sid.asUint16();
        for(int k=0; k<(int)best_triplet_segments.at(j).size(); k++){
          hitID[1] = best_triplet_segments.at(j).at(k).hitID;
          if(hitID[0] == hitID[1]) hitindex.push_back(ih);
        }
      }
      //fill
      TimeCluster new_tc;
      for(int k=0; k<(int)hitindex.size(); k++){
        int ih = hitindex[k];
        new_tc._strawHitIdxs.push_back(ordchcol[ih]);
      }
      initTimeCluster(new_tc);
      tccol1->push_back(new_tc);
      }

    }

    // Save the output time cluster collection for the diagnostics
    if (_diagLevel > 0) {
        _data._tccolnew = tccol1.get();
        //make Histograms
        _hmanager->fillHistograms(&_data);
    }

     Event.put(std::move(tccol1));

  }
}
//-----------------------------------------------------------------------------
// macro that makes this class a module.
//-----------------------------------------------------------------------------
DEFINE_ART_MODULE(mu2e::PhiZSeedFinder)
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------

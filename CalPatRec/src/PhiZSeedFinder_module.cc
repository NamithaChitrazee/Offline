////////////////////////////////////////////////////////////////////////////
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
// C++ Standard Library
#include <algorithm>
#include <cmath>
#include <vector>
#include <numeric>    // for std::iota
// ROOT
#include "TH1F.h"
#include "TH2F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "TEfficiency.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TLine.h"
#include "TEllipse.h"
#include "TPaveText.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TF1.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMarker.h"
// art Framework
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/make_tool.h"
#include "art_root_io/TFileService.h"
// fhiclcpp
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
// Offline - Config Tools
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
// Offline - Geometry
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
// Offline - Magnetic Field
#include "Offline/BFieldGeom/inc/BFieldManager.hh"
// Offline - Conditions
#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
// Offline - Data Products
#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/DataProducts/inc/Helicity.hh"
// Offline - RecoDataProducts
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StereoHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "Offline/RecoDataProducts/inc/RobustHelix.hh"
#include "Offline/RecoDataProducts/inc/HelixHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
// Offline - MCDataProducts
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
// Offline - Global Constants
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
// Offline - CalPatRec
#include "Offline/CalPatRec/inc/ChannelID.hh"
#include "Offline/CalPatRec/inc/PhiZSeedFinder_types.hh"
#include "Offline/CalPatRec/inc/PhiZSeedFinderAlg.hh"
// Offline - Utilities
#include "Offline/Mu2eUtilities/inc/LsqSums2.hh"
#include "Offline/Mu2eUtilities/inc/LsqSums4.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
// CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
using namespace std;
using CalPatRec::ChannelID;
namespace mu2e {
  using namespace PhiZSeedFinderTypes;
  class PhiZSeedFinder: public art::EDProducer {
  public:
    struct Config {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>   shCollTag              {Name("shCollTag"        )     , Comment("StrawHit collection tag" ) };
      fhicl::Atom<art::InputTag>   chCollTag              {Name("chCollTag"          )     , Comment("ComboHit collection tag"    ) };
      fhicl::Atom<art::InputTag>   tcCollTag              {Name("tcCollTag"          )     , Comment("time cluster collection tag") };
      fhicl::Atom<art::InputTag>   sdmcCollTag            {Name("sdmcCollTag"        )             , Comment("StrawDigiMC collection tag" ) };
      fhicl::Atom<int>             debugLevel             {Name("debugLevel"         )     , Comment("debug level"                ) };
      fhicl::Atom<int>             diagLevel              {Name("diagLevel"          )     , Comment("diag level"                  ) };
      fhicl::Atom<int>             printErrors            {Name("printErrors"        )     , Comment("print errors"                ) };
      fhicl::Atom<int>             writeFilteredComboHits {Name("writeFilteredComboHits")  , Comment("0: write all CH, 1: write filtered CH") };
      fhicl::Atom<int>             writeStrawHits         {Name("writeStrawHits"     )     , Comment("1: write all SH, new flags" ) };
      fhicl::Atom<int>             testOrder              {Name("testOrder"          )     , Comment("1: test order"              ) };
      fhicl::Atom<bool>             testHitMask           {Name("testHitMask"        )     , Comment("true: test hit mask"        ) };
      fhicl::Sequence<std::string> goodHitMask            {Name("goodHitMask"        )     , Comment("good hit mask"              ) };
      fhicl::Sequence<std::string> bkgHitMask             {Name("bkgHitMask"         )     , Comment("background hit mask"        ) };
      fhicl::Sequence<int> Helicities                     {Name("Helicities"         )     , Comment("Helicity values"        ) };
      fhicl::Atom<bool>doSingleOutput                     {Name("doSingleOutput"     )     , Comment("Create a single ouputput with both helicities") };
      fhicl::Sequence<std::string> SaveHelixFlag          {Name("SaveHelixFlag"     )      , Comment("Save Helix Flag, 'HelixOK'") };
      fhicl::Table<McUtilsToolBase::Config> mcUtils     {Name("mcUtils"   )       , Comment("get MC info if debugging"      )  };
      fhicl::Table<PhiZSeedFinderTypes::Config> diagPlugin      {Name("diagPlugin"      )  , Comment("Diag plugin"           ) };
      fhicl::Table<PhiZSeedFinderAlg::Config>    finderParameters{Name("finderParameters") , Comment("finder alg parameters" ) };
    };
//-----------------------------------------------------------------------------
// PhiZSeedFinder Constructors
//-----------------------------------------------------------------------------
     struct ev5_HitsInNthStation {
          int hitIndice;
          double phi;
          int strawhits;
          double x;
          double y;
          double z;
          int station;
          int plane;
          int face;
          int panel;
          int hitID;
          int segmentIndex; //index of segment group
          bool used; //whether or not hit is used in fits, default = -1
          double phiDiag;
          double circleError2;
          double helixPhi;
          double helixPhiError2;
          int nturn;
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
      struct tripletPoint {
        const XYZVectorF* pos;
        int hitIndice;
      };
      struct triplet {
        tripletPoint i;
        tripletPoint j;
        tripletPoint k;
      };
      struct cHit {
        int hitIndice; // index of point in _data.chcol
        double circleError2;
        double helixPhi;
        double helixPhiError2;
        int helixPhiCorrection;
        //for 2Pi correction
        int segmentIndice;
        double ambigPhi;
        bool inHelix;
        bool used; // whether or not hit is used in fits
        bool isolated;
        bool averagedOut;
        bool notOnLine;
        bool uselessTripletSeed;
        bool notOnSegment;
        bool debugParticle; // only filled in debug mode -- true if mc particle, false if background
        int station;
        int plane;
        int face;
        int panel;
        double x;
        double y;
        double z;
        double phi;
        int strawhits;
      };
      struct cleanup {
        int tcindex;
        int tcindice;
        double chi2ndf;
      };
      // another struct for debugging
      //type1
      struct mcInfo {
        int   pdg;
        int   simID;
        int   nStrawHits;
        int tcIndex;
        float pMin;
        float pMax;
        float mcX0;
        float mcY0;
        float mcRadius;
      };
      //type2
      struct mcInfoList {
        int   nStrawHits;
        int station;
        int plane;
        int face;
        int panel;
        double x;
        double y;
        double z;
        double phi;
        float mom;
        int   pdg;
        int   simID;
      };
      struct mcDiffR {
        int   nHits;
        int   simID;
      };
      //weight info used for circle fit
      struct weightinfo {
        double  chi2_default;
        double  deltaR;
        double  sigma_default;
        double  weight_default;
        int  nstrawhits;
        double  sigma_wire;//[mm]
        double  sigma_transverse; //[mm]
        int     intersection;
        double  chi2_new;
        double  residula_wire; //[mm]
        double  residula_transverse; //[mm]
        double  sigma_new;//[mm]
        double  weight_new;
        double nhit_slope_a;
        double nhit_slope_b;
        double nhit_slope_c;
        double ortho_nhit_slope_a;
        double ortho_nhit_slope_b;
        double ortho_nhit_slope_c;
      };
      // tracker geometric information data list
      struct trackerData {
        int station;
        int plane;
        int face;
        int panel;
        double z;
      };
      // Define Point structure
      struct Point {
        double x;
        double y;
      };
      // Define Line structure
      struct Line {
        double slope;
        double intercept;
        bool isVertical;
        double verticalX;  // Only valid if the line is vertical
      };
      struct HelixFinderData {
      };
  protected:
//-----------------------------------------------------------------------------
// talk-to parameters: input collections and algorithm parameters
//-----------------------------------------------------------------------------
    art::InputTag     _shCollTag;
    art::InputTag    _chCollTag;
    art::InputTag    _tcCollTag;                 // time cluster coll tag
    art::InputTag    _sdmcCollTag;
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
// collections
//-----------------------------------------------------------------------------
    //const ComboHitCollection*      _chColl;
    //const TimeClusterCollection*   _tcColl;
    //const CaloClusterCollection*   _ccColl;
//-----------------------------------------------------------------------------
// cache event/geometry objects
//-----------------------------------------------------------------------------
    //TimeClusterCollection* tccol1;
    const StrawHitCollection* _shColl;
    HelixSeedCollection* _hsColl;
    const Tracker*               _tracker;
    //const DiskCalorimeter*       _calorimeter;
    const mu2e::Calorimeter*       _calorimeter;
    PhiZSeedFinderTypes::Data_t  _data;               // all data used
    int                           _testOrderPrinted;
    PhiZSeedFinderAlg*           _finder;
    int run;
    int subrun;
    int eventNumber;
    ::LsqSums2 _lineFitter;
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
// helper functions for TimeCluster
//-----------------------------------------------------------------------------
  void clusterInfo(int Tc);
//-----------------------------------------------------------------------------
// helper functions for phizseed finder
//-----------------------------------------------------------------------------
    enum SegmentComp{unique=-1,first=0,second=1};
    SegmentComp compareSegments(const std::vector<ev5_HitsInNthStation>& seg1,
    const std::vector<ev5_HitsInNthStation>& seg2);
    void ev5_FillHitsInTimeCluster(int tc);
    void ev5_FillHitsInTimeCluster_ver2(int tc);
    void ev5_SegmentSearchInTriplet(std::vector<std::vector<ev5_Segment>>& diag_best_triplet_segments, double thre_residual, int tc);
    double ev5_DeltaPhi(double x1, double y1, double x2, double y2);
    double ev5_ParticleDirection(double x1, double y1, double x2, double y2);
    double ev5_ResidualDeltaPhi(double alpha, double beta, double z1, double phi1, double z2, double phi2);
    bool ev5_TripletQuality(const std::vector<ev5_Segment> diag_hit);
    void ev5_select_best_segments_step_01(const std::vector<std::vector<ev5_HitsInNthStation>>& segment_candidates, const std::vector<std::vector<ev5_Segment>>& diag_segment_candidates, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi);
    void ev5_select_best_segments_step_02(int station, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi);
    void ev5_select_best_segments_step_03(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int station, int nCH, double threshold_deltaphi);
    void ev5_select_best_segments_step_030(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag);
    void ev5_select_best_segments_step_03A(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi);
    void ev5_select_best_segments_step_031(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi);
    void ev5_select_best_segments_step_04(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi);
    void ev5_select_best_segments_step_05(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi);
    void ev5_select_best_segments_step_06A(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag);
    void ev5_select_best_segments_step_06(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi);
    void ev5_select_best_segments_cleanup(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi);
    void countHits(const std::vector<ev5_HitsInNthStation>& seg1, const std::vector<ev5_HitsInNthStation>& seg2, unsigned& nh1, unsigned& nh2, unsigned& nover);
    void findchisq(std::vector<ev5_HitsInNthStation> const& segment, double& chizphi) const;
    void ev5_select_best_segments_step_07(std::vector<std::vector<ev5_HitsInNthStation>>& all_BestSegmentInfo, std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi, int& NumberOfSegments);
 void mergeSegmentsAll(
    std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment,
    std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag,
    double thre_residual);
    void ev5_select_best_segments_step_08(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi, int& NumberOfSegments);
    void ev5_fit_slope(const std::vector<ev5_Segment>& hit_diag, double& alpha, double& beta, double& chindf);
    void ev5_fit_slope_ver2(int i, double& alpha, double& beta, double& chindf);
    void ev5_fit_slope_ver3(int i, double& alpha, double& beta, double& chindf);
    void ev5_fit_slope_ver4(int index, double& alpha, double& alphaError, double& beta, double& betaError, double& chindf);
//-----------------------------------------------------------------------------
// helper functions for helix finder
//-----------------------------------------------------------------------------
    //void findHelix(int tc, int isegment, HelixSeedCollection& HSColl);
    void findHelix(int tc, int isegment, HelixSeedCollection& HSColl, HelixSeed& Temp_HSeed);
    void saveHelix(int tc, HelixSeed& Temp_HSeed);
    void segment_check(int tc, int isegment);
    void helix_check(int tc, int isegment);
    void get_diffrad(int tc, int isegment, double& r_diff);
    void initTriplet(triplet& trip, int& outcome);
    void initSeedCircle(int& outcome);
    void initHelixPhi();
    void computeHelixPhi(size_t& tcHitsIndex, double& xC, double& yC);
    void computeHelixPhi_ver2(int hitIndice, double& xC, double& yC, double& helixPhi, double& helixPhiError2);
    void tcHitsFill(int isegment);
    void tcHitsFill_Add(int isegment);
    void plot_PhiVsZ_OriginalTC(int tc);
    void plot_HelixPhiVsZ(int TC, int isegment);
    void plot_PhiVsZ_forSegment(int tc, int isegment);
    void plot_PhiVsZ_forSegment_ver2(int ith_segment, int jth_segment, double alpha, double beta, double Chi2NDF);
    void plot_PhiVsZ_forEachStep(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegmen, const char* filename, int tc, int loopIndex);
    void plot_CirclePhiVsZ_forSegment(int tc, int isegment);
    void plot_2PiAmbiguityPhiVsZ_forSegment(int tc, int isegment);
    void plot_2PiAmbiguityPhiVsZ_forSegment_mod(int tc, int isegment);
    void plot_RVsZ_forSegment(int tc, int isegment);
    double computeCircleResidual2(size_t& tcHitsIndex, double& xC, double& yC, double& rC);
    void computeCircleError2(size_t& tcHitsIndex, double& xC, double& yC, double& rC);
    double computeCircleError2_ver2(int hitIndice, int nStrawHits, double& xC, double& yC, double& rC);
    void mod_computeCircleError2(size_t& tcHitsIndex, double& xC, double& yC, double& rC);
    void plot_XVsY(int TC, int isegment, const char* filename, double& xC, double& yC, double& rC);
    void plot_XVsY_hit(int TC, int isegment, const char* filename, double& xC, double& yC, double& rC);
    int mcPreSelection(int tc);
    int tcPreSelection(int tc);
    void InitTrackGeometry(mu2e::RobustHelix track, std::vector<trackerData>& tracker_data);
    void calculateLineEquation(const Point& p1, const Point& p2, Line& line);
    void findIntersection(const Line& line1, const Line& line2, Point& intersection);
    double triangleArea(double x1, double y1, double x2, double y2, double x3, double y3);
    bool isPointInTriangle(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3);
//-----------------------------------------------------------------------------
// need to use mcUtils if in debug mode
//-----------------------------------------------------------------------------
    std::unique_ptr<McUtilsToolBase> _mcUtils;
//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------
    art::Event*             _event;
//-----------------------------------------------------------------------------
// debug
//-----------------------------------------------------------------------------
  std::vector<weightinfo> _printweight;
//-----------------------------------------------------------------------------
// stuff for doing segment search
//-----------------------------------------------------------------------------
   vector<ev5_HitsInNthStation> _HitsInCluster;
//-----------------------------------------------------------------------------
// stuff for doing helix search
//-----------------------------------------------------------------------------
    std::vector<std::vector<ev5_HitsInNthStation>> all_BestSegmentInfo;
    std::vector<std::vector<ev5_HitsInNthStation>> _segmentHits;//ComboHits and some helper variables for segments
    std::vector<cHit> _tcHits;
    ::LsqSums4 _circleFitter;
    float               _bz0;
    double _dphidz;
    double _fz0;
    std::vector<double> _FitRadius;
    std::vector<double> _DiffRadius;
    std::vector<Helicity> _hels; // helicity values to fit
    bool                                _doSingleOutput;
    TrkFitFlag        _saveflag; // write out all helices that satisfy these flags
//-----------------------------------------------------------------------------
// data members specifically for when doing debugging
//-----------------------------------------------------------------------------
    std::vector<std::vector<mcInfo>>   _simIDsPerTC; // filled once per TC
    std::vector<std::vector<mcInfoList>>  _simInfoPerTC; // filled once per TC
    //std::vector<std::vector<mcInfoList>>  _simPerTC; // filled once per SimID[/TC]
    int                                _tcIndex;
    int                                _simID;
    float                              _mcRadius;
    float                              _mcX0;
    float                              _mcY0;
    size_t                             _bestPlotIndex;
    size_t                             _bestLineSegment;
    int                                _mcParticleInTC;
//-----------------------------------------------------------------------------
// constants
//-----------------------------------------------------------------------------
  static constexpr float mmTconversion = CLHEP::c_light/1000.0;
//-----------------------------------------------------------------------------
// function
//-----------------------------------------------------------------------------
  void initTimeCluster(TimeCluster& tc);
//-----------------------------------------------------------------------------
// functions for debug mode and runDisplay mode
//-----------------------------------------------------------------------------
  void initDebugMode();
  void findBestTC();
  };
//-----------------------------------------------------------------------------
  PhiZSeedFinder::PhiZSeedFinder(const art::EDProducer::Table<Config>& config):
    art::EDProducer{config},
    _shCollTag             (config().shCollTag()         ),
    _chCollTag             (config().chCollTag()         ),
    _tcCollTag             (config().tcCollTag()         ),
    _sdmcCollTag           (config().sdmcCollTag()       ),
    _debugLevel            (config().debugLevel()        ),
    _diagLevel             (config().diagLevel()         ),
    _printErrors           (config().printErrors()       ),
    _testOrder             (config().testOrder()         ),
    _bkgHitMask            (config().bkgHitMask()        ),
    _doSingleOutput        (config().doSingleOutput()    ),
    _saveflag              (config().SaveHelixFlag()     )
    {
      std::vector<int> helvals = config().Helicities();
      for(auto hv : helvals) {
        Helicity hel(hv);
        _hels.push_back(hel);
      }
      if (_doSingleOutput){
        produces<HelixSeedCollection>();
      }else {
        std::vector<int> helvals = config().Helicities();
        for(auto hel : _hels) {
          produces<HelixSeedCollection>(Helicity::name(hel));
        }
      }
    consumes<TimeClusterCollection>(_tcCollTag);
    consumes<ComboHitCollection>   (_chCollTag);
    //produces<TimeClusterCollection>();
    //produces<HelixSeedCollection>();
    _finder = new PhiZSeedFinderAlg(config().finderParameters,&_data);
    _testOrderPrinted = 0;
    if (_diagLevel != 0) _hmanager = art::make_tool  <ModuleHistToolBase>(config().diagPlugin,"diagPlugin");
    else                 _hmanager = std::make_unique<ModuleHistToolBase>();
    if (_diagLevel != 0) _mcUtils = art::make_tool  <McUtilsToolBase>(config().mcUtils,"mcUtils");
    else              _mcUtils = std::make_unique<McUtilsToolBase>();
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
    GeomHandle<mu2e::Calorimeter> ch;
    _calorimeter = ch.get();
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(CLHEP::Hep3Vector(0.0, 0.0, 0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();
  }
//-----------------------------------------------------------------------------
  bool PhiZSeedFinder::findData(const art::Event& Evt) {
    auto tccH    = Evt.getValidHandle<mu2e::TimeClusterCollection>(_tcCollTag);
    _data.tccol = tccH.product();
    auto chcH = Evt.getValidHandle<mu2e::ComboHitCollection>(_chCollTag);
    _data.chcol = chcH.product();
    if (_diagLevel == 1) {
      auto sdmccH = Evt.getValidHandle<StrawDigiMCCollection>(_sdmcCollTag);
      _data.sdmcColl   = sdmccH.product();
      auto shcH = Evt.getValidHandle<mu2e::StrawHitCollection>(_shCollTag);
      _shColl   = shcH.product();
    }
    if(_diagLevel == 1) return ((_data.tccol != nullptr) and (_data.chcol != nullptr) and (_data.sdmcColl != nullptr));
    else return (_data.tccol != nullptr) and (_data.chcol != nullptr);
  }
//-----------------------------------------------------------------------------
    void PhiZSeedFinder::ev5_FillHitsInTimeCluster(int tc) {
    // loop over ComboHits in a TimeCluster
    for (size_t i = 0; i < _data.tccol->at(tc)._strawHitIdxs.size(); i++) {
      int hitIndice = _data.tccol->at(tc)._strawHitIdxs[i];
      std::vector<StrawDigiIndex> shids;
      _data.chcol->fillStrawDigiIndices(hitIndice, shids);
      ev5_HitsInNthStation hitsincluster;
      hitsincluster.hitIndice = hitIndice;
      hitsincluster.hitID = hitIndice;
      hitsincluster.phi = _data.chcol->at(hitIndice).pos().phi();
      hitsincluster.x = _data.chcol->at(hitIndice).pos().x();
      hitsincluster.y = _data.chcol->at(hitIndice).pos().y();
      hitsincluster.z = _data.chcol->at(hitIndice).pos().z();
      hitsincluster.strawhits = (int)shids.size();
      hitsincluster.station =  _data.chcol->at(hitIndice).strawId().station();
      hitsincluster.plane   =  _data.chcol->at(hitIndice).strawId().plane();
      hitsincluster.face    =  _data.chcol->at(hitIndice).strawId().face();
      hitsincluster.panel   =  _data.chcol->at(hitIndice).strawId().panel();
      hitsincluster.used    = false;
      hitsincluster.phiDiag = 0.0;
      hitsincluster.segmentIndex = 0;
      hitsincluster.circleError2 = 0.0;
      hitsincluster.helixPhi = 0.0;
      hitsincluster.helixPhiError2 = 0.0;
      hitsincluster.nturn = 0;
      _HitsInCluster.push_back(hitsincluster);
    }
    //sort the vector in ascending order of the z-coordinate
    std::sort(_HitsInCluster.begin(), _HitsInCluster.end(), [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b) { return a.z < b.z; } );
}
//-----------------------------------------------------------------------------
    void PhiZSeedFinder::ev5_FillHitsInTimeCluster_ver2(int tc) {
   /* int kstation = 18;
    int kplane = 0;//o to 35
    int kface = 2;//0 or 1
    int kpanel = 3;// 0-2-4 or 1-3-5
    //loop over station
    for (int istation=17; istation >= 0; istation--){
      //loop over plane
      for (int iplane=35; iplane >= 0; istation--){
        //loop over face
        for (int iface=0; iface < 2; istation++){
          hitsincluster.panel = _data.chcol->at(hitIndice).strawId().plane();
          push_back(hitsincluste);
        }
        kplane++;
      }
    }
*/
    // loop over ComboHits in a TimeCluster
    for (size_t i = 0; i < _data.tccol->at(tc)._strawHitIdxs.size(); i++) {
      int hitIndice = _data.tccol->at(tc)._strawHitIdxs[i];
      std::vector<StrawDigiIndex> shids;
      _data.chcol->fillStrawDigiIndices(hitIndice, shids);
      ev5_HitsInNthStation hitsincluster;
      hitsincluster.hitIndice = hitIndice;
      hitsincluster.hitID = hitIndice;
      hitsincluster.phi = _data.chcol->at(hitIndice).pos().phi();
      hitsincluster.x = _data.chcol->at(hitIndice).pos().x();
      hitsincluster.y = _data.chcol->at(hitIndice).pos().y();
      hitsincluster.z = _data.chcol->at(hitIndice).pos().z();
      hitsincluster.strawhits = (int)shids.size();
      hitsincluster.station =  _data.chcol->at(hitIndice).strawId().station();
      hitsincluster.plane   =  _data.chcol->at(hitIndice).strawId().plane();
      hitsincluster.face    =  _data.chcol->at(hitIndice).strawId().face();
      hitsincluster.panel   =  _data.chcol->at(hitIndice).strawId().panel();
      hitsincluster.used    = false;
      hitsincluster.phiDiag = 0.0;
      hitsincluster.segmentIndex = 0;
      hitsincluster.circleError2 = 0.0;
      hitsincluster.helixPhi = 0.0;
      hitsincluster.helixPhiError2 = 0.0;
      hitsincluster.nturn = 0;
      _HitsInCluster.push_back(hitsincluster);
    }
    //sort the vector in ascending order of the z-coordinate
    std::sort(_HitsInCluster.begin(), _HitsInCluster.end(), [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b) { return a.z < b.z; } );
}
//-----------------------------------------------------------------------------
    void PhiZSeedFinder::clusterInfo(int Tc){
    std::cout << ">>> INFORMATION in PhiZFinder::clusterInfo: " << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << " Time Cluster " << Tc << std::endl;
    std::cout << " Total Combo Hits: " << _HitsInCluster.size() << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << std::left << std::setw(4) << "No. |";
    std::cout << std::left << std::setw(4) << "hitIndice |";
    std::cout << std::left << std::setw(4) << " station |";
    std::cout << std::left << std::setw(10) << " x [mm] |";
    std::cout << std::left << std::setw(10) << " y [mm] |";
    std::cout << std::left << std::setw(10) << " z [mm] |";
    std::cout << std::left << std::setw(10) << " R [mm] |";
    std::cout << std::left << std::setw(10) << " phi [rad] |";
    std::cout << std::endl;
    for(int j=0; j<(int)_HitsInCluster.size(); j++){
      //std::cout << "j = " << j << std::endl;
      std::cout << std::left << std::setw(4)  << j;
      std::cout << std::left << std::setw(4)  << _HitsInCluster.at(j).hitIndice;
      std::cout << std::left << std::setw(4)  << _HitsInCluster.at(j).station;
      std::cout << std::left << std::setw(10) << _HitsInCluster.at(j).x;
      std::cout << std::left << std::setw(10) << _HitsInCluster.at(j).y;
      std::cout << std::left << std::setw(10) << _HitsInCluster.at(j).z;
      std::cout << std::left << std::setw(10) << sqrt(_HitsInCluster.at(j).x*_HitsInCluster.at(j).x + _HitsInCluster.at(j).y*_HitsInCluster.at(j).y);
      std::cout << std::left << std::setw(10) << _HitsInCluster.at(j).phi;
      std::cout << std::endl;
    }
}
//-----------------------------------------------------------------------------
      void PhiZSeedFinder::ev5_SegmentSearchInTriplet(std::vector<std::vector<ev5_Segment>>& diag_best_triplet_segments, double thre_residual, int tc) {
      std::vector<std::vector<ev5_HitsInNthStation>> all_ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>> all_ThisIsBestSegment_Diag;
      all_ThisIsBestSegment.clear();
      all_ThisIsBestSegment_Diag.clear();
      // number of stations
      // search slope in 3 consecutive stations
      // Loop from 0 to 16 station
      ///int nstation = 18;
      for(int n=17; n>1; n--){
         //if(n==6) break;
          //----------------------------------------------------
          // Find triplet
          //----------------------------------------------------
          // Take combo hits in n-th, (n+1)-th, (n+2)-th stations
          std::vector<ev5_HitsInNthStation> Hits_In_Station[3];
          for(int k=0; k<3;k++) Hits_In_Station[k].clear();
          for(int j=0; j<(int)_HitsInCluster.size(); j++){
              ev5_HitsInNthStation hitsin_nthstation;
              hitsin_nthstation.hitIndice    = _HitsInCluster.at(j).hitIndice;
              hitsin_nthstation.phi          = _HitsInCluster.at(j).phi;
              hitsin_nthstation.strawhits    = _HitsInCluster.at(j).strawhits;
              hitsin_nthstation.x            = _HitsInCluster.at(j).x;
              hitsin_nthstation.y            = _HitsInCluster.at(j).y;
              hitsin_nthstation.z            = _HitsInCluster.at(j).z;
              hitsin_nthstation.station      = _HitsInCluster.at(j).station;
              hitsin_nthstation.plane        = _HitsInCluster.at(j).plane;
              hitsin_nthstation.face         = _HitsInCluster.at(j).face;
              hitsin_nthstation.panel        = _HitsInCluster.at(j).panel;
              hitsin_nthstation.hitID        = _HitsInCluster.at(j).hitID;
              if(_HitsInCluster.at(j).used == true) continue;
              if(n == _HitsInCluster.at(j).station) Hits_In_Station[0].push_back(hitsin_nthstation);
              if(n-1 == _HitsInCluster.at(j).station) Hits_In_Station[1].push_back(hitsin_nthstation);
              if(n-2 == _HitsInCluster.at(j).station) Hits_In_Station[2].push_back(hitsin_nthstation);
              //std::cout<<"hitsin_nthstation.strawhits = "<<hitsin_nthstation.strawhits<<std::endl;
          }
        //2 >= ComboHits [/station]
        size_t nCHsInStn_1 = Hits_In_Station[0].size();
        size_t nCHsInStn_2 = Hits_In_Station[1].size();
        size_t nCHsInStn_3 = Hits_In_Station[2].size();
        int nHitInStation = (int)nCHsInStn_1 + (int)nCHsInStn_2 + (int)nCHsInStn_3;
        //if(!(nHitInStation >= 5)) continue;
        // (1st, 2nd, 3rd) = (CH>=1, CH>=1, CH>=1)
        //if(!((int)nCHsInStn_1 >= 1)) continue;
        //if(!((int)nCHsInStn_2 >= 1)) continue;
        //if(!((int)nCHsInStn_3 >= 1)) continue;
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
        if(nHitInStation >= 5 and (int)nCHsInStn_1 >= 1 and (int)nCHsInStn_2 >= 1 and (int)nCHsInStn_3 >= 1){
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
        }
        //-------------------------------------------------------------
        //        End Find segments in 3 consecutive stations
        //-------------------------------------------------------------
      //---------------------------------------------------------------------
      // Select several best candidates in 3 consecutive stations
      //---------------------------------------------------------------------
      std::vector<std::vector<ev5_HitsInNthStation>> ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>> ThisIsBestSegment_Diag;
      ThisIsBestSegment.clear();
      ThisIsBestSegment_Diag.clear();
      if(segment_candidates.size() != 0 and diag_segment_candidates.size() != 0){
        // remove duplicate segments in 3 station based on hitID and fit each slope
        ev5_select_best_segments_step_01(segment_candidates, diag_segment_candidates, ThisIsBestSegment, ThisIsBestSegment_Diag, nHitInStation, thre_residual);
        //plot_PhiVsZ_forEachStep(ThisIsBestSegment, "step_01", tc, n);
        /*std::cout<<"kimo_step1"<<std::endl;
        for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
          std::cout<<"j = "<<j<<std::endl;
          for(int k=0; k<(int)ThisIsBestSegment.at(j).size(); k++){
          std::cout<<"k/strawhits/hitIndice = "<<k<<"/"<<ThisIsBestSegment.at(j).at(k).strawhits<<"/"<<ThisIsBestSegment.at(j).at(k).hitIndice<<std::endl;
          }
        }*/
        // collect remaining hits in 3 station and add those hits to the segment
        ev5_select_best_segments_step_02(n, ThisIsBestSegment, ThisIsBestSegment_Diag, nHitInStation, thre_residual);
        //plot_PhiVsZ_forEachStep(ThisIsBestSegment, "step_02", tc, n);
        /*std::cout<<"kimo_step2"<<std::endl;
        for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
          std::cout<<"j = "<<j<<std::endl;
          for(int k=0; k<(int)ThisIsBestSegment.at(j).size(); k++){
          std::cout<<"k/strawhits = "<<k<<"/"<<ThisIsBestSegment.at(j).at(k).strawhits<<std::endl;
          }
        }*/
        //---------------------------------------------------------------------
        // flag hits as "used" in  best candidates in 3 consecutive stations
        //---------------------------------------------------------------------
        for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
          for(int k=0; k<(int)ThisIsBestSegment.at(j).size(); k++){
            ThisIsBestSegment.at(j).at(k).used = true;
            for(int l=0; l<(int)_HitsInCluster.size(); l++){
              if(_HitsInCluster.at(l).hitIndice == ThisIsBestSegment.at(j).at(k).hitIndice) _HitsInCluster.at(l).used = true;
            }
          }
        }
      }
      //---------------------------------------------------------------------
      // Before extending the slope, flag hits already used in best candidates
      //---------------------------------------------------------------------
      for(int j=0; j<(int)all_ThisIsBestSegment.size(); j++){
        for(int k=0; k<(int)all_ThisIsBestSegment.at(j).size(); k++){
          all_ThisIsBestSegment.at(j).at(k).used = true;
          for(int l=0; l<(int)_HitsInCluster.size(); l++){
            if(_HitsInCluster.at(l).hitIndice == all_ThisIsBestSegment.at(j).at(k).hitIndice) _HitsInCluster.at(l).used = true;
          }
        }
      }
      //---------------------------------------------------------------------
      // extend the slope to only +/-1 neighboring station and then, add hits to the segment and fit the segment
      //---------------------------------------------------------------------
      std::cout<<" Hussain = "<<n<<std::endl;
      std::cout<<" ThisIsBestSegment.size() = "<<ThisIsBestSegment.size()<<std::endl;
      std::cout<<" all_ThisIsBestSegment.size() = "<<all_ThisIsBestSegment.size()<<std::endl;
      if (!all_ThisIsBestSegment.empty()) {
        for (auto &seg : all_ThisIsBestSegment) {
          ThisIsBestSegment.push_back(seg);
        }
      }
      if (!all_ThisIsBestSegment_Diag.empty()) {
        for (auto &seg_diag : all_ThisIsBestSegment_Diag) {
          ThisIsBestSegment_Diag.push_back(seg_diag);
        }
      }
      std::cout<<" After ThisIsBestSegment.size() = "<<ThisIsBestSegment.size()<<std::endl;
      if(ThisIsBestSegment.size() == 0) continue;
      if(ThisIsBestSegment_Diag.size() == 0) continue;
      ev5_select_best_segments_step_03(ThisIsBestSegment, ThisIsBestSegment_Diag, n, nHitInStation, thre_residual);
      /*std::cout<<"kimo_step3"<<std::endl;
      for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
        std::cout<<"j = "<<j<<std::endl;
        for(int k=0; k<(int)ThisIsBestSegment.at(j).size(); k++){
        std::cout<<"k/strawhits/hitIndice = "<<k<<"/"<<ThisIsBestSegment.at(j).at(k).strawhits<<"/"<<ThisIsBestSegment.at(j).at(k).hitIndice<<std::endl;
        }
      }*/
      // falg hits as used
      for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
        for(int k=0; k<(int)ThisIsBestSegment.at(j).size(); k++){
          ThisIsBestSegment.at(j).at(k).used = true;
          for(int l=0; l<(int)_HitsInCluster.size(); l++){
            if(_HitsInCluster.at(l).hitIndice == ThisIsBestSegment.at(j).at(k).hitIndice) _HitsInCluster.at(l).used = true;
          }
        }
      }
      // falg hits as used and removed used hits in other segments to protect hits already used for other station cycle
      //std::cout<<" Here ThisIsBestSegment.size() = "<<ThisIsBestSegment.size()<<std::endl;
      //ev5_select_best_segments_step_030(ThisIsBestSegment, ThisIsBestSegment_Diag);
      //std::cout<<"kimo_step30"<<std::endl;
      //std::cout<<"ThisIsBestSegment size = "<<ThisIsBestSegment.size()<<std::endl;
      /*for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
        std::cout<<"j = "<<j<<std::endl;
        for(int k=0; k<(int)ThisIsBestSegment.at(j).size(); k++){
        std::cout<<"k/strawhits/hitIndice = "<<k<<"/"<<ThisIsBestSegment.at(j).at(k).strawhits<<"/"<<ThisIsBestSegment.at(j).at(k).hitIndice<<std::endl;
        }
      }*/
      // extend the slope to neighboring and then, add hits to the segment and fit the segment
      //ev5_select_best_segments_step_03A(ThisIsBestSegment, ThisIsBestSegment_Diag, nHitInStation, thre_residual);
      _segmentHits.clear();
      _segmentHits = ThisIsBestSegment;
      //if(n == 8 or n == 7)
      //plot_PhiVsZ_forEachStep(_segmentHits, "step_030", tc, n);
      /*std::cout<<"kimo_step3A"<<std::endl;
      std::cout<<"ThisIsBestSegment size = "<<ThisIsBestSegment.size()<<std::endl;
      for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
        std::cout<<"j = "<<j<<std::endl;
        for(int k=0; k<(int)ThisIsBestSegment.at(j).size(); k++){
        std::cout<<"k/strawhits/hitIndice = "<<k<<"/"<<ThisIsBestSegment.at(j).at(k).strawhits<<"/"<<ThisIsBestSegment.at(j).at(k).hitIndice<<std::endl;
        }
      }*/
      //plot_PhiVsZ_forEachStep(ThisIsBestSegment, "step_03", tc, n);
      // extend the slope to neighboring, allowing 2 gapped station and then, add hits to the segment and fit the segment
      //ev5_select_best_segments_step_031(ThisIsBestSegment, ThisIsBestSegment_Diag, nHitInStation, thre_residual);
      //plot_PhiVsZ_forEachStep(ThisIsBestSegment, "step_031", tc, n);
      /*std::cout<<"kimo_step31"<<std::endl;
      for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
        std::cout<<"j = "<<j<<std::endl;
        for(int k=0; k<(int)ThisIsBestSegment.at(j).size(); k++){
        std::cout<<"k/strawhits/hitIndice = "<<k<<"/"<<ThisIsBestSegment.at(j).at(k).strawhits<<"/"<<ThisIsBestSegment.at(j).at(k).hitIndice<<std::endl;
        }
      }*/
      //if Segments >= 2, remove duplicate segments based on hitID
      //ev5_select_best_segments_step_04(ThisIsBestSegment, ThisIsBestSegment_Diag, nHitInStation, thre_residual);
      //plot_PhiVsZ_forEachStep(ThisIsBestSegment, "step_04", tc, n);
      /*std::cout<<"kimo_step4"<<std::endl;
      for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
        std::cout<<"j = "<<j<<std::endl;
        for(int k=0; k<(int)ThisIsBestSegment.at(j).size(); k++){
        std::cout<<"k/strawhits/hitIndice = "<<k<<"/"<<ThisIsBestSegment.at(j).at(k).strawhits<<"/"<<ThisIsBestSegment.at(j).at(k).hitIndice<<std::endl;
        }
      }*/
      //if Segments >= 2, remove duplicate segments based on slope value and fraction of overlapped hits
      //ev5_select_best_segments_step_05(ThisIsBestSegment, ThisIsBestSegment_Diag, nHitInStation, thre_residual);
      //plot_PhiVsZ_forEachStep(ThisIsBestSegment, "step_05", tc, n);
      _segmentHits = ThisIsBestSegment;
      /*std::cout<<"kimo_step5"<<std::endl;
      for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
        std::cout<<"j = "<<j<<std::endl;
        for(int k=0; k<(int)ThisIsBestSegment.at(j).size(); k++){
        std::cout<<"k/strawhits/hitIndice = "<<k<<"/"<<ThisIsBestSegment.at(j).at(k).strawhits<<"/"<<ThisIsBestSegment.at(j).at(k).hitIndice<<std::endl;
        }
      }*/
      //-------------------------------
      // Fill best candidate
      //-------------------------------
      if (!all_ThisIsBestSegment.empty()) all_ThisIsBestSegment.clear();
      if (!all_ThisIsBestSegment_Diag.empty()) all_ThisIsBestSegment_Diag.clear();
      for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
        //_segmentHits.push_back(ThisIsBestSegment.at(j));
        diag_best_triplet_segments.push_back(ThisIsBestSegment_Diag.at(j));
      }
      for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
        all_ThisIsBestSegment.push_back(ThisIsBestSegment.at(j));
        all_ThisIsBestSegment_Diag.push_back(ThisIsBestSegment_Diag.at(j));
      }
      std::cout<<"Fill best size = "<<all_ThisIsBestSegment.size()<<std::endl;
      for(int j=0; j<(int)all_ThisIsBestSegment.size(); j++){
        std::cout<<"j = "<<j<<std::endl;
        for(int k=0; k<(int)all_ThisIsBestSegment.at(j).size(); k++){
        //std::cout<<"k/strawhits/hitIndice = "<<k<<"/"<<all_ThisIsBestSegment.at(j).at(k).strawhits<<"/"<<all_ThisIsBestSegment.at(j).at(k).hitIndice<<std::endl;
        }
      }
      //break;
      }
      //-------------------------------
      // End loop on stations
      //-------------------------------
      _segmentHits.clear();
      _segmentHits = all_ThisIsBestSegment;
      std::cout<<"End loop on stations"<<std::endl;
      for(int j=0; j<(int)_segmentHits.size(); j++){
        std::cout<<"j = "<<j<<std::endl;
        for(int k=0; k<(int)_segmentHits.at(j).size(); k++){
        std::cout<<"k/hitIndice = "<<k<<"/"<<_segmentHits.at(j).at(k).hitIndice<<std::endl;
        }
      }
      //plot_PhiVsZ_forEachStep(_segmentHits, "step_06", tc, 999);
      /*std::cout<<"kimo_step6_before"<<std::endl;
      for(int j=0; j<(int)_segmentHits.size(); j++){
        std::cout<<"j = "<<j<<std::endl;
        for(int k=0; k<(int)_segmentHits.at(j).size(); k++){
        std::cout<<"k/hitIndice = "<<k<<"/"<<_segmentHits.at(j).at(k).hitIndice<<std::endl;
        }
      }*/
      // falg hits as used and removed used hits in other segments to protect hits already used for other station cycle
      //std::cout<<" Here ThisIsBestSegment.size() = "<<ThisIsBestSegment.size()<<std::endl;
      ev5_select_best_segments_step_06A(all_ThisIsBestSegment, all_ThisIsBestSegment_Diag);
      _segmentHits.clear();
      _segmentHits = all_ThisIsBestSegment;
      //std::cout<<"kimo_step6A"<<std::endl;
      //std::cout<<"ThisIsBestSegment size = "<<ThisIsBestSegment.size()<<std::endl;
      /*for(int j=0; j<(int)ThisIsBestSegment.size(); j++){
        std::cout<<"j = "<<j<<std::endl;
        for(int k=0; k<(int)ThisIsBestSegment.at(j).size(); k++){
        std::cout<<"k/strawhits/hitIndice = "<<k<<"/"<<ThisIsBestSegment.at(j).at(k).strawhits<<"/"<<ThisIsBestSegment.at(j).at(k).hitIndice<<std::endl;
        }
      }*/
      //plot_PhiVsZ_forEachStep(_segmentHits, "step_06", tc, 999);
      // remove duplicate segments based on hitID. remove segments based on slope value, fraction of overlapped hits, and Chi2/NDF
      ev5_select_best_segments_step_06(all_ThisIsBestSegment, all_ThisIsBestSegment_Diag, thre_residual);
      _segmentHits.clear();
      _segmentHits = all_ThisIsBestSegment;
      std::cout<<"ev5_select_best_segments_step_06 size = "<<_segmentHits.size()<<std::endl;
      for(int j=0; j<(int)_segmentHits.size(); j++){
        std::cout<<"size = "<<all_ThisIsBestSegment.at(j).size()<<std::endl;
        //for(int k=0; k<(int)_segmentHits.at(j).size(); k++){
        //std::cout<<"k/hitIndice = "<<k<<"/"<<_segmentHits.at(j).at(k).hitIndice<<std::endl;
        //}
      }
      for(int j=0; j<(int)all_ThisIsBestSegment_Diag.size(); j++){
        std::cout<<"size = "<<all_ThisIsBestSegment_Diag.at(j).size()<<std::endl;
        //for(int k=0; k<(int)all_ThisIsBestSegment_Diag.at(j).size(); k++){
        //std::cout<<"k/hitIndice = "<<k<<"/"<<_segmentHits.at(j).at(k).hitIndice<<std::endl;
        //}all_ThisIsBestSegment_Diag.at(i).j.
      }
      // ------------------------------------------------------------
      // Remove redundant or overlapping segments from the list of candidate "best segments".
      // If two segments share many of the same hits, keep only the better one (based on #hits or chi2).
      // This operates at the "segment" level, not the "hit" level.
      // Steps:
//   1. Loop over all pairs of segments.
//   2. Compare each pair with compareSegments().
//   3. If they are unique (little/no overlap) -> keep both.
//   4. If they overlap strongly:
//        - Prefer the one with more hits (if difference > threshold).
//        - If similar #hits, choose the one with lower chi2.
//   5. Remove the worse segment from the vectors.
// ------------------------------------------------------------
      ev5_select_best_segments_cleanup(all_ThisIsBestSegment, all_ThisIsBestSegment_Diag, thre_residual);
      _segmentHits.clear();
      _segmentHits = all_ThisIsBestSegment;
      //plot_PhiVsZ_forEachStep(_segmentHits, "step_06", tc, 999);
      //remove hits in diag that are not used in segments
      for (size_t i = 0; i < all_ThisIsBestSegment_Diag.size(); i++) {
    // Reference to the i-th group of diag segments
    std::vector<ev5_Segment> &diagSegs = all_ThisIsBestSegment_Diag[i];
    // Reference to the i-th group of hit segments
    std::vector<ev5_HitsInNthStation> &hitSegs = all_ThisIsBestSegment[i];
    // Loop through diagSegs with index j
    for (size_t j = 0; j < diagSegs.size(); ) {
        ev5_Segment dseg = diagSegs[j];
        bool foundMatch = false;
        // Compare against all hits in the same group
        for (size_t k = 0; k < hitSegs.size(); k++) {
            ev5_HitsInNthStation hseg = hitSegs[k];
            double diff_z = std::fabs(dseg.z - hseg.z);
            if (diff_z < 0.1) {
                foundMatch = true; // found a match
                break;
            }
        }
        if (!foundMatch) {
            // erase this diag segment if no match was found
            diagSegs.erase(diagSegs.begin() + j);
            // do not increment j, because elements shift left after erase
        } else {
            j++; // only increment if nothing was erased
        }
    }
}
      _segmentHits.clear();
      _segmentHits = all_ThisIsBestSegment;
      std::cout<<"merge segments based on slope value size = "<<_segmentHits.size()<<std::endl;
      for(int j=0; j<(int)_segmentHits.size(); j++){
        std::cout<<"size = "<<all_ThisIsBestSegment.at(j).size()<<std::endl;
        //for(int k=0; k<(int)_segmentHits.at(j).size(); k++){
        //std::cout<<"k/hitIndice = "<<k<<"/"<<_segmentHits.at(j).at(k).hitIndice<<std::endl;
        //}
      }
      for(int j=0; j<(int)all_ThisIsBestSegment_Diag.size(); j++){
        std::cout<<"size = "<<all_ThisIsBestSegment_Diag.at(j).size()<<std::endl;
        //for(int k=0; k<(int)all_ThisIsBestSegment_Diag.at(j).size(); k++){
        //std::cout<<"k/hitIndice = "<<k<<"/"<<_segmentHits.at(j).at(k).hitIndice<<std::endl;
        //}
      }
      // merge segments based on slope value, phi range. remove overlapped hits in the segments based on hitID
      //int number_of_merged_segments = 0;
      //ev5_select_best_segments_step_07(all_BestSegmentInfo, all_ThisIsBestSegment, all_ThisIsBestSegment_Diag, thre_residual, number_of_merged_segments);
      //_segmentHits.clear();
      //_segmentHits = all_ThisIsBestSegment;
      // merge segments based on slope value, phi range. remove overlapped hits in the segments based on hitID
      //int number_of_merged_segments = 0;
      //mergeSegmentsAll(all_BestSegmentInfo, all_ThisIsBestSegment, all_ThisIsBestSegment_Diag, thre_residual, number_of_merged_segments);
      double cut_thre_residual = 2.0;//[rad]
      mergeSegmentsAll(all_ThisIsBestSegment, all_ThisIsBestSegment_Diag, cut_thre_residual);
      _segmentHits.clear();
      _segmentHits = all_ThisIsBestSegment;
      //plot_PhiVsZ_forEachStep(_segmentHits, "step_07", tc, 999);
      std::cout<<"mergeSegmentsAll size = "<<_segmentHits.size()<<std::endl;
      for(int j=0; j<(int)_segmentHits.size(); j++){
        std::cout<<"j = "<<j<<std::endl;
        for(int k=0; k<(int)_segmentHits.at(j).size(); k++){
        std::cout<<"k/hitIndice = "<<k<<"/"<<_segmentHits.at(j).at(k).hitIndice<<std::endl;
        }
      }
      //std::cout<<"number_of_merged_segments = "<<number_of_merged_segments<<std::endl;
      //std::cout<<"number_of_merged_segments = "<<_segmentHits.size()<<std::endl;
/*      // remove duplicate merged segments (in case)
      std::cout<<"ev5_select_best_segments_step_07"<<std::endl;
      std::cout<<"all_ThisIsBestSegment.size() = "<<all_ThisIsBestSegment.size()<<std::endl;
      for(int j=0; j<(int)all_ThisIsBestSegment.size(); j++){
        for(int k=0; k<(int)all_ThisIsBestSegment.at(j).size(); k++){
          std::cout<<"all_ThisIsBestSegment.at(j).hitIndice = "<<all_ThisIsBestSegment.at(j).at(k).hitIndice<<std::endl;
        }
      }
      ev5_select_best_segments_step_08(all_ThisIsBestSegment, all_ThisIsBestSegment_Diag, thre_residual, number_of_merged_segments);
*/
      //fill
      _segmentHits.clear();
      _segmentHits = all_ThisIsBestSegment;
      diag_best_triplet_segments = all_ThisIsBestSegment_Diag;
}//end ev5_SegmentSearchInTriplet
//-----------------------------------------------------------------------------
 //PhiZSeedFinder::HelixComp PhiZSeedFinder::compareHelices(art::Event const& evt, HelixSeed const& h1, HelixSeed const& h2) {
 /*PhiZSeedFinder::HelixComp PhiZSeedFinder::compareHelices(HelixSeed const& h1, HelixSeed const& h2) {
    HelixComp retval(unique);
    unsigned nh1, nh2, nover;
    // count the StrawHit overlap between the helices
    countHits(evt,h1,h2, nh1, nh2, nover);
    unsigned minh = std::min(nh1, nh2);
    float chih1xy(0),chih1zphi(0),chih2xy(0),chih2zphi(0);
    //findchisq(evt,h1,chih1xy,chih1zphi);
    // Calculate the chi-sq of the helices
    //findchisq(evt,h2,chih2xy,chih2zphi);
    // overlapping helices: decide which is best
    if(nover >= _minnover && nover/float(minh) > _minoverfrac) {
      if(h1.caloCluster().isNonnull() && h2.caloCluster().isNull())
        retval = first;
      // Pick the one with a CaloCluster first
      else if( h2.caloCluster().isNonnull() && h1.caloCluster().isNull())
        retval = second;
      // then compare active StrawHit counts and if difference of the StrawHit counts greater than deltanh
      else if((nh1 > nh2) && (nh1-nh2) > _deltanh)
        retval = first;
      else if((nh2 > nh1) && (nh2-nh1) > _deltanh)
        retval = second;
      // finally compare chisquared: sum xy and fz
      else if(chih1xy+chih1zphi  < chih2xy+chih2zphi)
        retval = first;
      else
        retval = second;
    }
    return retval;
  }*/
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
  void PhiZSeedFinder::ev5_select_best_segments_step_01(const std::vector<std::vector<ev5_HitsInNthStation>>& segment_candidates, const std::vector<std::vector<ev5_Segment>>& diag_segment_candidates, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi){
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
      double slope_alpha = 0.0;//Slope (a)
      double slope_beta = 0.0;//Intercept (b)
      double ChiNDF = 0.0;
      ev5_fit_slope(diag_segment_candidates.at(i), slope_alpha, slope_beta, ChiNDF);
      //push back segment
      ThisIsBestSegment.push_back(segment_candidates.at(i));
      ThisIsBestSegment_Diag.push_back(diag_segment_candidates.at(i));
      int index = (int)ThisIsBestSegment_Diag.size()-1;
      for(int j=0; j<(int)ThisIsBestSegment_Diag.at(index).size(); j++){
        ThisIsBestSegment_Diag.at(index).at(j).alpha = slope_alpha;
        ThisIsBestSegment_Diag.at(index).at(j).beta = slope_beta;
        ThisIsBestSegment_Diag.at(index).at(j).chiNDF = ChiNDF;
      }
    }//end nSegmentInTriplet
  }//end ev5_select_best_segments_step_01
//-----------------------------------------------------------------------------
  void PhiZSeedFinder::ev5_select_best_segments_step_02(int station, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi){
      std::vector<std::vector<ev5_HitsInNthStation>> segments = ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>> diag_segments = ThisIsBestSegment_Diag;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << " ev5_select_best_segments_step_02  " << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //---------------------------------------------------------------------------------------------
      // Search remaining ComboHit candidates in 3 consecutive stations
      //---------------------------------------------------------------------------------------------
      //obtain the minimum station in the segment
      int max_station = station;
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
        for(int k=max_station; k>=max_station-2; k--){
        std::vector<ev5_HitsInNthStation> Hits_In_Station;
        Hits_In_Station.clear();
        for(int j=0; j<(int)_HitsInCluster.size(); j++){
            ev5_HitsInNthStation hitsin_nthstation;
            hitsin_nthstation.hitIndice    = _HitsInCluster.at(j).hitIndice;
            hitsin_nthstation.phi          = _HitsInCluster.at(j).phi;
            hitsin_nthstation.strawhits    = _HitsInCluster.at(j).strawhits;
            hitsin_nthstation.x            = _HitsInCluster.at(j).x;
            hitsin_nthstation.y            = _HitsInCluster.at(j).y;
            hitsin_nthstation.z            = _HitsInCluster.at(j).z;
            hitsin_nthstation.station      = _HitsInCluster.at(j).station;
            hitsin_nthstation.plane        = _HitsInCluster.at(j).plane;
            hitsin_nthstation.face        = _HitsInCluster.at(j).face;
            hitsin_nthstation.panel        = _HitsInCluster.at(j).panel;
            hitsin_nthstation.hitID        = _HitsInCluster.at(j).hitID;
            if(_HitsInCluster.at(j).used == true) continue;
            if(k != _HitsInCluster.at(j).station) continue;
            int flag_alreadyUsed = 0;
            for (const auto &element : hit_candidates) {
              if(element.hitID == _HitsInCluster.at(j).hitID) flag_alreadyUsed = 1;
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
void PhiZSeedFinder::ev5_select_best_segments_step_03(std::vector<std::vector<ev5_HitsInNthStation>> &ThisIsBestSegment, std::vector<std::vector<ev5_Segment>> &ThisIsBestSegment_Diag, int station, int nCH, double threshold_deltaphi){
    // ------------------------------------------------------------
    // Make local copies of the current "best" segments
    // ------------------------------------------------------------
    std::vector<std::vector<ev5_HitsInNthStation>> segments = ThisIsBestSegment;
    std::vector<std::vector<ev5_Segment>> diag_segments = ThisIsBestSegment_Diag;
    // For ev5_HitsInNthStation segments
for (auto &seg : segments) {
    std::sort(seg.begin(), seg.end(),
              [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b) {
                  return a.z < b.z; // ascending
              });
}
// For ev5_Segment segments
for (auto &seg : diag_segments) {
    std::sort(seg.begin(), seg.end(),
              [](const ev5_Segment& a, const ev5_Segment& b) {
                  return a.z < b.z; // ascending
              });
}
    // Containers for the new refined segments
    std::vector<std::vector<ev5_HitsInNthStation>> new_segment;
    std::vector<std::vector<ev5_Segment>> new_diag_segment;
    // Allowed station range
    int min_station = station - 3;
    int max_station = station + 1;
    std::cout<<"max_station/min_station = "<<max_station<<"/"<<min_station<<std::endl;
    // ------------------------------------------------------------
    // Loop over all current segment candidates
    // ------------------------------------------------------------
    for (int i = 0; i < (int)segments.size(); i++) {
        std::vector<ev5_HitsInNthStation> hit_candidates = segments.at(i);
        std::vector<ev5_Segment> hit_diag_candidates = diag_segments.at(i);
        // --------------------------------------------------------
        // Step 1: Check that there are hits in 3 consecutive stations:
        // (station, station-1, station-2)
        // --------------------------------------------------------
        bool hit_exist[3] = {false, false, false};
        for (auto &h : hit_candidates) {
            if (h.station == station)   hit_exist[0] = true;
            if (h.station == station-1) hit_exist[1] = true;
            if (h.station == station-2) hit_exist[2] = true;
        }
        if (!(hit_exist[0] && hit_exist[1] && hit_exist[2])) {
          new_segment.push_back(hit_candidates);
          new_diag_segment.push_back(hit_diag_candidates);
          continue;
        }
        // --------------------------------------------------------
        // Step 2: Fit slope using current diagnostics
        // --------------------------------------------------------
        //obtain the minimum station in the segment
        int max_station = 0;
        int reference_indice = 0; //smallest "z"
        double max_z = -9999; //smallest "z"
        double phi_ref = 0.0;  // reference hit phi
        for(int j=0; j<(int)hit_candidates.size(); j++){
            int station = hit_candidates.at(j).station;
            double x = hit_candidates.at(j).x;
            double y = hit_candidates.at(j).y;
            double z = hit_candidates.at(j).z;
            if(max_station < station) max_station = station;
            if(max_z < z) {
              reference_indice = hit_candidates.at(j).hitIndice;
              phi_ref = atan2(y, x);
              max_z = z;
            }
        }
            std::cout << "reference_indice/phi_ref/max_z " << reference_indice<<"/"<<phi_ref<<"/"<<max_z<<std::endl;
          //recalculate deltaphi in hit_diag_candidates
          for (int j = 0; j < (int)hit_candidates.size(); j++) {
            if(reference_indice == hit_candidates.at(j).hitIndice) {
              hit_diag_candidates.at(j).deltaphi = 0.0;
              std::cout<<"hit_candidates hitIndice/phi/z = "<<hit_candidates.at(j).hitIndice<<"/"<<hit_diag_candidates.at(j).deltaphi<<"/"<<hit_diag_candidates.at(j).z<<std::endl;
              continue;
            }
            double x2 = hit_candidates.at(j).x;
            double y2 = hit_candidates.at(j).y;
            double phi = atan2(y2, x2);
            double dphi = phi - phi_ref;
            // normalize into [-pi, pi]
            if (dphi > M_PI)  dphi -= 2*M_PI;
            if (dphi < -M_PI) dphi += 2*M_PI;
            std::cout << "Hit " << j
                << " phi=" << phi
                << " dphi (relative to ref)=" << dphi << "\n";
            hit_diag_candidates.at(j).deltaphi = dphi;
         }
        double slope_alpha = 0.0;
        double slope_beta  = 0.0;
        double ChiNDF      = 0.0;
        ev5_fit_slope(hit_diag_candidates, slope_alpha, slope_beta, ChiNDF);
        // --------------------------------------------------------
        // Step 3: Extend ONLY to min_station and max_station
        // Update slope every time a hit is added
        // --------------------------------------------------------
        // -- Extend to min_station
       // --------------------------------------------------------
// Step 3: Extend ONLY to min_station
// Update slope every time a hit is added
// --------------------------------------------------------
if (min_station >= 0) {
    std::vector<ev5_HitsInNthStation> Hits_In_Station;
    for (auto &h : _HitsInCluster) {
        if (h.used) continue;
        if (h.station != min_station) continue;
        Hits_In_Station.push_back(h);
    }
    std::cout << "Hits in min station " << min_station << ":\n";
for (const auto &h : Hits_In_Station) {
    std::cout << "hitIndice = " << h.hitIndice
              << ", station = " << h.station
              << ", phi = " << h.phi
              << "\n";
}
    std::cout << std::fixed << std::setprecision(10);
std::cout << "slope_alpha/beta = "
          << slope_alpha << " / " << slope_beta
          << std::endl;
    for (auto &candHit : Hits_In_Station) {
          //recalculate deltaphi in hit_diag_candidates
            double x2 = candHit.x;
            double y2 = candHit.y;
            double phi = atan2(y2, x2);// candHit.phi
            double dphi = phi - phi_ref;
            // normalize into [-pi, pi]
            if (dphi > M_PI)  dphi -= 2*M_PI;
            if (dphi < -M_PI) dphi += 2*M_PI;
            std::cout << " dphi (relative to ref)=" << dphi << "\n";
        ev5_Segment hit;
        hit.deltaphi = dphi;
        hit.z = candHit.z;
        hit.alpha = slope_alpha;
        hit.station = candHit.station;
        hit.usedforfit = false;
        double residual_phi = ev5_ResidualDeltaPhi(slope_alpha, slope_beta,
                                                   max_z, 0.0,
                                                   candHit.z, dphi);
        double thre_residual = 0.3;
         std::cout << " hitIndice/residual_phi/candHit.z = " << candHit.hitIndice<<"/"<<residual_phi << "/"<<candHit.z<< "\n";
        if (residual_phi < thre_residual) {
            // Add hit
            hit_candidates.push_back(candHit);
            hit_diag_candidates.push_back(hit);
            // --- Update slope immediately ---
            /*ev5_fit_slope(hit_diag_candidates, slope_alpha, slope_beta, ChiNDF);
            for (auto &seg2 : hit_diag_candidates) {
                seg2.alpha = slope_alpha;
                seg2.chiNDF = ChiNDF;
            }*/
        }
    }
}
        // -- Extend to max_station
        if (max_station <= 17) {
    std::vector<ev5_HitsInNthStation> Hits_In_Station;
    for (auto &h : _HitsInCluster) {
        if (h.used) continue;
        if (h.station != max_station) continue;
        Hits_In_Station.push_back(h);
    }
    std::cout << "Hits in max station " << max_station << ":\n";
for (const auto &h : Hits_In_Station) {
    std::cout << "hitIndice = " << h.hitIndice
              << ", station = " << h.station
              << ", phi = " << h.phi
              << "\n";
}
    std::cout << std::fixed << std::setprecision(10);
std::cout << "slope_alpha/beta = "
          << slope_alpha << " / " << slope_beta
          << std::endl;
    for (auto &candHit : Hits_In_Station) {
          //recalculate deltaphi in hit_diag_candidates
            double x2 = candHit.x;
            double y2 = candHit.y;
            double phi = atan2(y2, x2);// candHit.phi
            double dphi = phi - phi_ref;
            // normalize into [-pi, pi]
            if (dphi > M_PI)  dphi -= 2*M_PI;
            if (dphi < -M_PI) dphi += 2*M_PI;
            std::cout << " dphi (relative to ref)=" << dphi << "\n";
        ev5_Segment hit;
        hit.deltaphi = dphi;
        hit.z = candHit.z;
        hit.alpha = slope_alpha;
        hit.station = candHit.station;
        hit.usedforfit = false;
        double residual_phi = ev5_ResidualDeltaPhi(slope_alpha, slope_beta,
                                                   max_z, 0.0,
                                                   candHit.z, dphi);
        double thre_residual = 0.3;
        std::cout<<"residual_phi = "<<residual_phi<<std::endl;
        std::cout << "phi = " << phi<<std::endl;
        if (residual_phi < thre_residual) {
            // Add hit
            hit_candidates.push_back(candHit);
            hit_diag_candidates.push_back(hit);
            // --- Update slope immediately ---
            /*ev5_fit_slope(hit_diag_candidates, slope_alpha, slope_beta, ChiNDF);
            for (auto &seg2 : hit_diag_candidates) {
                seg2.alpha = slope_alpha;
                seg2.chiNDF = ChiNDF;
            }*/
        }
    }
        }
        new_segment.push_back(hit_candidates);
        new_diag_segment.push_back(hit_diag_candidates);
    }
    // ------------------------------------------------------------
    // Overwrite the input with the new refined segment lists
    // ------------------------------------------------------------
    ThisIsBestSegment      = new_segment;
    ThisIsBestSegment_Diag = new_diag_segment;
}
//--------------------------------------------------------------------------------//
  void PhiZSeedFinder::ev5_select_best_segments_step_03A(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi){
      std::vector<std::vector<ev5_HitsInNthStation>> segments = ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>> diag_segments = ThisIsBestSegment_Diag;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << " ev5_select_best_segments_step_03A  " << std::endl;
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
      //obtain the minimum station in the segment
      int min_station = 999;
      int max_station = 0;
      int reference_indice = 0; //smallest "z"
      double min_z = 9999; //smallest "z"
      double phi_ref = 0.0;  // reference hit phi
      for(int j=0; j<(int)hit_candidates.size(); j++){
          int station = hit_candidates.at(j).station;
          double x = hit_candidates.at(j).x;
          double y = hit_candidates.at(j).y;
          double z = hit_candidates.at(j).z;
          if(min_station > station) min_station = station;
          if(max_station < station) max_station = station;
          if(min_z > z) {
            reference_indice = hit_candidates.at(j).hitIndice;
            phi_ref = atan2(y, x);
          }
      }
        //recalculate deltaphi in hit_diag_candidates
        for (int j = 0; j < (int)hit_candidates.size(); j++) {
          if(reference_indice == hit_candidates.at(j).hitIndice) hit_diag_candidates.at(j).deltaphi = 0.0;
          double x2 = hit_candidates.at(j).x;
          double y2 = hit_candidates.at(j).y;
          double phi = atan2(y2, x2);
          double dphi = phi - phi_ref;
          // normalize into [-pi, pi]
          if (dphi > M_PI)  dphi -= 2*M_PI;
          if (dphi < -M_PI) dphi += 2*M_PI;
    std::cout << "Hit " << j
              << " phi=" << phi
              << " dphi (relative to ref)=" << dphi << "\n";
          hit_diag_candidates.at(j).deltaphi = dphi;
       }
        // Take combo hits in (n-i)-th stations
        int n = min_station;
        int loop = n;
        int count = 1;
        for(int k=1; 0<=loop-k; k++){
        if(count != k) break;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "            (n-i)-th               " << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "station = " << loop-k<< std::endl;
        double slope_alpha = 0.0;//Slope (a)
        double slope_beta = 0.0;//Intercept (b)
        double ChiNDF = 0.0;
        ev5_fit_slope(hit_diag_candidates, slope_alpha, slope_beta, ChiNDF);
        std::cout<<"slope_alpha/beta = "<<slope_alpha<<"/"<<slope_beta<<std::endl;
        std::vector<ev5_HitsInNthStation> Hits_In_Station;
        Hits_In_Station.clear();
        for(int j=0; j<(int)_HitsInCluster.size(); j++){
            ev5_HitsInNthStation hitsin_nthstation;
            hitsin_nthstation.hitIndice    = _HitsInCluster.at(j).hitIndice;
            hitsin_nthstation.phi          = _HitsInCluster.at(j).phi;
            hitsin_nthstation.strawhits    = _HitsInCluster.at(j).strawhits;
            hitsin_nthstation.x            = _HitsInCluster.at(j).x;
            hitsin_nthstation.y            = _HitsInCluster.at(j).y;
            hitsin_nthstation.z            = _HitsInCluster.at(j).z;
            hitsin_nthstation.station      = _HitsInCluster.at(j).station;
            hitsin_nthstation.plane        = _HitsInCluster.at(j).plane;
            hitsin_nthstation.face        = _HitsInCluster.at(j).face;
            hitsin_nthstation.panel        = _HitsInCluster.at(j).panel;
            hitsin_nthstation.hitID        = _HitsInCluster.at(j).hitID;
            if(_HitsInCluster.at(j).used == true) continue;
            if(loop-k == _HitsInCluster.at(j).station) Hits_In_Station.push_back(hitsin_nthstation);
        }
          size_t nCHsInStn_1 = Hits_In_Station.size();
          if(!((int)nCHsInStn_1 >= 1)) break;
          for(size_t p=0; p<Hits_In_Station.size(); p++){
            std::cout<<"phi/nStrawHits/station/plane/face/panel/x/y/z = "<<Hits_In_Station.at(p).phi<<"/"<<Hits_In_Station.at(p).strawhits<<"/"<<Hits_In_Station.at(p).station<<"/"<<Hits_In_Station.at(p).plane<<"/"<<Hits_In_Station.at(p).face<<"/"<<Hits_In_Station.at(p).panel<<"/"<<Hits_In_Station.at(p).x<<"/"<<Hits_In_Station.at(p).y<<"/"<<Hits_In_Station.at(p).z<<std::endl;
          }
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
              std::cout<<"residual_phi = "<<residual_phi<<std::endl;
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
        loop = max_station;
        int nstation = 18;
        count = 1;
        for(int k=1; loop+k < nstation; k++){
        if(count != k) break;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "            (n+i)-th               " << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "station = " << loop+k  << std::endl;
        double slope_alpha = 0.0;//Slope (a)
        double slope_beta = 0.0;//Intercept (b)
        double ChiNDF = 0.0;
        ev5_fit_slope(hit_diag_candidates, slope_alpha, slope_beta, ChiNDF);
        std::cout<<"slope_alpha/beta = "<<slope_alpha<<"/"<<slope_beta<<std::endl;
        std::vector<ev5_HitsInNthStation> Hits_In_Station;
        Hits_In_Station.clear();
        for(int j=0; j<(int)_HitsInCluster.size(); j++){
            ev5_HitsInNthStation hitsin_nthstation;
            hitsin_nthstation.hitIndice    = _HitsInCluster.at(j).hitIndice;
            hitsin_nthstation.phi          = _HitsInCluster.at(j).phi;
            hitsin_nthstation.strawhits    = _HitsInCluster.at(j).strawhits;
            hitsin_nthstation.x            = _HitsInCluster.at(j).x;
            hitsin_nthstation.y            = _HitsInCluster.at(j).y;
            hitsin_nthstation.z            = _HitsInCluster.at(j).z;
            hitsin_nthstation.station      = _HitsInCluster.at(j).station;
            hitsin_nthstation.plane        = _HitsInCluster.at(j).plane;
            hitsin_nthstation.face         = _HitsInCluster.at(j).face;
            hitsin_nthstation.panel        = _HitsInCluster.at(j).panel;
            hitsin_nthstation.hitID        = _HitsInCluster.at(j).hitID;
            if(_HitsInCluster.at(j).used == true) continue;
            if(loop+k == _HitsInCluster.at(j).station) Hits_In_Station.push_back(hitsin_nthstation);
        }
          size_t nCHsInStn_1 = Hits_In_Station.size();
          if(!((int)nCHsInStn_1 >= 1)) break;
          //std::cout << "loop/k " << loop<<"/"<<k<< std::endl;
          for(size_t p=0 ; p<Hits_In_Station.size(); p++){
            std::cout<<"phi/nStrawHits/station/plane/face/panel/x/y/z = "<<Hits_In_Station.at(p).phi<<"/"<<Hits_In_Station.at(p).strawhits<<"/"<<Hits_In_Station.at(p).station<<"/"<<Hits_In_Station.at(p).plane<<"/"<<Hits_In_Station.at(p).face<<"/"<<Hits_In_Station.at(p).panel<<"/"<<Hits_In_Station.at(p).x<<"/"<<Hits_In_Station.at(p).y<<"/"<<Hits_In_Station.at(p).z<<std::endl;
          }
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
  }//end ev5_select_best_segments_step_03A
//--------------------------------------------------------------------------------//
void PhiZSeedFinder::ev5_select_best_segments_step_030(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag){
      std::vector<std::vector<ev5_HitsInNthStation>> segments = ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>> diag_segments = ThisIsBestSegment_Diag;
    std::cout << "-----------------------------------" << std::endl;
    std::cout << " ev5_select_best_segments_step_030 " << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    // ---------------------------------------------------------------------
    // Step 1: Build a mapping from "hit index" -> "which segments contain it"
    // ---------------------------------------------------------------------
    std::unordered_map<int, std::set<int>> hitToSegments;
    std::cout << "Step 1: Building hit -> segment mapping\n";
    for (size_t i = 0; i < segments.size(); ++i) {
        std::cout << " Segment " << i << " contains hits: ";
        for (auto& h : segments[i]) {
            std::cout << h.hitIndice << " ";
            hitToSegments[h.hitIndice].insert(i);
        }
        std::cout << "\n";
    }
    std::cout << "Hit -> segment map:\n";
    for (auto& kv : hitToSegments) {
        std::cout << " Hit " << kv.first << " in segments: ";
        for (int segIdx : kv.second) std::cout << segIdx << " ";
        std::cout << "\n";
    }
    // ---------------------------------------------------------------------
    // Step 2: For each shared hit, decide ownership
    // ---------------------------------------------------------------------
    for (auto& kv : hitToSegments) {
        if (kv.second.size() < 2) continue;  // skip non-shared hits
        int hitIdx = kv.first;
        auto& segSet = kv.second;
        std::cout << "\nResolving shared hit " << hitIdx << " present in segments: ";
        for (int segIdx : segSet) std::cout << segIdx << " ";
        std::cout << "\n";
        double bestChi2Diff = -1e12; // start very negative
        int bestSeg = -1;
        for (int segIdx : segSet) {
            auto& seg = segments[segIdx];
            // --- find the hit inside this segment ---
            auto it = std::find_if(seg.begin(), seg.end(),
                                   [&](auto& h){ return h.hitIndice == hitIdx; });
            if (it == seg.end()) continue;
            auto backupHit = *it;
            // --- compute chi2 WITH hit ---
            double chi2ndf_with = 0.0;
            findchisq(seg, chi2ndf_with);
            // --- temporarily remove the hit ---
            seg.erase(it);
            // --- compute chi2 WITHOUT hit ---
            double chi2ndf_without = 0.0;
            findchisq(seg, chi2ndf_without);
            // --- restore hit ---
            seg.push_back(backupHit);
            double chi2diff = chi2ndf_without - chi2ndf_with;
            std::cout << "  Segment " << segIdx
                      << " chi2/ndf WITH hit=" << chi2ndf_with
                      << ", WITHOUT hit=" << chi2ndf_without
                      << ", diff=" << chi2diff << "\n";
            if (chi2diff > bestChi2Diff) {
                bestChi2Diff = chi2diff;
                bestSeg = segIdx;
            }
        }
        std::cout << " Best segment for hit " << hitIdx << " is segment " << bestSeg
                  << " (max chi2 improvement = " << bestChi2Diff << ")\n";
        // ------------------------------------------------------------------
        // Step 3: Remove the hit from all non-best segments
        // ------------------------------------------------------------------
        for (int segIdx : segSet) {
            if (segIdx == bestSeg) continue;
            auto& seg = segments[segIdx];
            auto& diag = diag_segments[segIdx];
            seg.erase(std::remove_if(seg.begin(), seg.end(),
                                     [&](auto& h){ return h.hitIndice == hitIdx; }),
                      seg.end());
            diag.erase(std::remove_if(diag.begin(), diag.end(),
                                      [&](auto& d){ return d.reference_point == hitIdx; }),
                       diag.end());
            std::cout << " Removed hit " << hitIdx << " from segment " << segIdx << "\n";
        }
    }
    // ------------------------------------------------------------------
    // Step 4:  Remove small or empty segments
    // ------------------------------------------------------------------
    for (int k = segments.size() - 1; k >= 0; --k) { // iterate backward safely
        if (segments[k].size() < 3) {
            segments.erase(segments.begin() + k);
            diag_segments.erase(diag_segments.begin() + k);
        }
    }
      std::cout<<"nsegment = "<<segments.size()<<std::endl;
      std::cout<<"n diag_segments = "<<diag_segments.size()<<std::endl;
    // ------------------------------------------------------------------
    // Step 5: Flag hits in _HitsInCluster that are already used in segments
    // ------------------------------------------------------------------
    /*for (size_t i = 0; i < segments.size(); ++i) {
    for (size_t j = 0; j < segments[i].size(); ++j) {
        int segHitIndex = segments[i][j].hitIndice;  // hit index in this segment
        for (size_t k = 0; k < _HitsInCluster.size(); ++k) {
            if (_HitsInCluster[k].hitIndice == segHitIndex) {
                _HitsInCluster[k].used = true; // mark as used
            }
        }
    }
    }*/
    std::cout << "-----------------------------------\n";
    std::cout << "         _HitsInCluster dump       \n";
    std::cout << "-----------------------------------\n";
    for (size_t i = 0; i < _HitsInCluster.size(); ++i) {
        const auto& h = _HitsInCluster[i];
        std::cout << "Hit #" << i
                  << " (hitIndice=" << h.hitIndice << ") "
                  << "phi=" << h.phi
                  << " xyz=(" << h.x << ", " << h.y << ", " << h.z << ") "
                  << "station=" << h.station
                  << " plane=" << h.plane
                  << " face=" << h.face
                  << " panel=" << h.panel
                  << " segmentIndex=" << h.segmentIndex
                  << " used=" << (h.used ? "true" : "false")
                  << "\n";
    }
    //Re-Fill
    ThisIsBestSegment.clear();
    ThisIsBestSegment_Diag.clear();
    ThisIsBestSegment = segments;
    ThisIsBestSegment_Diag = diag_segments;
    std::cout<<"ThisIsBestSegment = "<<ThisIsBestSegment.size()<<std::endl;
    std::cout<<"ThisIsBestSegment_Diag = "<<ThisIsBestSegment_Diag.size()<<std::endl;
  }//end ev5_select_best_segments_step_030
//--------------------------------------------------------------------------------//
void PhiZSeedFinder::ev5_select_best_segments_step_06A(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag){
      std::vector<std::vector<ev5_HitsInNthStation>> segments = all_ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>> diag_segments = all_ThisIsBestSegment_Diag;
    std::cout << "-----------------------------------" << std::endl;
    std::cout << " ev5_select_best_segments_step_06A " << std::endl;
    std::cout << "-----------------------------------" << std::endl;
// ---------------------------------------------------------------------
// Step 1: Build a mapping from "hit index" -> "which segments contain it"
// ---------------------------------------------------------------------
std::unordered_map<int, std::set<int>> hitToSegments;
std::cout << "Step 1: Building hit -> segment mapping\n";
for (size_t i = 0; i < segments.size(); ++i) {
    std::cout << " Segment " << i << " contains hits: ";
    for (auto& h : segments[i]) {
        std::cout << h.hitIndice << " ";
        hitToSegments[h.hitIndice].insert(i);
    }
    std::cout << "\n";
}
std::cout << "Hit -> segment map:\n";
for (auto& kv : hitToSegments) {
    std::cout << " Hit " << kv.first << " in segments: ";
    for (int segIdx : kv.second) std::cout << segIdx << " ";
    std::cout << "\n";
}
// ---------------------------------------------------------------------
// Step 2: For each shared hit, decide ownership
// ---------------------------------------------------------------------
for (auto& kv : hitToSegments) {
    if (kv.second.size() < 2) continue;  // skip non-shared hits
    int hitIdx = kv.first;
    auto& segSet = kv.second;
    std::cout << "\nResolving shared hit " << hitIdx << " present in segments: ";
    for (int segIdx : segSet) std::cout << segIdx << " ";
    std::cout << "\n";
    double bestChi2Diff = -1e12; // start very negative
    int bestSeg = -1;
    for (int segIdx : segSet) {
        auto& seg = segments[segIdx];
        // --- find the hit inside this segment ---
        auto it = std::find_if(seg.begin(), seg.end(),
                               [&](auto& h){ return h.hitIndice == hitIdx; });
        if (it == seg.end()) continue;
        auto backupHit = *it;
        // --- compute chi2 WITH hit ---
        double chi2ndf_with = 0.0;
        findchisq(seg, chi2ndf_with);
        // --- temporarily remove the hit ---
        seg.erase(it);
        // --- compute chi2 WITHOUT hit ---
        double chi2ndf_without = 0.0;
        findchisq(seg, chi2ndf_without);
        // --- restore hit ---
        seg.push_back(backupHit);
        double chi2diff = chi2ndf_without - chi2ndf_with;
        std::cout << "  Segment " << segIdx
                  << " chi2/ndf WITH hit=" << chi2ndf_with
                  << ", WITHOUT hit=" << chi2ndf_without
                  << ", diff=" << chi2diff << "\n";
        if (chi2diff > bestChi2Diff) {
            bestChi2Diff = chi2diff;
            bestSeg = segIdx;
        }
    }
    std::cout << " Best segment for hit " << hitIdx << " is segment " << bestSeg
              << " (max chi2 improvement = " << bestChi2Diff << ")\n";
    // ------------------------------------------------------------------
    // Step 3: Remove the hit from all non-best segments
    // ------------------------------------------------------------------
    for (int segIdx : segSet) {
        if (segIdx == bestSeg) continue;
        auto& seg = segments[segIdx];
        auto& diag = diag_segments[segIdx];
        seg.erase(std::remove_if(seg.begin(), seg.end(),
                                 [&](auto& h){ return h.hitIndice == hitIdx; }),
                  seg.end());
        diag.erase(std::remove_if(diag.begin(), diag.end(),
                                  [&](auto& d){ return d.reference_point == hitIdx; }),
                   diag.end());
        std::cout << " Removed hit " << hitIdx << " from segment " << segIdx << "\n";
    }
}
/*
    // ---------------------------------------------------------------------
    // Step 1: Build a mapping from "hit index" -> "which segments contain it"
    // ---------------------------------------------------------------------
    std::unordered_map<int, std::set<int>> hitToSegments;
    std::cout << "Step 1: Building hit -> segment mapping\n";
    for (size_t i = 0; i < segments.size(); ++i) {
        std::cout << " Segment " << i << " contains hits: ";
        for (auto& h : segments[i]) {
            std::cout << h.hitIndice << " ";
            hitToSegments[h.hitIndice].insert(i);
        }
        std::cout << "\n";
    }
    std::cout << "Hit -> segment map:\n";
    for (auto& kv : hitToSegments) {
        std::cout << " Hit " << kv.first << " in segments: ";
        for (int segIdx : kv.second) std::cout << segIdx << " ";
        std::cout << "\n";
    }
// ---------------------------------------------------------------------
// Step 2: For each shared hit, decide ownership
// ---------------------------------------------------------------------
double chi2diffThreshold = 0.005;  // <-- require chi² improvement > 3
for (std::unordered_map<int, std::set<int>>::iterator kv = hitToSegments.begin();
     kv != hitToSegments.end(); ++kv) {
    if (kv->second.size() < 2) continue;  // skip non-shared hits
    int hitIdx = kv->first;
    std::set<int>& segSet = kv->second;
    std::cout << "\nResolving shared hit " << hitIdx << " present in segments: ";
    for (std::set<int>::iterator it = segSet.begin(); it != segSet.end(); ++it) {
        std::cout << *it << " ";
    }
    std::cout << "\n";
    double bestChi2Diff = -1e12;
    int bestSeg = -1;
    // Step 2a: check chi² effect of removing this hit
    for (std::set<int>::iterator it = segSet.begin(); it != segSet.end(); ++it) {
        int segIdx = *it;
        std::vector<ev5_HitsInNthStation>& seg = segments[segIdx];
        // --- find the hit inside this segment ---
        std::vector<ev5_HitsInNthStation>::iterator hitIt =
            std::find_if(seg.begin(), seg.end(),
                         [hitIdx](const ev5_HitsInNthStation& h){ return h.hitIndice == hitIdx; });
        if (hitIt == seg.end()) continue;
        ev5_HitsInNthStation backupHit = *hitIt;
        // --- compute chi2 WITH hit ---
        double chi2ndf_with = 0.0;
        findchisq(seg, chi2ndf_with);
        // --- temporarily remove the hit ---
        seg.erase(hitIt);
        // --- compute chi2 WITHOUT hit ---
        double chi2ndf_without = 0.0;
        findchisq(seg, chi2ndf_without);
        // --- restore hit ---
        seg.push_back(backupHit);
        double chi2diff = chi2ndf_without - chi2ndf_with;
        std::cout << "  Segment " << segIdx
                  << " chi2/ndf WITH hit=" << chi2ndf_with
                  << ", WITHOUT hit=" << chi2ndf_without
                  << ", diff=" << chi2diff << "\n";
        if (chi2diff > bestChi2Diff) {
            bestChi2Diff = chi2diff;
            bestSeg = segIdx;
        }
    }
    // ------------------------------------------------------------------
    // Step 3: Remove the hit only if chi² diff is significant
    // ------------------------------------------------------------------
    if (bestChi2Diff > chi2diffThreshold && bestSeg >= 0) {
        std::cout << " Best segment for hit " << hitIdx << " is segment " << bestSeg
                  << " (max chi2 improvement = " << bestChi2Diff << ")\n";
        for (std::set<int>::iterator it = segSet.begin(); it != segSet.end(); ++it) {
            int segIdx = *it;
            if (segIdx == bestSeg) continue;
            std::vector<ev5_HitsInNthStation>& seg = segments[segIdx];
            std::vector<ev5_Segment>& diag = diag_segments[segIdx];
            seg.erase(std::remove_if(seg.begin(), seg.end(),
                                     [hitIdx](const ev5_HitsInNthStation& h){ return h.hitIndice == hitIdx; }),
                      seg.end());
            diag.erase(std::remove_if(diag.begin(), diag.end(),
                                      [hitIdx](const ev5_Segment& d){ return d.reference_point == hitIdx; }),
                       diag.end());
            std::cout << " Removed hit " << hitIdx << " from segment " << segIdx << "\n";
        }
    } else {
        std::cout << " Hit " << hitIdx
                  << " kept in all segments (chi2 improvement too small: "
                  << bestChi2Diff << ")\n";
    }
}
*/
    // ------------------------------------------------------------------
    // Step 4:  Remove small or empty segments
    // ------------------------------------------------------------------
    for (int k = segments.size() - 1; k >= 0; --k) { // iterate backward safely
        if (segments[k].size() < 3) {
            segments.erase(segments.begin() + k);
            diag_segments.erase(diag_segments.begin() + k);
        }
    }
      std::cout<<"nsegment = "<<segments.size()<<std::endl;
      std::cout<<"n diag_segments = "<<diag_segments.size()<<std::endl;
    // ------------------------------------------------------------------
    // Step 5: Flag hits in _HitsInCluster that are already used in segments
    // ------------------------------------------------------------------
    /*for (size_t i = 0; i < segments.size(); ++i) {
    for (size_t j = 0; j < segments[i].size(); ++j) {
        int segHitIndex = segments[i][j].hitIndice;  // hit index in this segment
        for (size_t k = 0; k < _HitsInCluster.size(); ++k) {
            if (_HitsInCluster[k].hitIndice == segHitIndex) {
                _HitsInCluster[k].used = true; // mark as used
            }
        }
    }
    }*/
    std::cout << "-----------------------------------\n";
    std::cout << "         _HitsInCluster dump       \n";
    std::cout << "-----------------------------------\n";
    for (size_t i = 0; i < _HitsInCluster.size(); ++i) {
        const auto& h = _HitsInCluster[i];
        std::cout << "Hit #" << i
                  << " (hitIndice=" << h.hitIndice << ") "
                  << "phi=" << h.phi
                  << " xyz=(" << h.x << ", " << h.y << ", " << h.z << ") "
                  << "station=" << h.station
                  << " plane=" << h.plane
                  << " face=" << h.face
                  << " panel=" << h.panel
                  << " segmentIndex=" << h.segmentIndex
                  << " used=" << (h.used ? "true" : "false")
                  << "\n";
    }
    //Re-Fill
    all_ThisIsBestSegment.clear();
    all_ThisIsBestSegment_Diag.clear();
    all_ThisIsBestSegment = segments;
    all_ThisIsBestSegment_Diag = diag_segments;
    std::cout<<"ThisIsBestSegment = "<<all_ThisIsBestSegment.size()<<std::endl;
    std::cout<<"ThisIsBestSegment_Diag = "<<all_ThisIsBestSegment_Diag.size()<<std::endl;
  }//end ev5_select_best_segments_step_030
//--------------------------------------------------------------------------------//
  void PhiZSeedFinder::ev5_select_best_segments_step_031(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi){
      std::vector<std::vector<ev5_HitsInNthStation>> segments = ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>> diag_segments = ThisIsBestSegment_Diag;
      std::cout << "-----------------------------------" << std::endl;
      std::cout << "-----------------------------------" << std::endl;
      std::cout << " ev5_select_best_segments_step_031  " << std::endl;
      std::cout << "-----------------------------------" << std::endl;
      std::cout << "-----------------------------------" << std::endl;
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
      int max_station = 0;
      /*for(int i=0; i<(int)segments.size(); i++){
        for(int j=0; j<(int)segments.at(i).size(); j++){
          int station = segments.at(i).at(j).station;
          if(min_station > station) min_station = station;
          if(max_station < station) max_station = station;
        }
      }*/
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
        min_station = 18;
        max_station = 0;
        for(int j=0; j<(int)segments.at(i).size(); j++){
          hit_candidates.push_back(segments.at(i).at(j));
          hit_diag_candidates.push_back(diag_segments.at(i).at(j));
          if (min_station > segments.at(i).at(j).station) min_station = segments.at(i).at(j).station;
          if (max_station < segments.at(i).at(j).station) max_station = segments.at(i).at(j).station;
        }
        // Take combo hits in (n-i)-th stations
        int n = min_station;
        int loop = n;
        int count = 1;
        std::cout<<"min_station = "<<min_station<<std::endl;
        std::cout<<"max_station = "<<max_station<<std::endl;
        for(int k=1; 0<=loop-k; k++){
        int hit_yes = 0;
        if(count == 2) break;//gap only 1
        //std::cout << "-----------------------------------" << std::endl;
        //std::cout << "            (n-i)-th               " << std::endl;
        //std::cout << "-----------------------------------" << std::endl;
        std::cout << "station = " << loop-k<<std::endl;
        double slope_alpha = 0.0;//Slope (a)
        double slope_beta = 0.0;//Intercept (b)
        double ChiNDF = 0.0;
        ev5_fit_slope(hit_diag_candidates, slope_alpha, slope_beta, ChiNDF);
        //std::cout<<"slope_alpha/beta = "<<slope_alpha<<"/"<<slope_beta<<std::endl;
        std::vector<ev5_HitsInNthStation> Hits_In_Station;
        Hits_In_Station.clear();
        for(int j=0; j<(int)_HitsInCluster.size(); j++){
            ev5_HitsInNthStation hitsin_nthstation;
            hitsin_nthstation.hitIndice    = _HitsInCluster.at(j).hitIndice;
            hitsin_nthstation.phi          = _HitsInCluster.at(j).phi;
            hitsin_nthstation.strawhits    = _HitsInCluster.at(j).strawhits;
            hitsin_nthstation.x            = _HitsInCluster.at(j).x;
            hitsin_nthstation.y            = _HitsInCluster.at(j).y;
            hitsin_nthstation.z            = _HitsInCluster.at(j).z;
            hitsin_nthstation.station      = _HitsInCluster.at(j).station;
            hitsin_nthstation.plane        = _HitsInCluster.at(j).plane;
            hitsin_nthstation.face        = _HitsInCluster.at(j).face;
            hitsin_nthstation.panel        = _HitsInCluster.at(j).panel;
            hitsin_nthstation.hitID        = _HitsInCluster.at(j).hitID;
            if(_HitsInCluster.at(j).used == true) continue;
            if(loop-k == _HitsInCluster.at(j).station) Hits_In_Station.push_back(hitsin_nthstation);
        }
        size_t nCHsInStn_1 = Hits_In_Station.size();
        if(nCHsInStn_1 == 0) {
            count++;
            continue;
        }
        std::cout<<"nCHsInStn_1 = "<<nCHsInStn_1<<std::endl;
          for(size_t p=0; p<Hits_In_Station.size(); p++){
            std::cout<<"phi/nStrawHits/station/plane/face/panel/x/y/z = "<<Hits_In_Station.at(p).phi<<"/"<<Hits_In_Station.at(p).strawhits<<"/"<<Hits_In_Station.at(p).station<<"/"<<Hits_In_Station.at(p).plane<<"/"<<Hits_In_Station.at(p).face<<"/"<<Hits_In_Station.at(p).panel<<"/"<<Hits_In_Station.at(p).x<<"/"<<Hits_In_Station.at(p).y<<"/"<<Hits_In_Station.at(p).z<<std::endl;
          }
          /*std::vector<ev5_HitsInNthStation> LeftHit = segments.at(i);
          // re-sort vector in ascending order of z-coordinate
          std::sort(LeftHit.begin(), LeftHit.end(), [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b){
            return a.z < b.z;
          });*/
          //Check the left most hits in the segment
          double first_Hit[4] = {hit_candidates.at(0).x, hit_candidates.at(0).y, hit_candidates.at(0).z, hit_candidates.at(0).phi};
          //Check hits in (n-i)th station, if they configure the segment
          //For (n-i)th station
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
          if(hit_yes == 1) count = 0;
          else count++;
        }//end (n-i)-th stations
        // Take combo hits in (n+i)-th stations
        loop = max_station;
        int nstation = 18;
        count = 0;
        for(int k=1; loop+k < nstation; k++){
        std::cout<<"k/count = "<<k<<"/"<<count<<std::endl;
        int hit_yes = 0;
        if(count == 2) break;//gap only 1
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "            (n+i)-th               " << std::endl;
        std::cout << "-----------------------------------" << std::endl;
        std::cout << "station = " << loop+k  << std::endl;
        double slope_alpha = 0.0;//Slope (a)
        double slope_beta = 0.0;//Intercept (b)
        double ChiNDF = 0.0;
        ev5_fit_slope(hit_diag_candidates, slope_alpha, slope_beta, ChiNDF);
        //std::cout<<"slope_alpha/beta = "<<slope_alpha<<"/"<<slope_beta<<std::endl;
        std::vector<ev5_HitsInNthStation> Hits_In_Station;
        Hits_In_Station.clear();
        for(int j=0; j<(int)_HitsInCluster.size(); j++){
            ev5_HitsInNthStation hitsin_nthstation;
            hitsin_nthstation.hitIndice    = _HitsInCluster.at(j).hitIndice;
            hitsin_nthstation.phi          = _HitsInCluster.at(j).phi;
            hitsin_nthstation.strawhits    = _HitsInCluster.at(j).strawhits;
            hitsin_nthstation.x            = _HitsInCluster.at(j).x;
            hitsin_nthstation.y            = _HitsInCluster.at(j).y;
            hitsin_nthstation.z            = _HitsInCluster.at(j).z;
            hitsin_nthstation.station      = _HitsInCluster.at(j).station;
            hitsin_nthstation.plane        = _HitsInCluster.at(j).plane;
            hitsin_nthstation.face        = _HitsInCluster.at(j).face;
            hitsin_nthstation.panel        = _HitsInCluster.at(j).panel;
            hitsin_nthstation.hitID        = _HitsInCluster.at(j).hitID;
            if(_HitsInCluster.at(j).used == true) continue;
            if(loop+k == _HitsInCluster.at(j).station) Hits_In_Station.push_back(hitsin_nthstation);
        }
          size_t nCHsInStn_1 = Hits_In_Station.size();
          if(nCHsInStn_1 == 0) {
            count++;
            continue;
          }
          //std::cout << "loop/k " << loop<<"/"<<k<< std::endl;
          /*for(size_t p=0 ; p<Hits_In_Station.size(); p++){
            std::cout<<"phi/nStrawHits/station/plane/face/panel/x/y/z = "<<Hits_In_Station.at(p).phi<<"/"<<Hits_In_Station.at(p).strawhits<<"/"<<Hits_In_Station.at(p).station<<"/"<<Hits_In_Station.at(p).plane<<"/"<<Hits_In_Station.at(p).face<<"/"<<Hits_In_Station.at(p).panel<<"/"<<Hits_In_Station.at(p).x<<"/"<<Hits_In_Station.at(p).y<<"/"<<Hits_In_Station.at(p).z<<std::endl;
          }*/
          /*std::vector<ev5_HitsInNthStation> RightHit = segments.at(i);
          // re-sort vector in ascending order of z-coordinate
          std::sort(RightHit.begin(), RightHit.end(), [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b){
            return a.z < b.z;
          });*/
          //Check the left most hits in the segment
          //int index = (int)RightHit.size()-1;
          //double first_Hit[4] = {RightHit.at(0).x, RightHit.at(0).y, RightHit.at(0).z, RightHit.at(0).phi};
          double first_Hit[4] = {hit_candidates.at(0).x, hit_candidates.at(0).y, hit_candidates.at(0).z, hit_candidates.at(0).phi};
          //Check hits in (n-i)th station, if they configure the segment
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
          if(hit_yes == 1) count = 0;
          else count++;
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
  }//end ev5_select_best_segments_step_031
//--------------------------------------------------------------------------------//
/*    void PhiZSeedFinder::ev5_fit_slope(const std::vector<ev5_Segment>& hit_diag, double& alpha, double& beta, double& chindf){
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
      //std::cout << "Fit Results:" << std::endl;
      //std::cout << "Slope (a): " << fitFunc->GetParameter(0) << " +/- " << fitFunc->GetParError(0) << std::endl;
      //std::cout << "Intercept (b): " << fitFunc->GetParameter(1) << " +/- " << fitFunc->GetParError(1) << std::endl;
      //std::cout << "Chi-square *100.: " << fitFunc->GetChisquare()*100.0 << std::endl;
      //std::cout << "NDF: " << fitFunc->GetNDF() << std::endl;
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
*/
//--------------------------------------------------------------------------------//
    void PhiZSeedFinder::ev5_fit_slope(const std::vector<ev5_Segment>& hit_diag, double& alpha, double& beta, double& chindf){
      ::LsqSums2 fitter;
      fitter.clear();
      for(int j=0; j<(int)hit_diag.size(); j++){
        // z as x-axis and phi as y-axis
        double z = hit_diag.at(j).z;
        double phi = hit_diag.at(j).deltaphi;
        double seedError2 = 0.1;
        double seedWeight = 1.0 / (seedError2);
        std::cout<<"j/z/phi = "<<j<<"/"<<z<<"/"<<phi<<std::endl;
        fitter.addPoint(z, phi, seedWeight);
      }
      //return fitting values
      double dphidz = fitter.dydx();
      alpha = dphidz;
      beta = fitter.y0();
      chindf = fitter.chi2Dof();
  }
//--------------------------------------------------------------------------------//
    void PhiZSeedFinder::ev5_fit_slope_ver2(int i, double& alpha, double& beta, double& chindf){
    //Get the slope value from the 1st segment
    //continue collecting hits until station gap is > 1
    int Reference_HitIndex = 0;// hit should be located in the upstream tracker
    int Reference_HitIndice = -1;// hit should be located in the upstream tracker
    //int Reference_station = 0;
    double Reference_Phi = -999.9;
    for(size_t j=0; j<_segmentHits.at(i).size(); j++){
      Reference_HitIndex = j;
      Reference_HitIndice = _segmentHits.at(i).at(j).hitIndice;
      //Reference_station = _segmentHits.at(i).at(j).station;
      Reference_Phi = _segmentHits.at(i).at(j).phi;
      _segmentHits.at(i).at(j).phiDiag = Reference_Phi;
      /*std::cout<<"best_segments[i][j].Phi = "<<_segmentHits.at(i).at(j).phi<<std::endl;
      std::cout<<"best_segments[i][j].ambigPhi = "<<_segmentHits.at(i).at(j).ambigPhi<<std::endl;
      std::cout<<"station_z = "<<_segmentHits.at(i).at(j).z<<std::endl;
      std::cout<<"station_Reference = "<<_segmentHits.at(i).at(j).station<<std::endl;
      std::cout<<"Reference_Index = "<<Reference_HitIndex<<std::endl;
      std::cout<<"Reference_Indice = "<<Reference_HitIndice<<std::endl;
      std::cout<<"Reference_Phi = "<<Reference_Phi<<std::endl;*/
      break;
    }
    ::LsqSums2 _lineFitter;
    _lineFitter.clear();
    for(size_t j=0; j<_segmentHits.at(i).size(); j++){
      double z = _segmentHits.at(i).at(j).z;
      double phiWeight = 0.1;
      //only add reference hit to the fiiting function
      if(_segmentHits.at(i).at(j).hitIndice == Reference_HitIndice) {
        _lineFitter.addPoint(z, Reference_Phi, phiWeight);
      }
      // get other hits
      if(Reference_HitIndex == (int)i) continue;
      if(_segmentHits.at(i).at(j).hitIndice == Reference_HitIndice) continue;
      //coorect helixphi and consider 2pi boundary
      float deltaPhi = _segmentHits.at(i).at(j).phi - Reference_Phi;
      /*std::cout<<"hitIndice = "<<_tcHits[i].hitIndice<<std::endl;
      std::cout<<"Helixphi = "<<_tcHits[i].helixPhi<<std::endl;
      std::cout<<"Z = "<<z<<std::endl;
      std::cout<<"station = "<<_tcHits[i].station<<std::endl;
      std::cout<<"phi = "<<_tcHits[i].helixPhi<<std::endl;
      std::cout<<"deltaPhi = "<<deltaPhi<<std::endl;
      */
      // If it turns more than pi then, consinder the 2pi boundary
      int turns = 0;
      if (deltaPhi > M_PI) turns--;
      if (deltaPhi < -M_PI) turns++;
      double phi = _segmentHits.at(i).at(j).phi + turns * 2 * M_PI;
      _segmentHits.at(i).at(j).phiDiag = phi;
      //std::cout<<"_tcHits[i].ambigPhi = "<<phi<<std::endl;
      // quality cut for the 1st segment
      if(_lineFitter.qn() <= 2) {
        _lineFitter.addPoint(z, _segmentHits.at(i).at(j).phiDiag, phiWeight);
        continue;
      }
      if(_lineFitter.qn() > 2) {
      if(turns != 0){
          // corss check
          double lineSlope = _lineFitter.dydx();
          double lineIntercept = _lineFitter.y0();
          // Predict phi from the line
          double predictedPhi = lineSlope * z + lineIntercept;
          // Compute the difference between prediction and actual
          double diffPhi[2] = {0.0};
          diffPhi[0] = predictedPhi - _segmentHits.at(i).at(j).phiDiag;
          diffPhi[1] = predictedPhi - _segmentHits.at(i).at(j).phi;
          // choose the nearest assumption
          if(abs(diffPhi[1]) < abs(diffPhi[0])) turns = 0;
          phi = _segmentHits.at(i).at(j).phi + turns * 2 * M_PI;
          _segmentHits.at(i).at(j).phiDiag = phi;
          // Round delta/2π to nearest integer for wrapping correction
          _lineFitter.addPoint(z, _segmentHits.at(i).at(j).phiDiag, phiWeight);
      }
        else {
          _lineFitter.addPoint(z, _segmentHits.at(i).at(j).phiDiag, phiWeight);
        }
      }
  }
      //return fitting values
      double dphidz = _lineFitter.dydx();
      alpha = dphidz;
      beta = _lineFitter.y0();
      chindf = _lineFitter.chi2Dof();
  }
//--------------------------------------------------------------------------------//
    void PhiZSeedFinder::ev5_fit_slope_ver3(int i, double& alpha, double& beta, double& chindf){
    //Get the slope value from the 1st segment
    //continue collecting hits until station gap is > 1
    int Reference_HitIndex = 0;// hit should be located in the upstream tracker
    //int Reference_HitIndice = -1;// hit should be located in the upstream tracker
    //int Reference_station = 0;
    double Reference_Phi = -999.9;
    for(size_t j=0; j<_segmentHits.at(i).size(); j++){
      Reference_HitIndex = j;
      //Reference_HitIndice = _segmentHits.at(i).at(j).hitIndice;
      //Reference_station = _segmentHits.at(i).at(j).station;
      Reference_Phi = _segmentHits.at(i).at(j).phi;
      _segmentHits.at(i).at(j).phiDiag = Reference_Phi;
      /*std::cout<<"best_segments[i][j].Phi = "<<_segmentHits.at(i).at(j).phi<<std::endl;
      std::cout<<"best_segments[i][j].ambigPhi = "<<_segmentHits.at(i).at(j).ambigPhi<<std::endl;
      std::cout<<"station_z = "<<_segmentHits.at(i).at(j).z<<std::endl;
      std::cout<<"station_Reference = "<<_segmentHits.at(i).at(j).station<<std::endl;
      std::cout<<"Reference_Index = "<<Reference_HitIndex<<std::endl;
      std::cout<<"Reference_Indice = "<<Reference_HitIndice<<std::endl;
      std::cout<<"Reference_Phi = "<<Reference_Phi<<std::endl;*/
      break;
    }
    ::LsqSums2 _lineFitter;
    _lineFitter.clear();
      //only add reference hit to the fiiting function
    _lineFitter.addPoint(_segmentHits.at(i).at(Reference_HitIndex).z, _segmentHits.at(i).at(Reference_HitIndex).phi, 0.1);
    for(size_t j=0; j<_segmentHits.at(i).size(); j++){
      double z = _segmentHits.at(i).at(j).z;
      double phiError2 = 0.01;
      for (size_t k = 0; k < _tcHits.size(); k++){
        if(_segmentHits.at(i).at(j).hitIndice != _tcHits[k].hitIndice) continue;
        phiError2 = _tcHits[k].helixPhiError2;
      }
      double phiWeight = 1.0/phiError2;
      // get other hits
      if(Reference_HitIndex == (int)j) continue;
      //coorect helixphi and consider 2pi boundary
      float deltaPhi = _segmentHits.at(i).at(j).phi - Reference_Phi;
      /*std::cout<<"hitIndice = "<<_tcHits[i].hitIndice<<std::endl;
      std::cout<<"Helixphi = "<<_tcHits[i].helixPhi<<std::endl;
      std::cout<<"Z = "<<z<<std::endl;
      std::cout<<"station = "<<_tcHits[i].station<<std::endl;
      std::cout<<"phi = "<<_tcHits[i].helixPhi<<std::endl;
      std::cout<<"deltaPhi = "<<deltaPhi<<std::endl;
      */
      // If it turns more than pi then, consinder the 2pi boundary
      int turns = 0;
      if (deltaPhi > M_PI) turns--;
      if (deltaPhi < -M_PI) turns++;
      double phi = _segmentHits.at(i).at(j).phi + turns * 2 * M_PI;
      _segmentHits.at(i).at(j).phiDiag = phi;
      //std::cout<<"_tcHits[i].ambigPhi = "<<phi<<std::endl;
      // quality cut for the 1st segment
      if(_lineFitter.qn() <= 2) {
        _lineFitter.addPoint(z, _segmentHits.at(i).at(j).phiDiag, phiWeight);
        continue;
      }
      if(_lineFitter.qn() > 2) {
      if(turns != 0){
          // corss check
          double lineSlope = _lineFitter.dydx();
          double lineIntercept = _lineFitter.y0();
          // Predict phi from the line
          double predictedPhi = lineSlope * z + lineIntercept;
          // Compute the difference between prediction and actual
          double diffPhi[2] = {0.0};
          diffPhi[0] = predictedPhi - _segmentHits.at(i).at(j).phiDiag;
          diffPhi[1] = predictedPhi - _segmentHits.at(i).at(j).phi;
          // choose the nearest assumption
          if(abs(diffPhi[1]) < abs(diffPhi[0])) turns = 0;
          phi = _segmentHits.at(i).at(j).phi + turns * 2 * M_PI;
          _segmentHits.at(i).at(j).phiDiag = phi;
          // Round delta/2π to nearest integer for wrapping correction
          _lineFitter.addPoint(z, _segmentHits.at(i).at(j).phiDiag, phiWeight);
      }
        else {
          _lineFitter.addPoint(z, _segmentHits.at(i).at(j).phiDiag, phiWeight);
        }
      }
  }
      //return fitting values
      double dphidz = _lineFitter.dydx();
      alpha = dphidz;
      beta = _lineFitter.y0();
      chindf = _lineFitter.chi2Dof();
  }
//--------------------------------------------------------------------------------//
    void PhiZSeedFinder::ev5_fit_slope_ver4(int index, double& alpha, double& alphaError, double& beta, double& betaError, double& chindf){
    std::cout<<"ev5_fit_slope_ver4"<<std::endl;
    std::cout<<"_segmentHits.at("<<index<<").size = "<<_segmentHits.at(index).size()<<std::endl;
          for (size_t j = 0; j < _segmentHits.at(index).size(); ++j) {
              const ev5_HitsInNthStation& hit = _segmentHits[index][j];
              std::cout << "   No. " << j
                        << " hitIndice = " << hit.hitIndice
                        << " station = " << hit.station
                        << " plane = "   << hit.plane
                        << " face ="    << hit.face
                        << " panel ="   << hit.panel
                        << " x ="       << hit.x
                        << " y ="       << hit.y
                        << " z ="       << hit.z
                        << " phi ="     << hit.phi
                        << " helixPhi ="     << hit.helixPhi
                        << " hitID ="   << hit.hitID
                        << " segmentIndex =" << hit.segmentIndex
                        << " used ="    << hit.used
                        << " nturn ="    << hit.nturn
                        << std::endl;
          }
    //Get the slope value from the 1st segment
    //continue collecting hits until station gap is > 1
    int Reference_HitIndex = 0;// hit should be located in the upstream tracker
    int Reference_segmentIndex = 0;// hit should be located in the upstream tracker
    //int Reference_HitIndice = -1;// hit should be located in the upstream tracker
    //int Reference_station = 0;
    double Reference_Phi = -999.9;
    for(size_t j=0; j<_segmentHits.at(index).size(); j++){
      Reference_HitIndex = j;
      Reference_segmentIndex = _segmentHits.at(index).at(j).segmentIndex;
      //Reference_HitIndice = _segmentHits.at(i).at(j).hitIndice;
      //Reference_station = _segmentHits.at(i).at(j).station;
      Reference_Phi = _segmentHits.at(index).at(j).helixPhi;
      _segmentHits.at(index).at(j).phiDiag = Reference_Phi;
      /*std::cout<<"best_segments[i][j].Phi = "<<_segmentHits.at(i).at(j).phi<<std::endl;
      std::cout<<"best_segments[i][j].ambigPhi = "<<_segmentHits.at(i).at(j).ambigPhi<<std::endl;
      std::cout<<"station_z = "<<_segmentHits.at(i).at(j).z<<std::endl;
      std::cout<<"station_Reference = "<<_segmentHits.at(i).at(j).station<<std::endl;
      std::cout<<"Reference_Index = "<<Reference_HitIndex<<std::endl;
      std::cout<<"Reference_Indice = "<<Reference_HitIndice<<std::endl;
      std::cout<<"Reference_Phi = "<<Reference_Phi<<std::endl;*/
      break;
    }
    ::LsqSums2 _lineFitter;
    _lineFitter.clear();
      //only add reference hit to the fiiting function
    std::cout<<"Reference segmentIndex/z/Phi = "<<Reference_segmentIndex<<"/"<<_segmentHits.at(index).at(Reference_HitIndex).z<<"/"<<Reference_Phi<<std::endl;
    _lineFitter.addPoint(_segmentHits.at(index).at(Reference_HitIndex).z, Reference_Phi, 0.1);
    // Step 1: find unique segmentIndex values
    std::set<int> uniqueSegmentIndices;
    for (const auto &hit : _segmentHits.at(index)) {
      uniqueSegmentIndices.insert(hit.segmentIndex);
    }
    std::cout << "Number of unique segmentIndex values = "
    << uniqueSegmentIndices.size() << std::endl;
    // Step 2: iterate over each unique segmentIndex and do something
    for (int segIdx : uniqueSegmentIndices) {
      std::cout << "=============================== " << std::endl;
      std::cout << "Processing segmentIndex = " << segIdx << std::endl;
      std::cout << "=============================== " << std::endl;
      for(size_t j=0; j<_segmentHits.at(index).size(); j++){
        std::cout << "No "<<j<<": segIdx/segment" << segIdx << "/"<<_segmentHits.at(index).at(j).segmentIndex<<std::endl;
          double z = _segmentHits.at(index).at(j).z;
          double phiError2 = _segmentHits.at(index).at(j).helixPhiError2;
          double phiWeight = 1.0/phiError2;
        if(segIdx == Reference_segmentIndex){
          double z = _segmentHits.at(index).at(j).z;
          double phiError2 = _segmentHits.at(index).at(j).helixPhiError2;
          double phiWeight = 1.0/phiError2;
          // get other hits
          std::cout<<"kitagawa0 = "<< j<<std::endl;
          if(Reference_HitIndex == (int)j) continue;
          //coorect helixphi and consider 2pi boundary
          float deltaPhi = _segmentHits.at(index).at(j).helixPhi - Reference_Phi;
          /*std::cout<<"hitIndice = "<<_tcHits[i].hitIndice<<std::endl;
          std::cout<<"Helixphi = "<<_tcHits[i].helixPhi<<std::endl;
          std::cout<<"Z = "<<z<<std::endl;
          std::cout<<"station = "<<_tcHits[i].station<<std::endl;
          std::cout<<"phi = "<<_tcHits[i].helixPhi<<std::endl;
          std::cout<<"deltaPhi = "<<deltaPhi<<std::endl;
          */
          // If it turns more than pi then, consinder the 2pi boundary
          int turns = 0;
          if (deltaPhi > M_PI) turns--;
          if (deltaPhi < -M_PI) turns++;
          double phi = _segmentHits.at(index).at(j).helixPhi + turns * 2 * M_PI;
          //std::cout<<"_tcHits[i].ambigPhi = "<<phi<<std::endl;
          // quality cut for the 1st segment
          if(_lineFitter.qn() < 2) {
            phi = phi + _segmentHits.at(index).at(j).nturn * 2 * M_PI;
            _lineFitter.addPoint(z, phi, phiWeight);
            std::cout<<"kitagawa1 = "<< j<<std::endl;
            std::cout<<"Fitter : "<<j<<" hitIndice/z/helixPhi/phi :"<<_segmentHits.at(index).at(j).hitIndice<<"/"<<z<<"/"<<_segmentHits.at(index).at(j).helixPhi<<"/"<<phi<<std::endl;
            continue;
          }
          if(_lineFitter.qn() >= 2) {
          std::cout<<"kitagawa2 = "<< j<<std::endl;
          if(turns != 0){
              phi = phi + _segmentHits.at(index).at(j).nturn * 2 * M_PI;
              double phiOrg = _segmentHits.at(index).at(j).helixPhi + _segmentHits.at(index).at(j).nturn * 2 * M_PI;
              // corss check
              double lineSlope = _lineFitter.dydx();
              double lineIntercept = _lineFitter.y0();
              // Predict phi from the line
              double predictedPhi = lineSlope * z + lineIntercept;
              // Compute the difference between prediction and actual
              double diffPhi[2] = {0.0};
              diffPhi[0] = predictedPhi - phi;
              diffPhi[1] = predictedPhi - phiOrg;
              // choose the nearest assumption
              if(abs(diffPhi[1]) < abs(diffPhi[0])) phi = phiOrg;
              // Round delta/2π to nearest integer for wrapping correction
              std::cout<<"kitagawa3 = "<< j<<std::endl;
              std::cout<<"Fitter : "<<j<<" hitIndice/z/helixPhi/phi :"<<_segmentHits.at(index).at(j).hitIndice<<"/"<<z<<"/"<<_segmentHits.at(index).at(j).helixPhi<<"/"<<phi<<std::endl;
              _lineFitter.addPoint(z, phi, phiWeight);
          }
            else {
double phiOrg = phi;
double candidates[3] = { phiOrg, phiOrg + 2*M_PI, phiOrg - 2*M_PI };
              // corss check
              double lineSlope = _lineFitter.dydx();
              double lineIntercept = _lineFitter.y0();
              // Predict phi from the line
              double predictedPhi = lineSlope * z + lineIntercept;
int best = 0;
double bestDiff = std::fabs(predictedPhi - candidates[0]);
for (int k = 1; k < 3; ++k) {
  double d = std::fabs(predictedPhi - candidates[k]);
  if (d < bestDiff) { bestDiff = d; best = k; }
}
phi = candidates[best];  // choose closest
phi += _segmentHits.at(index).at(j).nturn * 2*M_PI;
_lineFitter.addPoint(z, phi, phiWeight);
              std::cout<<"kitagawa4 = "<< j<<std::endl;
              std::cout<<"Fitter : "<<j<<" hitIndice/z/helixPhi/phi :"<<_segmentHits.at(index).at(j).hitIndice<<"/"<<z<<"/"<<_segmentHits.at(index).at(j).helixPhi<<"/"<<phi<<std::endl;
            }
          }
        } else {
          double phi = _segmentHits.at(index).at(j).helixPhi + _segmentHits.at(index).at(j).nturn * 2 * M_PI;
          _lineFitter.addPoint(z, phi, phiWeight);
          std::cout<<"kitagawa5 = "<< j<<std::endl;
          std::cout<<"Fitter : "<<j<<" hitIndice/z/helixPhi/phi :"<<_segmentHits.at(index).at(j).hitIndice<<"/"<<z<<"/"<<_segmentHits.at(index).at(j).helixPhi<<"/"<<phi<<std::endl;
        }
          std::cout<<"kitagawa6 = "<< j<<std::endl;
      }
    }
      //return fitting values
      double dphidz = _lineFitter.dydx();
      alpha = dphidz;
      alphaError = _lineFitter.dydxErr();
      beta = _lineFitter.y0();
      betaError = _lineFitter.y0Err();
      chindf = _lineFitter.chi2Dof();
  }
//--------------------------------------------------------------------------------//
  void PhiZSeedFinder::ev5_select_best_segments_step_04(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi){
      if(!(ThisIsBestSegment.size() >= 2)) return;
      if(!(ThisIsBestSegment_Diag.size() >= 2)) return;
      std::cout << "-----------------------------------" << std::endl;
      std::cout << "-----------------------------------" << std::endl;
      std::cout << " ev5_select_best_segments_step_04     " << std::endl;
      std::cout << "-----------------------------------" << std::endl;
      std::cout << "-----------------------------------" << std::endl;
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
  void PhiZSeedFinder::ev5_select_best_segments_step_05(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi){
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
      //std::cout << " ev5_select_best_segments_step_06  " << std::endl;
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
  }//end ev5_select_best_segments_step_06
//-----------------------------------------------------------------------------
  PhiZSeedFinder::SegmentComp PhiZSeedFinder::compareSegments(const std::vector<ev5_HitsInNthStation>& seg1, const std::vector<ev5_HitsInNthStation>& seg2) {
    std::cout<<"compareSegments"<<std::endl;
    SegmentComp retval(unique);
    unsigned nh1, nh2, nover;
    // count the StrawHit overlap between the helices
    countHits(seg1, seg2, nh1, nh2, nover);
    std::cout<<"nh1/nh2/nover = "<<nh1<<"/"<<nh2<<"/"<<nover<<std::endl;
    unsigned minh = std::min(nh1, nh2);
    double chih1zphi(0),chih2zphi(0);
    // Calculate the chi-sq of the segment
    findchisq(seg1, chih1zphi);
    findchisq(seg2, chih2zphi);
    std::cout<<"chih1zphi/chih2zphi = "<<chih1zphi<<"/"<<chih2zphi<<"/"<<std::endl;
    std::cout<<"nover/float(minh) = "<<nover/float(minh)<<nover<<std::endl;
    double _minnover = 10;
    double _minoverfrac = 0.5;
    double _deltanh = 5;
    // overlapping segments: decide which is best
    if(nover >= _minnover && nover/float(minh) > _minoverfrac) {
      //if(h1.caloCluster().isNonnull() && h2.caloCluster().isNull())
        //retval = first;
      // Pick the one with a CaloCluster first
      //else if( h2.caloCluster().isNonnull() && h1.caloCluster().isNull())
        //retval = second;
      // then compare active StrawHit counts and if difference of the StrawHit counts greater than deltanh
      if((nh1 > nh2) && (nh1-nh2) > _deltanh)
        retval = first;
      else if((nh2 > nh1) && (nh2-nh1) > _deltanh)
        retval = second;
      // finally compare chisquared: sum xy and fz
      else if(chih1zphi  < chih2zphi)
        retval = first;
      else
        retval = second;
    }
    // if it is still retval = unqiue
    if(retval == unique && nover/float(minh) > _minoverfrac) {
      if(chih1zphi  < chih2zphi)
        retval = first;
      else
        retval = second;
    }
    std::cout<<"Esto = "<<retval<<std::endl;
    return retval;
  }
//-----------------------------------------------------------------------------
void PhiZSeedFinder::countHits(
    const std::vector<ev5_HitsInNthStation>& seg1,
    const std::vector<ev5_HitsInNthStation>& seg2,
    unsigned& nh1,
    unsigned& nh2,
    unsigned& nover){
    nh1 = seg1.size();
    nh2 = seg2.size();
    nover = 0;
    std::unordered_set<int> h1_indices;
    for (const auto& hit : seg1) {
        h1_indices.insert(hit.hitIndice);
    }
    for (const auto& hit : seg2) {
        if (h1_indices.count(hit.hitIndice)) {
            nover++;
        }
    }
  }
//-----------------------------------------------------------------------------
void PhiZSeedFinder::findchisq(std::vector<ev5_HitsInNthStation> const& segment, double& chizphi) const{
    ::LsqSums2 fitter;
    fitter.clear();
    for (size_t j = 0; j < segment.size(); ++j) {
        const auto& hit = segment[j];
        double z = hit.z;
        double phi = hit.phi;
        std::cout << "z/phi = " << z << " / " << phi << std::endl;
        double Error2 = 0.1 * 0.1;
        double phiWeight = 1.0 / Error2; // Replace with actual error if available
        fitter.addPoint(z, phi, phiWeight);
    }
    //return fitting values
    chizphi = fitter.chi2Dof();
  }
//-----------------------------------------------------------------------------
  void PhiZSeedFinder::ev5_select_best_segments_cleanup(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi){
    std::cout << "-----------------------------------" << std::endl;
    std::cout << " ev5_select_best_segments_cleanup  " << std::endl;
    std::cout << "-----------------------------------" << std::endl;
/*    // Make a combined index vector
std::vector<size_t> indices(all_ThisIsBestSegment.size());
std::iota(indices.begin(), indices.end(), 0); // 0, 1, 2, ...
// Define a lambda that returns the max z in a segment
auto getMaxZ = [&](const std::vector<ev5_HitsInNthStation>& segment) {
    double maxZ = -1e9;
    for (const auto& hit : segment) {
        if (hit.z > maxZ) maxZ = hit.z;
    }
    return maxZ;
};
// Sort indices based on descending max Z
std::sort(indices.begin(), indices.end(),
          [&](size_t a, size_t b) {
              return getMaxZ(all_ThisIsBestSegment[a]) >
                     getMaxZ(all_ThisIsBestSegment[b]);
          });
// Apply the new order to both vectors
std::vector<std::vector<ev5_HitsInNthStation>> sortedSegments;
std::vector<std::vector<ev5_Segment>> sortedSegmentsDiag;
for (size_t idx : indices) {
    sortedSegments.push_back(all_ThisIsBestSegment[idx]);
    sortedSegmentsDiag.push_back(all_ThisIsBestSegment_Diag[idx]);
}
// Replace originals
all_ThisIsBestSegment.swap(sortedSegments);
all_ThisIsBestSegment_Diag.swap(sortedSegmentsDiag);
*/
    auto iseg = all_ThisIsBestSegment.begin();
    auto idiag = all_ThisIsBestSegment_Diag.begin();
    while (iseg != all_ThisIsBestSegment.end()) {
        auto jseg = iseg + 1;
        auto jdiag = idiag + 1;
        while (jseg != all_ThisIsBestSegment.end()) {
            auto hcomp = compareSegments(*iseg, *jseg);
            if (hcomp == unique) {
                ++jseg;
                ++jdiag;
            } else if (hcomp == first) {
                jseg = all_ThisIsBestSegment.erase(jseg);
                jdiag = all_ThisIsBestSegment_Diag.erase(jdiag);
            } else if (hcomp == second) {
                iseg = all_ThisIsBestSegment.erase(iseg);
                idiag = all_ThisIsBestSegment_Diag.erase(idiag);
                break;
            }
        }
        if (jseg == all_ThisIsBestSegment.end()) {
            ++iseg;
            ++idiag;
        }
    }
    std::cout << " end  " << std::endl;
    std::cout << "-----------------------------------" << std::endl;
}
// -------------------------------------------------------------------
// Merge segments: each segment as reference in turn, last is standalone
// Remove duplicates when merging
// -------------------------------------------------------------------
void PhiZSeedFinder::mergeSegmentsAll(
    std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment,
    std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag,
    double thre_residual
) {
    std::cout << "mergeSegmentsAll" << std::endl;
    //initialization
   for (size_t segIdx = 0; segIdx < all_ThisIsBestSegment.size(); ++segIdx) {
    for (auto &hit : all_ThisIsBestSegment[segIdx]) {
        hit.nturn = 0;
        hit.segmentIndex = segIdx;
    }
   }
    _segmentHits.clear();
    _segmentHits = all_ThisIsBestSegment;
    int loop = 0;
    int nSegmentsInTC = (int)all_ThisIsBestSegment.size();
    std::cout << "nSegmentsInTC = " << nSegmentsInTC << std::endl;
    if (nSegmentsInTC == 0) return;
    bool didMerge = true;
    while (didMerge) {
      std::cout<<"didMerge loop = "<<loop<<std::endl;
      loop++;
      didMerge = false;
      _segmentHits.clear();
      _segmentHits = all_ThisIsBestSegment;
      std::cout << "\n--- Dumping all_ThisIsBestSegment ---" << std::endl;
      std::cout << "\n all_ThisIsBestSegment size = " << all_ThisIsBestSegment.size() << std::endl;
      for (size_t i = 0; i < all_ThisIsBestSegment.size(); ++i) {
          std::cout << " Segment " << i << " (#hits = "
                    << all_ThisIsBestSegment[i].size() << ")" << std::endl;
          for (size_t j = 0; j < all_ThisIsBestSegment[i].size(); ++j) {
              const ev5_HitsInNthStation& hit = all_ThisIsBestSegment[i][j];
              std::cout << "   No. " << j
                        << " hitIndice = " << hit.hitIndice
                        << " station = " << hit.station
                        << " plane = "   << hit.plane
                        << " face ="    << hit.face
                        << " panel ="   << hit.panel
                        << " x ="       << hit.x
                        << " y ="       << hit.y
                        << " z ="       << hit.z
                        << " phi ="     << hit.phi
                        << " hitID ="   << hit.hitID
                        << " segmentIndex =" << hit.segmentIndex
                        << " used ="    << hit.used
                        << " nturn ="    << hit.nturn
                        << std::endl;
          }
      }
      std::cout << "\n--- Dumping all_ThisIsBestSegment_Diag ---" << std::endl;
      for (size_t i = 0; i < all_ThisIsBestSegment_Diag.size(); ++i) {
          std::cout << " SegmentDiag " << i << " (#entries = "
                    << all_ThisIsBestSegment_Diag[i].size() << ")" << std::endl;
          for (size_t j = 0; j < all_ThisIsBestSegment_Diag[i].size(); ++j) {
              const ev5_Segment& seg = all_ThisIsBestSegment_Diag[i][j];
              std::cout << "   No. " << j
                        << " station="   << seg.station
                        << " z="         << seg.z
                        << " deltaphi="  << seg.deltaphi
                        << " alpha="     << seg.alpha
                        << " beta="      << seg.beta
                        << " chiNDF="    << seg.chiNDF
                        << " ref_point=" << seg.reference_point
                        << " usedForFit="<< seg.usedforfit
                        << std::endl;
          }
      }
          for (size_t refIdx = 0; refIdx < all_ThisIsBestSegment.size(); ++refIdx) {
              std::cout << "\n[mergeSegmentsAll] Reference segment = " << refIdx
                        << " (#hits = " << all_ThisIsBestSegment[refIdx].size() << ")" << std::endl;
              for (size_t testIdx = refIdx + 1; testIdx < all_ThisIsBestSegment.size(); ++testIdx) {
                // === Decide whether to merge refIdx and testIdx ===
                bool canMerge[3] = {false, false, false};
                std::vector<ev5_HitsInNthStation> segmentHits;
              std::cout << "\n[mergeSegmentsAll] Test Reference segment = " << testIdx
                        << " (#hits = " << all_ThisIsBestSegment[testIdx].size() << ")" << std::endl;
              //---------------------------------------
              // circle fit: Step1 w/o correct weight
              //---------------------------------------
              _circleFitter.clear();
              for (size_t i = 0; i < all_ThisIsBestSegment[refIdx].size(); i++) {
                  double x  = all_ThisIsBestSegment[refIdx][i].x;
                  double y  = all_ThisIsBestSegment[refIdx][i].y;
                  double wP = 0.1;  // tentative value
                  _circleFitter.addPoint(x, y, wP);
              }
              for (size_t i = 0; i < all_ThisIsBestSegment[testIdx].size(); i++) {
                  double x  = all_ThisIsBestSegment[testIdx][i].x;
                  double y  = all_ThisIsBestSegment[testIdx][i].y;
                  double wP = 0.1;  // tentative value
                  _circleFitter.addPoint(x, y, wP);
              }
              double xC = _circleFitter.x0();
              double yC = _circleFitter.y0();
              double rC = _circleFitter.radius();
              //-------------------------------------------------------------------
              // circle fit: Step2 w/ correct weight
              //-------------------------------------------------------------------
              _circleFitter.clear();
              for (int idx : {refIdx, testIdx}) {
                for (size_t i = 0; i < all_ThisIsBestSegment[idx].size(); i++) {
                  int hitIndice = all_ThisIsBestSegment[idx][i].hitIndice;
                  //all_ThisIsBestSegment[idx][i].segmentIndex = idx;
                  int nStrawHits = all_ThisIsBestSegment[idx][i].strawhits;
                  double circleError2 = computeCircleError2_ver2(hitIndice, nStrawHits, xC, yC, rC);
                  all_ThisIsBestSegment[idx][i].circleError2 = circleError2;
                  double x = all_ThisIsBestSegment[idx][i].x;
                  double y = all_ThisIsBestSegment[idx][i].y;
                  double wP = 1.0 / (circleError2);
                  _circleFitter.addPoint(x, y, wP);
                }
                segmentHits.insert(segmentHits.end(), all_ThisIsBestSegment[idx].begin(), all_ThisIsBestSegment[idx].end());
              }
              for (auto &hit : segmentHits) hit.used = true;
              xC = _circleFitter.x0();
              yC = _circleFitter.y0();
              rC = _circleFitter.radius();
              //---------------------------------------
              // circle fit: Step3 ciclefit clean up
              //---------------------------------------
              std::cout<<"Before clean up "<<std::endl;
              std::cout<<"# of hits in Helix = "<<_circleFitter.qn()<<std::endl;
              std::cout<<"xC/yC = "<<_circleFitter.x0()<<"/"<<_circleFitter.y0()<<std::endl;
              std::cout<<"radius = "<<_circleFitter.radius()<<std::endl;
              std::cout<<"phi/dfdz/chi2DofC/chi2DofLineC = "<<_circleFitter.phi0()<<"/"<<_circleFitter.dfdz()<<"/"<<_circleFitter.chi2DofCircle()<<"/"<<_circleFitter.chi2DofLine()<<std::endl;
              if(_circleFitter.qn() > 10 and _circleFitter.chi2DofCircle() > 5.0){
              std::vector<cleanup> remove_hits;
              double chi2ndf = _circleFitter.chi2DofCircle();
              while(chi2ndf > 5.0){
                int remove_hitIndex;
                int remove_hitIndice;
                int find = 0;
                for(size_t i=0; i<segmentHits.size(); i++){
                  if(segmentHits[i].used == false) continue;
                  _circleFitter.clear();
                  for(size_t j=0; j<segmentHits.size(); j++){
                    if(segmentHits[j].used == false) continue;
                    if(i==j) continue;
                    int hitIndice = segmentHits[j].hitIndice;
                    int nStrawHits = segmentHits[j].strawhits;
                    double circleError2 = computeCircleError2_ver2(hitIndice, nStrawHits, xC, yC, rC);
                    segmentHits[j].circleError2 = circleError2;
                    double x = segmentHits[j].x;
                    double y = segmentHits[j].y;
                    double wP = 1.0 / (circleError2);
                    _circleFitter.addPoint(x, y, wP);
                  }
                  if(_circleFitter.chi2DofCircle() < chi2ndf) {
                    remove_hitIndex = i;
                    remove_hitIndice = segmentHits[i].hitIndice;
                    chi2ndf = _circleFitter.chi2DofCircle();
                    find++;
                  std::cout<<"Indice :"<<i<<"/"<<"chi2ndf: "<<_circleFitter.chi2DofCircle()<<std::endl;
                  std::cout<<"find :"<<find<<std::endl;
                  }
                }
              std::cout<<"find :"<<find<<std::endl;
              //check if chi2ndf is improved or not
              if(find > 0){
                cleanup circlefit;
                circlefit.tcindex = remove_hitIndex;
                circlefit.tcindice = remove_hitIndice;
                circlefit.chi2ndf = chi2ndf;
                remove_hits.push_back(circlefit);
                segmentHits[remove_hitIndex].used = false;
              //recalculate the circle parameter
              //Level 1:
              _circleFitter.clear();
              for(size_t i=0; i<segmentHits.size(); i++){
              if(segmentHits[i].used == false) continue;
                double x = segmentHits[i].x;
                double y = segmentHits[i].y;
                double wP = 1.0 / (segmentHits[i].circleError2);
                _circleFitter.addPoint(x, y, wP);
              }
              chi2ndf = _circleFitter.chi2DofCircle();
              xC = _circleFitter.x0();
              yC = _circleFitter.y0();
              rC = _circleFitter.radius();
              //std::cout<<"remove_hitIndex :"<< remove_hitIndex<<"/"<<"chi2ndf: "<<chi2ndf<<std::endl;
              if(chi2ndf < 5.0) break;
              if((segmentHits.size() - remove_hits.size()) <= 10) break;
              }else{
              break;
              }
              }
              //std::sort(remove_hits.begin(), remove_hits.end(), [](const cleanup& a, const cleanup& b) { return a.chi2ndf < b.chi2ndf; } );
              std::cout<<"==================================="<<std::endl;
              std::cout<<"clean-up"<<std::endl;
              std::cout<<"==================================="<<std::endl;
              for(size_t i=0; i<remove_hits.size(); i++){
              int tcindex = remove_hits.at(i).tcindex;
              int tcdice = remove_hits.at(i).tcindice;
              double chi2ndf_ = remove_hits.at(i).chi2ndf;
              std::cout<<"No: "<<tcindex<< ", Indice = "<< tcdice <<", chi2DofCircle (w/o this hit) = "<<chi2ndf_<<std::endl;
              }
              //Step4: iterate over selected combohits and recalculate the weight value and refit again
              //recalculate the weight
              std::cout<<"==================================="<<std::endl;
              std::cout<<"After clean-up"<<std::endl;
              std::cout<<"==================================="<<std::endl;
              _circleFitter.clear();
              int count = 0;
              for(size_t i=0; i<segmentHits.size(); i++){
              if(segmentHits[i].used == false) continue;
              double x = segmentHits[i].x;
              double y = segmentHits[i].y;
              double wP = 1.0 / (segmentHits[i].circleError2);
              _circleFitter.addPoint(x, y, wP);
              count++;
              }
              std::cout<<"Hits used in circle: "<<count<< ", chi2DofCircle = "<<_circleFitter.chi2DofCircle()<<std::endl;
              // xC = _circleFitter.x0();
              // yC = _circleFitter.y0();
              // _circleFitter.clear();
              // for(size_t i=0; i<nComboHitsInSegment; i++){
              //   if(_tcHits[i].used == false) continue;
              //   computeCircleError2(i, xC, yC);
              //   double x = _tcHits.at(i).x;
              //   double y = _tcHits.at(i).y;
              //   double wP = 1.0 / (_tcHits[i].circleError2);
              //   _circleFitter.addPoint(x, y, wP);
              //   _tcHits[i].used = true;
              // }
              }
              // === STEP 0: Backup original data ===
              auto backup_all_ThisIsBestSegment      = all_ThisIsBestSegment;
              auto backup_all_ThisIsBestSegment_Diag = all_ThisIsBestSegment_Diag;
              // === STEP 1: Remove hits ===
              for (int idx : {refIdx, testIdx}) {
                  auto& segment     = all_ThisIsBestSegment[idx];
                  auto& segmentDiag = all_ThisIsBestSegment_Diag[idx];
                  if (segment.size() != segmentDiag.size()) {
                      std::cerr << "Warning: segment size mismatch for idx "
                                << idx << " (" << segment.size() << " vs "
                                << segmentDiag.size() << ")\n";
                  }
                  auto itSeg  = segment.begin();
                  auto itDiag = segmentDiag.begin();
                  while (itSeg != segment.end() && itDiag != segmentDiag.end()) {
                      bool toRemove = false;
                      for (const auto& tcHit : segmentHits) {
                          if (!tcHit.used && tcHit.hitIndice == itSeg->hitIndice) {
                              std::cout << "Removing hitIndice = " << itSeg->hitIndice
                                        << " from segment " << idx << std::endl;
                              toRemove = true;
                              break;
                          }
                      }
                      if (toRemove) {
                          itSeg  = segment.erase(itSeg);
                          itDiag = segmentDiag.erase(itDiag);
                      } else {
                          ++itSeg;
                          ++itDiag;
                      }
                  }
              }
              _segmentHits.clear();
              _segmentHits = all_ThisIsBestSegment;
              //for only plot
              tcHitsFill(refIdx);
              tcHitsFill_Add(testIdx);
              //plot_XVsY(refIdx, testIdx, "mergeSegmentsAll_step1", xC, yC, rC);
              //---------------------------------------
              // refIdx hits: compute helix phi
              //---------------------------------------
              for (size_t i = 0; i < segmentHits.size(); i++) {
                  if(segmentHits[i].used == false) continue;
                  int hitIndice = segmentHits[i].hitIndice;
                  double helixPhi = 0.0;
                  double helixPhiError2 = 0.0;
                  computeHelixPhi_ver2(hitIndice, xC, yC, helixPhi, helixPhiError2);
                  segmentHits[i].helixPhi = helixPhi;
                  segmentHits[i].helixPhiError2 = helixPhiError2;
                    std::cout << "refIdx  = " << refIdx
                    << " | index = " << i
                    << " | hitIndice = " << segmentHits[i].hitIndice
                    << " | wireErr = " << _data.chcol->at(hitIndice).wireRes()
                    << " | z = " << segmentHits[i].x
                    << " | z = " << segmentHits[i].y
                    << " | z = " << segmentHits[i].z
                    << " | phi = " << segmentHits[i].phi
                    << " | helixPhi = " << segmentHits[i].helixPhi
                    << " | helixPhiError2= " << segmentHits[i].helixPhiError2
                    << " | sqrt(helixPhiError2) = " << sqrt(segmentHits[i].helixPhiError2)
                    << std::endl;
                   for (size_t j = 0; j < all_ThisIsBestSegment[refIdx].size(); j++) {
                      if(hitIndice != all_ThisIsBestSegment[refIdx][j].hitIndice) continue;
                      all_ThisIsBestSegment[refIdx][j].helixPhi = helixPhi;
                      all_ThisIsBestSegment[refIdx][j].helixPhiError2 = helixPhiError2;
                    }
                    for (size_t j = 0; j < all_ThisIsBestSegment[testIdx].size(); j++) {
                      if(hitIndice != all_ThisIsBestSegment[testIdx][j].hitIndice) continue;
                      all_ThisIsBestSegment[testIdx][j].helixPhi = helixPhi;
                      all_ThisIsBestSegment[testIdx][j].helixPhiError2 = helixPhiError2;
                    }
              }
              _segmentHits.clear();
              _segmentHits = all_ThisIsBestSegment;
              // slope value for reference
              double slope_alpha = 0.0;
              double slope_alphaError = 0.0;
              double slope_beta  = 0.0;
              double slope_betaError  = 0.0;
              double ChiNDF      = 0.0;
              std::cout<<"didMerge loop = "<<loop<<std::endl;
              std::cout<<"refIdx = "<<refIdx<<std::endl;
              ev5_fit_slope_ver4(refIdx, slope_alpha, slope_alphaError, slope_beta, slope_betaError, ChiNDF);
              //for only plot
              tcHitsFill(refIdx);
              std::cout<<"refIdx = "<<refIdx<<std::endl;
              //plot_PhiVsZ_forSegment_ver2(refIdx, testIdx, slope_alpha, slope_beta, ChiNDF);
              std::cout << std::fixed << std::setprecision(10);
              std::cout<<"refIdx = "<<refIdx<<std::endl;
              std::cout << "slope_alpha :"<< slope_alpha <<" +/- " << slope_alphaError <<std::endl;
              std::cout << "slope_beta :"<< slope_beta <<" +/- " << slope_betaError <<std::endl;
              std::cout << "ChiNDF : " << ChiNDF << std::endl;
              // slope value for test
              double test_slope_alpha = 0.0;
              double test_slope_alphaError = 0.0;
              double test_slope_beta  = 0.0;
              double test_slope_betaError  = 0.0;
              double test_ChiNDF      = 0.0;
              std::cout<<"didMerge loop = "<<loop<<std::endl;
              std::cout<<"testIdx = "<<testIdx<<std::endl;
              ev5_fit_slope_ver4(testIdx, test_slope_alpha, test_slope_alphaError, test_slope_beta, test_slope_betaError, test_ChiNDF);
              //for only plot
              tcHitsFill(testIdx);
              //plot_PhiVsZ_forSegment_ver2(testIdx, refIdx, test_slope_alpha, test_slope_beta, test_ChiNDF);
              std::cout<<"testIdx = "<<testIdx<<std::endl;
              std::cout << "slope_alpha :"<< test_slope_alpha <<" +/- " << test_slope_alphaError <<std::endl;
              std::cout << "slope_beta :"<< test_slope_beta <<" +/- " << test_slope_betaError <<std::endl;
              std::cout << "ChiNDF : " << test_ChiNDF << std::endl;
              //---------------------------------------
              // slope difference method
              //---------------------------------------
              /*double alpha_diff = 0.0;
              double alpha_sigma = 0.0003199;
              if (slope_alpha * test_slope_alpha > 0)
                alpha_diff = fabs(slope_alpha - test_slope_alpha);
              else{
                alpha_diff = fabs(slope_alpha) + fabs(test_slope_alpha);
              }
              std::cout << "alpha_diff = "<< alpha_diff <<std::endl;
              std::cout << "alpha_sigma*5 = "<< alpha_sigma * 5 <<std::endl;
              if (alpha_diff > alpha_sigma * 5) continue;
              */
              double alpha_sigma = 0.0003199;
              double alpha_diff = (slope_alpha * test_slope_alpha > 0)
                      ? std::fabs(slope_alpha - test_slope_alpha)
                      : std::fabs(slope_alpha) + std::fabs(test_slope_alpha);
              double combined_alpha_error = std::sqrt(
                  slope_alphaError * slope_alphaError +
                  test_slope_alphaError * test_slope_alphaError
              );
              // prevent division-by-zero / unrealistically small combined error
              const double min_combined_error = 1e-12;
              if (combined_alpha_error < min_combined_error) combined_alpha_error = min_combined_error;
              double empirical_threshold = alpha_sigma * 5.0;
              double statistical_threshold = 5.0 * combined_alpha_error;
              // pick the more conservative (larger) threshold
              double threshold = std::max(empirical_threshold, statistical_threshold);
              if (alpha_diff > threshold) {
                  // Not compatible -> skip/continue
                  //continue;
                  canMerge[0] = false;
              }
              else canMerge[0] = true;
              // ok to merge
              //Pure statistical test (recommended if fit errors are trusted)
              /*double combined_alpha_error = std::sqrt(
              slope_alphaError * slope_alphaError +
                  test_slope_alphaError * test_slope_alphaError
              );
              const double min_err = 1e-12;
              if (combined_alpha_error < min_err) combined_alpha_error = min_err;
              double Z = alpha_diff / combined_alpha_error; // # of sigma
              if (Z > 5.0) continue; // require compatibility within 5 sigma
              // optionally still check empirical threshold:
              // if (alpha_diff > alpha_sigma*5.0) continue;
              */
              std::cout << std::fixed << std::setprecision(6)
              << "alpha_diff=" << alpha_diff
              << " combined_err=" << combined_alpha_error
              << " stat_thresh=" << statistical_threshold
              << " empirical_thresh=" << empirical_threshold
              << " chosen_thresh=" << threshold
              << " Z=" << alpha_diff / combined_alpha_error
              << std::endl;
              //---------------------------------------
              // station overlap check
              //---------------------------------------
              bool station_overlap = false;
              for (auto& hitRef : all_ThisIsBestSegment[refIdx]) {
                for (auto& hitTest : all_ThisIsBestSegment[testIdx]) {
                       if (hitRef.station == hitTest.station) {
                           station_overlap = true;
                           break;
                       }
                }
                if (station_overlap) break;
              }
              if (station_overlap) canMerge[1] = false;
              else canMerge[1] = true;
              //---------------------------------------
              // slope intercept method: 1
              //---------------------------------------
                  // Compute difference after 2pi correction
             /*     double deltaPhi = slope_beta - test_slope_beta;
                  int deltaCorrection = std::round(deltaPhi / (2 * M_PI));
                  double beta_aligned = test_slope_beta + deltaCorrection * 2 * M_PI;
                  // Absolute beta difference
                  double beta_diff = std::fabs(slope_beta - beta_aligned);
                  // Combine uncertainties
                  double combined_beta_error = std::sqrt(
                      slope_betaError * slope_betaError +
                      test_slope_betaError * test_slope_betaError
                  );
                  // Avoid divide-by-zero
                  if (combined_beta_error < 1e-12)
                      combined_beta_error = 1e-12;
                  // Define thresholds
                  double empirical_threshold_beta = 0.4; //tentative value
                  double statistical_threshold_beta = 5.0 * combined_beta_error;
                  // Pick the more conservative (larger) threshold
                  double beta_threshold = std::max(empirical_threshold_beta, statistical_threshold_beta);
                  // Decision
                  canMerge = (beta_diff <= beta_threshold);
                  // Optional: print diagnostic info
                  std::cout << std::fixed << std::setprecision(6)
                            << "refIdx=" << refIdx
                            << " testIdx=" << testIdx
                            << " | beta_diff=" << beta_diff
                            << " | combined_err=" << combined_beta_error
                            << " | threshold=" << beta_threshold
                            << " | (emp=" << empirical_threshold_beta
                            << ", stat=" << statistical_threshold_beta << ")"
                            << " Z=" << beta_diff / combined_beta_error
                            << std::endl;
                  */
              //---------------------------------------
              // slope intercept method: 2
              //---------------------------------------
              //---------------------------------------
              // Predicted phi consistency check
              //---------------------------------------
              // Compute difference after 2pi correction
              double deltaPhi = slope_beta - test_slope_beta;
              int deltaCorrection = std::round(deltaPhi / (2 * M_PI));
              test_slope_beta = test_slope_beta + deltaCorrection * 2 * M_PI;
              double z_test = 0.0;
              for (size_t i = 0; i < all_ThisIsBestSegment[testIdx].size(); i++) {
                z_test = z_test + all_ThisIsBestSegment[testIdx][i].z;
              }
              z_test = z_test/(int)all_ThisIsBestSegment[testIdx].size();
              std::cout<<"z_test = "<< z_test <<std::endl;
              // Predicted phi from ref segment at z of test segment
              double phi_pred = slope_alpha * z_test + slope_beta;
              double phi_obs  = test_slope_alpha * z_test + test_slope_beta;
              // Normalize to [-pi, pi] to handle wrap-around
              double delta_phi = phi_pred - phi_obs;
              if (delta_phi > M_PI)  delta_phi -= 2 * M_PI;
              if (delta_phi < -M_PI) delta_phi += 2 * M_PI;
              // Uncertainty propagation
              double phi_pred_err = std::sqrt(
                  (z_test * z_test * slope_alphaError * slope_alphaError) +
                  (slope_betaError * slope_betaError)
              );
              double phi_obs_err = std::sqrt(
                  (z_test * z_test * test_slope_alphaError * test_slope_alphaError) +
                  (test_slope_betaError * test_slope_betaError)
              );
              double combined_phi_err = std::sqrt(phi_pred_err * phi_pred_err +
                                                 phi_obs_err * phi_obs_err);
              double Z = std::fabs(delta_phi) / (combined_phi_err + 1e-12);
              // Define thresholds
              double empirical_threshold_phi = 0.4;  // radians (example)
              double statistical_threshold_phi = 5.0 * combined_phi_err;
              double phi_threshold = std::max(empirical_threshold_phi, statistical_threshold_phi);
              // Decision
              canMerge[2] = (fabs(delta_phi) <= phi_threshold);
              // Diagnostic printout
              // Diagnostic printout
std::cout << "phi_pred = " << phi_pred << std::endl;
std::cout << "phi_obs  = " << phi_obs << std::endl;
std::cout << "delta_phi = " << delta_phi << std::endl;
std::cout << "phi_pred_err = " << phi_pred_err << std::endl;
std::cout << "phi_obs_err  = " << phi_obs_err << std::endl;
std::cout << "combined_phi_err = " << combined_phi_err << std::endl;
std::cout << "empirical_threshold_phi = " << empirical_threshold_phi << std::endl;
std::cout << "statistical_threshold_phi = " << statistical_threshold_phi << std::endl;
std::cout << "phi_threshold = " << phi_threshold << std::endl;
std::cout << "z(delta_phi/combined_err) = " << Z << std::endl;
std::cout << "canMerge[0] = " << (canMerge[0] ? "true" : "false") << std::endl;
std::cout << "canMerge[1] = " << (canMerge[1] ? "true" : "false") << std::endl;
std::cout << "canMerge[2] = " << (canMerge[2] ? "true" : "false") << std::endl;
                  //shift 2 pi for
                  if (1 == canMerge[0] * canMerge[1] * canMerge[2]) {
                    for (size_t i = 0; i < all_ThisIsBestSegment[testIdx].size(); i++) {
                      all_ThisIsBestSegment[testIdx][i].nturn = deltaCorrection;
                    }
                  }
                  // ===== merge testIdx into refIdx =====
                  if (1 == canMerge[0] * canMerge[1] * canMerge[2]) {
                    std::cout << "   --> Merged testIdx " << testIdx
                            << " into refIdx " << refIdx << std::endl;
                // merge directly into refIdx
                all_ThisIsBestSegment[refIdx].insert(
                    all_ThisIsBestSegment[refIdx].end(),
                    all_ThisIsBestSegment[testIdx].begin(),
                    all_ThisIsBestSegment[testIdx].end()
                );
                all_ThisIsBestSegment_Diag[refIdx].insert(
                    all_ThisIsBestSegment_Diag[refIdx].end(),
                    all_ThisIsBestSegment_Diag[testIdx].begin(),
                    all_ThisIsBestSegment_Diag[testIdx].end()
                );
                // remove merged segment
                all_ThisIsBestSegment.erase(all_ThisIsBestSegment.begin() + testIdx);
                all_ThisIsBestSegment_Diag.erase(all_ThisIsBestSegment_Diag.begin() + testIdx);
                std::cout << "   --> merged " << testIdx << " into " << refIdx << std::endl;
                didMerge = true;
                // restart because the vectors changed
                goto restartLoop;
              }
                // === STEP 3: If processing failed, restore backup ===
                if (0 == canMerge[0] * canMerge[1] * canMerge[2]) {
                  std::cout << "Selection failed. Restoring original hit segments...\n";
                  all_ThisIsBestSegment      = backup_all_ThisIsBestSegment;
                  all_ThisIsBestSegment_Diag = backup_all_ThisIsBestSegment_Diag;
                }
            }
          }
        // no merge found in this pass
        break;
        restartLoop:
          std::cout << "[mergeSegmentsAll] Done. Total candidates = "
              << all_ThisIsBestSegment.size() << std::endl;
            // For ev5_HitsInNthStation segments
        for (auto &seg : all_ThisIsBestSegment) {
            std::sort(seg.begin(), seg.end(),
                      [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b) {
                          return a.z < b.z; // ascending
                      });
        }
        // For ev5_Segment segments
        for (auto &seg : all_ThisIsBestSegment_Diag) {
            std::sort(seg.begin(), seg.end(),
                      [](const ev5_Segment& a, const ev5_Segment& b) {
                          return a.z < b.z; // ascending
                      });
        }
        //break;
        if(loop == 2) break;
        continue;
    }
    std::cout << "Final #segments = " << all_ThisIsBestSegment.size() << std::endl;
}
//-----------------------------------------------------------------------------
  void PhiZSeedFinder::ev5_select_best_segments_step_07(std::vector<std::vector<ev5_HitsInNthStation>>& all_BestSegmentInfo, std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi, int& NumberOfSegments){
      std::cout << "-----------------------------------" << std::endl;
      std::cout << "-----------------------------------" << std::endl;
      std::cout << " ev5_select_best_segments_step_07  " << std::endl;
      std::cout << "-----------------------------------" << std::endl;
      std::cout << "-----------------------------------" << std::endl;
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
    int nSegmentsInTC = (int)_segmentHits.size();
    std::vector<double> alpha;
    std::vector<std::vector<int>> nstations;
    std::vector<std::vector<double>> phi_window;
    for(int i=0; i<nSegmentsInTC; i++){
      double slope_alpha = 0.0;//Slope (a)
      double slope_beta = 0.0;//Intercept (b)
      double ChiNDF = 0.0;
      ev5_fit_slope_ver2(i, slope_alpha, slope_beta, ChiNDF);
      alpha.push_back(slope_alpha);
      std::vector<int> station;
      station.clear();
      std::vector<double> phi;
      phi.clear();
      double hit_phi[2] = {9999.9, -9999.9};//[0] = min, [1] = max
      for(int j=0; j<(int)_segmentHits.at(i).size(); j++){
        if(hit_phi[0] > _segmentHits.at(i).at(j).phi) hit_phi[0] = _segmentHits.at(i).at(j).phi;
        if(hit_phi[1] < _segmentHits.at(i).at(j).phi) hit_phi[1] = _segmentHits.at(i).at(j).phi;
        std::cout<<"phi = "<<_segmentHits.at(i).at(j).phi<<std::endl;
        station.push_back(_segmentHits.at(i).at(j).station);
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
        int phi_flag[4] = {0, 0, 0, 0};
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
        std::cout<<"Pass PhiWindow"<<flag<<std::endl;
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
    //kitagawa
    std::cout<<"all_ThisIsBestSegment.size() = "<<all_ThisIsBestSegment.size()<<std::endl;
    for(int j=0; j<(int)all_ThisIsBestSegment.size(); j++){
      for(int k=0; k<(int)all_ThisIsBestSegment.at(j).size(); k++){
        std::cout<<"all_ThisIsBestSegment.at(j).segmentIndex = "<<all_ThisIsBestSegment.at(j).at(k).segmentIndex<<std::endl;
      }
    }
    std::cout<<"static_cast<int>(ncandidate.size()) = "<< (int)ncandidate.size()<<std::endl;
    for(int i = 0; i < static_cast<int>(ncandidate.size()); i++){
      std::cout<<"i = "<< i<<std::endl;
      for(int j = 0; j < static_cast<int>(ncandidate.at(i).size()); j++){
      int index = ncandidate.at(i).at(j);
      std::cout<<"segment_index = "<<index<<std::endl;
      for(int k=0; k<(int)all_ThisIsBestSegment.size(); k++){
        if(index != k) continue;
        for(int l=0; l<(int)all_ThisIsBestSegment.at(k).size(); l++){
          all_ThisIsBestSegment.at(k).at(l).segmentIndex = i;
          //std::cout<<"all_ThisIsBestSegment.at(j).segmentIndex = "<<all_ThisIsBestSegment.at(k).at(l).segmentIndex<<std::endl;
        }
      }
      }
    }
    std::cout<<"after_fill"<<std::endl;
    for(int j=0; j<(int)all_ThisIsBestSegment.size(); j++){
      std::cout<<" j = "<<j<<std::endl;
      for(int k=0; k<(int)all_ThisIsBestSegment.at(j).size(); k++){
        std::cout<<"all_ThisIsBestSegment.at(j).segmentIndex = "<<all_ThisIsBestSegment.at(j).at(k).segmentIndex<<std::endl;
      }
    }
    all_BestSegmentInfo = all_ThisIsBestSegment;
/////////////////////kitagawa
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
    _segmentHits = hitsIn_best_segments;
    NumberOfSegments = (int)hitsIn_best_segments.size();
  }//end ev5_select_best_segments_step_07
//-----------------------------------------------------------------------------
  void PhiZSeedFinder::ev5_select_best_segments_step_08(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi, int& NumberOfSegments){
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << " ev5_select_best_segments_step_08  " << std::endl;
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
        int overlappedHits = 0;
        for (int k = 0; k <(int)hitID_list.at(j).size(); k++) {
          for (int l = 0; l <(int)hitID_list.at(i).size(); l++) {
            if(hitID_list.at(i).at(l) == hitID_list.at(j).at(k)) overlappedHits++;
          }
        }
        double fraction[2] = {0.0};
        std::cout<<"hitsA/hitsB/OverlapHits = "<<hitID_list.at(i).size()<<"/"<<hitID_list.at(j).size()<<"/"<<overlappedHits<<std::endl;
        fraction[0] = (double)overlappedHits/hitID_list.at(i).size();
        fraction[1] = (double)overlappedHits/hitID_list.at(j).size();
        std::cout<<"fraction[0]/fraction[1] = "<<fraction[0]<<"/"<<fraction[1]<<std::endl;
        if(fraction[0] > 0.7 and fraction[1] > 0.7){
          if(fraction[0] > fraction[1]) delete_index.push_back(i);
          else delete_index.push_back(j);
        }
        if(fraction[0] > 0.8 and fraction[1] < 0.3){
          delete_index.push_back(i);
        }
      }
    }
    // make a new hitID after removing the overlap event
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
    all_ThisIsBestSegment.clear();
    all_ThisIsBestSegment_Diag.clear();
    all_ThisIsBestSegment = temp_ThisIsBestSegment;
    all_ThisIsBestSegment_Diag = temp_ThisIsBestSegment_Diag;
    //std::cout<<"ThisIsBestSegment size = "<<all_ThisIsBestSegment.size()<<std::endl;
  }//end ev5_select_best_segments_step_08
//-----------------------------------------------------------------------------
// finding circle from triplet
//-----------------------------------------------------------------------------
void PhiZSeedFinder::initTriplet(triplet& trip, int& outcome) {
  _circleFitter.addPoint(trip.i.pos->x(), trip.i.pos->y());
  _circleFitter.addPoint(trip.j.pos->x(), trip.j.pos->y());
  _circleFitter.addPoint(trip.k.pos->x(), trip.k.pos->y());
  // check if circle is valid for search or if we should continue to next triplet
  //double radius = _circleFitter.radius();
  //double pt = computeHelixPerpMomentum(radius, _bz0);
  //if (pt < _minHelixPerpMomentum || pt > _maxHelixPerpMomentum) {
    //outcome = 0;
  //} else {
    //outcome = 1;
  //}
  outcome = 1;
}
//-----------------------------------------------------------------------------
// fill vector with hits of a segment to search for helix
//-----------------------------------------------------------------------------
void PhiZSeedFinder::tcHitsFill(int isegment) {
  _tcHits.clear();
  cHit hit;
  hit.circleError2 = 1.0;
  hit.helixPhi = 0.0;
  hit.helixPhiError2 = 0.0;
  hit.helixPhiCorrection = 0.0;
  hit.ambigPhi = 0.0;
  hit.segmentIndice = 0;
  hit.inHelix = false;
  hit.used = false;//default is false
  hit.isolated = false;
  hit.averagedOut = false;
  hit.notOnLine = true;
  hit.uselessTripletSeed = false;
  hit.notOnSegment = true;
  hit.station = 0;
  hit.plane = 0;
  hit.face = 0;
  hit.panel = 0;
  hit.x = 0.0;
  hit.y = 0.0;
  hit.z = 0.0;
  hit.phi = 0.0;
  hit.strawhits = 0;
  // fill hits from time cluster of i-th segment
  //std::cout<<"====== tcHitsFill ======"<<std::endl;
  //std::cout<< "run: " << run << " subRun: " << subrun << " event: " << eventNumber<<std::endl;
  //std::cout<<"isegment = "<<isegment<<std::endl;
  //std::cout<<"isegment.size = "<<_segmentHits.at(isegment).size()<<std::endl;
  for (size_t i = 0; i < _segmentHits.at(isegment).size(); i++) {
   //std::cout<<"i-th = "<<i<<std::endl;
  // std::cout<<"hitIndice = "<<_segmentHits.at(isegment).at(i).hitIndice<<std::endl;
   hit.hitIndice = _segmentHits.at(isegment).at(i).hitIndice;
   hit.station = _segmentHits.at(isegment).at(i).station;
   hit.plane = _segmentHits.at(isegment).at(i).plane;
   hit.face = _segmentHits.at(isegment).at(i).face;
   hit.panel = _segmentHits.at(isegment).at(i).panel;
   hit.x = _segmentHits.at(isegment).at(i).x;
   hit.y = _segmentHits.at(isegment).at(i).y;
   hit.z = _segmentHits.at(isegment).at(i).z;
   hit.phi = _segmentHits.at(isegment).at(i).phi;
   hit.segmentIndice = i;
   hit.helixPhi = _segmentHits.at(isegment).at(i).helixPhi;
   hit.helixPhiError2 = _segmentHits.at(isegment).at(i).helixPhiError2;
   hit.used = true;
   //hit.strawhits = _segmentHits.at(isegment).at(i).strawhits;
   const ComboHit* ch = &_data.chcol->at(_segmentHits.at(isegment).at(i).hitIndice);
   hit.strawhits = ch->_nsh;
   //std::cout<<"hit.x = "<<hit.x<<std::endl;
   //std::cout<<"hit.y = "<<hit.y<<std::endl;
   //std::cout<<"hit.z = "<<hit.z<<std::endl;
   //std::cout<<"hit.station = "<<hit.station<<std::endl;
   //std::cout<<"hit.plane = "<<hit.plane<<std::endl;
   //std::cout<<"hit.face = "<<hit.face<<std::endl;
   //std::cout<<"hit.panel = "<<hit.panel<<std::endl;
   //std::cout<<"hit.strawhits = "<<hit.strawhits<<std::endl;
   _tcHits.push_back(hit);
  }
   //std::cout<<"kitagawa_checkdayo = "<<std::endl;
  // After filling _tcHits, sort by z ascending:
  std::sort(_tcHits.begin(), _tcHits.end(),
          [](const cHit& a, const cHit& b) {
              return a.z < b.z;
          });
 // Print the sorted _tcHits to crosscheck
//std::cout << "Sorted _tcHits by z:\n";
//for (size_t i = 0; i < _tcHits.size(); ++i) {
//    std::cout << "Index " << i
//              << ": hitIndice = " << _tcHits[i].hitIndice
//              << ": z = " << _tcHits[i].z
//              << ", x = " << _tcHits[i].x
//              << ", y = " << _tcHits[i].y
//              << ", phi = " << _tcHits[i].phi
//              << ", station = " << _tcHits[i].station
//              << std::endl;
//}
}
void PhiZSeedFinder::tcHitsFill_Add(int isegment) {
  if(_tcHits.size() == 0)return;
  else{
  cHit hit;
  hit.circleError2 = 1.0;
  hit.helixPhi = 0.0;
  hit.helixPhiError2 = 0.0;
  hit.helixPhiCorrection = 0.0;
  hit.ambigPhi = 0.0;
  hit.inHelix = false;
  hit.used = false;//default is false
  hit.isolated = false;
  hit.averagedOut = false;
  hit.notOnLine = true;
  hit.uselessTripletSeed = false;
  hit.notOnSegment = true;
  hit.station = 0;
  hit.plane = 0;
  hit.face = 0;
  hit.panel = 0;
  hit.x = 0.0;
  hit.y = 0.0;
  hit.z = 0.0;
  hit.phi = 0.0;
  hit.strawhits = 0;
  // fill hits from time cluster of i-th segment
  //std::cout<<"====== tcHitsFill ======"<<std::endl;
  //std::cout<< "run: " << run << " subRun: " << subrun << " event: " << eventNumber<<std::endl;
  //std::cout<<"isegment = "<<isegment<<std::endl;
  //std::cout<<"isegment.size = "<<_segmentHits.at(isegment).size()<<std::endl;
  for (size_t i = 0; i < _segmentHits.at(isegment).size(); i++) {
   //std::cout<<"i-th = "<<i<<std::endl;
  // std::cout<<"hitIndice = "<<_segmentHits.at(isegment).at(i).hitIndice<<std::endl;
   hit.hitIndice = _segmentHits.at(isegment).at(i).hitIndice;
   hit.station = _segmentHits.at(isegment).at(i).station;
   hit.plane = _segmentHits.at(isegment).at(i).plane;
   hit.face = _segmentHits.at(isegment).at(i).face;
   hit.panel = _segmentHits.at(isegment).at(i).panel;
   hit.x = _segmentHits.at(isegment).at(i).x;
   hit.y = _segmentHits.at(isegment).at(i).y;
   hit.z = _segmentHits.at(isegment).at(i).z;
   hit.phi = _segmentHits.at(isegment).at(i).phi;
   hit.used = true;
   //hit.strawhits = _segmentHits.at(isegment).at(i).strawhits;
   const ComboHit* ch = &_data.chcol->at(_segmentHits.at(isegment).at(i).hitIndice);
   hit.strawhits = ch->_nsh;
   //std::cout<<"hit.x = "<<hit.x<<std::endl;
   //std::cout<<"hit.y = "<<hit.y<<std::endl;
   //std::cout<<"hit.z = "<<hit.z<<std::endl;
   //std::cout<<"hit.station = "<<hit.station<<std::endl;
   //std::cout<<"hit.plane = "<<hit.plane<<std::endl;
   //std::cout<<"hit.face = "<<hit.face<<std::endl;
   //std::cout<<"hit.panel = "<<hit.panel<<std::endl;
   //std::cout<<"hit.strawhits = "<<hit.strawhits<<std::endl;
   _tcHits.push_back(hit);
  }
  // After filling _tcHits, sort by z ascending:
  std::sort(_tcHits.begin(), _tcHits.end(),
          [](const cHit& a, const cHit& b) {
              return a.z < b.z;
          });
  }
}
//-----------------------------------------------------------------------------
// start with initial seed circle
//-----------------------------------------------------------------------------
void PhiZSeedFinder::initSeedCircle(int& outcome) {
  // get triplet circle parameters then clear fitter
  double xC = _circleFitter.x0();
  double yC = _circleFitter.y0();
  double rC = _circleFitter.radius();
  _circleFitter.clear();
  std::cout<<"xC/yC/rC = "<<xC<<"/"<<yC<<"/"<<rC<<std::endl;
  // project error bars onto the triplet circle found and add to fitter those within defined max
  // residual
}
//-----------------------------------------------------------------------------
// compute phi relative to helix center, and set helixPhiError
//-----------------------------------------------------------------------------
void PhiZSeedFinder::computeHelixPhi(size_t& tcHitsIndex, double& xC, double& yC) {
  int hitIndice = _tcHits[tcHitsIndex].hitIndice;
  //std::cout<<"x/y = "<<_tcHits.at(tcHitsIndex).x<<"/"<<_tcHits.at(tcHitsIndex).y<<std::endl;
  double X = _tcHits.at(tcHitsIndex).x - xC;
  double Y = _tcHits.at(tcHitsIndex).y - yC;
  //std::cout<<"X/Y = "<<X<<"/"<<Y<<std::endl;
  _tcHits[tcHitsIndex].helixPhi = polyAtan2(Y, X);
  //std::cout<<"polyAtan2(Y, X) = "<<polyAtan2(Y, X)<<std::endl;
  if (_tcHits[tcHitsIndex].helixPhi < 0) {
    //_tcHits[tcHitsIndex].helixPhi = _tcHits[tcHitsIndex].helixPhi + 2 * 3.14;
  }
  //std::cout<<"helixPhi = "<<_tcHits[tcHitsIndex].helixPhi<<std::endl;
  // find phi error and initialize it
  // find phi error by projecting errors onto vector tangent to circle
  double deltaS2(0);
  double tanVecX = Y / std::sqrt(X * X + Y * Y);
  double tanVecY = -X / std::sqrt(X * X + Y * Y);
  double wireErr = _data.chcol->at(hitIndice).wireRes();
  double wireVecX = _data.chcol->at(hitIndice).uDir().x();
  double wireVecY = _data.chcol->at(hitIndice).uDir().y();
  double projWireErr = wireErr * (wireVecX * tanVecX + wireVecY * tanVecY);
  double transErr = _data.chcol->at(hitIndice).transRes();
  double transVecX = _data.chcol->at(hitIndice).uDir().y();
  double transVecY = -_data.chcol->at(hitIndice).uDir().x();
  double projTransErr = transErr * (transVecX * tanVecX + transVecY * tanVecY);
  deltaS2 = projWireErr * projWireErr + projTransErr * projTransErr;
  _tcHits[tcHitsIndex].helixPhiError2 = deltaS2 / (X * X + Y * Y);
  double constant = 0.01;
  if(_tcHits[tcHitsIndex].helixPhiError2 < 0.01) _tcHits[tcHitsIndex].helixPhiError2 = _tcHits[tcHitsIndex].helixPhiError2 + constant;
  //std::cout<<"tcHitsIndex = "<<tcHitsIndex<<std::endl;
  //std::cout<<"Error2/weight = "<<_tcHits[tcHitsIndex].helixPhiError2<<"/"<<1.0/(_tcHits[tcHitsIndex].helixPhiError2)<<std::endl;
}
//-----------------------------------------------------------------------------
// compute phi relative to helix center, and set helixPhiError
//-----------------------------------------------------------------------------
void PhiZSeedFinder::computeHelixPhi_ver2(int hitIndice, double& xC, double& yC, double& helixPhi, double& helixPhiError2) {
  double X = _data.chcol->at(hitIndice).pos().x() - xC;
  double Y = _data.chcol->at(hitIndice).pos().y() - yC;
  helixPhi = polyAtan2(Y, X);
  // find phi error and initialize it
  // find phi error by projecting errors onto vector tangent to circle
  double deltaS2(0);
  double tanVecX = Y / std::sqrt(X * X + Y * Y);
  double tanVecY = -X / std::sqrt(X * X + Y * Y);
  double wireErr = _data.chcol->at(hitIndice).wireRes();
  double wireVecX = _data.chcol->at(hitIndice).uDir().x();
  double wireVecY = _data.chcol->at(hitIndice).uDir().y();
  double projWireErr = wireErr * (wireVecX * tanVecX + wireVecY * tanVecY);
  double transErr = _data.chcol->at(hitIndice).transRes();
  double transVecX = _data.chcol->at(hitIndice).uDir().y();
  double transVecY = -_data.chcol->at(hitIndice).uDir().x();
  double projTransErr = transErr * (transVecX * tanVecX + transVecY * tanVecY);
  deltaS2 = projWireErr * projWireErr + projTransErr * projTransErr;
  helixPhiError2 = deltaS2 / (X * X + Y * Y);
  double constant = 0.01;
  if(helixPhiError2 < 0.01) helixPhiError2 = helixPhiError2 + constant;
}
//-----------------------------------------------------------------------------
// function to initialize phi info relative to helix center in _tcHits
//-----------------------------------------------------------------------------
void PhiZSeedFinder::initHelixPhi() {
  double xC = _circleFitter.x0();
  double yC = _circleFitter.y0();
  size_t nComboHitsInSegment = _tcHits.size();
  for (size_t i = 0; i < nComboHitsInSegment; i++) {
    computeHelixPhi(i, xC, yC);
    // initialize phi data member relative to circle center
    _tcHits[i].helixPhiCorrection = 0;
  }
}
//-----------------------------------------------------------------------------
void PhiZSeedFinder::plot_PhiVsZ_OriginalTC(int tc){
  const int n = (int)_data.tccol->at(tc)._strawHitIdxs.size();
  TGraph *gr = new TGraph(n);
  gr->SetTitle("");
  gr->SetMarkerStyle(1);
  double phi, z;
  std::vector<int> marker_color;
  std::vector<int> marker_style;
  std::vector<double> marker_size;
  for (int i = 0; i < n; i++) {
    int hitIndice = _data.tccol->at(tc)._strawHitIdxs[i];
    std::vector<StrawDigiIndex> shids;
    _data.chcol->fillStrawDigiIndices(hitIndice, shids);
    z = _data.chcol->at(hitIndice).pos().z();
    phi = _data.chcol->at(hitIndice).pos().phi();
    std::cout<<"z/phi = "<<z<<"/"<<phi<<std::endl;
    //loop over StrawHits(shids)
    for (size_t j = 0; j < shids.size(); j++) {
      const mu2e::SimParticle* _simParticle;
      _simParticle = _mcUtils->getSimParticle(_event, shids[j]);
      //int SimID = _mcUtils->strawHitSimId(_event, shids[j]);
      int PdgID = _simParticle->pdgId();
      int   style(0), color(0);
      double size(0.);
      if (PdgID ==  11) {style = 20; size = 0.8; color = kRed;}
      else if (PdgID ==  -11) {style = 24; size = 0.8; color = kBlue;}
      else if (PdgID ==  13) {style = 20; size = 0.8; color = kGreen+2;}
      else if (PdgID ==  -13) {style = 20; size = 0.8; color = kGreen-2;}
      else if (PdgID ==  2212) {style = 20; size = 0.8; color = kBlue+2;}
      else if (PdgID ==  -211) {style = 20; size = 0.8; color = kPink+2;}
      else if (PdgID ==  211) {style = 20; size = 0.8; color = kPink-2;}
      else {style = 20; size = 0.8; color = kMagenta;}
      marker_color.push_back(color);
      marker_style.push_back(style);
      marker_size.push_back(size);
      break;
    }
    gr->SetPoint(i, z, phi);
   }
  //Draw
  TCanvas *canvas = new TCanvas("canvas", "", 800, 600);
  canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
   gr->Draw("AP");
   // Add multiple colors logic here
   std::cout<<"Draw "<<std::endl;
   TMarker *m;
   for (int i = 0; i < n; i++) {
    gr->GetPoint(i, z, phi);
    std::cout<<"z/phi = "<<z<<"/"<<phi<<std::endl;
    m = new TMarker(z, phi, 20);
    //m->SetMarkerColor(i + 1);// Setting marker color with different color for each point
    m->SetMarkerStyle(marker_style[i]);
    m->SetMarkerSize(marker_size[i]);
    m->SetMarkerColor(marker_color[i]);// Setting marker color with different color for each point
    m->Draw();// Draw marker with different color
   }
  gr->GetXaxis()->SetTitle("Z [mm]");
  gr->GetYaxis()->SetTitle("Helix Phi [rad]");
  gr->GetXaxis()->SetLimits(-1600, 1600);
  gr->GetYaxis()->SetRangeUser(-M_PI, M_PI);
  // Add title at the top
  TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
  std::stringstream eventStringStream;
  eventStringStream << "run: " << run << " subRun: " << subrun << " event: " << eventNumber;
  title->AddText(Form("Phi vs. Z (Run-subRun-Event, TC) = (%d-%d-%d, #%d) (pbar1b0)", run, subrun, eventNumber, tc));
  title->SetFillColor(0);
  title->SetTextAlign(22);
  title->Draw("same");
  //Draw vertical lines at specified X-coordinate of stations (18 stations)
  double stations_x[18] = {-1518.320, -1344.320, -1170.320, -996.320, -822.320, -648.320, -474.320, -300.320, -126.320, 47.680, 221.680, 395.680, 569.680, 743.680, 917.680, 1091.680, 1265.680, 1439.680};
  TLine *line[18];
  for(int j=0; j<18; j++){
    double x = stations_x[j];
    line[j] = new TLine(x, -M_PI, x, M_PI);
    line[j]->SetLineStyle(2);  // Dashed line style
    line[j]->SetLineWidth(1);
    line[j]->SetLineColor(kBlack);
    line[j]->Draw("same");
  }
   canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/oroginal_TC/pbar_PhiVsZ-%04d-%04d-%06d_TC-%d.pdf", run, subrun, eventNumber, tc));
   //canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/ce/PhiVsZ/oroginal_TC/pbar_PhiVsZ-%04d-%04d-%06d_TC-%d.pdf", run, subrun, eventNumber, tc));
  //delete
  delete canvas;
  delete gr;
  for(int j=0; j<18; j++) delete line[j];
}
//-----------------------------------------------------------------------------
void PhiZSeedFinder::plot_HelixPhiVsZ(int TC, int isegment){
  //step2: calculate the slope for other segment
  ::LsqSums2 _lineFitter;
  _lineFitter.clear();
  TMultiGraph* graph  = new TMultiGraph();
  const int ngraph = 2;
  TGraph *gr[ngraph];
  for(int j=0; j<ngraph; j++) {
    gr[j] = new TGraph();
    gr[j]->SetMarkerStyle(20);
    gr[j]->SetMarkerSize(0.5);
  }
  gr[0]->SetMarkerColor(kRed);
  gr[1]->SetMarkerColor(kGreen);
  //Helix Phi
  int index[ngraph] = {0};
  for(size_t j=0; j<_tcHits.size(); j++) {
   if(_tcHits[j].used == false) continue;
   double z = _tcHits[j].z;
   double phi = _tcHits[j].phi;
   double helixPhi = _tcHits[j].helixPhi;
   //std::cout<<"Helix Phi z/phi = "<<j<<"/"<<z<<"/"<<phi<<std::endl;
   gr[0]->SetPoint(index[0]++, z, helixPhi);
   gr[1]->SetPoint(index[1]++, z, phi);
   double xC = 0.0;
   double yC = 0.0;
   computeHelixPhi(j, xC, yC);
   double phiWeight = 1.0 / (_tcHits[j].helixPhiError2);
   _lineFitter.addPoint(z, phi, phiWeight);
  }
  //graph->Add(gr[0],"AP");
  graph->Add(gr[1],"AP");
  //Draw
  TCanvas *canvas = new TCanvas("canvas", "My TGraph", 800, 600);
  canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
  graph->Draw("AP");
  graph->GetXaxis()->SetTitle("Z [mm]");
  graph->GetYaxis()->SetTitle("Helix Phi [rad]");
  graph->GetXaxis()->SetLimits(-1600, 1600);
  graph->GetYaxis()->SetRangeUser(-M_PI, M_PI);
  // Add title at the top
  TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
  std::stringstream eventStringStream;
  eventStringStream << "run: " << run << " subRun: " << subrun << " event: " << eventNumber;
  title->AddText(Form("Phi vs. Z (Run-subRun-Event, TC, Cand) = (%d-%d-%d, #%d, #%d) (pbar1b0)", run, subrun, eventNumber, TC, isegment));
  title->SetFillColor(0);
  title->SetTextAlign(22);
  title->Draw("same");
  //Draw vertical lines at specified X-coordinate of stations (18 stations)
  double stations_x[18] = {-1518.320, -1344.320, -1170.320, -996.320, -822.320, -648.320, -474.320, -300.320, -126.320, 47.680, 221.680, 395.680, 569.680, 743.680, 917.680, 1091.680, 1265.680, 1439.680};
  TLine *line[18];
  for(int j=0; j<18; j++){
    double x = stations_x[j];
    line[j] = new TLine(x, -M_PI, x, M_PI);
    line[j]->SetLineStyle(2);  // Dashed line style
    line[j]->SetLineWidth(1);
    line[j]->SetLineColor(kBlack);
    line[j]->Draw("same");
  }
  // Draw the fitted slope line: y = dydx * x + y0 from x = -1600 to 1600
  double x_min = -1600;
  double x_max = 1600;
  double y_min = _lineFitter.dydx() * x_min + _lineFitter.y0();
  double y_max = _lineFitter.dydx() * x_max + _lineFitter.y0();
  TLine *fitLine = new TLine(x_min, y_min, x_max, y_max);
  fitLine->SetLineColor(kBlack);
  fitLine->SetLineStyle(2);  // Dashed line
  fitLine->SetLineWidth(2);
  fitLine->Draw("same");
// Draw manual fit parameters text
TLatex paramsText;
paramsText.SetTextSize(0.04);
paramsText.SetTextAlign(13);
double textX = 0.15;
double textY = 0.85;
paramsText.DrawLatexNDC(textX, textY,        Form("LsqSums2 fit dydx: %f", _lineFitter.dydx()));
paramsText.DrawLatexNDC(textX, textY - 0.05, Form("LsqSums2 fit y0: %f", _lineFitter.y0()));
paramsText.DrawLatexNDC(textX, textY - 0.10, Form("LsqSums2 fit #chi^{2}/ndf: %f", _lineFitter.chi2Dof()));
   canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/HelixPhi/pbar_PhiVsZ-%04d-%04d-%06d_TC-%d_Segemnt_%d.pdf", run, subrun, eventNumber, TC, isegment));
   //canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/ce/PhiVsZ/HelixPhi/pbar_PhiVsZ-%04d-%04d-%06d_TC-%d.pdf", run, subrun, eventNumber, TC));
  //delete
  delete canvas;
  for(int j=0; j<ngraph; j++) delete gr[j];
  for(int j=0; j<18; j++) delete line[j];
}
//-----------------------------------------------------------------------------
void PhiZSeedFinder::plot_PhiVsZ_forEachStep(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, const char* filename, int tc, int loopIndex){
      //std::cout<<"======================="<<std::endl;
      //std::cout<<"plot_PhiVsZ_forEachStep"<<std::endl;
      //std::cout<<"Loop Index = "<<loopIndex<<std::endl;
      //std::cout<< filename <<std::endl;
      //std::cout<<"======================="<<std::endl;
      std::vector<std::vector<ev5_HitsInNthStation>> segments = ThisIsBestSegment;
     // std::vector<std::vector<ev5_Segment>>& diag_segments = ThisIsBestSegment_Diag;
      //Take segment
      //std::cout<<"# of segments = "<<(int)segments.size()<<std::endl;
      //for(int i=0; i<(int)segments.size(); i++){
      //  for(int j=0; j<(int)segments.at(i).size(); j++){
      //    std::cout<<"segments.at["<<j<<"] = "<<segments.at(i).at(j).phi<<std::endl;
      //  }
      //}
 for(int i=0; i<(int)segments.size(); i++){
  TGraphErrors* gr = new TGraphErrors();
  gr->SetTitle("");
  gr->SetMarkerStyle(20);
  //std::cout<<"Fill "<<std::endl;
  tcHitsFill(i);
  ::LsqSums2 fitter;
  fitter.clear();
  int index = 0;
  for(size_t j=0; j<_tcHits.size(); j++) {
    if(_tcHits[j].used == false) continue;
    //double z = segments.at(i).at(j).z;
    //double phi = segments.at(i).at(j).phi;
    double z = _tcHits[j].z;
    double phi = _tcHits[j].phi;
    //std::cout<<"z/phi = "<<z<<"/"<<phi<<std::endl;
    double xC = 0.0;
    double yC = 0.0;
    computeHelixPhi(j, xC, yC);
    double phiWeight = 1.0 / (_tcHits[j].helixPhiError2);
    fitter.addPoint(z, phi, phiWeight);
    gr->SetPoint(index, z, phi);
    gr->SetPointError(index, 0.0, std::sqrt(_tcHits[j].helixPhiError2));
    index++;
   }
   index = 0;
   _data.h_lineFitter_chi2Dof[0].push_back(fitter.chi2Dof());
  //Draw
  TCanvas *canvas = new TCanvas("canvas", "", 800, 600);
  canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
  gr->SetMarkerSize(0.8);
  gr->SetMarkerColor(kRed);
  gr->Draw("AP");
  gr->GetXaxis()->SetTitle("Z [mm]");
  gr->GetYaxis()->SetTitle("Phi [rad]");
  gr->GetXaxis()->SetLimits(-1600, 1600);
  gr->GetYaxis()->SetRangeUser(-M_PI, M_PI);
  // Add title at the top
  TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
  std::stringstream eventStringStream;
  eventStringStream << "run: " << run << " subRun: " << subrun << " event: " << eventNumber;
  title->AddText(Form("Phi vs. Z (Run-subRun-Event, #cand) = (%d-%d-%d, #%d) (pbar1b0)", run, subrun, eventNumber, i));
  title->SetFillColor(0);
  title->SetTextAlign(22);
  title->Draw("same");
  //Draw vertical lines at specified X-coordinate of stations (18 stations)
  double stations_x[18] = {-1518.320, -1344.320, -1170.320, -996.320, -822.320, -648.320, -474.320, -300.320, -126.320, 47.680, 221.680, 395.680, 569.680, 743.680, 917.680, 1091.680, 1265.680, 1439.680};
  TLine *line[18];
  for(int j=0; j<18; j++){
    double x = stations_x[j];
    line[j] = new TLine(x, -M_PI, x, M_PI);
    line[j]->SetLineStyle(2);  // Dashed line style
    line[j]->SetLineWidth(1);
    line[j]->SetLineColor(kBlack);
    line[j]->Draw("same");
  }
  // Draw the fitted slope line: y = dydx * x + y0 from x = -1600 to 1600
  double x_min = -1600;
  double x_max = 1600;
  double y_min = fitter.dydx() * x_min + fitter.y0();
  double y_max = fitter.dydx() * x_max + fitter.y0();
  TLine *fitLine = new TLine(x_min, y_min, x_max, y_max);
  fitLine->SetLineColor(kBlack);
  fitLine->SetLineStyle(2);  // Dashed line
  fitLine->SetLineWidth(2);
  if (strcmp(filename, "step_07") != 0) fitLine->Draw("same");
// Draw manual fit parameters text
TLatex paramsText;
paramsText.SetTextSize(0.04);
paramsText.SetTextAlign(13);
double textX = 0.15;
double textY = 0.85;
if (strcmp(filename, "step_07") != 0) paramsText.DrawLatexNDC(textX, textY,        Form("LsqSums2 fit dydx: %f", fitter.dydx()));
if (strcmp(filename, "step_07") != 0) paramsText.DrawLatexNDC(textX, textY - 0.05, Form("LsqSums2 fit y0: %f", fitter.y0()));
if (strcmp(filename, "step_07") != 0) paramsText.DrawLatexNDC(textX, textY - 0.10, Form("LsqSums2 fit #chi^{2}/ndf: %f", fitter.chi2Dof()));
if (strcmp(filename, "step_07") != 0) paramsText.DrawLatexNDC(textX, textY - 0.15, Form("nHits: %f", fitter.qn()));
   canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/segment_check/%s/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d-%d.pdf", filename, run, subrun, eventNumber, tc, i, loopIndex));
   //canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/ce/PhiVsZ/segment_check/%s/pbar_PhiVsZ-%04d-%04d-%06d_-%d-%d.pdf", filename, run, subrun, eventNumber, i, loopIndex));
  //delete
  delete canvas;
  delete gr;
  for(int j=0; j<18; j++) delete line[j];
  }
}
//-----------------------------------------------------------------------------
void PhiZSeedFinder::plot_PhiVsZ_forSegment(int tc, int isegment){
  const int n = (int)_segmentHits.at(isegment).size();
  TGraph *gr = new TGraph(n);
  gr->SetTitle("");
  gr->SetMarkerStyle(1);
  double phi, z;
  std::vector<int> marker_color;
  std::vector<int> marker_style;
  std::vector<double> marker_size;
    std::cout<<"Fill "<<std::endl;
  for(int i=0; i<n; i++) {
   double z = _segmentHits.at(isegment).at(i).z;
   double phi = _segmentHits.at(isegment).at(i).phi;
   int index = _segmentHits.at(isegment).at(i).hitIndice;
   int alreadyfill = 0;
   for (size_t j = 0; j < _data.tccol->at(tc)._strawHitIdxs.size(); j++) {
    int hitIndice = _data.tccol->at(tc)._strawHitIdxs[j];
    if(index != hitIndice) continue;
    if(alreadyfill == 1) continue;
    std::vector<StrawDigiIndex> shids;
    _data.chcol->fillStrawDigiIndices(hitIndice, shids);
    //z = _data.chcol->at(hitIndice).pos().z();
    //phi = _data.chcol->at(hitIndice).pos().phi();
    std::cout<<"z/phi/station = "<<z<<"/"<<phi<<"/"<<_segmentHits.at(isegment).at(i).station<<std::endl;
    for (size_t k = 0; k < shids.size(); k++) {
      const mu2e::SimParticle* _simParticle;
      _simParticle = _mcUtils->getSimParticle(_event, shids[k]);
      //int SimID = _mcUtils->strawHitSimId(_event, shids[j]);
      int PdgID = _simParticle->pdgId();
      int   style(0), color(0);
      double size(0.);
      if (PdgID ==  11) {style = 20; size = 0.8; color = kRed;}
      else if (PdgID ==  -11) {style = 24; size = 0.8; color = kBlue;}
      else if (PdgID ==  13) {style = 20; size = 0.8; color = kGreen+2;}
      else if (PdgID ==  -13) {style = 20; size = 0.8; color = kGreen-2;}
      else if (PdgID ==  2212) {style = 20; size = 0.8; color = kBlue+2;}
      else if (PdgID ==  -211) {style = 20; size = 0.8; color = kPink+2;}
      else if (PdgID ==  211) {style = 20; size = 0.8; color = kPink-2;}
      else {style = 20; size = 0.8; color = kMagenta;}
      marker_color.push_back(color);
      marker_style.push_back(style);
      marker_size.push_back(size);
      alreadyfill = 1;
      break;
    }
   }
    gr->SetPoint(i, z, phi);
   }
  //Draw
  TCanvas *canvas = new TCanvas("canvas", "", 800, 600);
  canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
  gr->Draw("AP");
  TMarker *m;
    std::cout<<"Draw "<<std::endl;
  for (int i = 0; i < n; i++) {
    gr->GetPoint(i, z, phi);
    std::cout<<"z/phi = "<<z<<"/"<<phi<<std::endl;
    m = new TMarker(z, phi, 20);
    //m->SetMarkerColor(i + 1);// Setting marker color with different color for each point
    m->SetMarkerStyle(marker_style[i]);
    m->SetMarkerSize(marker_size[i]);
    m->SetMarkerColor(marker_color[i]);// Setting marker color with different color for each point
    m->Draw();// Draw marker with different color
  }
  gr->GetXaxis()->SetTitle("Z [mm]");
  gr->GetYaxis()->SetTitle("Phi [rad]");
  gr->GetXaxis()->SetLimits(-1600, 1600);
  gr->GetYaxis()->SetRangeUser(-M_PI, M_PI);
  // Add title at the top
  TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
  std::stringstream eventStringStream;
  eventStringStream << "run: " << run << " subRun: " << subrun << " event: " << eventNumber;
  title->AddText(Form("Phi vs. Z (Run-subRun-Event, TC, #cand) = (%d-%d-%d, %d, #%d) (pbar1b0)", run, subrun, eventNumber, tc, isegment));
  title->SetFillColor(0);
  title->SetTextAlign(22);
  title->Draw("same");
  //Draw vertical lines at specified X-coordinate of stations (18 stations)
  double stations_x[18] = {-1518.320, -1344.320, -1170.320, -996.320, -822.320, -648.320, -474.320, -300.320, -126.320, 47.680, 221.680, 395.680, 569.680, 743.680, 917.680, 1091.680, 1265.680, 1439.680};
  TLine *line[18];
  for(int j=0; j<18; j++){
    double x = stations_x[j];
    line[j] = new TLine(x, -M_PI, x, M_PI);
    line[j]->SetLineStyle(2);  // Dashed line style
    line[j]->SetLineWidth(1);
    line[j]->SetLineColor(kBlack);
    line[j]->Draw("same");
  }
   canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/segment/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));
   //canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/ce/PhiVsZ/segment/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));
  //delete
  delete canvas;
  delete gr;
  for(int j=0; j<18; j++) delete line[j];
}
//-----------------------------------------------------------------------------
void PhiZSeedFinder::plot_PhiVsZ_forSegment_ver2(int ith_segment, int jth_segment, double alpha, double beta, double Chi2NDF){
  TGraphErrors* gr = new TGraphErrors();
  gr->SetTitle("");
  gr->SetMarkerStyle(20);
  /*for(int i=0; i<(int)_tcHits.size(); i++) {
    double z = _tcHits[i].z;
    double phi = _tcHits[i].helixPhi;
    gr->SetPoint(i, z, phi);
    gr->SetPointError(i, 0.0, std::sqrt(_tcHits[i].helixPhiError2));
  }*/
    std::cout<<" plot_PhiVsZ_forSegment_ver2 "<<std::endl;
    for(size_t j=0; j<_segmentHits.at(ith_segment).size(); j++){
      double z = _segmentHits.at(ith_segment).at(j).z;
      int turns = _segmentHits.at(ith_segment).at(j).nturn;
      double phi = _segmentHits.at(ith_segment).at(j).helixPhi + turns * 2 * M_PI;
      std::cout<<"z/turns/phi = "<<z<<"/"<<turns<<"/"<<_segmentHits.at(ith_segment).at(j).helixPhi<<std::endl;
      gr->SetPoint(j, z, phi);
      gr->SetPointError(j, 0.0, std::sqrt(_segmentHits.at(ith_segment).at(j).helixPhiError2));
    }
  //Draw
  TCanvas *canvas = new TCanvas("canvas", "", 800, 600);
  canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
  gr->SetMarkerSize(0.8);
  gr->SetMarkerColor(kRed);
  gr->Draw("AP");
  gr->GetXaxis()->SetTitle("Z [mm]");
  gr->GetYaxis()->SetTitle("Phi [rad]");
  gr->GetXaxis()->SetLimits(-1600, 1600);
  //gr->GetYaxis()->SetRangeUser(-M_PI, M_PI);
  gr->GetYaxis()->SetRangeUser(-3*M_PI, 3*M_PI);
  // Add title at the top
  TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
  std::stringstream eventStringStream;
  eventStringStream << "run: " << run << " subRun: " << subrun << " event: " << eventNumber;
  title->AddText(Form("Phi vs. Z (Run-subRun-Event, #cand) = (%d-%d-%d, #%d - %d) (pbar1b0)", run, subrun, eventNumber, ith_segment, jth_segment));
  title->SetFillColor(0);
  title->SetTextAlign(22);
  title->Draw("same");
  //Draw vertical lines at specified X-coordinate of stations (18 stations)
  double stations_x[18] = {-1518.320, -1344.320, -1170.320, -996.320, -822.320, -648.320, -474.320, -300.320, -126.320, 47.680, 221.680, 395.680, 569.680, 743.680, 917.680, 1091.680, 1265.680, 1439.680};
  TLine *line[18];
  for(int j=0; j<18; j++){
    double x = stations_x[j];
    //line[j] = new TLine(x, -M_PI, x, M_PI);
    line[j] = new TLine(x, -3*M_PI, x, 3*M_PI);
    line[j]->SetLineStyle(2);  // Dashed line style
    line[j]->SetLineWidth(1);
    line[j]->SetLineColor(kBlack);
    line[j]->Draw("same");
  }
  // Draw the fitted slope line: y = dydx * x + y0 from x = -1600 to 1600
  double x_min = -1600;
  double x_max = 1600;
  double y_min = alpha * x_min + beta;
  double y_max = alpha * x_max + beta;
  TLine *fitLine = new TLine(x_min, y_min, x_max, y_max);
  fitLine->SetLineColor(kBlack);
  fitLine->SetLineStyle(2);  // Dashed line
  fitLine->SetLineWidth(2);
  fitLine->Draw("same");
  // Draw manual fit parameters text
  TLatex paramsText;
  paramsText.SetTextSize(0.04);
  paramsText.SetTextAlign(13);
  double textX = 0.15;
  double textY = 0.85;
  paramsText.DrawLatexNDC(textX, textY,        Form("LsqSums2 fit dydx: %f", alpha));
  paramsText.DrawLatexNDC(textX, textY - 0.05, Form("LsqSums2 fit y0: %f", beta));
  paramsText.DrawLatexNDC(textX, textY - 0.10, Form("LsqSums2 fit #chi^{2}/ndf: %f", Chi2NDF));
  paramsText.DrawLatexNDC(textX, textY - 0.15, Form("nHits: %d", (int)_tcHits.size()));
   canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/segment/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, ith_segment, jth_segment));
  //delete
  delete canvas;
  delete gr;
  for(int j=0; j<18; j++) delete line[j];
}
//-----------------------------------------------------------------------------
void PhiZSeedFinder::plot_CirclePhiVsZ_forSegment(int tc, int isegment){
    _circleFitter.clear();
    //Step1: fit circle of i-th segment with fixed weight
    size_t nComboHitsInSegment = _tcHits.size();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 0.1;//tentative value
      _circleFitter.addPoint(x, y, wP);
      //_tcHits[i].used = true;
    }
    double xC = _circleFitter.x0();
    double yC = _circleFitter.y0();
    double rC = _circleFitter.radius();
    //Step2: fit circle of i-th segment with correct weight
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    _circleFitter.clear();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      std::cout<<"i = "<<i<<std::endl;
      computeCircleError2(i, xC, yC, rC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
      //_tcHits[i].used = true;
    }
    //Step3: fit circle of i-th segment with correct weight
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    _circleFitter.clear();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      computeCircleError2(i, xC, yC, rC);
      computeHelixPhi(i, xC, yC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
      //_tcHits[i].used = true;
    }
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      computeHelixPhi(i, xC, yC);
    }
    _circleFitter.clear();
  TGraph *gr = new TGraph();
  gr->SetTitle("");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.8);
  gr->SetMarkerColor(kRed);
  for(int i=0; i<(int) _segmentHits.at(isegment).size(); i++) {
    double z = _segmentHits.at(isegment).at(i).z;
    double phi_ = _segmentHits.at(isegment).at(i).phi;
    double phi = _tcHits[i].helixPhi;
    std::cout<<"phi_0/phi = "<<phi_<<"/"<<phi<<std::endl;
    gr->SetPoint(i, z, phi);
  }
  //Draw
  TCanvas *canvas = new TCanvas("canvas", "", 800, 600);
  canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
  gr->Draw("AP");
  gr->GetXaxis()->SetTitle("Z [mm]");
  gr->GetYaxis()->SetTitle("Phi [rad]");
  gr->GetXaxis()->SetLimits(-1600, 1600);
  gr->GetYaxis()->SetRangeUser(-M_PI, M_PI);
  // Add title at the top
  TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
  std::stringstream eventStringStream;
  eventStringStream << "run: " << run << " subRun: " << subrun << " event: " << eventNumber;
  title->AddText(Form("Phi vs. Z (Run-subRun-Event, TC, #cand) = (%d-%d-%d, %d, #%d) (pbar1b0)", run, subrun, eventNumber, tc, isegment));
  title->SetFillColor(0);
  title->SetTextAlign(22);
  title->Draw("same");
  //Draw vertical lines at specified X-coordinate of stations (18 stations)
  double stations_x[18] = {-1518.320, -1344.320, -1170.320, -996.320, -822.320, -648.320, -474.320, -300.320, -126.320, 47.680, 221.680, 395.680, 569.680, 743.680, 917.680, 1091.680, 1265.680, 1439.680};
  TLine *line[18];
  for(int j=0; j<18; j++){
    double x = stations_x[j];
    line[j] = new TLine(x, -M_PI, x, M_PI);
    line[j]->SetLineStyle(2);  // Dashed line style
    line[j]->SetLineWidth(1);
    line[j]->SetLineColor(kBlack);
    line[j]->Draw("same");
  }
   canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/segment_circle/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));
   //canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/ce/PhiVsZ/segment_circle/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));
  //delete
  delete canvas;
  delete gr;
  for(int j=0; j<18; j++) delete line[j];
}
//-----------------------------------------------------------------------------
void PhiZSeedFinder::plot_2PiAmbiguityPhiVsZ_forSegment(int tc, int isegment){
    std::cout<<"---------------------------------"<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    std::cout<<" 2PiAmbiguityPhiVsZ_forSegment "<<isegment<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    _circleFitter.clear();
    //Step1: fit circle of i-th segment with fixed weight
    size_t nComboHitsInSegment = _tcHits.size();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 0.1;//tentative value
      _circleFitter.addPoint(x, y, wP);
      //_tcHits[i].used = true;
    }
    double xC = _circleFitter.x0();
    double yC = _circleFitter.y0();
    double rC = _circleFitter.radius();
    //Step2: fit circle of i-th segment with correct weight
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    _circleFitter.clear();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      //std::cout<<"i = "<<i<<std::endl;
      computeCircleError2(i, xC, yC, rC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
      //_tcHits[i].used = true;
    }
    //Step3: fit circle of i-th segment with correct weight
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    _circleFitter.clear();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      computeCircleError2(i, xC, yC, rC);
      computeHelixPhi(i, xC, yC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
      //_tcHits[i].used = true;
    }
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      computeHelixPhi(i, xC, yC);
    }
    _circleFitter.clear();
  TGraph *gr = new TGraph();
  gr->SetTitle("");
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(0.8);
  gr->SetMarkerColor(kRed);
  for(int i=0; i<(int) _segmentHits.at(isegment).size(); i++) {
    double z = _segmentHits.at(isegment).at(i).z;
    double phi_ = _segmentHits.at(isegment).at(i).phi;
    double phi = _tcHits[i].helixPhi;
    std::cout<<"i/ z/phi_0/phi = "<<i<<"/"<<z<<"/"<<phi_<<"/"<<phi<<std::endl;
    gr->SetPoint(i, z, phi);
  }
    // sort all_BestSegmentInfo in increasing order station
    for (int i = 0; i < (int)all_BestSegmentInfo.size(); i++) {
      std::sort(all_BestSegmentInfo.at(i).begin(), all_BestSegmentInfo.at(i).end(),
              [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b) {
                  return a.station < b.station; // Ascending order by 'station'
              });
    }
  //step1: calculate the slope from the fisrt segment
  int minimum_index = 0;
  int minimum_station = 17;
  for(int i=0; i<(int)all_BestSegmentInfo.size(); i++){
    for(int j=0; j<(int)all_BestSegmentInfo.at(i).size(); j++){
      if(isegment != all_BestSegmentInfo.at(i).at(j).segmentIndex) continue;
      int station = all_BestSegmentInfo.at(i).at(j).station;
      if(minimum_station > station){
        minimum_station = station;
        minimum_index = i;
      }
    }
  }
  std::cout<<"minimum_index = "<<minimum_index<<std::endl;
  ::LsqSums2 fitter;
  std::vector<double> fitter_phi;
  std::vector<double> fitter_z;
  fitter.clear();
  for(int i=0; i<(int)all_BestSegmentInfo.size(); i++){
    std::cout<<"i = "<<i<<std::endl;
    if(i != minimum_index) continue;
    int already = 0;
    //int already = 0;
    for(int j=0; j<(int)all_BestSegmentInfo.at(i).size(); j++){
      if(isegment != all_BestSegmentInfo.at(i).at(j).segmentIndex) continue;
      already = 1;
      // z as x-axis and phi as y-axis
      double z = all_BestSegmentInfo.at(i).at(j).z;
      double phi = 0.0;
      int hitIndice = all_BestSegmentInfo.at(i).at(j).hitIndice;
      //get helixPhi
      for(int k=0; k<(int)_tcHits.size(); k++){
        if(hitIndice != _tcHits[k].hitIndice) continue;
        phi = _tcHits[k].helixPhi;
      }
      double seedError2 = 0.1;
      double seedWeight = 1.0/(seedError2);
      fitter.addPoint(z, phi, seedWeight);
      fitter_z.push_back(z);
      fitter_phi.push_back(phi);
      std::cout<<"z/phi = "<<z<<"/"<<phi<<std::endl;
    }
    if(already == 1) break;
  }
      //return fitting values
      double dphidz = fitter.dydx();
      double alpha = dphidz;
      double beta = fitter.y0();
      double chindf = fitter.chi2Dof();
      _dphidz = fitter.dydx();
      _fz0 = fitter.y0();
      std::cout<<"dphidz = "<<dphidz<<std::endl;
      std::cout<<"alpha = "<<alpha<<std::endl;
      std::cout<<"beta = "<<beta<<std::endl;
      std::cout<<"chindf = "<<chindf<<std::endl;
  //step2: calculate the slope for other segment
  ::LsqSums2 fitter_2;
  fitter_2.clear();
  std::vector<int> index;
  std::vector<int> loop_index;
  for(int i=0; i<(int)all_BestSegmentInfo.size(); i++){
    //skip 1st segment
    if(i == minimum_index) continue;
    std::cout<<"i = "<<i<<std::endl;
    int best_loop = 0;
    int already = 0;
    double diff = 99999;
    int loop = 16;
    for(int k=0; k<=loop; k++){
    std::cout<<"loop = "<<k<<std::endl;
    fitter_2.clear();
    for(int l=0; l<(int)fitter_z.size(); l++){
        double seedError2 = 0.1;
        double seedWeight = 1.0/(seedError2);
        double z = fitter_z[l];
        double phi = fitter_phi[l];
       fitter_2.addPoint(z, phi, seedWeight);
    }
    std::cout<<"fitter_2.dydx() = "<<fitter_2.dydx()<<std::endl;
    for(int j=0; j<(int)all_BestSegmentInfo.at(i).size(); j++){
      if(isegment != all_BestSegmentInfo.at(i).at(j).segmentIndex) continue;
      already = 1;
      // z as x-axis and phi as y-axis
      double z = all_BestSegmentInfo.at(i).at(j).z;
      double phi = 0.0;
      int hitIndice = all_BestSegmentInfo.at(i).at(j).hitIndice;
      //get helixPhi
      for(int k=0; k<(int)_tcHits.size(); k++){
        if(hitIndice != _tcHits[k].hitIndice) continue;
        phi = _tcHits[k].helixPhi;
      }
      double seedError2 = 0.1;
      double seedWeight = 1.0/(seedError2);
      phi = phi + M_PI*(double)(k-8);
      std::cout<<"z/phi = "<<z<<"/"<<phi<<std::endl;
      fitter_2.addPoint(z, phi, seedWeight);
    }
    if(already == 0) continue;
    double dphidz_diff = abs(dphidz - fitter_2.dydx());
    std::cout<<"dphidz_diff/dphidz/fitter_2.dydx() = "<<dphidz_diff<<"/"<<dphidz<<"/"<<fitter_2.dydx()<<std::endl;
    if(dphidz_diff < diff) {
        best_loop = k-8;
        diff = dphidz_diff;
    }
    std::cout<<" best_loop = "<<best_loop<<std::endl;
    }//end loop
    //found segment
    if(already != 1)  continue;
    index.push_back(i);
    loop_index.push_back(best_loop);
  }
    //reset _tcHits
    for(size_t i=0; i<nComboHitsInSegment; i++){
      _tcHits[i].used = false;
    }
  index.push_back(minimum_index);
  loop_index.push_back(0);
  for(int i=0; i<(int)index.size(); i++){
    std::cout<<"index/best_loop = "<<index[i]<<"/"<<loop_index[i]<<std::endl;
  }
  int gr_index = 0;
  for(int m=0; m<(int)index.size(); m++){
  ::LsqSums2 fitter_3;
  fitter_3.clear();
  for(int i=0; i<(int)all_BestSegmentInfo.size(); i++){
    if(i != index[m]) continue;
    for(int j=0; j<(int)all_BestSegmentInfo.at(i).size(); j++){
      if(isegment != all_BestSegmentInfo.at(i).at(j).segmentIndex) continue;
      // z as x-axis and phi as y-axis
      double z = all_BestSegmentInfo.at(i).at(j).z;
      double phi = 0.0;
      int hitIndice = all_BestSegmentInfo.at(i).at(j).hitIndice;
      //get helixPhi
      for(int k=0; k<(int)_tcHits.size(); k++){
        if(hitIndice != _tcHits[k].hitIndice) continue;
        _tcHits[k].used = true;
        phi = _tcHits[k].helixPhi;
      }
      phi = phi + M_PI*(double)(loop_index[m]);
      gr->SetPoint(gr_index++, z, phi);
      double seedError2 = 0.1;
      double seedWeight = 1.0/(seedError2);
      fitter_3.addPoint(z, phi, seedWeight);
    }
  }
  }
  //Draw
  TCanvas *canvas = new TCanvas("canvas", "", 800, 600);
  canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
  gr->Draw("AP");
  gr->GetXaxis()->SetTitle("Z [mm]");
  gr->GetYaxis()->SetTitle("Phi [rad]");
  gr->GetXaxis()->SetLimits(-1600, 1600);
  gr->GetYaxis()->SetRangeUser(-8*M_PI, 8*M_PI);
  // Add title at the top
  TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
  std::stringstream eventStringStream;
  eventStringStream << "run: " << run << " subRun: " << subrun << " event: " << eventNumber;
  title->AddText(Form("Phi vs. Z (Run-subRun-Event, TC, #cand) = (%d-%d-%d, %d, #%d) (pbar1b0)", run, subrun, eventNumber, tc, isegment));
  title->SetFillColor(0);
  title->SetTextAlign(22);
  title->Draw("same");
  //Draw vertical lines at specified X-coordinate of stations (18 stations)
  double stations_x[18] = {-1518.320, -1344.320, -1170.320, -996.320, -822.320, -648.320, -474.320, -300.320, -126.320, 47.680, 221.680, 395.680, 569.680, 743.680, 917.680, 1091.680, 1265.680, 1439.680};
  TLine *line[18];
  for(int j=0; j<18; j++){
    double x = stations_x[j];
    line[j] = new TLine(x, -8*M_PI, x, 8*M_PI);
    line[j]->SetLineStyle(2);  // Dashed line style
    line[j]->SetLineWidth(1);
    line[j]->SetLineColor(kBlack);
    line[j]->Draw("same");
  }
  // Add text for the double value
  TLatex dphidx;
  dphidx.SetTextSize(0.04);
  dphidx.SetTextAlign(13);
  dphidx.DrawLatexNDC(0.15, 0.85, Form("#alpha : #beta : #chi^{2}/ndf: %f : %f : %f", fitter.dydx(), fitter.y0(), fitter.chi2Dof()));
   canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/2PiAmbiguityPhiVsZ_forSegment/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));
   //canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/ce/PhiVsZ/2PiAmbiguityPhiVsZ_forSegment/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));
  //delete
  delete canvas;
  delete gr;
  for(int j=0; j<18; j++) delete line[j];
}
void PhiZSeedFinder::plot_2PiAmbiguityPhiVsZ_forSegment_mod(int tc, int isegment){
    std::cout<<"---------------------------------"<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    std::cout<<" 2PiAmbiguityPhiVsZ_forSegment_mod "<<isegment<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    _circleFitter.clear();
    //Step1: fit circle of i-th segment with fixed weight
    for(size_t i=0; i<_tcHits.size(); i++){
      if(_tcHits[i].used == false) continue;
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 0.1;//tentative value
      _circleFitter.addPoint(x, y, wP);
    }
    double xC = _circleFitter.x0();
    double yC = _circleFitter.y0();
    double rC = _circleFitter.radius();
    //Step2: fit circle of i-th segment with correct weight
    _circleFitter.clear();
    for(size_t i=0; i<_tcHits.size(); i++){
      if(_tcHits[i].used == false) continue;
      computeCircleError2(i, xC, yC, rC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
    }
    //Step3: fit circle of i-th segment with correct weight
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    _circleFitter.clear();
    for(size_t i=0; i<_tcHits.size(); i++){
      if(_tcHits[i].used == false) continue;
      computeCircleError2(i, xC, yC, rC);
      //computeHelixPhi(i, xC, yC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
    }
std::cout << "Sorted _tcHits by z:\n";
for (size_t i = 0; i < _tcHits.size(); ++i) {
    std::cout << "Index " << i
              << ": hitIndice = " << _tcHits[i].hitIndice
              << ": z = " << _tcHits[i].z
              << ", x = " << _tcHits[i].x
              << ", y = " << _tcHits[i].y
              << ", phi = " << _tcHits[i].phi
              << ", helixphi = " << _tcHits[i].helixPhi
              << ", used/no = " << _tcHits[i].used
              << ", station = " << _tcHits[i].station
              << std::endl;
}
    //Get the slope value from the 1st segment
    //continue collecting hits until station gap is > 1
    int Reference_HitIndex = 0;// hit should be located in the upstream tracker
    int Reference_HitIndice = -1;// hit should be located in the upstream tracker
    int Reference_station = 0;
    int LastStation = 0;
    double Reference_Phi = -999.9;
    std::cout<<"_tcHits.size() = "<<_tcHits.size()<<std::endl;
    for(size_t i=0; i<_tcHits.size(); i++){
      std::cout<<"hitIndice = "<<_tcHits[i].hitIndice<<std::endl;
      std::cout<<"i = "<<i<<std::endl;
      if(_tcHits[i].used == false) continue;
      Reference_HitIndex = i;
      Reference_HitIndice = _tcHits[i].hitIndice;
      Reference_station = _tcHits[i].station;
      computeHelixPhi(i, xC, yC);
      Reference_Phi = _tcHits[i].helixPhi;
      _tcHits[i].ambigPhi = _tcHits[i].helixPhi;
      std::cout<<"_tcHits[i].helixPhi = "<<_tcHits[i].helixPhi<<std::endl;
      std::cout<<"_tcHits[i].ambigPhi = "<<_tcHits[i].ambigPhi<<std::endl;
      std::cout<<"station_z = "<<_tcHits[i].z<<std::endl;
      std::cout<<"station_Reference = "<<_tcHits[i].station<<std::endl;
      std::cout<<"Reference_Index = "<<Reference_HitIndex<<std::endl;
      std::cout<<"Reference_Indice = "<<Reference_HitIndice<<std::endl;
      std::cout<<"Reference_Phi = "<<Reference_Phi<<std::endl;
      break;
    }
    _lineFitter.clear();
    // Prepare graph with errors
    TGraphErrors* gr = new TGraphErrors();
    int graphIndex = 0;
    std::cout<<"========================================"<<std::endl;
    for(size_t i=0; i<_tcHits.size(); i++){
      std::cout<<"i = "<<i<<std::endl;
      double z = _tcHits[i].z;
      computeHelixPhi(i, xC, yC);
      double phiWeight = 1.0 / (_tcHits[i].helixPhiError2);
      //only add reference hit to the fiiting function
      if(Reference_HitIndex == (int)i and _tcHits[i].hitIndice == Reference_HitIndice) {
        _lineFitter.addPoint(z, _tcHits[i].ambigPhi, phiWeight);
      }
      // get other hits
      if(Reference_HitIndex == (int)i) continue;
      if(_tcHits[i].used == false) continue;
      if(_tcHits[i].hitIndice == Reference_HitIndice) continue;
      //coorect helixphi and consider 2pi boundary
      float deltaPhi = _tcHits[i].helixPhi - Reference_Phi;
      std::cout<<"hitIndice = "<<_tcHits[i].hitIndice<<std::endl;
      std::cout<<"Helixphi = "<<_tcHits[i].helixPhi<<std::endl;
      std::cout<<"Z = "<<z<<std::endl;
      std::cout<<"station = "<<_tcHits[i].station<<std::endl;
      std::cout<<"phi = "<<_tcHits[i].helixPhi<<std::endl;
      std::cout<<"deltaPhi = "<<deltaPhi<<std::endl;
      // If it turns more than pi then, consinder the 2pi boundary
      int turns = 0;
      if (deltaPhi > M_PI) turns--;
      if (deltaPhi < -M_PI) turns++;
      double phi = _tcHits[i].helixPhi + turns * 2 * M_PI;
      _tcHits[i].ambigPhi = phi;
      std::cout<<"_tcHits[i].ambigPhi = "<<phi<<std::endl;
      // quality cut for the 1st segment
      if(abs(_tcHits[i].station - Reference_station) >= 3 and _lineFitter.qn() > 5) continue;
      if(_lineFitter.qn() < 2) _lineFitter.addPoint(z, _tcHits[i].ambigPhi, phiWeight);
      else {
        if(turns != 0){
          // corss check
          double lineSlope = _lineFitter.dydx();
          double lineIntercept = _lineFitter.y0();
          // Predict phi from the line
          double predictedPhi = lineSlope * z + lineIntercept;
          // Compute the difference between prediction and actual
          double diffPhi[2] = {0.0};
          diffPhi[0] = predictedPhi - _tcHits[i].ambigPhi;
          diffPhi[1] = predictedPhi - _tcHits[i].helixPhi;
          // choose the nearest assumption
          if(abs(diffPhi[1]) < abs(diffPhi[0])) turns = 0;
          phi = _tcHits[i].helixPhi + turns * 2 * M_PI;
          _tcHits[i].ambigPhi = phi;
          // Round delta/2π to nearest integer for wrapping correction
          _lineFitter.addPoint(z, _tcHits[i].ambigPhi, phiWeight);
        }
        else {
          _lineFitter.addPoint(z, _tcHits[i].ambigPhi, phiWeight);
        }
      }
      std::cout<<"=======================================0= "<<std::endl;
      std::cout<<"addd "<<std::endl;
      std::cout<<"=======================================0= "<<std::endl;
      std::cout<<"i = "<<i<<std::endl;
      //_lineFitter.addPoint(z, _tcHits[i].ambigPhi, phiWeight);
      LastStation = _tcHits[i].station;
      // Set point with error
      //gr->SetPoint(graphIndex, z, phi);
      //gr->SetPointError(graphIndex, 0.0, std::sqrt(_tcHits[i].helixPhiError2));
      //graphIndex++;
    }
      std::cout<<"=======================================0= "<<std::endl;
      std::cout<<"Fini add "<<std::endl;
      std::cout<<"=======================================0= "<<std::endl;
      std::cout<<"LastStation = "<<LastStation<<std::endl;
    //Final
    double lineSlope = _lineFitter.dydx();
    double lineIntercept = _lineFitter.y0();
    std::cout<<"lineSlope = "<<lineSlope<<std::endl;
    std::cout<<"lineIntercept = "<<lineIntercept<<std::endl;
    std::cout<<"========================================== "<<std::endl;
    std::cout<<"00======================================== "<<std::endl;
    std::cout<<"========================================= "<<std::endl;
    _lineFitter.clear();
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if(_tcHits[i].used == false) continue;
      std::cout<<"No  = "<<graphIndex<<std::endl;
      double z = _tcHits[i].z;
      double phiWeight = 1.0 / (_tcHits[i].helixPhiError2);
      double phi = _tcHits[i].ambigPhi;
      double deltaPhi = lineSlope * z + lineIntercept - phi;
      std::cout<<"Z = "<<z<<std::endl;
      std::cout<<"phi = "<<phi<<std::endl;
      std::cout<<"deltaPhi = "<<deltaPhi<<std::endl;
      std::cout<<"deltaPhi / (2 * M_PI) = "<<deltaPhi/ (2 * M_PI)<<std::endl;
      //_tcHits[i].helixPhiCorrection = std::floor(deltaPhi / (2 * M_PI));
      _tcHits[i].helixPhiCorrection = std::round(deltaPhi / (2 * M_PI));
      std::cout<<"helixPhiCorrection  = "<<_tcHits[i].helixPhiCorrection<<std::endl;
      phi = phi + _tcHits[i].helixPhiCorrection * 2 * M_PI;
      std::cout<<"phi  = "<<phi<<std::endl;
      //std::cout<<"No. : z/deltaPhi/helixPhiCorrection/helixphi/phi = "<<graphIndex<<" :  "<<z<<"/"<<deltaPhi<<"/"<<_tcHits[i].helixPhiCorrection<<"/"<<_tcHits[i].helixPhi<<"/"<<phi<<std::endl;
      _lineFitter.addPoint(z, phi, phiWeight);
      // Set point with error
      gr->SetPoint(graphIndex, z, phi);
      gr->SetPointError(graphIndex, 0.0, std::sqrt(_tcHits[i].helixPhiError2));
      graphIndex++;
      //use updated slope value
      if(LastStation < _tcHits[i].station){
        lineSlope = _lineFitter.dydx();
        lineIntercept = _lineFitter.y0();
      }
   }
    /*double lineSlope = _lineFitter.dydx();
    double lineIntercept = _lineFitter.y0();
    int graphIndex = 0;
    _lineFitter.clear();
    std::cout<<"_tcHits.size() = "<<_tcHits.size()<<std::endl;
    std::cout<<"lineSlope = "<<lineSlope<<std::endl;
    std::cout<<"lineIntercept = "<<lineIntercept<<std::endl;
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if(_tcHits[i].used == false) continue;
      std::cout<<"No  = "<<graphIndex<<std::endl;
      double z = _tcHits[i].z;
      double phiWeight = 1.0 / (_tcHits[i].helixPhiError2);
      double phi = _tcHits[i].helixPhi;
      double deltaPhi = lineSlope * z + lineIntercept - phi;
      // If it turns more than pi then, consinder the 2pi boundary
      int turns = 0;
      if (deltaPhi > M_PI) turns--;
      if (deltaPhi < -M_PI) turns++;
      std::cout<<"Z = "<<z<<std::endl;
      std::cout<<"phi = "<<phi<<std::endl;
      std::cout<<"deltaPhi = "<<deltaPhi<<std::endl;
      std::cout<<"deltaPhi / (2 * M_PI) = "<<deltaPhi/ (2 * M_PI)<<std::endl;
      _tcHits[i].helixPhiCorrection = std::floor(deltaPhi / (2 * M_PI));
      std::cout<<"helixPhiCorrection  = "<<_tcHits[i].helixPhiCorrection<<std::endl;
      phi = phi + _tcHits[i].helixPhiCorrection * 2 * M_PI;
      std::cout<<"phi  = "<<phi<<std::endl;
      //std::cout<<"No. : z/deltaPhi/helixPhiCorrection/helixphi/phi = "<<graphIndex<<" :  "<<z<<"/"<<deltaPhi<<"/"<<_tcHits[i].helixPhiCorrection<<"/"<<_tcHits[i].helixPhi<<"/"<<phi<<std::endl;
      _lineFitter.addPoint(z, phi, phiWeight);
      gr->SetPoint(graphIndex++, z, phi);
   }*/
  //gr->SetMarkerStyle(20);
  //gr->SetMarkerSize(0.8);
  //gr->SetMarkerColor(kRed);
gr->SetTitle("");
gr->SetMarkerStyle(20);        // 'x' style marker. Use 1 for '+' style if preferred
gr->SetMarkerSize(0.4);       // Increase size a bit to make the cross clearer
gr->SetMarkerColor(kRed);
gr->SetLineColor(kRed);       // Set line color for the error bars
// Perform ROOT linear fit on graph with errors
//TF1* fitFunc = new TF1("fitFunc", "pol1", -1600, 1600);
//gr->Fit(fitFunc, "Q"); // Quiet fit
// Extract ROOT fit parameters
//double root_intercept = fitFunc->GetParameter(0);
//double root_slope = fitFunc->GetParameter(1);
//double root_intercept_err = fitFunc->GetParError(0);
//double root_slope_err = fitFunc->GetParError(1);
//double chi2 = fitFunc->GetChisquare();
//int ndf = fitFunc->GetNDF();
  //Draw
  TCanvas *canvas = new TCanvas("canvas", "", 800, 600);
  canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
  double ymax = M_PI*4;
  double ymin = -M_PI*4;
  gr->Draw("AP");
  gr->GetXaxis()->SetTitle("Z [mm]");
  gr->GetYaxis()->SetTitle("Phi [rad]");
  gr->GetXaxis()->SetLimits(-1600, 1600);
  //gr->GetYaxis()->SetRangeUser(-8*M_PI, 8*M_PI);
  //gr->GetYaxis()->SetRangeUser(0, 4*M_PI);
  gr->GetYaxis()->SetRangeUser(ymin, ymax);
  // Add title at the top
  TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
  std::stringstream eventStringStream;
  eventStringStream << "run: " << run << " subRun: " << subrun << " event: " << eventNumber;
  title->AddText(Form("Phi vs. Z (Run-subRun-Event, TC, #cand) = (%d-%d-%d, %d, #%d) (pbar1b0)", run, subrun, eventNumber, tc, isegment));
  title->SetFillColor(0);
  title->SetTextAlign(22);
  title->Draw("same");
  //Draw vertical lines at specified X-coordinate of stations (18 stations)
  double stations_x[18] = {-1518.320, -1344.320, -1170.320, -996.320, -822.320, -648.320, -474.320, -300.320, -126.320, 47.680, 221.680, 395.680, 569.680, 743.680, 917.680, 1091.680, 1265.680, 1439.680};
  TLine *line[18];
  for(int j=0; j<18; j++){
    double x = stations_x[j];
    line[j] = new TLine(x, ymin, x, ymax);
    line[j]->SetLineStyle(2);  // Dashed line style
    line[j]->SetLineWidth(1);
    line[j]->SetLineColor(kBlack);
    line[j]->Draw("same");
  }
  // Draw the fitted slope line: y = dydx * x + y0 from x = -1600 to 1600
  double x_min = -1600;
  double x_max = 1600;
  double y_min = _lineFitter.dydx() * x_min + _lineFitter.y0();
  double y_max = _lineFitter.dydx() * x_max + _lineFitter.y0();
  TLine *fitLine = new TLine(x_min, y_min, x_max, y_max);
  fitLine->SetLineColor(kBlack);
  fitLine->SetLineStyle(2);  // Dashed line
  fitLine->SetLineWidth(2);
  fitLine->Draw("same");
// Draw manual fit parameters text
TLatex paramsText;
paramsText.SetTextSize(0.04);
paramsText.SetTextAlign(13);
double textX = 0.15;
double textY = 0.85;
paramsText.DrawLatexNDC(textX, textY,        Form("LsqSums2 fit dydx: %f", _lineFitter.dydx()));
paramsText.DrawLatexNDC(textX, textY - 0.05, Form("LsqSums2 fit y0: %f", _lineFitter.y0()));
paramsText.DrawLatexNDC(textX, textY - 0.10, Form("LsqSums2 fit #chi^{2}/ndf: %f", _lineFitter.chi2Dof()));
// Draw ROOT fit parameters text below manual
//paramsText.DrawLatexNDC(textX, textY - 0.18, Form("ROOT fit slope: %.5f #pm %.5f", root_slope, root_slope_err));
//paramsText.DrawLatexNDC(textX, textY - 0.23, Form("ROOT fit intercept: %.5f #pm %.5f", root_intercept, root_intercept_err));
//paramsText.DrawLatexNDC(textX, textY - 0.28, Form("ROOT fit #chi^{2}/ndf: %.2f / %d", chi2, ndf));
   canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/2PiAmbiguityPhiVsZ_forSegment/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));
   //canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/ce/PhiVsZ/2PiAmbiguityPhiVsZ_forSegment/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));
  //delete
  delete canvas;
delete gr;
for (int j = 0; j < 18; j++) delete line[j];
delete fitLine;
//delete fitFunc;
}
/*void PhiZSeedFinder::plot_2PiAmbiguityPhiVsZ_forSegment_mod(int tc, int isegment){
    std::cout<<"---------------------------------"<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    std::cout<<" 2PiAmbiguityPhiVsZ_forSegment_mod "<<isegment<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    std::cout<<"---------------------------------"<<std::endl;
    _circleFitter.clear();
    //Step1: fit circle of i-th segment with fixed weight
    for(size_t i=0; i<_tcHits.size(); i++){
      if(_tcHits[i].used == false) continue;
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 0.1;//tentative value
      _circleFitter.addPoint(x, y, wP);
    }
    double xC = _circleFitter.x0();
    double yC = _circleFitter.y0();
    double rC = _circleFitter.radius();
    //Step2: fit circle of i-th segment with correct weight
    _circleFitter.clear();
    for(size_t i=0; i<_tcHits.size(); i++){
      if(_tcHits[i].used == false) continue;
      computeCircleError2(i, xC, yC, rC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
    }
    //Step3: fit circle of i-th segment with correct weight
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    _circleFitter.clear();
    for(size_t i=0; i<_tcHits.size(); i++){
      if(_tcHits[i].used == false) continue;
      computeCircleError2(i, xC, yC, rC);
      //computeHelixPhi(i, xC, yC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
    }
std::cout << "Sorted _tcHits by z:\n";
for (size_t i = 0; i < _tcHits.size(); ++i) {
    std::cout << "Index " << i
              << ": hitIndice = " << _tcHits[i].hitIndice
              << ": z = " << _tcHits[i].z
              << ", x = " << _tcHits[i].x
              << ", y = " << _tcHits[i].y
              << ", phi = " << _tcHits[i].phi
              << ", helixphi = " << _tcHits[i].helixPhi
              << ", used/no = " << _tcHits[i].used
              << ", station = " << _tcHits[i].station
              << std::endl;
}
    //Get the slope value from the 1st segment
    //continue collecting hits until station gap is > 1
    int Reference_HitIndex = 0;// hit should be located in the upstream tracker
    int Reference_HitIndice = -1;// hit should be located in the upstream tracker
    int Reference_station = 0;
    int LastStation = 0;
    double Reference_Phi = -999.9;
    std::cout<<"_tcHits.size() = "<<_tcHits.size()<<std::endl;
    for(size_t i=0; i<_tcHits.size(); i++){
      std::cout<<"hitIndice = "<<_tcHits[i].hitIndice<<std::endl;
      std::cout<<"i = "<<i<<std::endl;
      if(_tcHits[i].used == false) continue;
      Reference_HitIndex = i;
      Reference_HitIndice = _tcHits[i].hitIndice;
      Reference_station = _tcHits[i].station;
      computeHelixPhi(i, xC, yC);
      Reference_Phi = _tcHits[i].helixPhi;
      _tcHits[i].ambigPhi = _tcHits[i].helixPhi;
      std::cout<<"_tcHits[i].helixPhi = "<<_tcHits[i].helixPhi<<std::endl;
      std::cout<<"_tcHits[i].ambigPhi = "<<_tcHits[i].ambigPhi<<std::endl;
      std::cout<<"station_z = "<<_tcHits[i].z<<std::endl;
      std::cout<<"station_Reference = "<<_tcHits[i].station<<std::endl;
      std::cout<<"Reference_Index = "<<Reference_HitIndex<<std::endl;
      std::cout<<"Reference_Indice = "<<Reference_HitIndice<<std::endl;
      std::cout<<"Reference_Phi = "<<Reference_Phi<<std::endl;
      break;
    }
    _lineFitter.clear();
    // Prepare graph with errors
    TGraphErrors* gr = new TGraphErrors();
    int graphIndex = 0;
    std::cout<<"========================================"<<std::endl;
    for(size_t i=0; i<_tcHits.size(); i++){
      if(Reference_HitIndex == (int)i) continue;
      if(_tcHits[i].used == false) continue;
      if(_tcHits[i].hitIndice == Reference_HitIndice) continue;
      std::cout<<"i = "<<i<<std::endl;
      double z = _tcHits[i].z;
      computeHelixPhi(i, xC, yC);
      double phiWeight = 1.0 / (_tcHits[i].helixPhiError2);
      float deltaPhi = _tcHits[i].helixPhi - Reference_Phi;
      std::cout<<"hitIndice = "<<_tcHits[i].hitIndice<<std::endl;
      std::cout<<"Helixphi = "<<_tcHits[i].helixPhi<<std::endl;
      std::cout<<"Z = "<<z<<std::endl;
      std::cout<<"station = "<<_tcHits[i].station<<std::endl;
      std::cout<<"phi = "<<_tcHits[i].helixPhi<<std::endl;
      std::cout<<"deltaPhi = "<<deltaPhi<<std::endl;
      // If it turns more than pi then, consinder the 2pi boundary
      int turns = 0;
      if (deltaPhi > M_PI) turns--;
      if (deltaPhi < -M_PI) turns++;
      double phi = _tcHits[i].helixPhi + turns * 2 * M_PI;
      _tcHits[i].ambigPhi = phi;
      std::cout<<"_tcHits[i].ambigPhi = "<<phi<<std::endl;
      // quality cut for the 1st segment
      if(abs(_tcHits[i].station - Reference_station) >= 3 and _lineFitter.qn() > 5) continue;
      std::cout<<"============= "<<std::endl;
      std::cout<<"=======================================0= "<<std::endl;
      std::cout<<"addd "<<std::endl;
      std::cout<<"i = "<<i<<std::endl;
      _lineFitter.addPoint(z, _tcHits[i].ambigPhi, phiWeight);
      LastStation = _tcHits[i].station;
      // Set point with error
      //gr->SetPoint(graphIndex, z, phi);
      //gr->SetPointError(graphIndex, 0, phiError);
      //graphIndex++;
    }
    //Final
    double lineSlope = _lineFitter.dydx();
    double lineIntercept = _lineFitter.y0();
    std::cout<<"lineSlope = "<<lineSlope<<std::endl;
    std::cout<<"lineIntercept = "<<lineIntercept<<std::endl;
    std::cout<<"========================================== "<<std::endl;
    std::cout<<"00======================================== "<<std::endl;
    std::cout<<"========================================= "<<std::endl;
    _lineFitter.clear();
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if(_tcHits[i].used == false) continue;
      std::cout<<"No  = "<<graphIndex<<std::endl;
      double z = _tcHits[i].z;
      double phiWeight = 1.0 / (_tcHits[i].helixPhiError2);
      double phi = _tcHits[i].ambigPhi;
      double deltaPhi = lineSlope * z + lineIntercept - phi;
      std::cout<<"Z = "<<z<<std::endl;
      std::cout<<"phi = "<<phi<<std::endl;
      std::cout<<"deltaPhi = "<<deltaPhi<<std::endl;
      std::cout<<"deltaPhi / (2 * M_PI) = "<<deltaPhi/ (2 * M_PI)<<std::endl;
      //_tcHits[i].helixPhiCorrection = std::floor(deltaPhi / (2 * M_PI));
      _tcHits[i].helixPhiCorrection = std::round(deltaPhi / (2 * M_PI));
      std::cout<<"helixPhiCorrection  = "<<_tcHits[i].helixPhiCorrection<<std::endl;
      phi = phi + _tcHits[i].helixPhiCorrection * 2 * M_PI;
      std::cout<<"phi  = "<<phi<<std::endl;
      //std::cout<<"No. : z/deltaPhi/helixPhiCorrection/helixphi/phi = "<<graphIndex<<" :  "<<z<<"/"<<deltaPhi<<"/"<<_tcHits[i].helixPhiCorrection<<"/"<<_tcHits[i].helixPhi<<"/"<<phi<<std::endl;
      _lineFitter.addPoint(z, phi, phiWeight);
      // Set point with error
      gr->SetPoint(graphIndex, z, phi);
      gr->SetPointError(graphIndex, 0.0, std::sqrt(_tcHits[i].helixPhiError2));
      graphIndex++;
      //use updated slope value
      if(LastStation < _tcHits[i].station){
        lineSlope = _lineFitter.dydx();
        lineIntercept = _lineFitter.y0();
      }
   }
  //gr->SetMarkerStyle(20);
  //gr->SetMarkerSize(0.8);
  //gr->SetMarkerColor(kRed);
gr->SetTitle("");
gr->SetMarkerStyle(20);        // 'x' style marker. Use 1 for '+' style if preferred
gr->SetMarkerSize(0.4);       // Increase size a bit to make the cross clearer
gr->SetMarkerColor(kRed);
gr->SetLineColor(kRed);       // Set line color for the error bars
// Perform ROOT linear fit on graph with errors
//TF1* fitFunc = new TF1("fitFunc", "pol1", -1600, 1600);
//gr->Fit(fitFunc, "Q"); // Quiet fit
// Extract ROOT fit parameters
//double root_intercept = fitFunc->GetParameter(0);
//double root_slope = fitFunc->GetParameter(1);
//double root_intercept_err = fitFunc->GetParError(0);
//double root_slope_err = fitFunc->GetParError(1);
//double chi2 = fitFunc->GetChisquare();
//int ndf = fitFunc->GetNDF();
  //Draw
  TCanvas *canvas = new TCanvas("canvas", "", 800, 600);
  canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
  gr->Draw("AP");
  gr->GetXaxis()->SetTitle("Z [mm]");
  gr->GetYaxis()->SetTitle("Phi [rad]");
  gr->GetXaxis()->SetLimits(-1600, 1600);
  //gr->GetYaxis()->SetRangeUser(-8*M_PI, 8*M_PI);
  //gr->GetYaxis()->SetRangeUser(0, 4*M_PI);
  gr->GetYaxis()->SetRangeUser(-M_PI*3, M_PI*3);
  // Add title at the top
  TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
  std::stringstream eventStringStream;
  eventStringStream << "run: " << run << " subRun: " << subrun << " event: " << eventNumber;
  title->AddText(Form("Phi vs. Z (Run-subRun-Event, TC, #cand) = (%d-%d-%d, %d, #%d) (pbar1b0)", run, subrun, eventNumber, tc, isegment));
  title->SetFillColor(0);
  title->SetTextAlign(22);
  title->Draw("same");
  //Draw vertical lines at specified X-coordinate of stations (18 stations)
  double stations_x[18] = {-1518.320, -1344.320, -1170.320, -996.320, -822.320, -648.320, -474.320, -300.320, -126.320, 47.680, 221.680, 395.680, 569.680, 743.680, 917.680, 1091.680, 1265.680, 1439.680};
  TLine *line[18];
  for(int j=0; j<18; j++){
    double x = stations_x[j];
    line[j] = new TLine(x, -M_PI*3, x, M_PI*3);
    line[j]->SetLineStyle(2);  // Dashed line style
    line[j]->SetLineWidth(1);
    line[j]->SetLineColor(kBlack);
    line[j]->Draw("same");
  }
  // Draw the fitted slope line: y = dydx * x + y0 from x = -1600 to 1600
  double x_min = -1600;
  double x_max = 1600;
  double y_min = _lineFitter.dydx() * x_min + _lineFitter.y0();
  double y_max = _lineFitter.dydx() * x_max + _lineFitter.y0();
  TLine *fitLine = new TLine(x_min, y_min, x_max, y_max);
  fitLine->SetLineColor(kBlack);
  fitLine->SetLineStyle(2);  // Dashed line
  fitLine->SetLineWidth(2);
  fitLine->Draw("same");
// Draw manual fit parameters text
TLatex paramsText;
paramsText.SetTextSize(0.04);
paramsText.SetTextAlign(13);
double textX = 0.15;
double textY = 0.85;
paramsText.DrawLatexNDC(textX, textY,        Form("LsqSums2 fit dydx: %f", _lineFitter.dydx()));
paramsText.DrawLatexNDC(textX, textY - 0.05, Form("LsqSums2 fit y0: %f", _lineFitter.y0()));
paramsText.DrawLatexNDC(textX, textY - 0.10, Form("LsqSums2 fit #chi^{2}/ndf: %f", _lineFitter.chi2Dof()));
// Draw ROOT fit parameters text below manual
//paramsText.DrawLatexNDC(textX, textY - 0.18, Form("ROOT fit slope: %.5f #pm %.5f", root_slope, root_slope_err));
//paramsText.DrawLatexNDC(textX, textY - 0.23, Form("ROOT fit intercept: %.5f #pm %.5f", root_intercept, root_intercept_err));
//paramsText.DrawLatexNDC(textX, textY - 0.28, Form("ROOT fit #chi^{2}/ndf: %.2f / %d", chi2, ndf));
   canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/2PiAmbiguityPhiVsZ_forSegment/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));
   //canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/ce/PhiVsZ/2PiAmbiguityPhiVsZ_forSegment/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));
  //delete
  delete canvas;
delete gr;
for (int j = 0; j < 18; j++) delete line[j];
delete fitLine;
//delete fitFunc;
}
*/
//-----------------------------------------------------------------------------
void PhiZSeedFinder::plot_RVsZ_forSegment(int tc, int isegment){
  const int n = (int)_segmentHits.at(isegment).size();
  TGraph *gr = new TGraph(n);
  gr->SetTitle("");
  gr->SetMarkerStyle(1);
  double r, z;
  std::vector<int> marker_color;
  std::vector<int> marker_style;
  std::vector<double> marker_size;
    std::cout<<"Fill "<<std::endl;
  for(int i=0; i<n; i++) {
   double z = _segmentHits.at(isegment).at(i).z;
   double r = sqrt(_segmentHits.at(isegment).at(i).x * _segmentHits.at(isegment).at(i).x + _segmentHits.at(isegment).at(i).y * _segmentHits.at(isegment).at(i).y);
   int index = _segmentHits.at(isegment).at(i).hitIndice;
   int alreadyfill = 0;
   for (size_t j = 0; j < _data.tccol->at(tc)._strawHitIdxs.size(); j++) {
    int hitIndice = _data.tccol->at(tc)._strawHitIdxs[j];
    if(index != hitIndice) continue;
    if(alreadyfill == 1) continue;
    std::vector<StrawDigiIndex> shids;
    _data.chcol->fillStrawDigiIndices(hitIndice, shids);
    //z = _data.chcol->at(hitIndice).pos().z();
    //phi = _data.chcol->at(hitIndice).pos().phi();
    std::cout<<"z/r = "<<z<<"/"<<r<<std::endl;
    for (size_t k = 0; k < shids.size(); k++) {
      const mu2e::SimParticle* _simParticle;
      _simParticle = _mcUtils->getSimParticle(_event, shids[k]);
      //int SimID = _mcUtils->strawHitSimId(_event, shids[j]);
      int PdgID = _simParticle->pdgId();
      int   style(0), color(0);
      double size(0.);
      if (PdgID ==  11) {style = 20; size = 0.8; color = kRed;}
      else if (PdgID ==  -11) {style = 24; size = 0.8; color = kBlue;}
      else if (PdgID ==  13) {style = 20; size = 0.8; color = kGreen+2;}
      else if (PdgID ==  -13) {style = 20; size = 0.8; color = kGreen-2;}
      else if (PdgID ==  2212) {style = 20; size = 0.8; color = kBlue+2;}
      else if (PdgID ==  -211) {style = 20; size = 0.8; color = kPink+2;}
      else if (PdgID ==  211) {style = 20; size = 0.8; color = kPink-2;}
      else {style = 20; size = 0.8; color = kMagenta;}
      marker_color.push_back(color);
      marker_style.push_back(style);
      marker_size.push_back(size);
      alreadyfill = 1;
      break;
    }
   }
    gr->SetPoint(i, z, r);
   }
  //Draw
  TCanvas *canvas = new TCanvas("canvas", "", 800, 600);
  canvas->SetMargin(0.1, 0.1, 0.1, 0.1);
  gr->Draw("AP");
  TMarker *m;
  for (int i = 0; i < n; i++) {
    gr->GetPoint(i, z, r);
    std::cout<<"z/r = "<<z<<"/"<<r<<std::endl;
    m = new TMarker(z, r, 20);
    //m->SetMarkerColor(i + 1);// Setting marker color with different color for each point
    m->SetMarkerStyle(marker_style[i]);
    m->SetMarkerSize(marker_size[i]);
    m->SetMarkerColor(marker_color[i]);// Setting marker color with different color for each point
    m->Draw();// Draw marker with different color
  }
  gr->GetXaxis()->SetTitle("Z [mm]");
  gr->GetYaxis()->SetTitle("R [mm]");
  gr->GetXaxis()->SetLimits(-1600, 1600);
  //gr->GetYaxis()->SetLimits(350, 750);
  gr->GetYaxis()->SetRangeUser(350, 750);
  // Add title at the top
  TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
  std::stringstream eventStringStream;
  eventStringStream << "run: " << run << " subRun: " << subrun << " event: " << eventNumber;
  title->AddText(Form("R vs. Z (Run-subRun-Event, TC, #cand) = (%d-%d-%d, %d, #%d) (pbar1b0)", run, subrun, eventNumber, tc, isegment));
  title->SetFillColor(0);
  title->SetTextAlign(22);
  title->Draw("same");
  //Draw vertical lines at specified X-coordinate of stations (18 stations)
  double stations_x[18] = {-1518.320, -1344.320, -1170.320, -996.320, -822.320, -648.320, -474.320, -300.320, -126.320, 47.680, 221.680, 395.680, 569.680, 743.680, 917.680, 1091.680, 1265.680, 1439.680};
  TLine *line[18];
  for(int j=0; j<18; j++){
    double x = stations_x[j];
    line[j] = new TLine(x, 350, x, 750);
    line[j]->SetLineStyle(2);  // Dashed line style
    line[j]->SetLineWidth(1);
    line[j]->SetLineColor(kBlack);
    line[j]->Draw("same");
  }
   canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/RVsZ/segment/pbar_RVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));
   //canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/ce/RVsZ/segment/pbar_RVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));
  //delete
  delete canvas;
  delete gr;
  for(int j=0; j<18; j++) delete line[j];
}
//-----------------------------------------------------------------------------
void PhiZSeedFinder::plot_XVsY(int tc, int isegment, const char* filename, double& xC, double& yC, double& rC){
    std::cout<<"plot_XVsY"<<std::endl;
    //TGraph for Helix (can be removed in the future)
    TMultiGraph* graph  = new TMultiGraph();
    TGraph *gr1 = new TGraph(); gr1->SetMarkerStyle(20);
    TGraph *gr2 = new TGraph(); gr2->SetMarkerStyle(20);
    gr1->SetMarkerColor(kRed);
    gr1->SetMarkerSize(0.5);
    gr2->SetMarkerColor(kBlack);
    gr2->SetMarkerSize(0.5);
    size_t nComboHitsInSegment = _tcHits.size();
    int index[2] = {0};
    for(size_t j=0; j<nComboHitsInSegment; j++){
      if(_tcHits[j].used == false) continue;
      //if(j != 16) continue; //bokuno
      double x = _tcHits.at(j).x;
      double y = _tcHits.at(j).y;
      if(j >9999)gr1->SetPoint(index[0]++, x, y);
      gr2->SetPoint(index[1]++, x, y);
              std::cout
              << " | index = " << j
              << " | hitIndice = " << _tcHits[j].hitIndice
              << " | wireErr = " << _data.chcol->at(_tcHits[j].hitIndice).wireRes()
              << " | x = " << _tcHits[j].x
              << " | y = " << _tcHits[j].y
              << " | z = " << _tcHits[j].z
              << " | phi = " << _tcHits[j].phi
              << " | helixPhi (0) = " << _tcHits[j].helixPhi
              << " | helixPhiError2 = " << _tcHits[j].helixPhiError2
              << " | sqrt(helixPhiError2) = " << sqrt(_tcHits[j].helixPhiError2)
              << std::endl;
    }
    /*for(size_t j=0; j<nComboHitsInSegment; j++){
      if(_tcHits[j].used == false) continue;
      double x = _tcHits.at(j).x;
      double y = _tcHits.at(j).y;
      if(j >9999)gr1->SetPoint(index[0]++, x, y);
      gr2->SetPoint(index[1]++, x, y);
    }*/
    if(index[0] >= 1) graph->Add(gr1,"AP");
    if(index[1] >= 1) graph->Add(gr2,"AP");
    TCanvas *canvas4 = new TCanvas("canvas4", "My TGraph", 800, 800);
    canvas4->SetMargin(0.15, 0.1, 0.1, 0.1);
    graph->GetYaxis()->SetRangeUser(-900, 900);
    graph->GetXaxis()->SetLimits(-900, 900);
    graph->Draw("AP");
    graph->GetXaxis()->SetTitle("X [mm]");
    graph->GetYaxis()->SetTitle("Y [mm]");
    // Draw the first circle with radius 680
    TEllipse *circle_1 = new TEllipse(0, 0, 680);
    circle_1->SetLineColor(kGreen); // Set the line color to green
    circle_1->SetFillStyle(0);       // Set fill style to transparent
    circle_1->Draw("same");
    // Draw the second circle with radius 360
    TEllipse *circle_2 = new TEllipse(0, 0, 360);
    circle_2->SetLineColor(kGreen);  // Set the line color to green
    circle_2->SetFillStyle(0);       // Set fill style to transparent
    circle_2->Draw("same");          // Draw on the same canvas4
    // Draw the Helic: circle fit result
    TEllipse *circle_3 = new TEllipse(_circleFitter.x0(), _circleFitter.y0(), _circleFitter.radius());
    circle_3->SetLineColor(kBlack);  // Set the line color to green
    circle_3->SetFillStyle(0);       // Set fill style to transparent
    circle_3->Draw("same");          // Draw on the same canvas4
    // Draw the Helic: used for weight
    //TEllipse *circle4 = new TEllipse(xC, yC, rC);
    //circle4->SetLineColor(kBlue);  // Set the line color to green
    //circle4->SetFillStyle(0);       // Set fill style to transparent
    //circle4->Draw("same");          // Draw on the same canvas4
    // Draw the MC Helics
    //TEllipse *circle_5 = new TEllipse(_mcX0, _mcY0, _mcRadius);
    //circle_5->SetLineColor(kRed);  // Set the line color to green
    //circle_5->SetLineStyle(7);  // Set the line color to green
    //circle_5->SetLineWidth(3);  // Set the line color to green
    //circle_5->SetFillStyle(0);       // Set fill style to transparent
    //circle_5->Draw("same");          // Draw on the same canvas4
    std::cout<<"amanaki "<<std::endl;
    int num = 0;
    if(_mcParticleInTC == 1){
    for(size_t k=0; k<_simIDsPerTC.size(); k++) {
        for(size_t l = 0; l < _simIDsPerTC[k].size(); l++) {
          int TcIndex = _simIDsPerTC.at(k).at(l).tcIndex;
            if(TcIndex == tc) num++;
        }
     }
    }
    std::cout<<"amanaki "<<std::endl;
    std::cout<<"num = "<<num<<std::endl;
    TEllipse **circle = nullptr;
    //TEllipse **circle = new TEllipse*[num];  // Dynamically allocate an array of pointers to TEllipse
    if(num != 0) {
    circle = new TEllipse*[num];
    int index = 0;
    for(size_t k=0; k<_simIDsPerTC.size(); k++) {
        for(size_t l = 0; l < _simIDsPerTC[k].size(); l++) {
          int TcIndex = _simIDsPerTC.at(k).at(l).tcIndex;
          if(TcIndex != tc) continue;
            // Create a new TEllipse for each iteration
          circle[index] = new TEllipse(_simIDsPerTC.at(k).at(l).mcX0, _simIDsPerTC.at(k).at(l).mcY0, _simIDsPerTC.at(k).at(l).mcRadius);
          // Set properties for each circle
          circle[index]->SetLineColor(kRed);  // Set the line color to red
          circle[index]->SetLineStyle(7);     // Set the line style
          circle[index]->SetLineWidth(1);     // Set the line width
          circle[index]->SetFillStyle(0);     // Set fill style to transparent
          circle[index]->Draw("same");        // Draw on the same canvas
          index++;
        }
     }
    }
    std::cout<<"gomi1"<<std::endl;
    // Draw a cross line at the center (0, 0) in blue color
    TLine *crossLineX = new TLine(-50, 0, 50, 0);
    TLine *crossLineY = new TLine(0, -50, 0, 50);
    crossLineX->SetLineColor(kBlue);
    crossLineY->SetLineColor(kBlue);
    crossLineX->Draw("same");
    crossLineY->Draw("same");
    // Draw a cross line at the center of Helix in blue color
    TLine *crossHelixLineX = new TLine(-30+xC, yC, 30+xC, yC);
    TLine *crossHelixLineY = new TLine(xC, -30+yC, xC, 30+yC);
    crossHelixLineX->SetLineColor(kBlue);
    crossHelixLineY->SetLineColor(kBlue);
    crossHelixLineX->Draw("same");
    crossHelixLineY->Draw("same");
    // Draw a resolution on wire for each hit
    //TLine* LineW[nComboHitsInSegment];
    std::vector<TLine*> LineW(nComboHitsInSegment, nullptr);
    for(size_t j=0; j<nComboHitsInSegment; j++){
      if(_tcHits[j].used == false) continue;
      //if(j != 16) continue; //bokuno
      int hitIndice = _tcHits[j].hitIndice;
      double fSigW = _data.chcol->at(hitIndice).wireRes();
      //double fSigR  = 2.5;
      double DirX = _data.chcol->at(hitIndice).uDir().x();
      double DirY = _data.chcol->at(hitIndice).uDir().y();
      double LineW_xy[4] = {_tcHits.at(j).x-DirX*fSigW, _tcHits.at(j).y-DirY*fSigW, _tcHits.at(j).x+DirX*fSigW, _tcHits.at(j).y+DirY*fSigW};
      LineW[j] = new TLine(LineW_xy[0], LineW_xy[1], LineW_xy[2], LineW_xy[3]);
      //TLine *LineR = new TLine(_circleFitter.x0(), -30+_circleFitter.y0(), _circleFitter.x0(), 30+_circleFitter.y0());
      LineW[j]->SetLineColor(kBlack);
      LineW[j]->Draw("sames");
    }
    std::cout<<"gomi8"<<std::endl;
    // Add title at the top
    TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
    title->AddText(Form("X vs. Y (Run-subRun-Event, TC, #candidate) = (%d-%d-%d, %d, %d) (pbar1b0)", run, subrun, eventNumber, tc, isegment));
    title->SetFillColor(0);
    title->SetTextAlign(22);
    title->Draw("same");
    TLatex legend;
    legend.SetTextSize(0.03);
    //double w = 1.0/_tcHits[j].circleError2;
    //legend.DrawLatexNDC(0.18, 0.85, Form("%d CHs, #chi^{2}/ndf = %.3f, hit #%d, w = %.6f", (int)_circleFitter.qn(), _circleFitter.chi2DofCircle(), (int)j, w));
    legend.DrawLatexNDC(0.18, 0.85, Form("%d CHs, #chi^{2}/ndf = %.3f", (int)_circleFitter.qn(), _circleFitter.chi2DofCircle()));
    legend.DrawLatexNDC(0.18, 0.82, Form("(xC, yC, rC) = (%.1f, %.1f, %.1f)", _circleFitter.x0(), _circleFitter.y0(), _circleFitter.radius()));
    legend.DrawLatexNDC(0.18, 0.79, Form("%s", filename));
    //TLegend *legend1 = new TLegend(0.18,0.1,0.48,0.2);
    //legend1->AddEntry("circle_4","used for weight correction","l");
    //legend1->AddEntry("circle_3","fitting result","l");
    //legend1->Draw();
    canvas4->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/findHelix/%s/pbar_Helix_%04d-%04d-%04d_TC-%d_cand_%d.pdf", filename, run, subrun, eventNumber, tc, isegment));
    //canvas4->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/ce/findHelix/%s/pbar_Helix_%04d-%04d-%04d_TC-%d_cand_%d.pdf", filename, run, subrun, eventNumber, tc, isegment));
    delete canvas4;
    delete gr1;
    delete gr2;
    delete title;
    delete circle_1;
    delete circle_2;
    delete circle_3;
    //delete circle4;
    //delete circle_5;
    // Delete dynamically allocated memory
    std::cout<<"gomi8"<<std::endl;
    if(num != 0){
      for(int j = 0; j < num; j++) {
        delete circle[j];  // Delete each individual TEllipse object
      }
      delete[] circle;  // Delete the array of pointers
    }
    //delete legend1;
    delete crossLineX;
    delete crossLineY;
    delete crossHelixLineX;
    delete crossHelixLineY;
    //}
    std::cout<<"gomi14"<<std::endl;
}
//-----------------------------------------------------------------------------
void PhiZSeedFinder::plot_XVsY_hit(int tc, int isegment, const char* filename, double& xC, double& yC, double& rC){
    size_t nComboHitsInSegment = _tcHits.size();
    int index[2] = {0};
    for(size_t j=0; j<nComboHitsInSegment; j++){
    //TGraph for Helix (can be removed in the future)
    TMultiGraph* graph  = new TMultiGraph();
    TGraph *gr1 = new TGraph(); gr1->SetMarkerStyle(20);
    TGraph *gr2 = new TGraph(); gr2->SetMarkerStyle(20);
    gr1->SetMarkerColor(kRed);
    gr1->SetMarkerSize(0.5);
    gr2->SetMarkerColor(kBlack);
    gr2->SetMarkerSize(0.5);
      if(_tcHits[j].used == false) continue;
      //if(j != 16) continue; //bokuno
      double x = _tcHits.at(j).x;
      double y = _tcHits.at(j).y;
      if(j >9999)gr1->SetPoint(index[0]++, x, y);
      gr2->SetPoint(index[1]++, x, y);
    if(index[0] >= 1) graph->Add(gr1,"AP");
    if(index[1] >= 1) graph->Add(gr2,"AP");
    TCanvas *canvas4 = new TCanvas("canvas4", "My TGraph", 800, 800);
    canvas4->SetMargin(0.15, 0.1, 0.1, 0.1);
    graph->GetYaxis()->SetRangeUser(-900, 900);
    graph->GetXaxis()->SetLimits(-900, 900);
    graph->Draw("AP");
    graph->GetXaxis()->SetTitle("X [mm]");
    graph->GetYaxis()->SetTitle("Y [mm]");
    // Draw the first circle with radius 680
    TEllipse *circle1 = new TEllipse(0, 0, 680);
    circle1->SetLineColor(kGreen); // Set the line color to green
    circle1->SetFillStyle(0);       // Set fill style to transparent
    circle1->Draw("same");
    // Draw the second circle with radius 360
    TEllipse *circle2 = new TEllipse(0, 0, 360);
    circle2->SetLineColor(kGreen);  // Set the line color to green
    circle2->SetFillStyle(0);       // Set fill style to transparent
    circle2->Draw("same");          // Draw on the same canvas4
    // Draw the Helic: fitting result
    TEllipse *circle3 = new TEllipse(_circleFitter.x0(), _circleFitter.y0(), _circleFitter.radius());
    circle3->SetFillStyle(0);       // Set fill style to transparent
    circle3->SetLineColor(7);  // Set the line color to green
    circle3->Draw("same");          // Draw on the same canvas4
    // Draw the Helic: cirlc used for fitting
    TEllipse *circle4 = new TEllipse(xC, yC, rC);
    circle4->SetFillStyle(0);       // Set fill style to transparent
    circle4->SetLineColor(kBlue);  // Set the line color to green
    circle4->Draw("same");          // Draw on the same canvas4
    // Draw the MC Helics
    TEllipse *circle5 = new TEllipse(_mcX0, _mcY0, _mcRadius);
    circle5->SetLineStyle(7);  // Set the line color to green
    circle5->SetLineColor(kRed);  // Set the line color to green
    circle5->SetFillStyle(0);       // Set fill style to transparent
    circle5->Draw("same");          // Draw on the same canvas4
    // Draw a cross line at the center (0, 0) in blue color
    TLine *crossLineX = new TLine(-50, 0, 50, 0);
    TLine *crossLineY = new TLine(0, -50, 0, 50);
    crossLineX->SetLineColor(kBlue);
    crossLineY->SetLineColor(kBlue);
    crossLineX->Draw("same");
    crossLineY->Draw("same");
    // Draw a cross line at the center of Helix in blue color
    /*TLine *crossHelixLineX = new TLine(-30+_circleFitter.x0(), _circleFitter.y0(), 30+_circleFitter.x0(), _circleFitter.y0());
    TLine *crossHelixLineY = new TLine(_circleFitter.x0(), -30+_circleFitter.y0(), _circleFitter.x0(), 30+_circleFitter.y0());
    crossHelixLineX->SetLineColor(kBlue);
    crossHelixLineY->SetLineColor(kBlue);
    crossHelixLineX->Draw("same");
    crossHelixLineY->Draw("same");
*/
    // Draw a cross line at the center of Helix: used for fitting
    TLine *crossHelixLineX = new TLine(-30+xC, yC, 30+xC, yC);
    TLine *crossHelixLineY = new TLine(xC, -30+yC, xC, 30+yC);
    crossHelixLineX->SetLineColor(kBlue);
    crossHelixLineY->SetLineColor(kBlue);
    crossHelixLineX->Draw("same");
    crossHelixLineY->Draw("same");
    // Draw a resolution on wire for each hit
    //TLine* LineW[nComboHitsInSegment];
    std::vector<TLine*> LineW(nComboHitsInSegment, nullptr);
      if(_tcHits[j].used == false) continue;
      //if(j != 16) continue; //bokuno
      int hitIndice = _tcHits[j].hitIndice;
      double fSigW = _data.chcol->at(hitIndice).wireRes();
      //double fSigR  = 2.5;
      double DirX = _data.chcol->at(hitIndice).uDir().x();
      double DirY = _data.chcol->at(hitIndice).uDir().y();
      double LineW_xy[4] = {_tcHits.at(j).x-DirX*fSigW, _tcHits.at(j).y-DirY*fSigW, _tcHits.at(j).x+DirX*fSigW, _tcHits.at(j).y+DirY*fSigW};
      LineW[j] = new TLine(LineW_xy[0], LineW_xy[1], LineW_xy[2], LineW_xy[3]);
      //TLine *LineR = new TLine(_circleFitter.x0(), -30+_circleFitter.y0(), _circleFitter.x0(), 30+_circleFitter.y0());
      LineW[j]->SetLineColor(kBlack);
      LineW[j]->SetLineWidth(4);
      LineW[j]->Draw("sames");
    //Draw wire slope
    std::vector<TLine*> SlopeWire(nComboHitsInSegment, nullptr);
    double slope_x[2] = {_tcHits.at(j).x-300.0, _tcHits.at(j).x+300.0};
    double SlopeWire_xy[4] = {slope_x[0], (slope_x[0]*_printweight[j].nhit_slope_a + _printweight[j].nhit_slope_c)/(-_printweight[j].nhit_slope_b), slope_x[1], (slope_x[1]*_printweight[j].nhit_slope_a + _printweight[j].nhit_slope_c)/(-_printweight[j].nhit_slope_b)};
    SlopeWire[j] = new TLine(SlopeWire_xy[0], SlopeWire_xy[1], SlopeWire_xy[2], SlopeWire_xy[3]);
    SlopeWire[j]->SetLineStyle(2);  // Set the line color to green
    SlopeWire[j]->SetLineColor(kBlack);
    SlopeWire[j]->Draw("sames");
    //Draw slope orthogonal to wire
    std::vector<TLine*> orthogonal_SlopeWire(nComboHitsInSegment, nullptr);
    double orthogonal_slope_x[2] = {xC-300.0, xC+300.0};
    double slope_y[2] = {(orthogonal_slope_x[0]*_printweight[j].ortho_nhit_slope_a + _printweight[j].ortho_nhit_slope_c)/(-_printweight[j].ortho_nhit_slope_b), (orthogonal_slope_x[1]*_printweight[j].ortho_nhit_slope_a + _printweight[j].ortho_nhit_slope_c)/(-_printweight[j].ortho_nhit_slope_b)};
    std::cout<<"x1/y1"<<orthogonal_slope_x[0]<<"/"<<slope_y[0]<<std::endl;
    std::cout<<"x2/y2"<<orthogonal_slope_x[1]<<"/"<<slope_y[1]<<std::endl;
    double orthogonal_SlopeWire_xy[4] = {orthogonal_slope_x[0], slope_y[0], orthogonal_slope_x[1], slope_y[1]};
    orthogonal_SlopeWire[j] = new TLine(orthogonal_SlopeWire_xy[0], orthogonal_SlopeWire_xy[1], orthogonal_SlopeWire_xy[2], orthogonal_SlopeWire_xy[3]);
    orthogonal_SlopeWire[j]->SetLineStyle(2);  // Set the line color to green
    orthogonal_SlopeWire[j]->SetLineColor(kBlue);
    orthogonal_SlopeWire[j]->Draw("sames");
    // Add title at the top
    TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
    title->AddText(Form("X vs. Y (Run-subRun-Event, TC, #candidate) = (%d-%d-%d, %d, %d) (pbar1b0)", run, subrun, eventNumber, tc, isegment));
    title->SetFillColor(0);
    title->SetTextAlign(22);
    title->Draw("same");
    TLatex legend;
    legend.SetTextSize(0.03);
    double w = 1.0/_tcHits[j].circleError2;
    legend.DrawLatexNDC(0.18, 0.85, Form("%d CHs, #chi^{2}/ndf = %.3f, hit #%d, w = %.6f", (int)_circleFitter.qn(), _circleFitter.chi2DofCircle(), (int)j, w));
    //legend.DrawLatexNDC(0.18, 0.85, Form("%d CHs, #chi^{2}/ndf = %.3f", (int)_circleFitter.qn(), _circleFitter.chi2DofCircle()));
    legend.DrawLatexNDC(0.18, 0.82, Form("(xC, yC, rC) = (%.1f, %.1f, %.1f)", _circleFitter.x0(), _circleFitter.y0(), _circleFitter.radius()));
    legend.DrawLatexNDC(0.18, 0.79, Form("%s", filename));
    TLegend *legend1 = new TLegend(0.18,0.12,0.48,0.19);
    legend1->AddEntry("circle3","fitting result","l");
    legend1->AddEntry("circle4","used for weight correction","l");
    legend1->SetLineColor(kWhite);
    legend1->Draw("same");
    canvas4->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/findHelix/%s/pbar_Helix_%04d-%04d-%04d_TC-%d_cand_%d_hit_%d.pdf", filename, run, subrun, eventNumber, tc, isegment, (int)j));
    //canvas4->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/ce/findHelix/%s/pbar_Helix_%04d-%04d-%04d_TC-%d_cand_%d_hit_%d.pdf", filename, run, subrun, eventNumber, tc, isegment, (int)j));
    delete canvas4;
    delete gr1;
    delete gr2;
    delete title;
    delete circle1;
    delete circle2;
    delete circle3;
    delete circle4;
    delete circle5;
    delete legend1;
    delete crossLineX;
    delete crossLineY;
    delete crossHelixLineX;
    delete crossHelixLineY;
    }
}
//-----------------------------------------------------------------------------
// updates circleError of hit given some circle parameters
//-----------------------------------------------------------------------------
void PhiZSeedFinder::mod_computeCircleError2(size_t& tcHitsIndex, double& xC, double& yC, double& rC) {
      int hitIndice = _tcHits[tcHitsIndex].hitIndice;
      double transVar = _data.chcol->at(hitIndice).transVar();
      double x = _tcHits.at(tcHitsIndex).x;
      double y = _tcHits.at(tcHitsIndex).y;
      double dx = x - xC;
      double dy = y - yC;
      double dxn = dx * _data.chcol->at(hitIndice).vDir().x() + dy * _data.chcol->at(hitIndice).vDir().y();
      double costh2 = dxn * dxn / (dx * dx + dy * dy);
      double sinth2 = 1 - costh2;
      _tcHits[tcHitsIndex].circleError2 = _data.chcol->at(hitIndice).wireVar() * sinth2 + transVar * costh2;
      std::cout<<rC<<std::endl;
}
//-----------------------------------------------------------------------------
// updates circleError of hit given some circle parameters
//-----------------------------------------------------------------------------
void PhiZSeedFinder::computeCircleError2(size_t& tcHitsIndex, double& xC, double& yC, double& rC) {
  int hitIndice = _tcHits[tcHitsIndex].hitIndice;
  //default weight calculation
  double transVar = _data.chcol->at(hitIndice).transVar();
  double x = _tcHits.at(tcHitsIndex).x;
  double y = _tcHits.at(tcHitsIndex).y;
  double dx = x - xC;
  double dy = y - yC;
  double dxn = dx * _data.chcol->at(hitIndice).vDir().x() + dy * _data.chcol->at(hitIndice).vDir().y();
  double costh2 = dxn * dxn / (dx * dx + dy * dy);
  double sinth2 = 1 - costh2;
  _tcHits[tcHitsIndex].circleError2 = _data.chcol->at(hitIndice).wireVar() * sinth2 + transVar * costh2;
  //parameters for debug mode
  double deltaR = std::abs(rC - std::sqrt((x - xC) * (x - xC) + (y - yC) * (y - yC)));
  double wireErr = _data.chcol->at(hitIndice).wireRes();
  double wireVecX = _data.chcol->at(hitIndice).uDir().x();
  double wireVecY = _data.chcol->at(hitIndice).uDir().y();
  double transErr = _data.chcol->at(hitIndice).transRes();
  double chi2_default = deltaR*deltaR/_tcHits[tcHitsIndex].circleError2;
  double sigma_default = sqrt(_tcHits[tcHitsIndex].circleError2);
  double weight_default = 1.0/_tcHits[tcHitsIndex].circleError2;
  int nstrawhits = _tcHits[tcHitsIndex].strawhits;
  double sigma_wire = wireErr;
  double sigma_transverse = transErr;
  int intersection = 0;
  double chi2_new = 0.0;
  double residula_wire = 0.0;
  double residula_transverse = 0.0;
  double sigma_new = 0.0;
  double weight_new = 0.0;
  weightinfo hit;
  //method 2
  //liner equation: ax+by+c = 0
  //wire liner equation
  double wire_slope = wireVecY/wireVecX;
  double wire_a = wire_slope;
  double wire_b = -1.0;
  double wire_c = y - wire_a*x;
  //std::cout<<"wire_a/wire_b/wire_c = "<<wire_a<<"/"<<wire_b<<"/"<<wire_c<<std::endl;
  //distance between wire and helix center
  double distance_wire = fabs(wire_a*xC + wire_b*yC + wire_c)/sqrt(wire_a*wire_a + wire_b*wire_b);
  //std::cout<<"distance_wire  = "<<distance_wire <<std::endl;
  //orthogonal vector(unit vector) normal to wire vector
  //double trans_slope = -transVecY/transVecX;
  //double trans_a = trans_slope;
  //double trans_b = -1.0;
  //double trans_c = y + wire_a*x;
  double A = (wire_a*wire_a)/(wire_b*wire_b) + 1.0;
  double B = 2.0*((wire_a*wire_c)/(wire_b*wire_b) - xC + (wire_a/wire_b)*yC);
  double C = xC*xC + yC*yC + 2.0*wire_c*yC/wire_b + (wire_c*wire_c)/(wire_b*wire_b) - rC*rC;
  double D = B*B - 4.0*A*C;
  //if(distance_wire < rC) {
  if(D > 0.0) {
    std::cout<<"2 point "<<std::endl;
    //intersection: 2 points (x1, y1), (x2, y2)
    //std::cout<<"deltaDistance=  "<<deltaDistance<<std::endl;
    //std::cout<<"Default Chi2 = "<<deltaDistance*deltaDistance/_tcHits[tcHitsIndex].circleError2<<std::endl;
    double intersection_x[2] = {0.0};
    double intersection_y[2] = {0.0};
    //double d = fabs(wire_a*xC + wire_b*yC + wire_c);
    //std::cout<<"d = "<<d<<std::endl;
    //std::cout<<"A = "<<A<<std::endl;
    //std::cout<<"B = "<<B<<std::endl;
    //std::cout<<"C = "<<C<<std::endl;
    //std::cout<<"D = "<<D<<std::endl;
    //std::cout<<"wire_a = "<<wire_a<<std::endl;
    //std::cout<<"wire_b = "<<wire_b<<std::endl;
    //std::cout<<"wire_c = "<<wire_c<<std::endl;
    //std::cout<<"wire_a*d = "<<wire_a*d<<std::endl;
    //std::cout<<"wire_a*wire_a = "<<wire_a*wire_a<<std::endl;
    //std::cout<<"wire_b*wire_b = "<<wire_b*wire_b<<std::endl;
    //std::cout<<"rC*rC = "<<rC*rC<<std::endl;
    //std::cout<<"d*d = "<<d*d<<std::endl;
    //std::cout<<"sqrt((wire_a*wire_a + wire_b*wire_b)*rC*rC - d*d)) = "<<sqrt((wire_a*wire_a + wire_b*wire_b)*rC*rC - d*d)<<std::endl;
    //std::cout<<"wire_a*wire_a + wire_b*wire_b = "<<wire_a*wire_a + wire_b*wire_b<<std::endl;
    intersection_x[0] = (-B + sqrt(D))/(2.0*A);
    intersection_y[0] = wire_slope*intersection_x[0] + wire_c;
    intersection_x[1] = (-B - sqrt(D))/(2.0*A);
    intersection_y[1] = wire_slope*intersection_x[1] + wire_c;
    /*intersection_x[0] = (wire_a*d - wire_b*sqrt((wire_a*wire_a + wire_b*wire_b)*rC*rC - d*d))/(wire_a*wire_a + wire_b*wire_b) + xC;
    intersection_y[0] = (wire_b*d + wire_a*sqrt((wire_a*wire_a + wire_b*wire_b)*rC*rC - d*d))/(wire_a*wire_a + wire_b*wire_b) + yC;
    intersection_x[1] = (wire_a*d + wire_b*sqrt((wire_a*wire_a + wire_b*wire_b)*rC*rC - d*d))/(wire_a*wire_a + wire_b*wire_b) + xC;
    intersection_y[1] = (wire_b*d - wire_a*sqrt((wire_a*wire_a + wire_b*wire_b)*rC*rC - d*d))/(wire_a*wire_a + wire_b*wire_b) + yC;
    */
    //std::cout<<"intersection_x[0] = "<<intersection_x[0]<<std::endl;
    //std::cout<<"intersection_y[0] = "<<intersection_y[0]<<std::endl;
    //std::cout<<"intersection_x[1] = "<<intersection_x[1]<<std::endl;
    //std::cout<<"intersection_y[1] = "<<intersection_y[1]<<std::endl;
    double distance[4] = {0.0};
    distance[0] = sqrt((x - intersection_x[0]) * (x - intersection_x[0]) + (y - intersection_y[0]) * (y - intersection_y[0]));
    distance[1] = sqrt((x - intersection_x[1]) * (x - intersection_x[1]) + (y - intersection_y[1]) * (y - intersection_y[1]));
    double residual_distance = 0.0;
    if(distance[0] < distance[1]) residual_distance = distance[0];
    else residual_distance = distance[1];
    _tcHits[tcHitsIndex].circleError2 = (wireErr*wireErr) * (deltaR*deltaR) / (residual_distance*residual_distance);
    //std::cout<<"deltaDistance = "<<deltaDistance<<std::endl;
    //std::cout<<"distance[0] = "<<distance[0]<<std::endl;
    //std::cout<<"distance[1] = "<<distance[1]<<std::endl;
    std::cout<<"residual_distance = "<<residual_distance<<std::endl;
    //std::cout<<"circleError2 = "<<_tcHits[tcHitsIndex].circleError2<<std::endl;
    //std::cout<<"circleError = "<<sqrt(_tcHits[tcHitsIndex].circleError2)<<std::endl;
    //std::cout<<"1.0/circleError2 = "<<1.0/_tcHits[tcHitsIndex].circleError2<<std::endl;
    //std::cout<<"new Chi2 = "<<(residual_distance*residual_distance)/(wireErr*wireErr)<<std::endl;
    intersection = 2;
    chi2_new = (residual_distance*residual_distance)/(wireErr*wireErr);
    residula_wire = residual_distance;
    sigma_new = sqrt(_tcHits[tcHitsIndex].circleError2);
    weight_new = 1.0/_tcHits[tcHitsIndex].circleError2;
  //}else if(distance_wire == rC){
  }else if(D == 0.0){
    //std::cout<<"1 point"<<std::endl;
    intersection = 1;
  //intersection: 1 point (x1, y1)
   //std::cout<<"deltaDistance=  "<<deltaDistance<<std::endl;
  //std::cout<<"Default Chi2 = "<<deltaDistance*deltaDistance/_tcHits[tcHitsIndex].circleError2<<std::endl;
   double intersection_x = -B / (2.0*A);
   double intersection_y = wire_slope*intersection_x + wire_c;
   //double intersection_x = (wire_a*rC)/sqrt(wire_a*wire_a + wire_b*wire_b) + xC;
   //double intersection_y = (wire_b*rC)/sqrt(wire_a*wire_a + wire_b*wire_b) + yC;
   if(deltaR == 0) {
    _tcHits[tcHitsIndex].circleError2 = wireErr*wireErr;
    //std::cout<<"residual_distance = "<<deltaDistance<<std::endl;
    //std::cout<<"new Chi2 = "<<(deltaDistance*deltaDistance)/(wireErr*wireErr)<<std::endl;
    chi2_new = 0.0;
    residula_wire = deltaR;
    sigma_new = sqrt(_tcHits[tcHitsIndex].circleError2);
    weight_new = 1.0/_tcHits[tcHitsIndex].circleError2;
   }
   //if(deltaDistance == 0) _tcHits[tcHitsIndex].circleError2 = transErr*transErr;
   else {
    double residual_distance = sqrt((x - intersection_x) * (x - intersection_x) + (y - intersection_y) * (y - intersection_y));
    //std::cout<<"residual_distance = "<<residual_distance<<std::endl;
    _tcHits[tcHitsIndex].circleError2 = (wireErr*wireErr) * (deltaR*deltaR) / (residual_distance*residual_distance);
    //std::cout<<"new Chi2 = "<<(residual_distance*residual_distance)/(wireErr*wireErr)<<std::endl;
    chi2_new = (residual_distance*residual_distance)/(wireErr*wireErr);
    residula_wire = residual_distance;
    sigma_new = sqrt(_tcHits[tcHitsIndex].circleError2);
    weight_new = 1.0/_tcHits[tcHitsIndex].circleError2;
   }
   intersection = 1;
    //std::cout<<"circleError = "<<sqrt(_tcHits[tcHitsIndex].circleError2)<<std::endl;
    //std::cout<<"1.0/circleError2 = "<<1.0/_tcHits[tcHitsIndex].circleError2<<std::endl;
  }else{
    std::cout<<"else"<<std::endl;
  //intersection: 0 point
    //std::cout<<"deltaDistance=  "<<deltaDistance<<std::endl;
    //std::cout<<"Default Chi2 = "<<deltaDistance*deltaDistance/_tcHits[tcHitsIndex].circleError2<<std::endl;
    //distance between wire and helix radius
    double distance[2] = {0.0};//0:wire, 1:transverse
    //orthogonal vector(unit vector) normal to wire vector
    //double trans_a = wire_b;
    //double trans_b = -wire_a;
    //double trans_c = wire_b*xC - wire_a*yC;
    //double trans_a = -wire_b;
    //double trans_b = wire_a;
    double trans_a = wire_b/wire_a;
    double trans_b = -1.0;
    //double trans_c = wire_a*yC - wire_b*xC;
    //double trans_c = wire_b*xC + wire_a*yC;
    double trans_c = -wire_b/wire_a*xC + yC;
    if(_diagLevel > 0){
      hit.ortho_nhit_slope_a = trans_a;
      hit.ortho_nhit_slope_b = trans_b;
      hit.ortho_nhit_slope_c = trans_c;
    }
    std::cout<<"trans_a/trans_b/trans_c = "<<trans_a<<"/"<<trans_b<<"/"<<trans_c<<std::endl;
    std::cout<<"x/y = "<<x<<"/"<<y<<std::endl;
    std::cout<<"xC/yC = "<<xC<<"/"<<yC<<std::endl;
    std::cout<<"x1/y1 = "<<xC-300.0<<"/"<<(trans_a*(xC-300.0)+trans_c)/(-trans_b)<<std::endl;
    std::cout<<"x2/y2 = "<<xC+300.0<<"/"<<(trans_a*(xC+300.0)+trans_c)/(-trans_b)<<std::endl;
    //std::cout<<"OMG = "<<trans_a*wire_a + trans_b*wire_b<<std::endl;
    double intersection_x = (wire_b*trans_c - trans_b*wire_c)/(wire_a*trans_b - trans_a*wire_b);
    double intersection_y = (wire_c*trans_a - trans_c*wire_a)/(wire_a*trans_b - trans_a*wire_b);
    distance[0] = sqrt((x - intersection_x) * (x - intersection_x) + (y - intersection_y) * (y - intersection_y));
    distance[1] = fabs(rC - distance_wire);
    double chi2[2] = {0.0};
    chi2[0] = (distance[0]*distance[0]) / (wireErr*wireErr);
    chi2[1] = (distance[1]*distance[1]) / (transErr*transErr);
    _tcHits[tcHitsIndex].circleError2 = (deltaR*deltaR) / (chi2[0]+chi2[1]);
    //std::cout<<"deltaDistance = "<<deltaDistance<<std::endl;
    std::cout<<"intersection_x = "<<intersection_x<<std::endl;
    std::cout<<"intersection_y = "<<intersection_y<<std::endl;
    std::cout<<"distance[0] = "<<distance[0]<<std::endl;
    std::cout<<"distance[1] = "<<distance[1]<<std::endl;
    //std::cout<<"chi2[0] = "<<chi2[0]<<std::endl;
    //std::cout<<"chi2[1] = "<<chi2[1]<<std::endl;
    //std::cout<<"new Chi2 = "<<chi2[0]+chi2[1]<<std::endl;
    //std::cout<<"circleError1 = "<<sqrt(_tcHits[tcHitsIndex].circleError2)<<std::endl;
    //std::cout<<"1.0/circleError2 = "<<1.0/_tcHits[tcHitsIndex].circleError2<<std::endl;
    intersection = 0;
    chi2_new = chi2[0]+chi2[1];
    residula_wire = distance[0];//wire
    residula_transverse = distance[1];//transverse
    sigma_new = sqrt(_tcHits[tcHitsIndex].circleError2);
    weight_new = 1.0/_tcHits[tcHitsIndex].circleError2;
  }
  //std::cout<<" "<<std::endl;
  //std::cout<<"============================="<<std::endl;
  //std::cout<<" computeCircleError w/ correction "<<std::endl;
  //std::cout<<"============================="<<std::endl;
  //std::cout<<"radVecX = "<<radVecX<<std::endl;
  //std::cout<<"radVecY = "<<radVecY<<std::endl;
  //std::cout<<"wireErr = "<<wireErr<<std::endl;
  //std::cout<<"wireVecX = "<<wireVecX<<std::endl;
  //std::cout<<"wireVecY = "<<wireVecY<<std::endl;
  //std::cout<<"projWireErr = "<<projWireErr<<std::endl;
  //std::cout<<"transErr = "<<transErr<<std::endl;
  //std::cout<<"transVecX = "<<transVecX<<std::endl;
  //std::cout<<"transVecY = "<<transVecY<<std::endl;
  //std::cout<<"projTransErr = "<<projTransErr<<std::endl;
  //std::cout<<"circleError2 = "<<_tcHits[tcHitsIndex].circleError2<<std::endl;
  //std::cout<<"deltaDistance = "<<deltaDistance<<std::endl;
  //std::cout<<"1.0/circleError2 = "<<1.0/_tcHits[tcHitsIndex].circleError2<<std::endl;
  if(_diagLevel > 0){
    //weightinfo hit;
    hit.chi2_default = chi2_default;
    hit.deltaR = deltaR;
    hit.sigma_default = sigma_default;
    hit.weight_default = weight_default;
    hit.nstrawhits = nstrawhits;
    hit.sigma_wire = sigma_wire;//[mm]
    hit.sigma_transverse = sigma_transverse; //[mm]
    hit.intersection = intersection;
    hit.chi2_new = chi2_new;
    hit.residula_wire = residula_wire; //[mm]
    hit.residula_transverse = residula_transverse; //[mm]
    hit.sigma_new = sigma_new;//[mm]
    hit.weight_new = weight_new;
    hit.nhit_slope_a = wire_a;
    hit.nhit_slope_b = wire_b;
    hit.nhit_slope_c = wire_c;
    _printweight.push_back(hit);
    //std::cout<<"chi2_default = "<<chi2_default<<std::endl;
    //std::cout<<"deltaR = "<<deltaR<<std::endl;
    //std::cout<<"sigma_default = "<<sigma_default<<std::endl;
    //std::cout<<"weight_default = "<<weight_default<<std::endl;
    //std::cout<<"nstrawhits = "<<nstrawhits<<std::endl;
    //std::cout<<"sigma_wire = "<<sigma_wire<<std::endl;
    //std::cout<<"sigma_transverse = "<<sigma_transverse<<std::endl;
    //std::cout<<"intersection = "<<intersection<<std::endl;
    //std::cout<<"chi2_new = "<<chi2_new<<std::endl;
    //std::cout<<"residula_wire = "<<residula_wire<<std::endl;
    //std::cout<<"residula_transverse = "<<residula_transverse<<std::endl;
    //std::cout<<"sigma_new = "<<sigma_new<<std::endl;
    //std::cout<<"hit.ortho_nhit_slope_a = "<<hit.ortho_nhit_slope_a<<std::endl;
    //std::cout<<"hit.ortho_nhit_slope_b = "<<hit.ortho_nhit_slope_b<<std::endl;
    //std::cout<<"hit.ortho_nhit_slope_c = "<<hit.ortho_nhit_slope_c<<std::endl;
  }
   std::cout<<"KOMAMAMAMA= "<<std::endl;
}
//-----------------------------------------------------------------------------
// updates circleError of hit given some circle parameters
//-----------------------------------------------------------------------------
double PhiZSeedFinder::computeCircleError2_ver2(int hitIndice, int nStrawHits, double& xC, double& yC, double& rC) {
  //default weight calculation
  double transVar = _data.chcol->at(hitIndice).transVar();
  double x = _data.chcol->at(hitIndice).pos().x();
  double y = _data.chcol->at(hitIndice).pos().y();
  double dx = x - xC;
  double dy = y - yC;
  double dxn = dx * _data.chcol->at(hitIndice).vDir().x() + dy * _data.chcol->at(hitIndice).vDir().y();
  double costh2 = dxn * dxn / (dx * dx + dy * dy);
  double sinth2 = 1 - costh2;
  double circleError2 = _data.chcol->at(hitIndice).wireVar() * sinth2 + transVar * costh2;
  //parameters for debug mode
  double deltaR = std::abs(rC - std::sqrt((x - xC) * (x - xC) + (y - yC) * (y - yC)));
  double wireErr = _data.chcol->at(hitIndice).wireRes();
  double wireVecX = _data.chcol->at(hitIndice).uDir().x();
  double wireVecY = _data.chcol->at(hitIndice).uDir().y();
  double transErr = _data.chcol->at(hitIndice).transRes();
  double chi2_default = deltaR*deltaR/circleError2;
  double sigma_default = sqrt(circleError2);
  double weight_default = 1.0/circleError2;
  int nstrawhits = nStrawHits;
  double sigma_wire = wireErr;
  double sigma_transverse = transErr;
  int intersection = 0;
  double chi2_new = 0.0;
  double residula_wire = 0.0;
  double residula_transverse = 0.0;
  double sigma_new = 0.0;
  double weight_new = 0.0;
  weightinfo hit;
  //method 2
  //liner equation: ax+by+c = 0
  //wire liner equation
  double wire_slope = wireVecY/wireVecX;
  double wire_a = wire_slope;
  double wire_b = -1.0;
  double wire_c = y - wire_a*x;
  double distance_wire = fabs(wire_a*xC + wire_b*yC + wire_c)/sqrt(wire_a*wire_a + wire_b*wire_b);
  double A = (wire_a*wire_a)/(wire_b*wire_b) + 1.0;
  double B = 2.0*((wire_a*wire_c)/(wire_b*wire_b) - xC + (wire_a/wire_b)*yC);
  double C = xC*xC + yC*yC + 2.0*wire_c*yC/wire_b + (wire_c*wire_c)/(wire_b*wire_b) - rC*rC;
  double D = B*B - 4.0*A*C;
  if(D > 0.0) {
    double intersection_x[2] = {0.0};
    double intersection_y[2] = {0.0};
    intersection_x[0] = (-B + sqrt(D))/(2.0*A);
    intersection_y[0] = wire_slope*intersection_x[0] + wire_c;
    intersection_x[1] = (-B - sqrt(D))/(2.0*A);
    intersection_y[1] = wire_slope*intersection_x[1] + wire_c;
    double distance[4] = {0.0};
    distance[0] = sqrt((x - intersection_x[0]) * (x - intersection_x[0]) + (y - intersection_y[0]) * (y - intersection_y[0]));
    distance[1] = sqrt((x - intersection_x[1]) * (x - intersection_x[1]) + (y - intersection_y[1]) * (y - intersection_y[1]));
    double residual_distance = 0.0;
    if(distance[0] < distance[1]) residual_distance = distance[0];
    else residual_distance = distance[1];
    circleError2 = (wireErr*wireErr) * (deltaR*deltaR) / (residual_distance*residual_distance);
    intersection = 2;
    chi2_new = (residual_distance*residual_distance)/(wireErr*wireErr);
    residula_wire = residual_distance;
    sigma_new = sqrt(circleError2);
    weight_new = 1.0/circleError2;
  }else if(D == 0.0){
    intersection = 1;
   double intersection_x = -B / (2.0*A);
   double intersection_y = wire_slope*intersection_x + wire_c;
   if(deltaR == 0) {
    circleError2 = wireErr*wireErr;
    chi2_new = 0.0;
    residula_wire = deltaR;
    sigma_new = sqrt(circleError2);
    weight_new = 1.0/circleError2;
   }
   else {
    double residual_distance = sqrt((x - intersection_x) * (x - intersection_x) + (y - intersection_y) * (y - intersection_y));
    circleError2 = (wireErr*wireErr) * (deltaR*deltaR) / (residual_distance*residual_distance);
    chi2_new = (residual_distance*residual_distance)/(wireErr*wireErr);
    residula_wire = residual_distance;
    sigma_new = sqrt(circleError2);
    weight_new = 1.0/circleError2;
   }
   intersection = 1;
  }else{
    //distance between wire and helix radius
    double distance[2] = {0.0};//0:wire, 1:transverse
    //orthogonal vector(unit vector) normal to wire vector
    double trans_a = wire_b/wire_a;
    double trans_b = -1.0;
    double trans_c = -wire_b/wire_a*xC + yC;
    if(_diagLevel > 0){
      hit.ortho_nhit_slope_a = trans_a;
      hit.ortho_nhit_slope_b = trans_b;
      hit.ortho_nhit_slope_c = trans_c;
    }
    double intersection_x = (wire_b*trans_c - trans_b*wire_c)/(wire_a*trans_b - trans_a*wire_b);
    double intersection_y = (wire_c*trans_a - trans_c*wire_a)/(wire_a*trans_b - trans_a*wire_b);
    distance[0] = sqrt((x - intersection_x) * (x - intersection_x) + (y - intersection_y) * (y - intersection_y));
    distance[1] = fabs(rC - distance_wire);
    double chi2[2] = {0.0};
    chi2[0] = (distance[0]*distance[0]) / (wireErr*wireErr);
    chi2[1] = (distance[1]*distance[1]) / (transErr*transErr);
    circleError2 = (deltaR*deltaR) / (chi2[0]+chi2[1]);
    intersection = 0;
    chi2_new = chi2[0]+chi2[1];
    residula_wire = distance[0];//wire
    residula_transverse = distance[1];//transverse
    sigma_new = sqrt(circleError2);
    weight_new = 1.0/circleError2;
  }
  if(_diagLevel > 0){
    //weightinfo hit;
    hit.chi2_default = chi2_default;
    hit.deltaR = deltaR;
    hit.sigma_default = sigma_default;
    hit.weight_default = weight_default;
    hit.nstrawhits = nstrawhits;
    hit.sigma_wire = sigma_wire;//[mm]
    hit.sigma_transverse = sigma_transverse; //[mm]
    hit.intersection = intersection;
    hit.chi2_new = chi2_new;
    hit.residula_wire = residula_wire; //[mm]
    hit.residula_transverse = residula_transverse; //[mm]
    hit.sigma_new = sigma_new;//[mm]
    hit.weight_new = weight_new;
    hit.nhit_slope_a = wire_a;
    hit.nhit_slope_b = wire_b;
    hit.nhit_slope_c = wire_c;
    _printweight.push_back(hit);
  }
  return circleError2;
}
//-----------------------------------------------------------------------------
double PhiZSeedFinder::computeCircleResidual2(size_t& tcHitsIndex, double& xC, double& yC, double& rC) {
  double xP = _tcHits.at(tcHitsIndex).x;
  double yP = _tcHits.at(tcHitsIndex).y;
  double deltaDistance = std::abs(rC - std::sqrt((xP - xC) * (xP - xC) + (yP - yC) * (yP - yC)));
  //double circleSigma2 = _tcHits[tcHitsIndex].circleError2;
  //return deltaDistance * deltaDistance / circleSigma2;
  return deltaDistance;
}
//-----------------------------------------------------------------------------
  //void PhiZSeedFinder::findHelix(int tc, int isegment, HelixSeedCollection& HSColl){
  void PhiZSeedFinder::findHelix(int tc, int isegment, HelixSeedCollection& HSColl, HelixSeed& Temp_HSeed){
    _circleFitter.clear();
    //Step1: fit circle of i-th segment with fixed weight
    size_t nComboHitsInSegment = _tcHits.size();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 0.1;//tentative value
      _circleFitter.addPoint(x, y, wP);
      //_tcHits[i].used = true;
    }
    double xC = _circleFitter.x0();
    double yC = _circleFitter.y0();
    double rC = _circleFitter.radius();
    if (_diagLevel > 0) {
      std::cout<<"==================================="<<std::endl;
      std::cout<<"        Step1                      "<<std::endl;
      std::cout<<"==================================="<<std::endl;
      std::cout<<"# of hits in Helix = "<<_circleFitter.qn()<<std::endl;
      std::cout<<"xC/yC = "<<_circleFitter.x0()<<"/"<<_circleFitter.y0()<<std::endl;
      std::cout<<"radius = "<<_circleFitter.radius()<<std::endl;
      std::cout<<"phi/dfdz/chi2DofC/chi2DofLineC = "<<_circleFitter.phi0()<<"/"<<_circleFitter.dfdz()<<"/"<<_circleFitter.chi2DofCircle()<<"/"<<_circleFitter.chi2DofLine()<<std::endl;
      _data.h_circleFitter_chi2Dof[0].push_back(_circleFitter.chi2DofCircle());
      _data.h_circleFitter_nhits[0].push_back((int)_circleFitter.qn());
    }
    plot_XVsY(tc, isegment, "step1", xC, yC, rC);
    //plot_XVsY_hit(tc, isegment, "step1", xC, yC, rC);
    //Step2: fit circle of i-th segment with correct weight
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    _circleFitter.clear();
    _printweight.clear();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      std::cout<<"i = "<<i<<std::endl;
      computeCircleError2(i, xC, yC, rC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
      //_tcHits[i].used = true;
    }
    if (_diagLevel > 0) {
      std::cout<<"==================================="<<std::endl;
      std::cout<<"              Step2                "<<std::endl;
      std::cout<<"==================================="<<std::endl;
      std::cout<<"# of hits in Helix = "<<_circleFitter.qn()<<std::endl;
      std::cout<<"xC/yC = "<<_circleFitter.x0()<<"/"<<_circleFitter.y0()<<std::endl;
      std::cout<<"radius = "<<_circleFitter.radius()<<std::endl;
      std::cout<<"phi/dfdz/chi2DofC/chi2DofLineC = "<<_circleFitter.phi0()<<"/"<<_circleFitter.dfdz()<<"/"<<_circleFitter.chi2DofCircle()<<"/"<<_circleFitter.chi2DofLine()<<std::endl;
      _data.h_circleFitter_chi2Dof[1].push_back(_circleFitter.chi2DofCircle());
      _data.h_circleFitter_nhits[1].push_back((int)_circleFitter.qn());
      std::cout<<"chi2_default | deltaR | sigma_default | weight_default | nstrawhits | sigma_wire | sigma_transverse | intersection | chi2_new | residula_wire | residula_transverse | sigma_new | weight_new"<<std::endl;
      for(size_t i=0; i<_printweight.size(); i++){
      //printf("%4i %6.6f %6.6f %6.6f %6.6f %3i  %6.6f %6.6f %2i %6.6f %6.6f %6.6f %6.6f %6.6f\n",
        printf("%-5d %-8.6f %-10.6f %-13.6f %-10.6f %-5d %-10.6f %-10.6f %-5d %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f\n",
        (int)i,
       _printweight[i].chi2_default,
       _printweight[i].deltaR,
       _printweight[i].sigma_default,
       _printweight[i].weight_default,
       _printweight[i].nstrawhits,
       _printweight[i].sigma_wire,
       _printweight[i].sigma_transverse,
       _printweight[i].intersection,
       _printweight[i].chi2_new,
       _printweight[i].residula_wire,
       _printweight[i].residula_transverse,
       _printweight[i].sigma_new,
       _printweight[i].weight_new);
      }
    }
    plot_XVsY(tc, isegment, "step2", xC, yC, rC);
    //plot_XVsY_hit(tc, isegment, "step2", xC, yC, rC);
    // clean up hits in the circle
    //Step3: iterate over all combohits and find the best hits combination when it has a good chi2/ndf
    //int min_hit = 10;// minimum number of combohits to do a circle fit
    if(nComboHitsInSegment > 10 and _circleFitter.chi2DofCircle() > 5.0){
    //std::vector<cleanup> remove_hits;
   // double init_chi2ndf = _circleFitter.chi2DofCircle();
   // for(size_t i=0; i<nComboHitsInSegment; i++){
   //   _circleFitter.clear();
   //   for(size_t j=0; j<nComboHitsInSegment; j++){
   //     if(i==j) continue;
   //     double x = _tcHits.at(j).x;
   //     double y = _tcHits.at(j).y;
   //     double wP = 1.0 / (_tcHits[j].circleError2);
   //     _circleFitter.addPoint(x, y, wP);
   //   }
   //   if(_circleFitter.chi2DofCircle() < init_chi2ndf) {
   //     init_chi2ndf = _circleFitter.chi2DofCircle();
   //     cleanup circlefit;
   //     circlefit.tcindex = i;
   //     circlefit.chi2ndf = _circleFitter.chi2DofCircle();
   //     remove_hits.push_back(circlefit);
   //   }
   // }
   // std::sort(remove_hits.begin(), remove_hits.end(), [](const cleanup& a, const cleanup& b) { return a.chi2ndf < b.chi2ndf; } );
   //   std::cout<<"==================================="<<std::endl;
   //   std::cout<<"clean-up"<<std::endl;
   //   std::cout<<"==================================="<<std::endl;
   // for(size_t i=0; i<remove_hits.size(); i++){
   //   int tcindex = remove_hits.at(i).tcindex;
   //   double chi2ndf = remove_hits.at(i).chi2ndf;
   //   if(chi2ndf < 5.0) _tcHits[tcindex].used = false;
   //   std::cout<<"No: "<<tcindex<< ", chi2DofCircle = "<<remove_hits.at(i).chi2ndf<<std::endl;
   // }
    std::vector<cleanup> remove_hits;
    double chi2ndf = _circleFitter.chi2DofCircle();
    while(chi2ndf > 5.0){
      int remove_hitIndex;
      int find = 0;
      //Level0
      for(size_t i=0; i<nComboHitsInSegment; i++){
        if(_tcHits[i].used == false) continue;
        _circleFitter.clear();
        for(size_t j=0; j<nComboHitsInSegment; j++){
          if(_tcHits[j].used == false) continue;
          if(i==j) continue;
          double x = _tcHits.at(j).x;
          double y = _tcHits.at(j).y;
          double wP = 1.0 / (_tcHits[j].circleError2);
          _circleFitter.addPoint(x, y, wP);
        }
        if(_circleFitter.chi2DofCircle() < chi2ndf) {
          remove_hitIndex = i;
          chi2ndf = _circleFitter.chi2DofCircle();
          find++;
        std::cout<<"i :"<<i<<"/"<<"chi2ndf: "<<_circleFitter.chi2DofCircle()<<std::endl;
        std::cout<<"find :"<<find<<std::endl;
        }
      }
    std::cout<<"find :"<<find<<std::endl;
    //check if chi2ndf is improved or not
    if(find > 0){
      cleanup circlefit;
      circlefit.tcindex = remove_hitIndex;
      circlefit.chi2ndf = chi2ndf;
      remove_hits.push_back(circlefit);
    _tcHits[remove_hitIndex].used = false;
    //recalculate the circle parameter
    //Level 1:
    _circleFitter.clear();
    for(size_t i=0; i<nComboHitsInSegment; i++){
    if(_tcHits[i].used == false) continue;
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
    }
    chi2ndf = _circleFitter.chi2DofCircle();
    std::cout<<"remove_hitIndex :"<<remove_hitIndex<<"/"<<"chi2ndf: "<<chi2ndf<<std::endl;
    if(chi2ndf < 5.0) break;
    if((nComboHitsInSegment - (int)remove_hits.size()) <= 10) break;
    }else{
    break;
    }
    }
    //std::sort(remove_hits.begin(), remove_hits.end(), [](const cleanup& a, const cleanup& b) { return a.chi2ndf < b.chi2ndf; } );
    std::cout<<"==================================="<<std::endl;
    std::cout<<"clean-up"<<std::endl;
    std::cout<<"==================================="<<std::endl;
    for(size_t i=0; i<remove_hits.size(); i++){
    int tcindex = remove_hits.at(i).tcindex;
    double chi2ndf_ = remove_hits.at(i).chi2ndf;
    std::cout<<"No: "<<tcindex<< ", chi2DofCircle = "<<chi2ndf_<<std::endl;
    }
    //Step4: iterate over selected combohits and recalculate the weight value and refit again
    //recalculate the weight
    std::cout<<"==================================="<<std::endl;
    std::cout<<"After clean-up"<<std::endl;
    std::cout<<"==================================="<<std::endl;
    _circleFitter.clear();
    int count = 0;
    for(size_t i=0; i<nComboHitsInSegment; i++){
    if(_tcHits[i].used == false) continue;
    double x = _tcHits.at(i).x;
    double y = _tcHits.at(i).y;
    double wP = 1.0 / (_tcHits[i].circleError2);
    _circleFitter.addPoint(x, y, wP);
    count++;
    }
    std::cout<<"nhit used in circle:"<<count<< ", chi2DofCircle = "<<_circleFitter.chi2DofCircle()<<std::endl;
    // xC = _circleFitter.x0();
    // yC = _circleFitter.y0();
    // _circleFitter.clear();
    // for(size_t i=0; i<nComboHitsInSegment; i++){
    //   if(_tcHits[i].used == false) continue;
    //   computeCircleError2(i, xC, yC);
    //   double x = _tcHits.at(i).x;
    //   double y = _tcHits.at(i).y;
    //   double wP = 1.0 / (_tcHits[i].circleError2);
    //   _circleFitter.addPoint(x, y, wP);
    //   _tcHits[i].used = true;
    // }
    }
    //if(_circleFitter.chi2DofCircle() > 5) plot_XVsY(tc, isegment, "step3");
    //Step3: fit circle of i-th segment with correct weight
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    _circleFitter.clear();
    _printweight.clear();
    std::cout<<"==================================="<<std::endl;
    std::cout<<"            step3_start            "<<std::endl;
    std::cout<<"==================================="<<std::endl;
    for(size_t i=0; i<nComboHitsInSegment; i++){
    if(_tcHits[i].used == false) continue;
    computeCircleError2(i, xC, yC, rC);
    computeHelixPhi(i, xC, yC);
    double x = _tcHits.at(i).x;
    double y = _tcHits.at(i).y;
    //double z = _tcHits.at(i).z;
    //double phi = _tcHits.at(i).helixPhi;
    double wP = 1.0 / (_tcHits[i].circleError2);
    //double seedError2 = 0.1;
    //double seedWeight = 1.0/(seedError2);
    //if(i == 3)  wP = 0.007863;
    //if(i == 29)  wP = 0.4799;
    //if(i == 37) wP = 0.003453;
    //if(i == 40) wP = 0.461721;
    _circleFitter.addPoint(x, y, wP);
    //_tcHits[i].used = true;
    //hit = &helixData._chHitsToProcess[f];
    //ComboHit                hhit(*hit);
    //helixData._hseed._hhits.push_back(hhit);
    //hit = &_data.chcol->at(_tcHits[i].hitIndice);
    //ComboHit                hhit(*hit);
    //ComboHit hhit = _data.chcol->at(_tcHits[i].hitIndice);
    //Temp_HSeed._hhits.push_back(hhit);
    }
    plot_XVsY(tc, isegment, "step3", xC, yC, rC);
    _FitRadius.push_back(_circleFitter.radius());
    std::cout<<"==================================="<<std::endl;
    std::cout<<"Previous step3 "<<std::endl;
    std::cout<<"==================================="<<std::endl;
    std::cout<<"xC/yC = "<<xC<<"/"<<yC<<std::endl;
    //plot_XVsY_hit(tc, isegment, "step3", xC, yC, rC);
    if (_diagLevel > 0) {
    std::cout<<"==================================="<<std::endl;
    std::cout<<"Step3 "<<std::endl;
    std::cout<<"==================================="<<std::endl;
    std::cout<<"# of hits in Helix = "<<_circleFitter.qn()<<std::endl;
    std::cout<<"xC/yC = "<<_circleFitter.x0()<<"/"<<_circleFitter.y0()<<std::endl;
    std::cout<<"radius = "<<_circleFitter.radius()<<std::endl;
    std::cout<<"phi/dfdz/chi2DofCircle/chi2DofLine = "<<_circleFitter.phi0()<<"/"<<_circleFitter.dfdz()<<"/"<<_circleFitter.chi2DofCircle()<<"/"<<_circleFitter.chi2DofLine()<<std::endl;
    //if(_circleFitter.chi2DofCircle() > 5){
    _data.h_circleFitter_chi2Dof[2].push_back(_circleFitter.chi2DofCircle());
    _data.h_circleFitter_nhits[2].push_back((int)_circleFitter.qn());
    //std::cout<<"bad circle fit"<<std::endl;
    //std::cout<<"run subrun event :"<<run<<" "<<subrun<<" "<<eventNumber<<std::endl;
    std::cout<<"chi2_default | deltaR | sigma_default | weight_default | nstrawhits | sigma_wire | sigma_transverse | intersection | chi2_new | residula_wire | residula_transverse | sigma_new | weight_new"<<std::endl;
    //std::cout<<"_printweight.size() = "<<_printweight.size()<<std::endl;
    for(size_t i=0; i<_printweight.size(); i++){
    printf("%d %-8.6f %-10.6f %-13.6f %-10.6f %-5d %-10.6f %-10.6f %-5d %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f\n",
    //printf("%4i %6.6f %6.6f %6.6f %6.6f %3i  %6.6f %6.6f %2i %6.6f %6.6f %6.6f %6.6f %6.6f\n",
    (int)i,
    _printweight[i].chi2_default,
    _printweight[i].deltaR,
    _printweight[i].sigma_default,
    _printweight[i].weight_default,
    _printweight[i].nstrawhits,
    _printweight[i].sigma_wire,
    _printweight[i].sigma_transverse,
    _printweight[i].intersection,
    _printweight[i].chi2_new,
    _printweight[i].residula_wire,
    _printweight[i].residula_transverse,
    _printweight[i].sigma_new,
    _printweight[i].weight_new);
    }
    }
    /*xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    double rC = _circleFitter.radius();
    _circleFitter.clear();
    for(size_t i=0; i<nComboHitsInSegment; i++){
    std::cout<<"computeCircleError2 = "<<computeCircleResidual2(i, xC, yC, rC)<<std::endl;
    if(computeCircleResidual2(i, xC, yC, rC) < 40){
    double x = _tcHits.at(i).x;
    double y = _tcHits.at(i).y;
    double wP = 1.0 / (_tcHits[i].circleError2);
    _circleFitter.addPoint(x, y, wP);
    _tcHits[i].used = true;
    }else{
    _tcHits[i].used = false;
    }
    }*/
    /*if (_diagLevel > 0) {
    std::cout<<"Step3: after clean-up the cirlce "<<std::endl;
    std::cout<<"# of hits in Helix = "<<_circleFitter.qn()<<std::endl;
    std::cout<<"xC/yC = "<<_circleFitter.x0()<<"/"<<_circleFitter.y0()<<std::endl;
    std::cout<<"radius = "<<_circleFitter.radius()<<std::endl;
    std::cout<<"phi/dfdz/chi2DofC/chi2DofLineC = "<<_circleFitter.phi0()<<"/"<<_circleFitter.dfdz()<<"/"<<_circleFitter.chi2DofCircle()<<"/"<<_circleFitter.chi2DofLine()<<std::endl;
    _data.h_circleFitter_chi2Dof[2].push_back(_circleFitter.chi2DofCircle());
    }*/
    //int loopCondition;
    //initTriplet(tripletInfo, loopCondition);
    //initSeedCircle(loopCondition);
    std::cout<<"kitagawa_01"<<std::endl;
    //initHelixPhi();
    std::cout<<"kitagawa_02"<<std::endl;
    //if (_diagLevel > 0) plot_PhiVsZ(isegment);//(1)HelixPhi vs. Z, (2)Phi vs. Z, can be removed in the future
    //findSeedPhiLines();
    //resolve2PiAmbiguities();
    //initFinalSeed();
    //HelixSeed hseed;
    std::cout<<"kitagawa_03"<<std::endl;
    //hseed._t0 = tccol1->at(tc)._t0;
    //auto _tcCollH = _event->getValidHandle<TimeClusterCollection>(_tcLabel);
    //hseed._timeCluster = art::Ptr<mu2e::TimeCluster>(_tcCollH, tc);
    //hseed._hhits.setParent(_chColl->parent());
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    /*hseed._helix._radius = rC;
    hseed._helix._rcent  = sqrt(xC*xC + yC*yC);
    hseed._helix._fcent  = polyAtan2(yC, xC);
    hseed._helix._lambda = 1.0/_dphidz;
    hseed._helix._fz0    = fz0;
    hseed._helix._helicity = hseed._helix._lambda > 0 ? Helicity::poshel : Helicity::neghel;
    */
    Temp_HSeed._helix._radius = rC;
    Temp_HSeed._helix._rcent  = sqrt(xC*xC + yC*yC);
    Temp_HSeed._helix._fcent  = polyAtan2(yC, xC);
    //Temp_HSeed._helix._lambda = 1.0/_dphidz;
    //Temp_HSeed._helix._fz0    = fz0;
    Temp_HSeed._helix._helicity = Temp_HSeed._helix._lambda > 0 ? Helicity::poshel : Helicity::neghel;
    std::cout<<"xC/yC/rC = "<<xC<<"/"<<yC<<"/"<<rC<<std::endl;
    std::cout<<"r/f cent = "<<sqrt(xC*xC + yC*yC)<<"/"<<polyAtan2(yC, xC)<<std::endl;
    // auto _tcCollH = _event->getValidHandle<TimeClusterCollection>(_tcLabel);
    //hseed._timeCluster = art::Ptr<mu2e::TimeCluster>(_tcCollH, tc);
    //hseed._hhits.setParent(_chColl->parent());
    // include also the values of the chi2
    //hseed._helix._chi2dXY = _circleFitter.chi2DofCircle();
    Temp_HSeed._helix._chi2dXY = _circleFitter.chi2DofCircle();
    //Temp_HSeed._helix._chi2dZPhi = _lineFitter.chi2Dof();
    std::cout<<"_circleFitter.chi2DofCircle() = "<<_circleFitter.chi2DofCircle()<<std::endl;
    //Temp_HSeed = hseed;
    // push back the helix seed to the helix seed collection
    //HSColl.emplace_back(hseed);
    std::cout<<"kitagawa_04"<<std::endl;
    }
  void PhiZSeedFinder::saveHelix(int tc, HelixSeed& Temp_HSeed){
    std::cout << "saveHelix " << std::endl;
    std::cout << "Temp_HSeed._hhits.size() = " << Temp_HSeed._hhits.size() << std::endl;
for (size_t j = 0; j < Temp_HSeed._hhits.size(); ++j) {
  const mu2e::ComboHit& hit = Temp_HSeed._hhits.at(j);
  std::cout << "  [Hit " << j << "] pos = ("
            << hit.pos().x() << ", "
            << hit.pos().y() << ", "
            << hit.pos().z() << "), "
            << "nStrawHits = " << hit.nStrawHits() << ", "
            << "time = " << hit.time() << ", "
            << "correctedTime = " << hit.correctedTime() << ", "
            << "phi = " << hit.phi() << ", "
            << std::endl;
}
    //HelixSeed hseed;
    const TimeCluster* tculster = &_data.tccol->at(tc);
    Temp_HSeed._t0 = tculster->_t0;
    auto const& tccH    = _event->getValidHandle<mu2e::TimeClusterCollection>(_tcCollTag);
    Temp_HSeed._timeCluster = art::Ptr<TimeCluster>(tccH, tc);
    Temp_HSeed._hhits.setParent(_data.chcol->parent());
    // flag hits used in helix, and push to combo hit collection in helix seed
    // also add points to linear fitter to get t0
    std::cout << "Temp_HSeed._hhits.size() = " << Temp_HSeed._hhits.size() << std::endl;
for (size_t j = 0; j < Temp_HSeed._hhits.size(); ++j) {
  const mu2e::ComboHit& hit = Temp_HSeed._hhits.at(j);
  std::cout << "  [Hit " << j << "] pos = ("
            << hit.pos().x() << ", "
            << hit.pos().y() << ", "
            << hit.pos().z() << "), "
            << "nStrawHits = " << hit.nStrawHits() << ", "
            << "time = " << hit.time() << ", "
            << "correctedTime = " << hit.correctedTime() << ", "
            << "phi = " << hit.phi() << ", "
            << std::endl;
}
    ::LsqSums2 fitter;
    std::cout<<"_tcHits.size() = "<<_tcHits.size()<<std::endl;
    for (size_t i = 0; i < _tcHits.size(); i++) {
      if(_tcHits[i].used == false) continue;
      int hitIndice = _tcHits[i].hitIndice;
      const ComboHit* hit = &_data.chcol->at(hitIndice);
      fitter.addPoint(hit->pos().z(), hit->correctedTime(), 1 / (hit->timeRes() * hit->timeRes()));
      ComboHit hhit(*hit);
      hhit._hphi = _tcHits[i].helixPhi + _tcHits[i].helixPhiCorrection * 2 * M_PI;
      Temp_HSeed._hhits.push_back(hhit);
    }
    std::cout << "Temp_HSeed._hhits.size() = " << Temp_HSeed._hhits.size() << std::endl;
for (size_t j = 0; j < Temp_HSeed._hhits.size(); ++j) {
  const mu2e::ComboHit& hit = Temp_HSeed._hhits.at(j);
  std::cout << "  [Hit " << j << "] pos = ("
            << hit.pos().x() << ", "
            << hit.pos().y() << ", "
            << hit.pos().z() << "), "
            << "nStrawHits = " << hit.nStrawHits() << ", "
            << "time = " << hit.time() << ", "
            << "correctedTime = " << hit.correctedTime() << ", "
            << "phi = " << hit.phi() << ", "
            << std::endl;
}
    //float eDepAvg = Temp_HSeed._hhits.eDepAvg();
    Temp_HSeed._t0 = TrkT0(fitter.y0(), fitter.y0Err());
    Temp_HSeed._status.merge(TrkFitFlag::helixOK);
    Temp_HSeed._status.merge(TrkFitFlag::APRHelix);
    // take care of plotting if _diagLevel = 1
    if (_diagLevel == 1) {
      //hsInfo hsi;
      //hsi.eDepAvg = eDepAvg;
      //_diagInfo.helixSeedData.push_back(hsi);
    }
    //if (eDepAvg > _maxEDepAvg) return;
    // compute direction of propagation and make save decision
    /*HelixTool helTool(Temp_HSeed, _tracker);
    float tzSlope = 0.0;
    float tzSlopeErr = 0.0;
    float tzSlopeChi2 = 0.0;
    helTool.dirOfProp(tzSlope, tzSlopeErr, tzSlopeChi2);
    HelixRecoDir helDir(tzSlope, tzSlopeErr, tzSlopeChi2);
    hseed._recoDir = helDir;
    hseed._propDir = helDir.predictDirection(_tzSlopeSigThresh);
    if (!validHelixDirection(hseed._propDir)) return;
    */
    // push back the helix seed to the helix seed collection
    //HSColl.emplace_back(hseed);
    }
//-----------------------------------------------------------------------------
// calling logic that needs to be called to run debug mode
//-----------------------------------------------------------------------------
    void PhiZSeedFinder::initDebugMode() {
    /*
    art::Handle<mu2e::ComboHitCollection> shcH;
    const mu2e::ComboHitCollection*       shc(nullptr);
    fEvent->getByLabel(StrawHitCollTag,shcH);
    art::InputTag sdmc_tag = StrawDigiMCCollTag;
    if (sdmc_tag == "") sdmc_tag = fSdmcCollTag;
    art::Handle<mu2e::StrawDigiMCCollection> mcdH;
    fEvent->getByLabel<mu2e::StrawDigiMCCollection>(sdmc_tag,mcdH);
    const mu2e::StrawDigiMCCollection*  mcdigis(nullptr);
    if (mcdH.isValid())   mcdigis = mcdH.product();
    */
    // first find the TC we want to focus on
    findBestTC();
    }
    //-----------------------------------------------------------------------------
    // function to find the best time cluster in debug mode given particle of interest (set in fcl)
    //-----------------------------------------------------------------------------
    //  void PhiZSeedFinder::findBestTC() {
    //_simIDsPerTC.clear();
    //GlobalConstantsHandle<ParticleDataList> pdt;
    //int _debugPdgID = 11;
    //float qSign = pdt->particle(_debugPdgID).charge();
    /*   _data._nTimeClusters = _data.tccol->size();
    //int nsh = .size();
    //loop over TCs
    for (int i=0; i<_data._nTimeClusters; i++) {
    const TimeCluster* tc = &_data.tccol->at(i);
    _finder->run(tc);
    const std::vector<StrawHitIndex>& ordchcol = tc->hits();
    int nComboHitsInTC = ordchcol.size();
    std::cout<<"nComboHitsInTC = "<<nComboHitsInTC<<std::endl;
    for (int ih=0; ih<nComboHitsInTC; ih++) {
    int ind = ordchcol[ih];
    const ComboHit* ch = &_data.chcol->at(ind);
    //ev5_HitsInNthStation hitsincluster;
    //hitsincluster.hitIndice = ind;
    //hitsincluster.hitID = ch->_sid.asUint16();
    //hitsincluster.phi = ch->phi();
    //hitsincluster.x = ch->pos().x();
    //hitsincluster.y = ch->pos().y();
    //hitsincluster.z = ch->pos().z();
    //hitsincluster.strawhits = ch->_nsh;
    //hitsincluster.station = ch.strawId().station();
    //hitsincluster.station = ch->strawId().station();
    //ComboHitsInCluster.push_back(hitsincluster);
    const mu2e::SimParticle * sim (0);
    mu2e::GenId gen_id;
    int      pdg_id(-1), mother_pdg_id(-1), generator_id(-1), sim_id(-1);
    double   mc_mom(-1.), mc_pT(-1.), mc_pZ(0.);
    const mu2e::StrawDigiMC*  sdmc = &mcdigis->at(ind);
    const mu2e::StrawGasStep* step = sdmc->earlyStrawGasStep().get();
    if (Step) {
      art::Ptr<mu2e::SimParticle> const& simptr = Step->simParticle();
      art::Ptr<mu2e::SimParticle> mother        = simptr;
      while(mother->hasParent()) mother  = mother->parent();
      sim           = mother.operator ->();
      pdg_id        = simptr->pdgId();
      mother_pdg_id = sim->pdgId();
      if (simptr->fromGenerator()) generator_id = simptr->genParticle()->generatorId().id();
      else                         generator_id = -1;
      sim_id        = simptr->id().asInt();
      mc_mom        = Step->momvec().mag();
      mc_mom_z      = Step->momvec().z();
    }
    }
    }*/
    //}
    void PhiZSeedFinder::findBestTC() {
    _simIDsPerTC.clear();
    _simInfoPerTC.clear();
    GlobalConstantsHandle<ParticleDataList> pdt;
    //int _debugPdgID = 11;
    //float qSign = pdt->particle(_debugPdgID).charge();
    std::cout<<"================================"<<std::endl;
    std::cout<<"          findBestTC            "<<std::endl;
    std::cout<<"================================"<<std::endl;
    // loop over TCs to fill _simIDsPerTC
    for (size_t i = 0; i < _data.tccol->size(); i++) {
    std::vector<mcInfo> particlesInTC;
    std::vector<mcInfoList> particlesListInTC;
    const TimeCluster* tc = &_data.tccol->at(i);
    const std::vector<StrawHitIndex>& ordchcol = tc->hits();
    int nComboHitsInTC = ordchcol.size();
    std::cout<<"nComboHitsInTC "<<nComboHitsInTC<<std::endl;
    std::cout<<"_data->_chColl2->size() = "<<_data.chcol->size()<<std::endl;
    // loop over ComboHits in a TimeCluster: fill mcSimIDs data members simID, nHits, and pdgID
    //  for(size_t j=0; j<_data.chcol->size(); ++j){
    for (size_t j = 0; j < _data.tccol->at(i)._strawHitIdxs.size(); j++) {
      int hitIndice = _data.tccol->at(i)._strawHitIdxs[j];
      std::vector<StrawDigiIndex> shids;
      _data.chcol->fillStrawDigiIndices(hitIndice, shids);
      //int ind = ordchcol[];
      //int hitIndice = ordchcol[j];
      mcInfoList mchitInfo;
      mchitInfo.nStrawHits  = (int)shids.size();
      mchitInfo.station     = _data.chcol->at(hitIndice).strawId().station();
      mchitInfo.plane       = _data.chcol->at(hitIndice).strawId().plane();
      mchitInfo.face        = _data.chcol->at(hitIndice).strawId().face();
      mchitInfo.panel       = _data.chcol->at(hitIndice).strawId().panel();
      mchitInfo.x           = _data.chcol->at(hitIndice).pos().x();
      mchitInfo.y           = _data.chcol->at(hitIndice).pos().y();
      mchitInfo.z           = _data.chcol->at(hitIndice).pos().z();
      mchitInfo.phi         = _data.chcol->at(hitIndice).pos().phi();
      mchitInfo.mom         = 0.0;
      mchitInfo.pdg         = 0;
      mchitInfo.simID       = 0;
      //loop over StrawHits(shids)
      for (size_t k = 0; k < shids.size(); k++) {
        const mu2e::SimParticle* _simParticle;
        _simParticle = _mcUtils->getSimParticle(_event, shids[k]);
        //int _pdgID = _simParticle->pdgId();
        const XYZVectorF* simMomentum = _mcUtils->getMom(_event, shids[k]);
        float simXmomentum = simMomentum->x();
        float simYmomentum = simMomentum->y();
        float simZmomentum = simMomentum->z();
        float simPerpMomentum = std::sqrt(simXmomentum * simXmomentum + simYmomentum * simYmomentum);
        float simMomentumMag = std::sqrt(simPerpMomentum * simPerpMomentum + simZmomentum * simZmomentum);
        int SimID = _mcUtils->strawHitSimId(_event, shids[k]);
        bool particleAlreadyFound = false;
        mchitInfo.mom = (double)simMomentumMag;
        mchitInfo.pdg = _simParticle->pdgId();
        mchitInfo.simID = SimID;
        for (size_t n = 0; n < particlesInTC.size(); n++) {
          if (SimID == particlesInTC[n].simID) {
          particleAlreadyFound = true;
          particlesInTC[n].nStrawHits = particlesInTC[n].nStrawHits + 1;
          if (simMomentumMag > particlesInTC[n].pMax) {
            particlesInTC[n].pMax = simMomentumMag;
          }
          if (simMomentumMag < particlesInTC[n].pMin) {
            particlesInTC[n].pMin = simMomentumMag;
          }
            break;
          }
        }
        if (particleAlreadyFound) {
          continue;
        }
        mcInfo particle;
        particle.simID = SimID;
        particle.nStrawHits = 1;
        particle.pMax = simMomentumMag;
        particle.pMin = simMomentumMag;
        particle.tcIndex = i;
        particle.mcX0 = 0.;
        particle.mcY0 = 0.;
        particle.mcRadius = 0.;
        particlesInTC.push_back(particle);
      }//end loop StrawHits(shids)
      particlesListInTC.push_back(mchitInfo);
    }//end loop ComboHits in a TC
    //preselection of MC particle
    std::vector<mcInfo> particlesInTC_update;
    for (size_t j = 0; j < particlesInTC.size(); j++) {
      //int _debugStrawHitThresh = 12;
      int _debugStrawHitThresh = 20;
      if (particlesInTC[j].nStrawHits < _debugStrawHitThresh) {
        continue;
      }
      particlesInTC_update.push_back(particlesInTC[j]);
      //std::cout<<" pnStrawHits = "<<particlesInTC[j].nStrawHits<<std::endl;
    }
    if(0 != particlesInTC_update.size()) {
      _simIDsPerTC.push_back(particlesInTC_update);
    }
    _simInfoPerTC.push_back(particlesListInTC);
    }//end loop TimeClusters
    for (size_t i = 0; i < _simInfoPerTC.size(); i++) {
      std::cout << std::string(140, '=') << std::endl;
      std::cout << " Time Cluster " << i << std::endl;
      std::cout << " Total Combo Hits: " << _simInfoPerTC[i].size() << std::endl;
      std::cout << std::string(140, '=') << std::endl;
std::cout << std::right
          << std::setw(4)  << "No."
          << std::setw(12) << "nStrawHits"
          << std::setw(10) << "station"
          << std::setw(10) << "plane"
          << std::setw(10) << "face"
          << std::setw(10) << "panel"
          << std::setw(12) << "x [mm]"
          << std::setw(12) << "y [mm]"
          << std::setw(12) << "z [mm]"
          << std::setw(12) << "Phi [rad]"
          << std::setw(12) << "P [MeV/c]"
          << std::setw(10) << "pdg"
          << std::setw(10) << "simID"
          << std::endl;
          std::cout << std::string(136, '-') << std::endl;
for (size_t j = 0; j < _simInfoPerTC[i].size(); j++) {
    std::cout << std::right
              << std::setw(4)  << j
              << std::setw(12) << _simInfoPerTC[i][j].nStrawHits
              << std::setw(10) << _simInfoPerTC[i][j].station
              << std::setw(10) << _simInfoPerTC[i][j].plane
              << std::setw(10) << _simInfoPerTC[i][j].face
              << std::setw(10) << _simInfoPerTC[i][j].panel
              << std::fixed << std::setprecision(2)
              << std::setw(12) << _simInfoPerTC[i][j].x
              << std::setw(12) << _simInfoPerTC[i][j].y
              << std::setw(12) << _simInfoPerTC[i][j].z
              << std::setw(12) << _simInfoPerTC[i][j].phi
              << std::setw(12) << _simInfoPerTC[i][j].mom
              << std::setw(10) << _simInfoPerTC[i][j].pdg
              << std::setw(10) << _simInfoPerTC[i][j].simID
              << std::endl;
}
    }
    std::cout<<"_simIDsPerTC.size()"<<_simIDsPerTC.size()<<std::endl;
    for(size_t i = 0; i < _simIDsPerTC.size(); i++) {
    for(size_t j = 0; j < _simIDsPerTC[i].size(); j++) {
    int tcIndex = _simIDsPerTC[i].at(j).tcIndex;
    int simID = _simIDsPerTC.at(i).at(j).simID;
    std::cout<<" tcIndex/SimID = "<<tcIndex<<"/"<<simID<<std::endl;
    std::cout<<" nStrawHits = "<<_simIDsPerTC[i].at(j).nStrawHits<<std::endl;
    }
    }
    std::cout<<" compute helix MC circl  "<<std::endl;
    // compute helix MC circle parameters
    for(size_t i = 0; i < _simIDsPerTC.size(); i++) {
    for(size_t j = 0; j < _simIDsPerTC[i].size(); j++) {
    int tcIndex = _simIDsPerTC[i].at(j).tcIndex;
    int simID = _simIDsPerTC[i].at(j).simID;
    std::cout<<" tcIndex/SimID = "<<tcIndex<<"/"<<simID<<std::endl;
    std::cout<<" nStrawHits = "<<_simIDsPerTC[i].at(j).nStrawHits<<std::endl;
    for (size_t k = 0; k < _data.tccol->at(tcIndex)._strawHitIdxs.size(); k++) {
    int hitIndice = _data.tccol->at(tcIndex)._strawHitIdxs[k];
    std::vector<StrawDigiIndex> shids;
    _data.chcol->fillStrawDigiIndices(hitIndice, shids);
    bool foundParticle = false;
    for (size_t l = 0; l < shids.size(); l++) {
      if (_mcUtils->strawHitSimId(_event, shids[l]) != simID) {
        continue;
      } else {
        foundParticle = true;
        const XYZVectorF* simMomentum = _mcUtils->getMom(_event, shids[l]);
        const XYZVectorF* simPosition = _mcUtils->getPos(_event, shids[l]);
        float simXmomentum = simMomentum->x();
        float simYmomentum = simMomentum->y();
        float simPerpMomentum =
          std::sqrt(simXmomentum * simXmomentum + simYmomentum * simYmomentum);
        _mcRadius = (simPerpMomentum) / (_bz0 * mmTconversion);
        _simIDsPerTC[i].at(j).mcRadius = (simPerpMomentum) / (_bz0 * mmTconversion);
        float hitX = simPosition->x();
        float hitY = simPosition->y();
        std::cout<<"hitX/hitY/simID =  "<<hitX<<"/"<<hitY<<"/"<<_mcUtils->strawHitSimId(_event, shids[l])<<std::endl;
        const mu2e::SimParticle* simParticle;
        simParticle = _mcUtils->getSimParticle(_event, shids[l]);
        int pdgID = simParticle->pdgId();
        float qSign = pdt->particle(pdgID).charge();
        _mcX0 = hitX + qSign * simYmomentum * _mcRadius / simPerpMomentum;
        _mcY0 = hitY - qSign * simXmomentum * _mcRadius / simPerpMomentum;
        std::cout<<"_mcX0/_mcY0/_mcRadius =  "<<_mcX0<<"/"<<_mcY0<<"/"<<_mcRadius<<std::endl;
        _simIDsPerTC[i].at(j).mcX0 = hitX + qSign * simYmomentum * _mcRadius / simPerpMomentum;
        _simIDsPerTC[i].at(j).mcY0 = hitY - qSign * simXmomentum * _mcRadius / simPerpMomentum;
        std::cout<<"_mcX0/_mcY0/_mcRadius =  "<<_simIDsPerTC[i].at(j).mcX0<<"/"<<_simIDsPerTC[i].at(j).mcY0<<"/"<<_simIDsPerTC[i].at(j).mcRadius<<std::endl;
        break;
      }
    }
    if (foundParticle == true) {
      break;
    }
    }
    }
    }
    for(int i=0; i<(int)_simIDsPerTC.size(); i++){
    std::cout << "===========================================" << std::endl;
    std::cout << " MC Time Cluster  : #" <<i << std::endl;
    std::cout << " MC number of particles: " << (int)_simIDsPerTC.at(i).size() << std::endl;
    std::cout << "===========================================" << std::endl;
    std::cout << std::left << std::setw(4) << "No. |";
    std::cout << std::left << std::setw(10) << "tcIndex |";
    std::cout << std::left << std::setw(10) << " nStrawHits |";
    //std::cout << std::left << std::setw(4) << " station |";
    //std::cout << std::left << std::setw(10) << " x [mm] |";
    //std::cout << std::left << std::setw(10) << " y [mm] |";
    //std::cout << std::left << std::setw(10) << " z [mm] |";
    //std::cout << std::left << std::setw(10) << " R [mm] |";
    std::cout << std::left << std::setw(10) << " pmin |";
    std::cout << std::left << std::setw(10) << " pmax |";
    std::cout << std::left << std::setw(10) << " mcX0 |";
    std::cout << std::left << std::setw(10) << " mcY0 |";
    std::cout << std::left << std::setw(10) << " mcRadius |";
    std::cout << std::left << std::setw(10) << " simID |";
    //std::cout << std::left << std::setw(10) << " PDG ";
    std::cout << std::endl;
    for(int j=0; j<(int)_simIDsPerTC.at(i).size(); j++){
    //std::cout << "j = " << j << std::endl;
    std::cout << std::right << std::setw(4)  << j;
    std::cout << std::right << std::setw(4)  << _simIDsPerTC.at(i).at(j).tcIndex;
    //std::cout << std::left << std::setw(4)  << _simIDsPerTC.at(i).at(j).station;
    std::cout << std::right << std::setw(4) << _simIDsPerTC.at(i).at(j).nStrawHits;
    std::cout << std::right << std::setw(15) << _simIDsPerTC.at(i).at(j).pMin;
    std::cout << std::right << std::setw(15) << _simIDsPerTC.at(i).at(j).pMax;
    std::cout << std::right << std::setw(15) << _simIDsPerTC.at(i).at(j).mcX0;
    std::cout << std::right << std::setw(15) << _simIDsPerTC.at(i).at(j).mcY0;
    std::cout << std::right << std::setw(15) << _simIDsPerTC.at(i).at(j).mcRadius;
    std::cout << std::right << std::setw(10) << _simIDsPerTC.at(i).at(j).simID;
    //std::cout << std::left << std::setw(10) << _simIDsPerTC.at(i).at(j).phi;
    //std::cout << std::left << std::setw(10) << _simIDsPerTC.at(i).at(j).pdg;
    std::cout << std::endl;
    }
    }
    // loop over _simIDsPerTC to find best TC
    /*_tcIndex = -1;
    int mostStrawHits = 0;
    for (size_t i = 0; i < _simIDsPerTC.size(); i++) {
    for (size_t j = 0; j < _simIDsPerTC[i].size(); j++) {
    int _debugStrawHitThresh = 12;
    if (_simIDsPerTC[i].at(j).nStrawHits < _debugStrawHitThresh) {
      continue;
    }
    float momentumDiff = _simIDsPerTC[i].at(j).pMax - _simIDsPerTC[i].at(j).pMin;
    int _debugScatterThresh = 3;
    if (momentumDiff > _debugScatterThresh) {
      continue;
    }
    if (_simIDsPerTC[i].at(j).nStrawHits > mostStrawHits) {
      mostStrawHits = _simIDsPerTC[i].at(j).nStrawHits;
      _simID = _simIDsPerTC[i].at(j).simID;
      _tcIndex = (int)i;
    }
    }
    }
    // compute helix MC circle parameters
    if (_tcIndex != -1) {
    for (size_t j = 0; j < _data.tccol->at(_tcIndex)._strawHitIdxs.size(); j++) {
    int hitIndice = _data.tccol->at(_tcIndex)._strawHitIdxs[j];
    std::vector<StrawDigiIndex> shids;
    _data.chcol->fillStrawDigiIndices(hitIndice, shids);
    bool foundParticle = false;
    for (size_t k = 0; k < shids.size(); k++) {
      if (_mcUtils->strawHitSimId(_event, shids[k]) != _simID) {
        continue;
      } else {
        foundParticle = true;
        const XYZVectorF* simMomentum = _mcUtils->getMom(_event, shids[k]);
        const XYZVectorF* simPosition = _mcUtils->getPos(_event, shids[k]);
        float simXmomentum = simMomentum->x();
        float simYmomentum = simMomentum->y();
        float simPerpMomentum =
          std::sqrt(simXmomentum * simXmomentum + simYmomentum * simYmomentum);
        _mcRadius = (simPerpMomentum) / (_bz0 * mmTconversion);
        float hitX = simPosition->x();
        float hitY = simPosition->y();
        const mu2e::SimParticle* simParticle;
        simParticle = _mcUtils->getSimParticle(_event, shids[k]);
        int pdgID = simParticle->pdgId();
        float qSign = pdt->particle(pdgID).charge();
        _mcX0 = hitX + qSign * simYmomentum * _mcRadius / simPerpMomentum;
        _mcY0 = hitY - qSign * simXmomentum * _mcRadius / simPerpMomentum;
        std::cout<<"_mcX0/_mcY0/_mcRadius =  "<<_mcX0<<"/"<<_mcY0<<"/"<<_mcRadius<<std::endl;
        break;
      }
    }
    if (foundParticle == true) {
      break;
    }
    }
    }*/
    }
    //-----------------------------------------------------------------------------
    // tcpreselection
    //-----------------------------------------------------------------------------
    int PhiZSeedFinder::tcPreSelection(int tc){
    std::cout<<"tc = "<<tc<<std::endl;
    std::cout<<"======================"<<std::endl;
    std::cout<<"======================"<<std::endl;
    std::cout<<"  _tcPreSelection     "<<std::endl;
    std::cout<<"======================"<<std::endl;
    std::cout<<"======================"<<std::endl;
    // loop over ComboHits in a TimeCluster
    int flag = 1;
    int nstrawhits = 0;
    //int ncombohits = 0;
    for (size_t i = 0; i < _data.tccol->at(tc)._strawHitIdxs.size(); i++) {
      int hitIndice = _data.tccol->at(tc)._strawHitIdxs[i];
      std::vector<StrawDigiIndex> shids;
      _data.chcol->fillStrawDigiIndices(hitIndice, shids);
      nstrawhits = nstrawhits + (int)shids.size();
    }
    //if(nstrawhits >= 15) flag = 1;
    //if(ncombohits >= 10) flag = 1;
    std::cout<<"nstrawhits = "<<nstrawhits<<std::endl;
    return flag;
    }
    //-----------------------------------------------------------------------------
    // mcpreselection
    //-----------------------------------------------------------------------------
    int PhiZSeedFinder::mcPreSelection(int tc){
    std::cout<<"tc = "<<tc<<std::endl;
    std::cout<<"======================"<<std::endl;
    std::cout<<"======================"<<std::endl;
    std::cout<<"  _mcParticleInTC     "<<std::endl;
    std::cout<<"======================"<<std::endl;
    std::cout<<"======================"<<std::endl;
    int flag = 0;
    _mcParticleInTC = 0;
    //Condition 1: MC particle in a TC (Yes/No)
    for(size_t i = 0; i < _simIDsPerTC.size(); i++) {
    for(size_t j = 0; j < _simIDsPerTC[i].size(); j++) {
    int TcIndex = _simIDsPerTC.at(i).at(j).tcIndex;
    if(TcIndex == tc)  {
      _mcParticleInTC = 1;
      flag = 1;
    }
    if(_mcParticleInTC == 1) break;
    }
    }
    std::cout<<"_mcParticleInTC = "<<_mcParticleInTC<<std::endl;
    std::cout<<"flag = "<<flag<<std::endl;
    //Condition 2: MC circle is within the tracker
    for(size_t i = 0; i<_simIDsPerTC.size(); i++) {
    for(size_t j = 0; j < _simIDsPerTC[i].size(); j++) {
    int TcIndex = _simIDsPerTC.at(i).at(j).tcIndex;
    if(TcIndex != tc) continue;
    double distance = _simIDsPerTC.at(i).at(j).mcRadius + sqrt(_simIDsPerTC.at(i).at(j).mcX0*_simIDsPerTC.at(i).at(j).mcX0 + _simIDsPerTC.at(i).at(j).mcY0*_simIDsPerTC.at(i).at(j).mcY0);
    std::cout<<"distance = "<<distance<<std::endl;
    if(distance > 680) flag = 0;
    }
    }
    std::cout<<"flag = "<<flag<<std::endl;
    return flag;
    }
    //-----------------------------------------------------------------------------
    // Function to create the time clusters
    //-----------------------------------------------------------------------------
    void PhiZSeedFinder::initTimeCluster(TimeCluster& tc){
    int nstrs = tc._strawHitIdxs.size();
    tc._nsh = 0;
    double tacc(0),tacc2(0),xacc(0),yacc(0),zacc(0),weight(0);
    for (int i=0; i<nstrs; i++) {
    int loc = tc._strawHitIdxs[i];
    const ComboHit* ch = &_data.chcol->at(loc);
    double htime = ch->correctedTime();
    double hwt = ch->nStrawHits();
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
    // Function to calculate (MC radius - Helix radius of segment)
    //-----------------------------------------------------------------------------
void PhiZSeedFinder::get_diffrad(int tc, int isegment, double& r_diff){
    //double BestFraction = 0;
    std::vector<mcDiffR> particlesInTC;
    std::cout<<"get_diffrad"<<std::endl;
    //simID of particles in a isegment
    size_t nComboHitsInSegment = _tcHits.size();
    std::cout<<"nComboHitsInSegment = "<<nComboHitsInSegment<<std::endl;
    for(size_t j=0; j<nComboHitsInSegment; j++){
      int hitIndice = _tcHits.at(j).hitIndice;
      //const TimeCluster* TC = &_data.tccol->at(tc);
      //const std::vector<StrawHitIndex>& ordchcol = TC->hits();
      //int nParticle = 0;// number of particle for the specific simID
      //for(size_t k=0; k<ordchcol.size(); ++k){//nComboHitsInTC = ordchcol.size();
      for(size_t k=0; k<_data.tccol->at(tc)._strawHitIdxs.size(); ++k){//nComboHitsInTC = ordchcol.size();
        int _hitIndice = _data.tccol->at(tc)._strawHitIdxs[k];
        std::vector<StrawDigiIndex> shids;
        _data.chcol->fillStrawDigiIndices(_hitIndice, shids);
        //int ind = ordchcol[k];
        if(hitIndice != _hitIndice) continue;
        //const ComboHit* ch = &_data.chcol->at(ind);
        //loop over StrawHits(shids)
        int simID_SegmentHit = -1;
        for (size_t l = 0; l < shids.size(); l++) {
          //const mu2e::SimParticle* _simParticle;
          //_simParticle = _mcUtils->getSimParticle(_event, shids[l]);
          //int _pdgID = _simParticle->pdgId();
          simID_SegmentHit = _mcUtils->strawHitSimId(_event, shids[l]);
        }
            bool particleAlreadyFound = false;
            for (size_t n = 0; n < particlesInTC.size(); n++) {
              if (simID_SegmentHit == particlesInTC[n].simID) {
                particleAlreadyFound = true;
                particlesInTC[n].nHits = particlesInTC[n].nHits + 1;
                break;
              }
            }
            if (particleAlreadyFound) {
              continue;
            }
            mcDiffR particle;
            particle.nHits = 1;
            particle.simID = simID_SegmentHit;
            particlesInTC.push_back(particle);
      }
    }
    //nHits and simID of MC particle
    std::cout<<"nHits and simID of MC particle "<<std::endl;
    std::vector<mcDiffR> MCparticlesInTC;
    for(size_t i=0; i< _simIDsPerTC.size(); i++){//nMCParticlePerTC = _simIDsPerTC.at(tc).size();
      for(size_t j=0; j< _simIDsPerTC[i].size(); j++){//nMCParticlePerTC = _simIDsPerTC.at(tc).size();
      int TcIndex = _simIDsPerTC.at(i).at(j).tcIndex;
      if(TcIndex != tc) continue;
      int SimID = _simIDsPerTC.at(i).at(j).simID;
      //const TimeCluster* TC = &_data.tccol->at(tc);
      //const std::vector<StrawHitIndex>& ordchcol = TC->hits();
      //int nParticle = 0;// number of particle for the specific simID
      for(size_t k=0; k<_data.tccol->at(tc)._strawHitIdxs.size(); ++k){//nComboHitsInTC = ordchcol.size();
        int _hitIndice = _data.tccol->at(tc)._strawHitIdxs[k];
        std::vector<StrawDigiIndex> shids;
        _data.chcol->fillStrawDigiIndices(_hitIndice, shids);
        //int ind = ordchcol[ih];
        //const ComboHit* ch = &_data.chcol->at(ind);
        //loop over StrawHits(shids)
        int simID = -1;
        for (size_t l = 0; l < shids.size(); l++) {
          //const mu2e::SimParticle* _simParticle;
          //_simParticle = _mcUtils->getSimParticle(_event, shids[l]);
          //int _pdgID = _simParticle->pdgId();
          simID = _mcUtils->strawHitSimId(_event, shids[l]);
        }
        if(simID != SimID) continue;
            bool particleAlreadyFound = false;
            for (size_t n = 0; n < MCparticlesInTC.size(); n++) {
              if (simID == MCparticlesInTC[n].simID) {
                particleAlreadyFound = true;
                MCparticlesInTC[n].nHits = MCparticlesInTC[n].nHits + 1;
                break;
              }
            }
            if (particleAlreadyFound) {
              continue;
            }
            mcDiffR particle;
            particle.nHits = 1;
            particle.simID = simID;
            MCparticlesInTC.push_back(particle);
      }
      }
    }
    std::cout<<"======== particlesInTC ========"<<std::endl;
    double xC = _circleFitter.x0();
    double yC = _circleFitter.y0();
    double rC = _circleFitter.radius();
    std::cout<<"xC/yC/rC = "<<xC<<"/"<<yC<<"/"<<rC<<std::endl;
    for(size_t i = 0; i < particlesInTC.size(); i++) {
      int nhit = particlesInTC.at(i).nHits;
      int simID = particlesInTC.at(i).simID;
      std::cout<<" nhit/simID = "<<nhit<<"/"<<simID<<std::endl;
    }
    std::cout<<"======== MCparticlesInTC ========"<<std::endl;
    for(size_t i = 0; i < MCparticlesInTC.size(); i++) {
      int nhit = MCparticlesInTC.at(i).nHits;
      int simID = MCparticlesInTC.at(i).simID;
      std::cout<<" nhit/simID = "<<nhit<<"/"<<simID<<std::endl;
    }
    std::cout<<"======== CalculateDiffR ========"<<std::endl;
    for(size_t i = 0; i < MCparticlesInTC.size(); i++) {
      int nhit = MCparticlesInTC.at(i).nHits;
      int simID = MCparticlesInTC.at(i).simID;
      std::cout<<" nhit/simID = "<<nhit<<"/"<<simID<<std::endl;
      int max_nhit = 0;
      int max_simID = 0;
      for(size_t j = 0; j < particlesInTC.size(); j++) {
        if(max_nhit < particlesInTC.at(j).nHits){
          max_nhit = particlesInTC.at(j).nHits;
          max_simID = particlesInTC.at(j).simID;
        }
      }
      if(max_simID != simID) continue;
      double fraction = (double)max_nhit/(double)nhit;
      std::cout<<"fraction/max_nhit/nhit/max_simID = "<<fraction<<"/"<<max_nhit<<"/"<<nhit<<"/"<<max_simID<<std::endl;
      if(fraction < 0.5) continue;
      for(size_t k=0; k<_simIDsPerTC.size(); k++) {
        for(size_t l = 0; l < _simIDsPerTC[k].size(); l++) {
          int TcIndex = _simIDsPerTC.at(k).at(l).tcIndex;
          if(TcIndex != tc) continue;
          if(simID != _simIDsPerTC.at(k).at(l).simID) continue;
        double xC = _simIDsPerTC.at(k).at(l).mcX0;
        double yC = _simIDsPerTC.at(k).at(l).mcY0;
        double MCRadius = _simIDsPerTC.at(k).at(l).mcRadius;
        double Rdiff = MCRadius - rC;
        _data.h_diffradius.push_back(Rdiff);
        r_diff = Rdiff;
        if(fabs(Rdiff) > 200){
          std::cout<<"Ottamage!!"<<std::endl;
        }
        std::cout<<"run/subrun/eventNumber = "<<run<<"/"<<subrun<<"/"<<eventNumber<<std::endl;
        std::cout<<"xC/yC/MCRadius = "<<xC<<"/"<<yC<<"/"<<MCRadius<<std::endl;
        std::cout<<"Rdiff = "<<Rdiff<<std::endl;
        //_data.h_diffradius.push_back(Rdiff);
        }
      }
      std::cout<<" fraction = "<<fraction<<std::endl;
      break;
    }
  }
//-----------------------------------------------------------------------------
//
//-----------------------------------------------------------------------------
  void PhiZSeedFinder::segment_check(int tc, int isegment){
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "        segment_check              " << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    /*for(size_t i=0; i<nComboHitsInSegment; i++){
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 0.1;//tentative value
      _circleFitter.addPoint(x, y, wP);
      _tcHits[i].used = true;
    }*/
    // check number of gap satations
    int nTotGapStation = 0;
    // sort all_BestSegmentInfo in increasing order station
    for (int i = 0; i < (int)all_BestSegmentInfo.size(); i++) {
      std::sort(all_BestSegmentInfo.at(i).begin(), all_BestSegmentInfo.at(i).end(),
              [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b) {
                  return a.station < b.station; // Ascending order by 'station'
              });
    }
    std::vector<int> gap_station_list;
    gap_station_list.clear();
    std::cout<<"isegment = "<<isegment<<std::endl;
    for(int i=0; i<(int)all_BestSegmentInfo.size(); i++){
      int min_station = 0;
    //std::cout<<"i = "<<i<<std::endl;
      for(int j=0; j<(int)all_BestSegmentInfo.at(i).size(); j++){
        if(isegment != all_BestSegmentInfo.at(i).at(j).segmentIndex) break;
        if(j == 0) min_station = all_BestSegmentInfo.at(i).at(j).station;
        int station = all_BestSegmentInfo.at(i).at(j).station;
    //std::cout<<"station = "<<station<<std::endl;
        if(min_station < station){
          int gap = station - min_station;
          if(gap >= 2) {
            int nGapStation = gap - 1;
            nTotGapStation = nTotGapStation + nGapStation;
            for(int k=1;k<=nGapStation;k++){
              int nth_station = min_station + k;
              gap_station_list.push_back(nth_station);
            }
          }
         min_station = station;
        }
        // check if there is no hit in gap station, if Yes delete remove that station from gap_station_list
        for (int k = 0; k < static_cast<int>(gap_station_list.size()); k++) {
          if(gap_station_list.at(k) == station){
            gap_station_list.erase(std::remove(gap_station_list.begin(), gap_station_list.end(), station), gap_station_list.end());
            break;
          }
        }
        //std::cout<<"all_BestSegmentInfo["<<i<<"].segmentIndex = "<<all_BestSegmentInfo.at(i).at(j).segmentIndex<<std::endl;
      }
    //std::cout<<"nTotGapStation = "<<nTotGapStation<<std::endl;
    }
    std::cout<<"nTotGapStation = "<<nTotGapStation<<std::endl;
    //Sort gap_station_list in increasing order
    sort(gap_station_list.begin(), gap_station_list.end());
    // Delete duplicate
    //Example: before {1, 2, 3, 3, 4}
    //after {1, 2, 3, 3, 4}
    for (int i = 0; i < static_cast<int>(gap_station_list.size()); i++) {
        gap_station_list.erase(std::unique(gap_station_list.begin(), gap_station_list.end()), gap_station_list.end());
     }
    for (int i = 0; i < static_cast<int>(gap_station_list.size()); i++) {
      //std::cout<<"station = "<<gap_station_list.at(i)<<std::endl;
    }
    std::cout<<"gap_station_list.size() = "<<gap_station_list.size()<<std::endl;
    if(gap_station_list.size() >= 4){
    }
    // Calculate mean and standard deviation for circular data
    double sumX = 0.0, sumY = 0.0;
    int n = _tcHits.size();
    for (int i = 0; i < n; i++) {
      sumX += TMath::Cos(_tcHits.at(i).phi);
      sumY += TMath::Sin(_tcHits.at(i).phi);
    }
    double meanAngle = TMath::ATan2(sumY / n, sumX / n);
    double sumAngularDiffSq = 0.0;
    for (int i = 0; i < n; i++) {
      double angularDiff = _tcHits.at(i).phi - meanAngle;
      angularDiff = TMath::ATan2(TMath::Sin(angularDiff), TMath::Cos(angularDiff));
      sumAngularDiffSq += angularDiff * angularDiff;
    }
    double circularVariance = sumAngularDiffSq / n;
    double circularStdDev = TMath::Sqrt(circularVariance);
    // Handle wrapping of stdDev band around -pi and pi
    double lowerBound = meanAngle - circularStdDev;
    double upperBound = meanAngle + circularStdDev;
    for (int i = 0; i < n; i++) {
      int HitIsOk = 0;
      if (lowerBound < -TMath::Pi()) {
        double phi_lower[2] = {lowerBound + 2 * TMath::Pi(), -TMath::Pi()};
        double phi_upper[2] = {TMath::Pi(), upperBound};
        //check if phi is within the STD band
        if(phi_lower[0] < _tcHits.at(i).phi and _tcHits.at(i).phi < phi_upper[0]) HitIsOk = 1;
        if(phi_lower[1] < _tcHits.at(i).phi and _tcHits.at(i).phi < phi_upper[1]) HitIsOk = 1;
      } else if (upperBound > TMath::Pi()) {
        double phi_lower[2] = {lowerBound, -TMath::Pi()};
        double phi_upper[2] = {TMath::Pi(), upperBound - 2 * TMath::Pi()};
        //check if phi is within the STD band
        if(phi_lower[0] < _tcHits.at(i).phi and _tcHits.at(i).phi < phi_upper[0]) HitIsOk = 1;
        if(phi_lower[1] < _tcHits.at(i).phi and _tcHits.at(i).phi < phi_upper[1]) HitIsOk = 1;
      } else {
        double phi_lower = lowerBound;
        double phi_upper = upperBound;
        //check if phi is within the STD band
        if(phi_lower < _tcHits.at(i).phi and _tcHits.at(i).phi < phi_upper) HitIsOk = 1;
      }
      if(HitIsOk == 1) continue;
      HitIsOk = 1;
      //check if the hit has gap station in the neighbor
      int assumption_station[2] = {_tcHits.at(i).station - 1, _tcHits.at(i).station + 1};
      //std::cout<<"assumption_station[0]/[1] = "<<assumption_station[0]<<"/"<<assumption_station[1]<<std::endl;
      for (int j = 0; j < static_cast<int>(gap_station_list.size()); j++) {
        //std::cout<<"gap_station_list.at(j) = "<<gap_station_list.at(j)<<std::endl;
        if(assumption_station[0] == gap_station_list.at(j) or assumption_station[1] == gap_station_list.at(j)) HitIsOk = 0;
      }
      if(HitIsOk == 1) continue;
      _tcHits[i].used = false;
      std::cout<<"station/phi/z = "<<_tcHits.at(i).station<<"/"<<_tcHits.at(i).phi<<"/"<<_tcHits.at(i).z<<std::endl;
    }
/*
    _circleFitter.clear();
    //Step1: fit circle of i-th segment with fixed weight
    size_t nComboHitsInSegment = _tcHits.size();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 0.1;//tentative value
      _circleFitter.addPoint(x, y, wP);
      _tcHits[i].used = true;
    }
    double xC = _circleFitter.x0();
    double yC = _circleFitter.y0();
    double rC = _circleFitter.radius();
    //Step2: fit circle of i-th segment with correct weight
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    _circleFitter.clear();
    _printweight.clear();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      std::cout<<"i = "<<i<<std::endl;
      computeCircleError2(i, xC, yC, rC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
      _tcHits[i].used = true;
    }
    //double BestFraction = 0;
    std::vector<mcDiffR> particlesInTC;
    std::cout<<"get_diffrad"<<std::endl;
    //simID of particles in a isegment
    size_t nComboHitsInSegment = _tcHits.size();
    std::cout<<"nComboHitsInSegment = "<<nComboHitsInSegment<<std::endl;
    for(size_t j=0; j<nComboHitsInSegment; j++){
      int hitIndice = _tcHits.at(j).hitIndice;
      //const TimeCluster* TC = &_data.tccol->at(tc);
      //const std::vector<StrawHitIndex>& ordchcol = TC->hits();
      //int nParticle = 0;// number of particle for the specific simID
      //for(size_t k=0; k<ordchcol.size(); ++k){//nComboHitsInTC = ordchcol.size();
      for(size_t k=0; k<_data.tccol->at(tc)._strawHitIdxs.size(); ++k){//nComboHitsInTC = ordchcol.size();
        int _hitIndice = _data.tccol->at(tc)._strawHitIdxs[k];
        std::vector<StrawDigiIndex> shids;
        _data.chcol->fillStrawDigiIndices(_hitIndice, shids);
        //int ind = ordchcol[k];
        if(hitIndice != _hitIndice) continue;
        //const ComboHit* ch = &_data.chcol->at(ind);
        //loop over StrawHits(shids)
        int simID_SegmentHit = -1;
        for (size_t l = 0; l < shids.size(); l++) {
          //const mu2e::SimParticle* _simParticle;
          //_simParticle = _mcUtils->getSimParticle(_event, shids[l]);
          //int _pdgID = _simParticle->pdgId();
          simID_SegmentHit = _mcUtils->strawHitSimId(_event, shids[l]);
        }
            bool particleAlreadyFound = false;
            for (size_t n = 0; n < particlesInTC.size(); n++) {
              if (simID_SegmentHit == particlesInTC[n].simID) {
                particleAlreadyFound = true;
                particlesInTC[n].nHits = particlesInTC[n].nHits + 1;
                break;
              }
            }
            if (particleAlreadyFound) {
              continue;
            }
            mcDiffR particle;
            particle.nHits = 1;
            particle.simID = simID_SegmentHit;
            particlesInTC.push_back(particle);
      }
    }
    //nHits and simID of MC particle
    std::cout<<"nHits and simID of MC particle "<<std::endl;
    std::vector<mcDiffR> MCparticlesInTC;
    for(size_t i=0; i< _simIDsPerTC.size(); i++){//nMCParticlePerTC = _simIDsPerTC.at(tc).size();
      for(size_t j=0; j< _simIDsPerTC[i].size(); j++){//nMCParticlePerTC = _simIDsPerTC.at(tc).size();
      int TcIndex = _simIDsPerTC.at(i).at(j).tcIndex;
      if(TcIndex != tc) continue;
      int SimID = _simIDsPerTC.at(i).at(j).simID;
      //const TimeCluster* TC = &_data.tccol->at(tc);
      //const std::vector<StrawHitIndex>& ordchcol = TC->hits();
      //int nParticle = 0;// number of particle for the specific simID
      for(size_t k=0; k<_data.tccol->at(tc)._strawHitIdxs.size(); ++k){//nComboHitsInTC = ordchcol.size();
        int _hitIndice = _data.tccol->at(tc)._strawHitIdxs[k];
        std::vector<StrawDigiIndex> shids;
        _data.chcol->fillStrawDigiIndices(_hitIndice, shids);
        //int ind = ordchcol[ih];
        //const ComboHit* ch = &_data.chcol->at(ind);
        //loop over StrawHits(shids)
        int simID = -1;
        for (size_t l = 0; l < shids.size(); l++) {
          //const mu2e::SimParticle* _simParticle;
          //_simParticle = _mcUtils->getSimParticle(_event, shids[l]);
          //int _pdgID = _simParticle->pdgId();
          simID = _mcUtils->strawHitSimId(_event, shids[l]);
        }
        if(simID != SimID) continue;
            bool particleAlreadyFound = false;
            for (size_t n = 0; n < MCparticlesInTC.size(); n++) {
              if (simID == MCparticlesInTC[n].simID) {
                particleAlreadyFound = true;
                MCparticlesInTC[n].nHits = MCparticlesInTC[n].nHits + 1;
                break;
              }
            }
            if (particleAlreadyFound) {
              continue;
            }
            mcDiffR particle;
            particle.nHits = 1;
            particle.simID = simID;
            MCparticlesInTC.push_back(particle);
      }
      }
    }
    std::cout<<"======== particlesInTC ========"<<std::endl;
    double xC = _circleFitter.x0();
    double yC = _circleFitter.y0();
    double rC = _circleFitter.radius();
    std::cout<<"xC/yC/rC = "<<xC<<"/"<<yC<<"/"<<rC<<std::endl;
    for(size_t i = 0; i < particlesInTC.size(); i++) {
      int nhit = particlesInTC.at(i).nHits;
      int simID = particlesInTC.at(i).simID;
      std::cout<<" nhit/simID = "<<nhit<<"/"<<simID<<std::endl;
    }
    std::cout<<"======== MCparticlesInTC ========"<<std::endl;
    for(size_t i = 0; i < MCparticlesInTC.size(); i++) {
      int nhit = MCparticlesInTC.at(i).nHits;
      int simID = MCparticlesInTC.at(i).simID;
      std::cout<<" nhit/simID = "<<nhit<<"/"<<simID<<std::endl;
    }
    std::cout<<"======== CalculateDiffR ========"<<std::endl;
    for(size_t i = 0; i < MCparticlesInTC.size(); i++) {
      int nhit = MCparticlesInTC.at(i).nHits;
      int simID = MCparticlesInTC.at(i).simID;
      std::cout<<" nhit/simID = "<<nhit<<"/"<<simID<<std::endl;
      int max_nhit = 0;
      int max_simID = 0;
      for(size_t j = 0; j < particlesInTC.size(); j++) {
        if(max_nhit < particlesInTC.at(j).nHits){
          max_nhit = particlesInTC.at(j).nHits;
          max_simID = particlesInTC.at(j).simID;
        }
      }
      if(max_simID != simID) continue;
      double fraction = (double)max_nhit/(double)nhit;
      std::cout<<"fraction/max_nhit/nhit/max_simID = "<<fraction<<"/"<<max_nhit<<"/"<<nhit<<"/"<<max_simID<<std::endl;
      if(fraction < 0.5) continue;
      for(size_t k=0; k<_simIDsPerTC.size(); k++) {
        for(size_t l = 0; l < _simIDsPerTC[k].size(); l++) {
          int TcIndex = _simIDsPerTC.at(k).at(l).tcIndex;
          if(TcIndex != tc) continue;
          if(simID != _simIDsPerTC.at(k).at(l).simID) continue;
        double xC = _simIDsPerTC.at(k).at(l).mcX0;
        double yC = _simIDsPerTC.at(k).at(l).mcY0;
        double MCRadius = _simIDsPerTC.at(k).at(l).mcRadius;
        double Rdiff = MCRadius - rC;
        _data.h_diffradius.push_back(Rdiff);
        if(fabs(Rdiff) > 200){
          std::cout<<"Ottamage!!"<<std::endl;
        }
        std::cout<<"run/subrun/eventNumber = "<<run<<"/"<<subrun<<"/"<<eventNumber<<std::endl;
        std::cout<<"xC/yC/MCRadius = "<<xC<<"/"<<yC<<"/"<<MCRadius<<std::endl;
        std::cout<<"Rdiff = "<<Rdiff<<std::endl;
        //_data.h_diffradius.push_back(Rdiff);
        }
      }
      std::cout<<" fraction = "<<fraction<<std::endl;
      break;
    }
*/
    //----------------------------
    // 2PiAmbiguity PhiVsZ
    //----------------------------
  }
//-----------------------------------------------------------------------
  void PhiZSeedFinder::calculateLineEquation(const Point& p1, const Point& p2, Line& line) {
        // Handle the case where the line is vertical (same x values for both points)
        if (p1.x == p2.x) {
            line.isVertical = true;
            line.verticalX = p1.x; // The line equation is x = constant
        } else {
            line.isVertical = false;
            line.slope = (p2.y - p1.y) / (p2.x - p1.x);  // y = mx + b
            line.intercept = p1.y - line.slope * p1.x;  // b = y - mx
        }
        std::cout<<"p1 x/y = "<<p1.x<<"/"<<p1.y<<std::endl;
        std::cout<<"p2 x/y = "<<p2.x<<"/"<<p2.y<<std::endl;
        std::cout<<"line slope/intercept = "<<line.slope<<"/"<<line.intercept<<std::endl;
    }
//-----------------------------------------------------------------------
// Function to find the intersection of two lines
    void PhiZSeedFinder::findIntersection(const Line& line1, const Line& line2, Point& intersection) {
        if (line1.isVertical && line2.isVertical) {
            // Two vertical lines never intersect
            intersection.x = intersection.y = std::numeric_limits<double>::infinity();
            return;
        }
        if (line1.isVertical) {
            // Line 1 is vertical, use its x to find y from line2
            intersection.x = line1.verticalX;
            intersection.y = line2.slope * intersection.x + line2.intercept;
        } else if (line2.isVertical) {
            // Line 2 is vertical, use its x to find y from line1
            intersection.x = line2.verticalX;
            intersection.y = line1.slope * intersection.x + line1.intercept;
        } else {
            // General case: both lines are not vertical
            double denom = line1.slope - line2.slope;
            if (denom == 0) {
                // Parallel lines do not intersect
                intersection.x = intersection.y = std::numeric_limits<double>::infinity();
            } else {
                intersection.x = (line2.intercept - line1.intercept) / denom;
                intersection.y = line1.slope * intersection.x + line1.intercept;
            }
        }
    }
//-----------------------------------------------------------------------
double PhiZSeedFinder::triangleArea(double x1, double y1, double x2, double y2, double x3, double y3) {
    return std::abs(x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0;
}
//-----------------------------------------------------------------------
  bool PhiZSeedFinder::isPointInTriangle(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3) {
      std::cout<<"=========================="<<std::endl;
      std::cout<<"=========================="<<std::endl;
      std::cout<<"    isPointInTriangle     "<<std::endl;
      std::cout<<"=========================="<<std::endl;
      std::cout<<"=========================="<<std::endl;
      std::cout<<"x/y = "<<x<<"/"<<y<<std::endl;
      std::cout<<"x1/y1 = "<<x1<<"/"<<y1<<std::endl;
      std::cout<<"x2/y2 = "<<x2<<"/"<<y2<<std::endl;
      std::cout<<"x3/y3 = "<<x3<<"/"<<y3<<std::endl;
    double S_ABC = triangleArea(x1, y1, x2, y2, x3, y3);
    double S_PBC = triangleArea(x, y, x2, y2, x3, y3);
    double S_PCA = triangleArea(x1, y1, x, y, x3, y3);
    double S_PAB = triangleArea(x1, y1, x2, y2, x, y);
    std::cout<<" isPointInTriangle = "<<abs(S_ABC - (S_PBC + S_PCA + S_PAB))<<std::endl;
    return std::abs(S_ABC - (S_PBC + S_PCA + S_PAB)) < 1e-9;
}
//-----------------------------------------------------------------------
  void PhiZSeedFinder::InitTrackGeometry(mu2e::RobustHelix track, std::vector<trackerData>& tracker_data){
      std::cout<<"=========================="<<std::endl;
      std::cout<<"=========================="<<std::endl;
      std::cout<<"    InitTrackGeometry     "<<std::endl;
      std::cout<<"=========================="<<std::endl;
      std::cout<<"=========================="<<std::endl;
      mu2e::GeomHandle<mu2e::Tracker> tH;
      _tracker     = tH.get();
      ChannelID cx, co;
      int npl = _tracker->nPlanes();
      for (int ipl=0; ipl< npl; ipl++) {
        const Plane* pln = &_tracker->getPlane(ipl);
        int  ist = ipl/2;
        std::vector<double> panel_x0;
        std::vector<double> panel_y0;
        std::vector<double> panel_x1;
        std::vector<double> panel_y1;
        double panel_z0 = 0.0;
        double panel_z1 = 0.0;
        std::vector<int> v_station;
        std::vector<int> v_plane;
        std::vector<int> v_face;
        std::vector<int> v_panel;
        for (unsigned ipn=0; ipn<pln->nPanels(); ipn++) {
          const Panel* panel = &pln->getPanel(ipn);
          int face;
          if (panel->id().getPanel() % 2 == 0) face = 0;
          else                                 face = 1;
          cx.Station   = ist;
          cx.Plane     = ipl % 2;
          cx.Face      = face;
          cx.Panel     = ipn;
                  int _station = static_cast<int>(ist);
                  int _Plane = static_cast<int>(ipl);
                  int _Face = static_cast<int>(face);
                  int _Panel = static_cast<int>(ipn);
          ChannelID::orderID (&cx, &co);
          //FaceZ_t* fz;
          //fz = _data.fFaceData[co.Station][co.Face];
          //Pzz_t*    pz = fz->Panel(co.Panel);
          int fID      = 3*co.Face+co.Panel;
          double wx  = panel->straw0Direction().x();
          double wy  = panel->straw0Direction().y();
          double x  = panel->origin().x();
          double y  = panel->origin().y();
          double z  = panel->origin().z();
          //double x  = panel->wirePosition.x();
          //double y  = panel->wirePosition.y();
          //double z  = panel->wirePosition.z();
          std::cout<<"=================================="<<std::endl;
          std::cout<<"Station/Plane/Face/Panel = "<<cx.Station<<"/"<<cx.Plane<<"/"<<cx.Face<<"/"<<cx.Panel<<std::endl;
          std::cout<<"x/y/z = "<<x<<"/"<<y<<"/"<<z<<std::endl;
          std::cout<<"ipl/ipn/wx/wy/fID = "<<ipl<<"/"<<ipn<<"/"<<wx<<"/"<<wy<<"/"<<fID<<std::endl;
           //Access all straws in the panel
            auto const& straws = panel->straws();
            for (size_t strawIdx = 0; strawIdx < straws.size(); ++strawIdx) {
                const auto* straw = straws[strawIdx]; // Access the straw pointer
                //auto wirePos = straw->wirePosition(); // Get the wire position (CLHEP::Hep3Vector)
                auto wirePos = straw->strawPosition(); // Get the wire position (CLHEP::Hep3Vector)
                double x = wirePos.x();
                double y = wirePos.y();
                double z = wirePos.z();
                double half_wireL = straw->halfLength();
                // Get positions of both ends
                auto calEnd = straw->strawEnd(mu2e::StrawEnd::cal);
                auto hvEnd = straw->strawEnd(mu2e::StrawEnd::hv);
                if(cx.Face == 0){
                  panel_x0.push_back(calEnd.x());
                  panel_y0.push_back(calEnd.y());
                  panel_x0.push_back(hvEnd.x());
                  panel_y0.push_back(hvEnd.y());
                  panel_z0 = hvEnd.z();
                  v_station.push_back(_station);
                  v_plane.push_back(_Plane);
                  v_face.push_back(_Face);
                  v_panel.push_back(_Panel);
                }
                if(cx.Face == 1){
                  panel_x1.push_back(calEnd.x());
                  panel_y1.push_back(calEnd.y());
                  panel_x1.push_back(hvEnd.x());
                  panel_y1.push_back(hvEnd.y());
                  panel_z1 = hvEnd.z();
                  v_station.push_back(_station);
                  v_plane.push_back(_Plane);
                  v_face.push_back(_Face);
                  v_panel.push_back(_Panel);
                }
                std::cout << "Straw = "
                          << strawIdx << std::endl;
                std::cout << "Wire Position (x, y, z) = ("
                          << x << ", " << y << ", " << z << ")" << std::endl;
                std::cout << "half_wireL = "
                          << half_wireL << std::endl;
                std::cout << "Cal End (x, y, z) = ("
                          << calEnd.x() << ", " << calEnd.y() << ", " << calEnd.z() << ")" << std::endl;
                std::cout << "HV End (x, y, z) = ("
                          << hvEnd.x() << ", " << hvEnd.y() << ", " << hvEnd.z() << ")" << std::endl;
                        break;
            }
        }
        //calculate no-coverage region where no panel covers in the tracker
        for(int k=0;k<(int)panel_x0.size(); k++){
          std::cout << "face0_panel_x = "<<panel_x0[k]<<std::endl;
          std::cout << "face0_panel_y = "<<panel_y0[k]<<std::endl;
        }
        for(int k=0;k<(int)panel_x1.size(); k++){
          std::cout << "face1_panel_x1 = "<<panel_x1[k]<<std::endl;
          std::cout << "face1_panel_y1 = "<<panel_y1[k]<<std::endl;
        }
        //face0
        // Calculate line equations
        Line line[6];
        int index = 0;
        for(int k=0;k<(int)panel_x0.size(); k++){
          Point line1 = {panel_x0[k], panel_y0[k]};
          Point line2 = {panel_x0[k+1], panel_y0[k+1]};
          if(!(k % 2 == 0)) continue;
          calculateLineEquation(line1, line2, line[index]);
          index++;
        }
        //Find intersection point
        Point intersection[6];
        findIntersection(line[0], line[1], intersection[0]);
        findIntersection(line[0], line[2], intersection[1]);
        findIntersection(line[1], line[2], intersection[2]);
        //face1
        for(int k=0;k<(int)panel_x1.size(); k++){
          Point line1 = {panel_x1[k], panel_y1[k]};
          Point line2 = {panel_x1[k+1], panel_y1[k+1]};
          if(!(k % 2 == 0)) continue;
          calculateLineEquation(line1, line2, line[index]);
          index++;
        }
        //Find intersection point
        findIntersection(line[3], line[4], intersection[3]);
        findIntersection(line[3], line[5], intersection[4]);
        findIntersection(line[4], line[5], intersection[5]);
        //expect hit in face 0
        auto position = track.position(panel_z0);
        double px = position.x(), py = position.y();
        std::cout<<"px/py/z  = "<<px<<"/"<<py<<"/"<<panel_z0<<std::endl;
        bool inside = 0;
        inside = isPointInTriangle(px, py, intersection[0].x, intersection[0].y, intersection[1].x, intersection[1].y, intersection[2].x, intersection[2].y);
        if(inside != 1){
          for(int m=0;m<(int)v_station.size();m++){
            if(v_face[m] != 0) continue;
            trackerData filltrack;
            filltrack.station = v_station[m];
            filltrack.plane = v_plane[m];
            filltrack.face = v_face[m];
            filltrack.panel = v_panel[m];
            filltrack.z = panel_z0;
            tracker_data.push_back(filltrack);
            break;
          }
        }
        if(inside == 0) std::cout<<"Face_0 NOOOOOOOOOOO = "<<panel_z0<<std::endl;
        if(inside == 1) std::cout<<"Face_0 YESSSSSSSSSS = "<<panel_z0<<std::endl;
        //expect hit:face 1
        auto _position = track.position(panel_z1);
        px = _position.x(), py = _position.y();
        inside = isPointInTriangle(px, py, intersection[3].x, intersection[3].y, intersection[4].x, intersection[4].y, intersection[5].x, intersection[5].y);
        if(inside != 1){
          for(int m=0;m<(int)v_station.size();m++){
            if(v_face[m] != 1) continue;
            trackerData filltrack;
            filltrack.station = v_station[m];
            filltrack.plane = v_plane[m];
            filltrack.face = v_face[m];
            filltrack.panel = v_panel[m];
            filltrack.z = panel_z1;
            tracker_data.push_back(filltrack);
            break;
          }
        }
        if(inside == 0) std::cout<<"Face_1 NOOOOOOOOOOO = "<<panel_z1<<std::endl;
        if(inside == 1) std::cout<<"Face_1 YESSSSSSSSSS = "<<panel_z1<<std::endl;
        for(int k=0;k<6; k++){
          std::cout << "Intersection Point: (" << intersection[k].x << ", " << intersection[k].y << ")" << std::endl;
        }
      }
  }
//-----------------------------------------------------------------------
  void PhiZSeedFinder::helix_check(int tc, int isegment){
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "        helix_check                " << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    //Step1: fit circle of i-th segment with fixed weight
    size_t nComboHitsInSegment = _tcHits.size();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 0.1;//tentative value
      _circleFitter.addPoint(x, y, wP);
      //_tcHits[i].used = true;
    }
    double xC = _circleFitter.x0();
    double yC = _circleFitter.y0();
    double rC = _circleFitter.radius();
    //Step2: fit circle of i-th segment with correct weight
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    _circleFitter.clear();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      //std::cout<<"i = "<<i<<std::endl;
      computeCircleError2(i, xC, yC, rC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
      //_tcHits[i].used = true;
    }
    //Step3: fit circle of i-th segment with correct weight
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    _circleFitter.clear();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      computeCircleError2(i, xC, yC, rC);
      computeHelixPhi(i, xC, yC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);
      //_tcHits[i].used = true;
    }
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    std::cout<<"xC/yC/rC = "<<xC<<"/"<<yC<<"/"<<rC<<std::endl;
    //return x, y of the on the circumference
    //point(x, y) -> return(st, pl, face, panel) -> 1 hit
    //number of missing hits >= 10
    //trajectory calculation to calculate the expected number o f hits
    //------------------------------------
    // count missing hits in the tracker
    //------------------------------------
    int nhits_miss = 0;
    int nhits_expect = 0;
    double fraction = 0.0; // nhits_miss/nhits_expect
    std::vector<trackerData> tracker_data;
    tracker_data.clear();
    //RobustHelix: helix trajectory, this function returns (x, y) for each z
    mu2e::RobustHelix track;
    double rcent = sqrt(xC*xC + yC*yC);
    double fcent = polyAtan2(yC, xC);
    double radius = rC;
    double lambda = 1.0/_dphidz;
    double fz0 = _fz0;
    if (fz0 > M_PI) {
      fz0 = fz0 - 2 * M_PI;
    }
    if (fz0 < -M_PI) {
      fz0 = fz0 + 2 * M_PI;
    }
    track = mu2e::RobustHelix(rcent, fcent, radius, lambda, fz0);
    std::cout<<" track.momentum() = "<<track.momentum()<<std::endl;
    std::cout<<" rcent = "<<rcent<<std::endl;
    std::cout<<" fcent = "<<fcent<<std::endl;
    std::cout<<" lambda = "<<lambda<<std::endl;
    std::cout<<" fz0 = "<<fz0<<std::endl;
    InitTrackGeometry(track, tracker_data);
    std::cout<<"======= tracker_data-print = "<<(int)tracker_data.size()<<std::endl;
    for (size_t i = 0; i < tracker_data.size(); ++i) {
      int track_station = tracker_data[i].station;
      int track_plane = tracker_data[i].plane;
      int track_face = tracker_data[i].face;
      std::cout<<"i: station/plane/face = "<<i<<" / "<<track_station<<"/"<<track_plane<<"/"<<track_face<<std::endl;
    }
    std::cout<<"====== nComboHitsInSegment-print"<<std::endl;
      for(size_t j=0; j<nComboHitsInSegment; j++){
        if(_tcHits[j].used == false) continue;
          std::cout<<"station/plane/face/panel/z = "<<_tcHits[j].station<<"/"<<_tcHits[j].plane<<"/"<<_tcHits[j].face<<"/"<<_tcHits[j].panel<<"/"<<_tcHits[j].z<<std::endl;
      }
    std::cout<<"===== nhits_expect-calculation"<<std::endl;
    for (size_t i = 0; i < tracker_data.size(); ++i) {
      nhits_expect++;
      // check if expected hits can make a hit in the panel or not
      int track_station = tracker_data[i].station;
      int track_plane = tracker_data[i].plane;
      int track_face = tracker_data[i].face;
      //int track_panel = tracker_data[i].panel;
      int hitfound = 0;
      for(size_t j=0; j<nComboHitsInSegment; j++){
        if(_tcHits[j].used == false) continue;
        if(track_station == _tcHits[j].station and track_plane == _tcHits[j].plane and track_face == _tcHits[j].face) {
          hitfound = 1;
          std::cout<<"station/plane/face/panel/z = "<<_tcHits[j].station<<"/"<<_tcHits[j].plane<<"/"<<_tcHits[j].face<<"/"<<_tcHits[j].panel<<"/"<<_tcHits[j].z<<std::endl;
        }
      }
      if(hitfound == 0) nhits_miss++;
    }
/*    for (size_t i = 0; i < tracker_data.size(); ++i) {
      float zpos = tracker_data[i].z;
      auto position = track.position(zpos);
      double rhit = sqrt(position.x()*position.x() + position.y()*position.y());
      //hits which can produce hit in a panel
      int already = 0;
      // check if expected hits can make a hit in the panel or not
      already = InitTrackGeometry(position.x(), position.y());// Yes: 1 No: 0
      if(360 < rhit and rhit < 680){
       nhits_expect++;
       already = 1;
      }
      if(already == 0) continue;
      double wP = 0.1;
      _circleFitter.addPoint(position.x(), position.y(), wP);
      std::cout <<i<<":"<< "Position at z = " << zpos << " : "<< "(x, y, z, r) = ("<< position.x() << ", "<< position.y() << ", "<< position.z() << ", "<<sqrt(position.x()*position.x() + position.y()*position.y()) <<")"<< std::endl;
      int track_station = tracker_data[i].station;
      int track_plane = tracker_data[i].plane;
      int track_face = tracker_data[i].face;
      //int track_panel = tracker_data[i].panel;
      already = 0;
      for(size_t j=0; j<nComboHitsInSegment; j++){
        if(_tcHits[j].used == false) continue;
        if(track_station == _tcHits[j].station and track_plane == _tcHits[j].plane and track_face == _tcHits[j].face) {
          already = 1;
          std::cout<<"station/plane/face/panel/z = "<<_tcHits[j].station<<"/"<<_tcHits[i].plane<<"/"<<_tcHits[j].face<<"/"<<_tcHits[j].panel<<"/"<<_tcHits[j].z<<std::endl;
        }
      }
      if(already == 0) nhits_miss++;
    }
*/
    fraction = (double)nhits_miss/(double)nhits_expect;
    std::cout<<"nhit = "<<nComboHitsInSegment<<std::endl;
    std::cout<<"nhits_miss = "<<nhits_miss<<std::endl;
    std::cout<<"nhits_expect = "<<nhits_expect<<std::endl;
    std::cout<<"fraction_helix = "<<fraction<<std::endl;
    if(fraction < 0.3) {
      std::cout<<"okashii"<<std::endl;
      std::cout<<"run/subrun/eventNumber = "<<run<<"/"<<subrun<<"/"<<eventNumber<<std::endl;
    }
    _data.h_fraction.push_back(fraction);
    //for (size_t i = 0; i < tracker_data.size(); ++i) {
      //float zpos = tracker_data[i].z;
      //auto position = track.position(zpos);
      //hits which can produce hit in a panel
      //std::cout <<i<<":"<< "Position at z = " << zpos << " : "<< "(x, y, z, r) = ("<< position.x() << ", "<< position.y() << ", "<< position.z() << ", "<<sqrt(position.x()*position.x() + position.y()*position.y()) <<")"<< std::endl;
    //}
  }
//-----------------------------------------------------------------------------
// event entry point
//-----------------------------------------------------------------------------
  void PhiZSeedFinder::produce(art::Event& Event) {
    if (_debugLevel) printf("* >>> PhiZSeedFinder::produce  event number: %10i\n",Event.event());
    //-----------------------------------------------------------------------------
    // Run-SubRun-EventNumber
    //-----------------------------------------------------------------------------
    run = Event.id().run();
    subrun = Event.id().subRun();
    eventNumber = Event.id().event();
    _event = &Event;
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
    // MC info
    if (_diagLevel == 1) {
      initDebugMode();
    }
    // prepare diagnostic tool data members
    if (_diagLevel == 1) {
      _data.h_diffradius.clear();
      _data.h_fraction.clear();
      //_data._nTimeClusters.clear();
      //MC info
      _data.h_MCnStrawHitsPerParticle.clear();
      _data.h_MCnComboHitsPerParticle.clear();
      _data.h_MCMom.clear();
      _data.h_MCTransverseMom.clear();
      _data.h_MCTanLambda.clear();
    }
    //-----------------------------------------------------------------------------
    // produce a new "HelixSeedCollection" as an output
    //-----------------------------------------------------------------------------
    std::unique_ptr<HelixSeedCollection> hsColl(new HelixSeedCollection);
    //_hsColl = hsColl.get();
    _hsColl = nullptr;
    // create output: separate by helicity
    std::map<Helicity,unique_ptr<HelixSeedCollection>> helcols;
    int counter(0);
    if (!_doSingleOutput)  {
      for( auto const& hel : _hels) {
        helcols[hel] = std::unique_ptr<HelixSeedCollection>(new HelixSeedCollection());
        //_data.nseeds [counter] = 0;
        ++counter;
      }
    }else {
      helcols[0] = std::unique_ptr<HelixSeedCollection>(new HelixSeedCollection());
      //_data.nseeds [counter] = 0;
    }
    //-----------------------------------------------------------------------------
    // run PhiZ finder to search segment
    //-----------------------------------------------------------------------------
    _data._nTimeClusters = _data.tccol->size();
    //loop over TimeClusters (TCs)
    //std::cout<<"Loop_start"<<std::endl;
    //std::cout<<"_data._nTimeClusters = "<<_data._nTimeClusters<<std::endl;
    for (int i=0; i<_data._nTimeClusters; i++) {
      const TimeCluster* tc = &_data.tccol->at(i);
      //const TimeCluster& tc = &_data.tccol->at(i);
      _data._nComboHits  = _data.chcol->size();
      _data._nStrawHits   = -1;
    //std::cout<<"_data._nComboHits = "<<_data._nComboHits<<std::endl;
    //std::cout<<"_data._nStrawHits = "<<_data._nStrawHits<<std::endl;
      _finder->run(tc);
 /*         int nseeds = _data.nSeeds();
    int count = 0;
    //std::cout << "Total seeds found: " << nseeds << "\n";
for (int i = 0; i < _data.nSeeds(); ++i) {
    PhiZSeed* seed = _data.seed(i);
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
}*/
      //--------------------------------------------
      // initialize structure for each time cluster
      //--------------------------------------------
      //HelixSeed hseed;
      //hseed._status.merge(TrkFitFlag::TPRHelix);
      //_ev_timeCluster = NULL;
       //set variables used for searching the helix candidate
      //_ev_hseed              = hseed;
      //_ev_timeCluster        = tc;
      //_ev_hseed._hhits.setParent(chcol.parent());
      //_ev_hseed._t0          = tc->_t0;
      //_ev_hseed._timeCluster = art::Ptr<TimeCluster>(*_data.tccol, i);
      //_ev_hseed._timeCluster = art::Ptr<TimeCluster>(tccH, i);
      //_ev_hseed._status.merge(TrkFitFlag::hitsOK);
      //std::cout<<"_ev_hseed._t0 = "<<_ev_hseed._t0<<std::endl;
      //std::cout << "  Number of hits: " << _ev_hseed._hhits.size() << std::endl;
      //std::cout<<"_ev_hseed_timeCluster = "<< _ev_hseed._timeCluster->nhits()<<std::endl;
      //std::cout<<"_ev_hseed_timeCluster = "<< _ev_hseed._timeCluster->_nsh<<std::endl;
      //_ev_hseed._timeCluster->_nsh = 70;
      //--------------------
      // for diagnostic
      //--------------------
      if(_diagLevel > 0) plot_PhiVsZ_OriginalTC(i);//(1)HelixPhi vs. Z, (2)Phi vs. Z, can be removed in the future
      //--------------------
      // Pre-selection on strawhits
      //--------------------
      int tcflag = 0;
      tcflag = tcPreSelection(i);
      if(tcflag == 0) continue;
      //--------------------
      // for MC
      //--------------------
      if(_diagLevel > 0){
        //int flag = 0;
        //flag = mcPreSelection(i);
        //if(flag == 0) continue;
      }
      //clear functions
      all_BestSegmentInfo.clear();
      //--------------------
      // PhiZ finder start
      //--------------------
      _HitsInCluster.clear();
      ev5_FillHitsInTimeCluster(i);
      //ev5_FillHitsInTimeCluster_ver2(i);
      //clusterInfo(i);
      double thre_residual = 0.4;// +/-0.4[rad] delta phi window for the slope
      //-----------------------------------------
      // find segment in 3 consecutive stations
      //-----------------------------------------
      std::vector<std::vector<ev5_Segment>>  diag_best_triplet_segments;
      diag_best_triplet_segments.clear();
      ev5_SegmentSearchInTriplet(diag_best_triplet_segments, thre_residual, i);
      //-----------------------------------------------------------------------------
      // make a new time cluster (maybe not needed, can be removed in the future)
      //-----------------------------------------------------------------------------
      for(int j=0; j<(int)_segmentHits.size(); j++){
      std::vector<int> hitindex;
      for(size_t k=0; k<_data.tccol->at(i)._strawHitIdxs.size(); k++){
        int hitIndice = _data.tccol->at(i)._strawHitIdxs[k];
        std::vector<StrawDigiIndex> shids;
        _data.chcol->fillStrawDigiIndices(hitIndice, shids);
        int hitID[2];
        hitID[0] = hitIndice;
        for(int l=0; l<(int)_segmentHits.at(j).size(); l++){
          hitID[1] = _segmentHits.at(j).at(l).hitID;
          if(hitID[0] == hitID[1]) hitindex.push_back(l);
        }
      }
      //std::cout<<"kitagawa_9092"<<std::endl;
      //fill to the new TimeCluster
      TimeCluster new_tc;
      //std::cout<<"(int)hitindex.size() = "<<(int)hitindex.size()<<std::endl;
      //std::cout<<"kitagawa_180"<<std::endl;
      for(int k=0; k<(int)hitindex.size(); k++){
        int ih = hitindex[k];
        //new_tc._strawHitIdxs.push_back(ordchcol[ih]);
        new_tc._strawHitIdxs.push_back(StrawHitIndex(ih));
      }
        initTimeCluster(new_tc);
        //tccol1->push_back(new_tc);
      }
      //-------------------------------
      // find Helix for each segment
      //-------------------------------
      //if((int)_segmentHits.size() == 0) continue;
      std::vector<HelixSeed>          helix_seed_vec;
      _FitRadius.clear();
      int nSegments = _segmentHits.size();
      //nSegments = 0;//kitagawa 2025/08/26
      if(nSegments == 0) continue;
      //if(nSegments != 2) continue;
      //if(_ev_hseed._timeCluster->_nsh < 112) continue;
      //if(tc->_nsh < 112) continue;
      std::cout<<"LESSSSSGOOOO"<<std::endl;
      std::cout<<"nSegments = "<<nSegments<<std::endl;
      std::cout<<"nsh = "<<tc->_nsh<<std::endl;
      for(int j=0; j<nSegments; j++){
        int hit_count = 0;
        for(int k=0; k<(int)_segmentHits.at(j).size(); k++){
          hit_count =  _segmentHits.at(j).at(k).strawhits + hit_count;
        }
        std::cout<<"Segment #/ComboHits/StrawHits = "<<j<<"/"<<(int)_segmentHits.at(j).size()<<"/"<<hit_count<<std::endl;
      }
      for(int j=0; j<nSegments; j++){
        tcHitsFill(j);
      //-----------------------------------------------
      // Pre-selection on hits of helix candidate
      //-----------------------------------------------
      std::cout<<"Pre-selection"<<std::endl;
      int thre_strawhits = 0;
      std::cout<<"segment j/combohits = "<<j<<"/"<<_segmentHits.at(j).size()<<"/"<<std::endl;
      for(int l=0; l<(int)_segmentHits.at(j).size(); l++){
        std::cout<<"strawhits = "<<_segmentHits.at(j).at(l).strawhits<<std::endl;
        thre_strawhits = thre_strawhits + _segmentHits.at(j).at(l).strawhits;
        std::cout<<"thre_strawhits = "<<thre_strawhits<<std::endl;
      }
      if(thre_strawhits < 15) continue;
        //define Helix
        //HelixSeed temp_hseed;
        //findHelix(i, j, *_hsColl);
        //findHelix(i, j, *hsColl);
        //findHelix(i, j, *hsColl, temp_hseed);
        //findHelix(i, j, *_hsColl, _ev_hseed);
        //if(_mcParticleInTC != 1) continue;
        //double r_diff = 111;
        //get_diffrad(i, j, r_diff);
        //if(fabs(r_diff) > 10.0) continue;
        //std::cout<<"hakken "<<std::endl;
        //std::cout<<"run/subrun/eventNumber = "<<run<<"/"<<subrun<<"/"<<eventNumber<<std::endl;
        //segment_check(i, j);//remove unwanted hits here from the segment
        //helix_check(i, j);//remove unwanted hits here from the segment
        //if (_diagLevel > 0) {
          //plot_PhiVsZ_forSegment(i, j);
          //plot_CirclePhiVsZ_forSegment(i, j);
          //plot_2PiAmbiguityPhiVsZ_forSegment(i, j);
          //plot_RVsZ_forSegment(i, j);
        //}
        //kitagawa
        std::cout<<"kitagawa_05 "<<std::endl;
        //for(auto const& hel : _hels ) {
          // tentatively put a copy with the specified helicity in the appropriate output vector
          //_ev_hseed._helix._helicity = hel;
      //std::cout<<"_ev_hseed_timeCluster = "<< _ev_hseed._timeCluster->nhits()<<std::endl;
      //std::cout<<"_ev_hseed._t0 = "<< _ev_hseed._t0 <<std::endl;
      //std::cout<<"_ev_hseed_helix = "<< _ev_hseed._helix._helicity <<std::endl;
          HelixSeed   _ev_hseed;
          findHelix(i, j, *_hsColl, _ev_hseed);//_hsColl is not filled, make sure if parametrs are added into _hsColl
        std::cout<<"CROSSCHECK111"<<std::endl;
        if (_diagLevel > 0) {
          //plot_PhiVsZ_forSegment(i, j);
          //plot_CirclePhiVsZ_forSegment(i, j);
          //plot_2PiAmbiguityPhiVsZ_forSegment(i, j);
          plot_2PiAmbiguityPhiVsZ_forSegment_mod(i, j);
          //plot_RVsZ_forSegment(i, j);
          plot_HelixPhiVsZ(i, j);
        }
        std::cout<<"CROSSCHECK222"<<std::endl;
        _dphidz = _lineFitter.dydx();
        _ev_hseed._helix._fz0    = _lineFitter.y0();
        if (_ev_hseed._helix._fz0 > M_PI) { _ev_hseed._helix._fz0 = _ev_hseed._helix._fz0 - 2 * M_PI; }
    if (_ev_hseed._helix._fz0 < -M_PI) { _ev_hseed._helix._fz0 = _ev_hseed._helix._fz0 + 2 * M_PI; }
        _ev_hseed._helix._lambda = 1.0/_dphidz;
        _ev_hseed._helix._chi2dZPhi = _lineFitter.chi2Dof();
          saveHelix(i, _ev_hseed);//_hsColl is not filled, make sure if parametrs are added into _hsColl
          std::cout<<"CROSSCHECK = "<<std::endl;
          size_t nComboHitsInSegment = _tcHits.size();
          int nstrawhits__ = 0;
          for (size_t k = 0; k < nComboHitsInSegment; k++) {
            if(_tcHits[k].used == false) continue;
            //std::cout<<"_tcHits["<<k<<"].hitIndice = "<<_tcHits[k].hitIndice<<std::endl;
            std::cout<<"_tcHits["<<k<<"].hitIndice = "<<_tcHits[k].hitIndice<<std::endl;
            nstrawhits__ = nstrawhits__ + _tcHits[k].strawhits;
          }
          std::cout<<"nstrawhits = "<<nstrawhits__<<std::endl;
      //-----------------------------------------------
      // Pre-selection on hits of helix candidate
      //-----------------------------------------------
      double      radius   = _ev_hseed._helix._radius;
      double      lambda   = _ev_hseed._helix._lambda;
      double      tanDip   = lambda/radius;
      // assuming B=1T
      double      mm2MeV   = 3/10.;
      double      pT       = radius*mm2MeV;
      double      p        = pT/std::cos( std::atan(tanDip));
      if(p < 50) continue;
    if (abs(lambda) <  100.) std::cout<<"0_CHECKKK__lambda<100 "<<std::endl;
    if (abs(lambda) >  1500.) std::cout<<"1_CHECKKK__lambda>1500 "<<std::endl;
    if (pT <  40.) std::cout<<"2_CHECKKK__pT<40 "<<std::endl;
    if (nstrawhits__ > 100) std::cout<<"3_CHECKKK__nhits>100 "<<std::endl;
    if (_ev_hseed._helix._chi2dXY > 5) std::cout<<"4_CHECKKK__chi2dxy "<<std::endl;
    if (_ev_hseed._helix._chi2dZPhi > 5) std::cout<<"5_CHECKKK__chi2dZPhi "<<std::endl;
    if (_ev_hseed._helix._chi2dZPhi > 5 and _ev_hseed._helix._chi2dXY > 5) std::cout<<"6_CHECKKK__chi2dXY_chi2dZPhi"<<std::endl;
    if (pT > 150) std::cout<<"7_CHECKKK__pT150 "<<std::endl;
    if (p > 500) std::cout<<"8_CHECKKK__p500 "<<std::endl;
    std::cout<<"========================================="<<std::endl;
    std::cout<<"run/subrun/eventNumber = "<<run<<"/"<<subrun<<"/"<<eventNumber<<std::endl;
    //std::cout<<"run_subrun_eventNumber_"<<run<<"_"<<subrun<<"_"<<eventNumber<<std::endl;
    std::cout<<"nTimeClusters_/_#_:"<<_data._nTimeClusters<<" / "<<i<<std::endl;
    std::cout<<"nSegments_/_#_:"<<nSegments<<" / "<<j<<std::endl;
    std::cout<<"Combohits/Strawhits_:"<<_ev_hseed._hhits.size()<<" / "<<nstrawhits__<<std::endl;
    std::cout << "timeclsuter_nChits = " << _data.tccol->at(i)._strawHitIdxs.size() << std::endl;
    std::cout << "timeclsuter_nShits = " << tc->_nsh << std::endl;
    std::cout << "radius = " << radius << std::endl;
    std::cout << "lambda = " << lambda << std::endl;
    std::cout << "pT     = " << pT << std::endl;
    std::cout << "p      = " << p  << std::endl;
    std::cout << "chi2dXY = " << _ev_hseed._helix._chi2dXY << std::endl;
    std::cout << "chi2dZPhi = " << _ev_hseed._helix._chi2dZPhi << std::endl;
    std::cout<<"========================================="<<std::endl;
        std::cout<<"kitagawa_06 "<<std::endl;
          //_ev_hseed._helix._helicity = hel;
            std::cout<<"kitagawa_07 "<<std::endl;
          //if (_ev_hseed.status().hasAnyProperty(_saveflag)){
            std::cout<<"kitagawa_08 "<<std::endl;
            helix_seed_vec.push_back(_ev_hseed);
            std::cout<<"kitagawa_09 "<<std::endl;
      std::cout<<"_ev_hseed_radius = "<< _ev_hseed._helix._radius<<std::endl;
      std::cout<<"_ev_hseed_timeCluster = "<< _ev_hseed._timeCluster->nhits()<<std::endl;
      //std::cout<<"_ev_hseed._t0 = "<< _ev_hseed._t0 <<std::endl;
      //std::cout<<"_ev_hseed_helix = "<< _ev_hseed._helix._helicity <<std::endl;
              //kitagawa
            for (int m = 0; m < (int) helix_seed_vec.size(); ++m) {
    const HelixSeed& seed = helix_seed_vec.at(m); // Get the current HelixSeed
    std::cout << "HelixSeed [" << m << "]:" << std::endl;
    // Print time cluster information
    if (seed._timeCluster.isNonnull()) {
        std::cout << "  TimeCluster nhits: " << seed._timeCluster->nhits() << std::endl;
    } else {
        std::cout << "  TimeCluster: NULL PTR" << std::endl;
    }
    if (seed._timeCluster.isNonnull()) {
    std::cout << "  TimeCluster nhits: " << seed._timeCluster->nhits() << std::endl;
    // Try printing t0 using available methods
    std::cout << "  t0: " << seed._timeCluster->t0().t0() << std::endl;
    std::cout << "  t0 error: " << seed._timeCluster->t0().t0Err() << std::endl;
    // Print basic properties
    std::cout << "  Status: " << seed._status << std::endl;
    // Print helix parameters
    //std::cout << "  Helicity: " << seed._helix._helicity << std::endl;
    std::cout << "  Radius: " << seed._helix._radius << std::endl;
    std::cout << "  chi2ndf: " << seed._helix._chi2dXY<< std::endl;
    std::cout << "  chi2dZPhi: " << seed._helix._chi2dZPhi<< std::endl;
    //std::cout << "  Phi0: " << seed._helix._phi0 << std::endl;
    //std::cout << "  Lambda: " << seed._helix._lambda << std::endl;
    // Print hit information
    std::cout << "  Number of hits: " << seed._hhits.size() << std::endl;
} else {
    std::cout << "  TimeCluster: NULL PTR" << std::endl;
}
    std::cout << "--------------------------------------" << std::endl;
}
          //}
        //}//end loop over the helicity
        //int    index_best(-1);
    /*    int    index_best(0);
        //pickBestHelix(helix_seed_vec, index_best);
        std::cout<<"kitagawa_10 "<<std::endl;
        if ( (index_best>=0) && (index_best < 2) ){
          Helicity              hel_best = helix_seed_vec[index_best]._helix._helicity;
        std::cout<<"kitagawa_11 "<<std::endl;
          //Helicity              hel_best = temp_hseed._helix._helicity;
          if (_doSingleOutput) {
            hel_best = 0;
          }
          std::cout<<"hel_best = "<<hel_best<<std::endl;//kitagawa
          HelixSeedCollection*  hcol     = helcols[hel_best].get();
          hcol->push_back(helix_seed_vec[index_best]);
          //_hsColl
          //_hsColl     = helcols[hel_best].get();
          std::cout<<"kitagawa_111 "<<std::endl;
          //_hsColl->push_back(helix_seed_vec[index_best]);
          _hsColl = hcol;
          std::cout << "_hsColl_size = " << _hsColl->size() << std::endl;
          //helcols
          //helcols[hel_best]->push_back(helix_seed_vec[index_best]);
          std::cout<<"kitagawa_12 "<<std::endl;
          std::cout << "_hsColl_size = " << _hsColl->size() << std::endl;
          //std::cout << "_helcols = " << helcols->size() << std::endl;
       std::cout<<"kitagawa_12"<<std::endl;
    if (_hsColl && !_hsColl->empty()) {
    for (int m = 0; m < (int)_hsColl->size(); ++m) {
        if (_hsColl->at(m)._timeCluster.isNonnull()) {
            std::cout << "_hsColl_nhits[" << m << "] = " << _hsColl->at(m)._timeCluster->nhits() << std::endl;
        } else {
            std::cout << "_hsColl_nhits[" << m << "] = NULL PTR" << std::endl;
        }
    }
} else {
    std::cout << "_hsColl is nullptr or empty" << std::endl;
}
       std::cout<<"kitagawa_1"<<std::endl;
          //hcol->push_back(temp_hseed);
        } else if (index_best == 2){//both helices need to be saved
        std::cout<<"kitagawa_12_5 "<<std::endl;
*/
          /*for (unsigned k=0; k<_hels.size(); ++k){
            Helicity              hel   = helix_seed_vec[k]._helix._helicity;
            if (_doSingleOutput) {
              hel = 0;
            }
            HelixSeedCollection*  hcol  = helcols[hel].get();
            hcol->push_back(helix_seed_vec[k]);
          }*/
        //}
        std::cout<<"kitagawa_13 "<<std::endl;
        //std::cout << "_hsColl_size = " << _hsColl->size() << std::endl;
        //break;//kitagawa
      }// end loop over nSegments
      std::cout<<"kitagawa_14"<<std::endl;
      int    index_best(0);
      if ( (index_best>=0) && (index_best < 2) ){
      std::cout<<"(int)helix_seed_vec.size() = "<<(int)helix_seed_vec.size()<<std::endl;
      for(int s=0;s<(int)helix_seed_vec.size();s++){
        index_best = s;
        Helicity hel_best = helix_seed_vec[index_best]._helix._helicity;
          //Helicity              hel_best = temp_hseed._helix._helicity;
          if (_doSingleOutput) {
            hel_best = 0;
          }
          //std::cout<<"hel_best = "<<hel_best<<std::endl;//kitagawa
          HelixSeedCollection*  hcol     = helcols[hel_best].get();
          hcol->push_back(helix_seed_vec[index_best]);
          //hsColl->push_back(helix_seed_vec[index_best]);
          //_hsColl = hsColl;
      }
      }
      //Fill
      /*if(nSegments == 0) continue;
      double rmax = 0.0;
      for(int j=0; j<(int)_FitRadius.size(); j++){
        if(_FitRadius[j] > rmax) rmax = _FitRadius[j];
      }
      double diff = _mcRadius - rmax;
      if(_mcRadius > 10.0) _data.h_diffradius.push_back(diff);
      */
      //_data.h_diffradius.push_back(50.0);
      //break;
    }//end loop over TimeClusters (TCs)
    //-----------------------------------------------------------------------------
    // set diagnostic tool data members
    //-----------------------------------------------------------------------------
    if (_diagLevel > 0) {
        //clear
        //_data._tccolnew->clear();
        //fill
        //_data._tccolnew = tccol1.get();
        //_data._hscolnew = hsColl.get(); // helix histo table
        std::cout<<"kitagawa_"<<std::endl;
        _hmanager->fillHistograms(&_data);
        std::cout<<"kitagawa_"<<std::endl;
    }
    //-----------------------------------------------------------------------------
    // put helix seed collection into the event record
    //-----------------------------------------------------------------------------
    std::cout<<"kitagawa_15"<<std::endl;
    //std::cout<<"hsColl = "<<hsColl->size()<<std::endl;
    //Event.put(std::move(tccol1));
    //Event.put(std::move(hsColl));
    // Print out all helicity values before moving hsColl
    //    std::cout<<"hsColl/hel = "<<"/"<<hel<<std::endl;
    //    Event.put(std::move(hsColl),Helicity::name(hel));
   // }
       // Print out all helicity values before moving hsColl
    for (auto const& hel : _hels) {
        std::cout << "Helicity: " << Helicity::name(hel) << std::endl;
    }
       std::cout<<"kitagawa_16"<<std::endl;
    if (_hsColl && !_hsColl->empty()) {
    for (int j = 0; j < (int)_hsColl->size(); ++j) {
        if (_hsColl->at(j)._timeCluster.isNonnull()) {
            std::cout << "_hsColl_timeclsuter_nhits[" << j << "] = " << _hsColl->at(j)._timeCluster->nhits() << std::endl;
            std::cout << "_hsColl_hits[" << j << "] = " << _hsColl->at(j)._hhits.size() << std::endl;
            if((int)_hsColl->at(j)._hhits.size() < 10 ) std::cout<<"akanyaro"<<std::endl;
            std::cout << "_hsColl_radius[" << j << "] = " << _hsColl->at(j)._helix._radius << std::endl;
            std::cout << "_hsColl_chi2dXY[" << j << "] = " << _hsColl->at(j)._helix._chi2dXY << std::endl;
            std::cout << "_hsColl_chi2dXY[" << j << "] = " << _hsColl->at(j)._helix._chi2dZPhi << std::endl;
            //std::cout << "_hsColl_mom[" << j << "] = " << _hsColl->at(j)._helix._momentum < std::endl;
        } else {
            std::cout << "_hsColl_nhits[" << j << "] = NULL PTR" << std::endl;
        }
    }
    } else {
      std::cout << "_hsColl is nullptr or empty" << std::endl;
    }
       std::cout<<"kitagawa_17"<<std::endl;
// Geometry handle once (not inside loops)
// hoist geometry handle once (outside hot loops)
mu2e::GeomHandle<mu2e::Tracker> tH;
const mu2e::Tracker* tracker = tH.get();
  const mu2e::HelixSeedCollection* hcol0 = helcols[0].get();
  printf("=========== StrawHits from helcols[0] (pl pn ly st  time  energyDep  dt  TOT  x  y  z) ==========\n");
  printf("#Seeds=%zu\n", hcol0->size());
for (size_t i = 0; i < hcol0->size(); ++i) {
  const mu2e::HelixSeed& hseed = hcol0->at(i);
  const mu2e::ComboHitCollection& hits = hseed.hits();
  printf("\n[Seed %zu] comboHits=%zu  (t0=%.3f)\n",
         i, hits.size(), hseed.t0().t0());
  printf("No Stn Pln Pnl Str     time     energyDep       dt       TOT          x          y          z\n");
  int counthits = 0;
  for (size_t j = 0; j < hits.size(); ++j) {
    const mu2e::ComboHit& ch = hits.at(j);
    const int nsh = ch.nStrawHits();
    for (int ish = 0; ish < nsh; ++ish) {
      // No StrawHitIndex typedef here—use auto or int:
      auto lsh = ch.index(ish);                // returns a small integer type
      size_t idx = static_cast<size_t>(lsh);   // promote for bounds check
      if (!_shColl || idx >= _shColl->size()) continue;
      const mu2e::StrawHit& sh = _shColl->at(idx);
      const mu2e::Straw&    straw = tracker->getStraw(sh.strawId());
      const auto&           pos   = straw.getMidPoint();
      printf("%2d %2d %2d %2d %2d  %8.3f   %10.4f   %8.3f  %8.3f   %10.3f %10.3f %10.3f\n",
             counthits++,
             straw.id().getStation(),
             straw.id().getPlane(),
             straw.id().getPanel(),
             //straw.id().getLayer(),
             straw.id().getStraw(),
             sh.time(),
             sh.energyDep(),
             sh.dt(),
             sh.TOT(),
             pos.x(), pos.y(), pos.z());
    }
  }
}
  printf("==================================================================================================\n");
    // Store a copy of hsColl for each helicity
    /*for (auto const& hel : _hels) {
        std::cout << "hsColl/hel = " << Helicity::name(hel) << std::endl;
        Event.put(std::make_unique<HelixSeedCollection>(*hsColl), Helicity::name(hel));
    }*/
     if (_doSingleOutput) {
      std::cout<<"_doSingleOutput"<<std::endl;
      Event.put(std::move(helcols[0]));
      //Event.put(std::move(hsColl));
    } else{
      std::cout<<"kitagawa_hel"<<std::endl;
      std::cout<<"_hels.size() = "<<_hels.size()<<std::endl;
      for(auto const& hel : _hels ) {
        //std::cout<<"hel = "<<hel<<std::endl;
        //Event.put(std::move(helcols[hel]),Helicity::name(hel));
        Event.put(std::move(helcols[hel]),Helicity::name(hel));
      }
    }
    std::cout << "kitagawa_18" << std::endl;
  }//end event entry point
}
//-----------------------------------------------------------------------------
// macro that makes this class a module.
//-----------------------------------------------------------------------------
DEFINE_ART_MODULE(mu2e::PhiZSeedFinder)
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------

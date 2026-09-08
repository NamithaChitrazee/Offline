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
#include <string>     // for std::string  (debug messages)
#include <sstream>    // for std::ostringstream (debug messages)
#include <iomanip>    // for std::setw / std::setprecision
#include <iostream>
#include <map>        // for std::map (duplicate-hit bookkeeping in ev5_fit_slope_ver4)
#include <tuple>      // for std::tuple (debug point dump in ev5_fit_slope_ver4)
#include <limits>     // for std::numeric_limits (turn-number scan in ev5_fit_slope_ver5)
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
        int nturn;
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
      // ---------------------------------------------------------
      // New Structures inspired by DeltaFinderTypes::FaceZ_t
      // ---------------------------------------------------------
      // Holds all hits belonging to one specific Face (Z-layer)
      struct PhiZFace {
        // Indices pointing to the original _HitsInCluster or _chcol vector
        std::vector<int> fHitIndices;

        // Helper to clear for next event
        void clear() { fHitIndices.clear(); }
      };

      // Holds the 4 faces (Z-layers) for one Station
      struct PhiZStation {
        // 2 Planes * 2 Faces = 4 unique Z-layers per station
        PhiZFace fFaces[4];

        void clear() {
            for(int i=0; i<4; ++i) fFaces[i].clear();
        }
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
      // Result of evaluating one candidate merge pair
      struct ev5_MergeCandidate {
        int    refIdx      = -1;
        int    testIdx     = -1;
        bool   valid       = false;   // passed all vetoes

        // --- ranking quantities, filled only for a mergeable pair ---
        // primary   : chi2/ndf of the merged phi-z fit after the n*2pi shift
        // secondary : chi2/ndf of the circle fit, used only to break ties
        double chi2_circle = 0.0;   // circle fit chi2 of the merged hypothesis (total)
        double chi2_phiz   = 0.0;   // phi-z fit chi2 of the merged hypothesis (total)
        double ndf_circle  = 0.0;
        double ndf_phiz    = 0.0;
        double Zalpha      = 0.0;   // pull of the slope difference (VETO-1)
        double Zphi        = 0.0;   // pull of the extrapolated phi difference (VETO-2)
        double fDropped    = 0.0;   // fraction of hits removed by the clean-up

        // --- helix parameters ---
        double xC = 0.0, yC = 0.0, rC = 0.0;
        double alpha = 0.0, alphaError = 0.0;
        double beta  = 0.0, betaError  = 0.0;
        int    deltaCorrection = 0;   // 2pi turn correction applied to the test segment

        // --- bookkeeping filled in for the debug summary table ---
        std::string rejectReason;      // empty when the pair is accepted
        int    nHitRef        = 0;     // hits in ref  after the circle clean-up
        int    nHitTest       = 0;     // hits in test after the circle clean-up
        int    nHitRefIn      = 0;     // hits in ref  before the clean-up
        int    nHitTestIn     = 0;     // hits in test before the clean-up
        int    nDropped       = 0;
        double chi2ndf_circle = 0.0;
        double chi2ndf_phiz   = 0.0;
        double alphaRef       = 0.0, alphaRefErr  = 0.0, chi2ndfRef  = 0.0;
        double alphaTest      = 0.0, alphaTestErr = 0.0, chi2ndfTest = 0.0;
        double alphaDiff      = 0.0, alphaThr     = 0.0;
        double dPhi           = 0.0, dPhiThr      = 0.0;
        bool   okAlpha        = false;
        bool   okPhi          = false;
        bool   okFit          = false;

        // --- the merge product itself (used directly when committing) ---
        std::vector<ev5_HitsInNthStation> mergedHits;
        std::vector<ev5_Segment>          mergedDiag;
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
    void organizeHitsByFace();
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
    void findchisq_ver2(const std::vector<ev5_HitsInNthStation>& segment);
    void ev5_select_best_segments_step_07(std::vector<std::vector<ev5_HitsInNthStation>>& all_BestSegmentInfo, std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi, int& NumberOfSegments);
 void mergeSegmentsAll(
    std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment,
    std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag,
    double thre_residual);
    // Side-effect free evaluation of one merge pair.
    // Returns true if the pair passes all vetoes and is a merge candidate.
    bool evaluateMergePair(const std::vector<ev5_HitsInNthStation>& segRef,
                           const std::vector<ev5_Segment>&          diagRef,
                           const std::vector<ev5_HitsInNthStation>& segTest,
                           const std::vector<ev5_Segment>&          diagTest,
                           int refIdx, int testIdx,
                           ev5_MergeCandidate& cand);
    void ev5_select_best_segments_step_08(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& all_ThisIsBestSegment_Diag, double threshold_deltaphi, int& NumberOfSegments);
    void ev5_fit_slope(const std::vector<ev5_Segment>& hit_diag, double& alpha, double& beta, double& chindf);
    void ev5_fit_slope_ver2(int i, double& alpha, double& beta, double& chindf);
    void ev5_fit_slope_ver3(int i, double& alpha, double& beta, double& chindf);
    void ev5_fit_slope_ver4(int index, double& alpha, double& alphaError, double& beta, double& betaError, double& chindf);
    void ev5_fit_slope_ver5(int index, double& alpha, double& alphaError, double& beta, double& betaError, double& chindf);
//-----------------------------------------------------------------------------
// helper functions for helix finder
//-----------------------------------------------------------------------------
    //void findHelix(int tc, int isegment, HelixSeedCollection& HSColl);
    void findHelix(int tc, int isegment, HelixSeedCollection& HSColl, HelixSeed& Temp_HSeed);
    void findHelix_ver2(int tc, int isegment, HelixSeedCollection& HSColl, HelixSeed& Temp_HSeed);
    void findHelix_ver3(int tc, int isegment);
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
    void plot_PhiVsZ_forSegment_ver3(int tc, int isegment);
    void plot_PhiVsZ_forEachStep(std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegmen, const char* filename, int tc, int loopIndex);
    void plot_PhiVsZ_alignment_step(int tc, int isegment, int stepIdx, int currentSegIdx, double slope, double intercept);
    void plot_PhiVsZ_RawStep(const std::vector<std::vector<ev5_HitsInNthStation>>& segments,
                         const std::string& stepName,
                         int tc,
                         int station);
    void plot_CirclePhiVsZ_forSegment(int tc, int isegment);
    void plot_2PiAmbiguityPhiVsZ_forSegment(int tc, int isegment);
    void plot_2PiAmbiguityPhiVsZ_forSegment_mod(int tc, int isegment);
    void plot_RVsZ_forSegment(int tc, int isegment);
    void plot_PhiVsZ_forSegment_debug(int ith_segment, int jth_segment,double alpha,double beta,double Chi2NDF,int plotID,const std::string& stage, bool merged);
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
    std::vector<std::vector<ev5_HitsInNthStation>> _SSegmentHits;//ComboHits and some helper variables for segments
    std::vector<cHit> _tcHits;
    ::LsqSums4 _circleFitter;
    // thresholds for mergeSegmentsAll (can be promoted to fhicl params)
    double _mergeMaxChi2NDF  = 5.0;    // target chi2/ndf for the circle clean-up
    double _mergeMaxZalpha   = 5.0;    // upper limit on the slope pull (VETO-1)
    double _mergeMaxZphi     = 5.0;    // upper limit on the phase pull (VETO-2)
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
    // The main container: 18 Stations
    PhiZStation _stationData[18];
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
void PhiZSeedFinder::organizeHitsByFace() {
    // 1. Clear previous event data
    for(int s=0; s<18; ++s) _stationData[s].clear();

    // 2. Loop over your collected hits (Assuming _HitsInCluster is already filled)
    for (size_t i = 0; i < _HitsInCluster.size(); ++i) {
        const auto& hit = _HitsInCluster[i];

        // 3. Calculate the unique Face Index (0 to 3) within the station
        // Logic: Plane (0 or 1) * 2 + Face (0 or 1)
        // Plane 0 Face 0 -> Index 0
        // Plane 0 Face 1 -> Index 1
        // Plane 1 Face 0 -> Index 2
        // Plane 1 Face 1 -> Index 3

        // Note: hit.plane is usually 0-35 globally, or 0-1 locally.
        // Ensure you use the local plane index (0 or 1).
        int localPlane = hit.plane % 2;
        int faceIndex  = (localPlane * 2) + hit.face;

        // 4. Store the index of this hit in the structured container
        if (hit.station < 18 && faceIndex < 4) {
            _stationData[hit.station].fFaces[faceIndex].fHitIndices.push_back(i);
        }
    }
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

      // 1. Organize hits into the [Station][Face] grid ONCE
      organizeHitsByFace();


      // number of stations
      // search slope in 3 consecutive stations
      // Loop from 0 to 16 station
      //int nstation = 18;
      //-------------------------------
      // start loop on stations
      //-------------------------------
      for(int n=17; n>1; n--){
         //if(n==6) break;
          //----------------------------------------------------
          // Find triplet
          //----------------------------------------------------
          // Take combo hits in n-th, (n+1)-th, (n+2)-th stations
          /*std::vector<ev5_HitsInNthStation> Hits_In_Station[3];
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
        *///if(!(nHitInStation >= 5)) continue;
        // (1st, 2nd, 3rd) = (CH>=1, CH>=1, CH>=1)
        //if(!((int)nCHsInStn_1 >= 1)) continue;
        //if(!((int)nCHsInStn_2 >= 1)) continue;
        //if(!((int)nCHsInStn_3 >= 1)) continue;

        // ----------------------------------------------------
          // Find triplet (Optimized using _stationData)
          // ----------------------------------------------------

          // 1. Clear local vectors for the 3 stations in this triplet
          std::vector<ev5_HitsInNthStation> Hits_In_Station[3];

          // 2. Define which stations we are looking at: n, n-1, n-2
          int targetStations[3] = {n, n - 1, n - 2};

          // 3. Loop over the 3 target stations (0=n, 1=n-1, 2=n-2)
          for (int k = 0; k < 3; ++k) {
              int stnIdx = targetStations[k];

              // Safety check: ensure station index is valid (0 to 17)
              if (stnIdx < 0 || stnIdx >= 18) continue;

              // 4. Loop over the 4 Faces in this station (Plane 0/1, Face 0/1)
              for (int f = 0; f < 4; ++f) {
                  // Access the pre-sorted list of indices for this face
                  const auto& faceHits = _stationData[stnIdx].fFaces[f].fHitIndices;

                  // 5. Retrieve the actual hit objects using the indices
                  for (int hitIndex : faceHits) {

                      // Access the master list directly by index
                      // NOTE: We make a COPY here so we can add it to the local list.
                      // If you want to modify 'used' flags later, you rely on the index match.
                      ev5_HitsInNthStation hitObj = _HitsInCluster.at(hitIndex);

                      // CRITICAL: Skip if this hit is already used by a previous segment
                      if (hitObj.used) continue;

                      // Add to the local list for this station
                      Hits_In_Station[k].push_back(hitObj);
                  }
              }
          }

          // ----------------------------------------------------
          // End of optimized collection
          // ----------------------------------------------------

          // [The rest of your code follows unchanged...]
          // 2 >= ComboHits [/station]
          size_t nCHsInStn_1 = Hits_In_Station[0].size();
          size_t nCHsInStn_2 = Hits_In_Station[1].size();
          size_t nCHsInStn_3 = Hits_In_Station[2].size();
          int nHitInStation = (int)nCHsInStn_1 + (int)nCHsInStn_2 + (int)nCHsInStn_3;
          // ...
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
        // ... (Outer loops j and k remain the same) ...
        if(nHitInStation >= 5 and (int)nCHsInStn_1 >= 1 and (int)nCHsInStn_2 >= 1 and (int)nCHsInStn_3 >= 1){
          for(int j=0; j<(int)nCHsInStn_1; j++){
            for(int k=0; k<(int)nCHsInStn_3; k++){

              // -----------------------------------------------------------
              // 1. Setup the Seed (Slope Calculation)
              // -----------------------------------------------------------
              int flag_hit[3] = {0};
              double first_Hit[4] = {Hits_In_Station[0].at(j).x, Hits_In_Station[0].at(j).y, Hits_In_Station[0].at(j).z, Hits_In_Station[0].at(j).phi};
              double last_Hit[4]  = {Hits_In_Station[2].at(k).x, Hits_In_Station[2].at(k).y, Hits_In_Station[2].at(k).z, Hits_In_Station[2].at(k).phi};

              double DeltaPhi = ev5_DeltaPhi(first_Hit[0], first_Hit[1], last_Hit[0], last_Hit[1]);
              double ref_sign = ev5_ParticleDirection(first_Hit[0], first_Hit[1], last_Hit[0], last_Hit[1]);

              if(DeltaPhi > 2.0) continue;

              double phi[3] = {0.0, 0.0, ref_sign*DeltaPhi};
              double z[3]   = {first_Hit[2], 0.0, last_Hit[2]};

              // Linear Fit for Slope (Alpha/Beta)
              double mean_phi = (phi[0] + phi[2])/2.0;
              double mean_z   = (z[0] + z[2])/2.0;
              double m_n[3] = {0.0};
              double m_d[3] = {0.0};
              m_n[0] = (z[0] - mean_z)*(phi[0] - mean_phi);
              m_n[2] = (z[2] - mean_z)*(phi[2] - mean_phi);
              m_d[0] = pow(z[0] - mean_z, 2);
              m_d[2] = pow(z[2] - mean_z, 2);
              double slope_alpha = (m_n[0] + m_n[2])/(m_d[0] + m_d[2]);
              double slope_beta  = phi[0] - slope_alpha*z[0];

              // -----------------------------------------------------------
              // 2. Prepare Candidate Containers
              // -----------------------------------------------------------
              std::vector<ev5_HitsInNthStation> hit_candidates;
              std::vector<ev5_Segment> SegmentInTripletStation;

              // Add Seed Hit 1 (Station 0)
              hit_candidates.push_back(Hits_In_Station[0].at(j));
              ev5_Segment hit_first;
              hit_first.deltaphi = 0.0;
              hit_first.z = Hits_In_Station[0].at(j).z;
              hit_first.alpha = slope_alpha;
              hit_first.beta = slope_beta;
              hit_first.station = Hits_In_Station[0].at(j).station;
              hit_first.usedforfit = true;
              SegmentInTripletStation.push_back(hit_first);
              flag_hit[0] = 1;

              // Add Seed Hit 3 (Station 2)
              hit_candidates.push_back(Hits_In_Station[2].at(k));
              ev5_Segment hit_last;
              hit_last.deltaphi = phi[2];
              hit_last.z = Hits_In_Station[2].at(k).z;
              hit_last.alpha = slope_alpha;
              hit_last.beta = slope_beta;
              hit_last.station = Hits_In_Station[2].at(k).station;
              hit_last.usedforfit = true;
              SegmentInTripletStation.push_back(hit_last);
              flag_hit[2] = 1;

              // Identify faces of the seed hits to avoid duplicates
              // Formula: FaceIndex = (Plane % 2) * 2 + Face
              int seedFace0 = (Hits_In_Station[0].at(j).plane % 2) * 2 + Hits_In_Station[0].at(j).face;
              int seedFace2 = (Hits_In_Station[2].at(k).plane % 2) * 2 + Hits_In_Station[2].at(k).face;

              // -----------------------------------------------------------
              // 3. Search Loop: ONE HIT PER FACE Logic
              // -----------------------------------------------------------
              // We need arrays to store the BEST hit index for each of the 4 faces in each station
              // Initialize with -1 (no hit found)
              int bestIdx_St0[4] = {-1, -1, -1, -1}; double bestRes_St0[4] = {999., 999., 999., 999.};
              int bestIdx_St1[4] = {-1, -1, -1, -1}; double bestRes_St1[4] = {999., 999., 999., 999.};
              int bestIdx_St2[4] = {-1, -1, -1, -1}; double bestRes_St2[4] = {999., 999., 999., 999.};

              // --- Station 0 (n) ---
              for(int l=0; l<(int)nCHsInStn_1; l++){
                if(l==j) continue; // Skip seed hit itself
                int f = (Hits_In_Station[0].at(l).plane % 2) * 2 + Hits_In_Station[0].at(l).face;
                if(f == seedFace0) continue; // Don't pick another hit from the seed's face

                double middle_Hit[3] = {Hits_In_Station[0].at(l).x, Hits_In_Station[0].at(l).y, Hits_In_Station[0].at(l).z};
                double dPhi = ev5_DeltaPhi(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1]);
                double sign = ev5_ParticleDirection(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1]);
                double this_phi = sign*dPhi;
                double residual = ev5_ResidualDeltaPhi(slope_alpha, slope_beta, first_Hit[2], 0.0, middle_Hit[2], this_phi);

                // If better than current best for this face, store it
                if(residual < thre_residual && residual < bestRes_St0[f]){
                   bestRes_St0[f] = residual;
                   bestIdx_St0[f] = l;
                }
              }

              // --- Station 1 (n-1) ---
              for(int l=0; l<(int)nCHsInStn_2; l++){
                int f = (Hits_In_Station[1].at(l).plane % 2) * 2 + Hits_In_Station[1].at(l).face;

                double middle_Hit[3] = {Hits_In_Station[1].at(l).x, Hits_In_Station[1].at(l).y, Hits_In_Station[1].at(l).z};
                double dPhi = ev5_DeltaPhi(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1]);
                double sign = ev5_ParticleDirection(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1]);
                double this_phi = sign*dPhi;
                double residual = ev5_ResidualDeltaPhi(slope_alpha, slope_beta, first_Hit[2], 0.0, middle_Hit[2], this_phi);

                if(residual < thre_residual && residual < bestRes_St1[f]){
                   bestRes_St1[f] = residual;
                   bestIdx_St1[f] = l;
                   // Mark flag true if we found at least one good hit in Station 1
                   flag_hit[1] = 1;
                }
              }

              // --- Station 2 (n-2) ---
              for(int l=0; l<(int)nCHsInStn_3; l++){
                if(l==k) continue; // Skip seed hit itself
                int f = (Hits_In_Station[2].at(l).plane % 2) * 2 + Hits_In_Station[2].at(l).face;
                if(f == seedFace2) continue; // Don't pick another hit from the seed's face

                double middle_Hit[3] = {Hits_In_Station[2].at(l).x, Hits_In_Station[2].at(l).y, Hits_In_Station[2].at(l).z};
                double dPhi = ev5_DeltaPhi(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1]);
                double sign = ev5_ParticleDirection(first_Hit[0], first_Hit[1], middle_Hit[0], middle_Hit[1]);
                double this_phi = sign*dPhi;
                double residual = ev5_ResidualDeltaPhi(slope_alpha, slope_beta, first_Hit[2], 0.0, middle_Hit[2], this_phi);

                if(residual < thre_residual && residual < bestRes_St2[f]){
                   bestRes_St2[f] = residual;
                   bestIdx_St2[f] = l;
                }
              }

              // -----------------------------------------------------------
              // 4. Fill the Vectors with the Winners
              // -----------------------------------------------------------

              // Helper lambda to add a hit (simplifies code)
              auto addHit = [&](const ev5_HitsInNthStation& h, double res) {
                  hit_candidates.push_back(h);

                  // Re-calculate phi for the ev5_Segment struct
                  double mid[3] = {h.x, h.y, h.z};
                  double dP = ev5_DeltaPhi(first_Hit[0], first_Hit[1], mid[0], mid[1]);
                  double sg = ev5_ParticleDirection(first_Hit[0], first_Hit[1], mid[0], mid[1]);

                  ev5_Segment seg;
                  seg.deltaphi = sg * dP;
                  seg.z = h.z;
                  seg.alpha = slope_alpha;
                  seg.beta = slope_beta;
                  seg.station = h.station;
                  seg.usedforfit = false;
                  SegmentInTripletStation.push_back(seg);
              };

              // Collect Station 0 Winners
              for(int f=0; f<4; ++f) {
                  if(bestIdx_St0[f] != -1) {
                      addHit(Hits_In_Station[0].at(bestIdx_St0[f]), bestRes_St0[f]);
                      flag_hit[0] = 1; // Mark strictly if we added extras, though seed is already there
                  }
              }
              // Collect Station 1 Winners
              for(int f=0; f<4; ++f) {
                  if(bestIdx_St1[f] != -1) {
                      addHit(Hits_In_Station[1].at(bestIdx_St1[f]), bestRes_St1[f]);
                  }
              }
              // Collect Station 2 Winners
              for(int f=0; f<4; ++f) {
                  if(bestIdx_St2[f] != -1) {
                      addHit(Hits_In_Station[2].at(bestIdx_St2[f]), bestRes_St2[f]);
                      flag_hit[2] = 1;
                  }
              }

              // -----------------------------------------------------------
              // 5. Final Quality Check
              // -----------------------------------------------------------
              int hitInSegment = (int)hit_candidates.size();
              if(hitInSegment < 5) continue;
              if(flag_hit[0] == 0 or flag_hit[1] == 0 or flag_hit[2] == 0) continue;

              bool segment_quality = ev5_TripletQuality(SegmentInTripletStation);
              if(segment_quality == 0) continue;

              segment_candidates.push_back(hit_candidates);
              diag_segment_candidates.push_back(SegmentInTripletStation);

            } // end loop k
          } // end loop j
        } // end if check

        // ==========================================================
        // INSERT PLOTTING HERE
        // This executes once per Station 'n', showing all candidates
        // found in the triplet (n, n-1, n-2)
        // ==========================================================
        std::cout<<"step_00_raw = "<<segment_candidates.size()<<std::endl;
        if (!segment_candidates.empty()) {
             // Pass "step_00_raw" or similar as the stage name
             // 'tc' is the TimeCluster index, 'n' is the current Station index
             plot_PhiVsZ_RawStep(segment_candidates, "step_00_raw", tc, n);
        }

 /*       if(nHitInStation >= 5 and (int)nCHsInStn_1 >= 1 and (int)nCHsInStn_2 >= 1 and (int)nCHsInStn_3 >= 1){
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
        */
        //-------------------------------------------------------------
        //        End Find segments in 3 consecutive stations
        //-------------------------------------------------------------




      //---------------------------------------------------------------------
      // Select several best candidates in 3 consecutive stations
      // Remove duplicat segments
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
        std::cout<<"step_01 = "<<ThisIsBestSegment.size()<<std::endl;
        // ========================================================
        // DEBUG PRINTOUT: Print hits inside ThisIsBestSegment
        // ========================================================
        std::cout << "\n========================================================================\n";
        std::cout << " [DEBUG] HITS IN ThisIsBestSegment (Station: " << n << ")\n";
        std::cout << " Total Segments: " << ThisIsBestSegment.size() << "\n";
        std::cout << "------------------------------------------------------------------------\n";
        std::cout << std::left
                  << std::setw(10) << "SegIdx"
                  << std::setw(10) << "HitIdx"
                  << std::setw(10) << "Station"
                  << std::setw(12) << "X [mm]"
                  << std::setw(12) << "Y [mm]"
                  << std::setw(12) << "Z [mm]"
                  << std::setw(12) << "Phi [rad]" << "\n";
        std::cout << "------------------------------------------------------------------------\n";

        // Outer loop: iterate through each segment
        for (size_t seg_idx = 0; seg_idx < ThisIsBestSegment.size(); ++seg_idx) {
            // Inner loop: iterate through the hits within this specific segment
            for (const auto& h : ThisIsBestSegment[seg_idx]) {
                double phi = std::atan2(h.y, h.x); // Calculate phi for the printout
                std::cout << std::left
                          << std::setw(10) << seg_idx
                          << std::setw(10) << h.hitIndice
                          << std::setw(10) << h.station
                          << std::setw(12) << std::fixed << std::setprecision(3) << h.x
                          << std::setw(12) << std::fixed << std::setprecision(3) << h.y
                          << std::setw(12) << std::fixed << std::setprecision(3) << h.z
                          << std::setw(12) << std::fixed << std::setprecision(4) << phi << "\n";
            }
        }
        std::cout << "========================================================================\n";
        // ========================================================
        // ========================================================
        std::cout<<"statin = "<<n<<std::endl;
        if (!ThisIsBestSegment.empty()) {
             plot_PhiVsZ_RawStep(ThisIsBestSegment, "step_01", tc, n);
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
      std::cout<<"step_03 = "<<ThisIsBestSegment.size()<<std::endl;
        if (!ThisIsBestSegment.empty()) {
             plot_PhiVsZ_RawStep(ThisIsBestSegment, "step_03", tc, n);
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
      plot_PhiVsZ_forEachStep(_segmentHits, "step_05A", tc, 999);
      std::cout<<"ev5_select_best_segments_step_05A size = "<<_segmentHits.size()<<std::endl;
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
      plot_PhiVsZ_forEachStep(_segmentHits, "step_06A", tc, 999);
      std::cout<<"ev5_select_best_segments_step_06A size = "<<_segmentHits.size()<<std::endl;
      // remove duplicate segments based on hitID. remove segments based on slope value, fraction of overlapped hits, and Chi2/NDF
      ev5_select_best_segments_step_06(all_ThisIsBestSegment, all_ThisIsBestSegment_Diag, thre_residual);
      _segmentHits.clear();
      _segmentHits = all_ThisIsBestSegment;
      plot_PhiVsZ_forEachStep(_segmentHits, "step_06", tc, 999);
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

      //ckeck if segment has hits > 5
      const size_t minHits = 5;

        // (optional) sanity check
        if (all_ThisIsBestSegment.size() != all_ThisIsBestSegment_Diag.size()) {
          std::cerr << "[WARN] size mismatch: seg=" << all_ThisIsBestSegment.size()
                    << " diag=" << all_ThisIsBestSegment_Diag.size() << "\n";
        }

        for (int j = (int)all_ThisIsBestSegment.size() - 1; j >= 0; --j) {
          if (all_ThisIsBestSegment.at(j).size() < minHits) {
            all_ThisIsBestSegment.erase(all_ThisIsBestSegment.begin() + j);
            // delete corresponding diag index (guard in case of mismatch)
            if (j < (int)all_ThisIsBestSegment_Diag.size()) {
              all_ThisIsBestSegment_Diag.erase(all_ThisIsBestSegment_Diag.begin() + j);
            }
          }
        }
        _segmentHits.clear();
              _segmentHits = all_ThisIsBestSegment;


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
      plot_PhiVsZ_forEachStep(_segmentHits, "step_07", tc, 999);
      std::cout<<"mergeSegmentsAll size = "<<_segmentHits.size()<<std::endl;


// Print header once
// Column widths (adjust as needed)
const int W1 = 12;
const int W = 14;

// Print header once
std::cout
  << "-------------------------------------------------------------------------------------------------------------\n"
  << std::left
  << std::setw(W1) << "#Segment"
  << std::setw(W1) << "k"
  << std::setw(W)  << "hitIndice"
  << std::setw(W)  << "hitID"
  << std::setw(W)  << "station"
  << std::setw(W)  << "plane"
  << std::setw(W)  << "face"
  << std::setw(W)  << "panel"
  << std::setw(W)  << "segmentIndex"
  << std::setw(W)  << "used"
  << std::setw(W)  << "nturn"
  << std::setw(W)  << "strawhits"
  << std::setw(W)  << "x"
  << std::setw(W)  << "y"
  << std::setw(W)  << "z"
  << std::setw(W)  << "phi"
  << std::setw(W)  << "phiDiag"
  << std::setw(W)  << "helixPhi"
  << std::setw(W)  << "circleError2"
  << std::setw(W)  << "helixPhiError2"
  << "\n"
  << "-------------------------------------------------------------------------------------------------------------\n";

for (size_t j = 0; j < _segmentHits.size(); ++j) {
  for (size_t k = 0; k < _segmentHits[j].size(); ++k) {

    const auto& h = _segmentHits[j][k];

    // For x, y, z: print with 2 decimal places
std::cout << std::left
    << std::setw(W1) << j
    << std::setw(W1) << k
    << std::setw(W)  << h.hitIndice
    << std::setw(W)  << h.hitID
    << std::setw(W)  << h.station
    << std::setw(W)  << h.plane
    << std::setw(W)  << h.face
    << std::setw(W)  << h.panel
    << std::setw(W)  << h.segmentIndex
    << std::setw(W)  << h.used
    << std::setw(W)  << h.nturn
    << std::setw(W)  << h.strawhits
    << std::setw(W)  << std::fixed << std::setprecision(2) << h.x
    << std::setw(W)  << std::fixed << std::setprecision(2) << h.y
    << std::setw(W)  << std::fixed << std::setprecision(2) << h.z
    << std::setw(W)  << h.phi
    << std::setw(W)  << h.phiDiag
    << std::setw(W)  << h.helixPhi
    << std::setw(W)  << h.circleError2
    << std::setw(W)  << h.helixPhiError2
    << "\n";
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

      //print



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
    std::cout << " ev5_select_best_segments_step_01     " << std::endl;
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

// Make sure to include this at the top of your file
// #include <format>
// #include <iostream>
if (_debugLevel) {
    // --- 1. Print Total Collection Size ---
    std::cout << "\n" << std::string(180, '=') << "\n";
    std::cout << std::format(" TOTAL: {} Segments found.\n", ThisIsBestSegment.size());
    std::cout << std::string(180, '=') << "\n";

    // --- 2. Loop over each Segment ---
    for (size_t i = 0; i < ThisIsBestSegment.size(); ++i) {

        // --- A. Print Summary for THIS Segment ---
        std::cout << "\n" << std::string(180, '=') << "\n";
        std::cout << std::format(" Segment Index: {} | Number of Hits: {}\n", i, ThisIsBestSegment[i].size());
        std::cout << std::string(180, '=') << "\n";

        // --- B. Print Table Header ---
        // Using {:<N} for Left Align, width N
        std::cout << std::format(
            "{:<8} {:<8} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} "
            "{:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n",
            "#Seg", "Hit#", "HitInd", "HitID", "Stn", "Pln", "Fce", "Pnl", "SegIdx", "Used", "nTurn", "Straws",
            "X", "Y", "Z", "Phi", "PhiDiag", "HelixPhi", "CircErr2", "HPhiErr2"
        );
        std::cout << std::string(180, '-') << "\n";

        // --- C. Loop over Hits in this Segment ---
        for (size_t j = 0; j < ThisIsBestSegment[i].size(); ++j) {
            const auto& h = ThisIsBestSegment[i][j];

            // .2f = 2 decimal places, .3f = 3 decimal places
            std::cout << std::format(
                "{:<8} {:<8} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} {:<10} "
                "{:<12.2f} {:<12.2f} {:<12.2f} {:<12.3f} {:<12.3f} {:<12.3f} {:<12.3f} {:<12.3f}\n",
                i, j,
                h.hitIndice, h.hitID, h.station, h.plane, h.face, h.panel,
                h.segmentIndex, h.used, h.nturn, h.strawhits,
                h.x, h.y, h.z, h.phi, h.phiDiag, h.helixPhi, h.circleError2, h.helixPhiError2
            );
        }
    }
    // Final closing line
    std::cout << std::string(180, '-') << "\n";
}

  }//end ev5_select_best_segments_step_01
//-----------------------------------------------------------------------------
  void PhiZSeedFinder::ev5_select_best_segments_step_02(int station, std::vector<std::vector<ev5_HitsInNthStation>>& ThisIsBestSegment, std::vector<std::vector<ev5_Segment>>& ThisIsBestSegment_Diag, int nCH, double threshold_deltaphi){
      std::vector<std::vector<ev5_HitsInNthStation>> segments = ThisIsBestSegment;
      std::vector<std::vector<ev5_Segment>> diag_segments = ThisIsBestSegment_Diag;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      std::cout << " ev5_select_best_segments_step_02  " << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      //std::cout << "-----------------------------------" << std::endl;
      std::cout<<"station = "<<station<<std::endl;
      //---------------------------------------------------------------------------------------------
      // Search remaining ComboHit candidates in 3 consecutive stations
      //---------------------------------------------------------------------------------------------
      int max_station = station;
      std::vector<std::vector<ev5_HitsInNthStation>> new_segments;
      std::vector<std::vector<ev5_Segment>>  new_diag_segments;
      new_segments.clear();
      new_diag_segments.clear();

      // Take segment
      for(int i = 0; i < (int)segments.size(); i++){
        std::vector<ev5_HitsInNthStation> hit_candidates;
        std::vector<ev5_Segment> hit_diag_candidates;

        for(int j = 0; j < (int)segments.at(i).size(); j++){
          hit_candidates.push_back(segments.at(i).at(j));
          hit_diag_candidates.push_back(diag_segments.at(i).at(j));
        }

        // ========================================================
        // DEBUG PRINTOUT: Print Segment and Hit Info
        // ========================================================
        std::cout << "\n========================================================================\n";
        std::cout << " [DEBUG] PROCESSING SEGMENT INDEX: " << i << "\n";
        std::cout << " Total Hits in this Segment: " << hit_candidates.size() << "\n";
        std::cout << "------------------------------------------------------------------------\n";
        std::cout << std::left
                  << std::setw(10) << "HitIdx"
                  << std::setw(10) << "Station"
                  << std::setw(12) << "X [mm]"
                  << std::setw(12) << "Y [mm]"
                  << std::setw(12) << "Z [mm]"
                  << std::setw(12) << "Phi [rad]" << "\n";
        std::cout << "------------------------------------------------------------------------\n";

        for (const auto& h : hit_candidates) {
            double phi = std::atan2(h.y, h.x); // Calculate phi for the printout
            std::cout << std::left
                      << std::setw(10) << h.hitIndice
                      << std::setw(10) << h.station
                      << std::setw(12) << std::fixed << std::setprecision(3) << h.x
                      << std::setw(12) << std::fixed << std::setprecision(3) << h.y
                      << std::setw(12) << std::fixed << std::setprecision(3) << h.z
                      << std::setw(12) << std::fixed << std::setprecision(4) << phi << "\n";
        }
        std::cout << "========================================================================\n";

        // ========================================================
        // 1. Sort hit_candidates and hit_diag_candidates by Z (Descending: max Z to min Z)
        // ========================================================
        std::vector<int> idx(hit_candidates.size());
        std::iota(idx.begin(), idx.end(), 0); // Fill with 0, 1, 2...
        std::sort(idx.begin(), idx.end(), [&](int a, int b) {
            return hit_candidates[a].z > hit_candidates[b].z;
        });

        std::vector<ev5_HitsInNthStation> sorted_hits;
        std::vector<ev5_Segment> sorted_diags;
        for (int index : idx) {
            sorted_hits.push_back(hit_candidates[index]);
            sorted_diags.push_back(hit_diag_candidates[index]);
        }
        hit_candidates = sorted_hits;
        hit_diag_candidates = sorted_diags;

        // ========================================================
        // 2. Setup Reference Hit & Calculate Linear Fit (Phi vs Z)
        // ========================================================
        size_t refIdx = 0; // After sorting, 0 is the max Z
        int reference_indice = hit_candidates[refIdx].hitIndice; // <--- ADD THIS LINE
        double refZ   = hit_candidates[refIdx].z;
        double refPhi = std::atan2(hit_candidates[refIdx].y, hit_candidates[refIdx].x);
        double defaultWeight = 1.0 / 0.01; // From your previous snippet

        // Clear the fitter for this new segment
        _lineFitter.clear();

        for (auto& h : hit_candidates) {
            double raw_phi = std::atan2(h.y, h.x);

            // Calculate absolute shortest delta phi to avoid +/- Pi wrap-around
            double dPhi = raw_phi - refPhi;
            while (dPhi > M_PI)  dPhi -= 2.0 * M_PI;
            while (dPhi < -M_PI) dPhi += 2.0 * M_PI;

            // This creates a smooth continuous Phi value for the linear fit
            double unwrapped_phi = refPhi + dPhi;

            // Add the point to your LineFitter (X = Z, Y = unwrapped_phi)
            // Note: Adjust the method name "addPoint" if your class uses something else!
            _lineFitter.addPoint(h.z, unwrapped_phi, defaultWeight);
        }

        // Execute the fit (if your class requires an explicit fit command)
        // _lineFitter.fit();

        // Extract the parameters from the fitter
        // Note: Adjust "slope()" and "intercept()" to match your class methods
        double lineSlope     = _lineFitter.dydx(); // The Slope (dphi/dz)
        double lineIntercept = _lineFitter.y0();   // The Intercept (phi at z=0)

        if (_debugLevel > 0) {
            std::cout << "\n[DEBUG] --- Segment " << i << " _lineFitter Results ---\n";
            std::cout << "Ref Hit Z: " << refZ << " | Ref Phi: " << refPhi << "\n";
            std::cout << "Slope (dPhi/dZ): " << lineSlope << " | Intercept: " << lineIntercept << "\n";
        }



        // ========================================================
        // 3. Scan remaining ComboHits to see if they fit the slope
        // ========================================================
        for(int k = max_station; k >= max_station - 2; k--) {
            std::vector<ev5_HitsInNthStation> Hits_In_Station;

            for(int j = 0; j < (int)_HitsInCluster.size(); j++) {
                if(_HitsInCluster.at(j).used == true) continue;
                if(k != _HitsInCluster.at(j).station) continue;

                // Check if already in hit_candidates
                bool flag_alreadyUsed = false;
                for (const auto &element : hit_candidates) {
                  if(element.hitID == _HitsInCluster.at(j).hitID) {
                      flag_alreadyUsed = true;
                      break;
                  }
                }
                if(flag_alreadyUsed) continue;

                // Grab coordinates
                double x = _HitsInCluster.at(j).x;
                double y = _HitsInCluster.at(j).y;
                double z = _HitsInCluster.at(j).z;

                double raw_phi = std::atan2(y, x);
                double predictedPhi = lineSlope * z + lineIntercept;

                // Residual distance: Actual Phi - Predicted Phi
                double residual_phi = raw_phi - predictedPhi;

                // CRITICAL: Handle the +/- Pi wrap boundary for the residual!
                while (residual_phi > M_PI)  residual_phi -= 2.0 * M_PI;
                while (residual_phi < -M_PI) residual_phi += 2.0 * M_PI;

                if (_debugLevel > 0) {
                    std::cout << "  -> Testing HitID " << _HitsInCluster.at(j).hitID
                              << " at Z: " << z
                              << " | RawPhi: " << raw_phi
                              << " | PredPhi: " << predictedPhi
                              << " | Residual: " << std::abs(residual_phi) << "\n";
                }

                // If within threshold, accept it!
                if(std::abs(residual_phi) < threshold_deltaphi) {

                    // 1. Create and push the standard hit
                    ev5_HitsInNthStation new_hit = _HitsInCluster.at(j);
                    hit_candidates.push_back(new_hit);

                    // 2. Create and push the diagnostic hit using your exact struct
                    ev5_Segment new_diag;
                    new_diag.deltaphi   = residual_phi;       // The delta phi we just calculated
                    new_diag.z          = new_hit.z;          // Z coordinate
                    new_diag.station    = new_hit.station;    // Station number

                    // Fill in the rest of the struct with the fit info
                    new_diag.alpha      = lineSlope;          // The slope we used
                    new_diag.beta       = lineIntercept;      // The intercept we used
                    new_diag.chiNDF     = 0.0;                // Default/placeholder
                    new_diag.reference_point = reference_indice; // From your reference hit
                    new_diag.usedforfit = false;              // Set to false since it was added AFTER the fit

                    hit_diag_candidates.push_back(new_diag);

                    if (_debugLevel > 0) {
                        std::cout << "     *** HIT ACCEPTED! ***\n";
                    }
                }
            } // end loop on cluster hits
        } // end loop on stations

        new_segments.push_back(hit_candidates);
        new_diag_segments.push_back(hit_diag_candidates);

      } // end loop on segment

      // Re-Fill
      ThisIsBestSegment.clear();
      ThisIsBestSegment_Diag.clear();
      ThisIsBestSegment = new_segments;
      ThisIsBestSegment_Diag = new_diag_segments;
        // ========================================================
        // DEBUG PRINTOUT: Print hits inside ThisIsBestSegment
        // ========================================================
        std::cout << "end ev5_select_best_segments_step_02\n";
        std::cout << "\n========================================================================\n";
        std::cout << " [DEBUG] HITS IN ThisIsBestSegment (Station: " << station << ")\n";
        std::cout << " Total Segments: " << ThisIsBestSegment.size() << "\n";
        std::cout << "------------------------------------------------------------------------\n";
        std::cout << std::left
                  << std::setw(10) << "SegIdx"
                  << std::setw(10) << "HitIdx"
                  << std::setw(10) << "Station"
                  << std::setw(12) << "X [mm]"
                  << std::setw(12) << "Y [mm]"
                  << std::setw(12) << "Z [mm]"
                  << std::setw(12) << "Phi [rad]" << "\n";
        std::cout << "------------------------------------------------------------------------\n";

        // Outer loop: iterate through each segment
        for (size_t seg_idx = 0; seg_idx < ThisIsBestSegment.size(); ++seg_idx) {
            // Inner loop: iterate through the hits within this specific segment
            for (const auto& h : ThisIsBestSegment[seg_idx]) {
                double phi = std::atan2(h.y, h.x); // Calculate phi for the printout
                std::cout << std::left
                          << std::setw(10) << seg_idx
                          << std::setw(10) << h.hitIndice
                          << std::setw(10) << h.station
                          << std::setw(12) << std::fixed << std::setprecision(3) << h.x
                          << std::setw(12) << std::fixed << std::setprecision(3) << h.y
                          << std::setw(12) << std::fixed << std::setprecision(3) << h.z
                          << std::setw(12) << std::fixed << std::setprecision(4) << phi << "\n";
            }
        }
        std::cout << "========================================================================\n";
  }//end ev5_select_best_segments_step_02
//--------------------------------------------------------------------------------//
void PhiZSeedFinder::ev5_select_best_segments_step_03(std::vector<std::vector<ev5_HitsInNthStation>> &ThisIsBestSegment, std::vector<std::vector<ev5_Segment>> &ThisIsBestSegment_Diag, int station, int nCH, double threshold_deltaphi){

    std::cout<<"================================"<<std::endl;
    std::cout<<"ev5_select_best_segments_step_03"<<std::endl;
    std::cout<<"================================"<<std::endl;
      std::cout<<"station = "<<station<<std::endl;
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
        // ========================================================
        // DEBUG PRINTOUT: Print Segment and Hit Info
        // ========================================================
        std::cout << "\n========================================================================\n";
        std::cout << " [DEBUG] PROCESSING SEGMENT INDEX: " << i << "\n";
        std::cout << " Total Hits in this Segment: " << hit_candidates.size() << "\n";
        std::cout << "------------------------------------------------------------------------\n";
        std::cout << std::left
                  << std::setw(10) << "HitIdx"
                  << std::setw(10) << "Station"
                  << std::setw(12) << "X [mm]"
                  << std::setw(12) << "Y [mm]"
                  << std::setw(12) << "Z [mm]"
                  << std::setw(12) << "Phi [rad]" << "\n";
        std::cout << "------------------------------------------------------------------------\n";

        for (const auto& h : hit_candidates) {
            double phi = std::atan2(h.y, h.x); // Calculate phi for the printout
            std::cout << std::left
                      << std::setw(10) << h.hitIndice
                      << std::setw(10) << h.station
                      << std::setw(12) << std::fixed << std::setprecision(3) << h.x
                      << std::setw(12) << std::fixed << std::setprecision(3) << h.y
                      << std::setw(12) << std::fixed << std::setprecision(3) << h.z
                      << std::setw(12) << std::fixed << std::setprecision(4) << phi << "\n";
        }
        std::cout << "========================================================================\n";
        // ========================================================

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
            //std::cout << std::fixed << std::setprecision(2);
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
    //std::cout << std::fixed << std::setprecision(10);
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
            //std::cout << std::fixed << std::setprecision(2);
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

      if (_debugLevel) {
        std::cout << "--------------------------------------------------------------------------------\n";
        std::cout << " ev5_select_best_segments_step_06A (Resolve Shared Hits)\n";
        std::cout << "--------------------------------------------------------------------------------\n";
    }

    // =====================================================================
    // DEFINITIONS
    // =====================================================================
    const double CHI2_NDF_CUT = 5.0;  // User suggested ~5.0
    const double RATIO_CUT    = 3.0;  // The best match must be 3x closer to win exclusive rights

    // ---------------------------------------------------------------------
    // Step 1: Build a mapping from "hit index" -> "which segments contain it"
    // ---------------------------------------------------------------------
    std::unordered_map<int, std::set<int>> hitToSegments;

    if (_debugLevel) std::cout << "\n[SharedHit] --- ALL SEGMENTS AND HITS ---\n";

    for (size_t i = 0; i < segments.size(); ++i) {
        if (_debugLevel) {
            std::cout << std::format("\nSegment {} contains {} hits:\n", i, segments[i].size());
            // Print the header row just once per segment
            std::cout << "  Index | X        | Y        | Z         | Phi\n";
            std::cout << "  --------------------------------------------------\n";
        }

        for (auto& h : segments[i]) {
            hitToSegments[h.hitIndice].insert(i);

            // Print only the values, aligned to match the header
            if (_debugLevel) {
                std::cout << std::format("  {:<5} | {:<8.2f} | {:<8.2f} | {:<9.2f} | {:.3f}\n",
                                         h.hitIndice, h.x, h.y, h.z, h.phi);
            }
        }
    }

    // --- Print Summary of Shared Hits ---
    if (_debugLevel) {
        std::cout << "\n[SharedHit] --- SHARED HITS SUMMARY ---\n";

        // Table header for shared hits
        std::cout << "  Index | Segments      | X        | Y        | Z         | Phi\n";
        std::cout << "  -----------------------------------------------------------------------\n";

        bool foundShared = false;

        for (const auto& kv : hitToSegments) {
            if (kv.second.size() > 1) {
                foundShared = true;
                int hitIdx = kv.first;

                // Group the segment indices into a single string to fit nicely in the column
                std::string segList;
                for (int segIdx : kv.second) {
                    segList += std::to_string(segIdx) + " ";
                }

                // Grab coordinates from the first segment that owns this hit
                int firstSegIdx = *kv.second.begin();
                auto it = std::find_if(segments[firstSegIdx].begin(), segments[firstSegIdx].end(),
                                       [&](const auto& h){ return h.hitIndice == hitIdx; });

                if (it != segments[firstSegIdx].end()) {
                    // Print the tabular row
                    std::cout << std::format("  {:<5} | {:<13} | {:<8.2f} | {:<8.2f} | {:<9.2f} | {:.3f}\n",
                                             hitIdx, segList, it->x, it->y, it->z, it->phi);
                }
            }
        }
        if (!foundShared) {
            std::cout << "  No shared hits found in this event.\n";
        }
        std::cout << "--------------------------------------------------------------------------------\n";
    }

// ---------------------------------------------------------------------
// Step 2: Resolve Shared Hits (Corrected)
// ---------------------------------------------------------------------
for (auto& kv : hitToSegments) {
    // kv.first is the HitIndex
    // kv.second is the Set of Segment Indices containing this hit
    if (kv.second.size() < 2) continue; // If hit is unique (not shared), skip.

    int hitIdx = kv.first;
    std::set<int>& segSet = kv.second;

    struct SegCandidate {
        int segIdx;
        double chi2ndf;
        double distance; // Geometric distance (residual)
    };
    std::vector<SegCandidate> goodCandidates;

    // Print the header for this specific hit's evaluation
    if (_debugLevel) {
      std::cout << std::format("\n  Evaluating Hit {} shared by {} segments:\n", hitIdx, segSet.size());
      std::cout << "  Seg | Chi2/NDF | Alpha    | Beta     | Hit Phi | Pred Phi | Resid (rad)\n";
      std::cout << "  ---------------------------------------------------------------------------\n";
    }

    // --- SUB-STEP A: Check Chi2/NDF Quality ---
    for (int segIdx : segSet) {
        auto& seg = segments[segIdx];

        // 1. Run the fit on this segment
        findchisq_ver2(seg);

        // 2. Retrieve values directly from the fitter class
        double this_chi2 = _lineFitter.chi2Dof();
        double alpha     = _lineFitter.dydx(); // The Slope (dphi/dz)
        double beta      = _lineFitter.y0();   // The Intercept (phi at z=0)

        // 3. Threshold Check
        // If the segment is "bad" (chi2 > 5), we ignore it as a candidate for this hit.
        if (this_chi2 < CHI2_NDF_CUT) {

            // 4. Calculate the distance (residual) of THIS hit from the line
            // We need to find the specific hit object to get its Z and Phi coordinates
            auto it = std::find_if(seg.begin(), seg.end(),
                                   [&](auto& h){ return h.hitIndice == hitIdx; });

            if (it != seg.end()) {
                double hitZ   = it->z;
                double hitPhi = it->phi;

                // Mathematical Distance: | Measured - Predicted |
                // Predicted Phi = alpha * Z + beta
                double prediction = alpha * hitZ + beta;
                double residual   = std::abs(hitPhi - prediction);

                // Store this segment as a valid candidate
                goodCandidates.push_back({segIdx, this_chi2, residual});
            }
        }
    }

    // --- SUB-STEP B: Decision Making ---

    // Case 1: All segments failed the Chi2 check (All are > 5.0)
    if (goodCandidates.empty()) {
        std::cout << "Hit " << hitIdx << ": All sharing segments have poor Chi2. Removing from ALL.\n";
        // Remove this hit from *every* segment in the original set
        for (int segIdx : segSet) {
            auto& seg = segments[segIdx];
            auto& diag = diag_segments[segIdx];
            // Remove hit from hits vector
            seg.erase(std::remove_if(seg.begin(), seg.end(),
                [&](auto& h){ return h.hitIndice == hitIdx; }), seg.end());
            // Remove hit metadata from diag vector (if aligned)
            // Note: Adjust criteria if diag stores differently
             diag.erase(std::remove_if(diag.begin(), diag.end(),
                [&](auto& d){ return d.reference_point == hitIdx; }), diag.end());
        }
        continue; // Done with this hit
    }

    // Case 2: We have at least one good segment. Let's compare distances.
    // Sort candidates by distance (closest/smallest residual first)
    std::sort(goodCandidates.begin(), goodCandidates.end(),
        [](const SegCandidate& a, const SegCandidate& b) {
            return a.distance < b.distance;
        });

    int bestSegIdx = goodCandidates[0].segIdx;
    double bestDist = goodCandidates[0].distance;

    // Determine if we should be Exclusive or Inclusive (Shared)
    bool makeExclusive = false;

    if (goodCandidates.size() > 1) {
        double secondBestDist = goodCandidates[1].distance;

        // Ratio Check: Is the best one significantly better?
        // e.g. if Best=1.0 and Second=1.5 (Ratio 1.5), Keep Shared.
        //      if Best=1.0 and Second=5.0 (Ratio 5.0), Give to Best.
        if (secondBestDist > (RATIO_CUT * bestDist)) {
            makeExclusive = true;
            std::cout << "Hit " << hitIdx << ": Winner found (Ratio " << secondBestDist/bestDist << " > " << RATIO_CUT << "). Assigned to Seg " << bestSegIdx << "\n";
        } else {
            std::cout << "Hit " << hitIdx << ": Ambiguous (Ratio " << secondBestDist/bestDist << "). Keeping shared among good segments.\n";
        }
    } else {
        // Only one segment passed the Chi2 cut, so it wins by default
        makeExclusive = true;
        std::cout << "Hit " << hitIdx << ": Only one valid segment (Seg " << bestSegIdx << "). Assigned exclusively.\n";
    }

    // --- SUB-STEP C: Execute Removal ---

    // We iterate over the ORIGINAL set of segments that claimed this hit
    for (int segIdx : segSet) {
        bool keepHit = false;

        if (makeExclusive) {
            // If exclusive, keep ONLY in the best segment
            if (segIdx == bestSegIdx) keepHit = true;
        } else {
            // If shared/ambiguous, keep in ANY segment that passed the Chi2 cut
            for (const auto& cand : goodCandidates) {
                if (cand.segIdx == segIdx) {
                    keepHit = true;
                    break;
                }
            }
        }

        if (!keepHit) {
            // Remove the hit from this segment
            auto& seg = segments[segIdx];
            auto& diag = diag_segments[segIdx];
            seg.erase(std::remove_if(seg.begin(), seg.end(),
                [&](auto& h){ return h.hitIndice == hitIdx; }), seg.end());
            // Adjust diag removal logic as per your struct details
             diag.erase(std::remove_if(diag.begin(), diag.end(),
                [&](auto& d){ return d.reference_point == hitIdx; }), diag.end());
        }
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
  }//end ev5_select_best_segments_step_06A

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
//--------------------------------------------------------------------------------//
// ev5_fit_slope_ver4
//
// Fits a single straight line phi = alpha*z + beta through the hits stored in
// _segmentHits.at(index), resolving the 2*pi ambiguity of helixPhi as it goes.
//
// 2*pi disambiguation policy:
//   Hits close to the reference hit (in station number) are NOT used to
//   predict phi from the running line fit, because a line built from very
//   few hits spanning a short station range is numerically unreliable and
//   can send the 2*pi assignment into a runaway error. Those hits are added
//   to the fitter using only their own stored nturn (no line-based
//   correction). Line-prediction-based 2*pi resolution is only allowed once
//   a hit satisfies BOTH of the following:
//     (a) the hit's station is at least kMinStationGapForTrustedFit stations
//         downstream of the reference hit's station
//     (b) the fitter already holds at least kMinHitsForTrustedPrediction points
//   Condition (a) is the primary criterion for this detector geometry;
//   condition (b) is kept as an additional safety net for segments that
//   happen to have very few hits per station.
//
// Debug levels (controlled by _debugLevel):
//   > 0 : one-line summary of the input, the reference hit, and the final
//         fit result
//   > 1 : per-hit trace - trust decision (station gap / hit count), branch
//         taken, and every phi candidate considered
//   > 2 : full dump of every (z,phi,weight) point handed to _lineFitter, plus
//         a duplicate-hit check (an ALERT line if the same hitIndice is added
//         to the fitter more than once)
//--------------------------------------------------------------------------------//
void PhiZSeedFinder::ev5_fit_slope_ver4(int index, double& alpha, double& alphaError, double& beta, double& betaError, double& chindf) {

  const bool dbg  = (_debugLevel > 0);
  const bool dbg2 = (_debugLevel > 0);
  const bool dbg3 = (_debugLevel > 0);

  const std::string T = "[fit_slope_ver4 idx=" + std::to_string(index) + "]";

  auto& hits = _segmentHits.at(index);

  // --- tunable thresholds for trusting the line-based 2pi resolution ---
  const int kMinStationGapForTrustedFit  = 2;   // hit.station - Reference_station
  const int kMinHitsForTrustedPrediction = 4;   // safety net: minimum points already in the fitter

  if (dbg) {
    std::cout << "\n" << T << " ---------------------------------------------" << std::endl;
    std::cout << T << " nHit = " << hits.size() << std::endl;
  }
  if (dbg2) {
    std::cout << T << " input hits (z-sorted order expected):" << std::endl;
    std::cout << T << "   j  hitIndice  station  segIdx  nturn      z      helixPhi" << std::endl;
    for (size_t j = 0; j < hits.size(); ++j) {
      const auto& h = hits[j];
      std::cout << std::fixed
                << T << "  " << std::setw(3) << j
                << std::setw(10) << h.hitIndice
                << std::setw(9)  << h.station
                << std::setw(8)  << h.segmentIndex
                << std::setw(7)  << h.nturn
                << std::setw(10) << std::setprecision(1) << h.z
                << std::setw(12) << std::setprecision(4) << h.helixPhi
                << std::endl;
    }
  }

  //---------------------------------------------------------------------------
  // Reference hit: the first entry in the array is the anchor for the whole
  // fit. This assumes the array is already sorted by z (most upstream hit
  // first) - the caller is responsible for that ordering.
  //---------------------------------------------------------------------------
  const int    Reference_HitIndex     = 0;
  const int    Reference_segmentIndex = hits.at(Reference_HitIndex).segmentIndex;
  const int    Reference_station      = hits.at(Reference_HitIndex).station;
  const double Reference_Phi          = hits.at(Reference_HitIndex).helixPhi;
  const double Reference_z            = hits.at(Reference_HitIndex).z;

  hits.at(Reference_HitIndex).phiDiag = Reference_Phi;

  if (dbg) {
    std::cout << std::fixed << std::setprecision(4)
              << T << " reference hit : j=0  hitIndice=" << hits.at(0).hitIndice
              << "  segmentIndex=" << Reference_segmentIndex
              << "  station=" << Reference_station
              << "  z=" << std::setprecision(1) << Reference_z
              << "  helixPhi=" << std::setprecision(4) << Reference_Phi << std::endl;
    std::cout << T << " trust thresholds : station gap >= " << kMinStationGapForTrustedFit
              << "  AND  fitter points >= " << kMinHitsForTrustedPrediction << std::endl;
  }

  ::LsqSums2 _lineFitter;
  _lineFitter.clear();
  _lineFitter.addPoint(Reference_z, Reference_Phi, 0.1);

  // bookkeeping used only for the debug dump / duplicate check below;
  // does not affect the fit itself.
  std::vector<std::tuple<int,double,double,double>> addedPoints;  // (hitIndice, z, phi, weight)
  std::map<int,int> addCount;                                    // hitIndice -> number of times added
  addedPoints.emplace_back(hits.at(Reference_HitIndex).hitIndice, Reference_z, Reference_Phi, 0.1);
  addCount[hits.at(Reference_HitIndex).hitIndice]++;

  //---------------------------------------------------------------------------
  // Unique segmentIndex values present in this (possibly merged) array.
  //---------------------------------------------------------------------------
  std::set<int> uniqueSegmentIndices;
  for (const auto& h : hits) uniqueSegmentIndices.insert(h.segmentIndex);

  if (dbg) {
    std::cout << T << " unique segmentIndex values (" << uniqueSegmentIndices.size() << "): {";
    bool first = true;
    for (int s : uniqueSegmentIndices) { std::cout << (first ? "" : ",") << s; first = false; }
    std::cout << "}   Reference_segmentIndex=" << Reference_segmentIndex << std::endl;
  }

  for (int segIdx : uniqueSegmentIndices) {

    if (dbg2) {
      std::cout << T << " =============================================" << std::endl;
      std::cout << T << " outer loop segIdx = " << segIdx
                << (segIdx == Reference_segmentIndex ? "  (== Reference_segmentIndex)"
                                                     : "  (!= Reference_segmentIndex)")
                << std::endl;
      std::cout << T << " =============================================" << std::endl;
    }

    for (size_t j = 0; j < hits.size(); ++j) {

      const auto& h          = hits.at(j);
      const double z         = h.z;
      const double phiError2 = h.helixPhiError2;
      const double phiWeight = 1.0 / phiError2;

      std::string branch;
      double      phiChosen = 0.0;
      bool        added     = false;

      if (segIdx == Reference_segmentIndex) {

        if (Reference_HitIndex == (int)j) {
          if (dbg2) std::cout << T << "  j=" << j << "  hitIndice=" << h.hitIndice
                              << "  -> SKIP (reference hit itself)" << std::endl;
          continue;
        }

        // Small correction for the 2*pi boundary relative to the reference
        // phi. This is independent of the trust decision below: it only
        // keeps deltaPhi within [-pi,pi] before any line-based reasoning.
        float deltaPhi = h.helixPhi - Reference_Phi;
        int   turns    = 0;
        if (deltaPhi >  M_PI) turns--;
        if (deltaPhi < -M_PI) turns++;
        double phi = h.helixPhi + turns * 2 * M_PI;

        const int  stationGap      = h.station - Reference_station;
        const bool stationOK       = (stationGap >= kMinStationGapForTrustedFit);
        const bool hitCountOK      = (_lineFitter.qn() >= kMinHitsForTrustedPrediction);
        const bool trustPrediction = stationOK && hitCountOK;

        if (dbg2)
          std::cout << std::fixed << std::setprecision(4)
                    << T << "  j=" << j << "  hitIndice=" << h.hitIndice
                    << "  station=" << h.station << " (gap=" << stationGap << ")"
                    << "  qn=" << _lineFitter.qn()
                    << "  trust=" << (trustPrediction ? "YES" : "NO")
                    << " (stationOK=" << stationOK << ", hitCountOK=" << hitCountOK << ")"
                    << std::endl;

        if (!trustPrediction) {
          // Station gap and/or hit count threshold not met yet: do not use
          // the line to resolve the 2pi ambiguity. Use the hit's own nturn
          // directly. The hit is still added to the fitter so that qn() and
          // the station gap can grow toward the trusted regime.
          branch    = "UNTRUSTED (station/hit-count threshold not met)";
          phi       = phi + h.nturn * 2 * M_PI;
          phiChosen = phi;
          _lineFitter.addPoint(z, phi, phiWeight);
          added = true;

        } else {
          // Both thresholds satisfied: the running line fit is now
          // numerically trustworthy enough to be used for 2pi disambiguation.
          double lineSlope     = _lineFitter.dydx();
          double lineIntercept = _lineFitter.y0();
          double predictedPhi  = lineSlope * z + lineIntercept;

          if (turns != 0) {
            branch = "WRAP (turns=" + std::to_string(turns) + ")";
            phi = phi + h.nturn * 2 * M_PI;
            double phiOrg = h.helixPhi + h.nturn * 2 * M_PI;
            double diffPhi[2] = { predictedPhi - phi, predictedPhi - phiOrg };
            if (dbg2)
              std::cout << std::fixed << std::setprecision(4)
                        << T << "    predictedPhi=" << predictedPhi
                        << "  candidate(turns)=" << phi << " diff=" << diffPhi[0]
                        << "  candidate(orig)=" << phiOrg << " diff=" << diffPhi[1]
                        << std::endl;
            if (std::abs(diffPhi[1]) < std::abs(diffPhi[0])) {
              phi = phiOrg;
              branch += " -> chose ORIGINAL";
            } else {
              branch += " -> chose TURNS-CORRECTED";
            }
            phiChosen = phi;
            _lineFitter.addPoint(z, phi, phiWeight);
            added = true;

          } else {
            branch = "NO-WRAP (candidate search)";
            double phiOrg = phi;
            double candidates[3] = { phiOrg, phiOrg + 2 * M_PI, phiOrg - 2 * M_PI };
            int    best     = 0;
            double bestDiff = std::fabs(predictedPhi - candidates[0]);
            for (int k = 1; k < 3; ++k) {
              double d = std::fabs(predictedPhi - candidates[k]);
              if (d < bestDiff) { bestDiff = d; best = k; }
            }
            if (dbg2)
              std::cout << std::fixed << std::setprecision(4)
                        << T << "    predictedPhi=" << predictedPhi
                        << "  candidates=[" << candidates[0] << "," << candidates[1]
                        << "," << candidates[2] << "]  chosen index=" << best
                        << " (bestDiff=" << bestDiff << ")" << std::endl;
            phi = candidates[best];
            phi += h.nturn * 2 * M_PI;
            phiChosen = phi;
            _lineFitter.addPoint(z, phi, phiWeight);
            added = true;
            branch += "  best=" + std::to_string(best);
          }
        }

      } else {
        // Hit belongs to a different segmentIndex than the reference: use
        // its own nturn directly, no line-based disambiguation.
        branch     = "OTHER-SEGMENT (direct nturn)";
        double phi = h.helixPhi + h.nturn * 2 * M_PI;
        phiChosen  = phi;
        _lineFitter.addPoint(z, phi, phiWeight);
        added = true;
      }

      if (added) {
        addCount[h.hitIndice]++;
        addedPoints.emplace_back(h.hitIndice, z, phiChosen, phiWeight);

        if (dbg2) {
          std::cout << std::fixed << std::setprecision(4)
                    << T << "    ADD  hitIndice=" << h.hitIndice
                    << "  z=" << std::setprecision(1) << z
                    << "  helixPhi=" << std::setprecision(4) << h.helixPhi
                    << "  phi(used)=" << phiChosen
                    << "  branch=[" << branch << "]"
                    << "  qn_after=" << _lineFitter.qn() << std::endl;
        }
        if (addCount[h.hitIndice] > 1) {
          std::cout << T << "  !!! ALERT !!! hitIndice=" << h.hitIndice
                    << " has been added to _lineFitter " << addCount[h.hitIndice]
                    << " times (outer segIdx=" << segIdx
                    << ", own segmentIndex=" << h.segmentIndex
                    << ") - this hit is being double-counted in the fit." << std::endl;
        }
      }
    }
  }

  //---------------------------------------------------------------------------
  // Final fit result
  //---------------------------------------------------------------------------
  alpha      = _lineFitter.dydx();
  alphaError = _lineFitter.dydxErr();
  beta       = _lineFitter.y0();
  betaError  = _lineFitter.y0Err();
  chindf     = _lineFitter.chi2Dof();

  if (dbg3) {
    std::cout << T << " points actually given to the fitter (in add order):" << std::endl;
    std::cout << T << "   #   hitIndice        z        phi     weight" << std::endl;
    for (size_t k = 0; k < addedPoints.size(); ++k) {
      const auto& [hid, zAdded, phiAdded, w] = addedPoints[k];
      std::cout << std::fixed
                << T << "  " << std::setw(3) << k
                << std::setw(10) << hid
                << std::setw(11) << std::setprecision(1) << zAdded
                << std::setw(12) << std::setprecision(4) << phiAdded
                << std::setw(11) << std::setprecision(2) << w << std::endl;
    }
  }

  const int nPointsAdded = _lineFitter.qn();
  if (dbg) {
    std::cout << T << " nInputHits=" << hits.size()
              << "  nPointsInFitter=" << nPointsAdded
              << (uniqueSegmentIndices.size() == 1 && nPointsAdded != (int)hits.size()
                    ? "   <-- unexpected mismatch (single segment but counts differ)"
                    : "") << std::endl;
  }
  if (dbg) {
    std::cout << std::scientific << std::setprecision(4)
              << T << " RESULT alpha=" << alpha << " +- " << alphaError
              << "  beta=" << beta << " +- " << betaError
              << std::fixed << std::setprecision(3)
              << "  chi2/ndf=" << chindf
              << "  (nPointsInFitter=" << nPointsAdded << ")" << std::endl;
   }

}//end ev5_fit_slope_ver4

//----------------------------------------------------------------------------
//--------------------------------------------------------------------------------//
// ev5_fit_slope_ver5
//
// Fits a single straight line phi = alpha*z + beta through the hits stored in
// _segmentHits.at(index), resolving the 2*n*pi ambiguity of helixPhi as it goes.
// This routine is used ONLY to judge whether two PhiZ segments can be merged;
// nothing it computes is meant to survive into the final helix parameters, and
// it is deliberately free of side effects on anything outside _segmentHits.
//
// Physics model
// -------------
// A particle spirals through the tracker, so its trajectory is a straight line
// in the (z, helixPhi) plane. helixPhi comes from atan2 and is therefore folded
// into (-pi,pi], while the true phi keeps growing with z. Two distinct effects
// have to be undone, and they are NOT the same thing:
//
//   m_i  (intra-segment winding)
//        Even inside one segment the trajectory keeps turning, so as z grows
//        the folded helixPhi wraps past +/-pi one or more times. Because the
//        segment is continuously observed, m_i follows unambiguously from the
//        hit ordering - it varies from hit to hit and is NOT constant within a
//        segment.
//
//   N_S  (per-segment offset)
//        While the particle crosses the hollow central region of the tracker no
//        hits are produced, and it reappears one or more full turns later as a
//        separate segment with a different segmentIndex. That invisible number
//        of turns is unknown, but it is a single integer shared by every hit of
//        the segment, since the segment itself is continuously observed.
//
// So the total turn number of hit i in segment S is
//
//        n_i = m_i + N_S
//
// and only the N_S part is common to a segment. Treating the whole segment as
// one common n (as an earlier version did) is wrong; resolving every hit
// independently is also wrong, because a single noisy hit could then be pushed
// a full turn away from its own segment.
//
// Dependency order
// ----------------
//   alpha (trusted slope)  ->  m_i  ->  N_S
//
// m_i is measured relative to the segment's OWN first hit, so it is a purely
// relative quantity and does not need N_S. It also uses only alpha, never the
// intercept beta: early in the fit beta is the unstable parameter, and keeping
// it out of this step is what prevents the runaway seen in earlier versions.
// N_S is then fixed by comparing the m-corrected segment against the global
// line, and therefore has to come second.
//
// Algorithm
// ---------
//   PHASE 1 (seed, reference segment)
//     Reference hit = smallest z in the whole array. Hits are added with only a
//     [-pi,pi] wrap relative to the reference phi - no line prediction - until
//     the fit is trustworthy: station gap >= kMinStationGapForTrustedFit AND
//     fitter points >= kMinHitsForTrustedPrediction. Over such a short lever
//     arm the trajectory cannot have turned by more than pi, so the simple wrap
//     is safe here. Seed hits are added to the fitter so it grows toward the
//     trusted regime.
//
//   PHASE 2 (reference segment, m_i)
//     Once a trusted alpha exists, m_i is recomputed for EVERY hit of the
//     reference segment - including the seed hits - and the fitter is rebuilt
//     from scratch. This repairs any seed hit that the naive wrap got wrong.
//     N_S is 0 for the reference segment by definition: it anchors the fit.
//
//   PHASE 3 (remaining segments, m_i then N_S)
//     Segments are visited in order of increasing z. For each one, m_i is
//     derived from the current alpha, then a single N_S is chosen by scanning
//     candidate offsets and keeping the one with the smallest weighted residual
//     against the current line. All hits are then added with n_i = m_i + N_S,
//     after which alpha and beta update for the next segment.
//
// Note on the running fit: LsqSums2 accumulates sums, so dydx() and y0() change
// after every addPoint(). Predictions therefore always reflect every hit added
// so far - which is also why a hit added with a wrong turn number immediately
// contaminates the next prediction.
//
// Debug levels (_debugLevel):
//   > 0 : input summary, reference hit, phase transitions, per-segment N_S,
//         final result
//   > 1 : per-hit trace with trust decision, m_i and n_i
//   > 2 : full N_S scan per segment and the complete list of fitted points
//--------------------------------------------------------------------------------//
void PhiZSeedFinder::ev5_fit_slope_ver5(int index, double& alpha, double& alphaError,
                                        double& beta, double& betaError, double& chindf) {

  const bool dbg  = (_debugLevel > 0);
  const bool dbg2 = (_debugLevel > 1);
  const bool dbg3 = (_debugLevel > 2);

  const std::string T = "[fit_slope_ver5 idx=" + std::to_string(index) + "]";

  auto& hits = _segmentHits.at(index);

  // --- tunable thresholds ---
  const int kMinStationGapForTrustedFit  = 2;   // |hit.station - Reference_station|
  const int kMinHitsForTrustedPrediction = 4;   // minimum points already in the fitter
  const int kMaxTurnSearch               = 6;   // range of N_S offsets scanned

  alpha = alphaError = beta = betaError = chindf = 0.0;
  if (hits.empty()) {
    if (dbg) std::cout << T << " empty segment, nothing to fit" << std::endl;
    return;
  }

  //---------------------------------------------------------------------------
  // Order hits by z without disturbing the caller's array: work on indices.
  //---------------------------------------------------------------------------
  std::vector<size_t> order(hits.size());
  for (size_t k = 0; k < hits.size(); ++k) order[k] = k;
  std::sort(order.begin(), order.end(),
            [&hits](size_t a, size_t b) { return hits[a].z < hits[b].z; });

  const size_t Reference_HitIndex     = order.front();   // smallest z overall
  const int    Reference_segmentIndex = hits.at(Reference_HitIndex).segmentIndex;
  const int    Reference_station      = hits.at(Reference_HitIndex).station;
  const double Reference_Phi          = hits.at(Reference_HitIndex).helixPhi;
  const double Reference_z            = hits.at(Reference_HitIndex).z;

  //---------------------------------------------------------------------------
  // Group hits by segmentIndex (each group z-ordered); visit groups by min z.
  //---------------------------------------------------------------------------
  std::map<int, std::vector<size_t>> hitsBySegment;
  for (size_t p : order) hitsBySegment[hits[p].segmentIndex].push_back(p);

  std::vector<int> segmentOrder;
  for (const auto& kv : hitsBySegment) segmentOrder.push_back(kv.first);
  std::sort(segmentOrder.begin(), segmentOrder.end(),
            [&](int a, int b) {
              return hits[hitsBySegment[a].front()].z < hits[hitsBySegment[b].front()].z;
            });

  if (dbg) {
    std::cout << "\n" << T << " ---------------------------------------------" << std::endl;
    std::cout << T << " nHit = " << hits.size()
              << " , nSegment = " << segmentOrder.size() << std::endl;
    std::cout << std::fixed << std::setprecision(4)
              << T << " reference hit (smallest z) : arrayIdx=" << Reference_HitIndex
              << "  hitIndice=" << hits.at(Reference_HitIndex).hitIndice
              << "  segmentIndex=" << Reference_segmentIndex
              << "  station=" << Reference_station
              << "  z=" << std::setprecision(1) << Reference_z
              << "  helixPhi=" << std::setprecision(4) << Reference_Phi << std::endl;
    std::cout << T << " trust thresholds : station gap >= " << kMinStationGapForTrustedFit
              << "  AND  fitter points >= " << kMinHitsForTrustedPrediction << std::endl;
    for (int s : segmentOrder) {
      const auto& g = hitsBySegment[s];
      std::cout << std::fixed << std::setprecision(1)
                << T << "   segmentIndex=" << s
                << "  nHit=" << g.size()
                << "  z [" << hits[g.front()].z << "," << hits[g.back()].z << "]"
                << "  station [" << hits[g.front()].station << "," << hits[g.back()].station << "]"
                << (s == Reference_segmentIndex ? "   <== reference segment" : "")
                << std::endl;
    }
  }

  ::LsqSums2 _lineFitter;
  _lineFitter.clear();

  std::vector<std::tuple<int,double,double,int>> fitted;   // (hitIndice, z, phi, n)

  // Add one hit with turn number n; also records it for the debug dump.
  auto addHit = [&](size_t arrayIdx, int n, const char* branch) {
    auto& h = hits.at(arrayIdx);
    const double phi = h.helixPhi + n * 2.0 * M_PI;
    _lineFitter.addPoint(h.z, phi, 1.0 / h.helixPhiError2);
    h.nturn   = n;      // diagnostic only; not relied on downstream
    h.phiDiag = phi;    // diagnostic only; the phi actually fitted
    fitted.emplace_back(h.hitIndice, h.z, phi, n);
    if (dbg2)
      std::cout << std::fixed << std::setprecision(4)
                << T << "    ADD hitIndice=" << std::setw(4) << h.hitIndice
                << "  station=" << std::setw(2) << h.station
                << "  z=" << std::setw(9) << std::setprecision(1) << h.z
                << "  helixPhi=" << std::setw(8) << std::setprecision(4) << h.helixPhi
                << "  n=" << std::setw(3) << n
                << "  phi=" << std::setw(9) << phi
                << "  [" << branch << "]"
                << "  qn_after=" << _lineFitter.qn() << std::endl;
  };

  // m_i for one hit: winding relative to the segment's own first hit, using
  // only the slope. Independent of N_S and of the intercept.
  auto intraSegmentWinding = [&](size_t arrayIdx, size_t firstIdx, double slope) {
    const auto& h  = hits.at(arrayIdx);
    const auto& h0 = hits.at(firstIdx);
    const double dPhiExpected = slope * (h.z - h0.z);
    const double dPhiObserved = h.helixPhi - h0.helixPhi;
    return (int)std::round((dPhiExpected - dPhiObserved) / (2.0 * M_PI));
  };

  const std::vector<size_t>& refGroup = hitsBySegment[Reference_segmentIndex];

  //---------------------------------------------------------------------------
  // PHASE 1 : seed the fit from the reference segment
  //---------------------------------------------------------------------------
  if (dbg) std::cout << T << " --- PHASE 1 : seed (reference segment "
                     << Reference_segmentIndex << ") ---" << std::endl;

  addHit(Reference_HitIndex, 0, "PHASE1 reference hit");

  bool trustReached = false;
  for (size_t k = 0; k < refGroup.size(); ++k) {
    const size_t p = refGroup[k];
    if (p == Reference_HitIndex) continue;

    const auto& h = hits.at(p);
    const int  stationGap = std::abs(h.station - Reference_station);
    const bool stationOK  = (stationGap >= kMinStationGapForTrustedFit);
    const bool countOK    = (_lineFitter.qn() >= kMinHitsForTrustedPrediction);

    if (stationOK && countOK) {
      trustReached = true;
      if (dbg)
        std::cout << std::scientific << std::setprecision(4)
                  << T << "   trust reached at hitIndice=" << h.hitIndice
                  << " (stationGap=" << stationGap << ", qn=" << _lineFitter.qn() << ")"
                  << "  seed slope=" << _lineFitter.dydx()
                  << std::fixed << std::endl;
      break;
    }

    // Short lever arm: the trajectory cannot have turned by more than pi here,
    // so a plain [-pi,pi] wrap against the reference phi is sufficient.
    const double deltaPhi = h.helixPhi - Reference_Phi;
    int n = 0;
    if (deltaPhi >  M_PI) n = -1;
    if (deltaPhi < -M_PI) n = +1;
    addHit(p, n, "PHASE1 seed wrap");
  }

  //---------------------------------------------------------------------------
  // PHASE 2 : reference segment, proper m_i for every hit (N_S = 0)
  //---------------------------------------------------------------------------
  if (trustReached) {

    const double alphaSeed = _lineFitter.dydx();

    if (dbg)
      std::cout << std::scientific << std::setprecision(4)
                << T << " --- PHASE 2 : reference segment m_i with alphaSeed="
                << alphaSeed << std::fixed << " ---" << std::endl;

    // Rebuild the fit from scratch so that seed hits wrapped by the naive rule
    // are corrected too.
    _lineFitter.clear();
    fitted.clear();

    for (size_t p : refGroup) {
      const int m = intraSegmentWinding(p, Reference_HitIndex, alphaSeed);
      addHit(p, m, "PHASE2 m_i (N_S=0)");
    }

    if (dbg)
      std::cout << std::scientific << std::setprecision(4)
                << T << "   after PHASE 2 : alpha=" << _lineFitter.dydx()
                << " beta=" << _lineFitter.y0()
                << std::fixed << std::setprecision(3)
                << " chi2/ndf=" << _lineFitter.chi2Dof()
                << " (qn=" << _lineFitter.qn() << ")" << std::endl;

  } else if (dbg) {
    std::cout << T << " --- PHASE 2 skipped : trust never reached, keeping seed wraps"
              << " (qn=" << _lineFitter.qn() << ") ---" << std::endl;
  }

  //---------------------------------------------------------------------------
  // PHASE 3 : remaining segments, m_i then a single N_S each
  //---------------------------------------------------------------------------
  for (int segIdx : segmentOrder) {
    if (segIdx == Reference_segmentIndex) continue;

    const std::vector<size_t>& group = hitsBySegment[segIdx];

    if (dbg) std::cout << T << " --- PHASE 3 : segment " << segIdx
                       << " (" << group.size() << " hits) ---" << std::endl;

    if (_lineFitter.qn() < 2) {
      if (dbg) std::cout << T << "   no usable line (qn=" << _lineFitter.qn()
                         << "), keeping stored nturn" << std::endl;
      for (size_t p : group) addHit(p, hits.at(p).nturn, "PHASE3 fallback");
      continue;
    }

    const double lineSlope     = _lineFitter.dydx();
    const double lineIntercept = _lineFitter.y0();
    const size_t firstIdx      = group.front();     // smallest z in this segment

    // --- step A: m_i inside this segment, from the slope only ---
    std::map<size_t,int> m;
    for (size_t p : group) m[p] = intraSegmentWinding(p, firstIdx, lineSlope);

    if (dbg2) {
      std::cout << std::fixed << std::setprecision(4)
                << T << "   line: slope=" << lineSlope << " intercept=" << lineIntercept
                << std::endl;
      std::cout << T << "   intra-segment winding m_i:" << std::endl;
      std::cout << T << "     hitIndice  station        z   helixPhi   m_i" << std::endl;
      for (size_t p : group) {
        const auto& h = hits.at(p);
        std::cout << std::fixed
                  << T << "     " << std::setw(9) << h.hitIndice
                  << std::setw(9)  << h.station
                  << std::setw(10) << std::setprecision(1) << h.z
                  << std::setw(11) << std::setprecision(4) << h.helixPhi
                  << std::setw(6)  << m[p] << std::endl;
      }
    }

    // --- step B: one common offset N_S for the whole segment ---
    int    bestN    = 0;
    double bestCost = std::numeric_limits<double>::max();

    for (int N = -kMaxTurnSearch; N <= kMaxTurnSearch; ++N) {
      double cost = 0.0;
      for (size_t p : group) {
        const auto& h = hits.at(p);
        const double phi   = h.helixPhi + (m[p] + N) * 2.0 * M_PI;
        const double resid = phi - (lineSlope * h.z + lineIntercept);
        cost += (1.0 / h.helixPhiError2) * resid * resid;
      }
      if (dbg3)
        std::cout << T << "     trial N_S=" << std::setw(3) << N
                  << "  weighted cost=" << std::scientific << std::setprecision(4)
                  << cost << std::fixed << std::endl;
      if (cost < bestCost) { bestCost = cost; bestN = N; }
    }

    if (dbg)
      std::cout << T << "   chosen N_S=" << bestN
                << "  (weighted cost=" << std::scientific << std::setprecision(4)
                << bestCost << std::fixed << ")" << std::endl;

    for (size_t p : group) addHit(p, m[p] + bestN, "PHASE3 m_i + N_S");

    if (dbg)
      std::cout << std::scientific << std::setprecision(4)
                << T << "   after segment " << segIdx << " : alpha=" << _lineFitter.dydx()
                << " beta=" << _lineFitter.y0()
                << std::fixed << std::setprecision(3)
                << " chi2/ndf=" << _lineFitter.chi2Dof()
                << " (qn=" << _lineFitter.qn() << ")" << std::endl;
  }

  //---------------------------------------------------------------------------
  // Final fit result
  //---------------------------------------------------------------------------
  alpha      = _lineFitter.dydx();
  alphaError = _lineFitter.dydxErr();
  beta       = _lineFitter.y0();
  betaError  = _lineFitter.y0Err();
  chindf     = _lineFitter.chi2Dof();

  if (dbg3) {
    std::cout << T << " points given to the fitter (in add order):" << std::endl;
    std::cout << T << "   #   hitIndice        z        phi    n     resid" << std::endl;
    for (size_t k = 0; k < fitted.size(); ++k) {
      const auto& [hid, zf, phif, nf] = fitted[k];
      const double resid = phif - (alpha * zf + beta);
      std::cout << std::fixed
                << T << "  " << std::setw(3) << k
                << std::setw(10) << hid
                << std::setw(11) << std::setprecision(1) << zf
                << std::setw(12) << std::setprecision(4) << phif
                << std::setw(5)  << nf
                << std::setw(10) << std::setprecision(4) << resid << std::endl;
    }
  }

  if (dbg) {
    std::cout << T << " nInputHits=" << hits.size()
              << "  nPointsInFitter=" << _lineFitter.qn()
              << (_lineFitter.qn() != (int)hits.size()
                    ? "   <-- mismatch: some hits were not fitted"
                    : "") << std::endl;
    std::cout << std::scientific << std::setprecision(4)
              << T << " RESULT alpha=" << alpha << " +- " << alphaError
              << "  beta=" << beta << " +- " << betaError
              << std::fixed << std::setprecision(3)
              << "  chi2/ndf=" << chindf << std::endl;
  }

}//end ev5_fit_slope_ver5



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

    // --- Helper Lambda for Printing ---
    auto printCurrentSegments = [&](const std::string& label, const std::vector<std::vector<ev5_HitsInNthStation>>& currentSegs) {
        if (!_debugLevel) return; // Only print if debug is on
        std::cout << "\n================================================================================\n";
        std::cout << " [Step 06 State] " << label << " (Count: " << currentSegs.size() << ")\n";
        std::cout << "================================================================================\n";

        for (size_t i = 0; i < currentSegs.size(); ++i) {
             std::cout << std::format("\nSegment {} contains {} hits:\n", i, currentSegs[i].size());
             std::cout << "  HitID | X        | Y        | Z         | Phi\n";
             std::cout << "  --------------------------------------------------\n";
             for (const auto& h : currentSegs[i]) {
                  // Assuming h.hitID is the unique identifier you are filtering by
                  std::cout << std::format("  {:<5} | {:<8.2f} | {:<8.2f} | {:<9.2f} | {:.3f}\n",
                                           h.hitID, h.x, h.y, h.z, h.phi);
             }
        }
        std::cout << "--------------------------------------------------------------------------------\n";
    };

    // 1. PRINT INITIAL STATE
    printCurrentSegments("INITIAL INPUT", all_ThisIsBestSegment);

    //-----------------------------------------------------------------
    //    Remove identical segments
    //-----------------------------------------------------------------

    // (1). Create Sorted Hit Lists for Uniqueness Check
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

    // Remove duplicates from the list of IDs
    std::sort(hitID_list.begin(), hitID_list.end());
    hitID_list.erase(std::unique(hitID_list.begin(), hitID_list.end()), hitID_list.end());

    //---------------------------------------------------------------------------------------------------------
    // (2). Remove Subsets (e.g., if 0-1-2 and 0-1-2-4 exist, remove 0-1-2)
    //---------------------------------------------------------------------------------------------------------
    std::vector<int> delete_index;
    delete_index.clear();
    for (int i = 0; i < (int)hitID_list.size(); i++) {
        for (int j = 0; j <(int)hitID_list.size(); j++) {
            if(i == j) continue;
            // If lists are identical, skip (handled by unique above, but safe to keep)
            if(hitID_list.at(j) == hitID_list.at(i)) continue;

            int size = 0;
            // Check if List I is fully contained in List J
            for (int k = 0; k <(int)hitID_list.at(j).size(); k++) {
                for (int l = 0; l <(int)hitID_list.at(i).size(); l++) {
                    if(hitID_list.at(i).at(l) == hitID_list.at(j).at(k)) size++;
                }
            }
            // If all hits in I are found in J, mark I for deletion
            if(size == (int)hitID_list.at(i).size()) delete_index.push_back(i);
        }
    }

    // Make a new hitID list after removing the subset events
    std::vector<std::vector<int>> new_hitID_list;
    for(int i = 0; i <(int)hitID_list.size(); i++){
        int go = 1;
        for(int j=0; j<(int)delete_index.size(); j++){
            if(delete_index.at(j) == i) go = 0;
        }
        if(go == 1) new_hitID_list.push_back(hitID_list.at(i));
    }

    //---------------------------------------------------------------------------------------------------------
    // (3) Map the "New HitIDs" back to Segment Objects
    //---------------------------------------------------------------------------------------------------------
    std::vector<int> find_index;
    find_index.clear();
    for (int i = 0; i < (int)new_hitID_list.size(); i++) {
        bool flag = 0;
        int index_for_segment = 0;

        // Find the first segment in the original list that matches this HitID pattern
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

    // Create the temporary "Cleaned" lists
    std::vector<std::vector<ev5_HitsInNthStation>> temp_ThisIsBestSegment;
    std::vector<std::vector<ev5_Segment>> temp_ThisIsBestSegment_Diag;
    temp_ThisIsBestSegment.clear();
    temp_ThisIsBestSegment_Diag.clear();

    for(int i=0; i<(int)all_ThisIsBestSegment.size(); i++){
        bool go = 0;
        for(int j=0; j<(int)find_index.size(); j++){
            if(find_index[j] == i) go = 1;
        }
        if(go != 1) continue;
        temp_ThisIsBestSegment.push_back(all_ThisIsBestSegment.at(i));
        temp_ThisIsBestSegment_Diag.push_back(all_ThisIsBestSegment_Diag.at(i));
    }

    // 2. PRINT AFTER SUBSET/DUPLICATE REMOVAL
    printCurrentSegments("AFTER SUBSET/DUPLICATE REMOVAL", temp_ThisIsBestSegment);

    //---------------------------------------------------------------------------------------------------------
    // Select good segments and remove unwanted segments based on Overlap/Alpha
    //---------------------------------------------------------------------------------------------------------
    /*std::vector<int> delete_segmentID;
    delete_segmentID.clear();

    for(int i=0; i<(int)temp_ThisIsBestSegment.size(); i++){
        for(int j=0; j<(int)temp_ThisIsBestSegment.size(); j++){
            if(i==j) continue;

            double alpha[2] = {temp_ThisIsBestSegment_Diag.at(i).at(0).alpha, temp_ThisIsBestSegment_Diag.at(j).at(0).alpha};
            double alpha_diff = 0.0;
            double alpha_sigma = 0.0003199; // Obtained from "ev5_plot_AlphaDiagBestSegment"

            if(alpha[0]*alpha[1] > 0){ // Both slopes same sign
                alpha_diff = fabs(alpha[0] - alpha[1]);
            }
            if(alpha[0]*alpha[1] < 0){ // Opposite signs
                alpha[0] = fabs(alpha[0]);
                alpha[1] = fabs(alpha[1]);
                alpha_diff = alpha[0] + alpha[1];
            }

            double chiNDF[2] = {temp_ThisIsBestSegment_Diag.at(i).at(0).chiNDF, temp_ThisIsBestSegment_Diag.at(j).at(0).chiNDF};
            int TotalHitinSegments[2] = {(int)temp_ThisIsBestSegment_Diag.at(i).size(), (int)temp_ThisIsBestSegment_Diag.at(j).size()};

            int overlappedHits = 0;
            for(int k=0; k<(int)temp_ThisIsBestSegment.at(i).size(); k++){
                for(int l=0; l<(int)temp_ThisIsBestSegment.at(j).size(); l++){
                    // Using hitID to check overlap
                    if(temp_ThisIsBestSegment.at(i).at(k).hitID == temp_ThisIsBestSegment.at(j).at(l).hitID) overlappedHits++;
                }
            }

            double fraction[2] = {0.0};
            if(TotalHitinSegments[0] > 0) fraction[0] = (double)overlappedHits/TotalHitinSegments[0];
            if(TotalHitinSegments[1] > 0) fraction[1] = (double)overlappedHits/TotalHitinSegments[1];

            if(alpha_diff > alpha_sigma*5) continue; // Slopes are too different, likely distinct tracks

            // (1) 2 segments are not identical but ComboHits are overlapped -> Remove the one with less unique content?
            if(fraction[0] <= 0.3 and fraction[1] >= 0.6){
                delete_segmentID.push_back(i);
            }
            // (2) 2 segments are identical/highly overlapped -> Select the one with better Chi2
            if(fraction[0] >= 0.7 and fraction[1] >= 0.7){
                if(chiNDF[0] < chiNDF[1]) delete_segmentID.push_back(j); // Delete J if I is better
                else delete_segmentID.push_back(i);                      // Delete I if J is better
            }
        }
    }*/
    //---------------------------------------------------------------------------------------------------------
    // Select good segments and remove unwanted segments based on Overlap/Alpha
    //---------------------------------------------------------------------------------------------------------
    std::vector<int> delete_segmentID;
    delete_segmentID.clear();

    for(int i=0; i<(int)temp_ThisIsBestSegment.size(); i++){
        for(int j=0; j<(int)temp_ThisIsBestSegment.size(); j++){
            if(i == j) continue;

            // --- 1. Slope Difference Check ---
            double alpha[2] = {temp_ThisIsBestSegment_Diag.at(i).at(0).alpha, temp_ThisIsBestSegment_Diag.at(j).at(0).alpha};
            double alpha_diff = 0.0;
            double alpha_sigma = 0.0003199; // Obtained from "ev5_plot_AlphaDiagBestSegment"

            if(alpha[0]*alpha[1] > 0){ // Both slopes same sign
                alpha_diff = fabs(alpha[0] - alpha[1]);
            }
            else { // Opposite signs
                alpha_diff = fabs(fabs(alpha[0]) + fabs(alpha[1])); // Corrected logic for sign flip
            }

            // If slopes are too different, they are likely different tracks crossing, so don't merge/delete.
            if(alpha_diff > alpha_sigma * 5) continue;

            // --- 2. Calculate Overlap ---
            double chiNDF[2] = {temp_ThisIsBestSegment_Diag.at(i).at(0).chiNDF, temp_ThisIsBestSegment_Diag.at(j).at(0).chiNDF};

            // Note: Ensure we divide by the correct size.
            // Using the size of the hit vector is safer for the fraction calculation.
            int size_i = (int)temp_ThisIsBestSegment.at(i).size();
            int size_j = (int)temp_ThisIsBestSegment.at(j).size();

            int overlappedHits = 0;
            for(int k=0; k < size_i; k++){
                for(int l=0; l < size_j; l++){
                    if(temp_ThisIsBestSegment.at(i).at(k).hitID == temp_ThisIsBestSegment.at(j).at(l).hitID) {
                        overlappedHits++;
                    }
                }
            }

            double frac_i = (size_i > 0) ? (double)overlappedHits / size_i : 0.0;
            double frac_j = (size_j > 0) ? (double)overlappedHits / size_j : 0.0;

            // --- 3. The "More than 50% Mutual Overlap" Logic ---
            if (frac_i > 0.5 && frac_j > 0.5) {

                // Case A: Segment i has a higher shared fraction (it is "more covered" / shorter). Remove i.
                if (frac_i > frac_j) {
                    delete_segmentID.push_back(i);
                    if(_debugLevel) std::cout << std::format("  -> Removing Seg {} (Frac {:.2f}) vs Seg {} (Frac {:.2f})\n", i, frac_i, j, frac_j);
                }
                // Case B: Segment j has a higher shared fraction. Remove j.
                else if (frac_j > frac_i) {
                    delete_segmentID.push_back(j);
                     if(_debugLevel) std::cout << std::format("  -> Removing Seg {} (Frac {:.2f}) vs Seg {} (Frac {:.2f})\n", j, frac_j, i, frac_i);
                }
                // Case C: Exact same overlap fraction (likely same length). Use Chi2 tie-breaker.
                else {
                    if (chiNDF[0] > chiNDF[1]) { // i has worse (larger) Chi2
                        delete_segmentID.push_back(i);
                    } else { // j has worse Chi2
                        delete_segmentID.push_back(j);
                    }
                }
            }

        }
    }

    std::sort(delete_segmentID.begin(), delete_segmentID.end());
    delete_segmentID.erase(std::unique(delete_segmentID.begin(), delete_segmentID.end()), delete_segmentID.end());

    // Construct the Final Output List
    all_ThisIsBestSegment.clear();
    all_ThisIsBestSegment_Diag.clear();

    for(int i=0; i<(int)temp_ThisIsBestSegment.size(); i++){
        int flag = 0;
        for(int j=0; j<(int)delete_segmentID.size(); j++){
            if(delete_segmentID[j] == i) flag = 1;
        }
        if(flag == 1) continue;

        all_ThisIsBestSegment.push_back(temp_ThisIsBestSegment.at(i));
        all_ThisIsBestSegment_Diag.push_back(temp_ThisIsBestSegment_Diag.at(i));
    }

    // 3. PRINT FINAL STATE
    printCurrentSegments("FINAL OUTPUT (AFTER OVERLAP REMOVAL)", all_ThisIsBestSegment);

    // =====================================================================
    // STEP 4: Resolve Shared Hits (The 06A Logic)
    // =====================================================================
    // This handles hits shared by distinct tracks that weren't merged/deleted above.

    const double CHI2_NDF_CUT_FINAL = 5.0;
    const double RATIO_CUT_FINAL    = 3.0;

    // 1. Build the mapping for the survivors
    std::unordered_map<int, std::set<int>> hitToSegments;
    for (size_t i = 0; i < all_ThisIsBestSegment.size(); ++i) {
        for (auto& h : all_ThisIsBestSegment[i]) {
            hitToSegments[h.hitID].insert(i);
        }
    }

    if (_debugLevel) {
        std::cout << "\n[SharedHit-Fine] Starting Final Hit-Level Resolution...\n";
    }

    for (auto& kv : hitToSegments) {
        if (kv.second.size() < 2) continue; // Only process hits shared by 2+ survivors

        int hitIdx = kv.first;
        std::set<int>& segSet = kv.second;

        struct FinalCandidate {
            int segIdx;
            double chi2ndf;
            double residual;
        };
        std::vector<FinalCandidate> goodCandidates;

        if (_debugLevel) {
            std::cout << std::format("\n  Refining Shared Hit {}:\n", hitIdx);
            std::cout << "  Seg | Chi2/NDF | Resid (rad)\n";
            std::cout << "  ----------------------------\n";
        }

        for (int segIdx : segSet) {
            auto& seg = all_ThisIsBestSegment[segIdx];

            // Perform fit to get current Alpha/Beta
            findchisq_ver2(seg);
            double this_chi2 = _lineFitter.chi2Dof();
            double alpha     = _lineFitter.dydx();
            double beta      = _lineFitter.y0();

            auto it = std::find_if(seg.begin(), seg.end(),
                                   [&](auto& h){ return h.hitID == hitIdx; });

            if (it != seg.end()) {
                double residual = std::abs(it->phi - (alpha * it->z + beta));

                if (_debugLevel) {
                    std::cout << std::format("  {:<3} | {:<8.2f} | {:.4f} {}\n",
                                             segIdx, this_chi2, residual,
                                             (this_chi2 >= CHI2_NDF_CUT_FINAL ? "(FAIL)" : ""));
                }

                if (this_chi2 < CHI2_NDF_CUT_FINAL) {
                    goodCandidates.push_back({segIdx, this_chi2, residual});
                }
            }
        }

        // --- Decision Logic ---
        if (goodCandidates.empty()) {
            if (_debugLevel) {
                std::cout << std::format("  RESULT for Hit {}: REMOVED from all segments (Zero candidates passed Chi2 cut).\n", hitIdx);
            }
            // Remove hit from all segments if none are good fits
            for (int segIdx : segSet) {
                auto& s = all_ThisIsBestSegment[segIdx];
                auto& d = all_ThisIsBestSegment_Diag[segIdx];
                s.erase(std::remove_if(s.begin(), s.end(), [&](auto& h){ return h.hitID == hitIdx; }), s.end());
                d.erase(std::remove_if(d.begin(), d.end(), [&](auto& x){ return x.reference_point == hitIdx; }), d.end());
            }
            continue;
        }

        // Sort candidates by residual (best fit first)
        // Sort candidates by residual (best fit first)
        std::sort(goodCandidates.begin(), goodCandidates.end(), [](auto& a, auto& b){ return a.residual < b.residual; });

        // --- New Robust Decision Logic (N > 2 Support) ---

        // 1. Determine which segments are "Good Enough" to keep the hit
        std::set<int> segmentsToKeep;

        // The best candidate ALWAYS keeps the hit
        int bestSegIdx = goodCandidates[0].segIdx;
        double r0 = goodCandidates[0].residual;
        segmentsToKeep.insert(bestSegIdx);

        if (_debugLevel) {
             std::cout << std::format("  RESULT for Hit {}: Winner is Seg {} (Resid {:.4f})\n", hitIdx, bestSegIdx, r0);
        }

        // Check the runners-up (indices 1 to N)
        for (size_t i = 1; i < goodCandidates.size(); ++i) {
            int currSegIdx = goodCandidates[i].segIdx;
            double r_curr  = goodCandidates[i].residual;

            // Calculate ratio relative to the winner
            // Protect against divide-by-zero if perfect fit
            double ratio = (r0 > 1e-9) ? (r_curr / r0) : 999.0;

            if (ratio <= RATIO_CUT_FINAL) {
                // This segment is close enough to the winner (Ambiguous). Keep it.
                segmentsToKeep.insert(currSegIdx);
                if (_debugLevel) {
                    std::cout << std::format("    -> Seg {} is ambiguous (Ratio {:.2f} <= {}). Keeping hit.\n",
                                             currSegIdx, ratio, RATIO_CUT_FINAL);
                }
            } else {
                // This segment is significantly worse. Drop it.
                // (Note: It is not added to segmentsToKeep)
                if (_debugLevel) {
                    std::cout << std::format("    -> Seg {} is too poor (Ratio {:.2f} > {}). Marking for removal.\n",
                                             currSegIdx, ratio, RATIO_CUT_FINAL);
                }
            }
        }

        // 2. Execution of removal
        // Iterate over ALL segments that originally claimed this hit (segSet).
        // If a segment is NOT in 'segmentsToKeep', we remove the hit.
        // This handles:
        //   a. Segments that passed Chi2 but failed the Ratio cut.
        //   b. Segments that failed the initial Chi2 cut (never made it to goodCandidates).
        for (int segIdx : segSet) {
            if (segmentsToKeep.find(segIdx) == segmentsToKeep.end()) {

                // Remove hit from this segment
                auto& s = all_ThisIsBestSegment[segIdx];
                auto& d = all_ThisIsBestSegment_Diag[segIdx];

                s.erase(std::remove_if(s.begin(), s.end(), [&](auto& h){ return h.hitID == hitIdx; }), s.end());
                d.erase(std::remove_if(d.begin(), d.end(), [&](auto& x){ return x.reference_point == hitIdx; }), d.end());

                if (_debugLevel) {
                     // Optional verbose logging for removal
                     // std::cout << std::format("       [Removed hit {} from Seg {}]\n", hitIdx, segIdx);
                }
            }
        }

        /*if (goodCandidates.empty()) {
            // Remove hit from all segments if none are good fits
            for (int segIdx : segSet) {
                auto& s = all_ThisIsBestSegment[segIdx];
                auto& d = all_ThisIsBestSegment_Diag[segIdx];
                s.erase(std::remove_if(s.begin(), s.end(), [&](auto& h){ return h.hitID == hitIdx; }), s.end());
                d.erase(std::remove_if(d.begin(), d.end(), [&](auto& x){ return x.reference_point == hitIdx; }), d.end());
            }
            continue;
        }

        std::sort(goodCandidates.begin(), goodCandidates.end(), [](auto& a, auto& b){ return a.residual < b.residual; });

        int bestSegIdx = goodCandidates[0].segIdx;
        bool makeExclusive = (goodCandidates.size() == 1) || (goodCandidates[1].residual / goodCandidates[0].residual > RATIO_CUT_FINAL);

        for (int segIdx : segSet) {
            if (makeExclusive && segIdx != bestSegIdx) {
                // Remove from losers
                auto& s = all_ThisIsBestSegment[segIdx];
                auto& d = all_ThisIsBestSegment_Diag[segIdx];
                s.erase(std::remove_if(s.begin(), s.end(), [&](auto& h){ return h.hitID == hitIdx; }), s.end());
                d.erase(std::remove_if(d.begin(), d.end(), [&](auto& x){ return x.reference_point == hitIdx; }), d.end());
            }
        }*/

    }

    if (_debugLevel) {
        printCurrentSegments("06 FINAL STATE (AFTER MACRO & MICRO CLEANING)", all_ThisIsBestSegment);
    }

} // end ev5_select_best_segments_step_06

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
void PhiZSeedFinder::findchisq_ver2(const std::vector<ev5_HitsInNthStation>& segment) {
    // Use the class member fitter
    _lineFitter.clear();

    for (const auto& hit : segment) {
        double z = hit.z;
        double phi = hit.phi;

        // Define error/weight (currently fixed at 0.1, but should ideally come from hit resolution)
        double sigma = 0.1;
        double weight = 1.0 / (sigma * sigma);

        _lineFitter.addPoint(z, phi, weight);
    }
    // No return needed.
    // The alpha, beta, and chi2 are now stored in _lineFitter state.
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
// mergeSegmentsAll : "greedy first-match" -> "global best-pair" version
// All pairs are evaluated first, then only the most consistent pair is
// merged; the whole procedure is repeated until no pair is left.
// -------------------------------------------------------------------

//-----------------------------------------------------------------------------
// Pair evaluation: never modifies all_ThisIsBestSegment (side-effect free).
// Returns true if the pair is a valid merge candidate.
//
// Structure (one step per stage, in this order):
//   1. VETO-0  station overlap
//   2. circle fit of the two segments together
//   3. clean-up: drop hits until the circle chi2/ndf is acceptable
//   4. refit the circle, recompute helixPhi about the common centre
//   5. fit dphi/dz separately for ref and test
//   6. VETO-1  are the two slopes compatible?
//   7. VETO-2  extrapolate the ref line to the centre of the test segment and
//              check that the test segment, shifted by n*2pi as a whole, lands
//              on it within the combined error
//   8. RANKING (only reached when the pair is mergeable): refit the merged
//              segment as one line and store the numbers used to rank
//              candidates against each other
//
// There is deliberately no absolute cut on the merged chi2/ndf: VETO-1 and
// VETO-2 decide mergeability, and the merged chi2/ndf is used only to order
// the surviving candidates.
//
// Debug output is controlled by _debugLevel:
//    >0 : one summary line per pair, plus the reason for every rejection
//    >1 : full step-by-step dump (circle fit, clean-up, both slope fits,
//         every veto with its numbers and thresholds)
//    >2 : per-hit dump of the merged segment
//-----------------------------------------------------------------------------
bool PhiZSeedFinder::evaluateMergePair(const std::vector<ev5_HitsInNthStation>& segRef,
                                       const std::vector<ev5_Segment>&          diagRef,
                                       const std::vector<ev5_HitsInNthStation>& segTest,
                                       const std::vector<ev5_Segment>&          diagTest,
                                       int refIdx, int testIdx,
                                       ev5_MergeCandidate& cand) {

  cand = ev5_MergeCandidate();
  cand.refIdx     = refIdx;
  cand.testIdx    = testIdx;
  cand.nHitRefIn  = (int)segRef.size();
  cand.nHitTestIn = (int)segTest.size();

  const bool dbg  = (_debugLevel > 0);
  const bool dbg2 = (_debugLevel > 1);
  const bool dbg3 = (_debugLevel > 2);

  std::ostringstream tag;
  tag << "[evalPair " << std::setw(2) << refIdx << "-" << std::setw(2) << testIdx << "]";
  const std::string T = tag.str();

  auto reject = [&](const std::string& why) {
    cand.valid        = false;
    cand.rejectReason = why;
    if (dbg) std::cout << T << "  REJECT : " << why << std::endl;
    return false;
  };

  if (dbg2) {
    std::cout << "\n" << T << " ======================================================" << std::endl;
    std::cout << T << " nHit ref = " << segRef.size()
              << " , nHit test = " << segTest.size() << std::endl;
  }

  if (segRef.empty() || segTest.empty())
    return reject("empty segment");

  //---------------------------------------------------------------------------
  // Station / z ranges (context for the later messages)
  //---------------------------------------------------------------------------
  int    stRefMin = 9999, stRefMax = -1, stTstMin = 9999, stTstMax = -1;
  double zRefLo = 1e9, zRefHi = -1e9, zTstLo = 1e9, zTstHi = -1e9;
  for (const auto& h : segRef) {
    stRefMin = std::min(stRefMin, h.station);  stRefMax = std::max(stRefMax, h.station);
    zRefLo   = std::min(zRefLo,   h.z);        zRefHi   = std::max(zRefHi,   h.z);
  }
  for (const auto& h : segTest) {
    stTstMin = std::min(stTstMin, h.station);  stTstMax = std::max(stTstMax, h.station);
    zTstLo   = std::min(zTstLo,   h.z);        zTstHi   = std::max(zTstHi,   h.z);
  }
  if (dbg2) {
    std::cout << std::fixed << std::setprecision(1)
              << T << " ref  : station [" << stRefMin << "," << stRefMax << "]"
              << "  z [" << zRefLo << "," << zRefHi << "]" << std::endl;
    std::cout << T << " test : station [" << stTstMin << "," << stTstMax << "]"
              << "  z [" << zTstLo << "," << zTstHi << "]" << std::endl;
  }

  //===========================================================================
  // STEP 1 : VETO-0 - station overlap
  //          Two segments sharing a station cannot come from the same track.
  //===========================================================================
  {
    std::vector<int> shared;
    for (const auto& hR : segRef)
      for (const auto& hT : segTest)
        if (hR.station == hT.station &&
            std::find(shared.begin(), shared.end(), hR.station) == shared.end())
          shared.push_back(hR.station);

    if (!shared.empty()) {
      std::ostringstream m;
      m << "VETO-0 station overlap, n=" << shared.size() << " station(s) {";
      for (size_t k = 0; k < shared.size(); ++k) m << (k ? "," : "") << shared[k];
      m << "}";
      return reject(m.str());
    }
    if (dbg2) std::cout << T << " STEP1 VETO-0 station overlap : PASS" << std::endl;
  }

  // Local copies: everything below operates on temporaries only.
  std::vector<ev5_HitsInNthStation> locRef  = segRef;
  std::vector<ev5_HitsInNthStation> locTest = segTest;
  std::vector<ev5_Segment>          locDRef = diagRef;
  std::vector<ev5_Segment>          locDTest= diagTest;

  for (auto& h : locRef ) { h.used = true; h.nturn = 0; }
  for (auto& h : locTest) { h.used = true; h.nturn = 0; }

  //===========================================================================
  // STEP 2 : circle fit of both segments together
  //          (flat weights first, then proper weights)
  //===========================================================================
  _circleFitter.clear();
  for (const auto& h : locRef ) _circleFitter.addPoint(h.x, h.y, 0.1);
  for (const auto& h : locTest) _circleFitter.addPoint(h.x, h.y, 0.1);
  double xC = _circleFitter.x0();
  double yC = _circleFitter.y0();
  double rC = _circleFitter.radius();

  if (dbg2)
    std::cout << std::fixed << std::setprecision(3)
              << T << " STEP2 circle (flat weights)   : xC=" << xC << " yC=" << yC
              << " rC=" << rC << " chi2/ndf=" << _circleFitter.chi2DofCircle() << std::endl;

  std::vector<ev5_HitsInNthStation> allHits;
  allHits.insert(allHits.end(), locRef.begin(),  locRef.end());
  allHits.insert(allHits.end(), locTest.begin(), locTest.end());

  _circleFitter.clear();
  for (auto& h : allHits) {
    h.circleError2 = computeCircleError2_ver2(h.hitIndice, h.strawhits, xC, yC, rC);
    _circleFitter.addPoint(h.x, h.y, 1.0 / h.circleError2);
  }
  xC = _circleFitter.x0();
  yC = _circleFitter.y0();
  rC = _circleFitter.radius();

  if (dbg2)
    std::cout << std::fixed << std::setprecision(3)
              << T << " STEP2 circle (proper weights) : xC=" << xC << " yC=" << yC
              << " rC=" << rC << " chi2/ndf=" << _circleFitter.chi2DofCircle()
              << " nPoint=" << _circleFitter.qn() << std::endl;

  //===========================================================================
  // STEP 3 : clean-up - drop hits one at a time while the circle chi2/ndf
  //          keeps improving
  //===========================================================================
  int nDropped = 0;
  if (_circleFitter.qn() > 10 && _circleFitter.chi2DofCircle() > _mergeMaxChi2NDF) {
    double chi2ndf = _circleFitter.chi2DofCircle();
    if (dbg2)
      std::cout << T << " STEP3 clean-up starts : chi2/ndf=" << chi2ndf
                << " > " << _mergeMaxChi2NDF << std::endl;
    while (chi2ndf > _mergeMaxChi2NDF) {
      int    bestRemove = -1;
      double bestChi2   = chi2ndf;
      for (size_t i = 0; i < allHits.size(); ++i) {
        if (!allHits[i].used) continue;
        _circleFitter.clear();
        for (size_t j = 0; j < allHits.size(); ++j) {
          if (!allHits[j].used || i == j) continue;
          _circleFitter.addPoint(allHits[j].x, allHits[j].y, 1.0 / allHits[j].circleError2);
        }
        if (_circleFitter.chi2DofCircle() < bestChi2) {
          bestChi2   = _circleFitter.chi2DofCircle();
          bestRemove = (int)i;
        }
      }
      if (bestRemove < 0) {
        if (dbg2)
          std::cout << T << " STEP3 clean-up stops : no single hit improves chi2/ndf ("
                    << chi2ndf << ")" << std::endl;
        break;
      }
      allHits[bestRemove].used = false;
      ++nDropped;
      if (dbg2)
        std::cout << std::fixed << std::setprecision(3)
                  << T << " STEP3   drop hit " << allHits[bestRemove].hitIndice
                  << " (station " << allHits[bestRemove].station
                  << ", from " << (bestRemove < (int)locRef.size() ? "ref" : "test")
                  << ") chi2/ndf " << chi2ndf << " -> " << bestChi2 << std::endl;
      _circleFitter.clear();
      for (auto& h : allHits) {
        if (!h.used) continue;
        h.circleError2 = computeCircleError2_ver2(h.hitIndice, h.strawhits, xC, yC, rC);
        _circleFitter.addPoint(h.x, h.y, 1.0 / h.circleError2);
      }
      chi2ndf = _circleFitter.chi2DofCircle();
      xC = _circleFitter.x0();
      yC = _circleFitter.y0();
      rC = _circleFitter.radius();
      if (chi2ndf < _mergeMaxChi2NDF) {
        if (dbg2) std::cout << T << " STEP3 clean-up done : chi2/ndf=" << chi2ndf << std::endl;
        break;
      }
      if (_circleFitter.qn() <= 10) {
        if (dbg2) std::cout << T << " STEP3 clean-up stops : only " << _circleFitter.qn()
                            << " hits left" << std::endl;
        break;
      }
    }
  }
  else if (dbg2) {
    std::cout << T << " STEP3 clean-up skipped (nPoint=" << _circleFitter.qn()
              << ", chi2/ndf=" << _circleFitter.chi2DofCircle() << ")" << std::endl;
  }

  const int nUsed = (int)_circleFitter.qn();
  cand.nDropped   = nDropped;
  if (nUsed < 5) {
    std::ostringstream m;
    m << "too few hits survive the circle clean-up (" << nUsed << " < 5, dropped " << nDropped << ")";
    return reject(m.str());
  }

  //===========================================================================
  // STEP 4 : final circle parameters, then recompute helixPhi about the
  //          common circle centre
  //===========================================================================
  cand.xC = xC;  cand.yC = yC;  cand.rC = rC;
  cand.ndf_circle     = std::max(1.0, (double)(nUsed - 3));
  cand.chi2ndf_circle = _circleFitter.chi2DofCircle();
  cand.chi2_circle    = cand.chi2ndf_circle * cand.ndf_circle;
  cand.fDropped       = (double)nDropped / (double)allHits.size();

  if (dbg2)
    std::cout << std::fixed << std::setprecision(3)
              << T << " STEP4 circle result : nUsed=" << nUsed << " nDropped=" << nDropped
              << " chi2/ndf_circle=" << cand.chi2ndf_circle
              << " R=" << rC << std::endl;

  // Remove the dropped hits from the local segments.
  auto pruneSegment = [&](std::vector<ev5_HitsInNthStation>& seg,
                          std::vector<ev5_Segment>& diag) {
    auto itS = seg.begin();
    auto itD = diag.begin();
    while (itS != seg.end()) {
      bool drop = false;
      for (const auto& h : allHits)
        if (!h.used && h.hitIndice == itS->hitIndice) { drop = true; break; }
      if (drop) {
        itS = seg.erase(itS);
        if (itD != diag.end()) itD = diag.erase(itD);
      } else {
        ++itS;
        if (itD != diag.end()) ++itD;
      }
    }
  };
  pruneSegment(locRef,  locDRef);
  pruneSegment(locTest, locDTest);

  cand.nHitRef  = (int)locRef.size();
  cand.nHitTest = (int)locTest.size();

  if (locRef.size() < 2 || locTest.size() < 2) {
    std::ostringstream m;
    m << "segment shrank below 2 hits after clean-up (ref=" << locRef.size()
      << ", test=" << locTest.size() << ")";
    return reject(m.str());
  }

  auto recomputeHelixPhi = [&](std::vector<ev5_HitsInNthStation>& seg) {
    for (auto& h : seg) {
      double hPhi = 0.0, hPhiErr2 = 0.0;
      computeHelixPhi_ver2(h.hitIndice, xC, yC, hPhi, hPhiErr2);
      h.helixPhi       = hPhi;
      h.helixPhiError2 = hPhiErr2;
    }
  };
  recomputeHelixPhi(locRef);
  recomputeHelixPhi(locTest);

  if (dbg2)
    std::cout << std::fixed << std::setprecision(4)
              << T << " STEP4 helixPhi recomputed about (" << xC << "," << yC << ")"
              << std::endl;

  //===========================================================================
  // STEP 5 : fit dphi/dz separately for ref and test
  //          ev5_fit_slope_ver5 reads the member _segmentHits by index, so a
  //          scratch container is swapped in and restored afterwards.
  //===========================================================================
  std::vector<std::vector<ev5_HitsInNthStation>> saved_segmentHits = _segmentHits;

  auto sortByZ = [](std::vector<ev5_HitsInNthStation>& v) {
    std::sort(v.begin(), v.end(),
              [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b) { return a.z < b.z; });
  };
  sortByZ(locRef);
  sortByZ(locTest);

  _segmentHits.clear();
  _segmentHits.push_back(locRef);    // index 0
  _segmentHits.push_back(locTest);   // index 1

  double aR = 0, aRe = 0, bR = 0, bRe = 0, cR = 0;
  double aT = 0, aTe = 0, bT = 0, bTe = 0, cT = 0;
  ev5_fit_slope_ver5(0, aR, aRe, bR, bRe, cR);
  ev5_fit_slope_ver5(1, aT, aTe, bT, bTe, cT);

  cand.alphaRef  = aR;  cand.alphaRefErr  = aRe;  cand.chi2ndfRef  = cR;
  cand.alphaTest = aT;  cand.alphaTestErr = aTe;  cand.chi2ndfTest = cT;

  if (dbg2) {
    std::cout << std::scientific << std::setprecision(4)
              << T << " STEP5 ref  fit : alpha=" << aR << " +- " << aRe
              << "  beta=" << bR << " +- " << bRe
              << std::fixed << std::setprecision(3)
              << "  chi2/ndf=" << cR << std::endl;
    std::cout << std::scientific << std::setprecision(4)
              << T << " STEP5 test fit : alpha=" << aT << " +- " << aTe
              << "  beta=" << bT << " +- " << bTe
              << std::fixed << std::setprecision(3)
              << "  chi2/ndf=" << cT << std::endl;
  }

  //===========================================================================
  // STEP 6 : VETO-1 - are the two slopes compatible?
  //===========================================================================
  const double alpha_sigma = 0.0003199;                 // empirical value
  double alpha_diff = (aR * aT > 0) ? std::fabs(aR - aT)
                                    : std::fabs(aR) + std::fabs(aT);
  double combErrA = std::sqrt(aRe * aRe + aTe * aTe);
  if (combErrA < 1e-12) combErrA = 1e-12;
  double thrA = std::max(alpha_sigma * _mergeMaxZalpha, _mergeMaxZalpha * combErrA);
  cand.Zalpha    = alpha_diff / combErrA;
  cand.alphaDiff = alpha_diff;
  cand.alphaThr  = thrA;
  bool okAlpha   = (alpha_diff <= thrA);
  cand.okAlpha   = okAlpha;

  if (dbg2)
    std::cout << std::scientific << std::setprecision(4)
              << T << " STEP6 VETO-1 slope : |dAlpha|=" << alpha_diff
              << " thr=" << thrA
              << " (sigma-term=" << alpha_sigma * _mergeMaxZalpha
              << " , err-term=" << _mergeMaxZalpha * combErrA << ")"
              << std::fixed << std::setprecision(2)
              << "  pull=" << cand.Zalpha
              << "  -> " << (okAlpha ? "PASS" : "FAIL")
              << (aR * aT > 0 ? "" : "  [opposite sign slopes]") << std::endl;

  if (!okAlpha) {
    _segmentHits = saved_segmentHits;
    std::ostringstream m;
    m << "VETO-1 slope mismatch: |dAlpha|=" << std::scientific << std::setprecision(3)
      << alpha_diff << " > " << thrA << ", pull=" << std::fixed << std::setprecision(2)
      << cand.Zalpha;
    return reject(m.str());
  }

  //===========================================================================
  // STEP 7 : VETO-2 - extrapolate the ref line to the centre of the test
  //          segment and check that the test segment, shifted as a whole by
  //          n*2pi, lands on that line within the combined error.
  //
  //          zTest is the hit-density-weighted centre of the test segment, so
  //          the comparison is made where the test segment is best determined
  //          and the extrapolation from ref is shortest.
  //          Both sides are evaluated from their own line fits rather than
  //          from a single hit, so hit-level noise does not drive the test.
  //===========================================================================
  double zTest = 0.0;
  for (const auto& h : locTest) zTest += h.z;
  zTest /= (double)locTest.size();

  const double phiPred   = aR * zTest + bR;   // ref line extrapolated to zTest
  const double phiObsRaw = aT * zTest + bT;   // test line evaluated at zTest

  // Rigid n*2pi shift of the whole test segment onto the extended ref line.
  const int    nShift  = (int)std::round((phiPred - phiObsRaw) / (2.0 * M_PI));
  const double phiObs  = phiObsRaw + nShift * 2.0 * M_PI;
  const double dPhi    = phiPred - phiObs;

  const double errPred  = std::sqrt(zTest * zTest * aRe * aRe + bRe * bRe);
  const double errObs   = std::sqrt(zTest * zTest * aTe * aTe + bTe * bTe);
  double combErrP = std::sqrt(errPred * errPred + errObs * errObs);
  if (combErrP < 1e-12) combErrP = 1e-12;
  const double thrP = std::max(0.4, _mergeMaxZphi * combErrP);

  cand.Zphi            = std::fabs(dPhi) / combErrP;
  cand.dPhi            = dPhi;
  cand.dPhiThr         = thrP;
  bool okPhi           = (std::fabs(dPhi) <= thrP);
  cand.okPhi           = okPhi;
  cand.deltaCorrection = nShift;

  if (dbg2) {
    std::cout << std::fixed << std::setprecision(4)
              << T << " STEP7 VETO-2 phase : zTest=" << std::setprecision(1) << zTest
              << std::setprecision(4)
              << "  phiPred=" << phiPred << " (+-" << errPred << ")"
              << "  phiObsRaw=" << phiObsRaw << " (+-" << errObs << ")"
              << "  n=" << nShift << std::endl;
    std::cout << std::fixed << std::setprecision(4)
              << T << "                     phiObs(shifted)=" << phiObs
              << "  dPhi=" << dPhi
              << " thr=" << thrP
              << " (floor=0.4 , err-term=" << _mergeMaxZphi * combErrP << ")"
              << std::setprecision(2) << "  pull=" << cand.Zphi
              << "  -> " << (okPhi ? "PASS" : "FAIL") << std::endl;
  }

  if (!okPhi) {
    _segmentHits = saved_segmentHits;
    std::ostringstream m;
    m << "VETO-2 phase mismatch: |dPhi|=" << std::fixed << std::setprecision(4)
      << std::fabs(dPhi) << " > " << thrP << ", pull=" << std::setprecision(2)
      << cand.Zphi << ", n=" << nShift;
    return reject(m.str());
  }

  //===========================================================================
  // STEP 8 : the pair IS mergeable. Build the merged segment and produce the
  //          numbers used to rank this candidate against the other mergeable
  //          pairs. Ranking order:
  //            primary   - chi2/ndf of the merged phi-z fit after the n*2pi
  //                        shift (smaller is better)
  //            secondary - chi2/ndf of the circle fit, used only when the
  //                        primary values are effectively equal
  //          No cut is applied here; this stage only scores.
  //===========================================================================
  std::vector<ev5_HitsInNthStation> merged;
  merged.reserve(locRef.size() + locTest.size());
  for (auto h : locRef ) { h.nturn = 0;      merged.push_back(h); }
  for (auto h : locTest) { h.nturn = nShift; merged.push_back(h); }
  sortByZ(merged);                        // most upstream hit becomes the reference

  _segmentHits.clear();
  _segmentHits.push_back(merged);         // index 0
  double aM = 0, aMe = 0, bM = 0, bMe = 0, cM = 0;
  ev5_fit_slope_ver5(0, aM, aMe, bM, bMe, cM);

  _segmentHits = saved_segmentHits;       // restore the member

  cand.alpha        = aM;  cand.alphaError = aMe;
  cand.beta         = bM;  cand.betaError  = bMe;
  cand.ndf_phiz     = std::max(1.0, (double)merged.size() - 2.0);
  cand.chi2ndf_phiz = cM;
  cand.chi2_phiz    = cM * cand.ndf_phiz;

  if (dbg2)
    std::cout << std::scientific << std::setprecision(4)
              << T << " STEP8 merged fit : alpha=" << aM << " +- " << aMe
              << "  beta=" << bM << " +- " << bMe
              << std::fixed << std::setprecision(3)
              << "  chi2/ndf=" << cM << " (nHit=" << merged.size() << ")" << std::endl;

  if (dbg3) {
    std::cout << T << " STEP8 merged hits :" << std::endl;
    std::cout << T << "   idx  station        z     helixPhi   nturn    resid" << std::endl;
    for (const auto& h : merged) {
      const double phiFit = aM * h.z + bM;
      const double resid  = (h.helixPhi + 2 * M_PI * h.nturn) - phiFit;
      std::cout << std::fixed
                << T << "  " << std::setw(5) << h.hitIndice
                << std::setw(9)  << h.station
                << std::setw(10) << std::setprecision(1) << h.z
                << std::setw(12) << std::setprecision(4) << h.helixPhi
                << std::setw(8)  << h.nturn
                << std::setw(10) << std::setprecision(4) << resid << std::endl;
    }
  }

  cand.okFit = true;     // no absolute quality cut at this stage
  cand.valid = true;

  cand.mergedHits = merged;
  cand.mergedDiag = locDRef;
  cand.mergedDiag.insert(cand.mergedDiag.end(), locDTest.begin(), locDTest.end());
  std::sort(cand.mergedDiag.begin(), cand.mergedDiag.end(),
            [](const ev5_Segment& a, const ev5_Segment& b) { return a.z < b.z; });

  if (dbg)
    std::cout << std::fixed << std::setprecision(4)
              << T << "  ACCEPT : chi2/ndf_phiz=" << cand.chi2ndf_phiz
              << " (primary)  chi2/ndf_circle=" << cand.chi2ndf_circle
              << " (secondary)  Zalpha=" << std::setprecision(2) << cand.Zalpha
              << " Zphi=" << cand.Zphi
              << " n=" << cand.deltaCorrection
              << " nHit=" << cand.mergedHits.size() << std::endl;

  return true;
}


//-----------------------------------------------------------------------------
// (3) New mergeSegmentsAll: evaluate all pairs -> merge only the best one ->
//     repeat until no mergeable pair is left.
//-----------------------------------------------------------------------------
void PhiZSeedFinder::mergeSegmentsAll(std::vector<std::vector<ev5_HitsInNthStation>>& all_ThisIsBestSegment,
                                      std::vector<std::vector<ev5_Segment>>&          all_ThisIsBestSegment_Diag,
                                      double thre_residual) {

  (void)thre_residual;   // kept for signature compatibility, not used any more

  const bool dbg  = 1;

  std::cout << "\n############################################################" << std::endl;
  std::cout << "# mergeSegmentsAll (best-pair version)  : start with "
            << all_ThisIsBestSegment.size() << " segment(s)" << std::endl;
  std::cout << "############################################################" << std::endl;

  // Initialisation: stamp the original segment number on every hit, since
  // ev5_fit_slope_ver4 groups hits by segmentIndex.
  for (size_t s = 0; s < all_ThisIsBestSegment.size(); ++s) {
    for (auto& h : all_ThisIsBestSegment[s]) {
      h.nturn        = 0;
      h.segmentIndex = (int)s;
      h.used         = true;
    }
  }

  if (dbg) {
    std::cout << "# input segments:" << std::endl;
    for (size_t s = 0; s < all_ThisIsBestSegment.size(); ++s) {
      int    stLo = 9999, stHi = -1;
      double zLo = 1e9, zHi = -1e9;
      for (const auto& h : all_ThisIsBestSegment[s]) {
        stLo = std::min(stLo, h.station);  stHi = std::max(stHi, h.station);
        zLo  = std::min(zLo,  h.z);        zHi  = std::max(zHi,  h.z);
      }
      std::cout << std::fixed << std::setprecision(1)
                << "#   seg " << std::setw(2) << s
                << " : nHit=" << std::setw(3) << all_ThisIsBestSegment[s].size()
                << "  station [" << stLo << "," << stHi << "]"
                << "  z [" << zLo << "," << zHi << "]" << std::endl;
    }
  }

  _segmentHits = all_ThisIsBestSegment;

  if (all_ThisIsBestSegment.size() < 2) {
    std::cout << "# only one segment, nothing to merge" << std::endl;
    std::cout << "Final #segments = " << all_ThisIsBestSegment.size() << std::endl;
    return;
  }

  int pass       = 0;
  int nMergeDone = 0;

  while (all_ThisIsBestSegment.size() >= 2) {

    ++pass;
    const int nSeg  = (int)all_ThisIsBestSegment.size();
    const int nPair = nSeg * (nSeg - 1) / 2;
    std::cout << "\n=========== merge pass " << pass
              << " : " << nSeg << " segment(s), " << nPair << " pair(s) to test"
              << " ===========" << std::endl;

    //-------------------------------------------------------------------------
    // PHASE A : evaluate every pair (nothing is modified at this point)
    //-------------------------------------------------------------------------
    std::vector<ev5_MergeCandidate> evaluated;   // every pair, accepted or not
    evaluated.reserve(nPair);

    for (size_t i = 0; i + 1 < all_ThisIsBestSegment.size(); ++i) {
      for (size_t j = i + 1; j < all_ThisIsBestSegment.size(); ++j) {
        ev5_MergeCandidate cand;
        evaluateMergePair(all_ThisIsBestSegment[i], all_ThisIsBestSegment_Diag[i],
                          all_ThisIsBestSegment[j], all_ThisIsBestSegment_Diag[j],
                          (int)i, (int)j, cand);
        evaluated.push_back(cand);
      }
    }

    //-------------------------------------------------------------------------
    // Summary table of the whole pass: one line per pair, accepted or rejected
    //-------------------------------------------------------------------------
    int nAccept = 0;
    for (const auto& c : evaluated) if (c.valid) ++nAccept;

    std::cout << "--- pass " << pass << " summary : "
              << nAccept << " accepted / " << evaluated.size() << " pairs ---" << std::endl;
    std::cout << "  pair    nHit(r,t)  drop  chi2/ndf_phiz  chi2/ndf_cir   Zalpha    Zphi   n   verdict"
              << std::endl;
    for (const auto& c : evaluated) {
      std::cout << std::fixed
                << "  " << std::setw(2) << c.refIdx << "-" << std::setw(2) << c.testIdx
                << "   (" << std::setw(3) << c.nHitRefIn << "," << std::setw(3) << c.nHitTestIn << ")"
                << std::setw(6) << c.nDropped;
      if (c.valid) {
        std::cout << std::setprecision(3)
                  << std::setw(15) << c.chi2ndf_phiz
                  << std::setw(14) << c.chi2ndf_circle
                  << std::setprecision(2)
                  << std::setw(9)  << c.Zalpha
                  << std::setw(8)  << c.Zphi
                  << std::setw(4)  << c.deltaCorrection;
      } else {
        std::cout << std::setw(15) << "-" << std::setw(14) << "-"
                  << std::setw(9)  << "-" << std::setw(8) << "-" << std::setw(4) << "-";
      }
      std::cout << "   " << (c.valid ? std::string("ACCEPT") : ("REJECT: " + c.rejectReason))
                << std::endl;
    }

    //-------------------------------------------------------------------------
    // PHASE B : pick the single most consistent pair
    //-------------------------------------------------------------------------
    std::vector<ev5_MergeCandidate> candidates;
    for (const auto& c : evaluated) if (c.valid) candidates.push_back(c);

    if (candidates.empty()) {
      std::cout << "[mergeSegmentsAll] no mergeable pair left in pass " << pass
                << " -> stop" << std::endl;
      break;
    }

    // Ranking: the merged phi-z chi2/ndf after the n*2pi shift decides.
    // The circle chi2/ndf is consulted only when the phi-z values are
    // effectively equal, so it acts purely as a tie-breaker.
    const double kChi2TieEpsilon = 1.0e-3;
    std::sort(candidates.begin(), candidates.end(),
              [kChi2TieEpsilon](const ev5_MergeCandidate& a, const ev5_MergeCandidate& b) {
                if (std::fabs(a.chi2ndf_phiz - b.chi2ndf_phiz) > kChi2TieEpsilon)
                  return a.chi2ndf_phiz < b.chi2ndf_phiz;
                return a.chi2ndf_circle < b.chi2ndf_circle;
              });

    std::cout << "--- accepted candidates ranked (primary: chi2/ndf_phiz, tie-break: chi2/ndf_circle)"
              << " (pass " << pass << ") ---" << std::endl;
    for (size_t k = 0; k < candidates.size(); ++k) {
      const auto& c = candidates[k];
      std::cout << std::fixed << std::setprecision(4)
                << "  #" << k << " (" << c.refIdx << "," << c.testIdx << ")"
                << " chi2ndf_phiz=" << std::setw(9) << c.chi2ndf_phiz
                << " chi2ndf_cir="  << std::setw(9) << c.chi2ndf_circle
                << " | Zalpha=" << std::setprecision(2) << std::setw(6) << c.Zalpha
                << " Zphi="     << std::setw(6) << c.Zphi
                << " n="        << std::setw(3) << c.deltaCorrection
                << " nHitMerged=" << std::setw(4) << c.mergedHits.size()
                << (k == 0 ? "   <== WINNER" : "")
                << std::endl;
    }

    const ev5_MergeCandidate& best = candidates.front();

    if (candidates.size() > 1) {
      const auto& second = candidates[1];
      std::cout << std::fixed << std::setprecision(4)
                << "  margin over runner-up (" << second.refIdx << "," << second.testIdx << ") : "
                << "d(chi2/ndf_phiz)=" << (second.chi2ndf_phiz - best.chi2ndf_phiz)
                << "  d(chi2/ndf_circle)=" << (second.chi2ndf_circle - best.chi2ndf_circle)
                << std::endl;
    }

    //-------------------------------------------------------------------------
    // PHASE C : commit only the winning pair
    //-------------------------------------------------------------------------
    std::cout << "  --> MERGE segment " << best.testIdx << " into " << best.refIdx
              << " : nHit " << best.nHitRefIn << " + " << best.nHitTestIn
              << " - " << best.nDropped << " dropped = " << best.mergedHits.size()
              << std::fixed << std::setprecision(4)
              << " , chi2/ndf_phiz=" << best.chi2ndf_phiz
              << " , chi2/ndf_circle=" << best.chi2ndf_circle
              << " , n=" << best.deltaCorrection
              << std::scientific << std::setprecision(4)
              << " , alpha=" << best.alpha << " +- " << best.alphaError
              << std::fixed << std::setprecision(2)
              << " , R=" << best.rC << std::endl;

    all_ThisIsBestSegment[best.refIdx]      = best.mergedHits;   // helixPhi/nturn already updated
    all_ThisIsBestSegment_Diag[best.refIdx] = best.mergedDiag;

    all_ThisIsBestSegment.erase     (all_ThisIsBestSegment.begin()      + best.testIdx);
    all_ThisIsBestSegment_Diag.erase(all_ThisIsBestSegment_Diag.begin() + best.testIdx);
    ++nMergeDone;

    // Keep everything ordered in increasing z
    for (auto& seg : all_ThisIsBestSegment)
      std::sort(seg.begin(), seg.end(),
                [](const ev5_HitsInNthStation& a, const ev5_HitsInNthStation& b) { return a.z < b.z; });
    for (auto& seg : all_ThisIsBestSegment_Diag)
      std::sort(seg.begin(), seg.end(),
                [](const ev5_Segment& a, const ev5_Segment& b) { return a.z < b.z; });

    _segmentHits = all_ThisIsBestSegment;

    if (dbg) {
      std::cout << "  segments after pass " << pass << ":" << std::endl;
      for (size_t s = 0; s < all_ThisIsBestSegment.size(); ++s) {
        int    stLo = 9999, stHi = -1;
        double zLo = 1e9, zHi = -1e9;
        for (const auto& h : all_ThisIsBestSegment[s]) {
          stLo = std::min(stLo, h.station);  stHi = std::max(stHi, h.station);
          zLo  = std::min(zLo,  h.z);        zHi  = std::max(zHi,  h.z);
        }
        std::cout << std::fixed << std::setprecision(1)
                  << "    seg " << std::setw(2) << s
                  << " : nHit=" << std::setw(3) << all_ThisIsBestSegment[s].size()
                  << "  station [" << stLo << "," << stHi << "]"
                  << "  z [" << zLo << "," << zHi << "]" << std::endl;
      }
    }

    // Next pass: re-evaluate all pairs with one fewer segment
  }

  _segmentHits = all_ThisIsBestSegment;

  std::cout << "\n############################################################" << std::endl;
  std::cout << "# mergeSegmentsAll done : " << nMergeDone << " merge(s) in "
            << pass << " pass(es)" << std::endl;
  std::cout << "Final #segments = " << all_ThisIsBestSegment.size() << std::endl;
  std::cout << "############################################################" << std::endl;
}


//----------------------------------------------------------------------------
void PhiZSeedFinder::plot_PhiVsZ_forSegment_debug(
    int ith_segment,
    int jth_segment,
    double alpha,
    double beta,
    double Chi2NDF,
    int plotID,
    const std::string& stage,
    bool merged
) {
    TGraphErrors* gr = new TGraphErrors();
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetMarkerColor(merged ? kGreen+2 : kRed);
    std::cout<<"plot_PhiVsZ_forSegment_debug"<<std::endl;
    std::cout<<"segment/plotID = "<<ith_segment<<"/"<<plotID<<std::endl;
    std::cout<<"segment.size = "<<_SSegmentHits.at(ith_segment).size()<<std::endl;

    for (size_t j = 0; j < _SSegmentHits.at(ith_segment).size(); j++) {
        const auto& hit = _SSegmentHits.at(ith_segment).at(j);
        double z   = hit.z;
        double phi = hit.helixPhi + hit.nturn * 2 * M_PI;
        gr->SetPoint(j, z, phi);
        std::cout<<"j: "<<j<<", z/phi = "<<z<<"/"<<phi<<std::endl;
        gr->SetPointError(j, 0.0, std::sqrt(hit.helixPhiError2));
    }

    /*for (size_t j = 0; j < _tcHits.size(); j++) {
        if(_tcHits[j].used == false) continue;
        double z   = _tcHits[j].z;
        double phi = _tcHits[j].helixPhi;
        gr->SetPoint(j, z, phi);
        gr->SetPointError(j, 0.0, std::sqrt(_tcHits[j].helixPhiError2));
    }*/

    TCanvas* canvas = new TCanvas(Form("c_%d", plotID), "", 800, 600);
    canvas->SetMargin(0.1, 0.1, 0.1, 0.1);

    gr->Draw("AP");
    gr->GetXaxis()->SetTitle("Z [mm]");
    gr->GetYaxis()->SetTitle("Phi [rad]");
    gr->GetXaxis()->SetLimits(-1600, 1600);
    gr->GetYaxis()->SetRangeUser(-5*M_PI, 5*M_PI);

    // ===== Title =====
    TPaveText* title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
    title->SetFillColor(0);
    title->SetTextAlign(22);
    title->AddText(Form(
        "Phi vs Z | %s | plotID=%d | seg %d vs %d | merged=%s",
        stage.c_str(),
        plotID,
        ith_segment,
        jth_segment,
        merged ? "YES" : "NO"
    ));
    title->Draw("same");

    // ===== Fit line =====
    double x_min = -1600, x_max = 1600;
    TLine* fitLine = new TLine(
        x_min, alpha * x_min + beta,
        x_max, alpha * x_max + beta
    );
    fitLine->SetLineStyle(2);
    fitLine->SetLineWidth(2);
    fitLine->Draw("same");

    // ===== Params =====
    TLatex text;
    text.SetTextSize(0.035);
    text.DrawLatexNDC(0.15, 0.85, Form("alpha = %.6e", alpha));
    text.DrawLatexNDC(0.15, 0.80, Form("beta  = %.6e", beta));
    text.DrawLatexNDC(0.15, 0.75, Form("#chi^{2}/ndf = %.3f", Chi2NDF));

    canvas->SaveAs(Form(
        "/exp/mu2e/data/users/kitagawa/output/20240424/"
        "PhiZSeedFinder/pbar/debug/PhiVsZ_debug_%06d.pdf",
        plotID
    ));

    delete fitLine;
    delete title;
    delete canvas;
    delete gr;
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
  _tcHits.reserve(_segmentHits.at(isegment).size());

  for (size_t i = 0; i < _segmentHits.at(isegment).size(); ++i) {
    const auto& sh = _segmentHits.at(isegment).at(i);

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
  hit.nturn = 0;

    // ---- Copy EVERYTHING that exists in _segmentHits element (ev5_HitsInNthStation) ----
    hit.hitIndice       = sh.hitIndice;
    hit.phi             = sh.phi;
    hit.strawhits        = sh.strawhits;

    hit.x               = sh.x;
    hit.y               = sh.y;
    hit.z               = sh.z;

    hit.station         = sh.station;
    hit.plane           = sh.plane;
    hit.face            = sh.face;
    hit.panel           = sh.panel;

    //hit.hitID           = sh.hitID;
    hit.segmentIndice   = sh.segmentIndex;   // from _segmentHits
    hit.used            = sh.used;
    hit.nturn           = sh.nturn;

    hit.circleError2    = sh.circleError2;

    hit.helixPhi        = sh.helixPhi;
    hit.helixPhiError2  = sh.helixPhiError2;

    hit.used = true;
   //hit.strawhits = _segmentHits.at(isegment).at(i).strawhits;
   const ComboHit* ch = &_data.chcol->at(_segmentHits.at(isegment).at(i).hitIndice);
   hit.strawhits = ch->_nsh;
    _tcHits.push_back(hit);
  }

  // Sort by z ascending (as you already do)
  std::sort(_tcHits.begin(), _tcHits.end(),
            [](const cHit& a, const cHit& b) { return a.z < b.z; });

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
//----------------------------------------------------------------------------
void PhiZSeedFinder::plot_PhiVsZ_RawStep(const std::vector<std::vector<ev5_HitsInNthStation>>& segments,
                                         const std::string& stepName,
                                         int tc,
                                         int station) {



    int maxSegs = segments.size();
    for (int i = 0; i < maxSegs; ++i) {
        const auto& seg = segments[i];
        if(seg.empty()) continue;

        TGraphErrors* gr = new TGraphErrors();
        gr->SetTitle("");
        gr->SetMarkerStyle(20);

        int index = 0;
        for (const auto& hit : seg) {
          gr->SetPoint(index, hit.z, hit.phi);
          index++;
        }

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
        title->AddText(Form("Phi vs. Z (Run-subRun-Event, #cand) = (%d-%d-%d, #%d) (pbar1b0)", run, subrun, eventNumber, tc));
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
        std::string filename = std::format("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/segment_check/{}/PhiZ_TC{:03d}_Stn{:02d}_{:02d}_{}.pdf", stepName, tc, station, i, stepName);
          canvas->SaveAs(filename.c_str());

          //delete
          delete canvas;
          delete gr;
          for(int j=0; j<18; j++) delete line[j];

        }


}

//-----------------------------------------------------------------------------
void PhiZSeedFinder::plot_HelixPhiVsZ(int TC, int isegment){
    // step2: calculate the slope for other segment
    ::LsqSums2 _lineFitter;
    _lineFitter.clear();

    // We only need one TGraph now, not an array or TMultiGraph
    TGraph* gr = new TGraph();
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.5);
    gr->SetMarkerColor(kRed);

    int index = 0;
    for(size_t j=0; j<_tcHits.size(); j++) {
        if(_tcHits[j].used == false) continue;

        double z = _tcHits[j].z;
        double phi = _tcHits[j].phi;
        double helixPhi = _tcHits[j].helixPhi;

        // Only set point for helixPhi
        gr->SetPoint(index++, z, helixPhi);

        double phiWeight = 1.0 / (_tcHits[j].helixPhiError2);
        _lineFitter.addPoint(z, phi, phiWeight);
    }

    // Draw
    TCanvas *canvas = new TCanvas("canvas", "My TGraph", 800, 600);
    canvas->SetMargin(0.1, 0.1, 0.1, 0.1);

    // Draw the single graph directly
    gr->Draw("AP");
    gr->GetXaxis()->SetTitle("Z [mm]");
    gr->GetYaxis()->SetTitle("Helix Phi [rad]");
    gr->GetXaxis()->SetLimits(-1600, 1600);
    gr->GetYaxis()->SetRangeUser(-4*M_PI, 4*M_PI);

    // Add title at the top
    TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
    std::stringstream eventStringStream;
    eventStringStream << "run: " << run << " subRun: " << subrun << " event: " << eventNumber;
    title->AddText(Form("Phi vs. Z (Run-subRun-Event, TC, Cand) = (%d-%d-%d, #%d, #%d) (pbar1b0)", run, subrun, eventNumber, TC, isegment));
    title->SetFillColor(0);
    title->SetTextAlign(22);
    title->Draw("same");

    // Draw vertical lines at specified X-coordinate of stations (18 stations)
    double stations_x[18] = {-1518.320, -1344.320, -1170.320, -996.320, -822.320, -648.320, -474.320, -300.320, -126.320, 47.680, 221.680, 395.680, 569.680, 743.680, 917.680, 1091.680, 1265.680, 1439.680};
    TLine *line[18];
    for(int j=0; j<18; j++){
        double x = stations_x[j];
        line[j] = new TLine(x, -4*M_PI, x, 4*M_PI);
        line[j]->SetLineStyle(2);  // Dashed line style
        line[j]->SetLineWidth(1);
        line[j]->SetLineColor(kBlack);
        line[j]->Draw("same");
    }

    // --- Save and Cleanup (Assuming you have a save block below this in your original code) ---
    // canvas->SaveAs(...);
    // delete canvas;
    // delete gr;
    // delete title;
    // for(int j=0; j<18; j++) delete line[j];

  // Draw the fitted slope line: y = dydx * x + y0 from x = -1600 to 1600
  /*double x_min = -1600;
  double x_max = 1600;
  double y_min = _lineFitter.dydx() * x_min + _lineFitter.y0();
  double y_max = _lineFitter.dydx() * x_max + _lineFitter.y0();
  TLine *fitLine = new TLine(x_min, y_min, x_max, y_max);
  fitLine->SetLineColor(kBlack);
  fitLine->SetLineStyle(2);  // Dashed line
  fitLine->SetLineWidth(2);
  fitLine->Draw("same");*/

// Draw manual fit parameters text
/*TLatex paramsText;
paramsText.SetTextSize(0.04);
paramsText.SetTextAlign(13);
double textX = 0.15;
double textY = 0.85;
paramsText.DrawLatexNDC(textX, textY,        Form("LsqSums2 fit dydx: %f", _lineFitter.dydx()));
paramsText.DrawLatexNDC(textX, textY - 0.05, Form("LsqSums2 fit y0: %f", _lineFitter.y0()));
paramsText.DrawLatexNDC(textX, textY - 0.10, Form("LsqSums2 fit #chi^{2}/ndf: %f", _lineFitter.chi2Dof()));
*/   canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/HelixPhi/pbar_PhiVsZ-%04d-%04d-%06d_TC-%d_Segemnt_%d.pdf", run, subrun, eventNumber, TC, isegment));
   //canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/ce/PhiVsZ/HelixPhi/pbar_PhiVsZ-%04d-%04d-%06d_TC-%d.pdf", run, subrun, eventNumber, TC));
  //delete
  delete canvas;
  delete gr;
  for(int j=0; j<18; j++) delete line[j];
}//end plot_HelixPhiVsZ

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

// Only Draw if filename is not "step_07"
if (strcmp(filename, "step_07") != 0) {
    // Using %.0f to show zero decimal places
    paramsText.DrawLatexNDC(textX, textY,         Form("LsqSums2 fit dydx: %.6f", fitter.dydx()));
    paramsText.DrawLatexNDC(textX, textY - 0.05, Form("LsqSums2 fit y0: %.6f",   fitter.y0()));
    paramsText.DrawLatexNDC(textX, textY - 0.10, Form("LsqSums2 fit #chi^{2}/ndf: %.2f", fitter.chi2Dof()));
    paramsText.DrawLatexNDC(textX, textY - 0.15, Form("nHits: %.0f",            fitter.qn()));
}

canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/segment_check/%s/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d-%d.pdf", filename, run, subrun, eventNumber, tc, i, loopIndex));


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
//----------------------------------------------------------------------------
void PhiZSeedFinder::plot_PhiVsZ_forSegment_ver3(int tc, int isegment){
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
   if(_tcHits[i].used == false) continue;
   //double z = _segmentHits.at(isegment).at(i).z;
   //double phi = _segmentHits.at(isegment).at(i).helixPhi;
   //int index = _segmentHits.at(isegment).at(i).hitIndice;
   double z = _tcHits[i].z;
   double phi = _tcHits[i].helixPhi;
   int index = _tcHits[i].hitIndice;
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
  gr->GetYaxis()->SetTitle("HelixPhi [rad]");
  gr->GetXaxis()->SetLimits(-1600, 1600);
  gr->GetYaxis()->SetRangeUser(-6*M_PI, 6*M_PI);

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
    line[j] = new TLine(x, -6*M_PI, x, 6*M_PI);
    line[j]->SetLineStyle(2);  // Dashed line style
    line[j]->SetLineWidth(1);
    line[j]->SetLineColor(kBlack);
    line[j]->Draw("same");
  }

  // Add text for the double value
  TLatex dphidx;
  dphidx.SetTextSize(0.04);
  dphidx.SetTextAlign(13);
  dphidx.DrawLatexNDC(0.15, 0.85, Form("#alpha : #beta : #chi^{2}/ndf: %f : %f : %f", _lineFitter.dydx(), _lineFitter.y0(), _lineFitter.chi2Dof()));

  // =======================================================================
  // NEW: Draw the fitted line
  // =======================================================================
  // Define a 1D function: [0]*x + [1] (which is slope*x + intercept)
  TF1 *fitLine = new TF1("fitLine", "[0]*x + [1]", -1600, 1600);
  fitLine->SetParameter(0, _lineFitter.dydx()); // Slope
  fitLine->SetParameter(1, _lineFitter.y0());   // Intercept

  fitLine->SetLineColor(kBlack);
  fitLine->SetLineStyle(2); // 2 is the dashed line style in ROOT
  fitLine->SetLineWidth(2); // Slightly thicker than station lines so it stands out
  fitLine->Draw("same");

  // Save Canvas
  canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/segment_check/step_10/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d.pdf", run, subrun, eventNumber, tc, isegment));

  // =======================================================================
  // Cleanup memory
  // =======================================================================
  delete canvas;
  delete gr;
  for(int j=0; j<18; j++) delete line[j];
  delete fitLine; // Delete the TF1 object

} //end plot_PhiVsZ_forSegment_ver3

//-----------------------------------------------------------------------------
void PhiZSeedFinder::plot_PhiVsZ_alignment_step(int tc, int isegment, int stepIdx, int currentSegIdx, double slope, double intercept) {
    // =======================================================================
    // 1. Create Canvas and Graph
    // =======================================================================
    TCanvas *canvas = new TCanvas("canvas", "Phi vs Z Alignment Step", 800, 600);
    TGraphErrors *gr = new TGraphErrors();

    int ptIdx = 0;
    for(size_t j = 0; j < _tcHits.size(); j++) {
        if(!_tcHits[j].used) continue;
        // Plot all valid hits with their current (possibly shifted) helixPhi
        gr->SetPoint(ptIdx, _tcHits[j].z, _tcHits[j].helixPhi);
        gr->SetPointError(ptIdx, 0.0, std::sqrt(_tcHits[j].helixPhiError2));
        ptIdx++;
    }

    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->Draw("AP");
    gr->GetXaxis()->SetTitle("Z [mm]");
    gr->GetYaxis()->SetTitle("HelixPhi [rad]");
    gr->GetXaxis()->SetLimits(-1600, 1600);
    gr->GetYaxis()->SetRangeUser(-6*M_PI, 6*M_PI);

    // =======================================================================
    // 2. Add Titles and Station Lines
    // =======================================================================
    TPaveText *title = new TPaveText(0.1, 0.92, 0.9, 0.98, "NDC");
    // Added Step Index and Segment Index to the title for clarity
    title->AddText(Form("Phi vs Z Align Step %d (Seg %d) | TC %d, #cand %d", stepIdx, currentSegIdx, tc, isegment));
    title->SetFillColor(0);
    title->SetTextAlign(22);
    title->Draw("same");

    double stations_x[18] = {-1518.320, -1344.320, -1170.320, -996.320, -822.320, -648.320, -474.320, -300.320, -126.320, 47.680, 221.680, 395.680, 569.680, 743.680, 917.680, 1091.680, 1265.680, 1439.680};
    TLine *line[18];
    for(int j=0; j<18; j++){
        double x = stations_x[j];
        line[j] = new TLine(x, -6*M_PI, x, 6*M_PI);
        line[j]->SetLineStyle(2);
        line[j]->SetLineWidth(1);
        line[j]->SetLineColor(kBlack);
        line[j]->Draw("same");
    }

    // =======================================================================
    // 3. Draw Predicted Fitting Line and Text
    // =======================================================================
    TLatex dphidx;
    dphidx.SetTextSize(0.04);
    dphidx.SetTextAlign(13);
    // Note: Chi2 is omitted here since this is the prediction from the reference segment
    dphidx.DrawLatexNDC(0.15, 0.85, Form("Ref Fit -> #alpha : #beta = %f : %f", slope, intercept));

    TF1 *fitLine = new TF1("fitLine", "[0]*x + [1]", -1600, 1600);
    fitLine->SetParameter(0, slope);
    fitLine->SetParameter(1, intercept);
    fitLine->SetLineColor(kBlack);
    fitLine->SetLineStyle(2);
    fitLine->SetLineWidth(2);
    fitLine->Draw("same");

    // =======================================================================
    // 4. Save and Cleanup (Notice the filename includes %02d for the step)
    // =======================================================================
    canvas->SaveAs(Form("/exp/mu2e/data/users/kitagawa/output/20240424/PhiZSeedFinder/pbar/PhiVsZ/segment_check/step_11/pbar_PhiVsZ-%04d-%04d-%06d_TC_%d-%d_AlignStep%02d.pdf", run, subrun, eventNumber, tc, isegment, stepIdx));

    delete title;
    delete fitLine;
    for(int j=0; j<18; j++) delete line[j];
    delete gr;
    delete canvas;
}//end plot_PhiVsZ_alignment_step

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
  void PhiZSeedFinder::findHelix(int tc, int isegment, HelixSeedCollection& HSColl, HelixSeed& Temp_HSeed) {
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
    } //end findHelix


//-----------------------------------------------------------------------------
  void PhiZSeedFinder::findHelix_ver2(int tc, int isegment, HelixSeedCollection& HSColl, HelixSeed& Temp_HSeed) {
    std::cout<<"findHelix_ver2"<<std::endl;
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

    //Step2: fit circle of i-th segment with correct weight
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    _circleFitter.clear();
    _lineFitter.clear();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      std::cout<<"i = "<<i<<std::endl;
      computeCircleError2(i, xC, yC, rC);
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 1.0 / (_tcHits[i].circleError2);
      _circleFitter.addPoint(x, y, wP);

      //_tcHits[i].used = true;
      int hitIndice = _tcHits[i].hitIndice;
      double helixPhi = 0.0;
      double helixPhiError2 = 0.0;
      computeHelixPhi_ver2(hitIndice, xC, yC, helixPhi, helixPhiError2);
      _tcHits[i].helixPhi = helixPhi;//recalculate helixphi from mergeSegments
      _tcHits[i].helixPhiError2 = helixPhiError2;

      double z = _tcHits[i].z;
      helixPhi = _tcHits[i].helixPhi;
      helixPhi = helixPhi + _tcHits[i].nturn * 2 * M_PI;
      double phiWeight = 1.0 / (_tcHits[i].helixPhiError2);
      _lineFitter.addPoint(z, helixPhi, phiWeight);
    }

    if (_diagLevel > 0) {
      std::cout<<"==================================="<<std::endl;
      std::cout<<"              Step2                "<<std::endl;
      std::cout<<"==================================="<<std::endl;
      std::cout<<"# of hits in Helix = "<<_circleFitter.qn()<<std::endl;
      std::cout<<"xC/yC = "<<_circleFitter.x0()<<"/"<<_circleFitter.y0()<<std::endl;
      std::cout<<"radius = "<<_circleFitter.radius()<<std::endl;
      std::cout<<"phi/dfdz/chi2DofC/chi2DofLineC = "<<_circleFitter.phi0()<<"/"<<_circleFitter.dfdz()<<"/"<<_circleFitter.chi2DofCircle()<<"/"<<_circleFitter.chi2DofLine()<<std::endl;
    }
    //hseed._hhits.setParent(_chColl->parent());
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();
    std::cout<<"r/f cent = "<<sqrt(xC*xC + yC*yC)<<"/"<<polyAtan2(yC, xC)<<std::endl;
    std::cout<<"_lineFitter.chi2D = "<< _lineFitter.chi2Dof() <<std::endl;

    Temp_HSeed._helix._radius = rC;
    Temp_HSeed._helix._rcent  = sqrt(xC*xC + yC*yC);
    Temp_HSeed._helix._fcent  = polyAtan2(yC, xC);
    //Temp_HSeed._helix._lambda = 1.0/_dphidz;
    //Temp_HSeed._helix._fz0    = fz0;
    Temp_HSeed._helix._helicity = Temp_HSeed._helix._lambda > 0 ? Helicity::poshel : Helicity::neghel;
    Temp_HSeed._helix._chi2dXY = _circleFitter.chi2DofCircle();
    //Temp_HSeed._helix._chi2dZPhi = _lineFitter.chi2Dof();
    std::cout.setf(std::ios::fixed);
std::cout << std::setprecision(4)
          << "[Temp_HSeed Helix]\n"
          << "  radius      = " << rC      << " [mm]\n"
          << "  rcent       = " << Temp_HSeed._helix._rcent       << " [mm]\n"
          << "  fcent(rad)  = " << Temp_HSeed._helix._fcent       << " [rad]\n"
          << "  center(x,y) = (" << xC << ", " << yC << ") [mm]\n"
          << "  dydx      = " <<  _lineFitter.dydx()  << " [mm/rad]\n"
          << "  lambda      = " <<  1.0/_lineFitter.dydx()  << " [mm/rad]\n"
          //<< "  fz0(rad)    = " << H._fz0         << " [rad]\n"
          //<< "  tanDip      = " << tanDip         << " (= lambda/|radius|)\n"
          //<< "  helicity    = " << (H._helicity==Helicity::poshel ? "poshel" : "neghel") << "\n"
          //<< "  chi2dXY     = " << H._chi2dXY
          << std::endl;
    } //end findHelix_ver2

//-----------------------------------------------------------------------------
  void PhiZSeedFinder::findHelix_ver3(int tc, int isegment) {
    std::cout<<"findHelix_ver3"<<std::endl;
    _circleFitter.clear();
    //Step1: fit circle of the segment with fixed weight
    size_t nComboHitsInSegment = _tcHits.size();
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      _tcHits[i].nturn = 0;
      double x = _tcHits.at(i).x;
      double y = _tcHits.at(i).y;
      double wP = 0.1;//tentative value
      _circleFitter.addPoint(x, y, wP);
    }
    double xC = _circleFitter.x0();
    double yC = _circleFitter.y0();
    double rC = _circleFitter.radius();

    // --- (1) Setup Debug Data Structure ---
    struct CircleDebugInfo {
        size_t hitIdx;
        double x;
        double y;
        double z;    // <-- Added Z coordinate
        double phi;
        double error2;
        double weight;
    };
    std::vector<CircleDebugInfo> debugData;
    if (_debugLevel > 0) {
        debugData.reserve(nComboHitsInSegment); // Pre-allocate memory for efficiency
    }

    //Step2: fit circle of the segment with corrected weight
    _circleFitter.clear();

    for(size_t i = 0; i < nComboHitsInSegment; i++) {
        if(_tcHits[i].used == false) continue;

        // This function prints its own debug info safely without mixing
        computeCircleError2(i, xC, yC, rC);

        double x = _tcHits.at(i).x;
        double y = _tcHits.at(i).y;
        double z = _tcHits.at(i).z;  // <-- Fetch Z coordinate
        double wP = 1.0 / (_tcHits[i].circleError2);

        _circleFitter.addPoint(x, y, wP);

        // Store debug info in the vector
        if (_debugLevel > 0) {
            double phi = polyAtan2(y, x);
            debugData.push_back({i, x, y, z, phi, _tcHits[i].circleError2, wP}); // <-- Store Z
        }
    }

    // --- Recalculate Parameters ---
    xC = _circleFitter.x0();
    yC = _circleFitter.y0();
    rC = _circleFitter.radius();

    // --- (2) Dump Formatted Table AFTER the Loop ---
    if (_debugLevel > 0) {
        std::cout << "\n" << std::string(110, '=') << "\n"; // Widened for the extra column
        std::cout << " [Circle Fitter] Weight & Hit Info (Step 2)\n";
        std::cout << std::string(110, '-') << "\n";

        // Header (Added Z)
        std::cout << std::format("{:<8} {:<12} {:<12} {:<12} {:<12} {:<15} {:<15}\n",
                                 "HitIdx", "X [mm]", "Y [mm]", "Z [mm]", "Phi [rad]", "Error2", "Weight (wP)");
        std::cout << std::string(110, '-') << "\n";

        // Rows (Added Z)
        for (const auto& data : debugData) {
            std::cout << std::format("{:<8} {:<12.4f} {:<12.4f} {:<12.4f} {:<12.4f} {:<15.6f} {:<15.6f}\n",
                                     data.hitIdx, data.x, data.y, data.z, data.phi, data.error2, data.weight);
        }

        // Summary
        std::cout << std::string(110, '-') << "\n";
        std::cout << std::format(" Final Circle Fit: xC = {:.4f}, yC = {:.4f}, rC = {:.4f}\n", xC, yC, rC);
        std::cout << std::format(" Chi2/DOF        : {:.4f}\n", _circleFitter.chi2DofCircle());
        std::cout << std::string(110, '=') << std::endl;
    }

    std::cout << "  rC = " << rC
    << "  rcent = " << sqrt(xC*xC + yC*yC)
    << "  fcent = " << polyAtan2(yC, xC)
    << "  chi2dXY = " << _circleFitter.chi2DofCircle()
    << std::endl;

    //Step3: calculate helixPhi and helixPhiError
    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;
      computeCircleError2(i, xC, yC, rC);//update circle error
      int hitIndice = _tcHits[i].hitIndice;
      double helixPhi = 0.0, helixPhiError2 = 0.0;
      computeHelixPhi_ver2(hitIndice, xC, yC, helixPhi, helixPhiError2);
      _tcHits[i].helixPhi = helixPhi;
      _tcHits[i].helixPhiError2 = helixPhiError2;
    }

    // --- (3) Dump Formatted Table AFTER Step 3 ---
    if (_debugLevel > 0) {
        std::cout << "\n" << std::string(60, '=') << "\n";
        std::cout << " [Helix Fitter] Computed Helix Phi (Step 3)\n";
        std::cout << std::string(60, '-') << "\n";

        // Header
        std::cout << std::format("{:<8} {:<12} {:<15} {:<15}\n",
                                 "HitIdx", "Z [mm]", "HelixPhi", "HelixPhiErr2");
        std::cout << std::string(60, '-') << "\n";

        // Rows
        for (size_t i = 0; i < nComboHitsInSegment; i++) {
            if (!_tcHits[i].used) continue;

            std::cout << std::format("{:<8} {:<12.4f} {:<15.6f} {:<15.6f}\n",
                                     i,
                                     _tcHits[i].z,
                                     _tcHits[i].helixPhi,
                                     _tcHits[i].helixPhiError2);
        }
        std::cout << std::string(60, '=') << std::endl;
    }


    _lineFitter.clear();

    // Step 1: find unique segmentIndex values
    std::set<int> uniqueSegmentIndices;
    for (const auto &hit : _tcHits) {
      uniqueSegmentIndices.insert(hit.segmentIndice);
    }

    std::cout << "Number of unique segmentIndex values = "
              << uniqueSegmentIndices.size() << std::endl;

    // Step 2: process each segmentIndex independently
    for (int segIdx : uniqueSegmentIndices) {

      // Step A: collect hit indices for this segmentIndex
      std::vector<size_t> hitIdx;
      for (size_t i = 0; i < nComboHitsInSegment; ++i) {
        if (!_tcHits[i].used) continue;
        if (_tcHits[i].segmentIndice == segIdx) {
          hitIdx.push_back(i);
        }
      }

      if (hitIdx.size() < 3) {
        if (_debugLevel > 0) {
          std::cout << "Segment " << segIdx
                    << ": not enough hits (" << hitIdx.size() << ")\n";
        }
        continue;
      }

      // Step B: sort hits by z
      std::sort(hitIdx.begin(), hitIdx.end(),
                [&](size_t a, size_t b) {
                  return _tcHits[a].z < _tcHits[b].z;
                });

      // Step C: Predictive Phi Unwrapping
      std::vector<double> phi_unwrapped;
      phi_unwrapped.reserve(hitIdx.size());

      _lineFitter.clear();

      // 1. Setup the Reference Hit (Most upstream hit, j=0)
      size_t refIdx = hitIdx[0];
      double refZ   = _tcHits[refIdx].z;
      double refPhi = _tcHits[refIdx].helixPhi;
      double refWeight = 1.0 / _tcHits[refIdx].helixPhiError2;

      phi_unwrapped.push_back(refPhi);
      _lineFitter.addPoint(refZ, refPhi, refWeight);

      if (_debugLevel > 0) {
        std::cout << "\n>>> STARTING PREDICTIVE UNWRAPPING TRACE <<<\n";
        std::cout << std::format("Init: Added Ref Hit[0] at Z = {:.4f}, rawPhi = {:.4f}\n",
                                 refZ, refPhi);
      }

      // 2. Loop over remaining hits to unwrap and fit dynamically
      for (size_t i = 1; i < hitIdx.size(); ++i) {
        size_t idx = hitIdx[i];

        double z      = _tcHits[idx].z;
        double rawPhi = _tcHits[idx].helixPhi;
        double weight = 1.0 / _tcHits[idx].helixPhiError2;

        double candidatePhi = rawPhi;

        if (_debugLevel > 0) {
          std::cout << std::string(60, '-') << "\n";
          std::cout << std::format("Hit i={}, Idx={}, Z={:.4f}, rawPhi={:.4f} | Fitter N={}\n",
                                   i, idx, z, rawPhi, _lineFitter.qn());
        }

        // If we have enough points in fitter, predict from line
        if (_lineFitter.qn() >= 3) {
          double lineSlope     = _lineFitter.dydx();
          double lineIntercept = _lineFitter.y0();

          double predictedPhi  = lineSlope * z + lineIntercept;

          double diff  = predictedPhi - rawPhi;
          int turns    = std::round(diff / (2.0 * M_PI));
          candidatePhi = rawPhi + turns * 2.0 * M_PI;

          if (_debugLevel > 0) {
            std::cout << std::format("  [Mode: PREDICT]\n");
            std::cout << std::format("  Slope = {:.6f}, Intercept = {:.6f}\n",
                                     lineSlope, lineIntercept);
            std::cout << std::format("  Predicted Phi = {:.4f}\n", predictedPhi);
            std::cout << std::format("  Diff = {:.4f} -> Turns = {}\n", diff, turns);
          }
        }
        else {
          // fallback: compare to previous hit
          double prevPhi = phi_unwrapped.back();
          double diff    = prevPhi - rawPhi;
          int turns      = std::round(diff / (2.0 * M_PI));
          candidatePhi   = rawPhi + turns * 2.0 * M_PI;

          if (_debugLevel > 0) {
            std::cout << std::format("  [Mode: FALLBACK (Prev Hit)]\n");
            std::cout << std::format("  Prev Phi = {:.4f}\n", prevPhi);
            std::cout << std::format("  Diff = {:.4f} -> Turns = {}\n", diff, turns);
          }
        }

        if (_debugLevel > 0) {
          std::cout << std::format("  => Chosen Candidate Phi = {:.4f}\n", candidatePhi);
        }

        phi_unwrapped.push_back(candidatePhi);
        _lineFitter.addPoint(z, candidatePhi, weight);
      }

      if (_debugLevel > 0) {
        std::cout << ">>> END OF UNWRAPPING TRACE <<<\n\n";
      }

      // Step E: phi-z fit using _lineFitter
      _lineFitter.clear();
      for (size_t i = 0; i < hitIdx.size(); ++i) {
        size_t idx = hitIdx[i];
        double z   = _tcHits[idx].z;
        double phi = phi_unwrapped[i];
        double w   = 1.0 / _tcHits[idx].helixPhiError2;

        _lineFitter.addPoint(z, phi, w);
      }

      // Step F: retrieve fit results
      double dphidz    = _lineFitter.dydx();
      double dphidzErr = _lineFitter.dydxErr();
      double phi0      = _lineFitter.y0();
      double phi0Err   = _lineFitter.y0Err();
      double chi2ndf   = _lineFitter.chi2Dof();

      std::string stageName = "Step G";

      if (_debugLevel > 0) {
        std::cout << "\n" << std::string(80, '=') << "\n";
        std::cout << std::format(" [findHelix_ver3] {}: Segment Index {}\n",
                                 stageName, segIdx);
        std::cout << std::string(80, '-') << "\n";

        std::cout << std::format("  {:<12} = {:<12.7f} +/- {:.7f}\n",
                                 "dPhi/dZ", dphidz, dphidzErr);
        std::cout << std::format("  {:<12} = {:<12.7f} +/- {:.7f}\n",
                                 "Phi0", phi0, phi0Err);
        std::cout << std::format("  {:<12} = {:<12.4f}\n",
                                 "Chi2/NDF", chi2ndf);
        std::cout << std::format("  {:<12} = {:<12}\n",
                                 "nHits", hitIdx.size());

        std::cout << std::string(80, '-') << "\n";
        std::cout << std::format("{:<8} {:<15} {:<15} {:<15}\n",
                                 "HitIdx", "Z [mm]", "HelixPhi", "PhiUnwrapped");
        std::cout << std::string(80, '-') << "\n";
      }

      for (size_t i = 0; i < hitIdx.size(); ++i) {
        size_t idx = hitIdx[i];

        if (_debugLevel > 0) {
          std::cout << std::format("{:<8} {:<15.4f} {:<15.7f} {:<15.7f}\n",
                                   idx, _tcHits[idx].z,
                                   _tcHits[idx].helixPhi, phi_unwrapped[i]);
        }

        // Critical update: store internally-unwrapped phi back to hit
        _tcHits[idx].helixPhi = phi_unwrapped[i];
      }

      if (_debugLevel > 0) {
        std::cout << std::string(80, '=') << std::endl;
      }
    }

    // Step 3: iteratively align and merge segments by ±2π
    // Step 3: iteratively align and merge segments by ±2π
    while (true) {

      // ------------------------------------------------------------
      // Step 3.1: find current unique segment indices
      // ------------------------------------------------------------
      std::set<int> currentSegments;
      for (size_t i = 0; i < nComboHitsInSegment; ++i) {
        if (!_tcHits[i].used) continue;
        currentSegments.insert(_tcHits[i].segmentIndice);
      }

      if (_debugLevel > 0) {
        std::cout << "\n=============================================\n";
        std::cout << "Current number of segments = " << currentSegments.size() << "\n";
        std::cout << "=============================================\n";
      }

      // No.1: if only 1 segment remains, stop
      if (currentSegments.size() <= 1) break;

      // ------------------------------------------------------------
      // Step 3.2: for each segment, compute fit information
      // ------------------------------------------------------------
      std::vector<int>    curSegIndexList;
      std::vector<int>    curSegNHits;
      std::vector<double> curSegZc;
      std::vector<double> curSegPhic;
      std::vector<double> curSegdPhidZ;
      std::vector<double> curSegPhi0;
      std::vector<double> curSegdPhidZErr;
      std::vector<double> curSegPhi0Err;
      std::vector<double> curSegChi2;
      std::vector<double> curSegPhiErr;

      for (int segIdx : currentSegments) {
        std::vector<size_t> hitIdx;

        for (size_t i = 0; i < nComboHitsInSegment; ++i) {
          if (!_tcHits[i].used) continue;
          if (_tcHits[i].segmentIndice == segIdx) hitIdx.push_back(i);
        }

        if (hitIdx.size() < 2) continue;

        std::sort(hitIdx.begin(), hitIdx.end(),
                  [&](size_t a, size_t b) {
                    return _tcHits[a].z < _tcHits[b].z;
                  });

        // fit current segment using already-unwrapped helixPhi
        _lineFitter.clear();
        double sumW = 0.0, sumZ = 0.0, sumP = 0.0;

        for (size_t k = 0; k < hitIdx.size(); ++k) {
          size_t idx = hitIdx[k];
          double z   = _tcHits[idx].z;
          double phi = _tcHits[idx].helixPhi;
          double w   = 1.0 / _tcHits[idx].helixPhiError2;

          _lineFitter.addPoint(z, phi, w);

          sumW += w;
          sumZ += w * z;
          sumP += w * phi;
        }

        double zCentroid   = sumZ / sumW;
        double phiCentroid = sumP / sumW;

        curSegIndexList.push_back(segIdx);
        curSegNHits.push_back((int)hitIdx.size());
        curSegZc.push_back(zCentroid);
        curSegPhic.push_back(phiCentroid);
        curSegdPhidZ.push_back(_lineFitter.dydx());
        curSegPhi0.push_back(_lineFitter.y0());
        curSegdPhidZErr.push_back(_lineFitter.dydxErr());
        curSegPhi0Err.push_back(_lineFitter.y0Err());
        curSegChi2.push_back(_lineFitter.chi2Dof());
        curSegPhiErr.push_back((sumW > 0.0) ? std::sqrt(1.0 / sumW) : 0.0);

        if (_debugLevel > 0) {
          std::cout << "Segment " << segIdx
                    << " nHits=" << hitIdx.size()
                    << " zC=" << zCentroid
                    << " phiC=" << phiCentroid
                    << " dphidz=" << _lineFitter.dydx()
                    << " phi0=" << _lineFitter.y0()
                    << " chi2=" << _lineFitter.chi2Dof()
                    << "\n";
        }
      }

      if (curSegIndexList.size() <= 1) break;

      // ------------------------------------------------------------
      // Step 3.3: choose reference segment = largest number of hits
      // ------------------------------------------------------------
      int iRef = 0;
      for (int i = 1; i < (int)curSegIndexList.size(); ++i) {
        if (curSegNHits[i] > curSegNHits[iRef]) iRef = i;
      }

      int refSegIdx = curSegIndexList[iRef];

      if (_debugLevel > 0) {
        std::cout << "\nReference segment = " << refSegIdx
                  << " with nHits = " << curSegNHits[iRef] << "\n";
      }

      // ------------------------------------------------------------
      // Step 3.4: find the best segment to merge into the reference
      // ------------------------------------------------------------
      bool foundMerge = false;
      int  bestSegArrayIndex = -1;
      int  bestNturn = 0;
      double bestResidual = 1.0e30;

      for (int i = 0; i < (int)curSegIndexList.size(); ++i) {
        if (i == iRef) continue;

        double zc   = curSegZc[i];
        double phiC = curSegPhic[i];

        // predict phi at this segment centroid using reference segment fit
        double predPhi = curSegdPhidZ[iRef] * zc + curSegPhi0[iRef];

        double deltaPhi = predPhi - phiC;
        int nturn = std::lround(deltaPhi / (2.0 * M_PI));

        double shiftedPhi = phiC + nturn * 2.0 * M_PI;
        double residual = std::abs(predPhi - shiftedPhi);

        if (_debugLevel > 0) {
          std::cout << " Test seg " << curSegIndexList[i]
                    << " : predPhi=" << predPhi
                    << " phiC=" << phiC
                    << " deltaPhi=" << deltaPhi
                    << " nturn=" << nturn
                    << " shiftedPhi=" << shiftedPhi
                    << " residual=" << residual
                    << "\n";
        }

        double mergeThreshold = 15;  // rad, tune if needed
        if (residual < mergeThreshold && residual < bestResidual) {
          bestResidual = residual;
          bestSegArrayIndex = i;
          bestNturn = nturn;
          foundMerge = true;
        }
      }

      // if nothing can be merged, stop
      if (!foundMerge) {
        if (_debugLevel > 0) {
          std::cout << "No more segments can be merged.\n";
        }
        break;
      }

      int mergeSegIdx = curSegIndexList[bestSegArrayIndex];

      if (_debugLevel > 0) {
        std::cout << "\nMerging segment " << mergeSegIdx
                  << " into reference segment " << refSegIdx
                  << " with nturn = " << bestNturn
                  << " residual = " << bestResidual << "\n";
      }

      // ------------------------------------------------------------
      // Step 3.5: shift all hits of that segment and renumber them
      // ------------------------------------------------------------
      for (size_t j = 0; j < nComboHitsInSegment; ++j) {
        if (!_tcHits[j].used) continue;
        if (_tcHits[j].segmentIndice != mergeSegIdx) continue;

        _tcHits[j].helixPhi += bestNturn * 2.0 * M_PI;
        _tcHits[j].segmentIndice = refSegIdx;
        _tcHits[j].nturn    += bestNturn;
      }

      // loop again, because now the reference segment has more hits
      // ------------------------------------------------------------
      // Step 3.7: optional debug plot after each merge
      // ------------------------------------------------------------
      if (_debugLevel > 0) {
        plot_PhiVsZ_alignment_step(tc, isegment,
                                   bestSegArrayIndex,
                                   mergeSegIdx,
                                   curSegdPhidZ[iRef],
                                   curSegPhi0[iRef]);
      }

      // then loop back and recompute everything from the new merged state
    }

      //fit phi-z with corrected helixPhi after 2pi shift
      _lineFitter.clear();

    if (_debugLevel > 0) {
        std::cout << "\n=======================================================\n";
        std::cout << "Final Phi-Z Fit (After 2pi Shifts Applied)\n";
        std::cout << "-------------------------------------------------------\n";
    }

    for(size_t i=0; i<nComboHitsInSegment; i++){
      if(_tcHits[i].used == false) continue;

      double z = _tcHits[i].z;
      double phi = _tcHits[i].helixPhi;
      double phiWeight = 1.0 / (_tcHits[i].helixPhiError2);

      _lineFitter.addPoint(z, phi, phiWeight);

      // Debug printout for each hit
      if (_debugLevel > 0) {
          // Back-calculate the original unshifted phi for the printout
          double originalPhi = phi - (_tcHits[i].nturn * 2.0 * M_PI);

          std::cout << "  HitIdx="    << std::setw(3) << i
                    << " | Seg="      << std::setw(2) << _tcHits[i].segmentIndice
                    << " | Station="  << std::setw(2) << _tcHits[i].station // Update this if your variable is named differently
                    << " | z="        << std::fixed << std::setprecision(3) << std::setw(8) << z
                    << " | Orig Phi=" << std::setw(8) << originalPhi
                    << " | n(turns)=" << std::setw(2) << _tcHits[i].nturn
                    << " | Shifted Phi=" << std::setw(8) << phi
                    << "\n";
      }
    }

    if (_debugLevel > 0) {
        std::cout << "-------------------------------------------------------\n";
        std::cout << "LineFitter Result:"
                  << " dPhidZ = " << _lineFitter.dydx()
                  << " | Phi0 = " << _lineFitter.y0()
                  << " | chi2dZPhi = " << _lineFitter.chi2Dof() << "\n";
        std::cout << "=======================================================\n\n";
    }



    } //end findHelix_ver3

//-----------------------------------------------------------------------------
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
      hhit._hphi = _tcHits[i].helixPhi + _tcHits[i].nturn * 2 * M_PI;
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
        //std::cout<<"Segment #/ComboHits/StrawHits = "<<j<<"/"<<(int)_segmentHits.at(j).size()<<"/"<<hit_count<<std::endl;
      }
      for(int j=0; j<nSegments; j++){
        tcHitsFill(j);

        //maybe we can delete this ver2 as well (no need anymore),_hsColl is not filled, make sure if parametrs are added into _hsColl
        findHelix_ver3(i, j);

        if (_debugLevel > 0 or _diagLevel > 0) {
        // --- Header ---
        std::cout << "\n" << std::string(130, '=') << "\n";
        std::cout << std::format(" [Segment Hits Detail] Loop i: {}, Segment j: {}\n", i, j);
        std::cout << std::string(130, '-') << "\n";
        std::cout << std::format("{:<8} {:<10} {:<10} {:<8} {:<8} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n",
                                 "Seg#", "HitIdx", "SimID", "Used", "nTurn", "X [mm]", "Y [mm]", "Z [mm]", "Phi [rad]", "HelixPhi", "Err2");
        std::cout << std::string(130, '-') << "\n";

        // --- Loop ---
        for (size_t k = 0; k < _tcHits.size(); ++k) {
            const auto& h = _tcHits[k];
            int simID = -999; // Default if MC is missing

            // --- MC Lookup Logic ---
            // We need both the ComboHit collection and the MC collection to exist
            if (_data.chcol && _data.sdmcColl) {
                // 1. Get the ComboHit using your stored index
                size_t chIndex = h.hitIndice;

                if (chIndex < _data.chcol->size()) {
                    // 2. Get the underlying StrawDigi index
                    // (ComboHits can combine multiple straws; .index(0) gets the first one)
                    size_t digiIndex = _data.chcol->at(chIndex).index(0);

                    // 3. Retrieve the SimID from the MC collection
                    if (digiIndex < _data.sdmcColl->size()) {
                        const auto& mcdigi = _data.sdmcColl->at(digiIndex);
                        if (mcdigi.earlyStrawGasStep().isNonnull()) {
                            simID = mcdigi.earlyStrawGasStep()->simParticle()->id().asInt();
                        }
                    }
                }
            }

            // --- Print Row ---
            std::cout << std::format("{:<8} {:<10} {:<10} {:<8} {:<8} {:<12.4f} {:<12.4f} {:<12.4f} {:<12.4f} {:<12.4f} {:<12.4f}\n",
                                     j,
                                     h.hitIndice,
                                     simID,               // <--- The retrieved SimID
                                     (h.used ? "YES" : "NO"),
                                     h.nturn,
                                     h.x,
                                     h.y,
                                     h.z,
                                     h.phi,
                                     h.helixPhi,
                                     h.helixPhiError2);
        }
        std::cout << std::string(130, '=') << std::endl;
    }

        // --- (2) Diagnostic Plotting ---
        plot_PhiVsZ_forSegment_ver3(i, j);
        double xC = _circleFitter.x0();
        double yC = _circleFitter.y0();
        double rC = _circleFitter.radius();
        plot_XVsY(i, j, "step10", xC, yC, rC);

        // --- (3) Pre-selection based on Hit count ---
        if (_debugLevel > 0) std::cout << "Pre-selection\n";
      int snhits = 0;
      for (size_t k = 0; k < _tcHits.size(); k++) {
            if(_tcHits[k].used == false) continue;
            snhits = snhits + _tcHits[k].strawhits;
      }
      if (_debugLevel > 0) {
        std::cout << std::format("  -> Pre-selection Check: snhits = {} (Threshold: 15)\n", snhits);
      }
      if (snhits < 15) {
        if (_debugLevel > 0) std::cout << "  -> Segment REJECTED (Too few hits)\n";
        //continue;
      }

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
          //findHelix(i, j, *_hsColl, _ev_hseed);//maybe we can delete this function(no need anymore), _hsColl is not filled, make sure if parametrs are added into _hsColl
          //findHelix_ver2(i, j, *_hsColl, _ev_hseed);////maybe we can delete this ver2 as well (no need anymore),_hsColl is not filled, make sure if parametrs are added into _hsColl
          // --- Calculate Helix Parameters ---

        _ev_hseed._helix._radius   = rC;
        _ev_hseed._helix._rcent    = sqrt(xC*xC + yC*yC);
        _ev_hseed._helix._fcent    = polyAtan2(yC, xC);

        // --- Z-Phi Line Fit Parameters ---
        _dphidz = _lineFitter.dydx();
        _ev_hseed._helix._fz0 = _lineFitter.y0();

        // Apply 2*Pi wrapping to fz0
        bool wrapped = false;
        if (_ev_hseed._helix._fz0 >  M_PI) { _ev_hseed._helix._fz0 -= 2 * M_PI; wrapped = true; }
        if (_ev_hseed._helix._fz0 < -M_PI) { _ev_hseed._helix._fz0 += 2 * M_PI; wrapped = true; }

        // _ev_hseed._helix._lambda = 1.0 / _lineFitter.dydx(); // Preserved comment
        // _ev_hseed._helix._fz0    = fz0;                     // Preserved comment

        _ev_hseed._helix._lambda    = 1.0 / _dphidz;
        _ev_hseed._helix._chi2dZPhi = _lineFitter.chi2Dof();
        _ev_hseed._helix._chi2dXY   = _circleFitter.chi2DofCircle();
        _ev_hseed._helix._helicity  = (_ev_hseed._helix._lambda > 0) ? Helicity::poshel : Helicity::neghel;
        // _ev_hseed._helix._chi2dZPhi = _lineFitter.chi2Dof(); // Preserved comment

        // --- Debug Printout ---
        // --- Debug Printout ---
        if (_debugLevel > 0) {
            std::cout << "\n" << std::string(60, '=') << "\n";
            std::cout << std::format(" [Helix Fit Results] Segment (i: {}, j: {})\n", i, j);
            std::cout << std::string(60, '-') << "\n";

            std::cout << std::fixed << std::setprecision(4);
            // Circle Fit Info
            std::cout << std::format("  {:<18} = {:>10.4f} [mm]\n", "Radius", rC);
            std::cout << std::format("  {:<18} = {:>10.4f} [mm]\n", "RCent", _ev_hseed._helix._rcent);
            std::cout << std::format("  {:<18} = ({:.2f}, {:.2f}) [mm]\n", "Center (X,Y)", xC, yC);

            // Line Fit Info
            std::cout << std::format("  {:<18} = {:>10.4f} [rad/mm]\n", "dPhi/dZ", _dphidz);
            std::cout << std::format("  {:<18} = {:>10.4f} [mm/rad]\n", "Lambda", _ev_hseed._helix._lambda);
            std::cout << std::format("  {:<18} = {:>10.4f} [rad] {}\n", "fz0", _ev_hseed._helix._fz0, (wrapped ? "(2Pi Wrapped)" : ""));

            // Quality Metrics
            std::cout << std::string(60, '-') << "\n";
            std::cout << std::format("  {:<18} = {:>10.4f}\n", "Chi2/DOF (XY)", _ev_hseed._helix._chi2dXY);
            std::cout << std::format("  {:<18} = {:>10.4f}\n", "Chi2/DOF (ZPhi)", _ev_hseed._helix._chi2dZPhi);
            std::cout << std::format("  {:<18} = {:>10}\n", "Helicity", (_ev_hseed._helix._helicity == Helicity::poshel ? "Positive" : "Negative"));
            std::cout << std::string(60, '=') << std::endl;

            // --- Future Use / Commented out variables ---
            /*
            std::cout << std::format("  {:<15} = {:>10.4f} [rad]\n", "fz0", _ev_hseed._helix._fz0);
            std::cout << std::format("  {:<15} = {:>10}\n", "Helicity", (_ev_hseed._helix._helicity == Helicity::poshel ? "poshel" : "neghel"));
            std::cout << std::format("  {:<15} = {:>10.4f}\n", "Chi2dXY", _ev_hseed._helix._chi2dXY);
            */
        }

        // --- Diagnostic Plots ---
        if (_diagLevel > 0) {
            // plot_PhiVsZ_forSegment(i, j);
            // plot_CirclePhiVsZ_forSegment(i, j);
            // plot_2PiAmbiguityPhiVsZ_forSegment(i, j);
            // plot_RVsZ_forSegment(i, j);

            // Active Plots
            // plot_2PiAmbiguityPhiVsZ_forSegment_mod(i, j);
            plot_HelixPhiVsZ(i, j);
        }

        std::cout<<"CROSSCHECK222"<<std::endl;
            // Print Header using std::format
            std::cout << std::format("{:<8} {:<10} {:<8} {:<8} {:<12} {:<12} {:<12} {:<12} {:<12} {:<12}\n",
                                     "Seg#", "HitID", "Used", "nTurn", "X [mm]", "Y [mm]", "Z [mm]", "Phi [rad]", "HelixPhi", "Err2");
            std::cout << std::string(115, '-') << "\n";

            // Print Rows using std::format
            _lineFitter.clear();
            for (const auto& h : _tcHits) {
                std::cout << std::format("{:<8} {:<10} {:<8} {:<8} {:<12.4f} {:<12.4f} {:<12.4f} {:<12.4f} {:<12.4f} {:<12.4f}\n",
                                         j,
                                         h.segmentIndice,
                                         (h.used ? "YES" : "NO"),
                                         h.nturn,
                                         h.x,
                                         h.y,
                                         h.z,
                                         h.phi,
                                         h.helixPhi,
                                         h.helixPhiError2);
      double phiWeight = 1.0 / (h.helixPhiError2);
      _lineFitter.addPoint(h.z, h.helixPhi, phiWeight);

            }
            std::cout << std::string(115, '=') << std::endl;
      std::cout << "LineFitter "
      << " dPhidZ =" << _lineFitter.dydx()
      << " Phi0 =" << _lineFitter.y0()
      << " chi2dZPhi =" << _lineFitter.chi2Dof() << std::endl;

        //saveHelix(i, _ev_hseed);//_hsColl is not filled, make sure if parametrs are added into _hsColl
    //HelixSeed hseed;
    const TimeCluster* tculster = &_data.tccol->at(i);
    _ev_hseed._t0 = tculster->_t0;
    auto const& tccH    = _event->getValidHandle<mu2e::TimeClusterCollection>(_tcCollTag);
    _ev_hseed._timeCluster = art::Ptr<TimeCluster>(tccH, i);
    _ev_hseed._hhits.setParent(_data.chcol->parent());
    ::LsqSums2 fitter;
    for (size_t m = 0; m < _tcHits.size(); m++) {
      if(_tcHits[m].used == false) continue;
      int hitIndice = _tcHits[m].hitIndice;
      const ComboHit* hit = &_data.chcol->at(hitIndice);
      fitter.addPoint(hit->pos().z(), hit->correctedTime(), 1 / (hit->timeRes() * hit->timeRes()));
      ComboHit hhit(*hit);
      //hhit._hphi = 33.33;
      hhit._hphi = _tcHits[m].helixPhi;
      _ev_hseed._hhits.push_back(hhit);
    }
    //float eDepAvg = _ev_hseed._hhits.eDepAvg();
    _ev_hseed._t0 = TrkT0(fitter.y0(), fitter.y0Err());
    _ev_hseed._status.merge(TrkFitFlag::helixOK);
    _ev_hseed._status.merge(TrkFitFlag::APRHelix);

          std::cout<<"CROSSCHECK = "<<std::endl;
          size_t nComboHitsInSegment = _tcHits.size();
          int nstrawhits__ = 0;
          for (size_t k = 0; k < nComboHitsInSegment; k++) {
            if(_tcHits[k].used == false) continue;
            //std::cout<<"_tcHits["<<k<<"].hitIndice = "<<_tcHits[k].hitIndice<<std::endl;
            std::cout<<"_tcHits["<<k<<"].hitIndice = "<<_tcHits[k].hitIndice<<std::endl;
            nstrawhits__ = nstrawhits__ + _tcHits[k].strawhits;
          }
          //std::cout<<"nstrawhits = "<<nstrawhits__<<std::endl;
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
      std::cout<<"p  = "<<p<<std::endl;
      std::cout<<"pT  = "<<pT<<std::endl;
      //if(p < 50) continue;
      //if(pT < 50) continue;
      if (_ev_hseed._helix._chi2dXY > 5) continue;
      if (_ev_hseed._helix._chi2dZPhi > 5) continue;
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
   // _ev_hseed._helix._chi2dZPhi = 100;
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

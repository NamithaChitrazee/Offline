#include "TH1.h"
#include "TH2.h"

#include <string.h>

#include "Offline/CalPatRec/inc/DeltaFinder_types.hh"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Utilities/ToolMacros.h"
#include "art/Utilities/make_tool.h"
#include "art/Framework/Principal/Event.h"

#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "Offline/CalPatRec/inc/HlPrint.hh"
#include "Offline/CalPatRec/inc/PhiZSeedFinderAlg.hh"
#include "Offline/CalPatRec/inc/McPart_t.hh"
#include "Offline/CalPatRec/inc/PhiZSeedFinder_types.hh"

#include "Offline/Mu2eUtilities/inc/ModuleHistToolBase.hh"
#include "Offline/Mu2eUtilities/inc/McUtilsToolBase.hh"
#include "Offline/DataProducts/inc/PDGCode.hh"

#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"

using namespace std;

namespace mu2e {

  using namespace PhiZSeedFinderTypes;

  class SimParticle;
  //class StrawDigiMCCollection;

  class PhiZSeedFinderDiag: public ModuleHistToolBase {

    enum {
      kNEventHistSets  =  10,
      kNSegmentHistSets  =  10,
      kNMCInfoHistSets  =  10,
    };

    struct EventHist_t {
      TH1F*  fEventNumber;
      TH1F*  fRunNumber;
      TH1F*  fDiffRadius;
      TH1F*  fFractionHits;
    };

    struct SegmentHist_t {
      TH1F*  h_Chi2dPhdZ;
    };

    struct MCInfoHist_t {
      TH1F*  h_MCnStrawHitsPerParticle;
      TH1F*  h_MCnComboHitsPerParticle;
      TH1F*  h_MCMom;
      TH1F*  h_MCTransverseMom;
      TH1F*  h_MCTanLambda;
    };

                                        // hits referred to here are the combo hits
    struct Hist_t {
      EventHist_t*    fEvent [kNEventHistSets ];
      SegmentHist_t*  _segmentHist [kNSegmentHistSets ];
      MCInfoHist_t*   _mcInfoHist [kNMCInfoHistSets ];
    };

  protected:

    bool           _mcDiag;

    int            _eventNumber;

    TObjArray      _listOfMcParticles; // list of particles with at least one ComboHit in the tracker

    Data_t*        _data;                 // shared data, passed from the caller
    Hist_t         _hist;

    std::unique_ptr<McUtilsToolBase>      _mcUtils;

    float          _ppi;                // proton pulse intensity

  public:

    PhiZSeedFinderDiag(const fhicl::Table<mu2e::PhiZSeedFinderTypes::Config>& config);
    ~PhiZSeedFinderDiag();
//-----------------------------------------------------------------------------
// accessors
//-----------------------------------------------------------------------------
    McPart_t*   mcPart   (int I) { return (McPart_t*) _listOfMcParticles.At(I); }
    //mu2e::PhiZSeedMcPart_t* mcPart(int I) {
    //  return static_cast<mu2e::PhiZSeedMcPart_t*>(_listOfMcParticles.At(I));
    //}
//-----------------------------------------------------------------------------
// other functions
//-----------------------------------------------------------------------------
    void        bookEventHistograms (EventHist_t*  Hist, art::TFileDirectory* Dir);
    void        bookSegmentHistograms (SegmentHist_t*  Hist, art::TFileDirectory* Dir);
    void        bookMCInfoHistograms (MCInfoHist_t*  Hist, art::TFileDirectory* Dir);

    void        fillEventHistograms (EventHist_t*  Hist);
    void        fillSegmentHistograms (SegmentHist_t*  Hist);
    void        fillMCInfoHistograms (MCInfoHist_t*  Hist);
//-----------------------------------------------------------------------------
// overriden virtual functions of the base class
//-----------------------------------------------------------------------------
  public:

    virtual int bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) override ;
    virtual int fillHistograms(void* Data, int Mode = -1 ) override ;
    virtual int debug         (void* Data, int Mode = -1 ) override ;
  };

//-----------------------------------------------------------------------------
// this routine is called once per job
//-----------------------------------------------------------------------------
  PhiZSeedFinderDiag::PhiZSeedFinderDiag(const fhicl::Table<mu2e::PhiZSeedFinderTypes::Config>& config):
    _mcDiag                (config().mcDiag()                )
    // _printOTracker         (config().printOTracker()         ),
    // _printComboHits        (config().printComboHits()        ),
    // _printGoodComboHits    (config().printGoodComboHits()    ),
    // _printShcol            (config().printShcol()            ),
  {
    printf(" PhiZSeedFinderDiag::PhiZSeedFinderDiag : HOORAY Config! \n");

    if (_mcDiag != 0) _mcUtils = art::make_tool  <McUtilsToolBase>(config().mcUtils,"mcUtils");
    else              _mcUtils = std::make_unique<McUtilsToolBase>();
  }

//-----------------------------------------------------------------------------
  PhiZSeedFinderDiag::~PhiZSeedFinderDiag() {
  }

//-----------------------------------------------------------------------------
  void PhiZSeedFinderDiag::bookEventHistograms(EventHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->fEventNumber     = Dir->make<TH1F>("event" , "Event Number", 100, 0., 100000.);
    Hist->fDiffRadius      = Dir->make<TH1F>("diff radius" , "MC Radius - Fit Radius", 400, -200.0, 200.0);
    Hist->fFractionHits    = Dir->make<TH1F>("fraction hits" , "Fraction hits", 20, 0, 1.0);
  }

//-----------------------------------------------------------------------------
  void PhiZSeedFinderDiag::bookSegmentHistograms(SegmentHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->h_Chi2dPhdZ     = Dir->make<TH1F>("Chi2dPhdZ" , "Chi2PhiZ per segment", 500, 0., 50.);
  }
//-----------------------------------------------------------------------------
  void PhiZSeedFinderDiag::bookMCInfoHistograms(MCInfoHist_t* Hist, art::TFileDirectory* Dir) {
    Hist->h_MCnStrawHitsPerParticle   = Dir->make<TH1F>("MCnStrawHitsPerParticle" , "number of StrawHits per MC particle", 100, 0., 100000.);
    Hist->h_MCnComboHitsPerParticle   = Dir->make<TH1F>("MCnComboHitsPerParticle" , "number of ComboHits per MC particle", 400, -200.0, 200.0);
    Hist->h_MCMom                     = Dir->make<TH1F>("MCMom" , "momentum", 400, -200.0, 200.0);
    Hist->h_MCTransverseMom           = Dir->make<TH1F>("MCTransverseMom" , "transverse momentum of MC particle", 400, -200.0, 200.0);
    Hist->h_MCTanLambda               = Dir->make<TH1F>("MCTanLambda" , "tan(lambda) of MC particle", 400, -200.0, 200.0);
  }

//-----------------------------------------------------------------------------
// this routine is called once per job (likely, from beginJob)
// TH1::AddDirectory makes sure one can have histograms with the same name
// in different subdirectories
//-----------------------------------------------------------------------------
  int PhiZSeedFinderDiag::bookHistograms(art::ServiceHandle<art::TFileService>& Tfs) {

    TH1::AddDirectory(0);
    char folder_name[20];
//-----------------------------------------------------------------------------
// book event-level histograms - fill them once per event
//-----------------------------------------------------------------------------
    int book_event_histset[kNEventHistSets];
    for (int i=0; i<kNEventHistSets; i++) book_event_histset[i] = 0;

    book_event_histset[ 0] = 1;                // all events

    for (int i=0; i<kNEventHistSets; i++) {
      if (book_event_histset[i] != 0) {
        sprintf(folder_name,"evt_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist.fEvent[i] = new EventHist_t;
        bookEventHistograms(_hist.fEvent[i],&dir);
      }
    }

//-----------------------------------------------------------------------------
// book Segment histograms - fill them once per event
//-----------------------------------------------------------------------------
    int book_segment_histset[kNSegmentHistSets];
    for (int i=0; i<kNSegmentHistSets; i++) book_segment_histset[i] = 0;

    book_segment_histset[ 0] = 1;                // all events

    for (int i=0; i<kNSegmentHistSets; i++) {
      if (book_segment_histset[i] != 0) {
        sprintf(folder_name,"evt_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist._segmentHist[i] = new SegmentHist_t;
        bookSegmentHistograms(_hist._segmentHist[i],&dir);
      }
    }
//-----------------------------------------------------------------------------
// book MC info histograms - fill them once per event
//-----------------------------------------------------------------------------
    int book_mcinfo_histset[kNMCInfoHistSets];
    for (int i=0; i<kNMCInfoHistSets; i++) book_mcinfo_histset[i] = 0;

    book_mcinfo_histset[ 0] = 1;                // all events

    for (int i=0; i<kNMCInfoHistSets; i++) {
      if (book_mcinfo_histset[i] != 0) {
        sprintf(folder_name,"evt_%i",i);
        art::TFileDirectory dir = Tfs->mkdir(folder_name);

        _hist._mcInfoHist[i] = new MCInfoHist_t;
        bookMCInfoHistograms(_hist._mcInfoHist[i],&dir);
      }
    }

//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    return 0;
  }

//-----------------------------------------------------------------------------
  void  PhiZSeedFinderDiag::fillEventHistograms(EventHist_t* Hist) {

    int event_number = _data->event->event();
    Hist->fEventNumber->Fill(event_number);

    for(int j=0; j<(int)_data->h_diffradius.size(); j++){
      double diffradius = _data->h_diffradius[j];
      Hist->fDiffRadius->Fill(diffradius);
    }
    //h_fraction
    for(int j=0; j<(int)_data->h_fraction.size(); j++){
      double fraction = _data->h_fraction[j];
      Hist->fFractionHits->Fill(fraction);
    }

  }

//-----------------------------------------------------------------------------
  void  PhiZSeedFinderDiag::fillSegmentHistograms(SegmentHist_t* Hist) {
  std::cout<<"fillSegmentHistograms "<<std::endl;
    for(int i=0; i<(int)_data->h_lineFitter_chi2Dof[0].size(); i++){
      Hist->h_Chi2dPhdZ->Fill(_data->h_lineFitter_chi2Dof[0].at(i));
      std::cout<<"_data->h_lineFitter_chi2Dof[0].at("<<i<<") = "<<_data->h_lineFitter_chi2Dof[0].at(i)<<std::endl;
    }

  }

//-----------------------------------------------------------------------------
  void  PhiZSeedFinderDiag::fillMCInfoHistograms(MCInfoHist_t* Hist) {

    for(int i=0; i<(int)_data->h_MCnStrawHitsPerParticle.size(); i++){
      Hist->h_MCnStrawHitsPerParticle->Fill(10);
    }
    for(int i=0; i<(int)_data->h_MCnComboHitsPerParticle.size(); i++){
      Hist->h_MCnComboHitsPerParticle->Fill(10);
    }
    for(int i=0; i<(int)_data->h_MCMom.size(); i++){
      Hist->h_MCMom->Fill(10);
    }
    for(int i=0; i<(int)_data->h_MCTransverseMom.size(); i++){
      Hist->h_MCTransverseMom->Fill(10);
    }
    for(int i=0; i<(int)_data->h_MCTanLambda.size(); i++){
      Hist->h_MCTanLambda->Fill(10);
    }


  }
//-----------------------------------------------------------------------------
// main fill histograms function called once per event
// 'Mode' not used
//-----------------------------------------------------------------------------
  int PhiZSeedFinderDiag::fillHistograms(void* Data, int Mode) {
//-----------------------------------------------------------------------------
// event histograms - just one set
//-----------------------------------------------------------------------------
    _data = (Data_t*) Data;
    fillEventHistograms(_hist.fEvent[0]);
    fillSegmentHistograms(_hist._segmentHist[0]);
    fillMCInfoHistograms(_hist._mcInfoHist[0]);
//-----------------------------------------------------------------------------
// done
//-----------------------------------------------------------------------------
    return 0;
  }

//-----------------------------------------------------------------------------
// debugLevel > 0: print seeds
// Mode = 1:
// Mode = 2: event
//-----------------------------------------------------------------------------
  int PhiZSeedFinderDiag::debug(void* Data, int Mode) {
    return 0;
  }

}

DEFINE_ART_CLASS_TOOL(mu2e::PhiZSeedFinderDiag)

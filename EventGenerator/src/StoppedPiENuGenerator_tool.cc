///////////////////////////////////////////////////////////////////////////////
// P.Murat: use 'ProcessCode::mu2ePienu' for pi--> e nu decay
//          also use GenId::piEplusNuGun
///////////////////////////////////////////////////////////////////////////////
#include "art/Utilities/ToolMacros.h"
#include <memory>

#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandGeneral.h"

#include "Offline/EventGenerator/inc/ParticleGeneratorTool.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/MCDataProducts/inc/GenId.hh"
#include "Offline/Mu2eUtilities/inc/RandomUnitSphere.hh"
#include "Offline/Mu2eUtilities/inc/BinnedSpectrum.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"

#include "fhiclcpp/types/DelegatedParameter.h"

namespace mu2e {
  class StoppedPiEnuGenerator : public ParticleGeneratorTool {

  private:
    PDGCode::type                       _pdgCode;
    ProcessCode                         _processCode;
    double                              _mass;
    BinnedSpectrum                      _spectrum;
    std::unique_ptr<RandomUnitSphere>   _randomUnitSphere;
    std::unique_ptr<CLHEP::RandGeneral> _randSpectrum;
    double                              _czmin;
    double                              _czmax;

    enum SpectrumVar  { TOTAL_ENERGY, KINETIC_ENERGY, MOMENTUM };
    SpectrumVar                         _var;

  public:
    struct PhysConfig {
      using Name   = fhicl::Name;
      using Comment= fhicl::Comment;

      fhicl::DelegatedParameter spectrum   {Name("spectrum"   ), Comment("BinnedSpectrum parameters)")   };
      fhicl::Atom<std::string>  processCode{Name("processCode"), Comment("Mu2e process code"         )   };
      fhicl::Atom<double>       czmin      {Name("czmin"      ), Comment("cos(theta) min"            ),-1};
      fhicl::Atom<double>       czmax      {Name("czmax"      ), Comment("cos(theta) max"            ), 1};
    };
    typedef art::ToolConfigTable<PhysConfig> Parameters;

    explicit StoppedPiEnuGenerator(Parameters const& conf) :
      _spectrum(BinnedSpectrum(conf().spectrum.get<fhicl::ParameterSet>())),
      _czmin       (conf().czmin       ()),
      _czmax       (conf().czmax       ())
   {
      _pdgCode     = PDGCode::type(conf().spectrum.get<fhicl::ParameterSet>().get<int>("pdgCode"));
      _processCode = ProcessCode::findByName(conf().processCode());
      _mass        = GlobalConstantsHandle<ParticleDataList>()->particle(_pdgCode).mass();

      std::string v = conf().spectrum.get<fhicl::ParameterSet>().get<std::string>("spectrumVariable");
      if      (v == "totalEnergy"  )  _var = TOTAL_ENERGY;
      else if (v == "kineticEnergy")  _var = KINETIC_ENERGY;
      else if (v == "momentum"     )  _var = MOMENTUM;
      else {
        throw cet::exception("BADCONFIG") << "StoppedPiEnuGenerator::generate: unknown var name "<< v << "\n";
      }
   }

    virtual std::vector<ParticleGeneratorTool::Kinematic> generate() override;

    virtual void generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) override;

    virtual ProcessCode   processCode() override { return _processCode; }

    virtual void finishInitialization(art::RandomNumberGenerator::base_engine_t& eng, const std::string&) override {
      _randomUnitSphere = std::make_unique<RandomUnitSphere>(eng);
      _randSpectrum     = std::make_unique<CLHEP::RandGeneral>(eng, _spectrum.getPDF(), _spectrum.getNbins());
    }
  };

//-----------------------------------------------------------------------------
  std::vector<ParticleGeneratorTool::Kinematic> StoppedPiEnuGenerator::generate() {
    std::vector<ParticleGeneratorTool::Kinematic>  res;

    double e(0), p(0);

    if (_var == TOTAL_ENERGY) {
      e = _spectrum.sample(_randSpectrum->fire());
      p = sqrt(e*e -_mass*_mass);
    }
    else if (_var == MOMENTUM) {
      p = _spectrum.sample(_randSpectrum->fire());
      e = sqrt(p*p +_mass*_mass);
    }
    else if (_var == KINETIC_ENERGY) {
      e = _mass + _spectrum.sample(_randSpectrum->fire());
      p = sqrt(e*e -_mass*_mass);
    }
//-----------------------------------------------------------------------------
// for now, just check the cos(theta)
//-----------------------------------------------------------------------------
    while (1) {
      CLHEP::HepLorentzVector fourmom(_randomUnitSphere->fire(p),e);
      double cz = fourmom.pz()/p;
      if ((cz >= _czmin) and (cz <= _czmax)) {
        Kinematic k{_pdgCode, ProcessCode::mu2ePienu, fourmom};
        res.emplace_back(k);
        break;
      }
    }

    return res;
  }

  void StoppedPiEnuGenerator::generate(std::unique_ptr<GenParticleCollection>& out, const IO::StoppedParticleF& stop) {
    const CLHEP::Hep3Vector pos(stop.x, stop.y, stop.z);
    const auto daughters = generate();
    for(const auto& d: daughters) {
      out->emplace_back(d.pdgId, GenId::piEplusNuGun, pos, d.fourmom, stop.t);
    }
  }

}
DEFINE_ART_CLASS_TOOL(mu2e::StoppedPiEnuGenerator)
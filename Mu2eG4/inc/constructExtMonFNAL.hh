// Andrei Gaponenko, 2012

#ifndef CONSTRUCTEXTMONFNAL_HH
#define CONSTRUCTEXTMONFNAL_HH

#include <string>

#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALPlaneStack.hh"
#include "Geant4/G4ThreeVector.hh"

namespace CLHEP { class HepRotation; }

namespace mu2e {

  class SimpleConfig;
  class ExtMonFNALMagnet;

  void constructExtMonFNAL(const VolumeInfo& collimator1Parent,
                           const CLHEP::HepRotation& collimator1ParentRotationInMu2e,
                           const VolumeInfo& mainParent,
                           const CLHEP::HepRotation& mainParentRotationInMu2e,
                           const SimpleConfig& config);

  void constructExtMonFNALBuilding(const VolumeInfo& collimator1Parent,
                                   const CLHEP::HepRotation& collimator1ParentRotationInMu2e,
                                   const VolumeInfo& mainParent,
                                   const CLHEP::HepRotation& mainParentRotationInMu2e,
                                   const SimpleConfig& config);

  void constructExtMonFNALDetector(const VolumeInfo& room, const SimpleConfig& config);

  void constructExtMonFNALDetector(const VolumeInfo& mainParent,
                                   const CLHEP::HepRotation& mainParentRotationInMu2e,
                                   const SimpleConfig& config);

  void constructExtMonFNALScintillators(const VolumeInfo& mother,
                                        const ExtMonFNALPlaneStack& stack,
                                        const std::string& volNameSuffix,
                                        const SimpleConfig& config,
                                        bool const forceAuxEdgeVisible,
                                        bool const doSurfaceCheck,
                                        bool const placePV
                                        );

}

#endif /* CONSTRUCTEXTMONFNAL_HH */

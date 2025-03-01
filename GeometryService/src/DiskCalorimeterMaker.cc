//
// Make a Disk Calorimeter.
//
// original authors  Bertrand Echenarrd
//
// Disk geometry
//
//  see Mu2eG4/src/constructDiskCalorimeter.cc for the geometry
//    the disk are split in three pieces: front palte, disk case and back plate
//
//  front plate has calibration pipes
//  middle case has crystal enveloped in wrapper with front cap (back of the crystal is free)
//  back plate has a solid plate with holes and readouts inside + FEE box outside
//
//    Coordinate system
//
//           - crystal:
//                      origin:      base of the crystal
//                      orientation: crystal oriented along the z direction (direction front face - back face)
//                                   the front (back) face is at z=0 (z=crystal length).
//           - disk:
//                      origin:      center of the disk
//                      orientation: along the z axis
//                      other:       extra coordinate system placed at the front of the crystal - used for track-calo matching
//
//           - calorimeter:
//                      origin:      center of mother volume
//
//         Note: The extra disk coordinate system is suffixed "FF" for FrontFace. The FF z origin is at the same position of the crystal z origin.
//
// For reference, git tag (ef94504f51edbbfeb54a5e63651856bdf5c0a60d) has generic placement of the disk origin.
//

#include "cetlib_except/exception.h"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/Disk.hh"
#include "Offline/CalorimeterGeom/inc/Crystal.hh"
#include "Offline/GeometryService/inc/DiskCalorimeterMaker.hh"

#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"




namespace mu2e {

    DiskCalorimeterMaker::DiskCalorimeterMaker(SimpleConfig const& config, double solenoidOffset)
    {

       calo_ = std::unique_ptr<DiskCalorimeter>(new DiskCalorimeter());

       verbosityLevel_ = config.getInt("calorimeter.verbosityLevel",0);

       calo_->caloInfo_.set("caloDiskRadiusIn",       config.getDouble("calorimeter.caloDiskRadiusIn") );    //inner radius of disk enveloppe
       calo_->caloInfo_.set("caloDiskRadiusOut",      config.getDouble("calorimeter.caloDiskRadiusOut") );   //outer radius of disk enveloppe
       calo_->caloInfo_.set("caloFEBRadiusOut",       config.getDouble("calorimeter.caloFEBRadiusOut") );    //outer radius of FEB enveloppe
       calo_->caloInfo_.set("caloMotherZ0",           config.getDouble("calorimeter.caloMotherZ0") );        //upstream Z of calorimeter enveloppe volume
       calo_->caloInfo_.set("caloMotherZ1",           config.getDouble("calorimeter.caloMotherZ1") );        //downstream Z of calorimeter enveloppe volume
       calo_->caloInfo_.set("vdThickness",            config.getDouble("calorimeter.vdThickness") );         //virtual detector thickness

       calo_->caloInfo_.set("diskInAlRingRIn",        config.getDouble("calorimeter.diskInAlRingRIn") );     //inner radius of Al ring inside disk
       calo_->caloInfo_.set("diskInAlRingZLength",    config.getDouble("calorimeter.diskInAlRingZLength") ); //Z length of Al ring inside disk
       calo_->caloInfo_.set("diskInCFRingRIn",        config.getDouble("calorimeter.diskInCFRingRIn") );     //inner radius of carbon fiber ring
       calo_->caloInfo_.set("diskInCFRingROut",       config.getDouble("calorimeter.diskInCFRingROut") );    //outer radius of carbon fiber ring
       calo_->caloInfo_.set("diskCaseZLength",        config.getDouble("calorimeter.crystalZLength")         //Z Length of crystal support structure  (large Al ring)
                                                    + config.getDouble("calorimeter.crystalCapZLength") );
       calo_->caloInfo_.set("diskCrystalRIn",         config.getDouble("calorimeter.diskCrystalRIn") );      //inner radius of volume containing crystals
       calo_->caloInfo_.set("diskCrystalROut",        config.getDouble("calorimeter.diskCrystalROut") );     //outer radius of volume containing crystals
       calo_->caloInfo_.set("diskCaseRingROut",       config.getDouble("calorimeter.diskCaseRingROut") );    //outer radius of crystal support structure (large Al ring)
       calo_->caloInfo_.set("diskOutRailZLength",     config.getDouble("calorimeter.diskOutRailZLength") );  //Z length of rail outside the crystal support structure
       calo_->caloInfo_.set("diskOutRailROut",        config.getDouble("calorimeter.diskOutRailROut") );     //outer radius of rail outside the crystal support structure
       calo_->caloInfo_.set("diskStepThickness",      config.getDouble("calorimeter.diskStepThickness") );   //thickness of carbon fiber steps glued to inner hole
       calo_->caloInfo_.set("diskZMotherShift",       getVDouble(config,"diskZMotherShift",CaloConst::_nDisk)); //Z offset betweenfront face of disk+FEB envelope and calo mother volume Z0

       calo_->caloInfo_.set("crystalXYLength",        config.getDouble("calorimeter.crystalXYLength") );     //crysal transverse dimension
       calo_->caloInfo_.set("crystalZLength",         config.getDouble("calorimeter.crystalZLength") );      //crystal Z length
       calo_->caloInfo_.set("crystalCapZLength",      config.getDouble("calorimeter.crystalCapZLength") );   //crystal cap Z length (only in front of crystal)
       calo_->caloInfo_.set("wrapperThickness",       config.getDouble("calorimeter.wrapperThickness") );    //crystal wrapper thickness
       calo_->caloInfo_.set("refractiveIndex",        config.getDouble("calorimeter.refractiveIndex") );     //crystal refractive index

       calo_->caloInfo_.set("nSiPMPerCrystal",        CaloConst::_nSiPMPerCrystal );                         //number of SiPM per crystal
       calo_->caloInfo_.set("readoutXLength",         config.getDouble("calorimeter.readoutXLength") );      //SiPM x length
       calo_->caloInfo_.set("readoutYLength",         config.getDouble("calorimeter.readoutYLength") );      //SiPM Y length
       calo_->caloInfo_.set("readoutZLength",         config.getDouble("calorimeter.readoutZLength") );      //SiPM Z length
       calo_->caloInfo_.set("shimStepsInRowId",       getVInt(config,"shimStepsInRowId") );                  //indices of crystal rows w.r.t center calorimeter with inner shims
       calo_->caloInfo_.set("shimStepsOutRowId",      getVInt(config,"shimStepsOutRowId") );                 //indices of crystal rows w.r.t center calorimeter with outer shims
       std::vector<int> tempInt;
       for (uint16_t i=0;i<CaloConst::_nCaphriCrystal;++i) tempInt.push_back(CaloConst::_caphriId[i]);
       calo_->caloInfo_.set("caphriCrystalId",tempInt);                                                      //LYSO crystal ids

       calo_->caloInfo_.set("hasCrates",              config.getBool  ("calorimeter.hasCrates") );           //include or exclude crates
       calo_->caloInfo_.set("hasFrontPanel",          config.getBool  ("calorimeter.hasFrontPanel") );       //include or exclude front panel with calibration pipes
       calo_->caloInfo_.set("hasBackPanel",           config.getBool  ("calorimeter.hasBackPanel") );        //include or exclude back panel with readouts and cooling pipes
       calo_->caloInfo_.set("FEEXLength",             config.getDouble("calorimeter.FEEXLength") );          //X length of FEE copper box
       calo_->caloInfo_.set("FEEYLength",             config.getDouble("calorimeter.FEEYLength") );          //Y length of FEE copper box
       calo_->caloInfo_.set("FEEZLength",             config.getDouble("calorimeter.FEEZLength") );          //Z length of FEE copper box
       calo_->caloInfo_.set("FEEBoxThickness",        config.getDouble("calorimeter.FEEBoxThickness") );     //FEE copper box thickness
       calo_->caloInfo_.set("BPStripThickness",       config.getDouble("calorimeter.BPStripThickness") );    //Back plate cooling strip thickness
       calo_->caloInfo_.set("BPHoleXLength",          config.getDouble("calorimeter.BPHoleXLength") );       //X length of hole in back plate
       calo_->caloInfo_.set("BPHoleYLength",          config.getDouble("calorimeter.BPHoleYLength") );       //Y length of hole in back plate
       calo_->caloInfo_.set("BPHoleZLength",          config.getDouble("calorimeter.BPHoleZLength") );       //Z length of hole in back plate
       calo_->caloInfo_.set("BPOuterRadius",          config.getDouble("calorimeter.BPOuterRadius") );       //Back plate outer radius
       calo_->caloInfo_.set("BPPipeRadiusHigh",       config.getDouble("calorimeter.BPPipeRadiusHigh") );    //Radius of large back plate outer pipe
       calo_->caloInfo_.set("BPPipeRadiusLow",        config.getDouble("calorimeter.BPPipeRadiusLow") );     //Radius of small back plate outer pipe
       calo_->caloInfo_.set("BPPipeThickness",        config.getDouble("calorimeter.BPPipeThickness") );     //Thickness of back plate outer pipe
       calo_->caloInfo_.set("BPPipeZOffset",          config.getDouble("calorimeter.BPPipeZOffset") );       //Z offset between outer pipe center and back plate

       calo_->caloInfo_.set("FPInnerRadius",          config.getDouble("calorimeter.FPInnerRadius") );       //Front plate inner radius
       calo_->caloInfo_.set("FPOuterRadius",          config.getDouble("calorimeter.FPOuterRadius") );       //Front plate outer radius
       calo_->caloInfo_.set("FPFoamZLength",          config.getDouble("calorimeter.FPFoamZLength") );       //Thickness of foam inside front plate
       calo_->caloInfo_.set("FPCarbonZLength",        config.getDouble("calorimeter.FPCarbonZLength") );     //Z length of carbon fiber panel on front plate
       calo_->caloInfo_.set("FPCoolPipeTorRadius",    config.getDouble("calorimeter.FPCoolPipeTorRadius") ); //Tor radius of front plate outer cooling pipe
       calo_->caloInfo_.set("FPCoolPipeRadius",       config.getDouble("calorimeter.FPCoolPipeRadius") );    //Inner radius of front plate outer cooling pipe
       calo_->caloInfo_.set("FPCoolPipeThickness",    config.getDouble("calorimeter.FPCoolPipeThickness") ); //Thickness of front plate outer cooling pipe
       calo_->caloInfo_.set("nPipes",                 config.getInt   ("calorimeter.nPipes") );              //Number of cooling pipes inside front plate
       calo_->caloInfo_.set("pipeRadius",             config.getDouble("calorimeter.pipeRadius") );          //Inner radius of cooling pipes
       calo_->caloInfo_.set("pipeThickness",          config.getDouble("calorimeter.pipeThickness") );       //Thickness of cooling pipes
       calo_->caloInfo_.set("pipeInitSeparation",     config.getDouble("calorimeter.pipeInitSeparation") );  //Cooling pipe characteristics (see constructDiskCalorimeter)
       calo_->caloInfo_.set("radSmTor",               config.getDouble("calorimeter.radSmTor") );            //Cooling pipe characteristics (see constructDiskCalorimeter)
       calo_->caloInfo_.set("xsmall",                 config.getDouble("calorimeter.xsmall") );              //Cooling pipe characteristics (see constructDiskCalorimeter)
       calo_->caloInfo_.set("xdistance",              config.getDouble("calorimeter.xdistance") );           //Cooling pipe characteristics (see constructDiskCalorimeter)
       calo_->caloInfo_.set("rInnerManifold",         config.getDouble("calorimeter.rInnerManifold") );      //Cooling pipe characteristics (see constructDiskCalorimeter)
       calo_->caloInfo_.set("pipeTorRadius",          getVDouble(config,"pipeTorRadius", calo_->caloInfo_.getInt("nPipes")) ); //Tor radius of cooling pipes and other
       calo_->caloInfo_.set("largeTorPhi",            getVDouble(config,"largeTorPhi",   calo_->caloInfo_.getInt("nPipes")) ); //cooling pipe characteristics
       calo_->caloInfo_.set("smallTorPhi",            getVDouble(config,"smallTorPhi",   calo_->caloInfo_.getInt("nPipes")) ); //(see constructDiskCalorimeter)
       calo_->caloInfo_.set("yposition",              getVDouble(config,"yposition",     calo_->caloInfo_.getInt("nPipes")) );
       calo_->caloInfo_.set("straightEndPhi",         getVDouble(config,"straightEndPhi",calo_->caloInfo_.getInt("nPipes")) );

       calo_->caloInfo_.set("nCrates",                config.getInt   ("calorimeter.nCrates") );              //number of FEB crates per disk
       calo_->caloInfo_.set("nBoards",                config.getInt   ("calorimeter.nBoards") );              //number of boards per crate
       calo_->caloInfo_.set("FEBToDiskZOffset",       config.getDouble("calorimeter.FEBToDiskZOffset") );     //Z offset between front face of FEB and front face of disk
       calo_->caloInfo_.set("crateXLength",           config.getDouble("calorimeter.crateXLength") );         //X length of crate
       calo_->caloInfo_.set("crateYLength",           config.getDouble("calorimeter.crateYLength") );         //Y length of crate
       calo_->caloInfo_.set("crateZLength",           config.getDouble("calorimeter.crateZLength") );         //Z length of crate
       calo_->caloInfo_.set("crateFShieldThickness",  config.getDouble("calorimeter.crateFShieldThickness") );//Front shield thickness
       calo_->caloInfo_.set("crateBShieldThickness",  config.getDouble("calorimeter.crateBShieldThickness") );//Bottom shield thickness
       calo_->caloInfo_.set("crateBShieldLength",     config.getDouble("calorimeter.crateBShieldLength") );   //Bottom shield length
       calo_->caloInfo_.set("crateTThickness",        config.getDouble("calorimeter.crateTThickness") );      //Thickness of top crate plate
       calo_->caloInfo_.set("crateSThickness",        config.getDouble("calorimeter.crateSThickness") );      //Thickness of side crate plate
       calo_->caloInfo_.set("crateFShieldYLength",    config.getDouble("calorimeter.crateFShieldYLength") );  //Y length front shield
       calo_->caloInfo_.set("crateFShieldDeltaZ",     config.getDouble("calorimeter.crateFShieldDeltaZ") );   //Z difference betrween front shield and crate
       calo_->caloInfo_.set("radiatorThickness",      config.getDouble("calorimeter.radiatorThickness") );    //radiator thickness
       calo_->caloInfo_.set("radiatorZLength",        config.getDouble("calorimeter.radiatorZLength") );      //radiator Z length
       calo_->caloInfo_.set("activeStripThickness",   config.getDouble("calorimeter.activeStripThickness") ); //active strip thickness
       calo_->caloInfo_.set("passiveStripThickness",  config.getDouble("calorimeter.passiveStripThickness") );//passive strip thickness
       calo_->caloInfo_.set("cratePhiAngles",         getVDouble(config,"cratePhiAngles",calo_->caloInfo_.getInt("nCrates")) ); //phi angle of the crates


       // CALORIMETER ORIGIN AND FRONT FACE (FF) COORDINATES
       double zTrackerCenter  =  config.getDouble("mu2e.detectorSystemZ0");
       double xOrigin         = -config.getDouble("mu2e.solenoidOffset");
       double zCaloStart      =  config.getDouble("calorimeter.caloMotherZ0");
       double zCaloEnd        =  config.getDouble("calorimeter.caloMotherZ1");

       calo_->geomUtil_.origin( CLHEP::Hep3Vector(xOrigin,0,0.5*(zCaloStart+zCaloEnd)) );
       calo_->geomUtil_.trackerCenter(CLHEP::Hep3Vector(xOrigin,0,zTrackerCenter));
       calo_->geomUtil_.crystalZLength(config.getDouble("calorimeter.crystalZLength"));

       checkIt();
       makeIt();
    }

    void DiskCalorimeterMaker::makeIt(void)
    {
        double vdThickness           = calo_->caloInfo_.getDouble("vdThickness");
        double caloMotherZ0          = calo_->caloInfo_.getDouble("caloMotherZ0");
        double caloMotherZ1          = calo_->caloInfo_.getDouble("caloMotherZ1");
        double crystalHalfXY         = calo_->caloInfo_.getDouble("crystalXYLength")/2.0;
        double wrapperHalfThick      = calo_->caloInfo_.getDouble("wrapperThickness")/2.0;
        double crystalHalfZLength    = calo_->caloInfo_.getDouble("crystalZLength")/2.0;
        double crystalCapHalfZLength = calo_->caloInfo_.getDouble("crystalCapZLength")/2.0;
        double crystalCellRadius     = crystalHalfXY + 2.0*wrapperHalfThick;
        double innerDiskRadius       = calo_->caloInfo_.getDouble("diskInAlRingRIn");
        double innerCrysRadius       = calo_->caloInfo_.getDouble("diskCrystalRIn");
        double outerCrysRadius       = calo_->caloInfo_.getDouble("diskCrystalROut");
        double outerCaseRadius       = calo_->caloInfo_.getDouble("diskCaseRingROut");
        double diskOutRailROut       = calo_->caloInfo_.getDouble("diskOutRailROut");
        auto   shimStepsInRowId      = calo_->caloInfo_.getVInt("shimStepsInRowId");
        auto   shimStepsOutRowId     = calo_->caloInfo_.getVInt("shimStepsOutRowId");
        double FPCarbonThick         = calo_->caloInfo_.getDouble("FPCarbonZLength");
        double FPFoamThick           = calo_->caloInfo_.getDouble("FPFoamZLength");
        double FPCoolPipeRadius      = calo_->caloInfo_.getDouble("FPCoolPipeRadius");
        double FPpipeRadius          = calo_->caloInfo_.getDouble("pipeRadius");
        double diskCaseHalfZLength   = crystalHalfZLength + crystalCapHalfZLength;
        double BPHoleHalfZ           = calo_->caloInfo_.getDouble("BPHoleZLength")/2.0;
        double FEEBoxHalfZ           = calo_->caloInfo_.getDouble("FEEZLength")/2.0;
        double FEEBoxThick           = calo_->caloInfo_.getDouble("FEEBoxThickness");
        double BPPipeRadiusHigh      = calo_->caloInfo_.getDouble("BPPipeRadiusHigh");
        double BPPipeHalfZOffset     = calo_->caloInfo_.getDouble("BPPipeZOffset")/2.0;
        double crateZLength          = calo_->caloInfo_.getDouble("crateZLength");
        double crateFShieldThick     = calo_->caloInfo_.getDouble("crateFShieldThickness");
        double crateFShieldDeltaZ    = calo_->caloInfo_.getDouble("crateFShieldDeltaZ");
        double FEBToDiskZOffset      = calo_->caloInfo_.getDouble("FEBToDiskZOffset");
        double motherHalfZ           = (caloMotherZ1 - caloMotherZ0)/2.0;


        // First, calculate the total z length of the disk and the feb since we are not allowed
        // to get this from constructDiskCalorimeter.cc (no dependency on MC simulation)
        // Make sure this matches the geometry implemented in constructDiskCalorimeter !!!
        //
        double FPHalfZLength     = (FPCarbonThick + FPFoamThick - FPpipeRadius + FPCoolPipeRadius)/2.0;
        double BPHalfZLength     = BPHoleHalfZ + FEEBoxHalfZ + 2.0*FEEBoxThick + BPPipeHalfZOffset + BPPipeRadiusHigh;
        double diskHalfZLength   = FPHalfZLength + diskCaseHalfZLength + BPHalfZLength + vdThickness;
        double FEBZLength        = (crateZLength + crateFShieldDeltaZ + crateFShieldThick + 2.0*vdThickness);

        // Offset between front face of FEB and front face of disk
        double crateToDiskDeltaZ = FEBToDiskZOffset + vdThickness;

        // Offsets between the disk and crystal cordinate systems, i.e. distance between center of disk and front face crystals
        double disp = -diskHalfZLength + vdThickness + 2*FPHalfZLength + diskCaseHalfZLength - crystalHalfZLength + crystalCapHalfZLength;
        auto   diskOriginToCrystalOrigin = CLHEP::Hep3Vector(0,0,disp);


        // Then create the disks
        for (int idisk=0; idisk<CaloConst::_nDisk; ++idisk) {
           size_t crystalOffset = calo_->fullCrystalList_.size();
           double separation    = calo_->caloInfo_.getVDouble("diskZMotherShift").at(idisk);
           double angleZ        = 0;

           double dR1 = innerDiskRadius;
           double dR2 = diskOutRailROut;
           double dZ  = 2.0*diskHalfZLength;

           CLHEP::Hep3Vector size(dR1,dR2,dZ) ;
           CLHEP::Hep3Vector originLocal(0, 0, -motherHalfZ + diskHalfZLength + separation + crateToDiskDeltaZ);

           CLHEP::Hep3Vector frontFaceCenter = calo_->geomUtil_.origin() + originLocal + diskOriginToCrystalOrigin;
           CLHEP::Hep3Vector backFaceCenter  = frontFaceCenter + CLHEP::Hep3Vector(0,0,2.0*crystalHalfZLength);
           CLHEP::HepRotation diskRotation   = CLHEP::HepRotation::IDENTITY*CLHEP::HepRotationZ(angleZ);

           auto thisDisk = std::make_shared<Disk>(idisk,innerCrysRadius,outerCrysRadius, 2.0*crystalCellRadius,
                                                  2.0*crystalHalfZLength, crystalOffset, diskOriginToCrystalOrigin);
           calo_->disks_.push_back(thisDisk);

           thisDisk->geomInfo().size(size);
           thisDisk->geomInfo().originLocal(originLocal);
           thisDisk->geomInfo().origin(calo_->geomUtil_.origin() + originLocal);
           thisDisk->geomInfo().originToCrystalOrigin(diskOriginToCrystalOrigin);
           thisDisk->geomInfo().rotation(diskRotation);
           thisDisk->geomInfo().frontFaceCenter(frontFaceCenter);
           thisDisk->geomInfo().backFaceCenter(backFaceCenter);
           thisDisk->geomInfo().FEBZOffset(crateToDiskDeltaZ);
           thisDisk->geomInfo().FEBZLength(FEBZLength);
           thisDisk->geomInfo().envelopeRad(dR1,dR2);
           thisDisk->geomInfo().crystalDirection(CLHEP::Hep3Vector(0,0,1));


           //fill the full Crystal List / diskId (direct access for performance optimization)
           for (unsigned icry=0;icry<thisDisk->nCrystals();++icry){
              Crystal& thisCrystal = thisDisk->crystal(icry);
              calo_->fullCrystalList_.push_back(&thisCrystal);

              //precompute the neighbors in the global frame
              thisCrystal.setNeighbors(calo_->neighborsByLevel(icry + crystalOffset,1));
              thisCrystal.setNextNeighbors(calo_->neighborsByLevel(icry + crystalOffset,2));

              //pre-compute the crystal position in the mu2e frame (aka global frame)
              CLHEP::Hep3Vector globalPosition = thisDisk->geomInfo().origin() + thisDisk->geomInfo().inverseRotation()*(thisCrystal.localPosition());
              thisCrystal.setPosition(globalPosition);
           }

           //calculate the position of the inner and outer steps, including the special shims to make the walls straight
           std::vector<double> par,stepsInX, stepsInY, stepsOutX, stepsOutY;
           int nrows = int(outerCaseRadius/2.0/crystalCellRadius)+1;
           for (int i=-nrows;i<nrows;++i) {
             thisDisk->boundingBoxes(i,par);
             if (par.empty()) continue;

             double p0(par[0]),p1(par[1]);
             if (std::find(shimStepsOutRowId.begin(),shimStepsOutRowId.end(),i) != shimStepsOutRowId.end()){
               p0 -= crystalCellRadius;
               p1 += crystalCellRadius;
             }
             stepsOutX.push_back(p0);
             stepsOutX.push_back(p1);
             stepsOutY.push_back(par[2]);
             stepsOutY.push_back(par[3]);

             if (par.size()==6) {
                double p4(par[4]),p5(par[5]);
                if (std::find(shimStepsInRowId.begin(),shimStepsInRowId.end(),i) != shimStepsInRowId.end()){
                  p4 += crystalCellRadius;
                  p5 -= crystalCellRadius;
                }
                stepsInX.push_back(p4);stepsInX.push_back(p5);
                stepsInY.push_back(par[2]);stepsInY.push_back(par[3]);
             }
           }
           calo_->caloInfo_.set("stepsInsideX",stepsInX);
           calo_->caloInfo_.set("stepsInsideY",stepsInY);
           calo_->caloInfo_.set("stepsOutsideX",stepsOutX);
           calo_->caloInfo_.set("stepsOutsideY",stepsOutY);

           if (verbosityLevel_) std::cout<<"Constructed Disk "<<thisDisk->id()<<":  Rin="<<thisDisk->geomInfo().innerEnvelopeR()
                                         <<"  Rout="<<thisDisk->geomInfo().outerEnvelopeR()
                                         <<" (X,Y,Z)="<<thisDisk->geomInfo().origin()<<"  local_(X,Y,Z)="<<thisDisk->geomInfo().originLocal()
                                         <<"  with "<<thisDisk->nCrystals()<<" crystals"<<std::endl;

           if (verbosityLevel_ > 1) thisDisk->print()                     ;
        }
    }



    void DiskCalorimeterMaker::checkIt(void)
    {

        //Consistency check between disk and CaloConst content
        for (const auto& diskPtr : calo_->diskPtrs()) {
           if (diskPtr->nCrystals() != CaloConst::_nCrystalPerDisk)
           throw cet::exception("DiskCaloGeom") << "The number of crystals ("<<diskPtr->nCrystals()
                                                <<")is differfent from "<<CaloConst::_nCrystalPerDisk<<"\n";
        }

        //check calorimeter fits inside mother envelope
        double diskRin       = calo_->caloInfo_.getDouble("diskInAlRingRIn");
        double diskRout      = calo_->caloInfo_.getDouble("diskCaseRingROut");
        double outerRingROut = calo_->caloInfo_.getDouble("diskOutRailROut");
        if (diskRin > diskRout)
              throw cet::exception("DiskCaloGeom") << "calorimeter.diskInnerRingIn > calorimeter.diskOuterRingOut \n";

        if (outerRingROut > calo_->caloInfo_.getDouble("caloDiskRadiusOut"))
              throw cet::exception("DiskCaloGeom") << "calorimeter outer radius larger than calorimeter mother \n";

        if (diskRin < calo_->caloInfo_.getDouble("caloDiskRadiusIn"))
              throw cet::exception("DiskCaloGeom") << "calorimeter inner radius smaller than calorimeter mother \n";


        //check that holes in back plate are smaller than crystal, RO smaller than holes and FEE boxes fit
        if (calo_->caloInfo_.getDouble("BPHoleXLength") > calo_->caloInfo_.getDouble("crystalXYLength") ||
            calo_->caloInfo_.getDouble("BPHoleYLength") > calo_->caloInfo_.getDouble("crystalXYLength") )
              throw cet::exception("DiskCaloGeom") << "calorimeter backplate hole greater than crystal dimensions in X or Y \n";

        if (calo_->caloInfo_.getDouble("readoutXLength") > calo_->caloInfo_.getDouble("BPHoleXLength")  ||
            calo_->caloInfo_.getDouble("readoutYLength") > calo_->caloInfo_.getDouble("BPHoleYLength"))
              throw cet::exception("DiskCaloGeom") << "calorimeter readout larger than hole in X or Y \n";

        if (calo_->caloInfo_.getDouble("readoutZLength") > calo_->caloInfo_.getDouble("BPHoleZLength"))
              throw cet::exception("DiskCaloGeom") << "calorimeter readout too thick to fit in hole \n";

        if (calo_->caloInfo_.getDouble("FEEXLength") > calo_->caloInfo_.getDouble("BPHoleXLength")/calo_->caloInfo_.getInt("nSiPMPerCrystal"))
              throw cet::exception("DiskCaloGeom") << "calorimeter FEE box does not fit in X direction \n";

        if (calo_->caloInfo_.getDouble("FEEYLength") > calo_->caloInfo_.getDouble("crystalXYLength")-2*calo_->caloInfo_.getDouble("FEEBoxThickness"))
              throw cet::exception("DiskCaloGeom") << "calorimeter FEE box does not fit in Y direction \n";


        //Check pipes
        for (int i=0;i<calo_->caloInfo().getInt("nPipes");++i){
           double pipeTorRadius     = calo_->caloInfo_.getVDouble("pipeTorRadius").at(i);
           double pipeRadius        = calo_->caloInfo_.getDouble("pipeRadius");
           double radiusIn          = calo_->caloInfo_.getDouble("FPInnerRadius");
           double radiusOut         = calo_->caloInfo_.getDouble("FPOuterRadius");
           double foamZLength       = calo_->caloInfo_.getDouble("FPFoamZLength");
           double carbonThick       = calo_->caloInfo_.getDouble("FPCarbonZLength");
           double foamThick         = calo_->caloInfo_.getDouble("FPFoamZLength");
           double coolPipeTorRadius = calo_->caloInfo_.getDouble("FPCoolPipeTorRadius");
           double coolPipeRadius    = calo_->caloInfo_.getDouble("FPCoolPipeRadius");

           if ( pipeTorRadius - pipeRadius < radiusIn)
                 throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is smaller than disk inner radius\n";

           if ( pipeTorRadius + pipeRadius > radiusOut)
                 throw cet::exception("DiskCaloGeom") << "element "<<i<<" of calorimeter.pipeTorRadius is larger than disk outer radius\n";

           if ( pipeRadius > foamZLength/2.0)
                 throw cet::exception("DiskCaloGeom") << "calorimeter pipe radius too large to fit inside Foam front panel\n";

           if ( carbonThick + foamThick - pipeRadius < coolPipeRadius)
                 throw cet::exception("DiskCaloGeom") << "calorimeter cooling pipe radius too large\n";

           if ( carbonThick + pipeRadius > coolPipeRadius)
                 throw cet::exception("DiskCaloGeom") << "calorimeter cooling pipe radius too small\n";

           if ( coolPipeTorRadius - coolPipeRadius < radiusOut)
                 throw cet::exception("DiskCaloGeom") << "cooling pipe too large, overlap with foam structure\n";
        }


        //Just a few checks on crates
        int nBoards           = calo_->caloInfo_.getInt("nBoards");
        double radiatorDY     = calo_->caloInfo_.getDouble("radiatorThickness")/2.0;
        double activeStripDY  = calo_->caloInfo_.getDouble("activeStripThickness")/2.0;
        double passiveStripDY = calo_->caloInfo_.getDouble("passiveStripThickness")/2.0;
        double crateXLength   = calo_->caloInfo_.getDouble("crateXLength");
        double crateYLength   = calo_->caloInfo_.getDouble("crateYLength");

        if ( nBoards*(radiatorDY+activeStripDY+passiveStripDY) > crateYLength)
              throw cet::exception("DiskCaloGeom") << "calorimeter FEB boards too thick\n";

        if (calo_->caloInfo_.getDouble("crateFShieldYLength") > crateYLength)
              throw cet::exception("DiskCaloGeom") << "calorimeter FEB front shile too long in Y direction\n";

        if (calo_->caloInfo_.getDouble("crateSThickness") > crateXLength)
              throw cet::exception("DiskCaloGeom") << "calorimeter FEB crate side too thick\n";

        if (calo_->caloInfo_.getDouble("crateTThickness") > crateYLength )
              throw cet::exception("DiskCaloGeom") << "calorimeter FEB crate top too thick\n";

    }


    std::vector<double> DiskCalorimeterMaker::getVDouble(const SimpleConfig& config, const std::string& key, int size)
    {
       std::vector<double> vec;
       if (size>0 ) config.getVectorDouble("calorimeter."+key, vec, size);
       else         config.getVectorDouble("calorimeter."+key, vec);
       return vec;
    }

    std::vector<int> DiskCalorimeterMaker::getVInt(const SimpleConfig& config, const std::string& key, int size)
    {
       std::vector<int> vec;
       if (size>0) config.getVectorInt("calorimeter."+key, vec, size);
       else        config.getVectorInt("calorimeter."+key, vec);
       return vec;
    }


}

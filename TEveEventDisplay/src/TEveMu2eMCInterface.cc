
#include "Offline/TEveEventDisplay/src/TEveMu2e_base_classes/TEveMu2eMCInterface.h"

using namespace mu2e;
namespace mu2e{

  /*------------Function to clear lists for new events:-------------*/
  template <typename T, typename U> void DataLists(T data, bool Redraw, bool accumulate, string title, TEveElementList **List3D, TEveElementList **List2DXY = 0, TEveElementList **List2DXZ = 0, U projection = 0){	
      if(data == 0 && Redraw){
        if (*List3D != 0){
          (*List3D)->DestroyElements();
        }

        if (*List2DXY != 0){
          (*List2DXY)->DestroyElements();

          projection->fXYMgr->ImportElements(*List2DXY, projection->fDetXYScene); 

        }if (*List2DXZ != 0){
          (*List2DXZ)->DestroyElements();

          projection->fRZMgr->ImportElements(*List2DXZ, projection->fDetRZScene);
        }
        gEve->AddElement(*List3D);
        gEve->Redraw3D(kTRUE); 
      } 
      if(data!=0){
        if (*List3D== 0) {
          *List3D = new TEveElementList((title + "3D").c_str());
          if(!accumulate){(*List3D)->DestroyElements();} if(!accumulate){(*List3D)->DestroyElements();}     
        }
        else {
          (*List3D)->DestroyElements();  
        }
        if (*List2DXY== 0) {
          *List2DXY = new TEveElementList((title + "2D").c_str());
          (*List2DXY)->IncDenyDestroy();     
        }
        else {
          if (!accumulate){(*List2DXY)->DestroyElements();}    
        }
      }
      if (*List2DXZ== 0) {
        *List2DXZ = new TEveElementList((title + "2D").c_str());
        (*List2DXZ)->IncDenyDestroy();     
      }
      else {
        if (!accumulate){(*List2DXZ)->DestroyElements();}    
      }
    
  }
  

  /*------------Function to add straight line MC Trajectory i.e. for Comsics in No field:-------------*/
  void TEveMu2eMCInterface::AddSimpleMCTrajectory(bool firstloop, const MCTrajectoryCollection *trajcol, TEveMu2e2DProjection *tracker2Dproj, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2){    
        DataLists<const MCTrajectoryCollection*, TEveMu2e2DProjection*>(trajcol, Redraw, accumulate, "MC Trajectory", &fTrackList3D, &fTrackList2DXY,&fTrackList2DXZ, tracker2Dproj);
    if(trajcol!=0){
      TEveElementList *HitList3D = new TEveElementList("MCtraj3D");
      TEveElementList *HitList2D = new TEveElementList("MCtraj2D");
      std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
      for(trajectoryIter=trajcol->begin(); trajectoryIter!=trajcol->end(); trajectoryIter++)
      {
        const std::vector<MCTrajectoryPoint> &points = trajectoryIter->second.points();
        string pdgId= to_string(trajectoryIter->first->pdgId());
        CLHEP::Hep3Vector StartHitPos(points[0].x(), points[0].y(), points[0].z());
        CLHEP::Hep3Vector EndHitPos(points[points.size()-1].x(), points[points.size()-1].y(), points[points.size()-1].z());
        TEveMu2eMCTraj *teve_hit3D = new TEveMu2eMCTraj();
        string energy = to_string(points[0].kineticEnergy());
        teve_hit3D->DrawSimpleLine("MCTraj PDG " + pdgId + "Energy = " + energy  + ", ",  StartHitPos, EndHitPos, HitList3D);

        fTrackList3D->AddElement(HitList3D);

        GeomHandle<DetectorSystem> det;
        StartHitPos = det->toMu2e(StartHitPos);
        EndHitPos = det->toMu2e(EndHitPos);
        TEveMu2eMCTraj *teve_hit2D = new TEveMu2eMCTraj();
        teve_hit2D->DrawSimpleLine("MCTraj PDG " + pdgId + "Energy = " + energy + ", ", StartHitPos, EndHitPos, HitList2D);
        fTrackList2DXZ->AddElement(HitList2D); 

      }
      
      tracker2Dproj->fXYMgr->ImportElements(fTrackList2DXZ, tracker2Dproj->fEvtXYScene); 
      tracker2Dproj->fRZMgr->ImportElements(fTrackList2DXZ, tracker2Dproj->fEvtRZScene);
      
      gEve->AddElement(fTrackList3D);
      gEve->Redraw3D(kTRUE);
      }
  
  }
  
  /*------------Function to help user select a list of PDG codes to display:-------------*/
  int TEveMu2eMCInterface::Contains(std::vector<int> v, int x)
  {
    return std::count(v.begin(), v.end(), x);
  }
  

    /*------------Function to make label :-------------*/
    TEveText *TEveMu2eMCInterface::GetLabel(int PDGCode, TEveMu2eCustomHelix *line, TEveMu2eCustomHelix *line_twoDXY, TEveMu2eCustomHelix *line_twoDXZ){
      const char* pid = "pid";
      auto t = new TEveText(pid);
      Color_t color;
      double posy = 0;
      double posz = 3000.0;
      t->SetFontSize(15);
      switch(PDGCode) {
          case 11:
              color = 2;
              pid = "electron -";
              posy = 1400.0;
              break;
          case -11:
              color = 5;
              pid = "positron +";
              posy = 1500.0;
              break;
          case 13:
              color = 3;
              pid = "muon - ";
              posy = 1600.0;
              break;
          case -13:
              color = 46;
              pid = "muon + ";
              posy = 1700.0;
              break;
          case -211:
              color = 6;
              pid = "pion -";
              posy = 1800.0;
              break;
          case 211:
              color = 9;
              pid = "pion +";
              posy = 1900.0;
              break;
          case 2212:
              color = 4;
              pid = "proton";
              posy = 2000.0;
              break;
          case 22:
              color = 8;
              pid = "gamma";
              posy = 2100.0;
              break;
          default:
              color = 7;
              pid = "other";
              posy = 2200.0;
              break;
          }
      t->SetText(pid);
      t->SetMainColor(color);
      line->SetLineColor(color);
      line_twoDXZ->SetLineColor(color);
      line_twoDXY->SetLineColor(color);
      t->RefMainTrans().SetPos(0.0,posy,posz);
      
      return t;
    }
    
      /*------------
  Function to display MCTracjories of any shape, these are made up of a series of TEveLines, in the same way as Reco Helices:      
  -------------*/
    void TEveMu2eMCInterface::AddFullMCTrajectory(bool firstloop, const MCTrajectoryCollection *trajcol, TEveMu2e2DProjection *tracker2Dproj, bool Redraw, bool accumulate, TEveProjectionManager *TXYMgr, TEveProjectionManager *TRZMgr, TEveScene *scene1, TEveScene *scene2, std::vector<int> particleIds){

        if(trajcol!=0){
          DataLists<const MCTrajectoryCollection*, TEveMu2e2DProjection*>(trajcol, Redraw, accumulate, "MC Trajectory", &fTrackList3D, &fTrackList2DXZ, &fTrackList2DXY, tracker2Dproj);
          std::map<art::Ptr<mu2e::SimParticle>,mu2e::MCTrajectory>::const_iterator trajectoryIter;
          for(trajectoryIter=trajcol->begin(); trajectoryIter!=trajcol->end(); trajectoryIter++)
          { 
            TEveMu2eCustomHelix *line = new TEveMu2eCustomHelix();
            TEveMu2eCustomHelix *line_twoDXZ = new TEveMu2eCustomHelix();
            TEveMu2eCustomHelix *line_twoDXY = new TEveMu2eCustomHelix();
            //check user defined list of particles to plot:
            int x = Contains(particleIds,trajectoryIter->first->pdgId()); 
            if(x == 1){
              
              const std::vector<MCTrajectoryPoint> &points = trajectoryIter->second.points();
              
              for(unsigned int i=0; i<points.size();i++){
                CLHEP::Hep3Vector Pos(points[i].x(), points[i].y(), points[i].z());
                GeomHandle<DetectorSystem> det;
                CLHEP::Hep3Vector HitPos2D = det->toDetector(Pos);
                                      
                if(i==0) {
                      line->SetPoint(i,(Pos.x()), (Pos.y()),(Pos.z())); 
                      line_twoDXZ->SetPoint(i,pointmmTocm(HitPos2D.x()), pointmmTocm(HitPos2D.y())+1000,pointmmTocm(HitPos2D.z()));
                      line_twoDXY->SetPoint(i,pointmmTocm(HitPos2D.x()), pointmmTocm(HitPos2D.y()),pointmmTocm(HitPos2D.z()));      
                      
                } else {
                    line->SetNextPoint((Pos.x()),(Pos.y()),(Pos.z()));
                    line_twoDXZ->SetNextPoint(pointmmTocm(HitPos2D.x()), pointmmTocm(HitPos2D.y())+1000,pointmmTocm(HitPos2D.z()));
                    line_twoDXY->SetNextPoint(pointmmTocm(HitPos2D.x()), pointmmTocm(HitPos2D.y()),pointmmTocm(HitPos2D.z()));
                }
              }
              
              string energy = to_string(points[0].kineticEnergy());
   
              const std::string title = " MCTrajectory "+ energy + " Creation code = " + to_string(trajectoryIter->first->creationCode()) + "Stopping code = " + to_string(trajectoryIter->first->stoppingCode()) + " End Global Time = " + to_string(trajectoryIter->first->endGlobalTime())  ;
              line->SetTitle(Form(title.c_str()));
              
              //Get PID label:
              TEveText *t = GetLabel(trajectoryIter->first->pdgId(), line, line_twoDXZ, line_twoDXY);
              line_twoDXZ->SetTitle(Form(title.c_str()));
              line_twoDXY->SetTitle(Form(title.c_str()));
              line->SetLineWidth(3);
              line->SetPickable(kTRUE);
              line_twoDXZ->SetLineWidth(3);
              line_twoDXZ->SetPickable(kTRUE);
              line_twoDXY->SetLineWidth(3);
              line_twoDXY->SetPickable(kTRUE);
              fTrackList2DXZ->AddElement(line_twoDXZ);
              fTrackList2DXY->AddElement(line_twoDXY);
              fTrackList3D->AddElement(line);
              fTrackList3D->AddElement(t);
              }
              
             else std::cout<<"Warning: No Particles of User-Specified Type In File "<<std::endl;
            }
            tracker2Dproj->fXYMgr->ImportElements(fTrackList2DXY, tracker2Dproj->fEvtXYScene); 
            tracker2Dproj->fRZMgr->ImportElements(fTrackList2DXZ, tracker2Dproj->fEvtRZScene);
            gEve->AddElement(fTrackList3D);
            gEve->Redraw3D(kTRUE);
        }
    }
}

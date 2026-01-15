#include "Offline/TrkHitReco/inc/DBSClusterer.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "TMath.h"
#include <algorithm>
#include <vector>
#include <queue>


namespace mu2e
{
  DBSClusterer::DBSClusterer(const std::optional<Config> config) :
    DBSminExpand_     (config.value().DBSminN()),
    deltaTime_        (config.value().hitDeltaTime()),
    deltaZ_           (config.value().hitDeltaZ()),
    deltaXY2_         (config.value().hitDeltaXY()*config.value().hitDeltaXY()),
    minClusterHits_   (config.value().minClusterHits()),
    bkgmask_          (config.value().bkgmsk()),
    sigmask_          (config.value().sigmsk()),
    testflag_         (config.value().testflag()),
    kerasW_           (config.value().kerasWeights()),
    diag_             (config.value().diag())
  {
  }


  //---------------------------------------------------------------------------------------
  void DBSClusterer::init() {
    //Add Classifier init here (see TNTClsuterer for example)
    /*ConfigFileLookupPolicy configFile;
     auto kerasWgtsFile = configFile(kerasW_);
     sofiePtr_          = std::make_shared<TMVA_SOFIE_TrainBkgDiag::Session>(kerasWgtsFile);*/
  }


  //----------------------------------------------------------------------------------------------------------
  void DBSClusterer::findClusters(BkgClusterCollection& clusters, const ComboHitCollection& chcol)
  {
     if (chcol.empty()) return;

     std::vector<unsigned> idx; //list of combo hit IDs
     idx.reserve(chcol.size());
     //std::cout<<"chcol size = "<<chcol.size()<<std::endl;
     for (size_t ich=0;ich<chcol.size();++ich) {
       if (testflag_ && (!chcol[ich].flag().hasAllProperties(sigmask_) || chcol[ich].flag().hasAnyProperty(bkgmask_))) continue;
       //std::cout<<"chcol "<<ich<<" time = "<<chcol[ich].correctedTime()<<" pos x = "<<chcol[ich].pos().x()<<" y = "<<chcol[ich].pos().y()
       //         <<" z = "<<chcol[ich].pos().z()<<std::endl;
       idx.emplace_back(ich);
     }
     //Sort combo hits which are not flagged as background according to their corrected Time
     std::sort(idx.begin(),idx.end(),[&chcol](auto i, auto j){return chcol[i].correctedTime() < chcol[j].correctedTime();});


     //--- DBScan algorithm
     const unsigned        noiseID(chcol.size()+1u);
     const unsigned        unprocessedID(chcol.size()+2u);
     unsigned              currentClusterID(0);
     std::queue<unsigned>  inspect;
     std::vector<unsigned> hitToCluster(idx.size(),unprocessedID); //number of combohits used for clustering, cluster ID
     std::vector<unsigned> neighbors;
     neighbors.reserve(32);
     for (size_t i=0;i<idx.size();++i) {

       // if a point has already been assigned to a cluster, continue
       int nNeighbors = 0;
       if ( hitToCluster[i] != unprocessedID) continue;

       // If the neighborhood is too sparse, assign it to noise
       nNeighbors = findNeighbors(i, idx, chcol, neighbors);
       //std::cout<<"1. nNeighbors = "<<nNeighbors<<" size = "<<neighbors.size()<<std::endl;
       //std::cout<<"nNeighbors = "<<nNeighbors<<std::endl;
       if (nNeighbors < DBSminExpand_) {
         hitToCluster[i] = noiseID;
         std::cout<<noiseID<<std::endl;
         continue;
       }

       hitToCluster[i] = currentClusterID;
       //std::cout<<"CurrentClusterID = "<<currentClusterID<<std::endl;
       BkgCluster thisCluster;
       thisCluster.addHit(idx[i]);
       // Extend the cluster by adding / expanding around neighbors
       for (const auto& j : neighbors) inspect.push(j);

       while (!inspect.empty()){
         auto j = inspect.front();
         inspect.pop();

         if (hitToCluster[j] == noiseID) {
           hitToCluster[j] = currentClusterID;
           thisCluster.addHit(idx[j]);
         }
         if (hitToCluster[j] != unprocessedID) continue;

         hitToCluster[j] = currentClusterID;
         thisCluster.addHit(idx[j]);

         nNeighbors = findNeighbors(j,idx,chcol,neighbors);
         //std::cout<<"2. nNeighbors = "<<nNeighbors<<" size = "<<neighbors.size()<<std::endl;
         //std::cout<<"nNeighbors in loop = "<<nNeighbors<<std::endl;
         if (nNeighbors < DBSminExpand_) continue;

         for (const auto& k : neighbors) {
           if (hitToCluster[k]==unprocessedID || hitToCluster[k]==noiseID) inspect.push(k);
         }
       }
       //std::cout<<"thisCluster hits size in loop = "<<thisCluster.hits().size()<<std::endl;
       if (thisCluster.hits().size() < minClusterHits_) continue; //minClusterHits_ = 1
       std::cout<<"currentClusterID = "<<currentClusterID<<" final neighbors = "<<nNeighbors<<" size = "<<neighbors.size()
                <<" i = "<<i<<" cluster hits size = "<<thisCluster.hits().size()<<std::endl;
       clusters.push_back(std::move(thisCluster));
       ++currentClusterID;
     }

     /*
     //Make noise hits into their own cluster if we need to
     for (size_t i=0;i<idx.size();++i) {
       if (hitToCluster[i] != noiseID) continue;
       clusters.emplace_back(BkgCluster());
       clusters.back().addHit(idx[i]);
     }
     */

     //Calculate the cluster properties
     for (auto& cluster : clusters) calculateCluster(cluster, chcol);

     /*for (size_t i = 0; i < idx.size(); ++i) {
       const auto& ch = chcol[idx[i]];
       std::cout << "i=" << i << " chIdx=" << idx[i] << " clu=" << hitToCluster[i] << " t=" << ch.correctedTime() << " x=" << ch.pos().x() << " y=" << ch.pos().y() << " z=" << ch.pos().z() << "\n";
     }*/

     //mergeClusters(clusters, chcol);
     //std::cout<<"End of calculate clusters"<<std::endl;
     if (diag_>1) dump(clusters);
  }



  //---------------------------------------------------------------------------------------
  // Find the neighbors of given a point - can use any suitable distance function
  int DBSClusterer::findNeighbors(unsigned ihit, const std::vector<unsigned>& idx, const ComboHitCollection& chcol, std::vector<unsigned>& neighbors)
  {
    // Define number of neighbors as number of straw hits in the neighboring point.
    // Commneted version is number of neighbor is 1 regardless of number of straw hits

    //unsigned nNeighbors(0);
    unsigned nNeighbors(chcol[idx[ihit]].nStrawHits()-1);
    neighbors.clear();
    float time0 = chcol[idx[ihit]].correctedTime();
    float x0    = chcol[idx[ihit]].pos().x();
    float y0    =  chcol[idx[ihit]].pos().y();
    float z0    = chcol[idx[ihit]].pos().z();
    float phi0  = chcol[idx[ihit]].pos().phi();
    //std::cout<<"time0 = "<<time0<<" x0 = "<<x0<<" y0  "<<y0<<" z0 "<<z0<<" phi0 = "<<phi0<<" nSH = "<<chcol[idx[ihit]].nStrawHits()<<std::endl;
    unsigned istart(ihit);
    /*for (size_t j=istart; j<idx.size(); ++j){
        if (j==ihit) continue;
        float deltat = chcol[idx[j]].correctedTime() - time0;     // deltaTime = 15
        float deltaz = chcol[idx[j]].pos().z()- z0;              //deltaZ = 800
        float dist = std::sqrt((chcol[idx[j]].pos().x()-x0)*(chcol[idx[j]].pos().x()-x0) +
                               (chcol[idx[j]].pos().y()-y0)*(chcol[idx[j]].pos().y()-y0));
        float dphi = chcol[idx[j]].pos().phi() - phi0;
        if(deltat < deltaTime_){
          std::cout<<"delta t = "<<deltat<<" delta z = "<<deltaz<<" delta XY = "<<dist<<" delta phi = "<<dphi<<" nSH = "<<chcol[idx[j]].nStrawHits()<<std::endl;
          std::cout<<"j hit time = "<<chcol[idx[j]].correctedTime()<<" x = "<<chcol[idx[j]].pos().x()<<" y = "<<chcol[idx[j]].pos().y()<<" z = "<<chcol[idx[j]].pos().z()<<" nSH = "<<chcol[idx[j]].nStrawHits()<<std::endl;
          }
    }*/
    //what is istart here? and why are we incrementing or decrementing it? I know it is the hit
    //we are comparing right now with all the other hits to look for neighbours.
    while (istart>0 && time0 - chcol[idx[istart]].correctedTime() < deltaTime_) --istart;
    if    (time0 - chcol[idx[istart]].correctedTime()             > deltaTime_) ++istart;
    //std::cout<<"idx size = "<<idx.size()<<std::endl;
    for (size_t j=istart; j<idx.size(); ++j){
      if (j==ihit) continue;
      if (chcol[idx[j]].correctedTime() - time0 > deltaTime_) break;     // deltaTime = 15
      if (std::abs(chcol[idx[j]].pos().z()- z0) > deltaZ_)    continue; //deltaZ = 800
      //std::cout<<"findNeighbors::deltaTime = "<<chcol[idx[j]].correctedTime() - time0<<" deltaZ = "<<std::abs(chcol[idx[j]].pos().z()- z0)<<std::endl;
      float dist = (chcol[idx[j]].pos().x()-x0)*(chcol[idx[j]].pos().x()-x0) +
                   (chcol[idx[j]].pos().y()-y0)*(chcol[idx[j]].pos().y()-y0);
      float dphi = chcol[idx[j]].pos().phi() - phi0;
      if (dphi > M_PI) dphi = 2*M_PI-dphi;
      //std::cout<<"delta t = "<<chcol[idx[j]].correctedTime() - time0<<" delta z = "<<chcol[idx[j]].pos().z()- z0<<" delta XY = "<<dist<<" dphi = "<<dphi<<" nSH = "<<chcol[idx[j]].nStrawHits()<<std::endl;
      //std::cout<<"findNeighbors::dist = "<<dist<<std::endl;
      if (dist > deltaXY2_) continue; //or std::abs(dphi) > 0.1) continue; //deltaXY2 = 2500
      //std::cout<<"Selected delta t = "<<chcol[idx[j]].correctedTime() - time0<<" delta z = "<<chcol[idx[j]].pos().z()- z0<<" delta XY = "<<dist<<" dphi = "<<dphi<<" nSH = "<<chcol[idx[j]].nStrawHits()<<" deltaXY2 = "<<deltaXY2_<<std::endl;
      //if(dist < deltaXY2_ or std::abs(dphi) < 0.4){
        //std::cout<<"Selected j hit time = "<<chcol[idx[j]].correctedTime()<<" x = "<<chcol[idx[j]].pos().x()<<" y = "<<chcol[idx[j]].pos().y()<<" z = "<<chcol[idx[j]].pos().z()<<" nSH = "<<chcol[idx[j]].nStrawHits()<<std::endl;
        neighbors.emplace_back(j);
        nNeighbors += chcol[idx[j]].nStrawHits();
        //}
      //++nNeighbors;
    }
    /*if(nNeighbors < 2) {
      std::cout<<"time0 = "<<time0<<" x0 = "<<x0<<" y0  "<<y0<<" z0 "<<z0<<std::endl;
      for (size_t j=istart; j<idx.size(); ++j){
        if (j==ihit) continue;
        float deltat = chcol[idx[j]].correctedTime() - time0;     // deltaTime = 15
        float deltaz = chcol[idx[j]].pos().z()- z0;              //deltaZ = 800
        float dist = (chcol[idx[j]].pos().x()-x0)*(chcol[idx[j]].pos().x()-x0) +
                   (chcol[idx[j]].pos().y()-y0)*(chcol[idx[j]].pos().y()-y0);
        std::cout<<"delta t = "<<deltat<<" delta z = "<<deltaz<<" delta XY = "<<dist<<std::endl;
      }
      }*/
    return nNeighbors;
  }


  //---------------------------------------------------------------------------------------
  // this is only used for diagnosis at this point
  float DBSClusterer::distance(const BkgCluster& cluster, const ComboHit& hit) const
  {
    float psep_x = hit.pos().x()-cluster.pos().x();
    float psep_y = hit.pos().y()-cluster.pos().y();
    return sqrt(psep_x*psep_x+psep_y*psep_y);
    // alterntively, could use the chi2 distance
    //return std::sqrt(cluster.points().dChi2(TwoDPoint(hit.pos(),hit.uDir(),hit.uVar(),hit.vVar())))
    //         - std::sqrt(cluster.points().chisquared());
  }



  //---------------------------------------------------------------------------------------
  void DBSClusterer::calculateCluster(BkgCluster& cluster, const ComboHitCollection& chcol)
  {
    if (cluster.hits().empty()) {cluster.time(0.0f);cluster.pos(XYZVectorF(0.0f,0.0f,0.0f));return;}
    //std::cout<<"calculateClsuter:: cluster hits size = "<<cluster.hits().size()<<std::endl;
    if (cluster.hits().size()==1) {
      int idx = cluster.hits().at(0);
      cluster.time(chcol[idx].correctedTime());
      cluster.edep(chcol[idx].energyDep());
      cluster.pos(XYZVectorF(chcol[idx].pos().x(),chcol[idx].pos().y(),chcol[idx].pos().z()));
      XYZVectorF hitpos(chcol[idx].pos().x(), chcol[idx].pos().y(), chcol[idx].pos().z());
      cluster.addHitPosition(hitpos);
      //std::cout<<"SingleCluster:: time = "<<cluster.time()<<" pos = "<<cluster.pos().x()<<"cy = "<<cluster.pos().y()<<" cz = "<<cluster.pos().z()<<std::endl;
      return;
    }

    float sumWeight(0),crho(0),ctime(0), cz(0), cedep(0), cphi(0);
    //std::vector<float> _phicluster;
    //std::vector<float> _z;
    //_phicluster.clear();
    std::cout<<"cluster hit size = "<<cluster.hits().size()<<std::endl;
    for (auto& idx : cluster.hits()) {
      float weight = chcol[idx].nStrawHits();
      float dt     = chcol[idx].correctedTime();
      float dr     = sqrtf(chcol[idx].pos().perp2());
      float dz     = chcol[idx].pos().z();
      float edep   = chcol[idx].energyDep();
      XYZVectorF hitpos(chcol[idx].pos().x(), chcol[idx].pos().y(), chcol[idx].pos().z());
      cluster.addHitPosition(hitpos);
      //std::cout<<"Cluster:: pos = "<<hitpos.x()<<" cy = "<<hitpos.y()<<" cz = "<<hitpos.z()<<" time = "<<chcol[idx].correctedTime()<<std::endl;
      float dp     = chcol[idx].phi();
      if (dp > M_PI)  dp -= 2*M_PI;
      if (dp < -M_PI) dp += 2*M_PI;
      //_phicluster.push_back(dp);
      //_z.push_back(chcol[idx].pos().z());
      ctime    += dt*weight;
      crho     += dr*weight;
      cphi     += dp*weight;
      cz       += dz*weight;
      cedep    += edep*weight;
      sumWeight += weight;
    }
    /*std::cout<<"phicluster size = "<<_phicluster.size()<<std::endl;
    for (size_t j = 0; j < _phicluster.size(); ++j) {
      std::cout <<"phi = "<<_phicluster[j] << " z = "<<_z[j]<<std::endl;
      }*/
    //float diff = _phicluster.back() - _phicluster.front();
    //float zdiff = _z.back() - _z.front();
    /*float phiDiff = 0.0;
    float zDiff = 0.0;
    if (_phicluster.size() > 2) {
      for (size_t i = 1; i < _phicluster.size(); ++i) {
        phiDiff += (_phicluster[i] - _phicluster[i - 1]);
      }
      for (size_t i = 1; i < _z.size(); ++i) {
        zDiff += (_z[i] - _z[i - 1]);
      }
    }
    phiDiff = phiDiff / static_cast<float>(_phicluster.size() - 1);
    zDiff = zDiff / static_cast<float>(_z.size() - 1);*/
    //if(_phicluster.size() > 2)
    //  std::cout<<"Phi diff = "<<phiDiff<<" z diff = "<<zDiff<<" cluster size = "<<_phicluster.size()<<std::endl;
    cphi  /= sumWeight;
    crho  /= sumWeight;
    ctime /= sumWeight;
    cz    /= sumWeight;
    cedep /= sumWeight;
    //std::cout<<"cedp = "<<cedep<<" size = "<<cluster.hits().size()<<std::endl;
    cluster.time(ctime);
    cluster.pos(XYZVectorF(crho*cos(cphi),crho*sin(cphi),cz));
    cluster.edep(cedep);
    //std::cout<<"calculateCluster:: time = "<<cluster.time()<<" pos = "<<cluster.pos().x()<<" cz = "<<cz<<std::endl;
  }

  void DBSClusterer::mergeClusters(std::vector<BkgCluster>& clusters, const ComboHitCollection& chcol){
    std::vector<size_t> singleClusters;
    int i(0);
    for(auto& cluster : clusters){
      if(cluster.hits().size() < 2)
        singleClusters.push_back(i);
      i++;
    }
    std::cout<<"chcol size = "<<chcol.size()<<std::endl;
  }


  //---------------------------------------------------------------------------------------
  void DBSClusterer::classifyCluster(BkgCluster& cluster, const ComboHitCollection& chcol){

    //code logic to classify cluster with MVA
    // count hits and planes
    /*std::array<int,StrawId::_nplanes> hitplanes{0};
    for (const auto& chit : cluster.hits()) {
    const ComboHit& ch = chcol[chit];
    hitplanes[ch.strawId().plane()] += ch.nStrawHits();
    }

    unsigned npexp(0),np(0),nhits(0);
    int ipmin(0),ipmax(StrawId::_nplanes-1);
    while (hitplanes[ipmin]==0 && ipmin<StrawId::_nplanes) ++ipmin;
    while (hitplanes[ipmax]==0 and ipmax>0)                --ipmax;
    int fp(ipmin),lp(ipmin-1),pgap(0);
    for (int ip = ipmin; ip <= ipmax; ++ip) {
    npexp++; // should use TTracker to see if plane is physically present FIXME!
    if (hitplanes[ip]> 0){
      ++np;
      if(lp > 0 && ip - lp -1 > pgap)pgap = ip - lp -1;
      if(ip > lp)lp = ip;
      if(ip < fp)fp = ip;
      lp = ip;
    }
    nhits += hitplanes[ip];
    }
    if(nhits < 1 || np < 2) return;
    // find averages
    double sumEdep(0.);
    double sqrSumDeltaTime(0.);
    double sqrSumDeltaX(0.);
    double sqrSumDeltaY(0.);
    double sqrSumQual(0.);
    double sumPitch(0.);
    double sumYaw(0.);
    double sumwPitch(0.);
    double sumwYaw(0.);
    double sumEcc(0.);
    double sumwEcc(0.);
    unsigned nsthits(0.);
    unsigned nchits = cluster.hits().size();
    for (const auto& chit : cluster.hits()) {
      sumEdep         +=  chcol[chit].energyDep()/chcol[chit].nStrawHits();
      sqrSumDeltaX    += std::pow(chcol[chit].pos().x() - cluster.pos().x(),2);
      sqrSumDeltaY    += std::pow(chcol[chit].pos().y() - cluster.pos().y(),2);
      sqrSumDeltaTime += std::pow(chcol[chit].time() - cluster.time(),2);
      auto hdir        = chcol[chit].hDir();
      auto wecc        = chcol[chit].nStrawHits();
      sumEcc          += std::sqrt(1-(chcol[chit].vVar()/chcol[chit].uVar()))*wecc;
      sumwEcc         += wecc;
      if (chcol[chit].flag().hasAllProperties(StrawHitFlag::sline)){

        //quality of SLine fit
        sqrSumQual += std::pow(chcol[chit].qual(),2);

        //angle with Mu2e-Y
        double varPitch = std::pow(TMath::ACos(std::sqrt(chcol[chit].hcostVar())),2);
        double wPitch = 1/varPitch;
        double signPitch = hdir.Y()/std::abs(hdir.Y());
        sumPitch += signPitch*wPitch*hdir.theta();
        sumwPitch += wPitch;

        ROOT::Math::XYZVectorF z = {0,0,1};
        ROOT::Math::XYZVectorF dxdz = {hdir.X(),0,hdir.Z()};
        float magdxdz = std::sqrt(dxdz.Mag2());

        //angle with Mu2e-Z
        double varYaw = std::sqrt(chcol[chit].hphiVar() + varPitch);
        double wYaw = 1/varYaw;
        double signYaw = hdir.X()/std::abs(hdir.X());
        sumYaw += signYaw*wYaw*TMath::ACos(dxdz.Dot(z)/magdxdz);
        sumwYaw += wYaw;

        // # of stereo hits with SLine
        nsthits++;
      }
    }

    // fill mva input variables
    std::array<float,12> kerasvars;
    kerasvars[0]  = cluster.pos().Rho(); // cluster rho, cyl coor
    kerasvars[1]  = fp;// first plane hit
    kerasvars[2]  = lp;// last plane hit
    kerasvars[3]  = pgap;// largest plane gap without hits between planes with hits
    kerasvars[4]  = np;// # of planes hit
    kerasvars[5]  =  static_cast<float>(np)/static_cast<float>(lp - fp);// fraction of planes hit between first and last plane
    kerasvars[6]  = nhits;// sum of straw hits
    kerasvars[7]  = std::sqrt((sqrSumDeltaX+sqrSumDeltaY)/nchits);  // RMS of cluster rho
    kerasvars[8]  = std::sqrt(sqrSumDeltaTime/nchits);// RMS of cluster time
    kerasvars[9]  = nsthits > 0 ? sumPitch/sumwPitch : 0.;
    kerasvars[10] = nsthits > 0 ? sumYaw/sumwYaw : 0.;
    kerasvars[11] = sumEcc/sumwEcc;

    std::vector<float> kerasout = sofiePtr_->infer(kerasvars.data());
    cluster.setKerasQ(kerasout[0]);

    if (diag_>0)std::cout << "kerasout = " << kerasout[0] << std::endl;*/
    cluster.setKerasQ(-1.0);
  }





  //-------------------------------------------------------------------------------------------
  void DBSClusterer::dump(const std::vector<BkgCluster>& clusters)
  {
    int iclu(0);
    for (auto& cluster: clusters) {
      //std::cout<<"Cluster "<<iclu<<" "<<cluster.pos()<<" "<<cluster.time()<<"  "<<cluster.hits().size()<<"  - ";
      for (auto& hit : cluster.hits()) std::cout<<hit<<" ";
      std::cout<<std::endl;
      ++iclu;
    }
  }

}

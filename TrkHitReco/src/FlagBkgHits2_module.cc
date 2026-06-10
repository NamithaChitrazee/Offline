#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/ConditionsService/inc/ConditionsHandle.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/BkgCluster.hh"
#include "Offline/RecoDataProducts/inc/BkgClusterHit.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/TrkHitReco/inc/TrainBkgDiag.hxx"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

namespace TMVA_SOFIE_TrainBkgDiag { class Session; }

namespace mu2e
{

  class FlagBkgHits2 : public art::EDProducer
  {
    public:

      struct Config
      {
        using Name    = fhicl::Name;
        using Comment = fhicl::Comment;

        // module-level parameters
        fhicl::Atom<art::InputTag>    comboHitCollection{ Name("ComboHitCollection"),    Comment("ComboHit collection name") };
        fhicl::Atom<bool>             filterHits{         Name("FilterHits"),            Comment("Produce filtered ComboHit collection") };
        fhicl::Sequence<std::string>  backgroundMask{     Name("BackgroundMask"),        Comment("Bkg hit selection mask for output filtering") };
        fhicl::Atom<bool>             saveBkgClusters{    Name("SaveBkgClusters"),       Comment("Save bkg clusters") };
        fhicl::Atom<bool>             countProtons{       Name("CountProtons"),          Comment("Count protons") };
        fhicl::Atom<float>            minEdep{            Name("MinEdep"),               Comment("Min Edep") };
        fhicl::Atom<std::string>      outputLevel{        Name("OutputLevel"),           Comment("Level of the output ComboHitCollection") };

        // DBScan clustering parameters
        fhicl::Atom<int>              DBSminExpand{       Name("DBSminExpand"),          Comment("Min number neighbors for DBScan algo") };
        fhicl::Atom<float>            hitDeltaTime{       Name("DeltaTime"),             Comment("Max time difference between hits") };
        fhicl::Atom<float>            hitDeltaZ{          Name("DeltaZ"),                Comment("Max Z difference between hits") };
        fhicl::Atom<float>            hitDeltaXY{         Name("DeltaXY"),               Comment("Max XY difference between hits") };
        fhicl::Atom<unsigned>         minClusterHits{     Name("MinClusterHits"),        Comment("Min number hits in cluster") };
        fhicl::Sequence<std::string>  clusterBkgMask{     Name("ClusterBackgroundMask"), Comment("Bkg hit mask for clustering") };
        fhicl::Sequence<std::string>  clusterSigMask{     Name("ClusterSignalMask"),     Comment("Signal hit mask for clustering") };
        fhicl::Atom<bool>             testflag{           Name("TestFlag"),              Comment("Test hit flags during clustering") };
        fhicl::Atom<int>              diag{               Name("Diag"),                  Comment("Diagnosis level"), 0 };
      };

      explicit FlagBkgHits2(const art::EDProducer::Table<Config>& config);
      void beginJob() override;
      void produce(art::Event& event) override;

    private:
      // module-level members
      const art::ProductToken<ComboHitCollection> chtoken_;
      bool                                        filter_;
      bool                                        savebkg_;
      bool                                        countprotons_;
      float                                       minedep_;
      StrawHitFlag                                bkgmsk_;
      StrawIdMask::Level                          level_;
      float                                       kerasQ_;

      // DBScan clustering members
      int                                         DBSminExpand_;
      float                                       deltaTime_;
      float                                       deltaZ_;
      float                                       deltaXY2_;
      unsigned                                    minClusterHits_;
      StrawHitFlag                                clusterBkgMask_;
      StrawHitFlag                                clusterSigMask_;
      bool                                        testflag_;
      int                                         diag_;

      // module-level methods
      void classifyCluster    (BkgClusterCollection& bkgccol, StrawHitFlagCollection& chfcol, const ComboHitCollection& chcol, std::vector<int>& hitToClusterMap, size_t iclStart = 0) const;
      void countProton        (BkgClusterCollection& bkgccol, StrawHitFlagCollection& chfcol, const ComboHitCollection& chcol) const;

      // DBScan clustering methods
      void  findClusters      (BkgClusterCollection& clusters, const ComboHitCollection& chcol, std::vector<unsigned>& idx, float deltaTime, float deltaZ, float deltaXY2);
      float distance          (const BkgCluster& cluster,     const ComboHit& hit) const;
      int   findNeighbors     (unsigned ihit, const std::vector<unsigned>& idx, const ComboHitCollection& chcol, std::vector<unsigned>& neighbors, float deltaTime, float deltaZ, float deltaXY2);
      void  calculateCluster  (BkgCluster& cluster,           const ComboHitCollection& chcol);
      void  dump              (const std::vector<BkgCluster>& clusters);
  };


  FlagBkgHits2::FlagBkgHits2(const art::EDProducer::Table<Config>& config) :
    art::EDProducer{config},
    chtoken_{       consumes<ComboHitCollection>(config().comboHitCollection()) },
    filter_(        config().filterHits()),
    savebkg_(       config().saveBkgClusters()),
    countprotons_(  config().countProtons()),
    minedep_(       config().minEdep()),
    bkgmsk_(        config().backgroundMask()),
    DBSminExpand_(  config().DBSminExpand()),
    deltaTime_(     config().hitDeltaTime()),
    deltaZ_(        config().hitDeltaZ()),
    deltaXY2_(      config().hitDeltaXY() * config().hitDeltaXY()),
    minClusterHits_(config().minClusterHits()),
    clusterBkgMask_(config().clusterBkgMask()),
    clusterSigMask_(config().clusterSigMask()),
    testflag_(      config().testflag()),
    diag_(          config().diag())
  {
    produces<ComboHitCollection>();
    if (savebkg_) {
      produces<BkgClusterHitCollection>();
      produces<BkgClusterCollection>();
    }
    StrawIdMask mask(config().outputLevel());
    level_ = mask.level();
  }


  void FlagBkgHits2::beginJob()
  {
  }


  //------------------------------------------------------------------------------------------
  void FlagBkgHits2::produce(art::Event& event)
  {
    auto chH = event.getValidHandle(chtoken_);
    const ComboHitCollection& chcol = *chH.product();
    unsigned nch = chcol.size();
    BkgClusterCollection bkgccol;
    BkgClusterHitCollection bkghitcol;
    bkgccol.reserve(nch/2);
    if (savebkg_) bkghitcol.reserve(nch);

    // round 1: cluster flag-filtered hits with tight parameters
    std::vector<unsigned> r1Idx;
    r1Idx.reserve(nch);
    for (size_t ich = 0; ich < nch; ++ich) {
      if (testflag_ && (!chcol[ich].flag().hasAllProperties(clusterSigMask_) || chcol[ich].flag().hasAnyProperty(clusterBkgMask_))) continue;
      r1Idx.emplace_back(ich);
    }
    findClusters(bkgccol, chcol, r1Idx, deltaTime_, deltaZ_, deltaXY2_);

    StrawHitFlagCollection chfcol(nch);
    std::vector<int> hitToClusterMap(nch, -1);
    classifyCluster(bkgccol, chfcol, chcol, hitToClusterMap);

    unsigned nUnclustered = 0;
    for (size_t ich = 0; ich < nch; ++ich)
      if (hitToClusterMap[ich] == -1) ++nUnclustered;
    if (countprotons_)
      countProton(bkgccol, chfcol, chcol);

    auto chcol_out = std::make_unique<ComboHitCollection>();
    if (chfcol.size() > 0) {
      if (level_ == chcol.level()) {
        chcol_out->setSameParent(chcol);
        chcol_out->reserve(nch);
        for (size_t ich = 0; ich < nch; ++ich) {
          StrawHitFlag const& flag = chfcol[ich];
          if (!filter_ || !flag.hasAnyProperty(bkgmsk_)) {
            chcol_out->push_back(chcol[ich]);
            chcol_out->back()._flag.merge(flag);
          }
        }
      } else {
        auto pptr = chcol.parent(level_);
        auto const& chcol_p = *pptr;
        chcol_out->setSameParent(chcol_p);
        ComboHitCollection::SHIV shiv;
        chcol_out->reserve(chcol_p.size()*2);
        for (size_t ich = 0; ich < nch; ++ich) {
          shiv.clear();
          StrawHitFlag const& flag = chfcol[ich];
          if (!filter_ || !flag.hasAnyProperty(bkgmsk_)) {
            if (&chcol_p == chcol.fillStrawHitIndices(ich, shiv, level_)) {
              for (auto ishi : shiv) {
                chcol_out->push_back(chcol_p[ishi]);
                chcol_out->back()._flag.merge(flag);
              }
            } else {
              throw cet::exception("RECO") << "FlagBkgHits2: inconsistent ComboHits" << std::endl;
            }
          }
        }
        if ((!filter_) && chcol_out->size() != chcol_p.size())
          throw cet::exception("RECO") << "FlagBkgHits2: inconsistent ComboHit output" << std::endl;
      }
    }

    if (savebkg_) {
      for (size_t ich = 0; ich < nch; ++ich) {
        int icl = hitToClusterMap[ich];
        if (icl > -1) bkghitcol.emplace_back(BkgClusterHit(distance(bkgccol[icl], chcol[ich]), chfcol[ich]));
        else          bkghitcol.emplace_back(BkgClusterHit(999.0, chfcol[ich]));
      }
    }
    event.put(std::move(chcol_out));
    if (savebkg_) {
      event.put(std::make_unique<BkgClusterHitCollection>(bkghitcol));
      event.put(std::make_unique<BkgClusterCollection>(bkgccol));
    }
  }


  //------------------------------------------------------------------------------------------
  void FlagBkgHits2::classifyCluster(BkgClusterCollection& bkgccol, StrawHitFlagCollection& chfcol, const ComboHitCollection& chcol, std::vector<int>& hitToClusterMap, size_t iclStart) const
  {
    float rmsphi(0.0);
    for (size_t icl = iclStart; icl < bkgccol.size(); ++icl) {
      auto& cluster = bkgccol[icl];
      // compute per-plane hit counts
      std::array<int,StrawId::_nplanes> hitplanes{0};
      for (const auto& chit : cluster.hits()) {
        const ComboHit& ch = chcol[chit];
        hitplanes[ch.strawId().plane()] += ch.nStrawHits();
      }
      int ipmin(0), ipmax(StrawId::_nplanes-1);
      while (hitplanes[ipmin] == 0 && ipmin < StrawId::_nplanes) ++ipmin;
      while (hitplanes[ipmax] == 0 && ipmax > 0)                --ipmax;
      unsigned np(0), nhits(0);
      for (int ip = ipmin; ip <= ipmax; ++ip) {
        if (hitplanes[ip] > 0) ++np;
        nhits += hitplanes[ip];
      }
      cluster.setKerasQ(0.1);
      if (nhits >= 5 && np >= 3) {
        float sqrSumDeltaPhi(0.f);
        float phimin =  std::numeric_limits<float>::max();
        float phimax = -std::numeric_limits<float>::max();
        float phiclust = cluster.pos().phi();
        if (phiclust >  M_PI) phiclust -= 2*M_PI;
        if (phiclust < -M_PI) phiclust += 2*M_PI;
        unsigned nchits = cluster.hits().size();
        for (const auto& chit : cluster.hits()) {
          const auto& hit = chcol[chit];
          float dphi_rel = hit.phi() - phiclust;
          if (dphi_rel >  M_PI) dphi_rel -= 2*M_PI;
          if (dphi_rel < -M_PI) dphi_rel += 2*M_PI;
          if (dphi_rel < phimin) phimin = dphi_rel;
          if (dphi_rel > phimax) phimax = dphi_rel;
          sqrSumDeltaPhi += dphi_rel*dphi_rel;
        }
        rmsphi = std::sqrt(sqrSumDeltaPhi / nchits);
        if(rmsphi < 0.12 or np > 20){
          cluster.setKerasQ(1.0);
          StrawHitFlag flag(StrawHitFlag::bkgclust);
          flag.merge(StrawHitFlag(StrawHitFlag::bkg));
          cluster._flag.merge(BkgClusterFlag::bkg);
          for (const auto& chit : cluster.hits()) {
            chfcol[chit].merge(flag);
            hitToClusterMap[chit] = icl;
          }
        }
      }
      //Count protons
      if (!cluster._flag.hasAllProperties(BkgClusterFlag::bkg)) {
        StrawHitFlag bkgFlag(StrawHitFlag::bkg);
        for (const auto& hitIdx : cluster.hits()) {
          if (chcol[hitIdx].energyDep() > minedep_)
            chfcol[hitIdx].merge(bkgFlag);
        }
      }
    }
  }

  //------------------------------------------------------------------------------------------
  void FlagBkgHits2::countProton(BkgClusterCollection& bkgccol, StrawHitFlagCollection& chfcol, const ComboHitCollection& chcol) const
  {
    for (size_t icl = 0; icl < bkgccol.size(); ++icl) {
      auto& cluster = bkgccol[icl];
      if (!cluster._flag.hasAllProperties(BkgClusterFlag::bkg)) {
        StrawHitFlag bkgFlag(StrawHitFlag::bkg);
        for (const auto& hitIdx : cluster.hits()) {
          if (chcol[hitIdx].energyDep() > minedep_)
            chfcol[hitIdx].merge(bkgFlag);
        }
      }
    }
  }

  //------------------------------------------------------------------------------------------
  void FlagBkgHits2::findClusters(BkgClusterCollection& clusters, const ComboHitCollection& chcol,
                                  std::vector<unsigned>& idx, float deltaTime, float deltaZ, float deltaXY2)
  {
    if (idx.empty()) return;

    std::sort(idx.begin(), idx.end(), [&chcol](auto i, auto j){
      return chcol[i].pos().z() > chcol[j].pos().z();
    });

    const unsigned        noiseID(chcol.size()+1u);
    const unsigned        unprocessedID(chcol.size()+2u);
    unsigned              currentClusterID(0);
    std::vector<unsigned> inspect;
    inspect.reserve(idx.size());
    std::vector<unsigned> hitToCluster(idx.size(), unprocessedID);
    std::vector<unsigned> neighbors;
    neighbors.reserve(256);
    size_t prevSize = clusters.size();
    clusters.reserve(prevSize + std::max(16UL, idx.size()/10));

    for (size_t i = 0; i < idx.size(); ++i) {
      if (hitToCluster[i] != unprocessedID) continue;

      int nNeighbors = findNeighbors(i, idx, chcol, neighbors, deltaTime, deltaZ, deltaXY2);
      if (nNeighbors < DBSminExpand_) {
        hitToCluster[i] = noiseID;
        continue;
      }

      hitToCluster[i] = currentClusterID;
      BkgCluster thisCluster;
      thisCluster.addHit(idx[i]);
      inspect.clear();
      for (const auto& j : neighbors) inspect.push_back(j);

      while (!inspect.empty()) {
        auto j = inspect.back();
        inspect.pop_back();

        if (hitToCluster[j] == noiseID) {
          hitToCluster[j] = currentClusterID;
          thisCluster.addHit(idx[j]);
        }
        if (hitToCluster[j] != unprocessedID) continue;

        hitToCluster[j] = currentClusterID;
        thisCluster.addHit(idx[j]);

        nNeighbors = findNeighbors(j, idx, chcol, neighbors, deltaTime, deltaZ, deltaXY2);
        if (nNeighbors >= DBSminExpand_) {
          for (const auto& k : neighbors) {
            if (hitToCluster[k] == unprocessedID || hitToCluster[k] == noiseID)
              inspect.push_back(k);
          }
        }
      }
      if (thisCluster.hits().size() >= minClusterHits_) {
        clusters.push_back(std::move(thisCluster));
        ++currentClusterID;
      }
    }
    for (size_t icl = prevSize; icl < clusters.size(); ++icl) calculateCluster(clusters[icl], chcol);
    if (diag_ > 1) dump(clusters);
  }


  //------------------------------------------------------------------------------------------
  int FlagBkgHits2::findNeighbors(unsigned ihit, const std::vector<unsigned>& idx, const ComboHitCollection& chcol, std::vector<unsigned>& neighbors,
                                  float deltaTime, float deltaZ, float deltaXY2)
  {
    neighbors.clear();
    const auto& hit0 = chcol[idx[ihit]];
    float time0 = hit0.correctedTime();
    float x0    = hit0.pos().x();
    float y0    = hit0.pos().y();
    float z0    = hit0.pos().z();
    unsigned nNeighbors = hit0.nStrawHits()-1;
    // idx is sorted z-descending; binary search for first hit with z <= z0+deltaZ
    float maxZ = z0 + deltaZ;
    auto it_start = std::lower_bound(idx.begin(), idx.end(), maxZ, [&chcol](unsigned i, float val){
      return chcol[i].pos().z() > val;
    });
    size_t istart = std::distance(idx.begin(), it_start);
    for (size_t j = istart; j < idx.size(); ++j) {
      if (j == ihit) continue;
      const auto& hitj = chcol[idx[j]];
      if ((z0 - hitj.pos().z()) > deltaZ) break;
      float dt = hitj.correctedTime() - time0;
      if (std::abs(dt) > deltaTime) continue;
      float dx = hitj.pos().x() - x0;
      float dy = hitj.pos().y() - y0;
      if ((dx*dx + dy*dy) <= deltaXY2) {
        neighbors.emplace_back(j);
        nNeighbors += hitj.nStrawHits();
      }
    }
    return nNeighbors;
  }


  //------------------------------------------------------------------------------------------
  float FlagBkgHits2::distance(const BkgCluster& cluster, const ComboHit& hit) const
  {
    float psep_x = hit.pos().x() - cluster.pos().x();
    float psep_y = hit.pos().y() - cluster.pos().y();
    return sqrtf(psep_x*psep_x + psep_y*psep_y);
  }


  //------------------------------------------------------------------------------------------
  void FlagBkgHits2::calculateCluster(BkgCluster& cluster, const ComboHitCollection& chcol)
  {
    if (cluster.hits().empty()) { cluster.time(0.0f); cluster.pos(XYZVectorF(0.0f,0.0f,0.0f)); return; }
    float sumWeight(0), cx(0), cy(0), ctime(0), cz(0), cedep(0);
    for (auto& idx : cluster.hits()) {
      float weight = chcol[idx].nStrawHits();
      XYZVectorF hitpos(chcol[idx].pos().x(), chcol[idx].pos().y(), chcol[idx].pos().z());
      cluster.addHitPosition(hitpos);
      ctime     += chcol[idx].correctedTime() * weight;
      cx        += chcol[idx].pos().x()       * weight;
      cy        += chcol[idx].pos().y()       * weight;
      cz        += chcol[idx].pos().z()       * weight;
      cedep     += chcol[idx].energyDep()     * weight;
      sumWeight += weight;
    }
    cluster.time(ctime / sumWeight);
    cluster.pos(XYZVectorF(cx/sumWeight, cy/sumWeight, cz/sumWeight));
    cluster.edep(cedep / sumWeight);
  }


  //------------------------------------------------------------------------------------------
  void FlagBkgHits2::dump(const std::vector<BkgCluster>& clusters)
  {
    for (auto& cluster : clusters) {
      for (auto& hit : cluster.hits()) std::cout << hit << " ";
      std::cout << std::endl;
    }
  }

}

DEFINE_ART_MODULE(mu2e::FlagBkgHits2)

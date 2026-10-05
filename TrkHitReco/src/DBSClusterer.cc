#include "Offline/TrkHitReco/inc/DBSClusterer.hh"
#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "Offline/TrkHitReco/inc/TrainBkgDiag.hxx"

#include <algorithm>
#include <cstring>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <limits>


namespace mu2e
{
  DBSClusterer::DBSClusterer(const std::optional<Config> config) :
    DBSminExpand_     (config.value().DBSminN()),
    deltaTime_        (config.value().hitDeltaTime()),
    deltaZ_           (config.value().hitDeltaZ()),
    deltaXY2_         (config.value().hitDeltaXY()*config.value().hitDeltaXY()),
    minClusterHits_   (config.value().minClusterHits()),
    minHitsInCluster_ (config.value().minHitsInCluster()),
    bkgmask_          (config.value().bkgmsk()),
    sigmask_          (config.value().sigmsk()),
    testflag_         (config.value().testflag()),
    kerasW_           (config.value().kerasWeights()),
    normFile_         (config.value().normFile()),
    diag_             (config.value().diag())
  {
  }


  //---------------------------------------------------------------------------------------
  void DBSClusterer::init() {
     ConfigFileLookupPolicy configFile;
     auto kerasWgtsFile = configFile(kerasW_);
     sofiePtr_          = std::make_shared<TMVA_SOFIE_TrainBkgDiag::Session>(kerasWgtsFile);

     std::ifstream fin(configFile(normFile_));
     if (!fin.is_open())
       throw std::runtime_error("DBSClusterer: cannot open normalization file: " + normFile_);
     std::string varName;
     for (int i = 0; i < 7; ++i) {
       if (!(fin >> varName >> normMeans_[i] >> normSigmas_[i]))
         throw std::runtime_error("DBSClusterer: failed reading normalization file at entry " + std::to_string(i));
     }
  }


  //----------------------------------------------------------------------------------------------------------
  void DBSClusterer::findClusters(BkgClusterCollection& clusters, const ComboHitCollection& chcol)
  {  if (chcol.empty()) return;

     std::vector<unsigned> idx; // list of combo hit IDs
     idx.reserve(chcol.size());
     for (size_t ich=0;ich<chcol.size();++ich) {
       if (testflag_ && (!chcol[ich].flag().hasAllProperties(sigmask_) || chcol[ich].flag().hasAnyProperty(bkgmask_))) continue;
       idx.emplace_back(ich);
     }
     // Sort combo hits by correctedTime using a 4-pass radix sort - O(N) vs O(N log N).
     // Adapted from DeltaFinderAlg::orderHits() / countingSortPass() in CalPatRec.
     // floatToSortableInt maps a float to a uint32_t that preserves sort order:
     // negative floats have all bits flipped; positive floats have only the sign bit flipped.
     auto floatToSortableInt = [](float f) -> uint32_t {
       uint32_t u;
       memcpy(&u, &f, sizeof(float));
       if (u & 0x80000000u) u = ~u;
       else                 u |= 0x80000000u;
       return u;
     };
     std::vector<unsigned> sort_tmp(idx.size());
     for (int shift = 0; shift < 32; shift += 8) {
       const bool even = (shift / 8) % 2 == 0;
       std::vector<unsigned>& v1 = ( even) ? idx      : sort_tmp;
       std::vector<unsigned>& v2 = (!even) ? idx      : sort_tmp;
       constexpr int N_BUCKETS = 256;
       uint32_t count[N_BUCKETS] = {};
       const size_t n = v1.size();
       for (size_t i = 0; i < n; i++) {
         const uint8_t byte = (floatToSortableInt(chcol[v1[i]].correctedTime()) >> shift) & 0xFF;
         count[byte]++;
       }
       for (int i = 1; i < N_BUCKETS; i++) count[i] += count[i-1];
       for (size_t i = n - 1; ; i--) {
         const uint8_t byte = (floatToSortableInt(chcol[v1[i]].correctedTime()) >> shift) & 0xFF;
         v2[count[byte] - 1] = v1[i];
         count[byte]--;
         if (i == 0) break; // guard: size_t wraps at 0
       }
     }

     // Build a compact, contiguous cache of the fields findNeighbors() needs.
     // Avoids repeated ComboHit dereferences in the hot loop.
     std::vector<HitData> hitCache;
     hitCache.reserve(idx.size());
     for (unsigned chIdx : idx) {
       const auto& h = chcol[chIdx];
       hitCache.push_back({h.correctedTime(), h.pos().x(), h.pos().y(), h.pos().z(),
                           h.nStrawHits(), chIdx});
     }

     // Precompute lower-bound start indices in O(N) via a two-pointer scan so
     // findNeighbors() no longer pays O(log N) per call.
     std::vector<size_t> lowerBound(hitCache.size(), 0);
     for (size_t i = 0, j = 0; i < hitCache.size(); ++i) {
       float minTime = hitCache[i].time - deltaTime_;
       while (j < hitCache.size() && hitCache[j].time < minTime) ++j;
       lowerBound[i] = j;
     }

     const unsigned        noiseID(chcol.size()+1u);
     const unsigned        unprocessedID(chcol.size()+2u);
     unsigned              currentClusterID(0);
     // Using DFS(stack) instead of BFS(queue)
     // For DBSCAN, cluster membership of core points is invariant,
     // but traversal order may affect assignment of border points
     // in rare ambiguous cases. This change was validated to not
     // impact physics performance while improving cache locality.
     std::vector<unsigned> inspect;
     inspect.reserve(idx.size());
     std::vector<unsigned> hitToCluster(idx.size(),unprocessedID); // number of combohits used for clustering, cluster ID
     std::vector<unsigned> neighbors;
     neighbors.reserve(256);
     clusters.reserve(std::max(16UL, idx.size()/10));
     unsigned nNeighbors = 0;
     for (size_t i=0;i<idx.size();++i) {
       // If a point has already been assigned to a cluster, continue
       if ( hitToCluster[i] != unprocessedID) continue;

       // If the neighborhood is too sparse, assign it to noise
       nNeighbors = findNeighbors(i, lowerBound[i], hitCache, neighbors);
       if (nNeighbors < DBSminExpand_) {
         hitToCluster[i] = noiseID;
         continue;
       }

       hitToCluster[i] = currentClusterID;
       BkgCluster thisCluster;
       thisCluster.addHit(hitCache[i].chIdx);
       // Extend the cluster by adding/expanding around neighbors
       inspect.clear();
       for (const auto& j : neighbors) inspect.push_back(j);

       while (!inspect.empty()){
         auto j = inspect.back();
         inspect.pop_back();

         if (hitToCluster[j] == noiseID) {
           hitToCluster[j] = currentClusterID;
           thisCluster.addHit(hitCache[j].chIdx);
         }
         if (hitToCluster[j] != unprocessedID) continue;

         hitToCluster[j] = currentClusterID;
         thisCluster.addHit(hitCache[j].chIdx);

         nNeighbors = findNeighbors(j, lowerBound[j], hitCache, neighbors);
         if (nNeighbors >= DBSminExpand_){
           for (const auto& k : neighbors) {
             if (hitToCluster[k] == unprocessedID || hitToCluster[k] == noiseID)
               inspect.push_back(k);
           }
         }
       }
       if (thisCluster.hits().size() >= minClusterHits_){
         clusters.push_back(std::move(thisCluster));
         ++currentClusterID;
       }
     }
  }

  //---------------------------------------------------------------------------------------
  // Find the neighbors of given a point - can use any suitable distance function
  unsigned DBSClusterer::findNeighbors(unsigned ihit, size_t istart, const std::vector<HitData>& hitCache, std::vector<unsigned>& neighbors)
  {
    neighbors.clear();
    const HitData& h0 = hitCache[ihit];
    unsigned nNeighbors = (h0.nsh > 0) ? h0.nsh - 1 : 0;
    for (size_t j = istart; j < hitCache.size(); ++j) {
      if (j == ihit) continue;
      const HitData& hj = hitCache[j];
      float dt = hj.time - h0.time;
      if (dt > deltaTime_) break;
      if (std::abs(hj.z - h0.z) > deltaZ_) continue;
      float dx = hj.x - h0.x;
      float dy = hj.y - h0.y;
      if ((dx*dx + dy*dy) <= deltaXY2_) {
        neighbors.emplace_back(j);
        nNeighbors += hj.nsh;
      }
    }
    return nNeighbors;
  }


  //---------------------------------------------------------------------------------------
  // This is only used for diagnosis at this point
  float DBSClusterer::distance(const BkgCluster& cluster, const ComboHit& hit) const
  {
    float psep_x = hit.pos().x()-cluster.pos().x();
    float psep_y = hit.pos().y()-cluster.pos().y();
    return sqrt(psep_x*psep_x+psep_y*psep_y);
  }


  //---------------------------------------------------------------------------------------
  // Compute cluster position/time and run the MVA classifier in a single call.
  // Clusters below minHitsInCluster_ get zeroed defaults and no MVA score.
  // Two loops are unavoidable: the second needs the cluster centre set by the first.
  void DBSClusterer::classifyCluster(BkgCluster& cluster, const ComboHitCollection& chcol)
  {
    if (cluster.hits().size() < minHitsInCluster_) {
      cluster.time(0.0f);
      cluster.pos(XYZVectorF(0.0f, 0.0f, 0.0f));
      cluster.setKerasQ(0.0);
      return;
    }

    // Loop 1: weighted cluster centre + quantities that don't need the centre
    float sumWeight(0), crho(0), ctime(0), cz(0), cedep(0), cphi(0);
    float phi_ref = chcol[cluster.hits().at(0)].phi();
    float zmin = std::numeric_limits<float>::max();
    float zmax = -std::numeric_limits<float>::max();
    unsigned nhits(0);
    unsigned nchits = cluster.hits().size();

    for (auto& hitIdx : cluster.hits()) {
      const auto& hit = chcol[hitIdx];
      float weight = hit.nStrawHits();
      float dt     = hit.correctedTime();
      float dr     = sqrtf(hit.pos().perp2());
      float dz     = hit.pos().z();
      float edep   = hit.energyDep();

      XYZVectorF hitpos = hit.pos();
      cluster.addHitPosition(hitpos);

      float dp   = hitpos.phi();
      float dphi = dp - phi_ref;
      if (dphi > M_PI)  dphi -= 2*M_PI;
      if (dphi < -M_PI) dphi += 2*M_PI;

      ctime     += dt*weight;
      crho      += dr*weight;
      cphi      += (phi_ref + dphi)*weight;
      cz        += dz*weight;
      cedep     += edep*weight;
      sumWeight += weight;

      nhits += hit.nStrawHits();
      if (dz < zmin) zmin = dz;
      if (dz > zmax) zmax = dz;
    }

    cphi  /= sumWeight;
    crho  /= sumWeight;
    ctime /= sumWeight;
    cz    /= sumWeight;
    cedep /= sumWeight;

    if (cphi > M_PI)  cphi -= 2*M_PI;
    if (cphi < -M_PI) cphi += 2*M_PI;

    cluster.time(ctime);
    cluster.pos(XYZVectorF(crho*cos(cphi), crho*sin(cphi), cz));
    cluster.edep(cedep);

    // Loop 2: MVA input variables that require the cluster centre
    double sqrSumDeltaTime(0.), sqrSumDeltaX(0.), sqrSumDeltaY(0.), sqrSumDeltaPhi(0.);
    float phimin = std::numeric_limits<float>::max();
    float phimax = -std::numeric_limits<float>::max();
    float phiclust = cluster.pos().phi();
    if (phiclust > M_PI)  phiclust -= 2*M_PI;
    if (phiclust < -M_PI) phiclust += 2*M_PI;

    for (const auto& chit : cluster.hits()) {
      const auto& hit = chcol[chit];
      float dx = hit.pos().x() - cluster.pos().x();
      float dy = hit.pos().y() - cluster.pos().y();
      float dt = hit.correctedTime() - cluster.time();
      sqrSumDeltaX    += dx*dx;
      sqrSumDeltaY    += dy*dy;
      sqrSumDeltaTime += dt*dt;
      float dphi_rel = hit.phi() - phiclust;
      if (dphi_rel > M_PI)  dphi_rel -= 2*M_PI;
      if (dphi_rel < -M_PI) dphi_rel += 2*M_PI;
      if (dphi_rel < phimin) phimin = dphi_rel;
      if (dphi_rel > phimax) phimax = dphi_rel;
      sqrSumDeltaPhi += dphi_rel*dphi_rel;
    }

    std::array<float,7> kerasvars;
    kerasvars[0] = cluster.pos().Rho();
    kerasvars[1] = zmax - zmin;
    kerasvars[2] = phimax - phimin;
    kerasvars[3] = nhits;
    kerasvars[4] = std::sqrt((sqrSumDeltaX+sqrSumDeltaY)/nchits);
    kerasvars[5] = std::sqrt(sqrSumDeltaTime/nchits);
    kerasvars[6] = std::sqrt(sqrSumDeltaPhi/nchits);
    for (int i = 0; i < 7; ++i)
      kerasvars[i] = (kerasvars[i] - normMeans_[i]) / normSigmas_[i];
    std::vector<float> kerasout = sofiePtr_->infer(kerasvars.data());
    cluster.setKerasQ(kerasout[0]);
  }

}

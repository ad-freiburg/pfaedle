// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_procs() 1
#endif

#include <limits>
#include <map>
#include <set>
#include <stack>
#include <unordered_map>
#include <utility>
#include <vector>

using util::graph::EDijkstra;

// _____________________________________________________________________________
template <typename TW>
std::map<size_t, EdgeListHops> RouterImpl<TW>::route(
    const TripTrie<pfaedle::gtfs::Trip>* trie, const EdgeCandMap& ecm,
    const RoutingOpts& rOpts, const osm::Restrictor& rest, HopCache* hopCache,
    bool noFastHops) const {
  std::map<size_t, EdgeListHops> ret;

  // the current node costs in our DAG
  CostsDAG costsDAG(trie->getNds().size());
  PredeDAG predeDAG(trie->getNds().size());
  std::vector<double> maxCosts(trie->getNds().size());

  // skip the root node, init all to inf
  for (size_t nid = 1; nid < trie->getNds().size(); nid++) {
    costsDAG[nid].resize(ecm.at(nid).size(), DBL_INF);
    predeDAG[nid].resize(ecm.at(nid).size(), NO_PREDE);
  }

  std::stack<size_t> st;

  // init cost of all first childs
  for (size_t cnid : trie->getNd(0).childs) {
    st.push(cnid);
    for (size_t frId = 0; frId < ecm.at(cnid).size(); frId++) {
      costsDAG[cnid][frId] = ecm.at(cnid)[frId].pen;
    }
  }

  while (!st.empty()) {
    size_t frTrNid = st.top();
    st.pop();
    const auto& frTrNd = trie->getNd(frTrNid);

    // determine the max speed for this hop
    double maxSpeed = 0;
    for (size_t nid = 0; nid < ecm.at(frTrNid).size(); nid++) {
      if (!ecm.at(frTrNid)[nid].e) continue;
      if (ecm.at(frTrNid)[nid].e->getFrom()->pl().getComp().maxSpeed > maxSpeed)
        maxSpeed = ecm.at(frTrNid)[nid].e->getFrom()->pl().getComp().maxSpeed;
    }

    for (size_t toTrNid : trie->getNd(frTrNid).childs) {
      CostMatrix costM, dists;
      const auto& toTrNd = trie->getNd(toTrNid);

      if (frTrNd.arr && !toTrNd.arr) {
        for (size_t toId = 0; toId < costsDAG[toTrNid].size(); toId++) {
          auto toCand = ecm.at(toTrNid)[toId];
          for (size_t frId : toCand.depPrede) {
            double newC = costsDAG[frTrNid][frId] + ecm.at(toTrNid)[toId].pen;
            if (newC < costsDAG[toTrNid][toId]) {
              costsDAG[toTrNid][toId] = newC;
              predeDAG[toTrNid][toId] = frId;
            }
          }
        }
        st.push(toTrNid);
        continue;
      }

      const double avgDepT = frTrNd.accTime / frTrNd.trips;
      const double avgArrT = toTrNd.accTime / toTrNd.trips;

      double hopDist = 0;

      hopDist = util::geo::haversine(frTrNd.lat, frTrNd.lng, toTrNd.lat,
                                     toTrNd.lng);

      double minTime = hopDist / maxSpeed;
      double hopTime = avgArrT - avgDepT;

      if (hopTime < minTime) hopTime = minTime;

      uint32_t newMaxCost = TW::maxCost(hopTime, rOpts);
      uint32_t maxCost = newMaxCost;

      bool found = false;
      int step = 0;

      while (!found && step <= MAX_ROUTE_COST_DOUBLING_STEPS) {
        maxCosts[toTrNid] = newMaxCost;
        maxCost = newMaxCost;

        // calculate n x n hops between layers
        if (noFastHops || !TW::ALLOWS_FAST_ROUTE) {
          hops(ecm.at(frTrNid), ecm.at(toTrNid), &costM, &dists, toTrNd.rAttrs,
               rOpts, rest, hopCache, maxCost);
        } else {
          hopsFast(ecm.at(frTrNid), ecm.at(toTrNid), costsDAG[frTrNid], &costM,
                   toTrNd.rAttrs, rOpts, rest, hopCache, maxCost);
        }

        for (size_t matrixI = 0; matrixI < costM.size(); matrixI++) {
          const auto& mVal = costM[matrixI];
          const size_t frId = mVal.first.first;
          const size_t toId = mVal.first.second;
          const uint32_t c = mVal.second;

          double mDist = 0;

          // the dists and the costM matrices have entries at exactly the same
          // loc
          if (TW::NEED_DIST) mDist = dists[matrixI].second;

          // calculate the transition weights
          const double depT = ecm.at(frTrNid)[frId].time;
          const double arrT = ecm.at(toTrNid)[toId].time;
          const double w = TW::weight(c, mDist, arrT - depT, hopDist, rOpts);

          // update costs to successors in next layer
          double newC = costsDAG[frTrNid][frId] + ecm.at(toTrNid)[toId].pen + w;
          if (newC < costsDAG[toTrNid][toId]) {
            costsDAG[toTrNid][toId] = newC;
            predeDAG[toTrNid][toId] = frId;
            found = true;
          }
        }

        if (newMaxCost <= std::numeric_limits<uint32_t>::max() / 2)
          newMaxCost *= 2;
        else
          newMaxCost = std::numeric_limits<uint32_t>::max();

        if (newMaxCost == maxCost) break;
        step++;
      }

      if (!found) {
        // write the cost for the NULL candidates as a fallback
        for (size_t frNid = 0; frNid < ecm.at(frTrNid).size(); frNid++) {
          double newC = costsDAG[frTrNid][frNid] + maxCost * 100;
          // in the time expanded case, there might be multiple null cands
          size_t nullCId = 0;
          while (nullCId < ecm.at(toTrNid).size() &&
                 !ecm.at(toTrNid)[nullCId].e) {
            if (newC < costsDAG[toTrNid][nullCId]) {
              predeDAG[toTrNid][nullCId] = frNid;
              costsDAG[toTrNid][nullCId] = newC;
            }
            nullCId++;
          }
        }

        // for the remaining, write dummy edges
        for (size_t frNid = 0; frNid < ecm.at(frTrNid).size(); frNid++) {
          // skip NULL candidates
          size_t toNid = 1;
          while (toNid < ecm.at(toTrNid).size() && !ecm.at(toTrNid)[toNid].e)
            toNid++;
          for (; toNid < ecm.at(toTrNid).size(); toNid++) {
            double newC = costsDAG[frTrNid][frNid] + ecm.at(toTrNid)[toNid].pen;
            if (newC < costsDAG[toTrNid][toNid]) {
              predeDAG[toTrNid][toNid] = frNid;
              costsDAG[toTrNid][toNid] = newC;
            }
          }
        }
      }

      st.push(toTrNid);
    }
  }

  // update sink costs
  std::unordered_map<size_t, double> sinkCosts;
  std::unordered_map<size_t, size_t> frontIds;
  for (auto leaf : trie->getNdTrips()) {
    sinkCosts[leaf.first] = DBL_INF;
    frontIds[leaf.first] = 0;

    for (size_t lastId = 0; lastId < ecm.at(leaf.first).size(); lastId++) {
      double nCost = costsDAG[leaf.first][lastId];
      if (nCost < sinkCosts[leaf.first]) {
        frontIds[leaf.first] = lastId;
        sinkCosts[leaf.first] = nCost;
      }
    }
  }

  // retrieve edges
  for (auto leaf : trie->getNdTrips()) {
    const auto leafNid = leaf.first;
    auto curTrieNid = leafNid;

    while (predeDAG[curTrieNid][frontIds[leafNid]] != NO_PREDE) {
      const auto curTrieParNid = trie->getNd(curTrieNid).parent;
      const auto frId = predeDAG[curTrieNid][frontIds[leafNid]];
      const auto toId = frontIds[leafNid];

      const auto frTrNd = trie->getNd(curTrieParNid);
      const auto toTrNd = trie->getNd(curTrieNid);

      // skip in-node hops
      if (frTrNd.arr && !toTrNd.arr) {
        frontIds[leafNid] = frId;
        curTrieNid = curTrieParNid;
        continue;
      }

      std::vector<trgraph::Edge*> edgs;

      const auto& fr = ecm.at(curTrieParNid)[frId];
      const auto& to = ecm.at(curTrieNid)[toId];

      // for subtracting and adding progression costs
      typename TW::CostFunc costPr(toTrNd.rAttrs, rOpts, rest, ROUTE_INF);

      if (fr.e && to.e) {
        // account for max progression start offset, do this exactly like
        // in the hops calculation to ensure that we can find the path again
        double maxProgrStart = 0;
        for (const auto& fr : ecm.at(curTrieParNid)) {
          if (!fr.e) continue;
          double progrStart = 0;
          if (fr.progr > 0) progrStart = costPr(fr.e, 0, 0) * fr.progr;
          if (progrStart > maxProgrStart) maxProgrStart = progrStart;
        }

        const double maxCostRt = maxCosts[curTrieNid] + maxProgrStart;
        uint32_t maxCostRtInt = maxCostRt;

        // avoid overflow
        if (maxCostRt >= std::numeric_limits<uint32_t>::max()) {
          maxCostRtInt = std::numeric_limits<uint32_t>::max();
        }

        typename TW::CostFunc cost(toTrNd.rAttrs, rOpts, rest, maxCostRtInt);
        typename TW::DistHeur distH(fr.e->getFrom()->pl().getComp().maxSpeed,
                                    rOpts, {to.e});

        const double c =
            EDijkstra::shortestPath(fr.e, to.e, cost, distH, &edgs);

        if (c < maxCostRtInt) {
          // a path was found, use it
          ret[leafNid].push_back(
              {edgs, fr.e, to.e, fr.progr, to.progr, {}, {}});
        } else {
          // no path was found, which is marked by an empty edge list
          ret[leafNid].push_back({{}, fr.e, to.e, fr.progr, to.progr, {}, {}});
        }
      } else {
        // fallback to the position given in candidate
        if (fr.e) {
          ret[leafNid].push_back({edgs, fr.e, 0, fr.progr, 0, {}, to.point});
        } else if (to.e) {
          ret[leafNid].push_back({edgs, 0, to.e, 0, to.progr, fr.point, {}});
        } else {
          ret[leafNid].push_back({edgs, 0, 0, 0, 0, fr.point, to.point});
        }
      }
      frontIds[leafNid] = frId;
      curTrieNid = curTrieParNid;
    }
  }

  return ret;
}

// _____________________________________________________________________________
template <typename TW>
void RouterImpl<TW>::hops(const EdgeCandGroup& froms, const EdgeCandGroup& tos,
                          CostMatrix* rCosts, CostMatrix* dists,
                          const RoutingAttrs& rAttrs, const RoutingOpts& rOpts,
                          const osm::Restrictor& rest, HopCache* hopCache,
                          uint32_t maxCost) const {
  // standard 1 -> n approach
  std::set<trgraph::Edge*> eFrs;
  for (const auto& from : froms) {
    if (!from.e) continue;
    eFrs.insert(from.e);
  }

  std::set<trgraph::Edge*> eTos;
  for (const auto& to : tos) {
    if (!to.e) continue;
    eTos.insert(to.e);
  }

  EdgeCostMatrix ecm;
  EdgeDistMatrix ecmDist;

  // account for max progression start offset
  double maxProgrStart = 0;
  typename TW::CostFunc cost(rAttrs, rOpts, rest, ROUTE_INF);
  for (const auto& fr : froms) {
    if (!fr.e) continue;
    double progrStart = 0;
    if (fr.progr > 0) progrStart = cost(fr.e, 0, 0) * fr.progr;
    if (progrStart > maxProgrStart) maxProgrStart = progrStart;
  }

  maxCost = addNonOverflow(maxCost, maxProgrStart);
  typename TW::CostFunc costF(rAttrs, rOpts, rest, maxCost);

  for (trgraph::Edge* eFrom : eFrs) {
    std::set<trgraph::Edge*> remTos;
    for (trgraph::Edge* eTo : eTos) {
      // init ecmDist
      ecmDist[eFrom][eTo] = ROUTE_INF;

      std::pair<uint32_t, bool> cached = {0, 0};
      if (hopCache) cached = hopCache->get(eFrom, eTo);

      // shortcut: if the nodes lie in two different connected components,
      // the distance between them is trivially infinite
      if (eFrom->getFrom()->pl().getCompId() !=
          eTo->getTo()->pl().getCompId()) {
        ecm[eFrom][eTo] = costF.inf();
      } else if (cached.second >= costF.inf()) {
        ecm[eFrom][eTo] = costF.inf();
      } else if (!TW::NEED_DIST && cached.second) {
        ecm[eFrom][eTo] = cached.first;
      } else {
        remTos.insert(eTo);
      }
    }

    if (remTos.size()) {
      typename TW::DistHeur distH(eFrom->getFrom()->pl().getComp().maxSpeed,
                                  rOpts, remTos);

      std::unordered_map<trgraph::Edge*, TrEList> paths;
      std::unordered_map<trgraph::Edge*, TrEList*> pathPtrs;
      for (auto to : tos) pathPtrs[to.e] = &paths[to.e];

      const auto& costs =
          EDijkstra::shortestPath(eFrom, remTos, costF, distH, pathPtrs);

      for (const auto& c : costs) {
        ecm[eFrom][c.first] = c.second;

        if (paths[c.first].size() == 0) {
          if (hopCache) hopCache->setMin(eFrom, c.first, maxCost);
          continue;  // no path found
        }

        if (hopCache) hopCache->setEx(eFrom, c.first, c.second);
      }

      if (TW::NEED_DIST) {
        for (const auto& c : costs) {
          if (!paths[c.first].size()) continue;
          double d = 0;
          // don't count last edge
          for (size_t i = paths[c.first].size() - 1; i > 0; i--) {
            d += paths[c.first][i]->pl().getLength();
          }
          ecmDist[eFrom][c.first] = d;
        }
      }
    }
  }

  // build return costs
  for (size_t frId = 0; frId < froms.size(); frId++) {
    auto fr = froms[frId];
    if (!fr.e) continue;
    auto costFr = costF(fr.e, 0, 0);
    for (size_t toId = 0; toId < tos.size(); toId++) {
      auto to = tos[toId];
      if (!to.e) continue;

      uint32_t c = ecm[fr.e][to.e];

      if (c >= maxCost) continue;

      double dist = 0;
      if (TW::NEED_DIST) dist = ecmDist[fr.e][to.e];

      if (fr.e == to.e) {
        if (fr.progr <= to.progr) {
          auto costTo = costF(to.e, 0, 0);
          const uint32_t progrCFr = costFr * fr.progr;
          const uint32_t progrCTo = costTo * to.progr;

          // calculate this in one step to avoid uint32_t underflow below
          c += progrCTo - progrCFr;
        } else {
          // trivial case we can ignore
          continue;
        }

      } else {
        // subtract progression cost on first edge
        if (fr.progr > 0) {
          const uint32_t progrCFr = costFr * fr.progr;
          c -= progrCFr;
          if (TW::NEED_DIST) dist -= fr.e->pl().getLength() * fr.progr;
        }

        // add progression cost on last edge
        if (to.progr > 0) {
          const auto costTo = costF(to.e, 0, 0);
          const uint32_t progrCTo = costTo * to.progr;
          c += progrCTo;
          if (TW::NEED_DIST) dist += to.e->pl().getLength() * to.progr;
        }
      }

      if (c < maxCost - maxProgrStart) {
        rCosts->push_back({{frId, toId}, c});
        if (TW::NEED_DIST)
          dists->push_back({{frId, toId}, static_cast<uint32_t>(dist)});
      }
    }
  }
}

// _____________________________________________________________________________
template <typename TW>
void RouterImpl<TW>::hopsFast(const EdgeCandGroup& froms,
                              const EdgeCandGroup& tos,
                              const LayerCostsDAG& rawInitCosts,
                              CostMatrix* rCosts, const RoutingAttrs& rAttrs,
                              const RoutingOpts& rOpts,
                              const osm::Restrictor& restr, HopCache* hopCache,
                              uint32_t maxCost) const {
  std::unordered_map<trgraph::Edge*, uint32_t> initCosts;

  std::set<trgraph::Edge*> eFrs, eTos;
  std::map<trgraph::Edge*, std::vector<size_t>> eFrCands, eToCands;

  double maxSpeed = 0;
  for (size_t frId = 0; frId < froms.size(); frId++) {
    if (rawInitCosts[frId] >= DBL_INF || !connected(froms[frId], tos)) continue;

    eFrs.insert(froms[frId].e);
    eFrCands[froms[frId].e].push_back(frId);

    if (froms[frId].e->getFrom()->pl().getComp().maxSpeed > maxSpeed)
      maxSpeed = froms[frId].e->getFrom()->pl().getComp().maxSpeed;
  }

  for (size_t toId = 0; toId < tos.size(); toId++) {
    if (!connected(froms, tos[toId]))
      continue;  // skip nodes not conn'ed to any <fr>

    if (hopCache && cacheDrop(hopCache, eFrs, tos[toId].e, maxCost))
      continue;  // skip nodes we have already encountered at higher cost

    eTos.insert(tos[toId].e);
    eToCands[tos[toId].e].push_back(toId);
  }

  if (eFrs.size() == 0 || eTos.size() == 0) return;

  // account for max progression start offset
  double maxProgrStart = 0;
  typename TW::CostFunc progrCostF(rAttrs, rOpts, restr, ROUTE_INF);
  for (const auto& fr : froms) {
    if (!fr.e) continue;
    double progrStart = 0;
    if (fr.progr > 0) progrStart = progrCostF(fr.e, 0, 0) * fr.progr;
    if (progrStart > maxProgrStart) maxProgrStart = progrStart;
  }

  // initialize init doubles
  LayerCostsDAG prepInitCosts(froms.size());
  for (size_t frId = 0; frId < froms.size(); frId++) {
    if (!froms[frId].e || rawInitCosts[frId] >= DBL_INF) continue;
    const auto& fr = froms[frId];
    // offset by progr start
    double progrStart = progrCostF(fr.e, 0, 0) * fr.progr;
    prepInitCosts[frId] =
        TW::invWeight(rawInitCosts[frId], rOpts) + maxProgrStart - progrStart;
  }

  // all init costs are inf
  for (const auto& fr : froms) initCosts[fr.e] = ROUTE_INF;

  // now chose the best offset cost
  for (size_t frId = 0; frId < froms.size(); frId++) {
    if (!froms[frId].e || rawInitCosts[frId] >= DBL_INF) continue;
    const auto& fr = froms[frId];
    if (prepInitCosts[frId] < initCosts[fr.e])
      initCosts[fr.e] = prepInitCosts[frId];
  }

  // get max init costs
  uint32_t maxInit = 0;
  uint32_t minInit = ROUTE_INF;
  for (const auto& c : initCosts) {
    if (!eFrs.count(c.first)) continue;
    if (c.second != ROUTE_INF && c.second > maxInit) maxInit = c.second;
    if (c.second < minInit) minInit = c.second;
  }

  for (auto& c : initCosts) c.second = c.second - minInit;

  // account for start offsets
  maxCost = addNonOverflow(maxCost, maxProgrStart);

  typename TW::CostFunc costF(rAttrs, rOpts, restr,
                              maxCost + (maxInit - minInit));

  std::unordered_map<trgraph::Edge*, TrEList> paths;
  std::unordered_map<trgraph::Edge*, TrEList*> pathPtrs;
  for (const auto& to : tos) pathPtrs[to.e] = &paths[to.e];

  typename TW::DistHeur distH(maxSpeed, rOpts, eTos);

  const auto& costs =
      EDijkstra::shortestPath(eFrs, eTos, initCosts, maxCost, costF, distH);

  for (const auto& c : costs) {
    auto toEdg = c.first;
    if (c.second.second >= costF.inf()) {
      if (hopCache) hopCache->setMin(eFrs, toEdg, maxCost);
      continue;  // no path found
    }
    auto fromEdg = c.second.first;
    uint32_t cost = c.second.second - initCosts[fromEdg];

    if (cost >= maxCost) continue;

    for (size_t frId : eFrCands.find(fromEdg)->second) {
      const auto& fr = froms[frId];
      auto costFr = costF(fr.e, 0, 0);

      for (size_t toId : eToCands.find(toEdg)->second) {
        const auto& to = tos[toId];
        uint32_t wrCost = cost;

        if (fr.e == to.e) {
          if (fr.progr <= to.progr) {
            const auto costTo = costF(to.e, 0, 0);
            const uint32_t progrCFr = costFr * fr.progr;
            const uint32_t progrCTo = costTo * to.progr;

            // calculate this in one step to avoid uint32_t underflow below
            wrCost += progrCTo - progrCFr;
          } else {
            // trivial case we can ignore
            continue;
          }
        } else {
          // subtract progression cost on first edge
          if (fr.progr > 0) {
            const uint32_t progrCFr = costFr * fr.progr;
            wrCost -= progrCFr;
          }

          // add progression cost on last edge
          if (to.progr > 0) {
            const auto costTo = costF(to.e, 0, 0);
            const uint32_t progrCTo = costTo * to.progr;
            wrCost += progrCTo;
          }
        }

        if (wrCost < maxCost - maxProgrStart) {
          rCosts->push_back({{frId, toId}, wrCost});
        }
      }
    }
  }
}

// _____________________________________________________________________________
template <typename TW>
bool RouterImpl<TW>::connected(const EdgeCand& fr,
                               const EdgeCandGroup& tos) const {
  if (!fr.e) return false;
  for (const auto& to : tos) {
    if (!to.e) continue;
    if (fr.e->getFrom()->pl().getCompId() == to.e->getFrom()->pl().getCompId())
      return true;
  }
  return false;
}

// _____________________________________________________________________________
template <typename TW>
bool RouterImpl<TW>::connected(const EdgeCandGroup& froms,
                               const EdgeCand& to) const {
  if (!to.e) return false;
  for (const auto& fr : froms) {
    if (!fr.e) continue;
    if (fr.e->getFrom()->pl().getCompId() == to.e->getFrom()->pl().getCompId())
      return true;
  }
  return false;
}

// _____________________________________________________________________________
template <typename TW>
bool RouterImpl<TW>::cacheDrop(HopCache* hopCache,
                               const std::set<trgraph::Edge*>& froms,
                               const trgraph::Edge* to,
                               uint32_t maxCost) const {
  for (auto fr : froms)
    if (hopCache->get(fr, to).first <= maxCost) return false;

  return true;
}

// _____________________________________________________________________________
template <typename TW>
uint32_t RouterImpl<TW>::addNonOverflow(uint32_t a, uint32_t b) const {
  if (a == std::numeric_limits<uint32_t>::max() ||
      b == std::numeric_limits<uint32_t>::max())
    return std::numeric_limits<uint32_t>::max();
  uint32_t res = a + b;
  if (res >= a && res >= b) return res;
  return std::numeric_limits<uint32_t>::max();
}

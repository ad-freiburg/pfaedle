// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_ROUTER_H_
#define PFAEDLE_ROUTER_ROUTER_H_

#include <limits>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <map>
#include <utility>
#include <vector>
#include "pfaedle/Def.h"
#include "pfaedle/osm/Restrictor.h"
#include "pfaedle/router/HopCache.h"
#include "pfaedle/router/Misc.h"
#include "pfaedle/router/RoutingAttrs.h"
#include "pfaedle/router/TripTrie.h"
#include "pfaedle/router/Weights.h"
#include "pfaedle/trgraph/Graph.h"
#include "util/Misc.h"
#include "util/geo/Geo.h"
#include "util/graph/EDijkstra.h"

namespace pfaedle {
namespace router {

constexpr static uint32_t ROUTE_INF = std::numeric_limits<uint32_t>::max();
constexpr static double DBL_INF = std::numeric_limits<double>::infinity();
constexpr static size_t NO_PREDE = std::numeric_limits<size_t>::max();

constexpr static int MAX_ROUTE_COST_DOUBLING_STEPS = 3;

typedef std::pair<size_t, size_t> HId;
typedef std::vector<double> LayerCostsDAG;
typedef std::vector<LayerCostsDAG> CostsDAG;
typedef std::vector<std::vector<size_t>> PredeDAG;

typedef std::unordered_map<const trgraph::Edge*,
                           std::unordered_map<const trgraph::Edge*, uint32_t>>
    EdgeCostMatrix;
typedef std::unordered_map<const trgraph::Edge*,
                           std::unordered_map<const trgraph::Edge*, double>>
    EdgeDistMatrix;
typedef util::graph::EDijkstra::EList<trgraph::NodePL, trgraph::EdgePL> TrEList;

typedef std::vector<std::pair<std::pair<size_t, size_t>, uint32_t>> CostMatrix;

class Router {
 public:
  virtual ~Router() = default;
  virtual std::map<size_t, EdgeListHops> route(const TripTrie* trie,
                                               const EdgeCandMap& ecm,
                                               const RoutingOpts& rOpts,
                                               const osm::Restrictor& rest,
                                               HopCache* hopCache,
                                               bool noFastHops) const = 0;
};

/*
 * Finds the most likely route of schedule-based vehicle between stops in a
 * physical transportation network
 */
template <typename TW>
class RouterImpl : public Router {
 public:
  // Find the most likely path through the graph for a trip trie.
  virtual std::map<size_t, EdgeListHops> route(
      const TripTrie* trie, const EdgeCandMap& ecm, const RoutingOpts& rOpts,
      const osm::Restrictor& rest, HopCache* hopCache, bool noFastHops) const;

 private:
  void hops(const EdgeCandGroup& from, const EdgeCandGroup& to,
            CostMatrix* rCosts, CostMatrix* dists, const RoutingAttrs& rAttrs,
            const RoutingOpts& rOpts, const osm::Restrictor& rest,
            HopCache* hopCache, uint32_t maxCost) const;

  void hopsFast(const EdgeCandGroup& from, const EdgeCandGroup& to,
                const LayerCostsDAG& initCosts, CostMatrix* rCosts,
                const RoutingAttrs& rAttrs, const RoutingOpts& rOpts,
                const osm::Restrictor& rest,

                HopCache* hopCache, uint32_t maxCost) const;

  bool connected(const EdgeCand& from, const EdgeCandGroup& tos) const;
  bool connected(const EdgeCandGroup& froms, const EdgeCand& to) const;

  bool cacheDrop(

      HopCache* hopCache, const std::set<trgraph::Edge*>& froms,
      const trgraph::Edge* to, uint32_t maxCost) const;

  uint32_t addNonOverflow(uint32_t a, uint32_t b) const;
};

#include "pfaedle/router/Router.tpp"
}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_ROUTER_H_

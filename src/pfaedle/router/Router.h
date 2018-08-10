// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_ROUTER_H_
#define PFAEDLE_ROUTER_ROUTER_H_

#include <mutex>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>
#include <set>
#include <limits>
#include <string>
#include "pfaedle/osm/Restrictor.h"
#include "pfaedle/router/Graph.h"
#include "pfaedle/router/Misc.h"
#include "pfaedle/router/RoutingAttrs.h"
#include "pfaedle/trgraph/Graph.h"
#include "util/graph/Dijkstra.h"
#include "util/graph/EDijkstra.h"
#include "util/geo/Geo.h"
#include "pfaedle/Def.h"

using util::graph::EDijkstra;
using util::graph::Dijkstra;

namespace pfaedle {
namespace router {

typedef std::unordered_map<const trgraph::Edge*, router::Node*> CombNodeMap;
typedef std::pair<size_t, size_t> HId;
typedef std::map<
    RoutingAttrs,
    std::unordered_map<const trgraph::Edge*,
                       std::unordered_map<const trgraph::Edge*,
                                          std::pair<EdgeCost, EdgeList> > > >
    Cache;

struct HopBand {
  double minD;
  double maxD;
  const trgraph::Edge* nearest;
  double maxInGrpDist;
};

struct CostFunc
    : public EDijkstra::CostFunc<trgraph::NodePL, trgraph::EdgePL, EdgeCost> {
  CostFunc(const RoutingAttrs& rAttrs, const RoutingOpts& rOpts,
           const osm::Restrictor& res, const trgraph::StatGroup* tgGrp,
           double max)
      : _rAttrs(rAttrs),
        _rOpts(rOpts),
        _res(res),
        _max(max),
        _tgGrp(tgGrp),
        _inf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, _max, 0) {}

  const RoutingAttrs& _rAttrs;
  const RoutingOpts& _rOpts;
  const osm::Restrictor& _res;
  double _max;
  const trgraph::StatGroup* _tgGrp;
  EdgeCost _inf;

  EdgeCost operator()(const trgraph::Edge* from, const trgraph::Node* n,
                      const trgraph::Edge* to) const;
  EdgeCost inf() const { return _inf; }

  double transitLineCmp(const trgraph::EdgePL& e) const;
};

struct NCostFunc
    : public Dijkstra::CostFunc<trgraph::NodePL, trgraph::EdgePL, EdgeCost> {
  NCostFunc(const RoutingAttrs& rAttrs, const RoutingOpts& rOpts,
            const osm::Restrictor& res, const trgraph::StatGroup* tgGrp)
      : _rAttrs(rAttrs),
        _rOpts(rOpts),
        _res(res),
        _tgGrp(tgGrp),
        _inf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             std::numeric_limits<double>::infinity(), 0) {}

  const RoutingAttrs& _rAttrs;
  const RoutingOpts& _rOpts;
  const osm::Restrictor& _res;
  const trgraph::StatGroup* _tgGrp;
  EdgeCost _inf;

  EdgeCost operator()(const trgraph::Node* from, const trgraph::Edge* e,
                      const trgraph::Node* to) const;
  EdgeCost inf() const { return _inf; }

  double transitLineCmp(const trgraph::EdgePL& e) const;
};

struct DistHeur
    : public EDijkstra::HeurFunc<trgraph::NodePL, trgraph::EdgePL, EdgeCost> {
  DistHeur(uint8_t minLvl, const RoutingOpts& rOpts,
           const std::set<trgraph::Edge*>& tos);

  const RoutingOpts& _rOpts;
  uint8_t _lvl;
  POINT _center;
  double _maxCentD;
  EdgeCost operator()(const trgraph::Edge* a,
                      const std::set<trgraph::Edge*>& b) const;
};

struct NDistHeur
    : public Dijkstra::HeurFunc<trgraph::NodePL, trgraph::EdgePL, EdgeCost> {
  NDistHeur(const RoutingOpts& rOpts, const std::set<trgraph::Node*>& tos);

  const RoutingOpts& _rOpts;
  POINT _center;
  double _maxCentD;
  EdgeCost operator()(const trgraph::Node* a,
                      const std::set<trgraph::Node*>& b) const;
};

struct CombCostFunc
    : public EDijkstra::CostFunc<router::NodePL, router::EdgePL, double> {
  explicit CombCostFunc(const RoutingOpts& rOpts) : _rOpts(rOpts) {}

  const RoutingOpts& _rOpts;

  double operator()(const router::Edge* from, const router::Node* n,
                    const router::Edge* to) const;
  double inf() const { return std::numeric_limits<double>::infinity(); }
};

/*
 * Finds the most likely route of schedule-based vehicle between stops in a
 * physical transportation network
 */
class Router {
 public:
  // Init this router with caches for numThreads threads
  explicit Router(size_t numThreads);
  ~Router();

  // Find the most likely path through the graph for a node candidate route.
  EdgeListHops route(const NodeCandRoute& route, const RoutingAttrs& rAttrs,
                     const RoutingOpts& rOpts, const osm::Restrictor& rest,
                     router::Graph* cgraph) const;

  // Find the most likely path through cgraph for a node candidate route, but
  // based on a greedy node to node approach
  EdgeListHops routeGreedy(const NodeCandRoute& route,
                           const RoutingAttrs& rAttrs, const RoutingOpts& rOpts,
                           const osm::Restrictor& rest) const;

  // Find the most likely path through cgraph for a node candidate route, but
  // based on a greedy node to node set approach
  EdgeListHops routeGreedy2(const NodeCandRoute& route,
                            const RoutingAttrs& rAttrs,
                            const RoutingOpts& rOpts,
                            const osm::Restrictor& rest) const;

  // Return the number of thread caches this router was initialized with
  size_t getCacheNumber() const;

 private:
  mutable std::vector<Cache*> _cache;
  HopBand getHopBand(const NodeCandGroup& a, const NodeCandGroup& b,
                     const RoutingAttrs& rAttrs, const RoutingOpts& rOpts,
                     const osm::Restrictor& rest) const;

  void hops(trgraph::Edge* from, const std::set<trgraph::Edge*>& froms,
            const std::set<trgraph::Edge*> to, const trgraph::StatGroup* tgGrp,
            const std::unordered_map<trgraph::Edge*, EdgeList*>& edgesRet,
            std::unordered_map<trgraph::Edge*, EdgeCost>* rCosts,
            const RoutingAttrs& rAttrs, const RoutingOpts& rOpts,
            const osm::Restrictor& rest, HopBand hopB) const;

  std::set<trgraph::Edge*> getCachedHops(
      trgraph::Edge* from, const std::set<trgraph::Edge*>& to,
      const std::unordered_map<trgraph::Edge*, EdgeList*>& edgesRet,
      std::unordered_map<trgraph::Edge*, EdgeCost>* rCosts,
      const RoutingAttrs& rAttrs) const;

  void cache(trgraph::Edge* from, trgraph::Edge* to, const EdgeCost& c,
             EdgeList* edges, const RoutingAttrs& rAttrs) const;

  void nestedCache(const EdgeList* el, const std::set<trgraph::Edge*>& froms,
                   const CostFunc& cost, const RoutingAttrs& rAttrs) const;

  bool compConned(const NodeCandGroup& a, const NodeCandGroup& b) const;
};
}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_ROUTER_H_

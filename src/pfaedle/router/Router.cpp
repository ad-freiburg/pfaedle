// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_procs() 1
#endif

#include <algorithm>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
#include "pfaedle/router/Comp.h"
#include "pfaedle/router/Router.h"
#include "pfaedle/router/RoutingAttrs.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/graph/Dijkstra.h"
#include "util/graph/EDijkstra.h"
#include "util/log/Log.h"

using pfaedle::router::Router;
using pfaedle::router::EdgeCost;
using pfaedle::router::CostFunc;
using pfaedle::router::DistHeur;
using pfaedle::router::NCostFunc;
using pfaedle::router::NDistHeur;
using pfaedle::router::CombCostFunc;
using pfaedle::router::EdgeListHop;
using pfaedle::router::EdgeListHops;
using pfaedle::router::RoutingOpts;
using pfaedle::router::RoutingAttrs;
using pfaedle::router::HopBand;
using pfaedle::router::NodeCandRoute;
using util::graph::EDijkstra;
using util::graph::Dijkstra;

// _____________________________________________________________________________
EdgeCost NCostFunc::operator()(const trgraph::Node* from,
                               const trgraph::Edge* e,
                               const trgraph::Node* to) const {
  UNUSED(to);
  if (!from) return EdgeCost();

  int oneway = e->pl().oneWay() == 2;
  int32_t stationSkip = 0;

  double transitLinePen = 0;  // transitLineCmp(e->pl());

  return EdgeCost(e->pl().lvl() == 0 ? e->pl().getLength() : 0,
                  e->pl().lvl() == 1 ? e->pl().getLength() : 0,
                  e->pl().lvl() == 2 ? e->pl().getLength() : 0,
                  e->pl().lvl() == 3 ? e->pl().getLength() : 0,
                  e->pl().lvl() == 4 ? e->pl().getLength() : 0,
                  e->pl().lvl() == 5 ? e->pl().getLength() : 0,
                  e->pl().lvl() == 6 ? e->pl().getLength() : 0,
                  e->pl().lvl() == 7 ? e->pl().getLength() : 0, 0, stationSkip,
                  e->pl().getLength() * oneway, oneway,
                  e->pl().getLength() * transitLinePen, 0, &_rOpts);
}

// _____________________________________________________________________________
EdgeCost CostFunc::operator()(const trgraph::Edge* from, const trgraph::Node* n,
                              const trgraph::Edge* to) const {
  if (!from) return EdgeCost();

  uint32_t fullTurns = 0;
  int oneway = from->pl().oneWay() == 2;
  int32_t stationSkip = 0;

  if (n) {
    if (from->getFrom() == to->getTo() && from->getTo() == to->getFrom()) {
      // trivial full turn
      fullTurns = 1;
    } else if (n->getDeg() > 2) {
      // otherwise, only intersection angles will be punished
      fullTurns = router::angSmaller(from->pl().backHop(), *n->pl().getGeom(),
                                     to->pl().frontHop(), _rOpts.fullTurnAngle);
    }

    if (from->pl().isRestricted() && !_res.may(from, to, n)) oneway = 1;

    // for debugging
    n->pl().setVisited();

    if (_tgGrp && n->pl().getSI() && n->pl().getSI()->getGroup() != _tgGrp)
      stationSkip = 1;
  }

  double transitLinePen = transitLineCmp(from->pl());

  return EdgeCost(from->pl().lvl() == 0 ? from->pl().getLength() : 0,
                  from->pl().lvl() == 1 ? from->pl().getLength() : 0,
                  from->pl().lvl() == 2 ? from->pl().getLength() : 0,
                  from->pl().lvl() == 3 ? from->pl().getLength() : 0,
                  from->pl().lvl() == 4 ? from->pl().getLength() : 0,
                  from->pl().lvl() == 5 ? from->pl().getLength() : 0,
                  from->pl().lvl() == 6 ? from->pl().getLength() : 0,
                  from->pl().lvl() == 7 ? from->pl().getLength() : 0, fullTurns,
                  stationSkip, from->pl().getLength() * oneway, oneway,
                  from->pl().getLength() * transitLinePen, 0, &_rOpts);
}

// _____________________________________________________________________________
double CostFunc::transitLineCmp(const trgraph::EdgePL& e) const {
  double best = 1;
  for (const auto* l : e.getLines()) {
    double cur = _rAttrs.simi(l);

    if (cur < 0.0001) return cur;
    if (cur < best) best = cur;
  }

  return best;
}

// _____________________________________________________________________________
NDistHeur::NDistHeur(const RoutingOpts& rOpts,
                     const std::set<trgraph::Node*>& tos)
    : _rOpts(rOpts), _maxCentD(0) {
  size_t c = 0;
  double x = 0, y = 0;
  for (auto to : tos) {
    x += to->pl().getGeom()->getX();
    y += to->pl().getGeom()->getY();
    c++;
  }

  x /= c;
  y /= c;
  _center = FPoint(x, y);

  for (auto to : tos) {
    double cur = static_cast<double>(static_cast<double>(
        util::geo::webMercMeterDist(*to->pl().getGeom(), _center)));
    if (cur > _maxCentD) _maxCentD = cur;
  }
}

// _____________________________________________________________________________
DistHeur::DistHeur(uint8_t minLvl, const RoutingOpts& rOpts,
                   const std::set<trgraph::Edge*>& tos)
    : _rOpts(rOpts), _lvl(minLvl), _maxCentD(0) {
  size_t c = 0;
  double x = 0, y = 0;
  for (auto to : tos) {
    x += to->getFrom()->pl().getGeom()->getX();
    y += to->getFrom()->pl().getGeom()->getY();
    c++;
  }

  x /= c;
  y /= c;
  _center = FPoint(x, y);

  for (auto to : tos) {
    double cur = static_cast<double>(static_cast<double>(
        util::geo::webMercMeterDist(*to->getFrom()->pl().getGeom(), _center) *
        _rOpts.levelPunish[_lvl]));
    if (cur > _maxCentD) _maxCentD = cur;
  }
}

// _____________________________________________________________________________
EdgeCost DistHeur::operator()(const trgraph::Edge* a,
                              const std::set<trgraph::Edge*>& b) const {
  double cur = static_cast<double>(static_cast<double>(
      util::geo::webMercMeterDist(*a->getTo()->pl().getGeom(), _center) *
      _rOpts.levelPunish[_lvl]));

  UNUSED(b);

  return EdgeCost(cur - _maxCentD, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

// _____________________________________________________________________________
EdgeCost NDistHeur::operator()(const trgraph::Node* a,
                               const std::set<trgraph::Node*>& b) const {
  double cur = static_cast<double>(static_cast<double>(
      util::geo::webMercMeterDist(*a->pl().getGeom(), _center)));

  UNUSED(b);

  return EdgeCost(cur - _maxCentD, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
}

// _____________________________________________________________________________
double CombCostFunc::operator()(const router::Edge* from, const router::Node* n,
                                const router::Edge* to) const {
  UNUSED(n);
  UNUSED(from);
  return to->pl().getCost().getValue();
}

// _____________________________________________________________________________
Router::Router(size_t numThreads) : _cache(numThreads) {
  for (size_t i = 0; i < numThreads; i++) {
    _cache[i] = new Cache();
  }
}

// _____________________________________________________________________________
Router::~Router() {
  for (size_t i = 0; i < _cache.size(); i++) {
    delete _cache[i];
  }
}

// _____________________________________________________________________________
bool Router::compConned(const NodeCandGroup& a, const NodeCandGroup& b) const {
  for (auto n1 : a) {
    for (auto n2 : b) {
      if (n1.nd->pl().getComp() == n2.nd->pl().getComp()) return true;
    }
  }

  return false;
}

// _____________________________________________________________________________
HopBand Router::getHopBand(const NodeCandGroup& a, const NodeCandGroup& b,
                           const RoutingAttrs& rAttrs, const RoutingOpts& rOpts,
                           const osm::Restrictor& rest) const {
  double pend = 0;
  for (size_t i = 0; i < a.size(); i++) {
    for (size_t j = 0; j < b.size(); j++) {
      double d = util::geo::webMercMeterDist(*a[i].nd->pl().getGeom(),
                                             *b[j].nd->pl().getGeom());
      if (d > pend) pend = d;
    }
  }

  LOG(VDEBUG) << "Pending max hop distance is " << pend << " meters";

  const trgraph::StatGroup* tgGrpTo = 0;
  if (b.begin()->nd->pl().getSI())
    tgGrpTo = b.begin()->nd->pl().getSI()->getGroup();

  CostFunc costF(rAttrs, rOpts, rest, tgGrpTo, pend * 50);

  std::set<trgraph::Edge *> from, to;

  // TODO(patrick): test if the two sets share a common connected component

  for (auto n : a)
    from.insert(n.nd->getAdjListOut().begin(), n.nd->getAdjListOut().end());

  for (auto n : b)
    to.insert(n.nd->getAdjListOut().begin(), n.nd->getAdjListOut().end());

  LOG(VDEBUG) << "Doing pilot run between " << from.size() << "->" << to.size()
              << " candidates";

  EdgeList el;
  EdgeCost ret = costF.inf();
  DistHeur distH(0, rOpts, to);

  if (compConned(a, b))
    ret = EDijkstra::shortestPath(from, to, costF, distH, &el);

  if (el.size() < 2 && costF.inf() <= ret) {
    LOG(VDEBUG) << "Pilot run: no connection between candidate groups,"
                << " setting max distance to 1";
    return HopBand{0, 1, 0, 0};
  }

  // cache the found path, will save a few dijkstra iterations
  nestedCache(&el, from, costF, rAttrs);

  auto na = el.back()->getFrom();
  auto nb = el.front()->getFrom();

  double maxStrD = 0;

  for (auto e : to) {
    double d = static_cast<double>(static_cast<double>(
        util::geo::webMercMeterDist(*el.front()->getFrom()->pl().getGeom(),
                                    *e->getTo()->pl().getGeom())));
    if (d > maxStrD) maxStrD = d;
  }

  // TODO(patrick): derive the punish level here automatically
  double maxD = std::max(ret.getValue(), pend * rOpts.levelPunish[2]) * 3;
  double minD = ret.getValue();

  LOG(VDEBUG) << "Pilot run: min distance between two groups is "
              << ret.getValue() << " (between nodes " << na << " and " << nb
              << "), using a max routing distance of " << maxD << ". The max"
              << " straight line distance from the pilot target to any other "
                 "target node was"
              << " " << maxStrD << ".";

  return HopBand{minD, maxD, el.front(), maxStrD};
}

// _____________________________________________________________________________
EdgeListHops Router::routeGreedy(const NodeCandRoute& route,
                                 const RoutingAttrs& rAttrs,
                                 const RoutingOpts& rOpts,
                                 const osm::Restrictor& rest) const {
  if (route.size() < 2) return EdgeListHops();
  EdgeListHops ret(route.size() - 1);

  for (size_t i = 0; i < route.size() - 1; i++) {
    const trgraph::StatGroup* tgGrp = 0;
    std::set<trgraph::Node *> from, to;
    for (auto c : route[i]) from.insert(c.nd);
    for (auto c : route[i + 1]) to.insert(c.nd);
    if (route[i + 1].begin()->nd->pl().getSI())
      tgGrp = route[i + 1].begin()->nd->pl().getSI()->getGroup();

    NCostFunc cost(rAttrs, rOpts, rest, tgGrp);
    NDistHeur dist(rOpts, to);

    NodeList nodesRet;
    EdgeListHop hop;
    Dijkstra::shortestPath(from, to, cost, dist, &hop.edges, &nodesRet);

    if (nodesRet.size() > 1) {
      // careful: nodesRet is reversed!
      hop.start = nodesRet.back();
      hop.end = nodesRet.front();
    } else {
      // just take the first candidate if no route could be found
      hop.start = *from.begin();
      hop.end = *to.begin();
    }

    ret[i] = hop;
  }

  return ret;
}

// _____________________________________________________________________________
EdgeListHops Router::routeGreedy2(const NodeCandRoute& route,
                                  const RoutingAttrs& rAttrs,
                                  const RoutingOpts& rOpts,
                                  const osm::Restrictor& rest) const {
  if (route.size() < 2) return EdgeListHops();
  EdgeListHops ret(route.size() - 1);

  for (size_t i = 0; i < route.size() - 1; i++) {
    const trgraph::StatGroup* tgGrp = 0;
    std::set<trgraph::Node *> from, to;

    if (i == 0)
      for (auto c : route[i]) from.insert(c.nd);
    else
      from.insert(const_cast<trgraph::Node*>(ret[i - 1].end));

    for (auto c : route[i + 1]) to.insert(c.nd);

    if (route[i + 1].begin()->nd->pl().getSI())
      tgGrp = route[i + 1].begin()->nd->pl().getSI()->getGroup();

    NCostFunc cost(rAttrs, rOpts, rest, tgGrp);
    NDistHeur dist(rOpts, to);

    NodeList nodesRet;
    EdgeListHop hop;
    Dijkstra::shortestPath(from, to, cost, dist, &hop.edges, &nodesRet);
    if (nodesRet.size() > 1) {
      // careful: nodesRet is reversed!
      hop.start = nodesRet.back();
      hop.end = nodesRet.front();
    } else {
      // just take the first candidate if no route could be found
      hop.start = *from.begin();
      hop.end = *to.begin();
    }

    ret[i] = hop;
  }

  return ret;
}

// _____________________________________________________________________________
EdgeListHops Router::route(const NodeCandRoute& route,
                           const RoutingAttrs& rAttrs, const RoutingOpts& rOpts,
                           const osm::Restrictor& rest,
                           router::Graph* cgraph) const {
  if (route.size() < 2) return EdgeListHops();
  EdgeListHops ret(route.size() - 1);

  CombCostFunc ccost(rOpts);
  router::Node* source = cgraph->addNd();
  router::Node* sink = cgraph->addNd();
  CombNodeMap nodes;
  CombNodeMap nextNodes;

  for (size_t i = 0; i < route[0].size(); i++) {
    for (const auto* e : route[0][i].nd->getAdjListOut()) {
      // we can be sure that each edge is exactly assigned to only one
      // node because the transitgraph is directed
      nodes[e] = cgraph->addNd(route[0][i].nd);
      cgraph->addEdg(source, nodes[e])
          ->pl()
          .setCost(EdgeCost(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                            route[0][i].pen, 0));
    }
  }

  size_t iters = EDijkstra::ITERS;
  double itPerSecTot = 0;
  size_t n = 0;
  for (size_t i = 0; i < route.size() - 1; i++) {
    nextNodes.clear();
    HopBand hopBand = getHopBand(route[i], route[i + 1], rAttrs, rOpts, rest);

    const trgraph::StatGroup* tgGrp = 0;
    if (route[i + 1].begin()->nd->pl().getSI())
      tgGrp = route[i + 1].begin()->nd->pl().getSI()->getGroup();

    std::set<trgraph::Edge*> froms;
    for (const auto& fr : route[i]) {
      froms.insert(fr.nd->getAdjListOut().begin(),
                   fr.nd->getAdjListOut().end());
    }

    for (auto eFr : froms) {
      router::Node* cNodeFr = nodes.find(eFr)->second;

      EdgeSet tos;
      std::map<trgraph::Edge*, router::Edge*> edges;
      std::map<trgraph::Edge*, double> pens;
      std::unordered_map<trgraph::Edge*, EdgeList*> edgeLists;
      std::unordered_map<trgraph::Edge*, EdgeCost> costs;

      assert(route[i + 1].size());

      for (const auto& to : route[i + 1]) {
        assert(to.nd->getAdjListOut().size());
        for (auto eTo : to.nd->getAdjListOut()) {
          tos.insert(eTo);
          if (!nextNodes.count(eTo)) nextNodes[eTo] = cgraph->addNd(to.nd);
          if (i == route.size() - 2) cgraph->addEdg(nextNodes[eTo], sink);

          auto* ce = cgraph->addEdg(cNodeFr, nextNodes[eTo]);
          edges[eTo] = ce;
          pens[eTo] = to.pen;

          edgeLists[eTo] = ce->pl().getEdges();
          ce->pl().setStartNode(eFr->getFrom());
          // for debugging
          ce->pl().setStartEdge(eFr);

          ce->pl().setEndNode(to.nd);
          // for debugging
          ce->pl().setEndEdge(eTo);
        }
      }

      size_t iters = EDijkstra::ITERS;
      auto t1 = TIME();

      assert(tos.size());
      assert(froms.size());

      hops(eFr, froms, tos, tgGrp, edgeLists, &costs, rAttrs, rOpts, rest,
           hopBand);
      double itPerSec =
          (static_cast<double>(EDijkstra::ITERS - iters)) / TOOK(t1, TIME());
      n++;
      itPerSecTot += itPerSec;

      LOG(VDEBUG) << "from " << eFr << ": 1-" << tos.size() << " ("
                  << route[i + 1].size() << " nodes) hop took "
                  << EDijkstra::ITERS - iters << " iterations, "
                  << TOOK(t1, TIME()) << "ms (tput: " << itPerSec << " its/ms)";
      for (auto& kv : edges) {
        kv.second->pl().setCost(
            EdgeCost(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pens[kv.first], 0) +
            costs[kv.first]);

        if (kv.second->pl().getEdges()->size()) {
          if (kv.second->pl().getEdges() &&
              kv.second->pl().getEdges()->size()) {
            // the reach edge is included, but we dont want it in the geometry
            kv.second->pl().getEdges()->erase(
                kv.second->pl().getEdges()->begin());
          }
        }
      }
    }

    std::swap(nodes, nextNodes);
  }

  LOG(VDEBUG) << "Hops took " << EDijkstra::ITERS - iters << " iterations,"
              << " average tput was " << (itPerSecTot / n) << " its/ms";

  iters = EDijkstra::ITERS;
  std::vector<router::Edge*> res;
  EDijkstra::shortestPath(source, sink, ccost, &res);
  size_t j = 0;

  LOG(VDEBUG) << "Optim graph solve took " << EDijkstra::ITERS - iters
              << " iterations.";

  for (auto i = res.rbegin(); i != res.rend(); i++) {
    const auto e = *i;
    if (e->getFrom() != source && e->getTo() != sink) {
      assert(e->pl().frontNode());
      assert(e->pl().backNode());

      ret[j] = EdgeListHop{std::move(*e->pl().getEdges()), e->pl().frontNode(),
                           e->pl().backNode()};
      j++;
    }
  }

  assert(ret.size() == j);
  return ret;
}

// _____________________________________________________________________________
void Router::hops(trgraph::Edge* from, const std::set<trgraph::Edge*>& froms,
                  const std::set<trgraph::Edge*> tos,
                  const trgraph::StatGroup* tgGrp,
                  const std::unordered_map<trgraph::Edge*, EdgeList*>& edgesRet,
                  std::unordered_map<trgraph::Edge*, EdgeCost>* rCosts,
                  const RoutingAttrs& rAttrs, const RoutingOpts& rOpts,
                  const osm::Restrictor& rest, HopBand hopB) const {
  std::set<trgraph::Edge*> rem;

  CostFunc cost(rAttrs, rOpts, rest, tgGrp, hopB.maxD);

  const auto& cached = getCachedHops(from, tos, edgesRet, rCosts, rAttrs);

  for (auto e : cached) {
    // shortcut: if the nodes lie in two different connected components,
    // the distance between them is trivially infinite
    if (e == from || e->getFrom() == from->getFrom() ||
        from->getFrom()->pl().getComp() != e->getTo()->pl().getComp() ||
        e->pl().oneWay() == 2 || from->pl().oneWay() == 2) {
      (*rCosts)[e] = cost.inf();
    } else {
      rem.insert(e);
    }
  }

  LOG(VDEBUG) << "From cache: " << tos.size() - rem.size()
              << ", have to cal: " << rem.size();

  if (rem.size()) {
    DistHeur dist(from->getFrom()->pl().getComp()->minEdgeLvl, rOpts, rem);
    const auto& ret = EDijkstra::shortestPath(from, rem, cost, dist, edgesRet);
    for (const auto& kv : ret) {
      nestedCache(edgesRet.at(kv.first), froms, cost, rAttrs);

      (*rCosts)[kv.first] = kv.second;
    }
  }
}

// _____________________________________________________________________________
void Router::nestedCache(const EdgeList* el,
                         const std::set<trgraph::Edge*>& froms,
                         const CostFunc& cost,
                         const RoutingAttrs& rAttrs) const {
  if (el->size() == 0) return;
  // iterate over result edges backwards
  EdgeList curEdges;
  EdgeCost curCost;

  size_t j = 0;

  for (auto i = el->begin(); i < el->end(); i++) {
    if (curEdges.size()) {
      curCost = curCost + cost(*i, (*i)->getTo(), curEdges.back());
    }

    curEdges.push_back(*i);

    if (froms.count(*i)) {
      EdgeCost startC = cost(0, 0, *i) + curCost;
      cache(*i, el->front(), startC, &curEdges, rAttrs);
      j++;
    }
  }
}

// _____________________________________________________________________________
std::set<pfaedle::trgraph::Edge*> Router::getCachedHops(
    trgraph::Edge* from, const std::set<trgraph::Edge*>& tos,
    const std::unordered_map<trgraph::Edge*, EdgeList*>& edgesRet,
    std::unordered_map<trgraph::Edge*, EdgeCost>* rCosts,
    const RoutingAttrs& rAttrs) const {
  std::set<trgraph::Edge*> ret;
  for (auto to : tos) {
    if ((*_cache[omp_get_thread_num()])[rAttrs][from].count(to)) {
      const auto& cv = (*_cache[omp_get_thread_num()])[rAttrs][from][to];
      (*rCosts)[to] = cv.first;
      *edgesRet.at(to) = cv.second;
    } else {
      ret.insert(to);
    }
  }

  return ret;
}

// _____________________________________________________________________________
void Router::cache(trgraph::Edge* from, trgraph::Edge* to, const EdgeCost& c,
                   EdgeList* edges, const RoutingAttrs& rAttrs) const {
  if (from == to) return;
  (*_cache[omp_get_thread_num()])[rAttrs][from][to] =
      std::pair<EdgeCost, EdgeList>(c, *edges);
}

// _____________________________________________________________________________
size_t Router::getCacheNumber() const { return _cache.size(); }

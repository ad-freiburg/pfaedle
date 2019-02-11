// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_procs() 1
#endif

#include <exception>
#include <map>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "pfaedle/Def.h"
#include "pfaedle/eval/Collector.h"
#include "pfaedle/gtfs/Feed.h"
#include "pfaedle/gtfs/StopTime.h"
#include "pfaedle/osm/OsmBuilder.h"
#include "pfaedle/router/ShapeBuilder.h"
#include "pfaedle/trgraph/StatGroup.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/geo/output/GeoJsonOutput.h"
#include "util/graph/EDijkstra.h"
#include "util/log/Log.h"

using util::geo::extendBox;
using util::geo::DBox;
using util::geo::DPoint;
using util::geo::minbox;

using util::geo::webMercMeterDist;
using util::geo::webMercToLatLng;
using util::geo::latLngToWebMerc;
using util::geo::output::GeoGraphJsonOutput;
using pfaedle::router::ShapeBuilder;
using pfaedle::router::FeedStops;
using pfaedle::router::NodeCandGroup;
using pfaedle::router::NodeCandRoute;
using pfaedle::router::RoutingAttrs;
using pfaedle::router::EdgeListHops;
using pfaedle::router::Clusters;
using pfaedle::osm::BBoxIdx;
using ad::cppgtfs::gtfs::Stop;
using pfaedle::gtfs::Trip;
using pfaedle::gtfs::Feed;
using pfaedle::gtfs::StopTime;
using ad::cppgtfs::gtfs::ShapePoint;

// _____________________________________________________________________________
ShapeBuilder::ShapeBuilder(Feed* feed, ad::cppgtfs::gtfs::Feed* evalFeed,
                           MOTs mots, const config::MotConfig& motCfg,
                           eval::Collector* ecoll, pfaedle::trgraph::Graph* g,
                           router::FeedStops* fStops, osm::Restrictor* restr,
                           const config::Config& cfg)
    : _feed(feed),
      _evalFeed(evalFeed),
      _mots(mots),
      _motCfg(motCfg),
      _ecoll(ecoll),
      _cfg(cfg),
      _g(g),
      _crouter(omp_get_num_procs(), cfg.useCaching),
      _stops(fStops),
      _curShpCnt(0),
      _restr(restr) {
  _numThreads = _crouter.getCacheNumber();
}

// _____________________________________________________________________________
const NodeCandGroup& ShapeBuilder::getNodeCands(const Stop* s) const {
  if (_stops->find(s) == _stops->end() || _stops->at(s) == 0) return _emptyNCG;
  return _stops->at(s)->pl().getSI()->getGroup()->getNodeCands(s);
}

// _____________________________________________________________________________
LINE ShapeBuilder::shapeL(const router::NodeCandRoute& ncr,
                          const router::RoutingAttrs& rAttrs) {
  try {
    const router::EdgeListHops& res = route(ncr, rAttrs);

    LINE l;
    for (const auto& hop : res) {
      const trgraph::Node* last = hop.start;
      if (hop.edges.size() == 0) {
        l.push_back(*hop.start->pl().getGeom());
        l.push_back(*hop.end->pl().getGeom());
      }
      for (auto i = hop.edges.rbegin(); i != hop.edges.rend(); i++) {
        const auto* e = *i;
        if ((e->getFrom() == last) ^ e->pl().isRev()) {
          l.insert(l.end(), e->pl().getGeom()->begin(),
                   e->pl().getGeom()->end());
        } else {
          l.insert(l.end(), e->pl().getGeom()->rbegin(),
                   e->pl().getGeom()->rend());
        }
        last = e->getOtherNd(last);
      }
    }

    return l;
  } catch (const std::runtime_error& e) {
    LOG(ERROR) << e.what();
    return LINE();
  }
}

// _____________________________________________________________________________
LINE ShapeBuilder::shapeL(Trip* trip) {
  return shapeL(getNCR(trip), getRAttrs(trip));
}

// _____________________________________________________________________________
EdgeListHops ShapeBuilder::route(const router::NodeCandRoute& ncr,
                                 const router::RoutingAttrs& rAttrs) const {
  router::Graph g;

  if (_cfg.solveMethod == "global") {
    const router::EdgeListHops& ret =
        _crouter.route(ncr, rAttrs, _motCfg.routingOpts, *_restr, &g);

    // write combination graph
    if (!_cfg.shapeTripId.empty() && _cfg.writeCombGraph) {
      LOG(INFO) << "Outputting combgraph.json...";
      std::ofstream pstr(_cfg.dbgOutputPath + "/combgraph.json");
      GeoGraphJsonOutput o;
      o.print(g, pstr);
    }

    return ret;
  } else if (_cfg.solveMethod == "greedy") {
    return _crouter.routeGreedy(ncr, rAttrs, _motCfg.routingOpts, *_restr);
  } else if (_cfg.solveMethod == "greedy2") {
    return _crouter.routeGreedy2(ncr, rAttrs, _motCfg.routingOpts, *_restr);
  } else {
    LOG(ERROR) << "Unknown solution method " << _cfg.solveMethod;
    exit(1);
  }

  return EdgeListHops();
}

// _____________________________________________________________________________
pfaedle::router::Shape ShapeBuilder::shape(Trip* trip) const {
  LOG(VDEBUG) << "Map-matching shape for trip #" << trip->getId() << " of mot "
              << trip->getRoute()->getType() << "(sn=" << trip->getShortname()
              << ", rsn=" << trip->getRoute()->getShortName()
              << ", rln=" << trip->getRoute()->getLongName() << ")";
  Shape ret;
  ret.hops = route(getNCR(trip), getRAttrs(trip));
  ret.avgHopDist = avgHopDist(trip);

  LOG(VDEBUG) << "Finished map-matching for #" << trip->getId();

  return ret;
}

// _____________________________________________________________________________
pfaedle::router::Shape ShapeBuilder::shape(Trip* trip) {
  LOG(VDEBUG) << "Map-matching shape for trip #" << trip->getId() << " of mot "
              << trip->getRoute()->getType() << "(sn=" << trip->getShortname()
              << ", rsn=" << trip->getRoute()->getShortName()
              << ", rln=" << trip->getRoute()->getLongName() << ")";

  Shape ret;
  ret.hops = route(getNCR(trip), getRAttrs(trip));
  ret.avgHopDist = avgHopDist(trip);

  LOG(VDEBUG) << "Finished map-matching for #" << trip->getId();

  return ret;
}

// _____________________________________________________________________________
void ShapeBuilder::shape(pfaedle::netgraph::Graph* ng) {
  TrGraphEdgs gtfsGraph;

  LOG(DEBUG) << "Clustering trips...";
  Clusters clusters = clusterTrips(_feed, _mots);
  LOG(DEBUG) << "Clustered trips into " << clusters.size() << " clusters.";

  std::map<std::string, size_t> shpUsage;
  for (auto t : _feed->getTrips()) {
    if (!t.getShape().empty()) shpUsage[t.getShape()]++;
  }

  // to avoid unfair load balance on threads
  std::random_shuffle(clusters.begin(), clusters.end());

  size_t iters = EDijkstra::ITERS;
  size_t totiters = EDijkstra::ITERS;
  size_t oiters = EDijkstra::ITERS;
  size_t j = 0;

  auto t1 = TIME();
  auto t2 = TIME();
  double totAvgDist = 0;
  size_t totNumTrips = 0;

#pragma omp parallel for num_threads(_numThreads)
  for (size_t i = 0; i < clusters.size(); i++) {
    j++;

    if (j % 10 == 0) {
#pragma omp critical
      {
        LOG(INFO) << "@ " << j << " / " << clusters.size() << " ("
                  << (static_cast<int>((j * 1.0) / clusters.size() * 100))
                  << "%, " << (EDijkstra::ITERS - oiters) << " iters, "
                  << "matching " << (10.0 / (TOOK(t1, TIME()) / 1000))
                  << " trips/sec)";

        oiters = EDijkstra::ITERS;
        t1 = TIME();
      }
    }

    // explicitly call const version of shape here for thread safety
    const Shape& cshp =
        const_cast<const ShapeBuilder&>(*this).shape(clusters[i][0]);
    totAvgDist += cshp.avgHopDist;

    if (_cfg.buildTransitGraph) {
#pragma omp critical
      { writeTransitGraph(cshp, &gtfsGraph, clusters[i]); }
    }

    std::vector<double> distances;
    const ad::cppgtfs::gtfs::Shape& shp =
        getGtfsShape(cshp, clusters[i][0], &distances);

    LOG(VDEBUG) << "Took " << EDijkstra::ITERS - iters << " iterations.";
    iters = EDijkstra::ITERS;

    totNumTrips += clusters[i].size();

    for (auto t : clusters[i]) {
      if (_cfg.evaluate && _evalFeed && _ecoll) {
        _ecoll->add(t, _evalFeed->getShapes().get(t->getShape()), shp,
                    distances);
      }

      if (!t->getShape().empty() && shpUsage[t->getShape()] > 0) {
        shpUsage[t->getShape()]--;
        if (shpUsage[t->getShape()] == 0) {
          _feed->getShapes().remove(t->getShape());
        }
      }
      setShape(t, shp, distances);
    }
  }

  LOG(INFO) << "Matched " << totNumTrips << " trips in " << clusters.size()
            << " clusters.";
  LOG(DEBUG) << "Took " << (EDijkstra::ITERS - totiters)
             << " iterations in total.";
  LOG(DEBUG) << "Took " << TOOK(t2, TIME()) << " ms in total.";
  LOG(DEBUG) << "Total avg. tput "
             << (static_cast<double>(EDijkstra::ITERS - totiters)) /
                    TOOK(t2, TIME())
             << " iters/sec";
  LOG(DEBUG) << "Total avg. trip tput "
             << (clusters.size() / (TOOK(t2, TIME()) / 1000)) << " trips/sec";
  LOG(DEBUG) << "Avg hop distance was "
             << (totAvgDist / static_cast<double>(clusters.size()))
             << " meters";

  if (_cfg.buildTransitGraph) {
    LOG(INFO) << "Building transit network graph...";
    buildTrGraph(&gtfsGraph, ng);
  }
}

// _____________________________________________________________________________
void ShapeBuilder::setShape(Trip* t, const ad::cppgtfs::gtfs::Shape& s,
                            const std::vector<double>& distances) {
  assert(distances.size() == t->getStopTimes().size());
  // set distances
  size_t i = 0;
  for (const auto& st : t->getStopTimes()) {
    const_cast<StopTime<Stop>&>(st).setShapeDistanceTravelled(distances[i]);
    i++;
  }

  std::lock_guard<std::mutex> guard(_shpMutex);
  // TODO(patrick):
  t->setShape(_feed->getShapes().add(s));
}

// _____________________________________________________________________________
ad::cppgtfs::gtfs::Shape ShapeBuilder::getGtfsShape(
    const Shape& shp, Trip* t, std::vector<double>* hopDists) {
  ad::cppgtfs::gtfs::Shape ret(getFreeShapeId(t));

  assert(shp.hops.size() == t->getStopTimes().size() - 1);

  size_t seq = 0;
  double dist = -1;
  double lastDist = -1;
  hopDists->push_back(0);
  POINT last(0, 0);
  for (const auto& hop : shp.hops) {
    const trgraph::Node* l = hop.start;
    if (hop.edges.size() == 0) {
      POINT ll = webMercToLatLng<PFAEDLE_PRECISION>(
          hop.start->pl().getGeom()->getX(), hop.start->pl().getGeom()->getY());

      if (dist > -0.5)
        dist += webMercMeterDist(last, *hop.start->pl().getGeom());
      else
        dist = 0;

      last = *hop.start->pl().getGeom();

      if (dist - lastDist > 0.01) {
        ret.addPoint(ShapePoint(ll.getY(), ll.getX(), dist, seq));
        seq++;
        lastDist = dist;
      }

      dist += webMercMeterDist(last, *hop.end->pl().getGeom());
      last = *hop.end->pl().getGeom();

      if (dist - lastDist > 0.01) {
        ll = webMercToLatLng<PFAEDLE_PRECISION>(
            hop.end->pl().getGeom()->getX(), hop.end->pl().getGeom()->getY());
        ret.addPoint(ShapePoint(ll.getY(), ll.getX(), dist, seq));
        seq++;
        lastDist = dist;
      }
    }

    for (auto i = hop.edges.rbegin(); i != hop.edges.rend(); i++) {
      const auto* e = *i;
      if ((e->getFrom() == l) ^ e->pl().isRev()) {
        for (size_t i = 0; i < e->pl().getGeom()->size(); i++) {
          const POINT& cur = (*e->pl().getGeom())[i];
          if (dist > -0.5)
            dist += webMercMeterDist(last, cur);
          else
            dist = 0;
          last = cur;
          if (dist - lastDist > 0.01) {
            POINT ll =
                webMercToLatLng<PFAEDLE_PRECISION>(cur.getX(), cur.getY());
            ret.addPoint(ShapePoint(ll.getY(), ll.getX(), dist, seq));
            seq++;
            lastDist = dist;
          }
        }
      } else {
        for (int64_t i = e->pl().getGeom()->size() - 1; i >= 0; i--) {
          const POINT& cur = (*e->pl().getGeom())[i];
          if (dist > -0.5)
            dist += webMercMeterDist(last, cur);
          else
            dist = 0;
          last = cur;
          if (dist - lastDist > 0.01) {
            POINT ll =
                webMercToLatLng<PFAEDLE_PRECISION>(cur.getX(), cur.getY());
            ret.addPoint(ShapePoint(ll.getY(), ll.getX(), dist, seq));
            seq++;
            lastDist = dist;
          }
        }
      }
      l = e->getOtherNd(l);
    }

    hopDists->push_back(lastDist);
  }

  return ret;
}

// _____________________________________________________________________________
std::string ShapeBuilder::getFreeShapeId(Trip* trip) {
  std::string ret;
  std::lock_guard<std::mutex> guard(_shpMutex);
  while (!ret.size() || _feed->getShapes().get(ret)) {
    _curShpCnt++;
    ret = "shp_";
    ret += std::to_string(trip->getRoute()->getType());
    ret += "_" + std::to_string(_curShpCnt);
  }

  return ret;
}

// _____________________________________________________________________________
const RoutingAttrs& ShapeBuilder::getRAttrs(const Trip* trip) {
  auto i = _rAttrs.find(trip);

  if (i == _rAttrs.end()) {
    router::RoutingAttrs ret;

    const auto& lnormzer = _motCfg.osmBuildOpts.lineNormzer;

    ret.shortName = lnormzer(trip->getRoute()->getShortName());

    if (ret.shortName.empty()) ret.shortName = lnormzer(trip->getShortname());

    if (ret.shortName.empty())
      ret.shortName = lnormzer(trip->getRoute()->getLongName());

    ret.fromString = _motCfg.osmBuildOpts.statNormzer(
        trip->getStopTimes().begin()->getStop()->getName());
    ret.toString = _motCfg.osmBuildOpts.statNormzer(
        (--trip->getStopTimes().end())->getStop()->getName());

    return _rAttrs
        .insert(std::pair<const Trip*, router::RoutingAttrs>(trip, ret))
        .first->second;
  } else {
    return i->second;
  }
}

// _____________________________________________________________________________
const RoutingAttrs& ShapeBuilder::getRAttrs(const Trip* trip) const {
  return _rAttrs.find(trip)->second;
}

// _____________________________________________________________________________
void ShapeBuilder::getGtfsBox(const Feed* feed, const MOTs& mots,
                              const std::string& tid, bool dropShapes,
                              osm::BBoxIdx* box) {
  for (const auto& t : feed->getTrips()) {
    if (!tid.empty() && t.getId() != tid) continue;
    if (tid.empty() && !t.getShape().empty() && !dropShapes) continue;
    if (t.getStopTimes().size() < 2) continue;

    if (mots.count(t.getRoute()->getType())) {
      DBox cur;
      for (const auto& st : t.getStopTimes()) {
        cur = extendBox(DPoint(st.getStop()->getLng(), st.getStop()->getLat()),
                        cur);
      }
      box->add(cur);
    }
  }
}

// _____________________________________________________________________________
NodeCandRoute ShapeBuilder::getNCR(Trip* trip) const {
  router::NodeCandRoute ncr(trip->getStopTimes().size());

  size_t i = 0;

  for (const auto& st : trip->getStopTimes()) {
    ncr[i] = getNodeCands(st.getStop());
    if (ncr[i].size() == 0) {
      throw std::runtime_error("No node candidate found for station '" +
                               st.getStop()->getName() + "' on trip '" +
                               trip->getId() + "'");
    }
    i++;
  }
  return ncr;
}

// _____________________________________________________________________________
double ShapeBuilder::avgHopDist(Trip* trip) const {
  size_t i = 0;
  double sum = 0;

  const Stop* prev = 0;

  for (const auto& st : trip->getStopTimes()) {
    if (!prev) {
      prev = st.getStop();
      continue;
    }
    auto a = util::geo::latLngToWebMerc<PFAEDLE_PRECISION>(prev->getLat(),
                                                           prev->getLng());
    auto b = util::geo::latLngToWebMerc<PFAEDLE_PRECISION>(
        st.getStop()->getLat(), st.getStop()->getLng());
    sum += util::geo::webMercMeterDist(a, b);

    prev = st.getStop();
    i++;
  }
  return sum / static_cast<double>(i);
}

// _____________________________________________________________________________
Clusters ShapeBuilder::clusterTrips(Feed* f, MOTs mots) {
  // building an index [start station, end station] -> [cluster]

  std::map<StopPair, std::vector<size_t>> clusterIdx;

  Clusters ret;
  for (auto& trip : f->getTrips()) {
    if (!trip.getShape().empty() && !_cfg.dropShapes) continue;
    if (trip.getStopTimes().size() < 2) continue;
    if (!mots.count(trip.getRoute()->getType()) ||
        !_motCfg.mots.count(trip.getRoute()->getType()))
      continue;

    bool found = false;
    auto spair = StopPair(trip.getStopTimes().begin()->getStop(),
                          trip.getStopTimes().rbegin()->getStop());
    const auto& c = clusterIdx[spair];

    for (size_t i = 0; i < c.size(); i++) {
      if (routingEqual(ret[c[i]][0], &trip)) {
        ret[c[i]].push_back(&trip);
        found = true;
        break;
      }
    }
    if (!found) {
      ret.push_back(Cluster{&trip});
      // explicit call to write render attrs to cache
      getRAttrs(&trip);
      clusterIdx[spair].push_back(ret.size() - 1);
    }
  }

  return ret;
}

// _____________________________________________________________________________
bool ShapeBuilder::routingEqual(const Stop* a, const Stop* b) {
  if (a == b) return true;  // trivial

  auto namea = _motCfg.osmBuildOpts.statNormzer(a->getName());
  auto nameb = _motCfg.osmBuildOpts.statNormzer(b->getName());
  if (namea != nameb) return false;

  auto tracka = _motCfg.osmBuildOpts.trackNormzer(a->getPlatformCode());
  auto trackb = _motCfg.osmBuildOpts.trackNormzer(b->getPlatformCode());
  if (tracka != trackb) return false;

  POINT ap =
      util::geo::latLngToWebMerc<PFAEDLE_PRECISION>(a->getLat(), a->getLng());
  POINT bp =
      util::geo::latLngToWebMerc<PFAEDLE_PRECISION>(b->getLat(), b->getLng());

  double d = util::geo::webMercMeterDist(ap, bp);

  if (d > 1) return false;

  return true;
}

// _____________________________________________________________________________
bool ShapeBuilder::routingEqual(Trip* a, Trip* b) {
  if (a->getStopTimes().size() != b->getStopTimes().size()) return false;
  if (getRAttrs(a) != getRAttrs(b)) return false;

  auto stb = b->getStopTimes().begin();
  for (const auto& sta : a->getStopTimes()) {
    if (!routingEqual(sta.getStop(), stb->getStop())) {
      return false;
    }
    stb++;
  }

  return true;
}

// _____________________________________________________________________________
const pfaedle::trgraph::Graph* ShapeBuilder::getGraph() const { return _g; }

// _____________________________________________________________________________
void ShapeBuilder::writeTransitGraph(const Shape& shp, TrGraphEdgs* edgs,
                                     const Cluster& cluster) const {
  for (auto hop : shp.hops) {
    for (const auto* e : hop.edges) {
      if (e->pl().isRev()) e = _g->getEdg(e->getTo(), e->getFrom());
      (*edgs)[e].insert(cluster.begin(), cluster.end());
    }
  }
}

// _____________________________________________________________________________
void ShapeBuilder::buildTrGraph(TrGraphEdgs* edgs,
                                pfaedle::netgraph::Graph* ng) const {
  std::unordered_map<trgraph::Node*, pfaedle::netgraph::Node*> nodes;

  for (auto ep : *edgs) {
    auto e = ep.first;
    pfaedle::netgraph::Node* from = 0;
    pfaedle::netgraph::Node* to = 0;
    if (nodes.count(e->getFrom())) from = nodes[e->getFrom()];
    if (nodes.count(e->getTo())) to = nodes[e->getTo()];
    if (!from) {
      from = ng->addNd(*e->getFrom()->pl().getGeom());
      nodes[e->getFrom()] = from;
    }
    if (!to) {
      to = ng->addNd(*e->getTo()->pl().getGeom());
      nodes[e->getTo()] = to;
    }

    ng->addEdg(from, to,
               pfaedle::netgraph::EdgePL(*e->pl().getGeom(), ep.second));
  }
}

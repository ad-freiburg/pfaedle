// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <atomic>
#include <cstdlib>
#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <stdexcept>
#include <thread>
#include <unordered_map>
#include <utility>

#include "ad/cppgtfs/gtfs/Feed.h"
#include "pfaedle/Def.h"
#include "pfaedle/gtfs/Feed.h"
#include "pfaedle/gtfs/StopTime.h"
#include "pfaedle/osm/OsmBuilder.h"
#include "pfaedle/router/ShapeBuilder.h"
#include "pfaedle/statsimi-classifier/StatsimiClassifier.h"
#include "util/geo/Geo.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/geo/output/GeoJsonOutput.h"
#include "util/graph/EDijkstra.h"
#include "util/log/Log.h"

using util::geo::DBox;
using util::geo::DPoint;
using util::geo::extendBox;
using util::geo::minbox;
using util::geo::PolyLine;

using ad::cppgtfs::gtfs::NO_COLOR;
using ad::cppgtfs::gtfs::ShapePoint;
using ad::cppgtfs::gtfs::Stop;
using pfaedle::gtfs::Feed;
using pfaedle::gtfs::StopTime;
using pfaedle::gtfs::Trip;
using pfaedle::osm::BBoxIdx;
using pfaedle::router::EdgeCandGroup;
using pfaedle::router::EdgeCandMap;
using pfaedle::router::EdgeListHops;
using pfaedle::router::FeedStops;
using pfaedle::router::RoutingAttrs;
using pfaedle::router::ShapeBuilder;
using pfaedle::router::Stats;
using pfaedle::router::TripForests;
using pfaedle::router::TripTrie;
using pfaedle::trgraph::EdgeGrid;
using pfaedle::trgraph::NodeGrid;
using util::geo::latLngToWebMerc;
using util::geo::M_PER_DEG;
using util::geo::output::GeoGraphJsonOutput;

// _____________________________________________________________________________
ShapeBuilder::ShapeBuilder(
    Feed* feed, MOTs mots, const config::MotConfig& motCfg,
    pfaedle::trgraph::Graph* g, router::FeedStops* fStops,
    osm::Restrictor* restr,
    const pfaedle::statsimiclassifier::StatsimiClassifier* classifier,
    router::Router* router, const config::Config& cfg)
    : _feed(feed),
      _mots(mots),
      _motCfg(motCfg),
      _cfg(cfg),
      _g(g),
      _stops(fStops),
      _curShpCnt(0),
      _restr(restr),
      _classifier(classifier),
      _router(router) {
  pfaedle::osm::BBoxIdx box(cfg.boxPadding);
  ShapeBuilder::getGtfsBox(feed, mots, cfg.shapeTripId, cfg.dropShapes, &box,
                           _motCfg.osmBuildOpts.maxSpeed, 0, cfg.verbosity);

  _eGrid = EdgeGrid(cfg.gridSize, cfg.gridSize, box.getFullBox(), false);
  _nGrid = NodeGrid(cfg.gridSize, cfg.gridSize, box.getFullBox(), false);

  LOG(DEBUG) << "Grid size of " << _nGrid.getXWidth() << "x"
             << _nGrid.getYHeight();

  buildIndex();
}

// _____________________________________________________________________________
void ShapeBuilder::buildIndex() {
  for (auto* n : _g->getNds()) {
    for (auto* e : n->getAdjListOut()) {
      if (e->pl().lvl() > _motCfg.osmBuildOpts.maxSnapLevel) continue;
      // don't snap to one way edges
      if (e->pl().oneWay() == 2) continue;

      _eGrid.add(*e->pl().getGeom(), e);
    }
  }

  for (auto* n : _g->getNds()) {
    // only station nodes
    if (n->pl().getSI()) {
      _nGrid.add(*n->pl().getGeom(), n);
    }
  }
}

// _____________________________________________________________________________
void ShapeBuilder::buildCandCache(const TripForests& forests) {
  std::set<const Stop*> stops;
  size_t count = 0;

  for (const auto& forest : forests) {
    for (const auto& trie : forest.second) {
      for (const auto& trips : trie.getNdTrips()) {
        for (const auto& st : trips.second[0]->getStopTimes()) {
          stops.insert(st.getStop());
        }
      }
    }
  }

  size_t numThreads = std::thread::hardware_concurrency();
  std::vector<std::thread> thrds(numThreads);
  std::vector<GrpCache> caches(numThreads);
  std::vector<std::vector<const Stop*>> threadStops(numThreads);

  size_t i = 0;
  for (auto stop : stops) {
    threadStops[i].push_back(stop);
    if (++i == numThreads) i = 0;
  }

  i = 0;
  for (auto& t : thrds) {
    t = std::thread(&ShapeBuilder::edgCandWorker, this, &threadStops[i],
                    &caches[i]);
    i++;
  }

  for (auto& thr : thrds) thr.join();

  // merge
  for (size_t i = 0; i < numThreads; i++) {
    for (const auto& c : caches[i]) {
      _grpCache[c.first] = c.second;
      count += c.second.size();
    }
  }

  if (_grpCache.size())
    LOG(DEBUG) << "Average candidate set size: "
               << ((count * 1.0) / _grpCache.size());
}

// _____________________________________________________________________________
EdgeCandGroup ShapeBuilder::getEdgCands(const Stop* s) const {
  auto cached = _grpCache.find(s);
  if (cached != _grpCache.end()) return cached->second;

  EdgeCandGroup ret;

  const auto& snormzer = _motCfg.osmBuildOpts.statNormzer;
  auto normedName = snormzer.norm(s->getName());

  // the first cand is a placeholder for the stop position itself, it is chosen
  // when no candidate yielded a feasible route
  auto pos = POINT(s->getLng(), s->getLat());
  ret.push_back({0, 0, 0, pos, 0, {}});

  double maxMDist = _motCfg.osmBuildOpts.maxStationCandDistance;

  double distor = util::geo::latLngDistFactor(pos);

  if (_cfg.gaussianNoise > 0) {
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine gen(seed);

    // the standard dev is given in meters, convert (roughly...) to degrees
    double standardDev = (_cfg.gaussianNoise / M_PER_DEG) / distor;

    // mean 0 (no movement), standard dev according to config
    std::normal_distribution<double> dist(0.0, standardDev);

    // add gaussian noise
    pos.setX(pos.getX() + dist(gen));
    pos.setY(pos.getY() + dist(gen));
  }

  std::set<trgraph::Node*> frNIdx;
  _nGrid.get(util::geo::pad(util::geo::getBoundingBox(pos),
                            (maxMDist / M_PER_DEG) / distor),
             &frNIdx);

  if (_motCfg.routingOpts.useStations) {
    for (auto nd : frNIdx) {
      assert(nd->pl().getSI());

      double mDist = util::geo::haversine(pos, *nd->pl().getGeom());
      if (mDist > maxMDist) continue;

      double nameMatchPunish = 0;
      double trackMatchPunish = 0;

      if (!_classifier->similar(normedName, pos, nd->pl().getSI()->getName(),
                                *nd->pl().getGeom())) {
        // stations do not match, punish
        nameMatchPunish = _motCfg.routingOpts.stationUnmatchedPen;
      }
      std::string platform = s->getPlatformCode();

      if (!platform.empty() && !nd->pl().getSI()->getTrack().empty() &&
          nd->pl().getSI()->getTrack() == platform) {
        trackMatchPunish = _motCfg.routingOpts.platformUnmatchedPen;
      }

      for (auto* e : nd->getAdjListOut()) {
        // don't snap to one way edges
        if (e->pl().oneWay() == 2) continue;
        ret.push_back({e,
                       emWeight(mDist) + nameMatchPunish + trackMatchPunish,
                       0,
                       {},
                       0,
                       {}});
      }
    }
  }

  maxMDist = _motCfg.osmBuildOpts.maxSnapDistance;

  std::set<trgraph::Edge*> frEIdx;
  _eGrid.get(util::geo::pad(util::geo::getBoundingBox(pos),
                            (maxMDist / M_PER_DEG) / distor),
             &frEIdx);

  std::set<trgraph::Edge*> selected;
  std::map<const trgraph::Edge*, double> scores;
  std::map<const trgraph::Edge*, double> progrs;

  for (auto edg : frEIdx) {
    if (selected.count(edg)) continue;

    auto reach = deg2reachable(edg, selected);

    double mDist = dist(pos, *edg->pl().getGeom()) * distor * M_PER_DEG;

    if (mDist > maxMDist) continue;

    if (!reach || mDist < scores[reach]) {
      if (reach) {
        selected.erase(selected.find(reach));
        scores.erase(scores.find(reach));
      }
      util::geo::PolyLine<double> pl(*edg->pl().getGeom());
      auto lp = pl.projectOn(pos);
      double progr = lp.totalPos;
      if (edg->pl().isRev()) progr = 1 - progr;
      selected.insert(edg);
      scores[edg] = mDist;
      progrs[edg] = progr;
    }
  }

  for (auto e : selected) {
    ret.push_back({e,
                   emWeight(scores[e]) + _motCfg.routingOpts.nonStationPen,
                   progrs[e],
                   {},
                   0,
                   {}});
  }

  if (ret.size() == 1 && _cfg.verbosity) {
    LOG(WARN) << "No snapping candidate found for stop '" << s->getName()
              << "' (" << s->getId() << ")";
  }

  return ret;
}

// _____________________________________________________________________________
pfaedle::trgraph::Edge* ShapeBuilder::deg2reachable(
    trgraph::Edge* e, std::set<trgraph::Edge*> edgs) const {
  trgraph::Edge* cur = e;

  // forward
  while (cur->getTo()->getDeg() == 2) {
    // dont allow backtracking on reverse edge
    auto next = e->getTo()->getAdjListOut().front()->getTo() == e->getFrom()
                    ? e->getTo()->getAdjListOut().back()
                    : e->getTo()->getAdjListOut().front();
    if (next == e || next == cur) break;  // avoid circles
    if (next->pl().oneWay() == 2) break;  // dont follow one way edges
    if (edgs.count(next)) return next;
    cur = next;
  }

  // backward
  while (cur->getFrom()->getDeg() == 2) {
    // dont allow backtracking on reverse edge
    auto next = e->getFrom()->getAdjListIn().front()->getFrom() == e->getTo()
                    ? e->getFrom()->getAdjListIn().back()
                    : e->getFrom()->getAdjListIn().front();
    if (next == e || next == cur) break;  // avoid circles
    if (next->pl().oneWay() == 2) break;  // dont follow one way edges
    if (edgs.count(cur)) return next;
    cur = next;
  }

  return 0;
}

// _____________________________________________________________________________
std::pair<std::vector<LINE>, Stats> ShapeBuilder::shapeL(Trip* trip) {
  Stats stats;
  try {
    T_START(t);
    EDijkstra::ITERS = 0;
    auto hops = shapeify(trip);
    stats.solveTime = T_STOP(t);
    stats.numTries = 1;
    stats.numTrieLeafs = 1;
    stats.totNumTrips = 1;
    stats.dijkstraIters = EDijkstra::ITERS;
    std::map<uint32_t, double> colors;
    LOG(INFO) << "Matched 1 trip in " << std::fixed << std::setprecision(2)
              << stats.solveTime << " ms.";
    // print to line
    return {getGeom(hops, getRAttrs(trip), &colors, trip, 1), stats};
  } catch (const std::runtime_error& e) {
    LOG(ERROR) << e.what();
    return {std::vector<LINE>(), stats};
  }
}

// _____________________________________________________________________________
std::map<size_t, pfaedle::router::EdgeListHops> ShapeBuilder::route(
    const TripTrie<pfaedle::gtfs::Trip>* trie, const EdgeCandMap& ecm,
    HopCache* hopCache) const {
  return _router->route(trie, ecm, _motCfg.routingOpts, *_restr, hopCache,
                        _cfg.noFastHops);
}

// _____________________________________________________________________________
std::map<size_t, EdgeListHops> ShapeBuilder::shapeify(
    const TripTrie<pfaedle::gtfs::Trip>* trie, HopCache* hopCache) const {
  LOG(VDEBUG) << "Map-matching trie " << trie;

  assert(trie->getNdTrips().size());
  assert(trie->getNdTrips().begin()->second.size());
  RoutingAttrs rAttrs = getRAttrs(trie->getNdTrips().begin()->second[0]);

  std::map<size_t, EdgeListHops> ret;

  const auto& routes = route(trie, getECM(trie), hopCache);

  for (const auto& route : routes) {
    ret[route.first] = route.second;
  }

  LOG(VDEBUG) << "Finished map-matching for trie " << trie;

  return ret;
}

// _____________________________________________________________________________
EdgeListHops ShapeBuilder::shapeify(Trip* trip) {
  LOG(VDEBUG) << "Map-matching shape for trip #" << trip->getId() << " of mot "
              << trip->getRoute()->getType() << "(sn=" << trip->getShortname()
              << ", rsn=" << trip->getRoute()->getShortName()
              << ", rln=" << trip->getRoute()->getLongName() << ")";
  TripTrie<pfaedle::gtfs::Trip> trie;
  trie.addTrip(trip, getRAttrs(trip),
               _motCfg.routingOpts.transPenMethod == "timenorm", false);
  const auto& routes = route(&trie, getECM(&trie), 0);

  return routes.begin()->second;
}

// _____________________________________________________________________________
Stats ShapeBuilder::shapeify(pfaedle::netgraph::Graph* outNg) {
  Stats stats;
  EDijkstra::ITERS = 0;

  T_START(cluster);
  LOG(DEBUG) << "Clustering trips...";
  const TripForests& forests = clusterTrips(_feed, _mots);
  for (const auto& forest : forests) {
    for (const auto& trie : forest.second) {
      stats.numTries++;
      stats.numTrieLeafs += trie.getNdTrips().size();
    }
  }
  LOG(DEBUG) << "Clustered trips into " << stats.numTries
             << " tries with a total of " << stats.numTrieLeafs << " leafs in "
             << T_STOP(cluster) << "ms";

  LOG(DEBUG) << "Building candidate cache...";
  buildCandCache(forests);
  LOG(DEBUG) << "Done.";

  std::map<std::string, size_t> shpUse;
  RouteRefColors refColors;

  for (auto t : _feed->getTrips()) {
    if (!t.getShape().empty()) shpUse[t.getShape()]++;

    // write the colors of trips we won't touch, but whose route we might
    if (t.getStopTimes().size() < 2) continue;
    if (!_mots.count(t.getRoute()->getType()) ||
        !_motCfg.mots.count(t.getRoute()->getType()))
      continue;

    if (!t.getShape().empty() && !_cfg.dropShapes) {
      refColors[t.getRoute()][t.getRoute()->getColor()].push_back(&t);
    }
  }

  // we implicitely cluster by routing attrs here. This ensures that now two
  // threads will access the same routing attrs later on, which safes us an
  // expensive locking mechanism later on for the hop cache
  std::vector<const TripForest*> tries;
  for (const auto& forest : forests) {
    tries.push_back(&(forest.second));
    for (const auto& trie : forest.second) {
      for (const auto& trips : trie.getNdTrips()) {
        stats.totNumTrips += trips.second.size();
      }
    }
  }

  auto tStart = TIME();
  std::atomic<size_t> at(0);

  size_t numThreads = std::thread::hardware_concurrency();
  std::vector<std::thread> thrds(numThreads);
  std::vector<RouteRefColors> colors(numThreads);
  std::vector<TrGraphEdgs> gtfsGraphs(numThreads);

  size_t i = 0;
  for (auto& t : thrds) {
    t = std::thread(&ShapeBuilder::shapeWorker, this, &tries, &at, &shpUse,
                    &colors[i], &gtfsGraphs[i]);
    i++;
  }

  for (auto& thr : thrds) thr.join();

  stats.solveTime = TOOK(tStart, TIME());

  LOG(INFO) << "Matched " << stats.totNumTrips << " trips in " << std::fixed
            << std::setprecision(2) << stats.solveTime << " ms.";

  // merge colors
  for (auto& cols : colors) {
    for (auto& route : cols) {
      for (auto& col : route.second) {
        refColors[route.first][col.first].insert(
            refColors[route.first][col.first].end(), col.second.begin(),
            col.second.end());
      }
    }
  }

  // update them in the routes, split routes if necessary
  updateRouteColors(refColors);

  if (_cfg.buildTransitGraph) {
    LOG(DEBUG) << "Building transit network graph...";

    // merge gtfsgraph from threads
    TrGraphEdgs gtfsGraph;

    for (auto& g : gtfsGraphs) {
      for (auto& ePair : g) {
        gtfsGraph[ePair.first].insert(gtfsGraph[ePair.first].begin(),
                                      ePair.second.begin(), ePair.second.end());
      }
    }
    buildNetGraph(&gtfsGraph, outNg);
  }

  stats.dijkstraIters = EDijkstra::ITERS;

  return stats;
}

// _____________________________________________________________________________
void ShapeBuilder::updateRouteColors(const RouteRefColors& refColors) {
  for (auto& route : refColors) {
    if (route.second.size() == 1) {
      // only one color found for this route, great!
      // update inplace...
      route.first->setColor(route.second.begin()->first);
      if (route.first->getColor() != NO_COLOR)
        route.first->setTextColor(getTextColor(route.first->getColor()));
    } else {
      // are there fare rules using this route?
      std::vector<
          std::pair<ad::cppgtfs::gtfs::Fare<ad::cppgtfs::gtfs::Route>*,
                    ad::cppgtfs::gtfs::FareRule<ad::cppgtfs::gtfs::Route>>>
          rules;

      for (auto& f : _feed->getFares()) {
        for (auto r : f.second->getFareRules()) {
          if (r.getRoute() == route.first) {
            rules.push_back({f.second, r});
          }
        }
      }

      // add new routes...
      for (auto& c : route.second) {
        // keep the original one intact
        if (c.first == route.first->getColor()) continue;

        auto routeCp = *route.first;

        // find free id
        std::string newId = route.first->getId() + "::1";
        size_t i = 1;
        while (_feed->getRoutes().get(newId)) {
          i++;
          newId = route.first->getId() + "::" + std::to_string(i);
        }

        routeCp.setId(newId);
        routeCp.setColor(c.first);
        routeCp.setTextColor(getTextColor(routeCp.getColor()));

        auto newRoute = _feed->getRoutes().add(routeCp);

        // update trips to use that route
        for (auto& t : c.second) t->setRoute(newRoute);

        // add new fare rules
        for (auto a : rules) {
          auto rule = a.second;
          rule.setRoute(newRoute);
          a.first->addFareRule(rule);
        }
      }
    }
  }
}

// _____________________________________________________________________________
void ShapeBuilder::setShape(Trip* t, const ad::cppgtfs::gtfs::Shape& s,
                            const std::vector<float>& distances) {
  assert(distances.size() == t->getStopTimes().size());
  // set distances
  size_t i = 0;
  for (const auto& st : t->getStopTimes()) {
    const_cast<StopTime<Stop>&>(st).setShapeDistanceTravelled(distances[i]);
    i++;
  }

  std::lock_guard<std::mutex> guard(_shpMutex);
  auto gtfsShp = _feed->getShapes().add(s);
  t->setShape(gtfsShp);
}

// _____________________________________________________________________________
ad::cppgtfs::gtfs::Shape ShapeBuilder::getGtfsShape(
    const EdgeListHops& hops, Trip* t, size_t numOthers,
    const RoutingAttrs& rAttrs, std::vector<float>* hopDists,
    uint32_t* bestColor) {
  ad::cppgtfs::gtfs::Shape ret(getFreeShapeId(t));

  assert(hops.size() == t->getStopTimes().size() - 1);

  std::map<uint32_t, double> colors;

  const std::vector<LINE>& gl = getGeom(hops, rAttrs, &colors, t, numOthers);
  const std::vector<float>& measures = getMeasure(gl);

  size_t seq = 0;
  hopDists->push_back(0);
  for (size_t i = 0; i < gl.size(); i++) {
    for (size_t j = 0; j < gl[i].size(); j++) {
      ret.addPoint(
          ShapePoint(gl[i][j].getY(), gl[i][j].getX(), measures[seq], seq));
      seq++;
    }
    hopDists->push_back(measures[seq - 1]);
  }

  // get most likely color
  double best = 0;
  *bestColor = NO_COLOR;
  for (const auto& c : colors) {
    double progr = c.second / measures.back();
    // TODO(patrick): make threshold configurable
    if (progr > 0.9 && progr > best) {
      best = progr;
      *bestColor = c.first;
    }
  }

  return ret;
}

// _____________________________________________________________________________
std::string ShapeBuilder::getFreeShapeId(Trip* trip) {
  std::string ret;
  std::lock_guard<std::mutex> guard(_shpMutex);
  while (!ret.size() || _feed->getShapes().has(ret)) {
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

    ret.classifier = _classifier;

    const auto& lnormzer = _motCfg.osmBuildOpts.lineNormzer;
    const auto& snormzer = _motCfg.osmBuildOpts.statNormzer;

    ret.shortName = lnormzer.norm(trip->getRoute()->getShortName());
    ret.lineFrom =
        snormzer.norm(trip->getStopTimes().front().getStop()->getName());
    ret.lineTo = {
        snormzer.norm(trip->getStopTimes().back().getStop()->getName())};

    // fallbacks for line name
    if (ret.shortName.empty())
      ret.shortName = lnormzer.norm(trip->getShortname());

    if (ret.shortName.empty())
      ret.shortName = lnormzer.norm(trip->getRoute()->getLongName());

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
                              osm::BBoxIdx* box, double maxSpeed,
                              std::vector<double>* hopDists,
                              uint8_t verbosity) {
  for (const auto& t : feed->getTrips()) {
    if (!tid.empty() && t.getId() != tid) continue;
    if (tid.empty() && !t.getShape().empty() && !dropShapes) continue;
    if (t.getStopTimes().size() < 2) continue;

    if (mots.count(t.getRoute()->getType())) {
      DBox cur;
      for (size_t i = 0; i < t.getStopTimes().size(); i++) {
        // skip outliers
        const auto& st = t.getStopTimes()[i];

        int toTime = std::numeric_limits<int>::max();
        double toD = 0;
        int fromTime = std::numeric_limits<int>::max();
        double fromD = 0;

        if (i > 0) {
          const auto& stPrev = t.getStopTimes()[i - 1];
          toTime = st.getArrivalTime().seconds() -
                   stPrev.getDepartureTime().seconds();
          toD = util::geo::haversine(
              st.getStop()->getLat(), st.getStop()->getLng(),
              stPrev.getStop()->getLat(), stPrev.getStop()->getLng());
          if (hopDists) hopDists->push_back(toD);
        }

        if (i < t.getStopTimes().size() - 1) {
          const auto& stNext = t.getStopTimes()[i + 1];
          fromTime = stNext.getArrivalTime().seconds() -
                     st.getDepartureTime().seconds();
          fromD = util::geo::haversine(
              st.getStop()->getLat(), st.getStop()->getLng(),
              stNext.getStop()->getLat(), stNext.getStop()->getLng());
        }

        const double reqToTime = toD / maxSpeed;
        const double reqFromTime = fromD / maxSpeed;

        const double BUFFER = 5 * 60;

        if (reqToTime > (BUFFER + toTime) * 3 * MAX_ROUTE_COST_DOUBLING_STEPS &&
            reqFromTime >
                (BUFFER + fromTime) * 3 * MAX_ROUTE_COST_DOUBLING_STEPS) {
          if (verbosity) {
            LOG(WARN)
                << "Skipping station '" << st.getStop()->getName() << "' ("
                << st.getStop()->getId() << ") @ " << st.getStop()->getLat()
                << ", " << st.getStop()->getLng()
                << " for bounding box as the vehicle cannot realistically "
                   "reach and leave it in the scheduled time";
          } else {
            LOG(DEBUG)
                << "Skipping station '" << st.getStop()->getName() << "' ("
                << st.getStop()->getId() << ") @ " << st.getStop()->getLat()
                << ", " << st.getStop()->getLng()
                << " for bounding box as the vehicle cannot realistically "
                   "reach and leave it in the scheduled time";
          }
          continue;
        }

        cur = extendBox(DPoint(st.getStop()->getLng(), st.getStop()->getLat()),
                        cur);
      }
      box->add(cur);
    }
  }
}

// _____________________________________________________________________________
std::vector<double> ShapeBuilder::getTransTimes(Trip* trip) const {
  std::vector<double> ret;

  for (size_t i = 0; i < trip->getStopTimes().size() - 1; i++) {
    auto cur = trip->getStopTimes()[i];
    auto next = trip->getStopTimes()[i + 1];

    int depTime = cur.getDepartureTime().seconds();
    int arrTime = next.getArrivalTime().seconds();

    int diff = arrTime - depTime;
    if (diff < 1) diff = 1;

    ret.push_back(diff);
    assert(ret.back() >= 0);
  }

  return ret;
}

// _____________________________________________________________________________
std::vector<double> ShapeBuilder::getTransDists(Trip* trip) const {
  std::vector<double> ret;

  for (size_t i = 0; i < trip->getStopTimes().size() - 1; i++) {
    auto cur = trip->getStopTimes()[i];
    auto next = trip->getStopTimes()[i + 1];

    double dist = util::geo::haversine(
        cur.getStop()->getLat(), cur.getStop()->getLng(),
        next.getStop()->getLat(), next.getStop()->getLng());

    ret.push_back(dist);
  }

  return ret;
}

// _____________________________________________________________________________
EdgeCandMap ShapeBuilder::getECM(
    const TripTrie<pfaedle::gtfs::Trip>* trie) const {
  EdgeCandMap ecm(trie->getNds().size());

  for (size_t nid = 1; nid < trie->getNds().size(); nid++) {
    auto trNd = trie->getNds()[nid];
    auto parentTrNd = trie->getNds()[trNd.parent];

    if (nid != 1 && !trNd.arr) continue;

    double avgT = 0;

    if (trNd.trips) avgT = trNd.accTime / trNd.trips;

    const auto& cands = getEdgCands(trNd.reprStop);
    ecm[nid].reserve(cands.size());

    for (auto& cand : cands) {
      const auto& timeExpCands = timeExpand(cand, avgT);
      assert(timeExpCands.size());

      for (size_t depChildId : trNd.childs) {
        if (nid == 1) break;
        auto chldTrNd = trie->getNds()[depChildId];
        double avgChildT = 0;
        if (chldTrNd.trips) avgChildT = chldTrNd.accTime / chldTrNd.trips;

        double timeDiff = avgChildT - avgT;
        if (timeDiff < 0) timeDiff = 0;

        for (size_t candId = 0; candId < timeExpCands.size(); candId++) {
          const auto& cand = timeExpCands[candId];
          ecm[depChildId].push_back(cand);
          ecm[depChildId].back().time += timeDiff;

          ecm[depChildId].back().pen = timePen(cand.time, avgChildT);

          for (size_t sucCandId = 0; sucCandId < timeExpCands.size();
               sucCandId++) {
            if (timeExpCands[sucCandId].time <= ecm[depChildId].back().time) {
              ecm[depChildId].back().depPrede.push_back(sucCandId +
                                                        ecm[nid].size());
            }
          }
          assert(ecm[depChildId].back().depPrede.size());
        }
      }
      ecm[nid].insert(ecm[nid].end(), timeExpCands.begin(), timeExpCands.end());
    }

    assert(ecm[nid].size() != 0);
  }

  return ecm;
}

// _____________________________________________________________________________
double ShapeBuilder::timePen(int candTime, int schedTime) const {
  // standard deviation of normal distribution
  double standarddev = 5 * 60;

  int diff = abs(candTime - schedTime);

  double cNorm = diff / standarddev;
  return cNorm * cNorm;
}

// _____________________________________________________________________________
EdgeCandGroup ShapeBuilder::timeExpand(const EdgeCand& ec, int time) const {
  EdgeCandGroup ret;
  // TODO(patrick): heuristic for time expansion variance, currently
  // unused
  for (int i = 0; i < 1; i++) {
    EdgeCand ecNew = ec;
    // in 30 sec steps
    ecNew.time = time + i * 30;
    ecNew.pen = ecNew.pen + timePen(ecNew.time, time);
    ret.push_back(ecNew);
  }

  return ret;
}

// _____________________________________________________________________________
TripForests ShapeBuilder::clusterTrips(Feed* f, MOTs mots) {
  TripForests forest;
  std::map<RoutingAttrs, std::vector<Trip*>> trips;

  // warm the stop name normalizer caches so a
  // multithreaded access later on will never write to the underlying cache
  for (auto& stop : f->getStops()) {
    const auto& snormzer = _motCfg.osmBuildOpts.statNormzer;
    auto normedName = snormzer.norm(stop.getName());
  }

  // cluster by routing attr for parallization later on
  for (auto& trip : f->getTrips()) {
    if (!_cfg.dropShapes && !trip.getShape().empty()) continue;
    if (trip.getStopTimes().size() < 2) continue;
    if (!mots.count(trip.getRoute()->getType()) ||
        !_motCfg.mots.count(trip.getRoute()->getType()))
      continue;

    // important: we are building the routing attributes here, so a
    // multithreaded access later on will never write to the underlying cache
    const auto& rAttrs = getRAttrs(&trip);

    trips[rAttrs].push_back(&trip);
    forest[rAttrs] = {};
  }

  size_t numThreads = std::thread::hardware_concurrency();
  std::vector<std::thread> thrds(numThreads);
  std::vector<std::vector<RoutingAttrs>> attrs(numThreads);

  size_t i = 0;
  for (auto it : trips) {
    attrs[i].push_back(it.first);
    if (++i == numThreads) i = 0;
  }

  i = 0;
  for (auto& t : thrds) {
    t = std::thread(&ShapeBuilder::clusterWorker, this, &attrs[i], &trips,
                    &forest);
    i++;
  }

  for (auto& thr : thrds) thr.join();

  return forest;
}

// _____________________________________________________________________________
void ShapeBuilder::clusterWorker(
    const std::vector<RoutingAttrs>* rAttrsVec,
    const std::map<RoutingAttrs, std::vector<Trip*>>* trips,
    TripForests* forest) {
  for (const auto& rAttrs : *rAttrsVec) {
    for (auto& trip : trips->at(rAttrs)) {
      bool ins = false;
      auto& subForest = forest->at(rAttrs);
      for (auto& trie : subForest) {
        if (trie.addTrip(trip, rAttrs,
                         _motCfg.routingOpts.transPenMethod == "timenorm",
                         _cfg.noTrie)) {
          ins = true;
          break;
        }
      }

      if (!ins) {
        subForest.resize(subForest.size() + 1);
        subForest.back().addTrip(
            trip, rAttrs, _motCfg.routingOpts.transPenMethod == "timenorm",
            false);
      }
    }
  }
}

// _____________________________________________________________________________
const pfaedle::trgraph::Graph* ShapeBuilder::getGraph() const { return _g; }

// _____________________________________________________________________________
void ShapeBuilder::writeTransitGraph(
    const router::EdgeListHops& hops, TrGraphEdgs* edgs,
    const std::vector<pfaedle::gtfs::Trip*>& trips) const {
  for (const auto& hop : hops) {
    for (const auto* e : hop.edges) {
      if (e->pl().isRev()) e = _g->getEdg(e->getTo(), e->getFrom());
      (*edgs)[e].insert((*edgs)[e].begin(), trips.begin(), trips.end());
    }
  }
}

// _____________________________________________________________________________
void ShapeBuilder::buildNetGraph(TrGraphEdgs* edgs,
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

// _____________________________________________________________________________
std::vector<LINE> ShapeBuilder::getGeom(const EdgeListHops& hops,
                                        const RoutingAttrs& rAttrs,
                                        std::map<uint32_t, double>* colors,
                                        Trip* t, size_t numOthers) const {
  std::vector<LINE> ret;

  for (size_t i = hops.size(); i > 0; i--) {
    const auto& hop = hops[i - 1];
    if (!hop.start || !hop.end) {
      // no hop was found, use the fallback geometry

      if (_cfg.verbosity) {
        const auto stopFr = t->getStopTimes()[hops.size() - i].getStop();
        const auto stopTo = t->getStopTimes()[hops.size() - i + 1].getStop();

        LOG(WARN) << "No viable hop found between stops '" << stopFr->getName()
                  << "' (" << stopFr->getId() << ") and '" << stopTo->getName()
                  << "' (" << stopTo->getId() << ") for trip " << t->getId()
                  << " of type '"
                  << ad::cppgtfs::gtfs::flat::Route::getTypeString(
                         t->getRoute()->getType())
                  << "'"
                  << (numOthers > 1 ? " (and " + std::to_string(numOthers) +
                                          " similar trips)"
                                    : "")
                  << ", falling back to straight line";
      }

      if (hop.start) {
        if (hop.progrStart > 0) {
          auto l = getLine(hop.start);
          PolyLine<PFDL_PREC> pl(l);
          const auto& seg = pl.getSegment(hop.progrStart, 1);
          ret.push_back({seg.getLine().front(), hop.pointEnd});
        } else {
          ret.push_back({*hop.start->getFrom()->pl().getGeom(), hop.pointEnd});
        }
      } else if (hop.end) {
        if (hop.progrEnd > 0) {
          auto l = getLine(hop.end);
          PolyLine<PFDL_PREC> pl(l);
          const auto& seg = pl.getSegment(0, hop.progrEnd);
          ret.push_back({hop.pointStart, seg.getLine().back()});
        } else {
          ret.push_back({hop.pointStart, *hop.end->getFrom()->pl().getGeom()});
        }
      } else {
        ret.push_back({hop.pointStart, hop.pointEnd});
      }
    } else {
      const auto& l = getLine(hop, rAttrs, colors);
      ret.push_back(l);
    }
  }

  return ret;
}

// _____________________________________________________________________________
LINE ShapeBuilder::getLine(const EdgeListHop& hop, const RoutingAttrs& rAttrs,
                           std::map<uint32_t, double>* colors) const {
  LINE l;

  const auto& curL = getLine(hop.start);

  if (hop.edges.size() == 0) {
    // draw direct line between positions on edges
    if (hop.progrStart > 0) {
      PolyLine<PFDL_PREC> pl(curL);
      const auto& seg = pl.getSegment(hop.progrStart, 1);
      l.push_back(seg.front());
    } else {
      l.push_back(curL.front());
    }

    if (hop.progrEnd > 0) {
      PolyLine<PFDL_PREC> pl(getLine(hop.end));
      const auto& seg = pl.getSegment(0, hop.progrEnd);
      l.push_back(seg.back());
    } else {
      l.push_back(*hop.end->getFrom()->pl().getGeom());
    }

    return l;
  }

  // special case: start and end are on the same edge!
  if (hop.edges.size() == 1 && hop.start == hop.end) {
    PolyLine<PFDL_PREC> pl(curL);
    const auto& seg = pl.getSegment(hop.progrStart, hop.progrEnd);
    l.insert(l.end(), seg.getLine().begin(), seg.getLine().end());

    for (const auto& color : getColorMatch(hop.start, rAttrs)) {
      (*colors)[color] += hop.start->pl().getLength();
    }

    return l;
  }

  auto from = hop.start->getFrom();

  if (hop.progrStart > 0) {
    PolyLine<PFDL_PREC> pl(curL);
    const auto& seg = pl.getSegment(hop.progrStart, 1);
    l.insert(l.end(), seg.getLine().begin(), seg.getLine().end());

    double l = hop.start->pl().getLength() * (1 - hop.progrStart);
    for (const auto& color : getColorMatch(hop.start, rAttrs)) {
      (*colors)[color] += l;
    }
  } else {
    l.insert(l.end(), curL.begin(), curL.end());

    double l = hop.start->pl().getLength();
    for (const auto& color : getColorMatch(hop.start, rAttrs)) {
      (*colors)[color] += l;
    }
  }

  from = hop.start->getOtherNd(from);

  if (hop.edges.size() > 1) {
    for (size_t j = hop.edges.size() - 2; j > 0; j--) {
      const auto* e = hop.edges[j];
      const auto& curL = getLine(e);
      l.insert(l.end(), curL.begin(), curL.end());
      from = e->getOtherNd(from);

      double l = e->pl().getLength();
      for (const auto& color : getColorMatch(e, rAttrs)) {
        (*colors)[color] += l;
      }
    }
  }

  if (hop.progrEnd > 0) {
    PolyLine<PFDL_PREC> pl(getLine(hop.end));
    const auto& seg = pl.getSegment(0, hop.progrEnd);
    l.insert(l.end(), seg.getLine().begin(), seg.getLine().end());

    double l = hop.end->pl().getLength() * hop.progrEnd;
    for (const auto& color : getColorMatch(hop.end, rAttrs)) {
      (*colors)[color] += l;
    }
  }

  if (l.size() > 1) return util::geo::simplify(l, 0.5 / M_PER_DEG);
  return l;
}

// _____________________________________________________________________________
LINE ShapeBuilder::getLine(const trgraph::Edge* e) const {
  LINE l;
  if (!e->pl().getGeom() || e->pl().getGeom()->size() == 0)
    return {*e->getFrom()->pl().getGeom(), *e->getTo()->pl().getGeom()};
  if (e->pl().isRev()) {
    l.insert(l.end(), e->pl().getGeom()->rbegin(), e->pl().getGeom()->rend());
  } else {
    l.insert(l.end(), e->pl().getGeom()->begin(), e->pl().getGeom()->end());
  }
  return l;
}

// _____________________________________________________________________________
std::vector<float> ShapeBuilder::getMeasure(
    const std::vector<LINE>& lines) const {
  assert(lines.size());
  assert(lines.front().size());
  std::vector<float> ret;
  POINT last = lines.front().front();

  for (const auto& l : lines) {
    for (size_t i = 0; i < l.size(); i++) {
      if (ret.size() == 0) {
        ret.push_back(0);
      } else {
        float v = ret.back() + util::geo::haversine(last, l[i]);
        assert(v >= ret.back());  // required by GTFS standard!
        ret.push_back(v);
      }
      last = l[i];
    }
  }

  return ret;
}

// _____________________________________________________________________________
void ShapeBuilder::shapeWorker(
    const std::vector<const TripForest*>* tries, std::atomic<size_t>* at,
    std::map<std::string, size_t>* shpUse,
    std::map<Route*, std::map<uint32_t, std::vector<gtfs::Trip*>>>* routeColors,
    TrGraphEdgs* gtfsGraph) {
  while (1) {
    size_t j = (*at)++;
    if (j >= tries->size()) return;

    int step = tries->size() < 10 ? tries->size() : 10;

    if (j % (tries->size() / step) == 0) {
      LOG(INFO) << "@ " << (static_cast<int>((j * 1.0) / tries->size() * 100))
                << "%";
      LOG(DEBUG) << "(@ trie forest " << j << "/" << tries->size() << ")";
    }

    const auto& forest = *((*tries)[j]);

    // hop cache per forest, thus per routing attributes
    HopCache hopCacheLoc;
    HopCache* hopCache = 0;

    if (!_cfg.noHopCache) hopCache = &hopCacheLoc;

    for (size_t i = 0; i < forest.size(); i++) {
      const TripTrie<pfaedle::gtfs::Trip>* trie = &(forest[i]);
      const auto& hops = shapeify(trie, hopCache);

      for (const auto& leaf : trie->getNdTrips()) {
        std::vector<float> distances;
        const RoutingAttrs& rAttrs = trie->getNd(leaf.first).rAttrs;

        uint32_t color;

        const ad::cppgtfs::gtfs::Shape& shp =
            getGtfsShape(hops.at(leaf.first), leaf.second[0],
                         leaf.second.size(), rAttrs, &distances, &color);

        if (_cfg.buildTransitGraph) {
          writeTransitGraph(hops.at(leaf.first), gtfsGraph, leaf.second);
        }

        for (auto t : leaf.second) {
          if (_cfg.writeColors && color != NO_COLOR &&
              t->getRoute()->getColor() == NO_COLOR &&
              t->getRoute()->getTextColor() == NO_COLOR) {
            (*routeColors)[t->getRoute()][color].push_back(t);
          } else {
            // else, use the original route color
            (*routeColors)[t->getRoute()][t->getRoute()->getColor()].push_back(
                t);
          }

          if (!t->getShape().empty() && (*shpUse)[t->getShape()] > 0) {
            (*shpUse)[t->getShape()]--;
            if ((*shpUse)[t->getShape()] == 0) {
              std::lock_guard<std::mutex> guard(_shpMutex);
              _feed->getShapes().remove(t->getShape());
            }
          }
          setShape(t, shp, distances);
        }
      }
    }
  }
}

// _____________________________________________________________________________
void ShapeBuilder::edgCandWorker(std::vector<const Stop*>* stops,
                                 GrpCache* cache) {
  for (auto stop : *stops) {
    (*cache)[stop] = getEdgCands(stop);
  }
}

// _____________________________________________________________________________
std::set<uint32_t> ShapeBuilder::getColorMatch(
    const trgraph::Edge* e, const RoutingAttrs& rAttrs) const {
  std::set<uint32_t> ret;
  for (const auto* l : e->pl().getLines()) {
    auto simi = rAttrs.simi(l);
    if (simi.nameSimilar && l->color != NO_COLOR) ret.insert(l->color);
  }

  return ret;
}

// _____________________________________________________________________________
uint32_t ShapeBuilder::getTextColor(uint32_t c) const {
  double r = (c & 0x00FF0000) >> 16;
  double g = (c & 0x0000FF00) >> 8;
  double b = (c & 0x000000FF);

  // gray value approx
  double a = sqrt((r * r + g * g + b * b) / 3);

  // below a certain gray value, use white, else black
  if (a < 140) return 0x00FFFFFF;
  return 0;
}

// _____________________________________________________________________________
double ShapeBuilder::emWeight(double mDist) const {
  if (_motCfg.routingOpts.emPenMethod == "exp") {
    return mDist * _motCfg.routingOpts.stationDistPenFactor;
  }

  if (_motCfg.routingOpts.emPenMethod == "norm") {
    double s = mDist * _motCfg.routingOpts.stationDistPenFactor;
    return 0.5 * s * s;
  }

  return mDist;
}

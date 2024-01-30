// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_SHAPEBUILDER_H_
#define PFAEDLE_ROUTER_SHAPEBUILDER_H_

#include <map>
#include <mutex>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "ad/cppgtfs/gtfs/Feed.h"
#include "pfaedle/Def.h"
#include "pfaedle/config/MotConfig.h"
#include "pfaedle/config/PfaedleConfig.h"
#include "pfaedle/gtfs/Feed.h"
#include "pfaedle/netgraph/Graph.h"
#include "pfaedle/osm/Restrictor.h"
#include "pfaedle/router/Misc.h"
#include "pfaedle/router/Router.h"
#include "pfaedle/router/Stats.h"
#include "pfaedle/router/TripTrie.h"
#include "pfaedle/statsimi-classifier/StatsimiClassifier.h"
#include "pfaedle/trgraph/Graph.h"
#include "util/geo/Geo.h"

namespace pfaedle {
namespace router {

typedef std::vector<TripTrie<pfaedle::gtfs::Trip>> TripForest;
typedef std::map<router::RoutingAttrs, TripForest> TripForests;
typedef std::pair<const ad::cppgtfs::gtfs::Stop*,
                  const ad::cppgtfs::gtfs::Stop*>
    StopPair;
typedef std::unordered_map<const pfaedle::gtfs::Trip*, router::RoutingAttrs>
    TripRAttrs;
typedef std::unordered_map<const trgraph::Edge*,
                           std::vector<const pfaedle::gtfs::Trip*>>
    TrGraphEdgs;
typedef std::map<Route*, std::map<uint32_t, std::vector<gtfs::Trip*>>>
    RouteRefColors;
typedef std::unordered_map<const ad::cppgtfs::gtfs::Stop*, EdgeCandGroup>
    GrpCache;

/*
 * Layer class for the router. Provides an interface for direct usage with
 * GTFS data
 */
class ShapeBuilder {
 public:
  ShapeBuilder(
      pfaedle::gtfs::Feed* feed, MOTs mots, const config::MotConfig& motCfg,
      trgraph::Graph* g, router::FeedStops* stops, osm::Restrictor* restr,
      const pfaedle::statsimiclassifier::StatsimiClassifier* classifier,
      router::Router* router, const config::Config& cfg);

  Stats shapeify(pfaedle::netgraph::Graph* outNg);

  router::FeedStops* getFeedStops();

  // shape single trip
  std::pair<std::vector<LINE>, Stats> shapeL(pfaedle::gtfs::Trip* trip);

  std::map<size_t, EdgeListHops> shapeify(
      const TripTrie<pfaedle::gtfs::Trip>* trie, HopCache* hopCache) const;
  EdgeListHops shapeify(pfaedle::gtfs::Trip* trip);

  const trgraph::Graph* getGraph() const;

  static void getGtfsBox(const pfaedle::gtfs::Feed* feed, const MOTs& mots,
                         const std::string& tid, bool dropShapes,
                         osm::BBoxIdx* box, double maxSpeed,
                         std::vector<double>* hopDists, uint8_t verbosity);

 private:
  pfaedle::gtfs::Feed* _feed;
  MOTs _mots;
  config::MotConfig _motCfg;
  config::Config _cfg;
  trgraph::Graph* _g;
  router::FeedStops* _stops;

  EdgeCandGroup _emptyNCG;

  size_t _curShpCnt;

  std::mutex _shpMutex;

  TripRAttrs _rAttrs;

  osm::Restrictor* _restr;
  const pfaedle::statsimiclassifier::StatsimiClassifier* _classifier;
  GrpCache _grpCache;

  router::Router* _router;

  TripForests clusterTrips(pfaedle::gtfs::Feed* f, MOTs mots);
  void buildNetGraph(TrGraphEdgs* edgs, pfaedle::netgraph::Graph* ng) const;

  std::string getFreeShapeId(pfaedle::gtfs::Trip* t);
  ad::cppgtfs::gtfs::Shape getGtfsShape(const EdgeListHops& shp,
                                        pfaedle::gtfs::Trip* t,
                                        size_t numOthers,
                                        const RoutingAttrs& rAttrs,
                                        std::vector<float>* hopDists,
                                        uint32_t* bestColor);

  void setShape(pfaedle::gtfs::Trip* t, const ad::cppgtfs::gtfs::Shape& s,
                const std::vector<float>& dists);

  EdgeCandGroup getEdgCands(const ad::cppgtfs::gtfs::Stop* s) const;

  router::EdgeCandMap getECM(const TripTrie<pfaedle::gtfs::Trip>* trie) const;
  std::vector<double> getTransTimes(pfaedle::gtfs::Trip* trip) const;
  std::vector<double> getTransDists(pfaedle::gtfs::Trip* trip) const;
  const router::RoutingAttrs& getRAttrs(const pfaedle::gtfs::Trip* trip) const;
  const router::RoutingAttrs& getRAttrs(const pfaedle::gtfs::Trip* trip);
  std::map<size_t, router::EdgeListHops> route(
      const TripTrie<pfaedle::gtfs::Trip>* trie, const EdgeCandMap& ecm,
      HopCache* hopCache) const;
  double emWeight(double mDist) const;

  void buildCandCache(const TripForests& clusters);
  void buildIndex();

  std::vector<LINE> getGeom(const EdgeListHops& shp, const RoutingAttrs& rAttrs,
                            std::map<uint32_t, double>* colors, Trip* t,
                            size_t numOthers) const;
  double timePen(int candTime, int schedTime) const;

  LINE getLine(const EdgeListHop& hop, const RoutingAttrs&,
               std::map<uint32_t, double>* colMap) const;
  LINE getLine(const trgraph::Edge* edg) const;
  std::vector<float> getMeasure(const std::vector<LINE>& lines) const;

  trgraph::Edge* deg2reachable(trgraph::Edge* e,
                               std::set<trgraph::Edge*> edgs) const;

  EdgeCandGroup timeExpand(const EdgeCand& ec, int time) const;

  std::set<uint32_t> getColorMatch(const trgraph::Edge* e,
                                   const RoutingAttrs& rAttrs) const;

  void updateRouteColors(const RouteRefColors& c);

  uint32_t getTextColor(uint32_t c) const;

  void writeTransitGraph(const router::EdgeListHops& shp, TrGraphEdgs* edgs,
                         const std::vector<pfaedle::gtfs::Trip*>& trips) const;

  void shapeWorker(
      const std::vector<const TripForest*>* tries, std::atomic<size_t>* at,
      std::map<std::string, size_t>* shpUsage,
      std::map<Route*, std::map<uint32_t, std::vector<gtfs::Trip*>>>*,
      TrGraphEdgs* gtfsGraph);

  void edgCandWorker(std::vector<const Stop*>* stops, GrpCache* cache);
  void clusterWorker(const std::vector<RoutingAttrs>* rAttrs,
                     const std::map<RoutingAttrs, std::vector<Trip*>>* trips,
                     TripForests* forest);

  pfaedle::trgraph::EdgeGrid _eGrid;
  pfaedle::trgraph::NodeGrid _nGrid;
};

}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_SHAPEBUILDER_H_

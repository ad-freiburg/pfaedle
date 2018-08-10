// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_SHAPEBUILDER_H_
#define PFAEDLE_ROUTER_SHAPEBUILDER_H_

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
#include "pfaedle/eval/Collector.h"
#include "pfaedle/netgraph/Graph.h"
#include "pfaedle/osm/Restrictor.h"
#include "pfaedle/router/Misc.h"
#include "pfaedle/router/Router.h"
#include "pfaedle/trgraph/Graph.h"
#include "util/geo/Geo.h"

namespace pfaedle {
namespace router {

using ad::cppgtfs::gtfs::Stop;
using ad::cppgtfs::gtfs::Trip;
using ad::cppgtfs::gtfs::Feed;

struct Shape {
  router::EdgeListHops hops;
  double avgHopDist;
};

typedef std::vector<Trip*> Cluster;
typedef std::vector<Cluster> Clusters;
typedef std::pair<const Stop*, const Stop*> StopPair;
typedef std::unordered_map<const Trip*, router::RoutingAttrs> TripRAttrs;
typedef std::unordered_map<const trgraph::Edge*, std::set<const Trip*>>
    TrGraphEdgs;

/*
 * Layer class for the router. Provides an interface for direct usage with
 * GTFS data
 */
class ShapeBuilder {
 public:
  ShapeBuilder(Feed* feed, MOTs mots, const config::MotConfig& motCfg,
               eval::Collector* ecoll, const config::Config& cfg);

  void shape(pfaedle::netgraph::Graph* ng);

  router::FeedStops* getFeedStops();

  const NodeCandGroup& getNodeCands(const Stop* s) const;

  LINE shapeL(const router::NodeCandRoute& ncr,
              const router::RoutingAttrs& rAttrs);
  LINE shapeL(Trip* trip);

  pfaedle::router::Shape shape(Trip* trip) const;
  pfaedle::router::Shape shape(Trip* trip);

  const trgraph::Graph* getGraph() const;

  static osm::BBoxIdx getPaddedGtfsBox(const Feed* feed, double pad,
                                       const MOTs& mots, const std::string& tid,
                                       bool dropShapes);

 private:
  Feed* _feed;
  MOTs _mots;
  config::MotConfig _motCfg;
  eval::Collector* _ecoll;
  config::Config _cfg;
  trgraph::Graph _g;
  router::Router _crouter;

  router::FeedStops _stops;

  NodeCandGroup _emptyNCG;

  size_t _curShpCnt, _numThreads;

  std::mutex _shpMutex;

  TripRAttrs _rAttrs;

  osm::Restrictor _restr;

  void writeMotStops();
  void buildGraph();

  Clusters clusterTrips(Feed* f, MOTs mots);
  void writeTransitGraph(const Shape& shp, TrGraphEdgs* edgs,
                         const Cluster& cluster) const;
  void buildTrGraph(TrGraphEdgs* edgs, pfaedle::netgraph::Graph* ng) const;

  std::string getFreeShapeId(Trip* t);

  ad::cppgtfs::gtfs::Shape* getGtfsShape(const Shape& shp, Trip* t,
                                         std::vector<double>* hopDists);

  void setShape(Trip* t, ad::cppgtfs::gtfs::Shape* s,
                const std::vector<double>& dists);

  router::NodeCandRoute getNCR(Trip* trip) const;
  double avgHopDist(Trip* trip) const;
  const router::RoutingAttrs& getRAttrs(const Trip* trip) const;
  const router::RoutingAttrs& getRAttrs(const Trip* trip);
  bool routingEqual(Trip* a, Trip* b);
  bool routingEqual(const Stop* a, const Stop* b);
  router::EdgeListHops route(const router::NodeCandRoute& ncr,
                             const router::RoutingAttrs& rAttrs) const;
};
}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_SHAPEBUILDER_H_

// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_MISC_H_
#define PFAEDLE_ROUTER_MISC_H_

#include <set>
#include <string>
#include <vector>
#include <unordered_map>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "ad/cppgtfs/gtfs/Route.h"
#include "pfaedle/trgraph/Graph.h"

using ad::cppgtfs::gtfs::Route;
using ad::cppgtfs::gtfs::Stop;

namespace pfaedle {
namespace router {

struct NodeCand {
  trgraph::Node* nd;
  double pen;
};

struct RoutingOpts {
  RoutingOpts()
      : fullTurnPunishFac(2000),
        fullTurnAngle(45),
        passThruStationsPunish(100),
        oneWayPunishFac(1),
        oneWayEdgePunish(0),
        lineUnmatchedPunishFact(0.5),
        platformUnmatchedPen(0),
        stationDistPenFactor(0) {}
  double fullTurnPunishFac;
  double fullTurnAngle;
  double passThruStationsPunish;
  double oneWayPunishFac;
  double oneWayEdgePunish;
  double lineUnmatchedPunishFact;
  double platformUnmatchedPen;
  double stationDistPenFactor;
  double nonOsmPen;
  double levelPunish[8];
};

inline bool operator==(const RoutingOpts& a, const RoutingOpts& b) {
  return fabs(a.fullTurnPunishFac - b.fullTurnPunishFac) < 0.01 &&
         fabs(a.fullTurnAngle - b.fullTurnAngle) < 0.01 &&
         fabs(a.passThruStationsPunish - b.passThruStationsPunish) < 0.01 &&
         fabs(a.oneWayPunishFac - b.oneWayPunishFac) < 0.01 &&
         fabs(a.oneWayEdgePunish - b.oneWayEdgePunish) < 0.01 &&
         fabs(a.lineUnmatchedPunishFact - b.lineUnmatchedPunishFact) < 0.01 &&
         fabs(a.platformUnmatchedPen - b.platformUnmatchedPen) < 0.01 &&
         fabs(a.stationDistPenFactor - b.stationDistPenFactor) < 0.01 &&
         fabs(a.nonOsmPen - b.nonOsmPen) < 0.01 &&
         fabs(a.levelPunish[0] - b.levelPunish[0]) < 0.01 &&
         fabs(a.levelPunish[1] - b.levelPunish[1]) < 0.01 &&
         fabs(a.levelPunish[2] - b.levelPunish[2]) < 0.01 &&
         fabs(a.levelPunish[3] - b.levelPunish[3]) < 0.01 &&
         fabs(a.levelPunish[4] - b.levelPunish[4]) < 0.01 &&
         fabs(a.levelPunish[5] - b.levelPunish[5]) < 0.01 &&
         fabs(a.levelPunish[6] - b.levelPunish[6]) < 0.01 &&
         fabs(a.levelPunish[7] - b.levelPunish[7]) < 0.01;
}

struct EdgeCost {
  EdgeCost()
      : meterDist(0),
        meterDistLvl1(0),
        meterDistLvl2(0),
        meterDistLvl3(0),
        meterDistLvl4(0),
        meterDistLvl5(0),
        meterDistLvl6(0),
        meterDistLvl7(0),
        fullTurns(0),
        passThruStations(0),
        oneWayMeters(0),
        oneWayEdges(0),
        lineUnmatchedMeters(0),
        reachPen(0),
        o(0) {}
  EdgeCost(double mDist, double mDistLvl1, double mDistLvl2, double mDistLvl3,
           double mDistLvl4, double mDistLvl5, double mDistLvl6,
           double mDistLvl7, uint32_t fullTurns, int32_t passThru,
           double oneWayMeters, size_t oneWayEdges, double lineUnmatchedMeters,
           double reachPen, const RoutingOpts* o)
      : meterDist(mDist),
        meterDistLvl1(mDistLvl1),
        meterDistLvl2(mDistLvl2),
        meterDistLvl3(mDistLvl3),
        meterDistLvl4(mDistLvl4),
        meterDistLvl5(mDistLvl5),
        meterDistLvl6(mDistLvl6),
        meterDistLvl7(mDistLvl7),
        fullTurns(fullTurns),
        passThruStations(passThru),
        oneWayMeters(oneWayMeters),
        oneWayEdges(oneWayEdges),
        lineUnmatchedMeters(lineUnmatchedMeters),
        reachPen(reachPen),
        o(o) {}
  double meterDist;
  double meterDistLvl1;
  double meterDistLvl2;
  double meterDistLvl3;
  double meterDistLvl4;
  double meterDistLvl5;
  double meterDistLvl6;
  double meterDistLvl7;
  uint32_t fullTurns;
  int32_t passThruStations;
  double oneWayMeters;
  size_t oneWayEdges;
  double lineUnmatchedMeters;
  double reachPen;
  const RoutingOpts* o;

  double getValue() const {
    if (!o) return meterDist + reachPen;
    return meterDist * o->levelPunish[0] + meterDistLvl1 * o->levelPunish[1] +
           meterDistLvl2 * o->levelPunish[2] +
           meterDistLvl3 * o->levelPunish[3] +
           meterDistLvl4 * o->levelPunish[4] +
           meterDistLvl5 * o->levelPunish[5] +
           meterDistLvl6 * o->levelPunish[6] +
           meterDistLvl7 * o->levelPunish[7] +
           oneWayMeters * o->oneWayPunishFac +
           oneWayEdges * o->oneWayEdgePunish +
           lineUnmatchedMeters * o->lineUnmatchedPunishFact +
           fullTurns * o->fullTurnPunishFac +
           passThruStations * o->passThruStationsPunish + reachPen;
  }

  double getTotalMeters() const {
    return meterDist + meterDistLvl1 + meterDistLvl2 + meterDistLvl3 +
           meterDistLvl4 + meterDistLvl5 + meterDistLvl6 + meterDistLvl7;
  }
};

inline EdgeCost operator+(const EdgeCost& a, const EdgeCost& b) {
  return EdgeCost(
      a.meterDist + b.meterDist, a.meterDistLvl1 + b.meterDistLvl1,
      a.meterDistLvl2 + b.meterDistLvl2, a.meterDistLvl3 + b.meterDistLvl3,
      a.meterDistLvl4 + b.meterDistLvl4, a.meterDistLvl5 + b.meterDistLvl5,
      a.meterDistLvl6 + b.meterDistLvl6, a.meterDistLvl7 + b.meterDistLvl7,
      a.fullTurns + b.fullTurns, a.passThruStations + b.passThruStations,
      a.oneWayMeters + b.oneWayMeters, a.oneWayEdges + b.oneWayEdges,
      a.lineUnmatchedMeters + b.lineUnmatchedMeters, a.reachPen + b.reachPen,
      a.o ? a.o : b.o);
}

inline bool operator<=(const EdgeCost& a, const EdgeCost& b) {
  return a.getValue() <= b.getValue();
}

inline bool operator==(const EdgeCost& a, const EdgeCost& b) {
  return a.getValue() == b.getValue();
}

inline bool operator>(const EdgeCost& a, const EdgeCost& b) {
  return a.getValue() > b.getValue();
}

template<typename F>
inline bool angSmaller(const Point<F>& f, const Point<F>& m, const Point<F>& t,
                      double ang) {
  if (util::geo::innerProd(m, f, t) < ang) return 1;
  return 0;
}

typedef std::set<trgraph::Node*> NodeSet;
typedef std::set<trgraph::Edge*> EdgeSet;
typedef std::unordered_map<const Stop*, trgraph::Node*> FeedStops;

typedef std::vector<NodeCand> NodeCandGroup;
typedef std::vector<NodeCandGroup> NodeCandRoute;

typedef std::vector<trgraph::Edge*> EdgeList;
typedef std::vector<trgraph::Node*> NodeList;

struct EdgeListHop {
  EdgeList edges;
  const trgraph::Node* start;
  const trgraph::Node* end;
};

typedef std::vector<EdgeListHop> EdgeListHops;

typedef std::set<Route::TYPE> MOTs;
}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_MISC_H_

// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_MISC_H_
#define PFAEDLE_ROUTER_MISC_H_

#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "ad/cppgtfs/gtfs/Feed.h"
#include "ad/cppgtfs/gtfs/Route.h"
#include "pfaedle/gtfs/Feed.h"
#include "pfaedle/trgraph/Graph.h"
#include "util/Nullable.h"

using ad::cppgtfs::gtfs::Route;
using ad::cppgtfs::gtfs::Stop;

namespace pfaedle {
namespace router {

extern double time;

struct EdgeCand {
  trgraph::Edge* e;
  double pen;
  double progr;
  POINT point;
  int time;
  std::vector<size_t> depPrede;
};

struct RoutingOpts {
  RoutingOpts()
      : fullTurnPunishFac(1000),
        fullTurnAngle(45),
        lineUnmatchedPunishFact(1),
        lineNameFromUnmatchedPunishFact(1),
        lineNameToUnmatchedPunishFact(1),
        noLinesPunishFact(1),
        platformUnmatchedPen(0),
        stationDistPenFactor(0),
        turnRestrCost(0),
        popReachEdge(true),
        noSelfHops(true) {}
  uint32_t fullTurnPunishFac;
  double fullTurnAngle;
  double lineUnmatchedPunishFact;
  double lineNameFromUnmatchedPunishFact;
  double lineNameToUnmatchedPunishFact;
  double noLinesPunishFact;
  double platformUnmatchedPen;
  double stationUnmatchedPen;
  double stationDistPenFactor;
  double nonStationPen;
  uint32_t turnRestrCost;
  bool popReachEdge;
  bool noSelfHops;
  bool useStations;
  double transitionPen;
  std::string transPenMethod;
  std::string emPenMethod;
  std::string statsimiMethod;
};

// _____________________________________________________________________________
inline bool operator==(const RoutingOpts& a, const RoutingOpts& b) {
  return a.fullTurnPunishFac == b.fullTurnPunishFac &&
         fabs(a.fullTurnAngle - b.fullTurnAngle) < 0.01 &&
         fabs(a.lineUnmatchedPunishFact - b.lineUnmatchedPunishFact) < 0.01 &&
         fabs(a.lineNameFromUnmatchedPunishFact -
              b.lineNameFromUnmatchedPunishFact) < 0.01 &&
         fabs(a.lineNameToUnmatchedPunishFact -
              b.lineNameToUnmatchedPunishFact) < 0.01 &&
         fabs(a.noLinesPunishFact - b.noLinesPunishFact) < 0.01 &&
         fabs(a.platformUnmatchedPen - b.platformUnmatchedPen) < 0.01 &&
         fabs(a.stationUnmatchedPen - b.stationUnmatchedPen) < 0.01 &&
         fabs(a.stationDistPenFactor - b.stationDistPenFactor) < 0.01 &&
         a.turnRestrCost == b.turnRestrCost &&
         fabs(a.transitionPen - b.transitionPen) < 0.01 &&
         fabs(a.nonStationPen - b.nonStationPen) < 0.01 &&
         a.transPenMethod == b.transPenMethod &&
         a.emPenMethod == b.emPenMethod &&
         a.statsimiMethod == b.statsimiMethod &&
         a.useStations == b.useStations && a.popReachEdge == b.popReachEdge &&
         a.noSelfHops == b.noSelfHops;
}

typedef std::set<trgraph::Node*> NodeSet;
typedef std::set<trgraph::Edge*> EdgeSet;
typedef std::unordered_map<const Stop*, trgraph::Node*> FeedStops;

typedef std::vector<EdgeCand> EdgeCandGroup;
typedef std::vector<EdgeCandGroup> EdgeCandMap;
typedef std::vector<EdgeCandGroup> EdgeCandRoute;

typedef std::vector<trgraph::Edge*> EdgeList;
typedef std::vector<trgraph::Node*> NodeList;

struct EdgeListHop {
  EdgeList edges;
  const trgraph::Edge* start;
  const trgraph::Edge* end;
  double progrStart;
  double progrEnd;
  POINT pointStart;
  POINT pointEnd;
};

typedef std::vector<EdgeListHop> EdgeListHops;

typedef std::set<Route::TYPE> MOTs;

// _____________________________________________________________________________
inline MOTs motISect(const MOTs& a, const MOTs& b) {
  MOTs ret;
  for (auto mot : a)
    if (b.count(mot)) ret.insert(mot);
  return ret;
}

// _____________________________________________________________________________
inline pfaedle::router::FeedStops writeMotStops(const pfaedle::gtfs::Feed* feed,
                                                const MOTs mots,
                                                const std::string& tid) {
  pfaedle::router::FeedStops ret;
  for (auto t : feed->getTrips()) {
    if (!tid.empty() && t.getId() != tid) continue;
    if (mots.count(t.getRoute()->getType())) {
      for (auto st : t.getStopTimes()) {
        // if the station has type STATION_ENTRANCE, use the parent
        // station for routing. Normally, this should not occur, as
        // this is not allowed in stop_times.txt
        if (st.getStop()->getLocationType() ==
                ad::cppgtfs::gtfs::flat::Stop::STATION_ENTRANCE &&
            st.getStop()->getParentStation()) {
          ret[st.getStop()->getParentStation()] = 0;
        } else {
          ret[st.getStop()] = 0;
        }
      }
    }
  }
  return ret;
}

// _____________________________________________________________________________
inline std::string getMotStr(const MOTs& mots) {
  MOTs tmp = mots;
  bool first = false;
  std::string motStr;

  std::string names[11] = {"tram",  "subway",     "rail",    "bus",
                           "ferry", "cablecar",   "gondola", "funicular",
                           "coach", "trolleybus", "monorail"};

  for (const auto& n : names) {
    const auto& types = ad::cppgtfs::gtfs::flat::Route::getTypesFromString(n);
    const auto& isect = motISect(tmp, types);

    if (isect.size() == types.size()) {
      if (first) motStr += ", ";
      motStr += "{" + n + "}";
      first = true;
      for (const auto& mot : isect) tmp.erase(mot);
    }
  }

  for (const auto& mot : tmp) {
    if (first) motStr += ", ";
    motStr += "<" + ad::cppgtfs::gtfs::flat::Route::getTypeString(mot) + ">";
    first = true;
  }

  return motStr;
}
}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_MISC_H_

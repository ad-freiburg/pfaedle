// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_OSM_OSMREADOPTS_H_
#define PFAEDLE_OSM_OSMREADOPTS_H_

#include <queue>
#include <unordered_set>
#include <string>
#include <unordered_map>
#include <map>
#include <utility>
#include <vector>
#include <set>
#include "pfaedle/osm/Osm.h"
#include "pfaedle/trgraph/Graph.h"
#include "pfaedle/trgraph/Normalizer.h"

namespace pfaedle {
namespace osm {

typedef std::unordered_set<std::string> AttrKeySet;
typedef std::unordered_map<osmid, trgraph::Node*> NIdMap;
typedef std::unordered_map<osmid, std::set<trgraph::Node*>> NIdMultMap;
typedef std::pair<double, trgraph::Edge*> EdgeCand;
typedef std::priority_queue<EdgeCand> EdgeCandPQ;
typedef std::unordered_map<osmid, std::vector<size_t>> RelMap;
typedef std::vector<AttrMap> RelVec;
typedef std::vector<std::string> AttrLst;

typedef std::pair<std::string, uint64_t> AttrFlagPair;
typedef std::unordered_map<std::string, std::map<std::string, uint64_t>>
    MultAttrMap;

typedef std::pair<std::string, std::string> KeyVal;
typedef std::set<size_t> FlatRels;

typedef std::unordered_map<const trgraph::Edge*, std::string> EdgTracks;

struct RelLst {
  RelVec rels;
  FlatRels flat;
};

enum FilterFlags : uint64_t {
  USE = 1,  // dummy value
  REL_NO_DOWN = 2,
  NO_RELATIONS = 4,
  NO_WAYS = 8,
  NO_NODES = 16,
  MULT_VAL_MATCH = 32
};

struct FilterRule {
  FilterRule() : kv(KeyVal("", "")) {}
  KeyVal kv;
  std::set<std::string> flags;
};

inline bool operator==(const FilterRule& a, const FilterRule& b) {
  return a.kv == b.kv && a.flags == b.flags;
}

struct DeepAttrRule {
  std::string attr;
  FilterRule relRule;
};

inline bool operator==(const DeepAttrRule& a, const DeepAttrRule& b) {
  return a.attr == b.attr && a.relRule == b.relRule;
}

typedef std::vector<DeepAttrRule> DeepAttrLst;

struct RelLineRules {
  AttrLst sNameRule;
  AttrLst fromNameRule;
  AttrLst toNameRule;
};

inline bool operator==(const RelLineRules& a, const RelLineRules& b) {
  return a.sNameRule == b.sNameRule && a.fromNameRule == b.fromNameRule &&
         a.toNameRule == b.toNameRule;
}

struct StationAttrRules {
  DeepAttrLst nameRule;
  DeepAttrLst platformRule;
};

inline bool operator==(const StationAttrRules& a, const StationAttrRules& b) {
  return a.nameRule == b.nameRule && a.platformRule == b.platformRule;
}

struct StatGroupNAttrRule {
  DeepAttrRule attr;
  double maxDist;
};

inline bool operator==(const StatGroupNAttrRule& a,
                       const StatGroupNAttrRule& b) {
  return a.attr == b.attr && a.maxDist == b.maxDist;
}

typedef std::unordered_map<
    std::string,
    std::unordered_map<std::string, std::vector<trgraph::StatGroup*>>>
    StAttrGroups;

struct OsmReadOpts {
  OsmReadOpts() {}

  MultAttrMap noHupFilter;
  MultAttrMap keepFilter;
  MultAttrMap levelFilters[8];
  MultAttrMap dropFilter;
  MultAttrMap oneWayFilter;
  MultAttrMap oneWayFilterRev;
  MultAttrMap twoWayFilter;
  MultAttrMap stationFilter;
  MultAttrMap stationBlockerFilter;
  std::vector<StatGroupNAttrRule> statGroupNAttrRules;

  trgraph::Normalizer statNormzer;
  trgraph::Normalizer lineNormzer;
  trgraph::Normalizer trackNormzer;

  RelLineRules relLinerules;
  StationAttrRules statAttrRules;

  DeepAttrLst edgePlatformRules;

  uint8_t maxSnapLevel;

  double maxAngleSnapReach;
  std::vector<double> maxSnapDistances;
  double maxSnapFallbackHeurDistance;
  double maxBlockDistance;

  double maxOsmStationDistance;

  // TODO(patrick): this is not implemented yet
  double levelSnapPunishFac[7] = {0, 0, 0, 0, 0, 0, 0};

  double fullTurnAngle;

  // restriction system
  MultAttrMap restrPosRestr;
  MultAttrMap restrNegRestr;
  MultAttrMap noRestrFilter;
};

inline bool operator==(const OsmReadOpts& a, const OsmReadOpts& b) {
  if (a.maxSnapDistances.size() != b.maxSnapDistances.size()) return false;
  for (size_t i = 0; i < a.maxSnapDistances.size(); i++) {
    if (fabs(a.maxSnapDistances[i] - b.maxSnapDistances[i]) >= 0.1)
      return false;
  }

  return a.noHupFilter == b.noHupFilter && a.keepFilter == b.keepFilter &&
         a.levelFilters[0] == b.levelFilters[0] &&
         a.levelFilters[1] == b.levelFilters[1] &&
         a.levelFilters[2] == b.levelFilters[2] &&
         a.levelFilters[3] == b.levelFilters[3] &&
         a.levelFilters[4] == b.levelFilters[4] &&
         a.levelFilters[5] == b.levelFilters[5] &&
         a.levelFilters[6] == b.levelFilters[6] &&
         a.dropFilter == b.dropFilter && a.oneWayFilter == b.oneWayFilter &&
         a.oneWayFilterRev == b.oneWayFilterRev &&
         a.twoWayFilter == b.twoWayFilter &&
         a.stationFilter == b.stationFilter &&
         a.stationBlockerFilter == b.stationBlockerFilter &&
         a.statGroupNAttrRules == b.statGroupNAttrRules &&
         a.statNormzer == b.statNormzer && a.lineNormzer == b.lineNormzer &&
         a.trackNormzer == b.trackNormzer && a.relLinerules == b.relLinerules &&
         a.statAttrRules == b.statAttrRules &&
         a.maxSnapLevel == b.maxSnapLevel &&
         fabs(a.maxAngleSnapReach - b.maxAngleSnapReach) < 0.1 &&
         fabs(a.maxOsmStationDistance - b.maxOsmStationDistance) < 0.1 &&
         fabs(a.maxSnapFallbackHeurDistance - b.maxSnapFallbackHeurDistance) <
             0.1 &&
         fabs(a.maxBlockDistance - b.maxBlockDistance) < 0.1 &&
         fabs(a.levelSnapPunishFac[0] - b.levelSnapPunishFac[0]) < 0.1 &&
         fabs(a.levelSnapPunishFac[1] - b.levelSnapPunishFac[1]) < 0.1 &&
         fabs(a.levelSnapPunishFac[2] - b.levelSnapPunishFac[2]) < 0.1 &&
         fabs(a.levelSnapPunishFac[3] - b.levelSnapPunishFac[3]) < 0.1 &&
         fabs(a.levelSnapPunishFac[4] - b.levelSnapPunishFac[4]) < 0.1 &&
         fabs(a.levelSnapPunishFac[5] - b.levelSnapPunishFac[5]) < 0.1 &&
         fabs(a.levelSnapPunishFac[6] - b.levelSnapPunishFac[6]) < 0.1 &&
         fabs(a.fullTurnAngle - b.fullTurnAngle) < 0.1 &&
         a.restrPosRestr == b.restrPosRestr &&
         a.restrNegRestr == b.restrNegRestr &&
         a.noRestrFilter == b.noRestrFilter;
}
}  // namespace osm
}  // namespace pfaedle
#endif  // PFAEDLE_OSM_OSMREADOPTS_H_

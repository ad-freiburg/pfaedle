// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_TRGRAPH_NODEPL_H_
#define PFAEDLE_TRGRAPH_NODEPL_H_

#include <map>
#include <string>
#include <unordered_map>
#include <vector>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "pfaedle/Def.h"
#include "pfaedle/trgraph/StatInfo.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"

using util::geograph::GeoNodePL;

namespace pfaedle {
namespace trgraph {

struct Component {
  float maxSpeed;
};

/*
 * A node payload class for the transit graph.
 */
class NodePL {
 public:
  NodePL();
  NodePL(const POINT& geom);  // NOLINT
  NodePL(const POINT& geom, const StatInfo& si);

  // Return the geometry of this node.
  const POINT* getGeom() const;
  void setGeom(const POINT& geom);

  // Fill obj with k/v pairs describing the parameters of this payload.
  util::json::Dict getAttrs() const;

  // Set the station info for this node
  void setSI(const StatInfo& si);

  // Return the station info for this node
  const StatInfo* getSI() const;
  StatInfo* getSI();

  // Delete the station info for this node
  void setNoStat();

  // Get the component of this node
  const Component& getComp() const;

  // Get the component of this node
  uint32_t getCompId() const;

  // Set the component of this node
  void setComp(uint32_t c);

  // Make this node a blocker
  void setBlocker();

  // Check if this node is a blocker
  bool isBlocker() const;

  // Make this node a turning cycle
  void setTurnCycle();

  // Check if this node is a blocker
  bool isTurnCycle() const;

  // Mark this node as visited (usefull for counting search space in Dijkstra)
  // (only works for DEBUG build type)
  void setVisited() const;

  static std::vector<Component> comps;

 private:
  POINT _geom;
  uint32_t _si;
  uint32_t _component;

#ifdef PFAEDLE_DBG
  mutable bool _vis;
#endif
  static std::vector<StatInfo> _statInfos;
};
}  // namespace trgraph
}  // namespace pfaedle

#endif  // PFAEDLE_TRGRAPH_NODEPL_H_

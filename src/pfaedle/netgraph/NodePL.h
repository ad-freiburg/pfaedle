// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_NETGRAPH_NODEPL_H_
#define PFAEDLE_NETGRAPH_NODEPL_H_

#include <map>
#include <string>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "util/geo/GeoGraph.h"

using util::geograph::GeoNodePL;


namespace pfaedle {
namespace netgraph {

/*
 * A payload class for edges on a network graph - that is a graph
 * that exactly represents a physical public transit network
 */
class NodePL : public GeoNodePL<PFAEDLE_PRECISION> {
 public:
  NodePL() {}
  NodePL(const POINT& geom) { _geom = geom; }  // NOLINT

  const POINT* getGeom() const { return &_geom; }
  util::json::Dict getAttrs() const {
    return util::json::Dict();
  }

 private:
  POINT _geom;
};
}  // namespace netgraph
}  // namespace pfaedle

#endif  // PFAEDLE_NETGRAPH_NODEPL_H_

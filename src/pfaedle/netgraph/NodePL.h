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
using util::geo::DPoint;

namespace pfaedle {
namespace netgraph {

/*
 * A payload class for edges on a network graph - that is a graph
 * that exactly represents a physical public transit network
 */
class NodePL : public GeoNodePL<double> {
 public:
  NodePL() {}
  NodePL(const util::geo::DPoint& geom) { _geom = geom; }  // NOLINT

  const DPoint* getGeom() const { return &_geom; }
  util::json::Dict getAttrs() const {
    return util::json::Dict();
  }

 private:
  DPoint _geom;
};
}  // namespace netgraph
}  // namespace pfaedle

#endif  // PFAEDLE_NETGRAPH_NODEPL_H_

// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_NODEPL_H_
#define PFAEDLE_ROUTER_NODEPL_H_

#include <map>
#include <string>
#include "pfaedle/trgraph/Graph.h"
#include "util/geo/GeoGraph.h"
#include "util/geo/Geo.h"
#include "pfaedle/Def.h"

using util::geograph::GeoNodePL;


namespace pfaedle {
namespace router {

class NodePL {
 public:
  NodePL() : _n(0) {}
  NodePL(const pfaedle::trgraph::Node* n) : _n(n) {}  // NOLINT

  const POINT* getGeom() const {
    return !_n ? 0 : _n->pl().getGeom();
  }
  util::json::Dict getAttrs() const {
    if (_n) return _n->pl().getAttrs();
    return util::json::Dict();
  }

 private:
  const pfaedle::trgraph::Node* _n;
};
}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_NODEPL_H_

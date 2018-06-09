// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_GRAPH_H_
#define PFAEDLE_ROUTER_GRAPH_H_

#include "pfaedle/trgraph/Graph.h"
#include "pfaedle/router/EdgePL.h"
#include "pfaedle/router/NodePL.h"
#include "util/graph/DirGraph.h"

using util::geo::Point;
using util::geo::Line;

namespace pfaedle {
namespace router {

typedef util::graph::Edge<router::NodePL, router::EdgePL> Edge;
typedef util::graph::Node<router::NodePL, router::EdgePL> Node;
typedef util::graph::DirGraph<router::NodePL, router::EdgePL> Graph;

}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_GRAPH_H_

// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_NETGRAPH_GRAPH_H_
#define PFAEDLE_NETGRAPH_GRAPH_H_

#include "pfaedle/netgraph/NodePL.h"
#include "pfaedle/netgraph/EdgePL.h"
#include "util/graph/UndirGraph.h"

using util::geo::Point;
using util::geo::Line;
using util::geo::FPoint;
using util::geo::FLine;

namespace pfaedle {
namespace netgraph {

/*
 * A payload class for edges on a network graph - that is a graph
 * that exactly represents a physical public transit network
 */
typedef util::graph::Edge<NodePL, EdgePL> Edge;
typedef util::graph::Node<NodePL, EdgePL> Node;
typedef util::graph::UndirGraph<NodePL, EdgePL> Graph;

}  // namespace netgraph
}  // namespace pfaedle

#endif  // PFAEDLE_NETGRAPH_GRAPH_H_

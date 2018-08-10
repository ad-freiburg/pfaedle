// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_TRGRAPH_GRAPH_H_
#define PFAEDLE_TRGRAPH_GRAPH_H_

#include "pfaedle/trgraph/NodePL.h"
#include "pfaedle/trgraph/EdgePL.h"
#include "util/graph/UndirGraph.h"
#include "util/graph/DirGraph.h"
#include "util/geo/Grid.h"

using util::geo::Grid;
using util::geo::Point;
using util::geo::Line;

namespace pfaedle {
namespace trgraph {

/*
 * A graph for physical transit networks
*/
typedef util::graph::Edge<NodePL, EdgePL> Edge;
typedef util::graph::Node<NodePL, EdgePL> Node;
typedef util::graph::DirGraph<NodePL, EdgePL> Graph;
typedef Grid<Node*, Point, double> NodeGrid;
typedef Grid<Edge*, Line, double> EdgeGrid;

}  // namespace trgraph
}  // namespace pfaedle

#endif  // PFAEDLE_TRGRAPH_GRAPH_H_

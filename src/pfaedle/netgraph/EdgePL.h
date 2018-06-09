// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_NETGRAPH_EDGEPL_H_
#define PFAEDLE_NETGRAPH_EDGEPL_H_

#include <set>
#include <string>
#include <map>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "util/geo/GeoGraph.h"

using util::geograph::GeoEdgePL;
using ad::cppgtfs::gtfs::Trip;

namespace pfaedle {
namespace netgraph {

/*
 * A payload class for edges on a network graph - that is a graph
 * that exactly represents a physical public transit network
 */
class EdgePL : public GeoEdgePL<float> {
 public:
  EdgePL() {}
  EdgePL(const util::geo::FLine& l, const std::set<const Trip*>& trips)
      : _l(l), _trips(trips) {}
  const util::geo::FLine* getGeom() const { return &_l; }
  void getAttrs(std::map<std::string, std::string>* obj) const {
    (*obj)["numtrips"] = std::to_string(_trips.size());
  }

 private:
  util::geo::FLine _l;
  std::set<const Trip*> _trips;
};
}  // namespace netgraph
}  // namespace pfaedle

#endif  // PFAEDLE_NETGRAPH_EDGEPL_H_

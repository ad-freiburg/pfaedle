// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_NETGRAPH_EDGEPL_H_
#define PFAEDLE_NETGRAPH_EDGEPL_H_

#include <set>
#include <string>
#include <vector>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "util/String.h"
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
      : _l(l), _trips(trips) {
    for (const auto t : _trips) {
      _routeShortNames.insert(t->getRoute()->getShortName());
      _tripShortNames.insert(t->getShortname());
    }
  }
  const util::geo::FLine* getGeom() const { return &_l; }
  util::json::Dict getAttrs() const {
    util::json::Dict obj;
    obj["num_trips"] = static_cast<int>(_trips.size());
    obj["route_short_names"] = util::json::Array(
        _routeShortNames.begin(), _routeShortNames.end());
    obj["trip_short_names"] = util::json::Array(_tripShortNames.begin(),
                                                       _tripShortNames.end());
    return obj;
  }

 private:
  util::geo::FLine _l;
  std::set<const Trip*> _trips;
  std::set<std::string> _routeShortNames;
  std::set<std::string> _tripShortNames;
};
}  // namespace netgraph
}  // namespace pfaedle

#endif  // PFAEDLE_NETGRAPH_EDGEPL_H_

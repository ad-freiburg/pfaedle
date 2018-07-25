// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_NETGRAPH_EDGEPL_H_
#define PFAEDLE_NETGRAPH_EDGEPL_H_

#include <map>
#include <set>
#include <string>
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
  void getAttrs(std::map<std::string, std::string>* obj) const {
    (*obj)["num_trips"] = std::to_string(_trips.size());
    (*obj)["route_short_names"] =
        util::implode(_routeShortNames.begin(), _routeShortNames.end(), ", ");
    (*obj)["trip_short_names"] =
        util::implode(_tripShortNames.begin(), _tripShortNames.end(), ", ");
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

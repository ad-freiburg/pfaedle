// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_ROUTINGATTRS_H_
#define PFAEDLE_ROUTER_ROUTINGATTRS_H_

#include <map>
#include <string>
#include "pfaedle/trgraph/EdgePL.h"

using pfaedle::trgraph::TransitEdgeLine;

namespace pfaedle {
namespace router {

struct RoutingAttrs {
  RoutingAttrs() : fromString(""), toString(""), shortName(""), _simiCache() {}
  std::string fromString;
  std::string toString;
  std::string shortName;

  mutable std::map<const TransitEdgeLine*, double> _simiCache;

  double simi(const TransitEdgeLine* line) const {
    auto i = _simiCache.find(line);
    if (i != _simiCache.end()) return i->second;

    double cur = 1;
    if (router::lineSimi(line->shortName, shortName) > 0.5) cur -= 0.33;

    if (line->toStr.empty() || router::statSimi(line->toStr, toString) > 0.5)
      cur -= 0.33;

    if (line->fromStr.empty() ||
        router::statSimi(line->fromStr, fromString) > 0.5)
      cur -= 0.33;

    _simiCache[line] = cur;

    return cur;
  }
};

inline bool operator==(const RoutingAttrs& a, const RoutingAttrs& b) {
  return a.shortName == b.shortName && a.toString == b.toString &&
         a.fromString == b.fromString;
}

inline bool operator!=(const RoutingAttrs& a, const RoutingAttrs& b) {
  return !(a == b);
}

inline bool operator<(const RoutingAttrs& a, const RoutingAttrs& b) {
  return a.fromString < b.fromString ||
         (a.fromString == b.fromString && a.toString < b.toString) ||
         (a.fromString == b.fromString && a.toString == b.toString &&
          a.shortName < b.shortName);
}

}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_ROUTINGATTRS_H_

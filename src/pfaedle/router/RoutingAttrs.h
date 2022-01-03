// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_ROUTINGATTRS_H_
#define PFAEDLE_ROUTER_ROUTINGATTRS_H_

#include <unordered_map>
#include <vector>
#include <string>
#include "pfaedle/statsimi-classifier/StatsimiClassifier.h"
#include "pfaedle/trgraph/EdgePL.h"

using pfaedle::trgraph::TransitEdgeLine;

namespace pfaedle {
namespace router {

struct LineSimilarity {
  bool nameSimilar : 1;
  bool fromSimilar : 1;
  bool toSimilar : 1;
};

inline bool operator<(const LineSimilarity& a, const LineSimilarity& b) {
  return (a.nameSimilar + a.fromSimilar + a.toSimilar) <
         (b.nameSimilar + b.fromSimilar + b.toSimilar);
}

struct RoutingAttrs {
  RoutingAttrs()
      : lineFrom(""), lineTo(), shortName(""), classifier(0), _simiCache() {}
  std::string lineFrom;
  std::vector<std::string> lineTo;
  std::string shortName;

  const pfaedle::statsimiclassifier::StatsimiClassifier* classifier;

  mutable std::unordered_map<const TransitEdgeLine*, LineSimilarity> _simiCache;

  LineSimilarity simi(const TransitEdgeLine* line) const {
    // shortcut, if we don't have a line information, classify as similar
    if (line->shortName.empty() && line->toStr.empty() && line->fromStr.empty())
      return {true, true, true};

    auto i = _simiCache.find(line);
    if (i != _simiCache.end()) return i->second;

    LineSimilarity ret{false, false, false};

    if (shortName.empty() || router::lineSimi(line->shortName, shortName) > 0.5)
      ret.nameSimilar = true;

    if (lineTo.size() == 0) {
      ret.toSimilar = true;
    } else {
      for (const auto& lTo : lineTo) {
        if (lTo.empty() || classifier->similar(line->toStr, lTo)) {
          ret.toSimilar = true;
          break;
        }
      }
    }

    if (lineFrom.empty() || classifier->similar(line->fromStr, lineFrom))
      ret.fromSimilar = true;

    _simiCache[line] = ret;

    return ret;
  }

  void merge(const RoutingAttrs& other) {
    assert(other.lineFrom == lineFrom);
    assert(other.shortName == shortName);

    for (const auto& l : other.lineTo) {
      auto i = std::lower_bound(lineTo.begin(), lineTo.end(), l);
      if (i != lineTo.end() && (*i) == l) continue;  // already present
      lineTo.insert(i, l);
    }
  }
};

inline bool operator==(const RoutingAttrs& a, const RoutingAttrs& b) {
  return a.shortName == b.shortName && a.lineFrom == b.lineFrom;
}

inline bool operator!=(const RoutingAttrs& a, const RoutingAttrs& b) {
  return !(a == b);
}

inline bool operator<(const RoutingAttrs& a, const RoutingAttrs& b) {
  return a.lineFrom < b.lineFrom ||
         (a.lineFrom == b.lineFrom && a.shortName < b.shortName);
}

}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_ROUTINGATTRS_H_

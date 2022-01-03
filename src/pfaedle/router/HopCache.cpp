// Copyright 2020, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <utility>
#include <set>
#include "pfaedle/router/HopCache.h"
#include "pfaedle/trgraph/Graph.h"
#include "util/Misc.h"

using pfaedle::router::HopCache;
using pfaedle::trgraph::Edge;

// _____________________________________________________________________________
void HopCache::setMin(const Edge* a, const Edge* b, uint32_t val) {
  _cache.set(a, b, val);
}

// _____________________________________________________________________________
void HopCache::setEx(const Edge* a, const Edge* b, uint32_t val) {
  int64_t v = val;
  _cache.set(a, b, -(v + 1));
}

// _____________________________________________________________________________
void HopCache::setMin(const Edge* a, const std::set<Edge*>& b, uint32_t val) {
  for (auto eb : b) _cache.set(a, eb, val);
}

// _____________________________________________________________________________
void HopCache::setMin(const std::set<Edge*>& a, const Edge* b, uint32_t val) {
  for (auto ea : a) _cache.set(ea, b, val);
}

// _____________________________________________________________________________
std::pair<uint32_t, bool> HopCache::get(const Edge* a, const Edge* b) const {
  int64_t v = _cache.get(a, b);
  if (v < 0) return {(-v) - 1, 1};
  return {v, 0};
}

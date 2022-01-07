// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_STATS_H_
#define PFAEDLE_ROUTER_STATS_H_

#include <algorithm>
#include <iostream>
#include <string>
#include "util/String.h"

namespace pfaedle {
namespace router {

struct Stats {
  Stats()
      : totNumTrips(0),
        numTries(0),
        numTrieLeafs(0),
        solveTime(0),
        dijkstraIters(0) {}
  size_t totNumTrips;
  size_t numTries;
  size_t numTrieLeafs;
  double solveTime;
  size_t dijkstraIters;
};

inline Stats operator+ (const Stats& c1, const Stats& c2) {
  Stats ret = c1;
  ret.totNumTrips += c2.totNumTrips;
  ret.numTries += c2.numTries;
  ret.numTrieLeafs += c2.numTrieLeafs;
  ret.solveTime += c2.solveTime;
  ret.dijkstraIters += c2.dijkstraIters;
  return ret;
}

inline Stats& operator+= (Stats& c1, const Stats& c2) {
  c1 = c1 + c2;
  return c1;
}

}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_STATS_H_

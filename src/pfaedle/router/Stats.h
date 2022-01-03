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

}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_STATS_H_

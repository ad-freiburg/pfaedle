// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_COMP_H_
#define PFAEDLE_ROUTER_COMP_H_

#include <iostream>
#include <algorithm>
#include <string>
#include "util/String.h"

namespace pfaedle {
namespace router {

using util::editDist;

// _____________________________________________________________________________
inline double statSimi(const std::string& a, const std::string& b) {
  if (a == b) return 1;

  if (a.empty() || b.empty()) return 0;

  if (a.size() > b.size() + 1) {
    // check if a begins with b
    if (a.compare(0, b.size() + 1, b + " ") == 0) {
      return 1;
    }

    // check if a ends with b
    if (a.compare(a.size() - (b.size() + 1), b.size() + 1, " " + b) == 0) {
      return 1;
    }
  }

  if (b.size() > a.size() + 1) {
    // check if b begins with a
    if (b.compare(0, a.size() + 1, a + " ") == 0) {
      return 1;
    }

    // check if b ends with a
    if (b.compare(b.size() - (a.size() + 1), a.size() + 1, " " + a) == 0) {
      return 1;
    }
  }

  if (static_cast<double>(editDist(a, b)) /
          (std::max(static_cast<double>(a.size()),
                    static_cast<double>(b.size()))) <
      0.05)
    return 1;

  return 0;
}

// _____________________________________________________________________________
inline double lineSimi(const std::string& a, const std::string& b) {
  if (a == b) return 1;

  if (a.empty() || b.empty()) return 0;

  // if one of the lines is completely contained in the other, return 1
  if (a.find(b) != std::string::npos) {
    return 1;
  } else if (b.find(a) != std::string::npos) {
    return 1;
  }

  return 0;
}
}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_COMP_H_

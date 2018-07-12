// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_LINE_H_
#define UTIL_GEO_LINE_H_

#include <vector>
#include "./Point.h"

namespace util {
namespace geo {

template<typename T>
using Line = std::vector<Point<T>>;

template<typename T>
using LineSegment = std::pair<Point<T>, Point<T>>;


}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_LINE_H_

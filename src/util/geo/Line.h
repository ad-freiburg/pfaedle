// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_LINE_H_
#define UTIL_GEO_LINE_H_

#include <vector>
#include "./Point.h"

namespace util {
namespace geon {

template<typename T>
using Line = std::vector<Point<T>>;


}  // namespace geon
}  // namespace util

#endif  // UTIL_GEO_LINE_H_

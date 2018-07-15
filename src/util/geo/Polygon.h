// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_POLYGON_H_
#define UTIL_GEO_POLYGON_H_

#include <vector>
#include "./Point.h"
#include "./Line.h"

namespace util {
namespace geo {

template <typename T>
class Polygon {
 public:
  // maximum inverse box as default value of box
  Polygon(const Line<T>& l) : _outer(l) {}

  const std::vector<Point<T>>& getOuter() const { return _outer; }
  std::vector<Point<T>>& getOuter() {return _outer; }

 private:
  std::vector<Point<T>> _outer;
};


}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_LINE_H_

// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_POLYGON_H_
#define UTIL_GEO_POLYGON_H_

#include <vector>
#include "./Box.h"
#include "./Line.h"
#include "./Point.h"

namespace util {
namespace geo {

template <typename T>
class Polygon {
 public:
  Polygon() {}

  Polygon(const Line<T>& l) : _outer(l) {}
  Polygon(const Box<T>& b)
      : _outer({b.getLowerLeft(),
                Point<T>(b.getUpperRight().getX(), b.getLowerLeft().getY()),
                b.getUpperRight(),
                Point<T>(b.getLowerLeft().getX(), b.getUpperRight().getY())}) {}

  const std::vector<Point<T>>& getOuter() const { return _outer; }
  std::vector<Point<T>>& getOuter() { return _outer; }

 private:
  std::vector<Point<T>> _outer;
};

}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_LINE_H_

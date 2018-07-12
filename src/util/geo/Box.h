// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Author: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef UTIL_GEO_BOX_H_
#define UTIL_GEO_BOX_H_

#include "./Point.h"

namespace util {
namespace geo {

template <typename T>
class Box {
 public:
  // maximum inverse box as default value of box
  Box()
      : _ll(std::numeric_limits<T>::max(), std::numeric_limits<T>::max()),
        _ur(std::numeric_limits<T>::min(), std::numeric_limits<T>::min()) {}
  Box(const Point<T>& ll, const Point<T>& ur) : _ll(ll), _ur(ur) {}
  const Point<T>& getLowerLeft() const { return _ll; }
  const Point<T>& getUpperRight() const { return _ur; }

  Point<T>& getLowerLeft() { return _ll; }
  Point<T>& getUpperRight() { return _ur; }

  void setLowerLeft(const Point<T>& ll) { _ll = ll; }
  void setUpperRight(const Point<T>& ur) { _ur = ur; }

  bool operator==(const Box<T>& b) const {
    return getLowerLeft() == b.getLowerLeft() &&
           getUpperRight == b.getUpperRight();
  }

  bool operator!=(const Box<T>& p) const { return !(*this == p); }

 private:
  Point<T> _ll, _ur;
};

}  // namespace geo
}  // namespace util

#endif  // UTIL_GEO_BOX_H_

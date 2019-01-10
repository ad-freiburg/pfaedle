// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_EVAL_RESULT_H_
#define PFAEDLE_EVAL_RESULT_H_

#include "pfaedle/gtfs/Feed.h"
#include "ad/cppgtfs/gtfs/Feed.h"

using pfaedle::gtfs::Trip;
using ad::cppgtfs::gtfs::Shape;

namespace pfaedle {
namespace eval {

/*
 * A single evaluation result.
 */
class Result {
 public:
  Result(const Trip* t, double dist) : _t(t), _dist(dist) {}

  double getDist() const { return _dist; }
  const Trip* getTrip() const { return _t; }

 private:
  const Trip* _t;
  double _dist;
};

inline bool operator<(const Result& lhs, const Result& rhs) {
  return lhs.getDist() < rhs.getDist() ||
         (lhs.getDist() == rhs.getDist() && lhs.getTrip() < rhs.getTrip());
}

}  // namespace eval
}  // namespace pfaedle

#endif  // PFAEDLE_EVAL_RESULT_H_

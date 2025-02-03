// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_EVAL_COLLECTOR_H_
#define PFAEDLE_EVAL_COLLECTOR_H_

#include <fstream>
#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "pfaedle/Def.h"
#include "shapevl/Result.h"
#include "util/geo/Geo.h"
#include "util/json/Writer.h"

using ad::cppgtfs::gtfs::Shape;
using ad::cppgtfs::gtfs::Trip;

namespace pfaedle {
namespace eval {

struct lineCmp {
  bool operator()(const LINE& a, const LINE& b) const {
    if (a.size() != b.size()) {
      return a.size() < b.size();
    }

    for (size_t i = 0; i < a.size(); i++) {
      if (util::geo::dist(a[i], b[i]) > .00001) {
        return (a[i].getX() < b[i].getX()) ||
               (a[i].getX() == b[i].getX() && a[i].getY() < b[i].getY());
        ;
      }
    }

    return false;
  }
};

/*
 * Collects routing results for evaluation
 */
class Collector {
 public:
  Collector(std::ostream* reportOut, int reportLevel)
      : _trips(0),
        _noOrigShp(0),
        _accFdSum(0),
        _unmatchedSegSum(0),
        _unmatchedSegLengthSum(0),
        _an0(0),
        _an5(0),
        _an10(0),
        _an30(0),
        _an50(0),
        _an70(0),
        _an90(0),
        _reportOut(reportOut), _reportLevel(reportLevel) {}

  // Add a shape found by our tool newS for a trip t with newly calculated
  // station dist values with the old shape oldS
  double add(const Trip* oldT, const Shape* oldS, size_t numOldTrips,
             const Trip* newT, const Shape* newS, double segLen);

  // Return the set of all Result objects
  const std::set<Result>& getResults() const;

  // Print general stats to os
  void printStats(std::ostream* os) const;

  // Print general stats to os
  void printShortStats(std::ostream* os) const;

  // Get JSON stats
  std::map<string, double> getStats();

  // Print a CSV for the results to os
  void printCsv(std::ostream* os, const std::set<Result>& result) const;

  // Return the averaged average frechet distance
  double getAvgDist() const;

  static LINE getLine(const Shape* s, std::vector<double>* dists);

  double getAcc() const;

 private:
  std::set<Result> _results;
  std::map<LINE, std::map<LINE, double, lineCmp>, lineCmp> _accFdCache;
  std::map<LINE, std::map<LINE, double, lineCmp>, lineCmp> _dACache;
  std::map<LINE, std::map<LINE, double, lineCmp>, lineCmp> _fdCache;
  std::map<LINE, std::map<LINE, double, lineCmp>, lineCmp> _dCache;
  std::map<LINE, std::map<LINE, double, lineCmp>, lineCmp> _lenDiffCache;

  size_t _trips;
  size_t _noOrigShp;

  std::vector<double> _distDiffs;
  std::vector<double> _hopDists;

  double _accFdSum;
  size_t _unmatchedSegSum;
  double _unmatchedSegLengthSum;

  size_t _an0 = 0;
  size_t _an5 = 0;
  size_t _an10 = 0;
  size_t _an20 = 0;
  size_t _an30 = 0;
  size_t _an50 = 0;
  size_t _an70 = 0;
  size_t _an90 = 0;

  std::ostream* _reportOut;
  int _reportLevel = 1;

  std::pair<size_t, double> getDa(const std::vector<LINE>& a,
                                  const std::vector<LINE>& b, double segLen);

  static std::vector<LINE> segmentize(
      const Trip* t, const LINE& shape, const std::vector<double>& dists,
      std::vector<std::pair<double, double>>& lenDist);
};

}  // namespace eval
}  // namespace pfaedle

#endif  // PFAEDLE_EVAL_COLLECTOR_H_

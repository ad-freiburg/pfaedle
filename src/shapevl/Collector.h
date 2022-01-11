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

/*
 * Collects routing results for evaluation
 */
class Collector {
 public:
  Collector(std::ostream* reportOut)
      : _trips(0),
        _noOrigShp(0),
        _fdSum(0),
        _unmatchedSegSum(0),
        _unmatchedSegLengthSum(0),
        _an0(0),
        _an5(0),
        _an10(0),
        _an30(0),
        _an50(0),
        _an70(0),
        _an90(0),
        _reportOut(reportOut) {}

  // Add a shape found by our tool newS for a trip t with newly calculated
  // station dist values with the old shape oldS
  double add(const Trip* oldT, const Shape* oldS, const Trip* newT,
             const Shape* newS);

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

  static LINE getWebMercLine(const Shape* s, std::vector<double>* dists);

  double getAcc() const;

 private:
  std::set<Result> _results;
  std::map<const Shape*, std::map<std::string, double> > _dCache;
  std::map<const Shape*, std::map<std::string, std::pair<size_t, double> > >
      _dACache;
  size_t _trips;
  size_t _noOrigShp;

  std::vector<double> _distDiffs;

  double _fdSum;
  size_t _unmatchedSegSum;
  double _unmatchedSegLengthSum;

  size_t _an0;
  size_t _an5;
  size_t _an10;
  size_t _an30;
  size_t _an50;
  size_t _an70;
  size_t _an90;

  std::ostream* _reportOut;

  static std::pair<size_t, double> getDa(const std::vector<LINE>& a,
                                         const std::vector<LINE>& b);

  static std::vector<LINE> segmentize(const Trip* t, const LINE& shape,
                                      const std::vector<double>& dists,
                                      std::vector<std::pair<double, double>>& lenDist);

  static std::vector<double> getBins(double mind, double maxd, size_t steps);
};

}  // namespace eval
}  // namespace pfaedle

#endif  // PFAEDLE_EVAL_COLLECTOR_H_

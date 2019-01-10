// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_EVAL_COLLECTOR_H_
#define PFAEDLE_EVAL_COLLECTOR_H_

#include <map>
#include <ostream>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "pfaedle/gtfs/Feed.h"
#include "pfaedle/Def.h"
#include "pfaedle/eval/Result.h"
#include "util/geo/Geo.h"

using pfaedle::gtfs::Trip;
using ad::cppgtfs::gtfs::Shape;

namespace pfaedle {
namespace eval {

/*
 * Collects routing results for evaluation
 */
class Collector {
 public:
  Collector(const std::string& evalOutPath, const std::vector<double>& dfBins)
      : _noOrigShp(0),
        _fdSum(0),
        _unmatchedSegSum(0),
        _unmatchedSegLengthSum(0),
        _evalOutPath(evalOutPath),
        _dfBins(dfBins) {}

  // Add a shape found by our tool newS for a trip t with newly calculated
  // station dist values with the old shape oldS
  double add(const Trip* t, const Shape* oldS, const Shape& newS,
             const std::vector<double>& newDists);

  // Return the set of all Result objects
  const std::set<Result>& getResults() const;

  // Print general stats to os
  void printStats(std::ostream* os) const;

  // Print histogramgs for the results to os
  void printHisto(std::ostream* os, const std::set<Result>& result,
                  const std::vector<double>& bins) const;

  // Print a CSV for the results to os
  void printCsv(std::ostream* os, const std::set<Result>& result,
                const std::vector<double>& bins) const;

  // Return the averaged average frechet distance
  double getAvgDist() const;

  static LINE getWebMercLine(const Shape* s, double from, double to);
  static LINE getWebMercLine(const Shape* s, double from, double to,
                             std::vector<double>* dists);

 private:
  std::set<Result> _results;
  std::set<Result> _resultsAN;
  std::set<Result> _resultsAL;
  std::map<const Shape*, std::map<std::string, double> > _dCache;
  std::map<const Shape*, std::map<std::string, std::pair<size_t, double> > >
      _dACache;
  size_t _noOrigShp;

  double _fdSum;
  size_t _unmatchedSegSum;
  double _unmatchedSegLengthSum;

  std::string _evalOutPath;

  std::vector<double> _dfBins;

  static std::pair<size_t, double> getDa(const std::vector<LINE>& a,
                                         const std::vector<LINE>& b);

  static std::vector<LINE> segmentize(const Trip* t, const LINE& shape,
                                      const std::vector<double>& dists,
                                      const std::vector<double>* newTripDists);

  static std::vector<double> getBins(double mind, double maxd, size_t steps);
};

}  // namespace eval
}  // namespace pfaedle

#endif  // PFAEDLE_EVAL_COLLECTOR_H_

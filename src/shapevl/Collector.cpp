// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <fstream>
#include <set>
#include <string>
#include <utility>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "pfaedle/Def.h"
#include "shapevl/Collector.h"
#include "shapevl/Result.h"
#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"
#include "util/geo/output/GeoJsonOutput.h"
#include "util/log/Log.h"

using util::geo::PolyLine;

using ad::cppgtfs::gtfs::Shape;
using ad::cppgtfs::gtfs::Trip;
using pfaedle::eval::Collector;
using pfaedle::eval::Result;
using util::geo::output::GeoJsonOutput;

// _____________________________________________________________________________
double Collector::add(const Trip* oldT, const Shape* oldS, const Trip* newT,
                      const Shape* newS) {
  // This adds a new trip with a new shape to our evaluation.

  _trips++;

  if (!oldS) {
    // If there is no original shape, we cannot compare them - abort!
    _noOrigShp++;
    return 0;
  }

  for (auto st : oldT->getStopTimes()) {
    if (st.getShapeDistanceTravelled() < 0) {
      // we cannot safely compare trips without shape dist travelled
      // information - abort!
      _noOrigShp++;
      return 0;
    }
  }

  for (auto st : newT->getStopTimes()) {
    if (st.getShapeDistanceTravelled() < 0) {
      // we cannot safely compare trips without shape dist travelled
      // information - abort!
      _noOrigShp++;
      return 0;
    }
  }

  double fd = 0;

  // A "segment" is a path from station s_i to station s_{i+1}

  size_t unmatchedSegments;        // number of unmatched segments
  double unmatchedSegmentsLength;  // total _an. length of unmatched segments

  std::vector<double> oldDists;
  LINE oldL = getWebMercLine(oldS, &oldDists);

  std::vector<double> newDists;
  LINE newL = getWebMercLine(newS, &newDists);

  auto oldSegs = segmentize(oldT, oldL, oldDists);
  auto newSegs = segmentize(newT, newL, newDists);

  // new lines build from cleaned-up shapes
  LINE oldLCut;
  LINE newLCut;

  for (auto oldL : oldSegs)
    oldLCut.insert(oldLCut.end(), oldL.begin(), oldL.end());

  for (auto newL : newSegs)
    newLCut.insert(newLCut.end(), newL.begin(), newL.end());

  // determine the scale factor between the distance in projected
  // coordinates and the real-world distance in meters
  auto avgY =
      (oldSegs.front().front().getY() + oldSegs.back().back().getY()) / 2;
  double fac = cos(2 * atan(exp(avgY / 6378137.0)) - 1.5707965);

  double SEGL = 10;

  if (_dCache.count(oldS) && _dCache.find(oldS)->second.count(newS->getId())) {
    fd = _dCache[oldS][newS->getId()];
  } else {
    fd = util::geo::accFrechetDistC(oldLCut, newLCut, SEGL / fac) * fac;
    _dCache[oldS][newS->getId()] = fd;
  }

  if (_dACache.count(oldS) &&
      _dACache.find(oldS)->second.count(newS->getId())) {
    unmatchedSegments = _dACache[oldS][newS->getId()].first;
    unmatchedSegmentsLength = _dACache[oldS][newS->getId()].second;
  } else {
    auto dA = getDa(oldSegs, newSegs);
    _dACache[oldS][newS->getId()] = dA;
    unmatchedSegments = dA.first;
    unmatchedSegmentsLength = dA.second;
  }

  double totL = 0;
  for (auto l : oldSegs) totL += util::geo::len(l) * fac;

  // filter out shapes with a length of under 5 meters - they are most likely
  // artifacts
  if (totL < 5) {
    _noOrigShp++;
    return 0;
  }

  _fdSum += fd / totL;
  _unmatchedSegSum += unmatchedSegments;
  _unmatchedSegLengthSum += unmatchedSegmentsLength;

  double avgFd = fd / totL;
  double AN = static_cast<double>(unmatchedSegments) /
              static_cast<double>(oldSegs.size());
  double AL = unmatchedSegmentsLength / totL;

  _results.insert(Result(oldT, avgFd));

  if (AN <= 0.0001) _an0++;
  if (AN <= 0.05) _an5++;
  if (AN <= 0.1) _an10++;
  if (AN <= 0.3) _an30++;
  if (AN <= 0.5) _an50++;
  if (AN <= 0.7) _an70++;
  if (AN <= 0.9) _an90++;

  LOG(VDEBUG) << "This result (" << oldT->getId()
              << "): A_N/N = " << unmatchedSegments << "/" << oldSegs.size()
              << " = " << AN << " A_L/L = " << unmatchedSegmentsLength << "/"
              << totL << " = " << AL << " d_f = " << avgFd;

  if (_reportOut) {
    (*_reportOut) << oldT->getId() << "\t" << AN << "\t" << AL << "\t" << avgFd
                  << "\t" << util::geo::getWKT(oldSegs) << "\t"
                  << util::geo::getWKT(newSegs) << "\n";
  }

  return avgFd;
}

// _____________________________________________________________________________
std::vector<LINE> Collector::segmentize(const Trip* t, const LINE& shape,
                                        const std::vector<double>& dists) {
  // The straightforward way to segmentize the shape would be to just cut it at
  // the exact measurements in stop_times.txt. We have tried that, but found
  // that it produces misleading results for the following reason:
  //
  // 1) The measurement specifies an exact position on the shape.
  // 2) Even if we consider correct rail or bus tracks, the "right" position
  //    where a vehicle may come to a halt is not a point - its a line segment,
  //    basically the entire track in railroad term
  // 3) The position point on the shape in real-world feeds may be either a) the
  //    position where a train comes to a halt, b) the position where a carriage
  //    comes to a halt, c) the beginning of the tracks line segment, d) the end
  //    of the tracks line segment, e) the center of the tracks line segment, f)
  //    ... (any position on the tracks line segment.
  // 4) The "correct" position is NOT well-defined.
  // 5) As tracks are often longer than 20 meters, this will dillute our AN
  //    measure, although the shape is CORRECT (because the ground truth uses
  //    a different position philosophy than the test data)
  // 6) To normalize this, we use the following approach:
  //      a) Get the exact progression of the measurment on the shape
  //      b) Extract a segment of 200 meters, with the measurement progress in
  //      the middle c) Project the GROUND TRUTH station coordinate to this
  //      segment d) The result is the cutting point
  // 7) If a completely wrong track was chosen, the frechet distance will still
  //    be greater than 20 meters and AN will measure an unmatch.
  // 8) TODO: implement this, explain this in diss
  std::vector<LINE> ret;

  if (t->getStopTimes().size() < 2) return ret;

  POLYLINE pl(shape);
  std::vector<std::pair<POINT, double> > cuts;

  size_t i = 0;
  for (auto st : t->getStopTimes()) {
    cuts.push_back(std::pair<POINT, double>(
        util::geo::latLngToWebMerc<PFDL_PREC>(st.getStop()->getLat(),
                                              st.getStop()->getLng()),
        st.getShapeDistanceTravelled()));
    i++;
  }

  // get first half of geometry, and search for start point there!
  size_t before = std::upper_bound(dists.begin(), dists.end(), cuts[1].second) -
                  dists.begin();
  if (before + 1 > shape.size()) before = shape.size() - 1;
  assert(shape.begin() + before + 1 <= shape.end());
  POLYLINE l(LINE(shape.begin(), shape.begin() + before + 1));
  auto lastLp = l.projectOn(cuts.front().first);

  for (size_t i = 1; i < cuts.size(); i++) {
    size_t before = shape.size();
    if (i < cuts.size() - 1 && cuts[i + 1].second > -0.5) {
      before =
          std::upper_bound(dists.begin(), dists.end(), cuts[i + 1].second) -
          dists.begin();
    }

    POLYLINE beforePl(LINE(shape.begin(), shape.begin() + before));

    auto curLp = beforePl.projectOnAfter(cuts[i].first, lastLp.lastIndex);

    ret.push_back(pl.getSegment(lastLp, curLp).getLine());
    lastLp = curLp;
  }

  return ret;
}

// _____________________________________________________________________________
LINE Collector::getWebMercLine(const Shape* s, std::vector<double>* dists) {
  LINE ret;

  for (size_t i = 0; i < s->getPoints().size(); i++) {
    ret.push_back(util::geo::latLngToWebMerc<PFDL_PREC>(s->getPoints()[i].lat,
                                                        s->getPoints()[i].lng));
    (*dists).push_back(s->getPoints()[i].travelDist);
  }

  return ret;
}

// _____________________________________________________________________________
const std::set<Result>& Collector::getResults() const { return _results; }

// _____________________________________________________________________________
double Collector::getAvgDist() const { return _fdSum / _results.size(); }

// _____________________________________________________________________________
std::vector<double> Collector::getBins(double mind, double maxd, size_t steps) {
  double bin = (maxd - mind) / steps;
  double curE = mind + bin;

  std::vector<double> ret;
  while (curE <= maxd) {
    ret.push_back(curE);
    curE += bin;
  }
  return ret;
}

// _____________________________________________________________________________
void Collector::printCsv(std::ostream* os,
                         const std::set<Result>& result) const {
  for (auto r : result) (*os) << r.getDist() << "\n";
}

// _____________________________________________________________________________
double Collector::getAcc() const {
  return static_cast<double>(_an0) / static_cast<double>(_results.size());
}

// _____________________________________________________________________________
void Collector::printShortStats(std::ostream* os) const {
  if (_results.size()) {
    (*os) << (static_cast<double>(_an0) /
              static_cast<double>(_results.size())) *
                 100
          << ",";
    (*os) << (static_cast<double>(_an5) /
              static_cast<double>(_results.size())) *
                 100
          << ",";
    (*os) << (static_cast<double>(_an10) /
              static_cast<double>(_results.size())) *
                 100
          << ",";
    (*os) << (static_cast<double>(_an30) /
              static_cast<double>(_results.size())) *
                 100
          << ",";
    (*os) << (static_cast<double>(_an50) /
              static_cast<double>(_results.size())) *
                 100
          << ",";
    (*os) << (static_cast<double>(_an70) /
              static_cast<double>(_results.size())) *
                 100
          << ",";
    (*os) << (static_cast<double>(_an90) /
              static_cast<double>(_results.size())) *
                 100;
  }
}

// _____________________________________________________________________________
void Collector::printStats(std::ostream* os) const {
  (*os) << std::setfill(' ') << std::setw(50) << "  # of trips: " << _trips
        << "\n";
  (*os) << std::setfill(' ') << std::setw(50)
        << "  # of trips new shapes were matched for: " << _results.size()
        << "\n";
  (*os) << std::setw(50) << "  # of trips without input shapes: " << _noOrigShp
        << "\n";

  if (_results.size()) {
    (*os) << std::setw(50) << "  highest avg frechet distance to input shapes: "
          << (--_results.end())->getDist() << " (on trip #"
          << (--_results.end())->getTrip()->getId() << ")\n";
    (*os) << std::setw(50) << "  lowest distance to input shapes: "
          << (_results.begin())->getDist() << " (on trip #"
          << (_results.begin())->getTrip()->getId() << ")\n";
    (*os) << std::setw(50)
          << "  averaged avg frechet distance: " << getAvgDist() << "\n";

    (*os) << "\n";
    (*os) << "  an-0: "
          << (static_cast<double>(_an0) /
              static_cast<double>(_results.size())) *
                 100
          << " %"
          << "\n";
    (*os) << "  an-5: "
          << (static_cast<double>(_an5) /
              static_cast<double>(_results.size())) *
                 100
          << " %"
          << "\n";
    (*os) << "  an-10: "
          << (static_cast<double>(_an10) /
              static_cast<double>(_results.size())) *
                 100
          << " %"
          << "\n";
    (*os) << "  acc-30: "
          << (static_cast<double>(_an30) /
              static_cast<double>(_results.size())) *
                 100
          << " %"
          << "\n";
    (*os) << "  acc-50: "
          << (static_cast<double>(_an50) /
              static_cast<double>(_results.size())) *
                 100
          << " %"
          << "\n";
    (*os) << "  acc-70: "
          << (static_cast<double>(_an70) /
              static_cast<double>(_results.size())) *
                 100
          << " %"
          << "\n";
    (*os) << "  acc-90: "
          << (static_cast<double>(_an90) /
              static_cast<double>(_results.size())) *
                 100
          << " %"
          << "\n";
  }

  (*os) << std::endl;
}

// _____________________________________________________________________________
std::pair<size_t, double> Collector::getDa(const std::vector<LINE>& a,
                                           const std::vector<LINE>& b) {
  assert(a.size() == b.size());
  std::pair<size_t, double> ret{0, 0};

  // euclidean distance on web mercator is in meters on equator,
  // and proportional to cos(lat) in both y directions
  double fac =
      cos(2 * atan(exp((a.front().front().getY() + a.back().back().getY()) /
                       6378137.0)) -
          1.5707965);

  for (size_t i = 0; i < a.size(); i++) {
    double fd = util::geo::frechetDist(a[i], b[i], 3 / fac) * fac;
    if (fd >= 20) {
      ret.first++;
      ret.second += util::geo::len(a[i]) * fac;
    }
  }

  return ret;
}

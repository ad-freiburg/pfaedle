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

  std::vector<std::pair<double, double>> newLenDists;
  std::vector<std::pair<double, double>> oldLenDists;

  auto oldSegs = segmentize(oldT, oldL, oldDists, newLenDists);
  auto newSegs = segmentize(newT, newL, newDists, oldLenDists);

  for (const auto& p : oldLenDists) {
    _distDiffs.push_back(fabs(p.first - p.second));
  }

  // new lines build from cleaned-up shapes
  LINE oldLCut;
  LINE newLCut;

  for (auto oldL : oldSegs)
    oldLCut.insert(oldLCut.end(), oldL.begin(), oldL.end());

  for (auto newL : newSegs) {
    newLCut.insert(newLCut.end(), newL.begin(), newL.end());
  }

  // determine the scale factor between the distance in projected
  // coordinates and the real-world distance in meters
  auto avgY =
      (oldSegs.front().front().getY() + oldSegs.back().back().getY()) / 2;
  double fac = cos(2 * atan(exp(avgY / 6378137.0)) - 1.5707965);

  double SEGL = 10;

  auto old = _dCache.find(oldLCut);
  if (old != _dCache.end()) {
    auto match = old->second.find(newLCut);
    if (match != old->second.end()) {
      fd = match->second;
    } else {
      fd = util::geo::accFrechetDistC(oldLCut, newLCut, SEGL / fac) * fac;
      _dCache[oldLCut][newLCut] = fd;
    }
  } else {
    fd = util::geo::accFrechetDistC(oldLCut, newLCut, SEGL / fac) * fac;
    _dCache[oldLCut][newLCut] = fd;
  }

  auto dA = getDa(oldSegs, newSegs);
  unmatchedSegments = dA.first;
  unmatchedSegmentsLength = dA.second;

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
std::vector<LINE> Collector::segmentize(
    const Trip* t, const LINE& shape, const std::vector<double>& dists,
    std::vector<std::pair<double, double>>& lenDist) {
  std::vector<LINE> ret;

  if (t->getStopTimes().size() < 2) return ret;

  POLYLINE pl(shape);
  std::vector<std::pair<POINT, double>> cuts;

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

    POLYLINE afterPl(LINE(shape.begin(), shape.begin() + before));

    auto curLp = afterPl.projectOnAfter(cuts[i].first, lastLp.lastIndex);

    auto curL = pl.getSegment(lastLp, curLp).getLine();

    double dist =
        util::geo::haversine(t->getStopTimes()[i - 1].getStop()->getLat(),
                             t->getStopTimes()[i - 1].getStop()->getLng(),
                             t->getStopTimes()[i].getStop()->getLat(),
                             t->getStopTimes()[i].getStop()->getLng());
    double len = util::geo::webMercLen(curL);
    lenDist.push_back({dist, len});

    ret.push_back(curL);
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
std::map<string, double> Collector::getStats() {
  std::map<string, double> stats;

  if (_distDiffs.size()) {
    auto i = _distDiffs.begin() + _distDiffs.size() / 2;
    std::nth_element(_distDiffs.begin(), i, _distDiffs.end());

    stats["median-dist-diff"] = *i;
  } else {
    stats["median-dist-diff"] = -1;
  }

  stats["num-trips"] = _trips;
  stats["num-trips-matched"] = _results.size();
  stats["num-trips-wo-shapes"] = _noOrigShp;
  stats["avg-fr"] = getAvgDist();
  if (_results.size()) {
    stats["max-avg-frech-dist"] = (--_results.end())->getDist();
  } else {
    stats["max-avg-frech-dist"] = -1;
  }
  stats["an-0"] =
      (static_cast<double>(_an0) / static_cast<double>(_results.size())) * 100;
  stats["an-5"] =
      (static_cast<double>(_an5) / static_cast<double>(_results.size())) * 100;
  stats["an-10"] =
      (static_cast<double>(_an10) / static_cast<double>(_results.size())) * 100;
  stats["an-30"] =
      (static_cast<double>(_an30) / static_cast<double>(_results.size())) * 100;
  stats["an-50"] =
      (static_cast<double>(_an50) / static_cast<double>(_results.size())) * 100;
  stats["an-70"] =
      (static_cast<double>(_an70) / static_cast<double>(_results.size())) * 100;
  stats["an-90"] =
      (static_cast<double>(_an90) / static_cast<double>(_results.size())) * 100;

  return stats;
}

// _____________________________________________________________________________
std::pair<size_t, double> Collector::getDa(const std::vector<LINE>& a,
                                           const std::vector<LINE>& b) {
  assert(a.size() == b.size());
  std::pair<size_t, double> ret{0, 0};

  // euclidean distance on web mercator is in meters on equator,
  // and proportional to cos(lat) in both y directions

  double fac = webMercDistFactor(a.front().front());

  double SEGL = 10;

  for (size_t i = 0; i < a.size(); i++) {
    double fd = 0;
    auto old = _dACache.find(a[i]);
    if (old != _dACache.end()) {
      auto match = old->second.find(b[i]);
      if (match != old->second.end()) {
        fd = match->second;
      } else {
        fd = util::geo::frechetDist(a[i], b[i], SEGL / fac) * fac;
        _dACache[a[i]][b[i]] = fd;
      }
    } else {
        fd = util::geo::frechetDist(a[i], b[i], SEGL / fac) * fac;
        _dACache[a[i]][b[i]] = fd;
    }

    if (fd >= 50) {
      ret.first++;
      ret.second += util::geo::webMercLen(a[i]) * 100;
    }
  }

  return ret;
}

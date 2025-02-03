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
double Collector::add(const Trip* oldT, const Shape* oldS, size_t numOldTrips,
                      const Trip* newT, const Shape* newS, double segLen) {
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

  double accFd = 0;
  double fd = 0;
  double d = 0;
  double lenDiff = 0;

  // A "segment" is a path from station s_i to station s_{i+1}

  size_t unmatchedSegments;        // number of unmatched segments
  double unmatchedSegmentsLength;  // total _an. length of unmatched segments

  std::vector<double> oldDists;
  LINE oldL = getLine(oldS, &oldDists);

  std::vector<double> newDists;
  LINE newL = getLine(newS, &newDists);

  // check dist between anchor points

  if ((util::geo::latLngLen(oldL) * 1.0) / (oldL.size() * 1.0) > 1000) {
    // most likely input with a degenerated shape - dont compare
    _noOrigShp++;
    return 0;
  }

  if ((util::geo::latLngLen(newL) * 1.0) / (newL.size() * 1.0) > 1000) {
    // most likely input with a degenerated shape - dont compare
    _noOrigShp++;
    return 0;
  }

  std::vector<std::pair<double, double>> newLenDists;
  std::vector<std::pair<double, double>> oldLenDists;

  auto oldSegs = segmentize(oldT, oldL, oldDists, newLenDists);
  auto newSegs = segmentize(newT, newL, newDists, oldLenDists);

  for (const auto& p : oldLenDists) {
    _distDiffs.push_back(fabs(p.first - p.second));
    _hopDists.push_back(p.first);
  }

  // new lines build from cleaned-up shapes
  LINE oldLCut;
  LINE newLCut;

  for (auto oldL : oldSegs)
    oldLCut.insert(oldLCut.end(), oldL.begin(), oldL.end());

  for (auto newL : newSegs) {
    newLCut.insert(newLCut.end(), newL.begin(), newL.end());
  }

  double f = util::geo::webMercDistFactor(oldLCut.front());

  // roughly half a meter
  auto oldLCutS =
      util::geo::simplify(oldLCut, f * (0.5 / util::geo::M_PER_DEG));
  auto newLCutS =
      util::geo::simplify(newLCut, f * (0.5 / util::geo::M_PER_DEG));

  auto old = _accFdCache.find(oldLCutS);
  if (old != _accFdCache.end()) {
    auto match = old->second.find(newLCutS);
    if (match != old->second.end()) {
      accFd = match->second;
    } else {
      accFd = util::geo::accFrechetDistCHav(oldLCutS, newLCutS, segLen);
      _accFdCache[oldLCutS][newLCutS] = accFd;
    }
  } else {
    accFd = util::geo::accFrechetDistCHav(oldLCutS, newLCutS, segLen);
    _accFdCache[oldLCutS][newLCutS] = accFd;
  }

  old = _fdCache.find(oldLCutS);
  if (old != _fdCache.end()) {
    auto match = old->second.find(newLCutS);
    if (match != old->second.end()) {
      fd = match->second;
    } else {
      fd = util::geo::frechetDistHav(oldLCutS, newLCutS, segLen);
      _fdCache[oldLCutS][newLCutS] = fd;
    }
  } else {
    fd = util::geo::frechetDistHav(oldLCutS, newLCutS, segLen);
    _fdCache[oldLCutS][newLCutS] = fd;
  }

  old = _dCache.find(oldLCutS);
  if (old != _dCache.end()) {
    auto match = old->second.find(newLCutS);
    if (match != old->second.end()) {
      d = match->second;
    } else {
      d = util::geo::dist(oldLCutS, newLCutS);
      _dCache[oldLCutS][newLCutS] = d;
    }
  } else {
    d = util::geo::dist(oldLCutS, newLCutS);
    _dCache[oldLCutS][newLCutS] = d;
  }

  old = _lenDiffCache.find(oldLCutS);
  if (old != _lenDiffCache.end()) {
    auto match = old->second.find(newLCutS);
    if (match != old->second.end()) {
      lenDiff = match->second;
    } else {
      lenDiff = util::geo::latLngLen(oldLCutS) - util::geo::latLngLen(newLCutS);
      _lenDiffCache[oldLCutS][newLCutS] = lenDiff;
    }
  } else {
    lenDiff = util::geo::latLngLen(oldLCutS) - util::geo::latLngLen(newLCutS);
    _lenDiffCache[oldLCutS][newLCutS] = lenDiff;
  }

  auto dA = getDa(oldSegs, newSegs, segLen);
  unmatchedSegments = dA.first;
  unmatchedSegmentsLength = dA.second;

  double totL = 0;
  for (auto l : oldSegs) totL += util::geo::latLngLen(l);

  // filter out shapes with a length of under 5 meters - they are most likely
  // artifacts
  if (totL < 5) {
    _noOrigShp++;
    return 0;
  }

  _accFdSum += accFd / totL;
  _unmatchedSegSum += unmatchedSegments;
  _unmatchedSegLengthSum += unmatchedSegmentsLength;

  double avgFd = accFd / totL;
  double AN = static_cast<double>(unmatchedSegments) /
              static_cast<double>(oldSegs.size());
  double AL = unmatchedSegmentsLength / totL;

  _results.insert(Result(oldT, avgFd));

  if (AN <= 0.0001) _an0++;
  if (AN <= 0.05) _an5++;
  if (AN <= 0.1) _an10++;
  if (AN <= 0.2) _an20++;
  if (AN <= 0.3) _an30++;
  if (AN <= 0.5) _an50++;
  if (AN <= 0.7) _an70++;
  if (AN <= 0.9) _an90++;

  LOG(VDEBUG) << "This result (" << oldT->getId()
              << "): A_N/N = " << unmatchedSegments << "/" << oldSegs.size()
              << " = " << AN << " A_L/L = " << unmatchedSegmentsLength << "/"
              << totL << " = " << AL << " d_f = " << avgFd;

  if (_reportOut) {
    if (_reportLevel == 0) {
      (*_reportOut) << std::fixed << std::setprecision(6);
      (*_reportOut) << oldT->getId() << "\t" << AN << "\t" << AL << "\t"
                    << avgFd << "\t" << oldT->getRoute()->getShortName()
                    << "\t";
    } else if (_reportLevel == 1) {
      (*_reportOut) << std::fixed << std::setprecision(6);
      (*_reportOut) << oldT->getId() << "\t" << AN << "\t" << AL << "\t"
                    << avgFd << "\t" << util::geo::getWKT(oldSegs) << "\t"
                    << util::geo::getWKT(newSegs) << "\t"
                    << oldT->getRoute()->getShortName() << "\t";

      for (const auto& st : oldT->getStopTimes()) {
        (*_reportOut) << st.getStop()->getName() << "\t"
                      << st.getStop()->getLat() << "\t"
                      << st.getStop()->getLng() << "\t";
      }
    } else if (_reportLevel == 2) {
      (*_reportOut) << std::fixed << std::setprecision(6);
      (*_reportOut) << oldT->getId() << "\t" << AN << "\t" << AL << "\t"
                    << avgFd << "\t" << fd << "\t"
                    << d
                    << "\t"
                    << lenDiff
                    << "\t" << util::geo::getWKT(oldSegs) << "\t"
                    << util::geo::getWKT(newSegs) << "\t"
                    << oldT->getRoute()->getShortName() << "\t"
                    << numOldTrips << "\t";

      for (const auto& st : oldT->getStopTimes()) {
        (*_reportOut) << st.getStop()->getName() << "\t"
                      << st.getStop()->getLat() << "\t"
                      << st.getStop()->getLng() << "\t";
      }
    }
    (*_reportOut) << "\n";
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
  std::vector<double> cuts;

  for (const auto& st : t->getStopTimes()) {
    cuts.push_back(st.getShapeDistanceTravelled());
  }

  size_t to =
      std::upper_bound(dists.begin(), dists.end(), cuts[0]) - dists.begin();

  POINT lastP;
  if (to >= dists.size()) {
    lastP = shape.back();
  } else if (to == 0) {
    lastP = shape.front();
  } else {
    double progr = (cuts[0] - dists[to - 1]) / (dists[to] - dists[to - 1]);
    lastP = shape[to - 1];
    lastP.setX(lastP.getX() +
               progr * (shape[to].getX() - shape[to - 1].getX()));
    lastP.setY(lastP.getY() +
               progr * (shape[to].getY() - shape[to - 1].getY()));
  }

  for (size_t i = 1; i < cuts.size(); i++) {
    size_t to =
        std::upper_bound(dists.begin(), dists.end(), cuts[i]) - dists.begin();

    POINT curP;
    if (to >= dists.size()) {
      curP = shape.back();
    } else if (to == 0) {
      curP = shape.front();
    } else {
      curP = shape[to - 1];
      double progr = (cuts[i] - dists[to - 1]) / (dists[to] - dists[to - 1]);
      curP.setX(curP.getX() +
                progr * (shape[to].getX() - shape[to - 1].getX()));
      curP.setY(curP.getY() +
                progr * (shape[to].getY() - shape[to - 1].getY()));
    }

    auto curL = pl.getSegment(lastP, curP).getLine();

    double dist =
        util::geo::haversine(t->getStopTimes()[i - 1].getStop()->getLat(),
                             t->getStopTimes()[i - 1].getStop()->getLng(),
                             t->getStopTimes()[i].getStop()->getLat(),
                             t->getStopTimes()[i].getStop()->getLng());
    double len = util::geo::latLngLen(curL);
    lenDist.push_back({dist, len});

    ret.push_back(curL);
    lastP = curP;
  }

  return ret;
}

// _____________________________________________________________________________
LINE Collector::getLine(const Shape* s, std::vector<double>* dists) {
  LINE ret;

  for (size_t i = 0; i < s->getPoints().size(); i++) {
    ret.push_back({s->getPoints()[i].lng, s->getPoints()[i].lat});
    (*dists).push_back(s->getPoints()[i].travelDist);
  }
  return ret;
}

// _____________________________________________________________________________
const std::set<Result>& Collector::getResults() const { return _results; }

// _____________________________________________________________________________
double Collector::getAvgDist() const { return _accFdSum / _results.size(); }

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
    (*os) << (static_cast<double>(_an20) /
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
    (*os) << "  an-20: "
          << (static_cast<double>(_an20) /
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

    // std::nth_element makes a partial sort of the first n elements
    std::nth_element(_distDiffs.begin(), i, _distDiffs.end());

    stats["median-dist-diff"] = *i;
  } else {
    stats["median-dist-diff"] = -1;
  }

  if (_hopDists.size()) {
    double s = 0;
    for (auto d : _hopDists) s += d;

    stats["avg-hop-dist"] = s / (_hopDists.size() * 1.0);
  } else {
    stats["avg-hop-dist"] = -1;
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
  stats["an-20"] =
      (static_cast<double>(_an20) / static_cast<double>(_results.size())) * 100;
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
                                           const std::vector<LINE>& b, double segLen) {
  assert(a.size() == b.size());
  std::pair<size_t, double> ret{0, 0};

  double MAX = 100;

  for (size_t i = 0; i < a.size(); i++) {
    double fdMeter = 0;

    double f = util::geo::webMercDistFactor(a[i].front());

    // roughly half a meter
    auto aSimpl = util::geo::simplify(a[i], f * (0.5 / util::geo::M_PER_DEG));
    auto bSimpl = util::geo::simplify(b[i], f * (0.5 / util::geo::M_PER_DEG));

    auto old = _dACache.find(aSimpl);
    if (old != _dACache.end()) {
      auto match = old->second.find(bSimpl);
      if (match != old->second.end()) {
        fdMeter = match->second;
      } else {
        fdMeter = util::geo::frechetDistHav(aSimpl, bSimpl, segLen);
        _dACache[aSimpl][bSimpl] = fdMeter;
      }
    } else {
      fdMeter = util::geo::frechetDistHav(aSimpl, bSimpl, segLen);
      _dACache[aSimpl][bSimpl] = fdMeter;
    }

    if (fdMeter >= MAX) {
      ret.first++;
      ret.second += util::geo::latLngLen(aSimpl);
    }
  }

  return ret;
}

// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <limits>
#include "pfaedle/router/Weights.h"

using pfaedle::router::DistDiffTransWeight;
using pfaedle::router::ExpoTransWeight;
using pfaedle::router::LineSimilarity;
using pfaedle::router::NormDistrTransWeight;
using util::geo::haversine;

// _____________________________________________________________________________
ExpoTransWeight::DistHeur::DistHeur(double maxV, const RoutingOpts& rOpts,
                                    const std::set<trgraph::Edge*>& tos)
    : _rOpts(rOpts), _maxV(maxV), _maxCentD(0), _lastE(0) {
  size_t c = 0;
  double x = 0, y = 0;

  for (const auto to : tos) {
    x += to->getFrom()->pl().getGeom()->getX();
    y += to->getFrom()->pl().getGeom()->getY();
    c++;
  }

  x /= c;
  y /= c;

  _center = POINT{x, y};

  for (const auto to : tos) {
    const double cur = haversine(*to->getFrom()->pl().getGeom(), _center);
    if (cur > _maxCentD) _maxCentD = cur;
  }

  _maxCentD /= _maxV;
}

// _____________________________________________________________________________
uint32_t ExpoTransWeight::DistHeur::operator()(
    const trgraph::Edge* a, const std::set<trgraph::Edge*>& b) const {
  UNUSED(b);

  // avoid repeated calculation for the same edge over and over again
  if (a == _lastE) return _lastC;

  _lastE = a;

  const double d = haversine(*a->getFrom()->pl().getGeom(), _center);
  const double heur = fmax(0, (d / _maxV - _maxCentD) * 10);

  // avoid overflow
  if (heur > std::numeric_limits<uint32_t>::max()) {
    _lastC = std::numeric_limits<uint32_t>::max();
    ;
    return _lastC;
  }

  _lastC = heur;
  return heur;
}

// _____________________________________________________________________________
uint32_t ExpoTransWeight::CostFunc::operator()(const trgraph::Edge* from,
                                               const trgraph::Node* n,
                                               const trgraph::Edge* to) const {
  if (!from) return 0;

  uint32_t c = from->pl().getCost();

  if (c == std::numeric_limits<uint32_t>::max()) return c;

  if (from == _lastFrom) {
    // the transit line simi calculation is independent of the "to" edge, so if
    // the last "from" edge was the same, skip it!
    c = _lastC;
  } else if (!_noLineSimiPen) {
    const auto& simi = transitLineSimi(from);

    if (!simi.nameSimilar) {
      if (_rOpts.lineUnmatchedPunishFact < 1) {
        c = std::ceil(static_cast<double>(c) * _rOpts.lineUnmatchedPunishFact);
      } else if (_rOpts.lineUnmatchedPunishFact > 1) {
        double a =
            std::round(static_cast<double>(c) * _rOpts.lineUnmatchedPunishFact);
        if (a > std::numeric_limits<uint32_t>::max())
          return std::numeric_limits<uint32_t>::max();
        c = a;
      }
    }

    if (!simi.fromSimilar) {
      if (_rOpts.lineNameFromUnmatchedPunishFact < 1) {
        c = std::ceil(static_cast<double>(c) *
                      _rOpts.lineNameFromUnmatchedPunishFact);
      } else if (_rOpts.lineNameFromUnmatchedPunishFact > 1) {
        double a = std::round(static_cast<double>(c) *
                              _rOpts.lineNameFromUnmatchedPunishFact);
        if (a > std::numeric_limits<uint32_t>::max())
          return std::numeric_limits<uint32_t>::max();
        c = a;
      }
    }

    if (!simi.toSimilar) {
      if (_rOpts.lineNameToUnmatchedPunishFact < 1) {
        c = std::ceil(static_cast<double>(c) *
                      _rOpts.lineNameToUnmatchedPunishFact);
      } else if (_rOpts.lineNameToUnmatchedPunishFact > 1) {
        double a = std::round(static_cast<double>(c) *
                              _rOpts.lineNameToUnmatchedPunishFact);
        if (a > std::numeric_limits<uint32_t>::max())
          return std::numeric_limits<uint32_t>::max();
        c = a;
      }
    }

    _lastC = c;
    _lastFrom = from;
  }

  uint32_t overflowCheck = c;

  if (n && !n->pl().isTurnCycle()) {
    if (_rOpts.fullTurnPunishFac != 0 && from->getFrom() == to->getTo() &&
        from->getTo() == to->getFrom()) {
      // trivial full turn
      c += _rOpts.fullTurnPunishFac;

      if (c <= overflowCheck) return std::numeric_limits<uint32_t>::max();
      overflowCheck = c;
    } else if (_rOpts.fullTurnPunishFac != 0 && n->getDeg() > 2) {
      // otherwise, only intersection angles will be punished

      double ang = util::geo::innerProd(
          *n->pl().getGeom(), from->pl().backHop(), to->pl().frontHop());

      if (ang < _rOpts.fullTurnAngle) {
        c += _rOpts.fullTurnPunishFac;
        if (c <= overflowCheck) return std::numeric_limits<uint32_t>::max();
        overflowCheck = c;
      }
    }

    // turn restriction cost
    if (_rOpts.turnRestrCost > 0 && from->pl().isRestricted() &&
        !_res.may(from, to, n)) {
      c += _rOpts.turnRestrCost;
      if (c <= overflowCheck) return std::numeric_limits<uint32_t>::max();
    }
  }

  return c;
}

// _____________________________________________________________________________
LineSimilarity ExpoTransWeight::CostFunc::transitLineSimi(
    const trgraph::Edge* e) const {
  if (_rAttrs.shortName.empty() && _rAttrs.lineFrom.empty() &&
      _rAttrs.lineTo.empty())
    return {true, true, true};

  LineSimilarity best = {false, false, false};
  for (const auto* l : e->pl().getLines()) {
    auto simi = _rAttrs.simi(l);
    if (simi.nameSimilar && simi.toSimilar && simi.fromSimilar) return simi;
    if (best < simi) best = simi;
  }

  return best;
}

// _____________________________________________________________________________
double ExpoTransWeight::weight(uint32_t c, double d, double t0, double d0,
                               const RoutingOpts& rOpts) {
  UNUSED(t0);
  UNUSED(d);
  UNUSED(d0);
  return rOpts.transitionPen * static_cast<double>(c) / 10.0;
}

// _____________________________________________________________________________
uint32_t ExpoTransWeight::invWeight(double c, const RoutingOpts& rOpts) {
  return std::round((c / rOpts.transitionPen) * 10.0);
}

// _____________________________________________________________________________
uint32_t ExpoTransWeight::maxCost(double tTime, const RoutingOpts& rOpts) {
  // abort after 3 times the scheduled time, but assume a min time of
  // 1 minute!
  return std::ceil(fmax(tTime, 60) * 3.0 * rOpts.lineUnmatchedPunishFact *
                   rOpts.lineNameToUnmatchedPunishFact *
                   rOpts.lineNameFromUnmatchedPunishFact * 10);
}

// _____________________________________________________________________________

// _____________________________________________________________________________
double NormDistrTransWeight::weight(uint32_t cs, double d, double t0, double d0,
                                    const RoutingOpts& rOpts) {
  UNUSED(d);
  UNUSED(d0);
  UNUSED(rOpts);

  double t = static_cast<double>(cs) / 10.0;

  // standard deviation of normal distribution
  double standarddev = 1;

  // no backwards time travel!
  if (t0 < 0) return std::numeric_limits<double>::infinity();

  // always assume it takes at least 10 seconds to travel
  t0 = fmax(10, t0);

  double cNorm = (t / t0 - 1) / standarddev;
  double normWeight = cNorm * cNorm;

  double expWeight = ExpoTransWeight::weight(cs, d, t0, d0, rOpts);

  return normWeight + expWeight;
}

// _____________________________________________________________________________
uint32_t NormDistrTransWeight::invWeight(double c, const RoutingOpts& rOpts) {
  UNUSED(rOpts);
  UNUSED(c);

  throw(std::runtime_error("Cannot apply inv weight to DistDiffTransWeight"));
}

// _____________________________________________________________________________

// _____________________________________________________________________________
double DistDiffTransWeight::weight(uint32_t c, double d, double t0, double d0,
                                   const RoutingOpts& rOpts) {
  UNUSED(t0);
  UNUSED(c);

  double w = fabs(d - d0);

  return rOpts.transitionPen * w;
}

// _____________________________________________________________________________
uint32_t DistDiffTransWeight::invWeight(double c, const RoutingOpts& rOpts) {
  UNUSED(rOpts);
  UNUSED(c);

  throw(std::runtime_error("Cannot apply inv weight to DistDiffTransWeight"));
}

// _____________________________________________________________________________
uint32_t DistDiffTransWeight::maxCost(double tTime, const RoutingOpts& rOpts) {
  UNUSED(tTime);
  UNUSED(rOpts);
  return std::numeric_limits<uint32_t>::max();
}

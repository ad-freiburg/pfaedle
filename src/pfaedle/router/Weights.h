// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_WEIGHTS_H_
#define PFAEDLE_ROUTER_WEIGHTS_H_

#include <set>
#include "pfaedle/osm/Restrictor.h"
#include "pfaedle/router/Misc.h"
#include "pfaedle/router/RoutingAttrs.h"
#include "pfaedle/trgraph/Graph.h"
#include "util/graph/EDijkstra.h"

namespace pfaedle {
namespace router {

typedef util::graph::EDijkstra::CostFunc<trgraph::NodePL, trgraph::EdgePL,
                                         uint32_t>
    RCostFunc;
typedef util::graph::EDijkstra::HeurFunc<trgraph::NodePL, trgraph::EdgePL,
                                         uint32_t>
    RHeurFunc;

class ExpoTransWeight {
 public:
  struct CostFunc : public RCostFunc {
    CostFunc(const RoutingAttrs& rAttrs, const RoutingOpts& rOpts,
             const osm::Restrictor& res, uint32_t max)
        : _rAttrs(rAttrs),
          _rOpts(rOpts),
          _res(res),
          _inf(max),
          _noLineSimiPen(false),
          _lastFrom(0) {
      if (_rAttrs.lineFrom.empty() && _rAttrs.lineTo.empty() &&
          _rAttrs.shortName.empty()) {
        _noLineSimiPen = true;
      }
      if (_rOpts.lineUnmatchedPunishFact == 1) {
        _noLineSimiPen = true;
      }
    }

    const RoutingAttrs& _rAttrs;
    const RoutingOpts& _rOpts;
    const osm::Restrictor& _res;
    uint32_t _inf;
    bool _noLineSimiPen;
    mutable const trgraph::Edge* _lastFrom;
    mutable uint32_t _lastC = 0;

    uint32_t operator()(const trgraph::Edge* from, const trgraph::Node* n,
                        const trgraph::Edge* to) const;
    uint32_t inf() const { return _inf; }

    LineSimilarity transitLineSimi(const trgraph::Edge* e) const;
  };

  struct DistHeur : RHeurFunc {
    DistHeur(double maxV, const RoutingOpts& rOpts,
             const std::set<trgraph::Edge*>& tos);

    const RoutingOpts& _rOpts;
    double _maxV;
    POINT _center;
    double _maxCentD;
    uint32_t operator()(const trgraph::Edge* a,
                        const std::set<trgraph::Edge*>& b) const;
    mutable const trgraph::Edge* _lastE;
    mutable uint32_t _lastC = 0;
  };

  static uint32_t maxCost(double tTime, const RoutingOpts& rOpts);
  static double weight(uint32_t c, double d, double t0, double d0,
                       const RoutingOpts& rOpts);
  static uint32_t invWeight(double cost, const RoutingOpts& rOpts);
  static const bool ALLOWS_FAST_ROUTE = true;
  static const bool NEED_DIST = false;
};

class ExpoTransWeightNoHeur : public ExpoTransWeight {
 public:
  struct DistHeur : RHeurFunc {
    DistHeur(double maxV, const RoutingOpts& rOpts,
             const std::set<trgraph::Edge*>& tos) {
      UNUSED(maxV);
      UNUSED(rOpts);
      UNUSED(tos);
    }

    uint32_t operator()(const trgraph::Edge* a,
                        const std::set<trgraph::Edge*>& b) const {
      UNUSED(a);
      UNUSED(b);
      return 0;
    }
  };
};

class NormDistrTransWeight : public ExpoTransWeight {
 public:
  static double weight(uint32_t c, double d, double t0, double d0,
                       const RoutingOpts& rOpts);
  static uint32_t invWeight(double cost, const RoutingOpts& rOpts);
  static const bool ALLOWS_FAST_ROUTE = false;
  static const bool NEED_DIST = false;
};

class NormDistrTransWeightNoHeur : public NormDistrTransWeight {
 public:
  struct DistHeur : RHeurFunc {
    DistHeur(double maxV, const RoutingOpts& rOpts,
             const std::set<trgraph::Edge*>& tos) {
      UNUSED(maxV);
      UNUSED(rOpts);
      UNUSED(tos);
    }

    uint32_t operator()(const trgraph::Edge* a,
                        const std::set<trgraph::Edge*>& b) const {
      UNUSED(a);
      UNUSED(b);
      return 0;
    }
  };
};

class DistDiffTransWeight : public ExpoTransWeight {
 public:
  static uint32_t maxCost(double tTime, const RoutingOpts& rOpts);
  static double weight(uint32_t c, double d, double t0, double d0,
                       const RoutingOpts& rOpts);
  static uint32_t invWeight(double cost, const RoutingOpts& rOpts);
  static const bool ALLOWS_FAST_ROUTE = false;
  static const bool NEED_DIST = true;
};

class DistDiffTransWeightNoHeur : public DistDiffTransWeight {
 public:
  struct DistHeur : RHeurFunc {
    DistHeur(double maxV, const RoutingOpts& rOpts,
             const std::set<trgraph::Edge*>& tos) {
      UNUSED(maxV);
      UNUSED(rOpts);
      UNUSED(tos);
    }

    uint32_t operator()(const trgraph::Edge* a,
                        const std::set<trgraph::Edge*>& b) const {
      UNUSED(a);
      UNUSED(b);
      return 0;
    }
  };
};

}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_WEIGHTS_H_

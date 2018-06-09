// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "pfaedle/router/EdgePL.h"
#include "pfaedle/router/Router.h"
#include "util/String.h"

using pfaedle::router::EdgePL;
using pfaedle::router::EdgeCost;
using pfaedle::router::EdgeList;
using pfaedle::trgraph::Node;

// _____________________________________________________________________________
EdgeList* EdgePL::getEdges() { return &_edges; }

// _____________________________________________________________________________
const EdgeList& EdgePL::getEdges() const { return _edges; }

// _____________________________________________________________________________
const FPoint& EdgePL::frontHop() const {
  if (!_edges.size()) return *_end->pl().getGeom();
  return _edges.back()->pl().frontHop();
}

// _____________________________________________________________________________
const FPoint& EdgePL::backHop() const {
  if (!_edges.size()) return *_start->pl().getGeom();
  return _edges.front()->pl().backHop();
}

// _____________________________________________________________________________
const Node* EdgePL::backNode() const { return _end; }

// _____________________________________________________________________________
const Node* EdgePL::frontNode() const { return _start; }

// _____________________________________________________________________________
const util::geo::FLine* EdgePL::getGeom() const {
  if (!_edges.size()) return 0;
  if (!_geom.size()) {
    const trgraph::Node* l = _start;
    for (auto i = _edges.rbegin(); i != _edges.rend(); i++) {
      const auto e = *i;
      if ((e->getFrom() == l) ^ e->pl().isRev()) {
        _geom.insert(_geom.end(), e->pl().getGeom()->begin(),
                     e->pl().getGeom()->end());
      } else {
        _geom.insert(_geom.end(), e->pl().getGeom()->rbegin(),
                     e->pl().getGeom()->rend());
      }
      l = e->getOtherNd(l);
    }
  }

  return &_geom;
}

// _____________________________________________________________________________
void EdgePL::setStartNode(const trgraph::Node* s) { _start = s; }

// _____________________________________________________________________________
void EdgePL::setEndNode(const trgraph::Node* e) { _end = e; }

// _____________________________________________________________________________
void EdgePL::setStartEdge(const trgraph::Edge* s) { _startE = s; }

// _____________________________________________________________________________
void EdgePL::setEndEdge(const trgraph::Edge* e) { _endE = e; }

// _____________________________________________________________________________
const EdgeCost& EdgePL::getCost() const { return _cost; }

// _____________________________________________________________________________
void EdgePL::setCost(const router::EdgeCost& c) { _cost = c; }


// _____________________________________________________________________________
void EdgePL::getAttrs(std::map<std::string, std::string>* obj) const {
  (*obj)["cost"] = std::to_string(_cost.getValue());
  (*obj)["from_edge"] = util::toString(_startE);
  (*obj)["to_edge"] = util::toString(_endE);
  (*obj)["cost_m_lvl1"] = std::to_string(_cost.meterDistLvl1);
  (*obj)["cost_m_lvl0"] = std::to_string(_cost.meterDist);
  (*obj)["cost_m_lvl1"] = std::to_string(_cost.meterDistLvl1);
  (*obj)["cost_m_lvl2"] = std::to_string(_cost.meterDistLvl2);
  (*obj)["cost_m_lvl3"] = std::to_string(_cost.meterDistLvl3);
  (*obj)["cost_m_lvl4"] = std::to_string(_cost.meterDistLvl4);
  (*obj)["cost_m_lvl5"] = std::to_string(_cost.meterDistLvl5);
  (*obj)["cost_m_lvl6"] = std::to_string(_cost.meterDistLvl6);
  (*obj)["cost_m_lvl7"] = std::to_string(_cost.meterDistLvl7);
  (*obj)["cost_fullturn"] = std::to_string(_cost.fullTurns);
  (*obj)["cost_st_passthru"] = std::to_string(_cost.passThruStations);
  (*obj)["cost_m_oneway"] = std::to_string(_cost.oneWayMeters);
  (*obj)["cost_m_lineunmatch"] = std::to_string(_cost.lineUnmatchedMeters);
  (*obj)["cost_reach_node_pen"] = std::to_string(_cost.reachPen);
  (*obj)["cost_oneway_event"] = std::to_string(_cost.oneWayEdges);
  (*obj)["dummy"] = _edges.size() ? "no" : "yes";
}

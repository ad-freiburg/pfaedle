// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <string>
#include <unordered_map>
#include "pfaedle/trgraph/NodePL.h"
#include "pfaedle/trgraph/StatGroup.h"
#include "pfaedle/trgraph/StatInfo.h"
#include "util/String.h"

using pfaedle::trgraph::StatInfo;
using pfaedle::trgraph::NodePL;
using pfaedle::trgraph::Component;

// we use the adress of this dummy station info as a special value
// of this node, meaning "is a station block". Re-using the _si field here
// saves some memory
StatInfo NodePL::_blockerSI = StatInfo();

std::unordered_map<const Component*, size_t> NodePL::_comps;

// _____________________________________________________________________________
NodePL::NodePL() : _geom(0, 0), _si(0), _component(0), _vis(0) {}

// _____________________________________________________________________________
NodePL::NodePL(const NodePL& pl)
    : _geom(pl._geom), _si(0), _component(pl._component), _vis(pl._vis) {
  if (pl._si) setSI(*(pl._si));
}

// _____________________________________________________________________________
NodePL::NodePL(const util::geo::FPoint& geom)
    : _geom(geom), _si(0), _component(0), _vis(0) {}

// _____________________________________________________________________________
NodePL::NodePL(const util::geo::FPoint& geom, const StatInfo& si)
    : _geom(geom), _si(0), _component(0), _vis(0) {
  setSI(si);
}

// _____________________________________________________________________________
NodePL::~NodePL() {
  if (getSI()) delete _si;
  if (_component) {
    _comps[_component]--;
    if (_comps[_component] == 0) {
      delete _component;
      _comps.erase(_comps.find(_component));
    }
  }
}

// _____________________________________________________________________________
void NodePL::setVisited() const { _vis = true; }

// _____________________________________________________________________________
void NodePL::setNoStat() { _si = 0; }

// _____________________________________________________________________________
const Component* NodePL::getComp() const { return _component; }

// _____________________________________________________________________________
void NodePL::setComp(const Component* c) {
  if (_component == c) return;
  _component = c;

  // NOT thread safe!
  if (!_comps.count(c))
    _comps[c] = 1;
  else
    _comps[c]++;
}

// _____________________________________________________________________________
const util::geo::FPoint* NodePL::getGeom() const { return &_geom; }

// _____________________________________________________________________________
void NodePL::setGeom(const util::geo::FPoint& geom) { _geom = geom; }

// _____________________________________________________________________________
void NodePL::getAttrs(std::map<std::string, std::string>* obj) const {
  (*obj)["component"] = std::to_string(reinterpret_cast<size_t>(_component));
  (*obj)["dijkstra_vis"] = _vis ? "yes" : "no";
  if (getSI()) {
    (*obj)["station_info_ptr"] = util::toString(_si);
    (*obj)["station_name"] = _si->getName();
    (*obj)["station_alt_names"] = util::implode(_si->getAltNames(), ",");
    (*obj)["from_osm"] = _si->isFromOsm() ? "yes" : "no";
    (*obj)["station_platform"] = _si->getTrack();
    (*obj)["station_group"] =
        std::to_string(reinterpret_cast<size_t>(_si->getGroup()));

    std::stringstream gtfsIds;
    if (_si->getGroup()) {
      for (auto* s : _si->getGroup()->getStops()) {
        gtfsIds << s->getId() << " (" << s->getName() << "),";
      }
    }

    (*obj)["station_group_stops"] = gtfsIds.str();
  }
}

// _____________________________________________________________________________
void NodePL::setSI(const StatInfo& si) { _si = new StatInfo(si); }

// _____________________________________________________________________________
const StatInfo* NodePL::getSI() const {
  if (isBlocker()) return 0;
  return _si;
}

// _____________________________________________________________________________
StatInfo* NodePL::getSI() {
  if (isBlocker()) return 0;
  return _si;
}

// _____________________________________________________________________________
void NodePL::setBlocker() { _si = &_blockerSI; }

// _____________________________________________________________________________
bool NodePL::isBlocker() const { return _si == &_blockerSI; }

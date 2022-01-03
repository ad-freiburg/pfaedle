// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <string>
#include <vector>
#include <limits>
#include <unordered_map>
#include "pfaedle/trgraph/NodePL.h"
#include "pfaedle/trgraph/StatInfo.h"
#include "util/String.h"

using pfaedle::trgraph::Component;
using pfaedle::trgraph::NodePL;
using pfaedle::trgraph::StatInfo;

std::vector<Component> NodePL::comps;
std::vector<StatInfo> NodePL::_statInfos;

// _____________________________________________________________________________
NodePL::NodePL()
    : _geom(0, 0),
      _si(0),
      _component(0)
#ifdef PFAEDLE_DBG
      ,
      _vis(0)
#endif
{
}

// _____________________________________________________________________________
NodePL::NodePL(const POINT& geom)
    : _geom(geom),
      _si(0),
      _component(0)
#ifdef PFAEDLE_DBG
      ,
      _vis(0)
#endif
{
}

// _____________________________________________________________________________
NodePL::NodePL(const POINT& geom, const StatInfo& si)
    : _geom(geom),
      _si(0),
      _component(0)
#ifdef PFAEDLE_DBG
      ,
      _vis(0)
#endif
{
  setSI(si);
}

// _____________________________________________________________________________
void NodePL::setVisited() const {
#ifdef PFAEDLE_DBG
  _vis = true;
#endif
}

// _____________________________________________________________________________
void NodePL::setNoStat() { _si = 0; }

// _____________________________________________________________________________
const Component& NodePL::getComp() const { return comps[_component - 1]; }

// _____________________________________________________________________________
uint32_t NodePL::getCompId() const { return _component; }

// _____________________________________________________________________________
void NodePL::setComp(uint32_t id) {
  _component = id;
}

// _____________________________________________________________________________
const POINT* NodePL::getGeom() const { return &_geom; }

// _____________________________________________________________________________
void NodePL::setGeom(const POINT& geom) { _geom = geom; }

// _____________________________________________________________________________
util::json::Dict NodePL::getAttrs() const {
  util::json::Dict obj;
  obj["component"] = std::to_string(_component);
#ifdef PFAEDLE_DBG
  obj["dijkstra_vis"] = _vis ? "yes" : "no";
#endif
  if (getSI()) {
    obj["station_info_ptr"] = util::toString(_si);
    obj["station_name"] = getSI()->getName();
    obj["station_alt_names"] =
        util::implode(getSI()->getAltNames(), ",");
    obj["station_platform"] = getSI()->getTrack();

#ifdef PFAEDLE_STATION_IDS
    // only print this in debug mode
    obj["station_id"] = getSI()->getId();
#endif
  }
  return obj;
}

// _____________________________________________________________________________
void NodePL::setSI(const StatInfo& si) {
  _statInfos.emplace_back(si);
  _si = _statInfos.size();
}

// _____________________________________________________________________________
const StatInfo* NodePL::getSI() const {
  if (isBlocker()) return 0;
  if (isTurnCycle()) return 0;
  if (_si == 0) return 0;
  return &_statInfos[_si - 1];
}

// _____________________________________________________________________________
StatInfo* NodePL::getSI() {
  if (isBlocker()) return 0;
  if (isTurnCycle()) return 0;
  if (_si == 0) return 0;
  return &_statInfos[_si - 1];
}

// _____________________________________________________________________________
void NodePL::setTurnCycle() { _si = std::numeric_limits<uint32_t>::max() - 1; }

// _____________________________________________________________________________
bool NodePL::isTurnCycle() const {
  return _si == (std::numeric_limits<uint32_t>::max() - 1);
}

// _____________________________________________________________________________
void NodePL::setBlocker() { _si = std::numeric_limits<uint32_t>::max(); }

// _____________________________________________________________________________
bool NodePL::isBlocker() const {
  return _si == std::numeric_limits<uint32_t>::max();
}

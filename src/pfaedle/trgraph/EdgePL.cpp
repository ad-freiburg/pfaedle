// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <map>
#include <string>
#include <vector>
#include "pfaedle/trgraph/EdgePL.h"

using pfaedle::trgraph::EdgePL;
using pfaedle::trgraph::TransitEdgeLine;

std::map<util::geo::FLine*, size_t> EdgePL::_flines;
std::map<const TransitEdgeLine*, size_t> EdgePL::_tlines;

// _____________________________________________________________________________
EdgePL::EdgePL()
    : _length(0), _oneWay(0), _hasRestr(false), _rev(false), _lvl(0) {
  _l = new util::geo::FLine();
  _flines[_l] = 1;
}

// _____________________________________________________________________________
EdgePL::EdgePL(const EdgePL& pl) : EdgePL(pl, false) {}

// _____________________________________________________________________________
EdgePL::EdgePL(const EdgePL& pl, bool geoflat)
    : _length(pl._length),
      _oneWay(pl._oneWay),
      _hasRestr(pl._hasRestr),
      _rev(pl._rev),
      _lvl(pl._lvl) {
  if (geoflat) {
    _l = pl._l;
  } else {
    _l = new util::geo::FLine(*pl._l);
  }
  _flines[_l]++;

  for (auto l : _lines) addLine(l);
}

// _____________________________________________________________________________
EdgePL::~EdgePL() {
  if (_l) {
    _flines[_l]--;
    if (_flines[_l] == 0) delete _l;
  }

  for (auto l : _lines) unRefTLine(l);
}

// _____________________________________________________________________________
void EdgePL::unRefTLine(const TransitEdgeLine* l) {
  if (l) {
    _tlines[l]--;
    if (_tlines[l] == 0) {
      delete l;
      _tlines.erase(_tlines.find(l));
    }
  }
}

// _____________________________________________________________________________
EdgePL EdgePL::revCopy() const {
  EdgePL ret(*this);
  ret.setRev();
  if (ret.oneWay() == 1)
    ret.setOneWay(2);
  else if (ret.oneWay() == 2)
    ret.setOneWay(1);
  return ret;
}

// _____________________________________________________________________________
void EdgePL::setLength(double d) { _length = d; }

// _____________________________________________________________________________
double EdgePL::getLength() const { return _length; }

// _____________________________________________________________________________
void EdgePL::addLine(const TransitEdgeLine* l) {
  if (_lines.insert(l).second) {
    if (_tlines.count(l))
      _tlines[l]++;
    else
      _tlines[l] = 1;
  }
}

// _____________________________________________________________________________
void EdgePL::addLines(const std::vector<TransitEdgeLine*>& l) {
  for (auto line : l) addLine(line);
}

// _____________________________________________________________________________
const std::set<const TransitEdgeLine*>& EdgePL::getLines() const {
  return _lines;
}

// _____________________________________________________________________________
void EdgePL::addPoint(const util::geo::FPoint& p) { _l->push_back(p); }

// _____________________________________________________________________________
const util::geo::FLine* EdgePL::getGeom() const { return _l; }

// _____________________________________________________________________________
util::geo::FLine* EdgePL::getGeom() { return _l; }

// _____________________________________________________________________________
util::json::Dict EdgePL::getAttrs() const {
  util::json::Dict obj;
  obj["m_length"] = std::to_string(_length);
  obj["oneway"] = std::to_string(static_cast<int>(_oneWay));
  obj["level"] = std::to_string(_lvl);
  obj["restriction"] = isRestricted() ? "yes" : "no";

  std::stringstream ss;
  bool first = false;

  for (auto* l : _lines) {
    if (first) ss << ",";
    ss << l->shortName;
    if (l->fromStr.size() || l->toStr.size()) {
      ss << "(" << l->fromStr;
      ss << "->" << l->toStr << ")";
    }
    first = true;
  }

  obj["lines"] = ss.str();
  return obj;
}

// _____________________________________________________________________________
void EdgePL::setRestricted() { _hasRestr = true; }

// _____________________________________________________________________________
bool EdgePL::isRestricted() const { return _hasRestr; }

// _____________________________________________________________________________
uint8_t EdgePL::oneWay() const { return _oneWay; }

// _____________________________________________________________________________
void EdgePL::setOneWay(uint8_t dir) { _oneWay = dir; }

// _____________________________________________________________________________
void EdgePL::setOneWay() { _oneWay = 1; }

// _____________________________________________________________________________
void EdgePL::setLvl(uint8_t lvl) { _lvl = lvl; }

// _____________________________________________________________________________
uint8_t EdgePL::lvl() const { return _lvl; }

// _____________________________________________________________________________
void EdgePL::setRev() { _rev = true; }

// _____________________________________________________________________________
bool EdgePL::isRev() const { return _rev; }

// _____________________________________________________________________________
const util::geo::FPoint& EdgePL::backHop() const {
  if (isRev()) {
    return *(++(getGeom()->cbegin()));
  }
  return *(++(getGeom()->crbegin()));
}

// _____________________________________________________________________________
const util::geo::FPoint& EdgePL::frontHop() const {
  if (isRev()) {
    return *(++(getGeom()->crbegin()));
  }
  return *(++(getGeom()->cbegin()));
}

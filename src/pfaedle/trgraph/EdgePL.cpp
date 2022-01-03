// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <map>
#include <string>
#include <vector>
#include "pfaedle/trgraph/EdgePL.h"
#include "util/geo/Geo.h"

using pfaedle::trgraph::EdgePL;
using pfaedle::trgraph::TransitEdgeLine;

std::map<LINE*, size_t> EdgePL::_flines;
std::map<const TransitEdgeLine*, size_t> EdgePL::_tlines;

// _____________________________________________________________________________
EdgePL::EdgePL()
    : _oneWay(0), _hasRestr(false), _rev(false), _lvl(0), _cost(0), _l(0) {
}

// _____________________________________________________________________________
EdgePL::EdgePL(const EdgePL& pl) : EdgePL(pl, false) {}

// _____________________________________________________________________________
EdgePL::EdgePL(const EdgePL& pl, bool geoflat)
    : _oneWay(pl._oneWay),
      _hasRestr(pl._hasRestr),
      _rev(pl._rev),
      _lvl(pl._lvl),
      _cost(pl._cost),
      _l(0) {
  if (pl._l) {
    if (geoflat) {
      _l = pl._l;
    } else {
      _l = new LINE(*pl._l);
    }
    _flines[_l]++;
  }

  for (auto l : pl._lines) addLine(l);
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
double EdgePL::getLength() const {
  double len = 0;

  for (size_t i = 1; i < _l->size(); i++) {
    len += haversine((*_l)[i-1], (*_l)[i]);
  }

  return len;
}

// _____________________________________________________________________________
void EdgePL::addLine(const TransitEdgeLine* l) {
  auto lb = std::lower_bound(_lines.begin(), _lines.end(), l);
  if (lb == _lines.end() || *lb != l) {
    _lines.reserve(_lines.size() + 1);
    lb = std::lower_bound(_lines.begin(), _lines.end(), l);
    _lines.insert(lb, l);
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
const std::vector<const TransitEdgeLine*>& EdgePL::getLines() const {
  return _lines;
}

// _____________________________________________________________________________
void EdgePL::addPoint(const POINT& p) {
  if (!_l) {
    _l = new LINE();
    _flines[_l] = 1;
  }
  _l->push_back(p);
}

// _____________________________________________________________________________
const LINE* EdgePL::getGeom() const { return _l; }

// _____________________________________________________________________________
LINE* EdgePL::getGeom() { return _l; }

// _____________________________________________________________________________
util::json::Dict EdgePL::getAttrs() const {
  util::json::Dict obj;
  obj["m_length"] = std::to_string(getLength());
  obj["oneway"] = std::to_string(static_cast<int>(_oneWay));
  obj["cost"] = std::to_string(static_cast<double>(_cost) / 10.0);
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
uint32_t EdgePL::getCost() const { return _cost; }

// _____________________________________________________________________________
void EdgePL::setCost(uint32_t c) { _cost = c; }

// _____________________________________________________________________________
void EdgePL::setRev() { _rev = true; }

// _____________________________________________________________________________
bool EdgePL::isRev() const { return _rev; }

// _____________________________________________________________________________
const POINT& EdgePL::backHop() const {
  if (isRev()) {
    return *(++(getGeom()->cbegin()));
  }
  return *(++(getGeom()->crbegin()));
}

// _____________________________________________________________________________
const POINT& EdgePL::frontHop() const {
  if (isRev()) {
    return *(++(getGeom()->crbegin()));
  }
  return *(++(getGeom()->cbegin()));
}

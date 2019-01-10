// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <iostream>
#include <sstream>
#include <string>
#include "pfaedle/osm/OsmFilter.h"

using pfaedle::osm::OsmFilter;

// _____________________________________________________________________________
OsmFilter::OsmFilter(const MultAttrMap& keep, const MultAttrMap& drop)
    : _keep(keep), _drop(drop) {}

// _____________________________________________________________________________
OsmFilter::OsmFilter(const OsmReadOpts& o)
    : _keep(o.keepFilter),
      _drop(o.dropFilter),
      _nohup(o.noHupFilter),
      _oneway(o.oneWayFilter),
      _onewayrev(o.oneWayFilterRev),
      _twoway(o.twoWayFilter),
      _station(o.stationFilter),
      _blocker(o.stationBlockerFilter),
      _posRestr(o.restrPosRestr),
      _negRestr(o.restrNegRestr),
      _noRestr(o.noRestrFilter),
      _levels(o.levelFilters) {}

// _____________________________________________________________________________
uint64_t OsmFilter::keep(const AttrMap& attrs, Type t) const {
  return contained(attrs, _keep, t);
}

// _____________________________________________________________________________
uint64_t OsmFilter::drop(const AttrMap& attrs, Type t) const {
  return contained(attrs, _drop, t);
}

// _____________________________________________________________________________
uint64_t OsmFilter::nohup(const char* key, const char* v) const {
  const auto& dkv = _nohup.find(key);
  if (dkv != _nohup.end()) {
    for (const auto& val : dkv->second) {
      if (valMatches(v, val.first)) return true;
    }
  }

  return false;
}

// _____________________________________________________________________________
uint64_t OsmFilter::oneway(const AttrMap& attrs) const {
  if (contained(attrs, _twoway, WAY)) return false;
  return contained(attrs, _oneway, WAY);
}

// _____________________________________________________________________________
uint64_t OsmFilter::onewayrev(const AttrMap& attrs) const {
  if (contained(attrs, _twoway, WAY)) return false;
  return contained(attrs, _onewayrev, WAY);
}

// _____________________________________________________________________________
uint64_t OsmFilter::station(const AttrMap& attrs) const {
  return contained(attrs, _station, NODE);
}

// _____________________________________________________________________________
uint64_t OsmFilter::blocker(const AttrMap& attrs) const {
  return contained(attrs, _blocker, NODE);
}

// _____________________________________________________________________________
uint64_t OsmFilter::contained(const AttrMap& attrs, const MultAttrMap& map,
                              Type t) {
  for (const auto& kv : attrs) {
    const auto& dkv = map.find(kv.first);

    if (dkv != map.end()) {
      for (const auto& val : dkv->second) {
        bool multValMatch = val.second & osm::MULT_VAL_MATCH;
        if (val.second & t) continue;
        if (valMatches(kv.second, val.first, multValMatch)) return val.second;
      }
    }
  }

  return 0;
}

// _____________________________________________________________________________
uint64_t OsmFilter::contained(const AttrMap& attrs, const Attr& attr) {
  for (const auto& kv : attrs) {
    if (kv.first == attr.first) return valMatches(kv.second, attr.second);
  }

  return 0;
}

// _____________________________________________________________________________
uint8_t OsmFilter::level(const AttrMap& attrs) const {
  // the best matching level is always returned
  for (int16_t i = 0; i < 8; i++) {
    for (const auto& kv : attrs) {
      const auto& lkv = (_levels + i)->find(kv.first);
      if (lkv != (_levels + i)->end()) {
        for (const auto& val : lkv->second) {
          if (valMatches(kv.second, val.first)) return i;
        }
      }
    }
  }

  return 0;
}

// _____________________________________________________________________________
bool OsmFilter::valMatches(const std::string& a, const std::string& b) {
  return valMatches(a, b, false);
}

// _____________________________________________________________________________
bool OsmFilter::valMatches(const std::string& a, const std::string& b, bool m) {
  if (b == "*") return true;

  if (m) {
    // search for occurances in semicolon separated list
    if (a.find(std::string(";") + b) != std::string::npos) return true;
    if (a.find(b + ";") != std::string::npos) return true;
    if (a.find(std::string("; ") + b) != std::string::npos) return true;
    if (a.find(b + " ;") != std::string::npos) return true;
  }

  return a == b;
}

// _____________________________________________________________________________
std::vector<std::string> OsmFilter::getAttrKeys() const {
  std::vector<std::string> ret;
  for (const auto& kv : _keep) {
    ret.push_back(kv.first);
  }
  for (const auto& kv : _drop) {
    ret.push_back(kv.first);
  }
  for (const auto& kv : _nohup) {
    ret.push_back(kv.first);
  }
  for (const auto& kv : _oneway) {
    ret.push_back(kv.first);
  }
  for (const auto& kv : _twoway) {
    ret.push_back(kv.first);
  }
  for (const auto& kv : _station) {
    ret.push_back(kv.first);
  }
  for (const auto& kv : _blocker) {
    ret.push_back(kv.first);
  }
  for (const auto& kv : _posRestr) {
    ret.push_back(kv.first);
  }
  for (const auto& kv : _negRestr) {
    ret.push_back(kv.first);
  }
  for (const auto& kv : _noRestr) {
    ret.push_back(kv.first);
  }
  for (uint8_t i = 0; i < 8; i++) {
    for (const auto& kv : *(_levels + i)) {
      ret.push_back(kv.first);
    }
  }

  return ret;
}

// _____________________________________________________________________________
OsmFilter OsmFilter::merge(const OsmFilter& other) const {
  MultAttrMap keep;
  MultAttrMap drop;

  for (const auto& kv : _keep) {
    keep[kv.first].insert(kv.second.begin(), kv.second.end());
  }

  for (const auto& kv : other._keep) {
    keep[kv.first].insert(kv.second.begin(), kv.second.end());
  }

  return OsmFilter(keep, drop);
}

// _____________________________________________________________________________
std::string OsmFilter::toString() const {
  std::stringstream ss;
  ss << "[KEEP]\n\n";

  for (const auto& kv : _keep) {
    ss << " " << kv.first << "=";
    bool first = false;
    for (const auto& v : kv.second) {
      if (first) ss << ",";
      first = true;
      ss << v.first;
    }
    ss << "\n";
  }

  ss << "\n[DROP]\n\n";

  for (const auto& kv : _drop) {
    ss << " " << kv.first << "=";
    bool first = false;
    for (const auto& v : kv.second) {
      if (first) ss << ",";
      first = true;
      ss << v.first;
    }
    ss << "\n";
  }

  return ss.str();
}

// _____________________________________________________________________________
uint64_t OsmFilter::negRestr(const AttrMap& attrs) const {
  if (contained(attrs, _noRestr, ALL)) return false;
  return (contained(attrs, _negRestr, ALL));
}

// _____________________________________________________________________________
uint64_t OsmFilter::posRestr(const AttrMap& attrs) const {
  if (contained(attrs, _noRestr, ALL)) return false;
  return (contained(attrs, _posRestr, ALL));
}

// _____________________________________________________________________________
const pfaedle::osm::MultAttrMap& OsmFilter::getKeepRules() const {
  return _keep;
}

// _____________________________________________________________________________
const pfaedle::osm::MultAttrMap& OsmFilter::getDropRules() const {
  return _drop;
}

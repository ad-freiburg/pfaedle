// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "pfaedle/router/Comp.h"
#include "pfaedle/trgraph/StatGroup.h"
#include "pfaedle/trgraph/StatInfo.h"

using pfaedle::trgraph::StatInfo;
using pfaedle::trgraph::StatGroup;

std::unordered_map<const StatGroup*, size_t> StatInfo::_groups;

// _____________________________________________________________________________
StatInfo::StatInfo() : _name(""), _track(""), _fromOsm(false), _group(0) {}

// _____________________________________________________________________________
StatInfo::StatInfo(const StatInfo& si)
    : _name(si._name),
      _altNames(si._altNames),
      _track(si._track),
      _fromOsm(si._fromOsm),
      _group(0) {
  setGroup(si._group);
#ifdef PFAEDLE_STATION_IDS
  _id = si._id;
#endif
}

// _____________________________________________________________________________
StatInfo::StatInfo(const std::string& name, const std::string& track,
                   bool fromOsm)
    : _name(name), _track(track), _fromOsm(fromOsm), _group(0) {}

// _____________________________________________________________________________
StatInfo::~StatInfo() { unRefGroup(_group); }

// _____________________________________________________________________________
void StatInfo::unRefGroup(StatGroup* g) {
  if (g) {
    _groups[g]--;
    if (_groups[g] == 0) {
      // std::cout << "Deleting " << g << std::endl;
      delete g;
      _groups.erase(_groups.find(g));
    }
  }
}

// _____________________________________________________________________________
void StatInfo::setGroup(StatGroup* g) {
  if (_group == g) return;
  unRefGroup(_group);

  _group = g;

  // NOT thread safe!
  if (!_groups.count(g))
    _groups[g] = 1;
  else
    _groups[g]++;
}

// _____________________________________________________________________________
StatGroup* StatInfo::getGroup() const { return _group; }

// _____________________________________________________________________________
const std::string& StatInfo::getName() const { return _name; }

// _____________________________________________________________________________
const std::string& StatInfo::getTrack() const { return _track; }

// _____________________________________________________________________________
bool StatInfo::isFromOsm() const { return _fromOsm; }

// _____________________________________________________________________________
void StatInfo::setIsFromOsm(bool is) { _fromOsm = is; }

// _____________________________________________________________________________
double StatInfo::simi(const StatInfo* other) const {
  if (!other) return 0;
  if (router::statSimi(_name, other->getName()) > 0.5) return 1;

  for (const auto& a : _altNames) {
    if (router::statSimi(a, other->getName()) > 0.5) return 1;
    for (const auto& b : other->getAltNames()) {
      if (router::statSimi(a, b) > 0.5) return 1;
    }
  }

  for (const auto& b : other->getAltNames()) {
    if (router::statSimi(_name, b) > 0.5) return 1;
  }

  return 0;
}

// _____________________________________________________________________________
const std::vector<std::string>& StatInfo::getAltNames() const {
  return _altNames;
}

// _____________________________________________________________________________
void StatInfo::addAltName(const std::string& name) {
  _altNames.push_back(name);
}

// _____________________________________________________________________________
void StatInfo::setTrack(const std::string& tr) { _track = tr; }

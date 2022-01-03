// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "pfaedle/router/Comp.h"
#include "pfaedle/trgraph/StatInfo.h"

using pfaedle::trgraph::StatInfo;

// _____________________________________________________________________________
StatInfo::StatInfo() : _name(""), _track("") {}

// _____________________________________________________________________________
StatInfo::StatInfo(const StatInfo& si)
    : _name(si._name), _altNames(si._altNames), _track(si._track) {
#ifdef PFAEDLE_STATION_IDS
  _id = si._id;
#endif
}

// _____________________________________________________________________________
StatInfo::StatInfo(const std::string& name, const std::string& track)
    : _name(name), _track(track) {}

// _____________________________________________________________________________
const std::string& StatInfo::getName() const { return _name; }

// _____________________________________________________________________________
const std::string& StatInfo::getTrack() const { return _track; }

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

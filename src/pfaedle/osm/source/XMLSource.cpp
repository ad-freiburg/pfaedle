// Copyright 2024, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "pfaedle/osm/source/XMLSource.h"
#include "pfxml/pfxml.h"
#include "util/Misc.h"

using pfaedle::osm::source::XMLSource;
using pfaedle::osm::source::OsmSourceNode;
using pfaedle::osm::source::OsmSourceWay;
using pfaedle::osm::source::OsmSourceRelation;
using pfaedle::osm::source::OsmSourceRelationMember;
using pfaedle::osm::source::OsmSourceAttr;

// _____________________________________________________________________________
XMLSource::XMLSource(const std::string& path) : _path(path), _xml(path) {
  // init states
  _start = _xml.state();
  _nodeBeg = _xml.state();
  _wayBeg = _xml.state();
  _relBeg = _xml.state();

  // seek node begin
  seekNodes();

  // set individual states to node begin
  _nodeBeg = _xml.state();
  _wayBeg = _xml.state();
  _relBeg = _xml.state();
}

// _____________________________________________________________________________
const OsmSourceNode* XMLSource::nextNode() {
  do {
    const pfxml::tag& cur = _xml.get();

    if (_xml.level() != 2) continue;
    if (!_inNodeBlock && strcmp(cur.name, "node") == 0) _inNodeBlock = true;

    if (_inNodeBlock) {
      // block ended
      if (strcmp(cur.name, "node")) {
        _wayBeg = _xml.state();
        return 0;
      }

      _curNode.lon = util::atof(cur.attr("lon"), 7);
      _curNode.lat = util::atof(cur.attr("lat"), 7);

      _curNode.id = util::atoul(cur.attr("id"));

      return &_curNode;
    }
  } while (_xml.next());

  return 0;
}

// _____________________________________________________________________________
void XMLSource::seekNodes() {
  _xml.set_state(_nodeBeg);
  while (strcmp(_xml.get().name, "node") && _xml.next()) {}
}

// _____________________________________________________________________________
void XMLSource::seekWays() {
  _xml.set_state(_wayBeg);
  while (strcmp(_xml.get().name, "way") && _xml.next()) {}
}

// _____________________________________________________________________________
void XMLSource::seekRels() {
  _xml.set_state(_relBeg);
  while (strcmp(_xml.get().name, "relation") && _xml.next()) {}
}

// _____________________________________________________________________________
bool XMLSource::cont() {
  return _xml.next();
}

// _____________________________________________________________________________
const OsmSourceWay* XMLSource::nextWay() {
  do {
    const pfxml::tag& cur = _xml.get();

    if (_xml.level() != 2) continue;
    if (!_inWayBlock && strcmp(cur.name, "way") == 0) _inWayBlock = true;

    if (_inWayBlock) {
      // block ended
      if (strcmp(cur.name, "way")) {
        _relBeg = _xml.state();
        return 0;
      }
      _curWay.id = util::atoul(cur.attr("id"));

      return &_curWay;
    }

  } while (_xml.next());

  return 0;
}

// _____________________________________________________________________________
const OsmSourceRelationMember* XMLSource::nextMember() {
  do {
    const pfxml::tag& cur = _xml.get();

    if (_xml.level() != 3) return 0;

    if (strcmp(cur.name, "member") == 0) {
      _curMember.id = util::atoul(cur.attr("ref"));
      _curMember.type = 0;
      if (strcmp(cur.attr("type"), "way") == 0) _curMember.type = 1;
      else if (strcmp(cur.attr("type"), "relation") == 0) _curMember.type = 2;
      _curMember.role = cur.attr("role");
      if (!_curMember.role) _curMember.role = "";
      return &_curMember;
    } else {
      return 0;
    }
  } while (_xml.next());

  return 0;
}

// _____________________________________________________________________________
uint64_t XMLSource::nextMemberNode() {
  do {
    const pfxml::tag& cur = _xml.get();

    if (_xml.level() != 3) return 0;

    if (strcmp(cur.name, "nd") == 0) {
      return util::atoul(cur.attr("ref"));
    } else {
      return 0;
    }
  } while (_xml.next());

  return 0;
}

// _____________________________________________________________________________
const OsmSourceRelation* XMLSource::nextRel() {
  do {
    const pfxml::tag& cur = _xml.get();

    if (_xml.level() != 2) continue;
    if (!_inRelBlock && strcmp(cur.name, "relation") == 0) _inRelBlock = true;

    if (_inRelBlock) {
      // block ended
      if (strcmp(cur.name, "relation")) return 0;
      _curRel.id = util::atoul(cur.attr("id"));

      return &_curRel;
    }

  } while (_xml.next());

  return 0;
}

// _____________________________________________________________________________
const OsmSourceAttr XMLSource::nextAttr() {
  do {
    const pfxml::tag& cur = _xml.get();

    if (_xml.level() != 3) return {0, 0};

    if (strcmp(cur.name, "tag") == 0) {
      return {cur.attr("k"), cur.attr("v")};
    } else {
      return {0, 0};
    }
  } while (_xml.next());

  return {0, 0};
}


// _____________________________________________________________________________
util::geo::Box<double> XMLSource::getBounds() {
  _xml.set_state(_start);

  while (_xml.next() && strcmp(_xml.get().name, "bounds")) {}

  const pfxml::tag& cur = _xml.get();

  if (strcmp(cur.name, "bounds") != 0) {
    throw pfxml::parse_exc(
        std::string("Could not find required <bounds> tag"), _path, 0, 0, 0);
  }

  if (!cur.attr("minlat")) {
    throw pfxml::parse_exc(
        std::string(
            "Could not find required attribute \"minlat\" for <bounds> tag"),
        _path, 0, 0, 0);
  }
  if (!cur.attr("minlon")) {
    throw pfxml::parse_exc(
        std::string(
            "Could not find required attribute \"minlon\" for <bounds> tag"),
        _path, 0, 0, 0);
  }
  if (!cur.attr("maxlat")) {
    throw pfxml::parse_exc(
        std::string(
            "Could not find required attribute \"maxlat\" for <bounds> tag"),
        _path, 0, 0, 0);
  }
  if (!cur.attr("maxlon")) {
    throw pfxml::parse_exc(
        std::string(
            "Could not find required attribute \"maxlon\" for <bounds> tag"),
        _path, 0, 0, 0);
  }

  double minlat = atof(cur.attr("minlat"));
  double minlon = atof(cur.attr("minlon"));
  double maxlat = atof(cur.attr("maxlat"));
  double maxlon = atof(cur.attr("maxlon"));

  return util::geo::Box<double>({minlon, minlat}, {maxlon, maxlat});
}

// _____________________________________________________________________________
std::string XMLSource::decode(const char* str) const {
  return pfxml::file::decode(str);
}

// _____________________________________________________________________________
std::string XMLSource::decode(const std::string& str) const {
  return pfxml::file::decode(str);
}

// Copyright 2024, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_OSM_SOURCE_OSMSOURCE_H_
#define PFAEDLE_OSM_SOURCE_OSMSOURCE_H_

#include <stdint.h>
#include <vector>
#include "util/geo/Geo.h"

namespace pfaedle {
namespace osm {
namespace source {

struct OsmSourceNode {
  uint64_t id;
  double lat;
  double lon;
};

struct OsmSourceAttr {
  const char* key;
  const char* value;
};

struct OsmSourceWay {
  uint64_t id;
};

struct OsmSourceRelation {
  uint64_t id;
};

struct OsmSourceRelationMember {
  uint64_t id;
  uint8_t type;
  const char* role;
};

class OsmSource {
 public:
  virtual const OsmSourceNode* nextNode() = 0;
  virtual const OsmSourceAttr nextAttr() = 0;
  virtual const OsmSourceWay* nextWay() = 0;
  virtual uint64_t nextMemberNode() = 0;
  virtual const OsmSourceRelationMember* nextMember() = 0;
  virtual const OsmSourceRelation* nextRel() = 0;
  virtual bool cont() = 0;

  virtual ~OsmSource() {};

  virtual util::geo::Box<double> getBounds() = 0;

  virtual void seekNodes() = 0;
  virtual void seekWays() = 0;
  virtual void seekRels() = 0;

  virtual std::string decode(const char* str) const = 0;
  virtual std::string decode(const std::string& str) const = 0;
};

}  // namespace source
}  // namespace osm
}  // namespace pfaedle

#endif

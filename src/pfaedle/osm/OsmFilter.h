// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_OSM_OSMFILTER_H_
#define PFAEDLE_OSM_OSMFILTER_H_

#include <stdint.h>

#include <string>
#include <vector>

#include "pfaedle/osm/Osm.h"
#include "pfaedle/osm/OsmReadOpts.h"

namespace pfaedle {
namespace osm {

class OsmFilter {
 public:
  enum Type : uint64_t { NODE = 16, WAY = 8, REL = 4, ALL = 0 };
  OsmFilter() {}
  OsmFilter(const MultAttrMap& keep, const MultAttrMap& drop);
  explicit OsmFilter(const OsmReadOpts& o);
  uint64_t keep(const AttrMap& attrs, Type t) const;
  uint64_t drop(const AttrMap& attrs, Type t) const;
  uint64_t nohup(const char* key, const char* val) const;
  uint8_t level(const AttrMap& attrs) const;
  uint64_t oneway(const AttrMap& attrs) const;
  uint64_t onewayrev(const AttrMap& attrs) const;
  uint64_t station(const AttrMap& attrs) const;
  uint64_t blocker(const AttrMap& attrs) const;
  uint64_t turnCycle(const AttrMap& attrs) const;
  uint64_t negRestr(const AttrMap& attrs) const;
  uint64_t posRestr(const AttrMap& attrs) const;
  std::vector<std::string> getAttrKeys() const;

  OsmFilter merge(const OsmFilter& other) const;

  const MultAttrMap& getKeepRules() const;
  const MultAttrMap& getDropRules() const;

  std::string toString() const;

  static bool valMatches(const std::string& a, const std::string& b, bool m);
  static bool valMatches(const std::string& a, const std::string& b);
  static uint64_t contained(const AttrMap& attrs, const MultAttrMap& map,
                            Type t);
  static uint64_t contained(const AttrMap& attrs, const Attr& map);

 private:
  MultAttrMap _keep, _drop, _nohup, _oneway, _onewayrev, _twoway, _station,
      _blocker, _posRestr, _negRestr, _noRestr, _turnCycle;
  const MultAttrMap* _levels;
};
}  // namespace osm
}  // namespace pfaedle
#endif  // PFAEDLE_OSM_OSMFILTER_H_

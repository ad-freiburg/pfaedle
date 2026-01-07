// Copyright 2025, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_OSM_WRITER_OSMWRITER_H_
#define PFAEDLE_OSM_WRITER_OSMWRITER_H_

#include <stdint.h>

#include <vector>

#include "pfaedle/osm/Osm.h"
#include "pfaedle/osm/source/OsmSource.h"
#include "util/geo/Geo.h"

namespace pfaedle {
namespace osm {
namespace output {

class OsmWriter {
 public:
  virtual void writeNode(const OsmNode& node) = 0;
  virtual void writeWay(const OsmWay& node) = 0;
  virtual void writeRel(const OsmRel& nd, const OsmIdList& nodes,
                        const OsmIdList& ways,
                        std::vector<const char*> nodeRoles,
                        std::vector<const char*> wayRoles) = 0;

  virtual ~OsmWriter(){};
};

}  // namespace output
}  // namespace osm
}  // namespace pfaedle

#endif

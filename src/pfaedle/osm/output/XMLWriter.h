// Copyright 2024, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_OSM_WRITER_XMLWRITER_H_
#define PFAEDLE_OSM_WRITER_XMLWRITER_H_

#include "pfaedle/osm/output/OsmWriter.h"
#include "pfaedle/osm/source/OsmSource.h"
#include "util/xml/XmlWriter.h"

namespace pfaedle {
namespace osm {
namespace output {

class XMLWriter : public OsmWriter {
 public:
  XMLWriter(const std::string& path, const util::geo::Box<double>& latLngBox,
            const source::OsmSource* source);
  ~XMLWriter();
  virtual void writeNode(const OsmNode& node);
  virtual void writeWay(const OsmWay& node);
  virtual void writeRel(const OsmRel& nd, const OsmIdList& nodes,
                        const OsmIdList& ways,
                        std::vector<const char*> nodeRoles,
                        std::vector<const char*> wayRoles);

 private:
  util::xml::XmlWriter _wr;
  const source::OsmSource* _source;
};

}  // namespace output
}  // namespace osm
}  // namespace pfaedle

#endif

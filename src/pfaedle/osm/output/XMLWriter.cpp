// Copyright 2024, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "pfaedle/osm/output/XMLWriter.h"
#include "util/Misc.h"
#include "pfaedle/_config.h"

using pfaedle::osm::output::XMLWriter;
using pfaedle::osm::source::OsmSource;

// _____________________________________________________________________________
XMLWriter::XMLWriter(const std::string& path,
                     const util::geo::Box<double>& latLngBox,
                     const OsmSource* source)
    : _wr(path, false, 0), _source(source) {
  _wr.put("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  _wr.openTag("osm", {{"version", "0.6"},
                     {"generator", std::string("pfaedle/") + VERSION_FULL}});
  _wr.openTag(
      "bounds",
      {{"minlat", util::formatFloat(latLngBox.getLowerLeft().getY(), 7)},
       {"minlon", util::formatFloat(latLngBox.getLowerLeft().getX(), 7)},
       {"maxlat",
       util::formatFloat(latLngBox.getUpperRight().getY(), 7)},
       {"maxlon",
       util::formatFloat(latLngBox.getUpperRight().getX(), 7)}});
  _wr.closeTag();
}

// _____________________________________________________________________________
XMLWriter::~XMLWriter() {
  _wr.closeTags();
}

// _____________________________________________________________________________
void XMLWriter::writeWay(const OsmWay& w) {
    _wr.openTag("way", "id", std::to_string(w.id));
    for (osmid nid : w.nodes) {
      _wr.openTag("nd", "ref", std::to_string(nid));
      _wr.closeTag();
    }
    for (const auto& kv : w.attrs) {
      std::map<std::string, std::string> attrs;
      attrs["k"] = kv.first;
      attrs["v"] = _source->decode(kv.second);
      _wr.openTag("tag", attrs);
      _wr.closeTag();
    }
    _wr.closeTag();
}

// _____________________________________________________________________________
void XMLWriter::writeNode(const OsmNode& nd) {
  _wr.openTag("node", {{"id", std::to_string(nd.id)},
                       {"lat", util::formatFloat(nd.lat, 7)},
                       {"lon", util::formatFloat(nd.lng, 7)}});
  for (const auto& kv : nd.attrs) {
    _wr.openTag("tag", {{"k", kv.first}, {"v", _source->decode(kv.second)}});
    _wr.closeTag();
  }
  _wr.closeTag();
}

// _____________________________________________________________________________
void XMLWriter::writeRel(const OsmRel& rel, const OsmIdList& nodes,
                         const OsmIdList& ways,
                         std::vector<const char*> nodeRoles,
                         std::vector<const char*> wayRoles) {
  _wr.openTag("relation", "id", std::to_string(rel.id));

  for (size_t j = 0; j < nodes.size(); j++) {
    osmid nid = nodes[j];
    std::map<std::string, std::string> attrs;
    attrs["type"] = "node";
    if (strlen(nodeRoles[j])) attrs["role"] = nodeRoles[j];
    attrs["ref"] = std::to_string(nid);
    _wr.openTag("member", attrs);
    _wr.closeTag();
  }

  for (size_t j = 0; j < ways.size(); j++) {
    osmid wid = ways[j];
    std::map<std::string, std::string> attrs;
    attrs["type"] = "way";
    if (strlen(wayRoles[j])) attrs["role"] = wayRoles[j];
    attrs["ref"] = std::to_string(wid);
    _wr.openTag("member", attrs);
    _wr.closeTag();
  }

  for (const auto& kv : rel.attrs) {
    std::map<std::string, std::string> attrs = {
        {"k", kv.first}, {"v", _source->decode(kv.second)}};
    _wr.openTag("tag", attrs);
    _wr.closeTag();
  }

  _wr.closeTag();
}

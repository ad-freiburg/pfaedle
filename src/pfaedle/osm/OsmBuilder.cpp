// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <float.h>

#include <algorithm>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include "pfaedle/Def.h"
#include "pfaedle/_config.h"
#include "pfaedle/osm/BBoxIdx.h"
#include "pfaedle/osm/Osm.h"
#include "pfaedle/osm/OsmBuilder.h"
#include "pfaedle/osm/OsmFilter.h"
#include "pfaedle/osm/Restrictor.h"
#include "util/Misc.h"
#include "util/Nullable.h"
#include "util/log/Log.h"
#include "pfxml/pfxml.h"

using ad::cppgtfs::gtfs::Stop;
using pfaedle::osm::BlockSearch;
using pfaedle::osm::EdgeGrid;
using pfaedle::osm::EqSearch;
using pfaedle::osm::NodeGrid;
using pfaedle::osm::OsmBuilder;
using pfaedle::osm::OsmNode;
using pfaedle::osm::OsmRel;
using pfaedle::osm::OsmWay;
using pfaedle::trgraph::Component;
using pfaedle::trgraph::Edge;
using pfaedle::trgraph::EdgePL;
using pfaedle::trgraph::Graph;
using pfaedle::trgraph::Node;
using pfaedle::trgraph::NodePL;
using pfaedle::trgraph::Normalizer;
using pfaedle::trgraph::StatInfo;
using pfaedle::trgraph::TransitEdgeLine;
using util::Nullable;
using util::geo::Box;
using util::geo::M_PER_DEG;

// _____________________________________________________________________________
bool EqSearch::operator()(const Node* cand, const StatInfo* si) const {
  return cand->pl().getSI() && cand->pl().getSI()->simi(si) > minSimi;
}

// _____________________________________________________________________________
OsmBuilder::OsmBuilder() {}

// _____________________________________________________________________________
void OsmBuilder::read(const std::string& path, const OsmReadOpts& opts,
                      Graph* g, const BBoxIdx& bbox, double gridSize,
                      Restrictor* res) {
  if (!bbox.size()) return;

  LOG(INFO) << "Reading OSM file " << path << " ... ";

  NodeSet orphanStations;
  EdgTracks eTracks;
  {
    OsmIdSet bboxNodes, noHupNodes;

    NIdMap nodes;
    NIdMultMap multNodes;
    RelLst intmRels;
    RelMap nodeRels, wayRels;

    Restrictions rawRests;

    AttrKeySet attrKeys[3] = {};
    getKeptAttrKeys(opts, attrKeys);

    OsmFilter filter(opts);

    pfxml::file xml(path);

    // we do four passes of the file here to be as memory creedy as possible:
    // - the first pass collects all node IDs which are
    //    * inside the given bounding box
    //    * (TODO: maybe more filtering?)
    //   these nodes are stored on the HD via OsmIdSet (which implements a
    //   simple bloom filter / base 256 encoded id store
    // - the second pass collects filtered relations
    // - the third pass collects filtered ways which contain one of the nodes
    //   from pass 1
    // - the forth pass collects filtered nodes which were
    //    * collected as node ids in pass 1
    //    * match the filter criteria
    //    * have been used in a way in pass 3

    LOG(DEBUG) << "Reading bounding box nodes...";
    skipUntil(&xml, "node");
    pfxml::parser_state nodeBeg = xml.state();
    pfxml::parser_state edgesBeg =
        readBBoxNds(&xml, &bboxNodes, &noHupNodes, filter, bbox);

    LOG(DEBUG) << "Reading relations...";
    skipUntil(&xml, "relation");
    readRels(&xml, &intmRels, &nodeRels, &wayRels, filter, attrKeys[2],
             &rawRests);

    LOG(DEBUG) << "Reading edges...";
    xml.set_state(edgesBeg);
    readEdges(&xml, g, intmRels, wayRels, filter, bboxNodes, &nodes, &multNodes,
              noHupNodes, attrKeys[1], rawRests, res, intmRels.flat, &eTracks,
              opts);

    LOG(DEBUG) << "Reading kept nodes...";
    xml.set_state(nodeBeg);
    readNodes(&xml, g, intmRels, nodeRels, filter, bboxNodes, &nodes,
              &multNodes, &orphanStations, attrKeys[0], intmRels.flat, opts);
  }

  LOG(DEBUG) << "OSM ID set lookups: " << osm::OsmIdSet::LOOKUPS
             << ", file lookups: " << osm::OsmIdSet::FLOOKUPS;

  LOG(DEBUG) << "Applying edge track numbers...";
  writeEdgeTracks(eTracks);
  eTracks.clear();

  {
    LOG(DEBUG) << "Fixing gaps...";
    NodeGrid ng = buildNodeIdx(g, gridSize, bbox.getFullBox(), false);
    LOG(DEBUG) << "Grid size of " << ng.getXWidth() << "x" << ng.getYHeight();
    fixGaps(g, &ng);
  }

  LOG(DEBUG) << "Snapping stations...";
  snapStats(opts, g, bbox, gridSize, res, orphanStations);

  LOG(DEBUG) << "Collapsing edges...";
  collapseEdges(g);

  LOG(DEBUG) << "Writing edge geoms...";
  writeGeoms(g, opts);

  LOG(DEBUG) << "Deleting orphan nodes...";
  deleteOrphNds(g, opts);

  LOG(DEBUG) << "Writing graph components...";
  // the restrictor is needed here to prevent connections in the graph
  // which are not possible in reality
  uint32_t comps = writeComps(g, opts);

  LOG(DEBUG) << "Simplifying geometries...";
  simplifyGeoms(g);

  LOG(DEBUG) << "Writing other-direction edges...";
  writeODirEdgs(g, res);

  LOG(DEBUG) << "Write wrong-direction costs...";
  writeOneWayPens(g, opts);

  if (opts.noLinesPunishFact != 1.0) {
    LOG(DEBUG) << "Write no-line pens...";
    writeNoLinePens(g, opts);
  }

  LOG(DEBUG) << "Write dummy node self-edges...";
  writeSelfEdgs(g);

  size_t numEdges = 0;

  for (auto* n : g->getNds()) {
    numEdges += n->getAdjListOut().size();
  }

  LOG(DEBUG) << "Graph has " << g->getNds().size() << " nodes, " << numEdges
             << " edges and " << comps
             << " connected component(s) with more than 1 node";
  LOG(DEBUG) << _lines.size() << " transit lines have been read.";
}

// _____________________________________________________________________________
void OsmBuilder::osmfilterRuleWrite(std::ostream* out,
                                    const std::vector<OsmReadOpts>& opts,
                                    const BBoxIdx& latLngBox) const {
  UNUSED(latLngBox);
  OsmIdSet bboxNodes, noHupNodes;
  MultAttrMap emptyF;

  RelLst rels;
  OsmIdList ways;
  RelMap nodeRels, wayRels;

  NIdMap nodes;

  OsmFilter filter;

  AttrKeySet attrKeys[3] = {};

  for (const OsmReadOpts& o : opts) {
    filter = filter.merge(OsmFilter(o.keepFilter, o.dropFilter));
    getKeptAttrKeys(o, attrKeys);
  }

  *out << "--keep=\n";

  for (auto r : filter.getKeepRules()) {
    for (auto val : r.second) {
      *out << r.first << "=";
      if (val.first != "*") *out << val.first;
      *out << "\n";
    }
  }

  *out << "\n";

  *out << "--keep-tags=\n";
  *out << "all\n";

  for (const auto& keys : attrKeys) {
    for (auto val : keys) {
      *out << val << "=\n";
    }
  }
}

// _____________________________________________________________________________
void OsmBuilder::overpassQryWrite(std::ostream* out,
                                  const std::vector<OsmReadOpts>& opts,
                                  const BBoxIdx& latLngBox) const {
  OsmIdSet bboxNodes, noHupNodes;
  MultAttrMap emptyF;

  RelLst rels;
  OsmIdList ways;
  RelMap nodeRels, wayRels;

  // TODO(patrick): not needed here!
  Restrictions rests;

  NIdMap nodes;

  // always empty
  NIdMultMap multNodes;
  util::xml::XmlWriter wr(out, true, 4);

  *out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
  wr.openComment();
  wr.writeText(" - written by pfaedle -");
  wr.closeTag();
  wr.openTag("osm-script");

  OsmFilter filter;

  for (const OsmReadOpts& o : opts) {
    filter = filter.merge(OsmFilter(o.keepFilter, o.dropFilter));
  }

  wr.openTag("union");
  size_t c = 0;
  for (auto box : latLngBox.getLeafs()) {
    if (box.getLowerLeft().getX() > box.getUpperRight().getX()) continue;
    c++;
    wr.openComment();
    wr.writeText(std::string("Bounding box #") + std::to_string(c) + " (" +
                 std::to_string(box.getLowerLeft().getY()) + ", " +
                 std::to_string(box.getLowerLeft().getX()) + ", " +
                 std::to_string(box.getUpperRight().getY()) + ", " +
                 std::to_string(box.getUpperRight().getX()) + ")");
    wr.closeTag();
    for (auto t : std::vector<std::string>{"way", "node", "relation"}) {
      for (auto r : filter.getKeepRules()) {
        for (auto val : r.second) {
          if (t == "way" && (val.second & OsmFilter::WAY)) continue;
          if (t == "relation" && (val.second & OsmFilter::REL)) continue;
          if (t == "node" && (val.second & OsmFilter::NODE)) continue;

          wr.openTag("query", {{"type", t}});
          if (val.first == "*")
            wr.openTag("has-kv", {{"k", r.first}});
          else
            wr.openTag("has-kv", {{"k", r.first}, {"v", val.first}});
          wr.closeTag();
          wr.openTag("bbox-query",
                     {{"s", std::to_string(box.getLowerLeft().getY())},
                      {"w", std::to_string(box.getLowerLeft().getX())},
                      {"n", std::to_string(box.getUpperRight().getY())},
                      {"e", std::to_string(box.getUpperRight().getX())}});
          wr.closeTag();
          wr.closeTag();
        }
      }
    }
  }

  wr.closeTag();

  wr.openTag("union");
  wr.openTag("item");
  wr.closeTag();
  wr.openTag("recurse", {{"type", "down"}});
  wr.closeTag();
  wr.closeTag();
  wr.openTag("print");

  wr.closeTags();
}

// _____________________________________________________________________________
void OsmBuilder::filterWrite(const std::string& in, const std::string& out,
                             const std::vector<OsmReadOpts>& opts,
                             const BBoxIdx& box) {
  OsmIdSet bboxNodes, noHupNodes;
  MultAttrMap emptyF;

  RelLst rels;
  OsmIdList ways;
  RelMap nodeRels, wayRels;

  // TODO(patrick): not needed here!
  Restrictions rests;

  NIdMap nodes;

  // always empty
  NIdMultMap multNodes;

  pfxml::file xml(in);

  BBoxIdx latLngBox = box;

  if (latLngBox.size() == 0) {
    skipUntil(&xml, "bounds");

    const pfxml::tag& cur = xml.get();

    if (strcmp(cur.name, "bounds") != 0) {
      throw pfxml::parse_exc(
          std::string("Could not find required <bounds> tag"), in, 0, 0, 0);
    }

    if (!cur.attr("minlat")) {
      throw pfxml::parse_exc(
          std::string(
              "Could not find required attribute \"minlat\" for <bounds> tag"),
          in, 0, 0, 0);
    }
    if (!cur.attr("minlon")) {
      throw pfxml::parse_exc(
          std::string(
              "Could not find required attribute \"minlon\" for <bounds> tag"),
          in, 0, 0, 0);
    }
    if (!cur.attr("maxlat")) {
      throw pfxml::parse_exc(
          std::string(
              "Could not find required attribute \"maxlat\" for <bounds> tag"),
          in, 0, 0, 0);
    }
    if (!cur.attr("maxlon")) {
      throw pfxml::parse_exc(
          std::string(
              "Could not find required attribute \"maxlon\" for <bounds> tag"),
          in, 0, 0, 0);
    }

    double minlat = atof(cur.attr("minlat"));
    double minlon = atof(cur.attr("minlon"));
    double maxlat = atof(cur.attr("maxlat"));
    double maxlon = atof(cur.attr("maxlon"));

    latLngBox.add(Box<double>({minlon, minlat}, {maxlon, maxlat}));
  }

  util::xml::XmlWriter wr(out, false, 0);

  wr.put("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
  wr.openTag("osm", {{"version", "0.6"},
                     {"generator", std::string("pfaedle/") + VERSION_FULL}});
  wr.openTag(
      "bounds",
      {{"minlat", std::to_string(latLngBox.getFullBox().getLowerLeft().getY())},
       {"minlon", std::to_string(latLngBox.getFullBox().getLowerLeft().getX())},
       {"maxlat",
        std::to_string(latLngBox.getFullBox().getUpperRight().getY())},
       {"maxlon",
        std::to_string(latLngBox.getFullBox().getUpperRight().getX())}});
  wr.closeTag();

  OsmFilter filter;
  AttrKeySet attrKeys[3] = {};

  for (const OsmReadOpts& o : opts) {
    getKeptAttrKeys(o, attrKeys);
    filter = filter.merge(OsmFilter(o.keepFilter, o.dropFilter));
  }

  skipUntil(&xml, "node");
  pfxml::parser_state nodeBeg = xml.state();
  pfxml::parser_state edgesBeg =
      readBBoxNds(&xml, &bboxNodes, &noHupNodes, filter, latLngBox);

  skipUntil(&xml, "relation");
  readRels(&xml, &rels, &nodeRels, &wayRels, filter, attrKeys[2], &rests);

  xml.set_state(edgesBeg);
  readEdges(&xml, wayRels, filter, bboxNodes, attrKeys[1], &ways, &nodes,
            rels.flat);

  xml.set_state(nodeBeg);

  readWriteNds(&xml, &wr, nodeRels, filter, bboxNodes, &nodes, attrKeys[0],
               rels.flat);
  readWriteWays(&xml, &wr, &ways, attrKeys[1]);

  std::sort(ways.begin(), ways.end());
  skipUntil(&xml, "relation");
  readWriteRels(&xml, &wr, &ways, &nodes, filter, attrKeys[2]);

  wr.closeTags();
}

// _____________________________________________________________________________
void OsmBuilder::readWriteRels(pfxml::file* i, util::xml::XmlWriter* o,
                               OsmIdList* ways, NIdMap* nodes,
                               const OsmFilter& filter,
                               const AttrKeySet& keepAttrs) {
  OsmRel rel;
  while ((rel = nextRel(i, filter, keepAttrs)).id) {
    OsmIdList realNodes, realWays;
    std::vector<const char*> realNodeRoles, realWayRoles;

    for (size_t j = 0; j < rel.ways.size(); j++) {
      osmid wid = rel.ways[j];
      const auto& i = std::lower_bound(ways->begin(), ways->end(), wid);
      if (i != ways->end() && *i == wid) {
        realWays.push_back(wid);
        realWayRoles.push_back(rel.wayRoles[j].c_str());
      }
    }

    for (size_t j = 0; j < rel.nodes.size(); j++) {
      osmid nid = rel.nodes[j];
      if (nodes->count(nid)) {
        realNodes.push_back(nid);
        realNodeRoles.push_back(rel.nodeRoles[j].c_str());
      }
    }

    if (realNodes.size() || realWays.size()) {
      o->openTag("relation", "id", std::to_string(rel.id));

      for (size_t j = 0; j < realNodes.size(); j++) {
        osmid nid = realNodes[j];
        std::map<std::string, std::string> attrs;
        attrs["type"] = "node";
        if (strlen(realNodeRoles[j])) attrs["role"] = realNodeRoles[j];
        attrs["ref"] = std::to_string(nid);
        o->openTag("member", attrs);
        o->closeTag();
      }

      for (size_t j = 0; j < realWays.size(); j++) {
        osmid wid = realWays[j];
        std::map<std::string, std::string> attrs;
        attrs["type"] = "way";
        if (strlen(realWayRoles[j])) attrs["role"] = realWayRoles[j];
        attrs["ref"] = std::to_string(wid);
        o->openTag("member", attrs);
        o->closeTag();
      }

      for (const auto& kv : rel.attrs) {
        std::map<std::string, std::string> attrs = {
            {"k", kv.first}, {"v", pfxml::file::decode(kv.second)}};
        o->openTag("tag", attrs);
        o->closeTag();
      }

      o->closeTag();
    }
  }
}

// _____________________________________________________________________________
void OsmBuilder::readWriteWays(pfxml::file* i, util::xml::XmlWriter* o,
                               OsmIdList* ways,
                               const AttrKeySet& keepAttrs) const {
  OsmWay w;
  NIdMultMap empty;
  for (auto wid : *ways) {
    w = nextWayWithId(i, wid, keepAttrs);
    assert(w.id);
    o->openTag("way", "id", std::to_string(wid));
    for (osmid nid : w.nodes) {
      o->openTag("nd", "ref", std::to_string(nid));
      o->closeTag();
    }
    for (const auto& kv : w.attrs) {
      std::map<std::string, std::string> attrs;
      attrs["k"] = kv.first;
      attrs["v"] = pfxml::file::decode(kv.second);
      o->openTag("tag", attrs);
      o->closeTag();
    }
    o->closeTag();
  }
}

// _____________________________________________________________________________
NodePL OsmBuilder::plFromGtfs(const Stop* s, const OsmReadOpts& ops) {
  NodePL ret({s->getLat(), s->getLng()},
             StatInfo(ops.statNormzer.norm(s->getName()),
                      ops.trackNormzer.norm(s->getPlatformCode())));

#ifdef PFAEDLE_STATION_IDS
  // debug feature, store station id from GTFS
  ret.getSI()->setId(s->getId());
#endif

  if (s->getParentStation()) {
    ret.getSI()->addAltName(
        ops.statNormzer.norm(s->getParentStation()->getName()));
  }

  return ret;
}

// _____________________________________________________________________________
pfxml::parser_state OsmBuilder::readBBoxNds(pfxml::file* xml, OsmIdSet* nodes,
                                            OsmIdSet* nohupNodes,
                                            const OsmFilter& filter,
                                            const BBoxIdx& bbox) const {
  bool inNodeBlock = false;
  uint64_t curId = 0;

  do {
    const pfxml::tag& cur = xml->get();

    if (inNodeBlock && xml->level() == 3 && curId &&
        strcmp(cur.name, "tag") == 0) {
      if (filter.nohup(cur.attr("k"), cur.attr("v"))) {
        nohupNodes->add(curId);
      }
    }

    if (xml->level() != 2) continue;
    if (!inNodeBlock && strcmp(cur.name, "node") == 0) inNodeBlock = true;

    if (inNodeBlock) {
      // block ended
      if (strcmp(cur.name, "node")) return xml->state();
      double y = util::atof(cur.attr("lat"), 7);
      double x = util::atof(cur.attr("lon"), 7);

      curId = util::atoul(cur.attr("id"));

      if (bbox.contains(Point<double>(x, y))) {
        nodes->add(curId);
      } else {
        nodes->nadd(curId);
      }
    }
  } while (xml->next());

  return xml->state();
}

// _____________________________________________________________________________
OsmWay OsmBuilder::nextWayWithId(pfxml::file* xml, osmid wid,
                                 const AttrKeySet& keepAttrs) const {
  OsmWay w;

  do {
    const pfxml::tag& cur = xml->get();
    if (xml->level() == 2 || xml->level() == 0) {
      if (w.id || strcmp(cur.name, "way")) return w;

      osmid id = util::atoul(cur.attr("id"));
      if (id == wid) w.id = id;
    }

    if (w.id && xml->level() == 3) {
      if (strcmp(cur.name, "nd") == 0) {
        w.nodes.push_back(util::atoul(cur.attr("ref")));
      } else if (strcmp(cur.name, "tag") == 0) {
        if (keepAttrs.count(cur.attr("k")))
          w.attrs[cur.attr("k")] = cur.attr("v");
      }
    }
  } while (xml->next());

  if (w.id) return w;

  return OsmWay();
}

// _____________________________________________________________________________
void OsmBuilder::skipUntil(pfxml::file* xml, const std::string& s) const {
  while (xml->next() && strcmp(xml->get().name, s.c_str())) {
  }
}

// _____________________________________________________________________________
bool OsmBuilder::relKeep(osmid id, const RelMap& rels,
                         const FlatRels& fl) const {
  auto it = rels.find(id);

  if (it == rels.end()) return false;

  for (osmid relId : it->second) {
    // as soon as any of this entities relations is not flat, return true
    if (!fl.count(relId)) return true;
  }

  return false;
}

// _____________________________________________________________________________
OsmWay OsmBuilder::nextWay(pfxml::file* xml, const RelMap& wayRels,
                           const OsmFilter& filter, const OsmIdSet& bBoxNodes,
                           const AttrKeySet& keepAttrs,
                           const FlatRels& fl) const {
  OsmWay w;

  do {
    const pfxml::tag& cur = xml->get();
    if (xml->level() == 2 || xml->level() == 0) {
      if (keepWay(w, wayRels, filter, bBoxNodes, fl)) return w;
      if (strcmp(cur.name, "way")) return OsmWay();

      w.id = util::atoul(cur.attr("id"));
      w.nodes.clear();
      w.attrs.clear();
    }

    if (w.id && xml->level() == 3) {
      if (strcmp(cur.name, "nd") == 0) {
        osmid nid = util::atoul(cur.attr("ref"));
        w.nodes.push_back(nid);
      } else if (strcmp(cur.name, "tag") == 0) {
        if (keepAttrs.count(cur.attr("k")))
          w.attrs[cur.attr("k")] = cur.attr("v");
      }
    }
  } while (xml->next());

  if (keepWay(w, wayRels, filter, bBoxNodes, fl)) return w;
  return OsmWay();
}

// _____________________________________________________________________________
bool OsmBuilder::keepWay(const OsmWay& w, const RelMap& wayRels,
                         const OsmFilter& filter, const OsmIdSet& bBoxNodes,
                         const FlatRels& fl) const {
  if (w.id && w.nodes.size() > 1 &&
      (relKeep(w.id, wayRels, fl) || filter.keep(w.attrs, OsmFilter::WAY)) &&
      !filter.drop(w.attrs, OsmFilter::WAY)) {
    for (osmid nid : w.nodes) {
      if (bBoxNodes.has(nid)) {
        return true;
      }
    }
  }

  return false;
}

// _____________________________________________________________________________
void OsmBuilder::readEdges(pfxml::file* xml, const RelMap& wayRels,
                           const OsmFilter& filter, const OsmIdSet& bBoxNodes,
                           const AttrKeySet& keepAttrs, OsmIdList* ret,
                           NIdMap* nodes, const FlatRels& flat) {
  OsmWay w;
  while ((w = nextWay(xml, wayRels, filter, bBoxNodes, keepAttrs, flat)).id) {
    ret->push_back(w.id);
    for (auto n : w.nodes) {
      (*nodes)[n] = 0;
    }
  }
}

// _____________________________________________________________________________
void OsmBuilder::readEdges(pfxml::file* xml, Graph* g, const RelLst& rels,
                           const RelMap& wayRels, const OsmFilter& filter,
                           const OsmIdSet& bBoxNodes, NIdMap* nodes,
                           NIdMultMap* multiNodes, const OsmIdSet& noHupNodes,
                           const AttrKeySet& keepAttrs,
                           const Restrictions& rawRests, Restrictor* restor,
                           const FlatRels& fl, EdgTracks* eTracks,
                           const OsmReadOpts& opts) {
  OsmWay w;
  while ((w = nextWay(xml, wayRels, filter, bBoxNodes, keepAttrs, fl)).id) {
    Node* last = 0;
    std::vector<TransitEdgeLine*> lines;
    if (wayRels.count(w.id)) {
      lines = getLines(wayRels.find(w.id)->second, rels, opts);
    }
    std::string track =
        getAttrByFirstMatch(opts.edgePlatformRules, w.id, w.attrs, wayRels,
                            rels, opts.trackNormzer);

    osmid lastnid = 0;
    for (osmid nid : w.nodes) {
      Node* n = 0;
      if (noHupNodes.has(nid)) {
        n = g->addNd();
        (*multiNodes)[nid].insert(n);
      } else if (!nodes->count(nid)) {
        if (!bBoxNodes.has(nid)) continue;
        n = g->addNd();
        (*nodes)[nid] = n;
      } else {
        n = (*nodes)[nid];
      }

      if (last) {
        auto e = g->addEdg(last, n, EdgePL());
        if (!e) continue;

        processRestr(nid, w.id, rawRests, e, n, restor);
        processRestr(lastnid, w.id, rawRests, e, last, restor);

        e->pl().addLines(lines);
        e->pl().setLvl(filter.level(w.attrs));
        if (!track.empty()) (*eTracks)[e] = track;

        if (filter.oneway(w.attrs)) e->pl().setOneWay(1);
        if (filter.onewayrev(w.attrs)) e->pl().setOneWay(2);
      }
      lastnid = nid;
      last = n;
    }
  }
}

// _____________________________________________________________________________
void OsmBuilder::processRestr(osmid nid, osmid wid,
                              const Restrictions& rawRests, Edge* e, Node* n,
                              Restrictor* restor) const {
  if (rawRests.pos.count(nid)) {
    for (const auto& kv : rawRests.pos.find(nid)->second) {
      if (kv.eFrom == wid) {
        e->pl().setRestricted();
        restor->add(e, kv.eTo, n, true);
      } else if (kv.eTo == wid) {
        e->pl().setRestricted();
        restor->relax(wid, n, e);
      }
    }
  }

  if (rawRests.neg.count(nid)) {
    for (const auto& kv : rawRests.neg.find(nid)->second) {
      if (kv.eFrom == wid) {
        e->pl().setRestricted();
        restor->add(e, kv.eTo, n, false);
      } else if (kv.eTo == wid) {
        e->pl().setRestricted();
        restor->relax(wid, n, e);
      }
    }
  }
}

// _____________________________________________________________________________
OsmNode OsmBuilder::nextNode(pfxml::file* xml, NIdMap* nodes,
                             NIdMultMap* multNodes, const RelMap& nodeRels,
                             const OsmFilter& filter, const OsmIdSet& bBoxNodes,
                             const AttrKeySet& keepAttrs,
                             const FlatRels& fl) const {
  OsmNode n;

  do {
    const pfxml::tag& cur = xml->get();
    if (xml->level() == 2 || xml->level() == 0) {
      if (keepNode(n, *nodes, *multNodes, nodeRels, bBoxNodes, filter, fl))
        return n;
      // block ended
      if (strcmp(cur.name, "node")) return OsmNode();

      n.attrs.clear();
      n.lat = util::atof(cur.attr("lat"), 7);
      n.lng = util::atof(cur.attr("lon"), 7);
      n.id = util::atoul(cur.attr("id"));
    }

    if (xml->level() == 3 && n.id && strcmp(cur.name, "tag") == 0) {
      if (keepAttrs.count(cur.attr("k")))
        n.attrs[cur.attr("k")] = cur.attr("v");
    }
  } while (xml->next());

  if (keepNode(n, *nodes, *multNodes, nodeRels, bBoxNodes, filter, fl))
    return n;
  return OsmNode();
}

// _____________________________________________________________________________
bool OsmBuilder::keepNode(const OsmNode& n, const NIdMap& nodes,
                          const NIdMultMap& multNodes, const RelMap& nodeRels,
                          const OsmIdSet& bBoxNodes, const OsmFilter& filter,
                          const FlatRels& fl) const {
  if (n.id &&
      (nodes.count(n.id) || multNodes.count(n.id) ||
       relKeep(n.id, nodeRels, fl) || filter.keep(n.attrs, OsmFilter::NODE)) &&
      (nodes.count(n.id) || bBoxNodes.has(n.id)) &&
      (nodes.count(n.id) || multNodes.count(n.id) ||
       !filter.drop(n.attrs, OsmFilter::NODE))) {
    return true;
  }

  return false;
}

// _____________________________________________________________________________
void OsmBuilder::readWriteNds(pfxml::file* i, util::xml::XmlWriter* o,
                              const RelMap& nRels, const OsmFilter& filter,
                              const OsmIdSet& bBoxNds, NIdMap* nds,
                              const AttrKeySet& keepAttrs,
                              const FlatRels& f) const {
  OsmNode nd;
  NIdMultMap empt;
  while (
      (nd = nextNode(i, nds, &empt, nRels, filter, bBoxNds, keepAttrs, f)).id) {
    (*nds)[nd.id] = 0;
    o->openTag("node", {{"id", std::to_string(nd.id)},
                        {"lat", std::to_string(nd.lat)},
                        {"lon", std::to_string(nd.lng)}});
    for (const auto& kv : nd.attrs) {
      o->openTag("tag",
                 {{"k", kv.first}, {"v", pfxml::file::decode(kv.second)}});
      o->closeTag();
    }
    o->closeTag();
  }
}

// _____________________________________________________________________________
void OsmBuilder::readNodes(pfxml::file* xml, Graph* g, const RelLst& rels,
                           const RelMap& nodeRels, const OsmFilter& filter,
                           const OsmIdSet& bBoxNodes, NIdMap* nodes,
                           NIdMultMap* multNodes, NodeSet* orphanStations,
                           const AttrKeySet& keepAttrs, const FlatRels& fl,
                           const OsmReadOpts& opts) const {
  OsmNode nd;
  while ((nd = nextNode(xml, nodes, multNodes, nodeRels, filter, bBoxNodes,
                        keepAttrs, fl))
             .id) {
    Node* n = 0;
    POINT pos = {nd.lng, nd.lat};
    if (nodes->count(nd.id)) {
      n = (*nodes)[nd.id];
      n->pl().setGeom(pos);
      if (filter.station(nd.attrs)) {
        auto si = getStatInfo(nd.id, nd.attrs, nodeRels, rels, opts);
        if (!si.isNull()) n->pl().setSI(si);
      } else if (filter.blocker(nd.attrs)) {
        n->pl().setBlocker();
      } else if (filter.turnCycle(nd.attrs)) {
        n->pl().setTurnCycle();
      }
    } else if ((*multNodes).count(nd.id)) {
      for (auto* n : (*multNodes)[nd.id]) {
        n->pl().setGeom(pos);
        if (filter.station(nd.attrs)) {
          auto si = getStatInfo(nd.id, nd.attrs, nodeRels, rels, opts);
          if (!si.isNull()) n->pl().setSI(si);
        } else if (filter.blocker(nd.attrs)) {
          n->pl().setBlocker();
        } else if (filter.turnCycle(nd.attrs)) {
          n->pl().setTurnCycle();
        }
      }
    } else {
      // these are nodes without any connected edges
      if (filter.station(nd.attrs)) {
        auto tmp = g->addNd(NodePL(pos));
        auto si = getStatInfo(nd.id, nd.attrs, nodeRels, rels, opts);
        if (!si.isNull()) tmp->pl().setSI(si);
        if (tmp->pl().getSI()) {
          orphanStations->insert(tmp);
        }
      }
    }
  }
}

// _____________________________________________________________________________
OsmRel OsmBuilder::nextRel(pfxml::file* xml, const OsmFilter& filter,
                           const AttrKeySet& keepAttrs) const {
  OsmRel rel;

  do {
    const pfxml::tag& cur = xml->get();
    if (xml->level() == 2 || xml->level() == 0) {
      uint64_t keepFlags = 0;
      uint64_t dropFlags = 0;
      if (rel.id && rel.attrs.size() &&
          (keepFlags = filter.keep(rel.attrs, OsmFilter::REL)) &&
          !(dropFlags = filter.drop(rel.attrs, OsmFilter::REL))) {
        rel.keepFlags = keepFlags;
        rel.dropFlags = dropFlags;
        return rel;
      }

      // block ended
      if (strcmp(cur.name, "relation")) return OsmRel();

      rel.attrs.clear();
      rel.nodes.clear();
      rel.ways.clear();
      rel.nodeRoles.clear();
      rel.wayRoles.clear();
      rel.keepFlags = 0;
      rel.dropFlags = 0;
      rel.id = util::atoul(cur.attr("id"));
    }

    if (xml->level() == 3 && rel.id) {
      if (strcmp(cur.name, "member") == 0) {
        if (strcmp(cur.attr("type"), "node") == 0) {
          osmid id = util::atoul(cur.attr("ref"));
          // TODO(patrick): no need to push IDs that have been filtered out by
          // the bounding box!!!!
          rel.nodes.push_back(id);
          if (cur.attr("role")) {
            rel.nodeRoles.push_back(cur.attr("role"));
          } else {
            rel.nodeRoles.push_back("");
          }
        }
        if (strcmp(cur.attr("type"), "way") == 0) {
          osmid id = util::atoul(cur.attr("ref"));
          rel.ways.push_back(id);
          if (cur.attr("role")) {
            rel.wayRoles.push_back(cur.attr("role"));
          } else {
            rel.wayRoles.push_back("");
          }
        }
      } else if (strcmp(cur.name, "tag") == 0) {
        if (keepAttrs.count(cur.attr("k")))
          rel.attrs[cur.attr("k")] = cur.attr("v");
      }
    }
  } while (xml->next());

  // dont forget last relation
  uint64_t keepFlags = 0;
  uint64_t dropFlags = 0;
  if (rel.id && rel.attrs.size() &&
      (keepFlags = filter.keep(rel.attrs, OsmFilter::REL)) &&
      !(dropFlags = filter.drop(rel.attrs, OsmFilter::REL))) {
    rel.keepFlags = keepFlags;
    rel.dropFlags = dropFlags;
    return rel;
  }

  return OsmRel();
}

// _____________________________________________________________________________
void OsmBuilder::readRels(pfxml::file* xml, RelLst* rels, RelMap* nodeRels,
                          RelMap* wayRels, const OsmFilter& filter,
                          const AttrKeySet& keepAttrs,
                          Restrictions* rests) const {
  OsmRel rel;
  while ((rel = nextRel(xml, filter, keepAttrs)).id) {
    rels->rels.push_back(rel.attrs);
    if (rel.keepFlags & osm::REL_NO_DOWN) {
      rels->flat.insert(rels->rels.size() - 1);
    }
    for (osmid id : rel.nodes) (*nodeRels)[id].push_back(rels->rels.size() - 1);
    for (osmid id : rel.ways) (*wayRels)[id].push_back(rels->rels.size() - 1);

    // TODO(patrick): this is not needed for the filtering - remove it here!
    readRestr(rel, rests, filter);
  }
}

// _____________________________________________________________________________
void OsmBuilder::readRestr(const OsmRel& rel, Restrictions* rests,
                           const OsmFilter& filter) const {
  if (!rel.attrs.count("type")) return;
  if (rel.attrs.find("type")->second != "restriction") return;

  bool pos = filter.posRestr(rel.attrs);
  bool neg = filter.negRestr(rel.attrs);

  if (!pos && !neg) return;

  osmid from = 0;
  osmid to = 0;
  osmid via = 0;

  for (size_t i = 0; i < rel.ways.size(); i++) {
    if (rel.wayRoles[i] == "from") {
      if (from) return;  // only one from member supported
      from = rel.ways[i];
    }
    if (rel.wayRoles[i] == "to") {
      if (to) return;  // only one to member supported
      to = rel.ways[i];
    }
  }

  for (size_t i = 0; i < rel.nodes.size(); i++) {
    if (rel.nodeRoles[i] == "via") {
      via = rel.nodes[i];
      break;
    }
  }

  if (from && to && via) {
    if (pos)
      rests->pos[via].push_back(Restriction{from, to});
    else if (neg)
      rests->neg[via].push_back(Restriction{from, to});
  }
}

// _____________________________________________________________________________
std::string OsmBuilder::getAttrByFirstMatch(const DeepAttrLst& rule, osmid id,
                                            const AttrMap& am,
                                            const RelMap& entRels,
                                            const RelLst& rels,
                                            const Normalizer& normzer) const {
  std::string ret;
  for (const auto& s : rule) {
    ret = normzer.norm(pfxml::file::decode(getAttr(s, id, am, entRels, rels)));
    if (!ret.empty()) return ret;
  }

  return ret;
}

// _____________________________________________________________________________
std::vector<std::string> OsmBuilder::getAttrMatchRanked(
    const DeepAttrLst& rule, osmid id, const AttrMap& am, const RelMap& entRels,
    const RelLst& rels, const Normalizer& norm) const {
  std::vector<std::string> ret;
  for (const auto& s : rule) {
    std::string tmp =
        norm.norm(pfxml::file::decode(getAttr(s, id, am, entRels, rels)));
    if (!tmp.empty()) ret.push_back(tmp);
  }

  return ret;
}

// _____________________________________________________________________________
std::string OsmBuilder::getAttr(const DeepAttrRule& s, osmid id,
                                const AttrMap& am, const RelMap& entRels,
                                const RelLst& rels) const {
  if (s.relRule.kv.first.empty()) {
    if (am.find(s.attr) != am.end()) {
      return am.find(s.attr)->second;
    }
  } else {
    if (entRels.count(id)) {
      for (const auto& relId : entRels.find(id)->second) {
        if (OsmFilter::contained(rels.rels[relId], s.relRule.kv)) {
          if (rels.rels[relId].count(s.attr)) {
            return rels.rels[relId].find(s.attr)->second;
          }
        }
      }
    }
  }
  return "";
}

// _____________________________________________________________________________
Nullable<StatInfo> OsmBuilder::getStatInfo(osmid nid, const AttrMap& m,
                                           const RelMap& nodeRels,
                                           const RelLst& rels,
                                           const OsmReadOpts& ops) const {
  std::string platform;
  std::vector<std::string> names;

  names = getAttrMatchRanked(ops.statAttrRules.nameRule, nid, m, nodeRels, rels,
                             ops.statNormzer);
  platform = getAttrByFirstMatch(ops.statAttrRules.platformRule, nid, m,
                                 nodeRels, rels, ops.trackNormzer);

  if (!names.size()) return Nullable<StatInfo>();

  auto ret = StatInfo(names[0], platform);

#ifdef PFAEDLE_STATION_IDS
  ret.setId(getAttrByFirstMatch(ops.statAttrRules.idRule, nid, m, nodeRels,
                                rels, ops.idNormzer));
#endif

  for (size_t i = 1; i < names.size(); i++) ret.addAltName(names[i]);

  return ret;
}

// _____________________________________________________________________________
double OsmBuilder::dist(const Node* a, const Node* b) {
  return util::geo::haversine(*(a->pl().getGeom()), *(b->pl().getGeom()));
}

// _____________________________________________________________________________
void OsmBuilder::writeGeoms(Graph* g, const OsmReadOpts& opts) {
  for (auto* n : g->getNds()) {
    for (auto* e : n->getAdjListOut()) {
      if (!e->pl().getGeom()) {
        e->pl().addPoint(*e->getFrom()->pl().getGeom());
        e->pl().addPoint(*e->getTo()->pl().getGeom());
      }

      e->pl().setCost(costToInt(e->pl().getLength() /
                                opts.levelDefSpeed[e->pl().lvl()]));
    }
  }
}

// _____________________________________________________________________________
void OsmBuilder::fixGaps(Graph* g, NodeGrid* ng) {
  double METER = 1;
  for (auto* n : g->getNds()) {
    if (n->getInDeg() + n->getOutDeg() == 1) {
      // get all nodes in distance
      std::set<Node*> ret;
      double distor = util::geo::latLngDistFactor(*n->pl().getGeom());
      ng->get(util::geo::pad(util::geo::getBoundingBox(*n->pl().getGeom()),
                             (METER / M_PER_DEG) / distor),
              &ret);
      for (auto* nb : ret) {
        if (nb != n && (nb->getInDeg() + nb->getOutDeg()) == 1 &&
            dist(nb, n) <= METER) {
          // special case: both nodes are non-stations, move
          // the end point nb to n and delete nb
          if (!nb->pl().getSI() && !n->pl().getSI()) {
            Node* otherN;
            if (nb->getOutDeg())
              otherN = (*nb->getAdjListOut().begin())->getOtherNd(nb);
            else
              otherN = (*nb->getAdjListIn().begin())->getOtherNd(nb);

            Edge* e;
            if (nb->getOutDeg())
              e = g->addEdg(otherN, n, (*nb->getAdjListOut().begin())->pl());
            else
              e = g->addEdg(otherN, n, (*nb->getAdjListIn().begin())->pl());
            if (e) {
              g->delNd(nb);
              ng->remove(nb);
            }
          } else {
            // if one of the nodes is a station, just add an edge between them
            if (nb->getOutDeg())
              g->addEdg(n, nb, (*nb->getAdjListOut().begin())->pl());
            else
              g->addEdg(n, nb, (*nb->getAdjListIn().begin())->pl());
          }
        }
      }
    }
  }
}

// _____________________________________________________________________________
EdgeGrid OsmBuilder::buildEdgeIdx(Graph* g, double size, const BOX& box) {
  EdgeGrid ret(size, size, box, false);
  for (auto* n : g->getNds()) {
    for (auto* e : n->getAdjListOut()) {
      auto llGeom =
          LINE{*e->getFrom()->pl().getGeom(), *e->getTo()->pl().getGeom()};
      ret.add(llGeom, e);
    }
  }
  return ret;
}

// _____________________________________________________________________________
NodeGrid OsmBuilder::buildNodeIdx(Graph* g, double size, const BOX& box,
                                  bool which) {
  NodeGrid ret(size, size, box, false);
  for (auto* n : g->getNds()) {
    // only orphan nodes
    if (!which && n->getInDeg() + n->getOutDeg() == 1)
      ret.add(*n->pl().getGeom(), n);
    // only station nodes
    else if (which && n->pl().getSI())
      ret.add(*n->pl().getGeom(), n);
  }
  return ret;
}

// _____________________________________________________________________________
Node* OsmBuilder::depthSearch(const Edge* e, const StatInfo* si, const POINT& p,
                              double maxD, int maxFullTurns, double minAngle,
                              const SearchFunc& sfunc) {
  // shortcuts
  double dFrom = haversine(*e->getFrom()->pl().getGeom(), p);
  double dTo = haversine(*e->getTo()->pl().getGeom(), p);
  if (dFrom > maxD && dTo > maxD) return 0;

  if (dFrom <= maxD && sfunc(e->getFrom(), si)) return e->getFrom();
  if (dTo <= maxD && sfunc(e->getTo(), si)) return e->getTo();

  NodeCandPQ pq;
  NodeSet closed;
  pq.push(NodeCand{dFrom, e->getFrom(), e, 0});
  if (e->getFrom() != e->getTo()) pq.push(NodeCand{dTo, e->getTo(), e, 0});

  while (!pq.empty()) {
    auto cur = pq.top();
    pq.pop();
    if (closed.count(cur.node)) continue;
    closed.insert(cur.node);

    for (size_t i = 0; i < cur.node->getInDeg() + cur.node->getOutDeg(); i++) {
      trgraph::Node* cand;
      trgraph::Edge* edg;

      if (i < cur.node->getInDeg()) {
        edg = cur.node->getAdjListIn()[i];
        cand = edg->getFrom();
      } else {
        edg = cur.node->getAdjListOut()[i - cur.node->getInDeg()];
        cand = edg->getTo();
      }

      if (cand == cur.node) continue;  // dont follow self edges

      int fullTurn = 0;

      if (cur.fromEdge && cur.node->getInDeg() + cur.node->getOutDeg() >
                              2) {  // only intersection angles
        const POINT& toP = *cand->pl().getGeom();
        const POINT& fromP =
            *cur.fromEdge->getOtherNd(cur.node)->pl().getGeom();
        const POINT& nodeP = *cur.node->pl().getGeom();

        if (util::geo::innerProd(nodeP, fromP, toP) < minAngle) fullTurn = 1;
      }

      double eLen = dist(edg->getFrom(), edg->getTo());

      if ((maxFullTurns < 0 || cur.fullTurns + fullTurn <= maxFullTurns) &&
          cur.dist + eLen < maxD && !closed.count(cand)) {
        if (sfunc(cand, si)) {
          return cand;
        } else {
          pq.push(
              NodeCand{cur.dist + eLen, cand, edg, cur.fullTurns + fullTurn});
        }
      }
    }
  }

  return 0;
}

// _____________________________________________________________________________
bool OsmBuilder::isBlocked(const Edge* e, const StatInfo* si, const POINT& p,
                           double maxD, int maxFullTurns, double minAngle) {
  return depthSearch(e, si, p, maxD, maxFullTurns, minAngle, BlockSearch());
}

// _____________________________________________________________________________
Node* OsmBuilder::eqStatReach(const Edge* e, const StatInfo* si, const POINT& p,
                              double maxD, int maxFullTurns, double minAngle) {
  return depthSearch(e, si, p, maxD, maxFullTurns, minAngle, EqSearch());
}

// _____________________________________________________________________________
void OsmBuilder::getEdgCands(const POINT& geom, EdgeCandPQ* ret, EdgeGrid* eg,
                             double d) {
  double distor = util::geo::latLngDistFactor(geom);
  std::set<Edge*> neighs;
  BOX box =
      util::geo::pad(util::geo::getBoundingBox(geom), (d / M_PER_DEG) / distor);
  eg->get(box, &neighs);

  for (auto* e : neighs) {
    double dist = util::geo::distToSegment(*e->getFrom()->pl().getGeom(),
                                           *e->getTo()->pl().getGeom(), geom);

    if (dist * distor * M_PER_DEG <= d) {
      ret->push(EdgeCand(-dist, e));
    }
  }
}

// _____________________________________________________________________________
void OsmBuilder::snapStation(Graph* g, NodePL* s, EdgeGrid* eg, NodeGrid* sng,
                             const OsmReadOpts& opts, Restrictor* restor,
                             double d) {
  assert(s->getSI());

  EdgeCandPQ pq;

  getEdgCands(*s->getGeom(), &pq, eg, d);

  while (!pq.empty()) {
    auto* e = pq.top().second;
    pq.pop();
    auto geom =
        util::geo::projectOn(*e->getFrom()->pl().getGeom(), *s->getGeom(),
                             *e->getTo()->pl().getGeom());

    Node* eq = 0;
    if (!(eq = eqStatReach(e, s->getSI(), geom, 2 * d, 0,
                           opts.maxAngleSnapReach))) {
      if (e->pl().lvl() > opts.maxSnapLevel) continue;
      if (isBlocked(e, s->getSI(), geom, opts.maxBlockDistance, 0,
                    opts.maxAngleSnapReach)) {
        continue;
      }

      // if the projected position is near (< 0.5 meters) the end point of this
      // way and the endpoint is not already a station, place the station there.
      if (!e->getFrom()->pl().getSI() &&
          haversine(geom, *e->getFrom()->pl().getGeom()) < .5) {
        e->getFrom()->pl().setSI(*s->getSI());
      } else if (!e->getTo()->pl().getSI() &&
                 haversine(geom, *e->getTo()->pl().getGeom()) < .5) {
        e->getTo()->pl().setSI(*s->getSI());
      } else {
        s->setGeom(geom);
        Node* n = g->addNd(*s);
        sng->add(geom, n);

        auto ne = g->addEdg(e->getFrom(), n, e->pl());
        ne->pl().setCost(costToInt(dist(e->getFrom(), n) /
                                   opts.levelDefSpeed[ne->pl().lvl()]));
        eg->add({*e->getFrom()->pl().getGeom(), *n->pl().getGeom()}, ne);

        auto nf = g->addEdg(n, e->getTo(), e->pl());
        nf->pl().setCost(costToInt(dist(n, e->getTo()) /
                                   opts.levelDefSpeed[nf->pl().lvl()]));
        eg->add({*n->pl().getGeom(), *e->getTo()->pl().getGeom()}, nf);

        // replace edge in restrictor
        restor->replaceEdge(e, ne, nf);

        g->delEdg(e->getFrom(), e->getTo());
        eg->remove(e);
      }
    } else {
      // if the snapped station is very near to the original OSM station
      // write additional info from this snap station to the equivalent stat
      if (haversine(*s->getGeom(), *eq->pl().getGeom()) < 5) {
        if (eq->pl().getSI()->getTrack().empty())
          eq->pl().getSI()->setTrack(s->getSI()->getTrack());
      }
    }
  }
}

// _____________________________________________________________________________
std::vector<TransitEdgeLine*> OsmBuilder::getLines(
    const std::vector<size_t>& edgeRels, const RelLst& rels,
    const OsmReadOpts& ops) {
  std::vector<TransitEdgeLine*> ret;
  for (size_t relId : edgeRels) {
    TransitEdgeLine* elp = 0;

    if (_relLines.count(relId)) {
      elp = _relLines[relId];
    } else {
      TransitEdgeLine el;
      el.color = ad::cppgtfs::gtfs::NO_COLOR;

      bool found = false;
      for (const auto& r : ops.relLinerules.sNameRule) {
        for (const auto& relAttr : rels.rels[relId]) {
          if (relAttr.first == r) {
            el.shortName =
                ops.lineNormzer.norm(pfxml::file::decode(relAttr.second));
            if (!el.shortName.empty()) found = true;
          }
        }
        if (found) break;
      }

      found = false;
      for (const auto& r : ops.relLinerules.fromNameRule) {
        for (const auto& relAttr : rels.rels[relId]) {
          if (relAttr.first == r) {
            el.fromStr =
                ops.statNormzer.norm(pfxml::file::decode(relAttr.second));
            if (!el.fromStr.empty()) found = true;
          }
        }
        if (found) break;
      }

      found = false;
      for (const auto& r : ops.relLinerules.toNameRule) {
        for (const auto& relAttr : rels.rels[relId]) {
          if (relAttr.first == r) {
            el.toStr =
                ops.statNormzer.norm(pfxml::file::decode(relAttr.second));
            if (!el.toStr.empty()) found = true;
          }
        }
        if (found) break;
      }

      found = false;
      for (const auto& r : ops.relLinerules.colorRule) {
        for (const auto& relAttr : rels.rels[relId]) {
          if (relAttr.first == r) {
            auto dec = pfxml::file::decode(relAttr.second);
            auto color = parseHexColor(dec);
            if (color == ad::cppgtfs::gtfs::NO_COLOR)
              color = parseHexColor(std::string("#") + dec);
            if (color != ad::cppgtfs::gtfs::NO_COLOR) {
              found = true;
              el.color = color;
            }
          }
        }
        if (found) break;
      }

      if (!el.shortName.size() && !el.fromStr.size() && !el.toStr.size())
        continue;

      if (_lines.count(el)) {
        elp = _lines[el];
        _relLines[relId] = elp;
      } else {
        elp = new TransitEdgeLine(el);
        _lines[el] = elp;
        _relLines[relId] = elp;
      }
    }
    ret.push_back(elp);
  }

  return ret;
}

// _____________________________________________________________________________
void OsmBuilder::getKeptAttrKeys(const OsmReadOpts& opts,
                                 AttrKeySet sets[3]) const {
  for (const auto& i : opts.keepFilter) {
    for (size_t j = 0; j < 3; j++) sets[j].insert(i.first);
  }

  for (const auto& i : opts.dropFilter) {
    for (size_t j = 0; j < 3; j++) sets[j].insert(i.first);
  }

  for (const auto& i : opts.noHupFilter) {
    sets[0].insert(i.first);
  }

  for (const auto& i : opts.oneWayFilter) {
    sets[1].insert(i.first);
  }

  for (const auto& i : opts.oneWayFilterRev) {
    sets[1].insert(i.first);
  }

  for (const auto& i : opts.twoWayFilter) {
    sets[1].insert(i.first);
  }

  for (const auto& i : opts.stationFilter) {
    sets[0].insert(i.first);
    sets[2].insert(i.first);
  }

  for (const auto& i : opts.stationBlockerFilter) {
    sets[0].insert(i.first);
  }

  for (const auto& i : opts.turnCycleFilter) {
    sets[0].insert(i.first);
  }

  for (uint8_t j = 0; j < 7; j++) {
    for (const auto& kv : *(opts.levelFilters + j)) {
      sets[1].insert(kv.first);
    }
  }

  // restriction system
  for (const auto& i : opts.restrPosRestr) {
    sets[2].insert(i.first);
  }
  for (const auto& i : opts.restrNegRestr) {
    sets[2].insert(i.first);
  }
  for (const auto& i : opts.noRestrFilter) {
    sets[2].insert(i.first);
  }

  sets[2].insert("from");
  sets[2].insert("via");
  sets[2].insert("to");

  sets[2].insert(opts.relLinerules.toNameRule.begin(),
                 opts.relLinerules.toNameRule.end());
  sets[2].insert(opts.relLinerules.fromNameRule.begin(),
                 opts.relLinerules.fromNameRule.end());
  sets[2].insert(opts.relLinerules.sNameRule.begin(),
                 opts.relLinerules.sNameRule.end());
  sets[2].insert(opts.relLinerules.colorRule.begin(),
                 opts.relLinerules.colorRule.end());

  for (const auto& i : opts.statAttrRules.nameRule) {
    if (i.relRule.kv.first.empty()) {
      sets[0].insert(i.attr);
    } else {
      sets[2].insert(i.relRule.kv.first);
      sets[2].insert(i.attr);
    }
  }

  for (const auto& i : opts.edgePlatformRules) {
    if (i.relRule.kv.first.empty()) {
      sets[1].insert(i.attr);
    } else {
      sets[2].insert(i.relRule.kv.first);
      sets[2].insert(i.attr);
    }
  }

  for (const auto& i : opts.statAttrRules.platformRule) {
    if (i.relRule.kv.first.empty()) {
      sets[0].insert(i.attr);
    } else {
      sets[2].insert(i.relRule.kv.first);
      sets[2].insert(i.attr);
    }
  }

  for (const auto& i : opts.statAttrRules.idRule) {
    if (i.relRule.kv.first.empty()) {
      sets[0].insert(i.attr);
    } else {
      sets[2].insert(i.relRule.kv.first);
      sets[2].insert(i.attr);
    }
  }
}

// _____________________________________________________________________________
void OsmBuilder::deleteOrphNds(Graph* g, const OsmReadOpts& opts) {
  UNUSED(opts);
  for (auto i = g->getNds().begin(); i != g->getNds().end();) {
    if ((*i)->getInDeg() + (*i)->getOutDeg() != 0 || (*i)->pl().getSI()) {
      ++i;
      continue;
    }

    i = g->delNd(*i);
  }
}

// _____________________________________________________________________________
bool OsmBuilder::edgesSim(const Edge* a, const Edge* b) {
  if (static_cast<bool>(a->pl().oneWay()) ^ static_cast<bool>(b->pl().oneWay()))
    return false;
  if (a->pl().lvl() != b->pl().lvl()) return false;
  if (a->pl().getLines().size() != b->pl().getLines().size()) return false;
  if (a->pl().oneWay() && b->pl().oneWay()) {
    if (a->getFrom() != b->getTo() && a->getTo() != b->getFrom()) return false;
  }
  if (a->pl().isRestricted() || b->pl().isRestricted()) return false;
  if (a->pl().getLines() != b->pl().getLines()) return false;

  return true;
}

// _____________________________________________________________________________
const EdgePL& OsmBuilder::mergeEdgePL(Edge* a, Edge* b) {
  const Node* n = 0;
  if (a->getFrom() == b->getFrom())
    n = a->getFrom();
  else if (a->getFrom() == b->getTo())
    n = a->getFrom();
  else
    n = a->getTo();

  if (a->pl().getGeom() == 0) {
    a->pl().addPoint(*a->getFrom()->pl().getGeom());
    a->pl().addPoint(*a->getTo()->pl().getGeom());
  }

  if (a->getTo() == n && b->getTo() == n) {
    // --> n <--
    if (b->pl().getGeom()) {
      a->pl().getGeom()->insert(a->pl().getGeom()->end(),
                                b->pl().getGeom()->rbegin(),
                                b->pl().getGeom()->rend());
    } else {
      a->pl().getGeom()->push_back(*b->getFrom()->pl().getGeom());
    }
  } else if (a->getTo() == n && b->getFrom() == n) {
    // --> n -->
    if (b->pl().getGeom()) {
      a->pl().getGeom()->insert(a->pl().getGeom()->end(),
                                b->pl().getGeom()->begin(),
                                b->pl().getGeom()->end());
    } else {
      a->pl().getGeom()->push_back(*b->getTo()->pl().getGeom());
    }
  } else if (a->getFrom() == n && b->getTo() == n) {
    // <-- n <--
    std::reverse(a->pl().getGeom()->begin(), a->pl().getGeom()->end());
    if (b->pl().getGeom()) {
      a->pl().getGeom()->insert(a->pl().getGeom()->end(),
                                b->pl().getGeom()->rbegin(),
                                b->pl().getGeom()->rend());
    } else {
      a->pl().getGeom()->push_back(*b->getFrom()->pl().getGeom());
    }
  } else {
    // <-- n -->
    std::reverse(a->pl().getGeom()->begin(), a->pl().getGeom()->end());
    if (b->pl().getGeom()) {
      a->pl().getGeom()->insert(a->pl().getGeom()->end(),
                                b->pl().getGeom()->begin(),
                                b->pl().getGeom()->end());
    } else {
      a->pl().getGeom()->push_back(*b->getTo()->pl().getGeom());
    }
  }

  return a->pl();
}

// _____________________________________________________________________________
void OsmBuilder::collapseEdges(Graph* g) {
  for (auto n : g->getNds()) {
    if (n->getOutDeg() + n->getInDeg() != 2 || n->pl().getSI() ||
        n->pl().isTurnCycle())
      continue;

    Edge* ea;
    Edge* eb;
    if (n->getOutDeg() == 2) {
      ea = *n->getAdjListOut().begin();
      eb = *n->getAdjListOut().rbegin();
    } else if (n->getInDeg() == 2) {
      ea = *n->getAdjListIn().begin();
      eb = *n->getAdjListIn().rbegin();
    } else {
      ea = *n->getAdjListOut().begin();
      eb = *n->getAdjListIn().begin();
    }

    // important, we don't have a multigraph! if the same edge
    // already exists, leave this node
    if (g->getEdg(ea->getOtherNd(n), eb->getOtherNd(n))) continue;
    if (g->getEdg(eb->getOtherNd(n), ea->getOtherNd(n))) continue;

    if (edgesSim(ea, eb)) {
      if (ea->pl().oneWay() && ea->getOtherNd(n) != ea->getFrom()) {
        g->addEdg(eb->getOtherNd(n), ea->getOtherNd(n), mergeEdgePL(eb, ea));
      } else {
        g->addEdg(ea->getOtherNd(n), eb->getOtherNd(n), mergeEdgePL(ea, eb));
      }

      g->delEdg(ea->getFrom(), ea->getTo());
      g->delEdg(eb->getFrom(), eb->getTo());
    }
  }
}

// _____________________________________________________________________________
void OsmBuilder::simplifyGeoms(Graph* g) {
  for (auto* n : g->getNds()) {
    for (auto* e : n->getAdjListOut()) {
      (*e->pl().getGeom()) =
          util::geo::simplify(*e->pl().getGeom(), 0.5 / M_PER_DEG);
    }
  }
}

// _____________________________________________________________________________
uint32_t OsmBuilder::writeComps(Graph* g, const OsmReadOpts& opts) {
  NodePL::comps.clear();
  NodePL::comps.emplace_back(Component{0});
  uint32_t numC = 0;
  uint64_t numNds = 0;

  double fac = opts.maxSpeedCorFac;

  for (auto* n : g->getNds()) {
    if (!n->pl().getCompId()) {
      std::stack<std::pair<Node*, Edge*>> q;
      q.push(std::pair<Node*, Edge*>(n, 0));
      while (!q.empty()) {
        std::pair<Node*, Edge*> cur = q.top();
        q.pop();

        cur.first->pl().setComp(NodePL::comps.size());
        numNds++;
        for (auto* e : cur.first->getAdjListOut()) {
          double speed = opts.levelDefSpeed[e->pl().lvl()] / fac;
          if (speed > NodePL::comps.back().maxSpeed)
            NodePL::comps.back().maxSpeed = speed;
          if (!e->getOtherNd(cur.first)->pl().getCompId())
            q.push(std::pair<Node*, Edge*>(e->getOtherNd(cur.first), e));
        }
        for (auto* e : cur.first->getAdjListIn()) {
          double speed = opts.levelDefSpeed[e->pl().lvl()] / fac;
          if (speed > NodePL::comps.back().maxSpeed)
            NodePL::comps.back().maxSpeed = speed;
          if (!e->getOtherNd(cur.first)->pl().getCompId())
            q.push(std::pair<Node*, Edge*>(e->getOtherNd(cur.first), e));
        }
      }

      if (numNds > 1) numC++;
      NodePL::comps.emplace_back(Component{0});
      numNds = 0;
    }
  }

  return numC;
}

// _____________________________________________________________________________
void OsmBuilder::writeEdgeTracks(const EdgTracks& tracks) {
  for (const auto& tr : tracks) {
    if (tr.first->getTo()->pl().getSI() &&
        tr.first->getTo()->pl().getSI()->getTrack().empty()) {
      tr.first->getTo()->pl().getSI()->setTrack(tr.second);
    }
    if (tr.first->getFrom()->pl().getSI() &&
        tr.first->getFrom()->pl().getSI()->getTrack().empty()) {
      tr.first->getFrom()->pl().getSI()->setTrack(tr.second);
    }
  }
}

// _____________________________________________________________________________
void OsmBuilder::writeODirEdgs(Graph* g, Restrictor* restor) {
  for (auto* n : g->getNds()) {
    for (auto* e : n->getAdjListOut()) {
      if (g->getEdg(e->getTo(), e->getFrom())) continue;
      auto newE = g->addEdg(e->getTo(), e->getFrom(), e->pl().revCopy());
      assert(newE->pl().getGeom());
      if (e->pl().isRestricted()) restor->duplicateEdge(e, newE);
    }
  }
}

// _____________________________________________________________________________
void OsmBuilder::writeSelfEdgs(Graph* g) {
  // if a station only has degree 1, there is no way to arrive at this station
  // without doing a full turn (because the outgoing candidate edge is always
  // the incoming edge). This is a problem at end-stations. We solve this by
  // adding self-edges with infinite costs - this still allows usage as
  // arrivals, does not punish bends (because the node degree is still only 2)
  // and prevents the usage of the edge to circumvent turn penalties
  for (auto* n : g->getNds()) {
    if (n->pl().getSI() && n->getAdjListOut().size() == 1) {
      auto e = g->addEdg(n, n);
      e->pl().setCost(std::numeric_limits<uint32_t>::max());
      e->pl().addPoint(*e->getFrom()->pl().getGeom());
      e->pl().addPoint(*e->getTo()->pl().getGeom());
    }
  }
}

// _____________________________________________________________________________
void OsmBuilder::writeNoLinePens(Graph* g, const OsmReadOpts& opts) {
  for (auto* n : g->getNds()) {
    for (auto* e : n->getAdjListOut()) {
      if (e->pl().getLines().size() == 0) {
        double c = e->pl().getCost();
        c = c / 10.0;  // convert into seconds
        e->pl().setCost(costToInt(c * opts.noLinesPunishFact));
      }
    }
  }
}

// _____________________________________________________________________________
void OsmBuilder::writeOneWayPens(Graph* g, const OsmReadOpts& opts) {
  for (auto* n : g->getNds()) {
    for (auto* e : n->getAdjListOut()) {
      if (e->pl().oneWay() == 2) {
        double c = e->pl().getCost();
        c = c / 10.0;  // convert into seconds
        e->pl().setCost(
            costToInt(c * opts.oneWaySpeedPen + opts.oneWayEntryCost));
      }
    }
  }
}

// _____________________________________________________________________________
bool OsmBuilder::keepFullTurn(const trgraph::Node* n, double ang) {
  if (n->getInDeg() + n->getOutDeg() != 1) return false;

  const trgraph::Edge* e = 0;
  if (n->getOutDeg())
    e = n->getAdjListOut().front();
  else
    e = n->getAdjListIn().front();

  auto other = e->getOtherNd(n);

  if (other->getInDeg() + other->getOutDeg() == 3) {
    const trgraph::Edge* a = 0;
    const trgraph::Edge* b = 0;
    for (auto f : other->getAdjListIn()) {
      if (f != e && !a)
        a = f;
      else if (f != e && !b)
        b = f;
    }

    for (auto f : other->getAdjListOut()) {
      if (f != e && !a)
        a = f;
      else if (f != e && !b)
        b = f;
    }

    POINT ap, bp;

    if (!a || !b) return false;

    if (a->pl().getGeom() && b->pl().getGeom()) {
      ap = a->pl().backHop();
      bp = b->pl().backHop();
      if (a->getTo() != other) ap = a->pl().frontHop();
      if (b->getTo() != other) bp = b->pl().frontHop();
    } else {
      assert(!a->pl().getGeom());
      assert(!b->pl().getGeom());
      ap = *a->getTo()->pl().getGeom();
      bp = *b->getTo()->pl().getGeom();
      if (a->getTo() != other) ap = *a->getFrom()->pl().getGeom();
      if (b->getTo() != other) bp = *b->getFrom()->pl().getGeom();
    }

    return util::geo::innerProd(*other->pl().getGeom(), ap, bp) > ang;
  }

  return false;
}

// _____________________________________________________________________________
void OsmBuilder::snapStats(const OsmReadOpts& opts, Graph* g,
                           const BBoxIdx& bbox, double gridSize,
                           Restrictor* res, const NodeSet& orphanStations) {
  NodeGrid sng = buildNodeIdx(g, gridSize, bbox.getFullBox(), true);
  EdgeGrid eg = buildEdgeIdx(g, gridSize, bbox.getFullBox());

  LOG(DEBUG) << "Grid size of " << sng.getXWidth() << "x" << sng.getYHeight();

  for (double d : opts.maxOsmStationDistances) {
    for (auto s : orphanStations) {
      NodePL pl = s->pl();
      snapStation(g, &pl, &eg, &sng, opts, res, d);
    }
  }
}

// _____________________________________________________________________________
uint32_t OsmBuilder::costToInt(double c) {
  // always round upwards, otherwise when combined with the heuristic which
  // is always rounded downwards the PQ monotonicity is not ensured anymore -
  // with a downward rounding, the rounding errors may sum up so high that the
  // path will get cheaper than the heuristic cost
  uint32_t val = std::ceil(c * 10);
  if (std::ceil(c * 10) > std::numeric_limits<uint32_t>::max()) {
    LOG(DEBUG) << "Cost " << c
               << " does not fit in unsigned 32 bit integer, defaulting to "
               << std::numeric_limits<uint32_t>::max() << ".";
    return std::numeric_limits<uint32_t>::max();
  }
  return val;
}

// _____________________________________________________________________________
uint32_t OsmBuilder::parseHexColor(std::string s) const {
  // TODO(patrick): not very nice
  size_t proced = 0;
  std::transform(s.begin(), s.end(), s.begin(), ::toupper);
  std::string ret = "      ";
  if (s.size() == 7 && s[0] == '#') {
    for (size_t i = 1; i < 7; i++) {
      if (isdigit(s[i]))
        ret[i - 1] = s[i];
      else if (isalpha(s[i]) && (s[i] > 64 && s[i] < 71))
        ret[i - 1] = s[i];
      else
        return ad::cppgtfs::gtfs::NO_COLOR;
    }

    return std::stoul("0x" + ret, &proced, 16);
  }

  if (s.size() == 4 && s[0] == '#') {
    for (size_t i = 1; i < 4; i++) {
      if (isdigit(s[i])) {
        ret[(i - 1) * 2] = s[i];
        ret[(i - 1) * 2 + 1] = s[i];
      } else if (isalpha(s[i]) && (s[i] > 64 && s[i] < 71)) {
        ret[(i - 1) * 2] = s[i];
        ret[(i - 1) * 2 + 1] = s[i];
      } else {
        return ad::cppgtfs::gtfs::NO_COLOR;
      }
    }
    return std::stoul("0x" + ret, &proced, 16);
  }

  if (s == "BLACK") return 0x00000000;
  if (s == "SILVER") return 0x00C0C0C0;
  if (s == "GRAY") return 0x00808080;
  if (s == "WHITE") return 0x00FFFFFF;
  if (s == "MAROON") return 0x00800000;
  if (s == "RED") return 0x00FF0000;
  if (s == "PURPLE") return 0x00800080;
  if (s == "FUCHSIA") return 0x00FF00FF;
  if (s == "GREEN") return 0x00008000;
  if (s == "LIME") return 0x0000FF00;
  if (s == "OLIVE") return 0x00808000;
  if (s == "YELLOW") return 0x00FFFF00;
  if (s == "NAVY") return 0x00000080;
  if (s == "BLUE") return 0x000000FF;
  if (s == "TEAL") return 0x00008080;
  if (s == "AQUA") return 0x0000FFFF;

  if (ret.empty()) return ad::cppgtfs::gtfs::NO_COLOR;
  return std::stoul("0x" + ret, &proced, 16);
}

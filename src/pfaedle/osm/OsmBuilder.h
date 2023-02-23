// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_OSM_OSMBUILDER_H_
#define PFAEDLE_OSM_OSMBUILDER_H_
#include <map>
#include <queue>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "pfaedle/Def.h"
#include "pfaedle/osm/BBoxIdx.h"
#include "pfaedle/osm/OsmFilter.h"
#include "pfaedle/osm/OsmIdSet.h"
#include "pfaedle/osm/OsmReadOpts.h"
#include "pfaedle/osm/Restrictor.h"
#include "pfaedle/router/Router.h"
#include "pfaedle/trgraph/Graph.h"
#include "pfaedle/trgraph/Normalizer.h"
#include "pfaedle/trgraph/StatInfo.h"
#include "util/Nullable.h"
#include "util/geo/Geo.h"
#include "util/xml/XmlWriter.h"
#include "pfxml/pfxml.h"

namespace pfaedle {
namespace osm {

using ad::cppgtfs::gtfs::Stop;
using pfaedle::router::NodeSet;
using pfaedle::trgraph::Component;
using pfaedle::trgraph::Edge;
using pfaedle::trgraph::EdgeGrid;
using pfaedle::trgraph::EdgePL;
using pfaedle::trgraph::Graph;
using pfaedle::trgraph::Node;
using pfaedle::trgraph::NodeGrid;
using pfaedle::trgraph::NodePL;
using pfaedle::trgraph::Normalizer;
using pfaedle::trgraph::StatInfo;
using pfaedle::trgraph::TransitEdgeLine;
using util::Nullable;

struct NodeCand {
  double dist;
  Node* node;
  const Edge* fromEdge;
  int fullTurns;
};

struct SearchFunc {
  virtual bool operator()(const Node* n, const StatInfo* si) const = 0;
};

struct EqSearch : public SearchFunc {
  EqSearch() {}
  double minSimi = 0.9;
  bool operator()(const Node* cand, const StatInfo* si) const;
};

struct BlockSearch : public SearchFunc {
  bool operator()(const Node* n, const StatInfo* si) const {
    if (n->pl().getSI() && n->pl().getSI()->simi(si) < 0.5) return true;
    return n->pl().isBlocker();
  }
};

inline bool operator<(const NodeCand& a, const NodeCand& b) {
  return a.fullTurns > b.fullTurns || a.dist > b.dist;
}

typedef std::priority_queue<NodeCand> NodeCandPQ;

/*
 * Builds a physical transit network graph from OSM data
 */
class OsmBuilder {
 public:
  OsmBuilder();

  // Read the OSM file at path, and write a graph to g. Only elements
  // inside the bounding box will be read
  void read(const std::string& path, const OsmReadOpts& opts, Graph* g,
            const BBoxIdx& box, double gridSize, Restrictor* res);

  // Based on the list of options, output an overpass XML query for getting
  // the data needed for routing
  void overpassQryWrite(std::ostream* out, const std::vector<OsmReadOpts>& opts,
                        const BBoxIdx& latLngBox) const;

  // Based on the list of options, output an osmfilter configuration file
  // to filter the data needed for routing
  void osmfilterRuleWrite(std::ostream* out,
                          const std::vector<OsmReadOpts>& opts,
                          const BBoxIdx& latLngBox) const;

  // Based on the list of options, read an OSM file from in and output an
  // OSM file to out which contains exactly the entities that are needed
  // from the file at in
  void filterWrite(const std::string& in, const std::string& out,
                   const std::vector<OsmReadOpts>& opts, const BBoxIdx& box);

 private:
  pfxml::parser_state readBBoxNds(pfxml::file* xml, OsmIdSet* nodes,
                                  OsmIdSet* noHupNodes, const OsmFilter& filter,
                                  const BBoxIdx& bbox) const;

  void readRels(pfxml::file* f, RelLst* rels, RelMap* nodeRels, RelMap* wayRels,
                const OsmFilter& filter, const AttrKeySet& keepAttrs,
                Restrictions* rests) const;

  void readRestr(const OsmRel& rel, Restrictions* rests,
                 const OsmFilter& filter) const;

  void readNodes(pfxml::file* f, Graph* g, const RelLst& rels,
                 const RelMap& nodeRels, const OsmFilter& filter,
                 const OsmIdSet& bBoxNodes, NIdMap* nodes,
                 NIdMultMap* multNodes, NodeSet* orphanStations,
                 const AttrKeySet& keepAttrs, const FlatRels& flatRels,
                 const OsmReadOpts& opts) const;

  void readWriteNds(pfxml::file* i, util::xml::XmlWriter* o,
                    const RelMap& nodeRels, const OsmFilter& filter,
                    const OsmIdSet& bBoxNodes, NIdMap* nodes,
                    const AttrKeySet& keepAttrs, const FlatRels& f) const;

  void readWriteWays(pfxml::file* i, util::xml::XmlWriter* o, OsmIdList* ways,
                     const AttrKeySet& keepAttrs) const;

  void readWriteRels(pfxml::file* i, util::xml::XmlWriter* o, OsmIdList* ways,
                     NIdMap* nodes, const OsmFilter& filter,
                     const AttrKeySet& keepAttrs);

  void readEdges(pfxml::file* xml, Graph* g, const RelLst& rels,
                 const RelMap& wayRels, const OsmFilter& filter,
                 const OsmIdSet& bBoxNodes, NIdMap* nodes,
                 NIdMultMap* multNodes, const OsmIdSet& noHupNodes,
                 const AttrKeySet& keepAttrs, const Restrictions& rest,
                 Restrictor* restor, const FlatRels& flatRels,
                 EdgTracks* etracks, const OsmReadOpts& opts);

  void readEdges(pfxml::file* xml, const RelMap& wayRels,
                 const OsmFilter& filter, const OsmIdSet& bBoxNodes,
                 const AttrKeySet& keepAttrs, OsmIdList* ret, NIdMap* nodes,
                 const FlatRels& flatRels);

  OsmWay nextWay(pfxml::file* xml, const RelMap& wayRels,
                 const OsmFilter& filter, const OsmIdSet& bBoxNodes,
                 const AttrKeySet& keepAttrs, const FlatRels& flatRels) const;

  bool keepWay(const OsmWay& w, const RelMap& wayRels, const OsmFilter& filter,
               const OsmIdSet& bBoxNodes, const FlatRels& fl) const;

  OsmWay nextWayWithId(pfxml::file* xml, osmid wid,
                       const AttrKeySet& keepAttrs) const;

  OsmNode nextNode(pfxml::file* xml, NIdMap* nodes, NIdMultMap* multNodes,
                   const RelMap& nodeRels, const OsmFilter& filter,
                   const OsmIdSet& bBoxNodes, const AttrKeySet& keepAttrs,
                   const FlatRels& flatRels) const;

  bool keepNode(const OsmNode& n, const NIdMap& nodes,
                const NIdMultMap& multNodes, const RelMap& nodeRels,
                const OsmIdSet& bBoxNodes, const OsmFilter& filter,
                const FlatRels& fl) const;

  OsmRel nextRel(pfxml::file* xml, const OsmFilter& filter,
                 const AttrKeySet& keepAttrs) const;

 protected:
  Nullable<StatInfo> getStatInfo(osmid nid, const AttrMap& m,
                                 const RelMap& nodeRels, const RelLst& rels,
                                 const OsmReadOpts& ops) const;

  static void snapStats(const OsmReadOpts& opts, Graph* g, const BBoxIdx& bbox,
                        double gridSize, Restrictor* res,
                        const NodeSet& orphanStations);
  static void writeGeoms(Graph* g, const OsmReadOpts& opts);
  static void deleteOrphEdgs(Graph* g, const OsmReadOpts& opts);
  static void deleteOrphNds(Graph* g, const OsmReadOpts& opts);
  static double dist(const Node* a, const Node* b);

  static NodeGrid buildNodeIdx(Graph* g, double size, const BOX& box,
                               bool which);

  static EdgeGrid buildEdgeIdx(Graph* g, double size, const BOX& box);

  static void fixGaps(Graph* g, NodeGrid* ng);
  static void collapseEdges(Graph* g);
  static void writeODirEdgs(Graph* g, Restrictor* restor);
  static void writeSelfEdgs(Graph* g);
  static void writeOneWayPens(Graph* g, const OsmReadOpts& opts);
  static void writeNoLinePens(Graph* g, const OsmReadOpts& opts);
  static void writeEdgeTracks(const EdgTracks& tracks);
  static void simplifyGeoms(Graph* g);
  static uint32_t writeComps(Graph* g, const OsmReadOpts& opts);
  static bool edgesSim(const Edge* a, const Edge* b);
  static const EdgePL& mergeEdgePL(Edge* a, Edge* b);
  static void getEdgCands(const POINT& s, EdgeCandPQ* ret, EdgeGrid* eg,
                          double d);

  static void snapStation(Graph* g, NodePL* s, EdgeGrid* eg, NodeGrid* sng,
                          const OsmReadOpts& opts, Restrictor* restor,
                          double maxD);

  // Checks if from the edge e, a station similar to si can be reach with less
  // than maxD distance and less or equal to "maxFullTurns" full turns. If
  // such a station exists, it is returned. If not, 0 is returned.
  static Node* eqStatReach(const Edge* e, const StatInfo* si, const POINT& p,
                           double maxD, int maxFullTurns, double maxAng);

  static Node* depthSearch(const Edge* e, const StatInfo* si, const POINT& p,
                           double maxD, int maxFullTurns, double minAngle,
                           const SearchFunc& sfunc);

  static bool isBlocked(const Edge* e, const StatInfo* si, const POINT& p,
                        double maxD, int maxFullTurns, double minAngle);
  static bool keepFullTurn(const trgraph::Node* n, double ang);

  static NodePL plFromGtfs(const Stop* s, const OsmReadOpts& ops);

  std::vector<TransitEdgeLine*> getLines(const std::vector<size_t>& edgeRels,
                                         const RelLst& rels,
                                         const OsmReadOpts& ops);

  void getKeptAttrKeys(const OsmReadOpts& opts, AttrKeySet sets[3]) const;

  void skipUntil(pfxml::file* xml, const std::string& s) const;

  void processRestr(osmid nid, osmid wid, const Restrictions& rawRests, Edge* e,
                    Node* n, Restrictor* restor) const;

  std::string getAttrByFirstMatch(const DeepAttrLst& rule, osmid id,
                                  const AttrMap& attrs, const RelMap& entRels,
                                  const RelLst& rels,
                                  const Normalizer& norm) const;

  std::vector<std::string> getAttrMatchRanked(const DeepAttrLst& rule, osmid id,
                                              const AttrMap& attrs,
                                              const RelMap& entRels,
                                              const RelLst& rels,
                                              const Normalizer& norm) const;

  std::string getAttr(const DeepAttrRule& s, osmid id, const AttrMap& attrs,
                      const RelMap& entRels, const RelLst& rels) const;

  bool relKeep(osmid id, const RelMap& rels, const FlatRels& fl) const;

  uint32_t parseHexColor(std::string) const;

  static uint32_t costToInt(double c);

  std::map<TransitEdgeLine, TransitEdgeLine*> _lines;
  std::map<size_t, TransitEdgeLine*> _relLines;
};
}  // namespace osm
}  // namespace pfaedle
#endif  // PFAEDLE_OSM_OSMBUILDER_H_

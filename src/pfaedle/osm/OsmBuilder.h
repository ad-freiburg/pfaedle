// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_OSM_OSMBUILDER_H_
#define PFAEDLE_OSM_OSMBUILDER_H_

#include <queue>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "ad/cppgtfs/gtfs/Feed.h"
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
#include "util/xml/XmlWriter.h"
#include "xml/File.h"

namespace pfaedle {
namespace osm {

using pfaedle::trgraph::EdgeGrid;
using pfaedle::trgraph::NodeGrid;
using pfaedle::trgraph::Normalizer;
using pfaedle::trgraph::Graph;
using pfaedle::trgraph::Node;
using pfaedle::trgraph::NodePL;
using pfaedle::trgraph::Edge;
using pfaedle::trgraph::EdgePL;
using pfaedle::trgraph::TransitEdgeLine;
using pfaedle::trgraph::StatInfo;
using pfaedle::trgraph::StatGroup;
using pfaedle::trgraph::Component;
using pfaedle::router::NodeSet;
using ad::cppgtfs::gtfs::Stop;
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
  double minSimi = 0.9;
  bool operator()(const Node* cand, const StatInfo* si) const {
    return cand->pl().getSI() && cand->pl().getSI()->simi(si) > minSimi;
  }
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
            const BBoxIdx& box, size_t gridSize, router::FeedStops* fs,
            Restrictor* res);


  // Based on the list of options, read an OSM file from in and output an
  // OSM file to out which contains exactly the entities that are needed
  // from the file at in
  void filterWrite(const std::string& in, const std::string& out,
                   const std::vector<OsmReadOpts>& opts, const BBoxIdx& box);

 private:
  xml::ParserState readBBoxNds(xml::File* xml, OsmIdSet* nodes,
                               OsmIdSet* noHupNodes, const OsmFilter& filter,
                               const BBoxIdx& bbox) const;

  void readRels(xml::File* f, RelLst* rels, RelMap* nodeRels, RelMap* wayRels,
                const OsmFilter& filter, const AttrKeySet& keepAttrs,
                Restrictions* rests) const;

  void readRestr(const OsmRel& rel, Restrictions* rests,
                 const OsmFilter& filter) const;

  void readNodes(xml::File* f, Graph* g, const RelLst& rels,
                 const RelMap& nodeRels, const OsmFilter& filter,
                 const OsmIdSet& bBoxNodes, NIdMap* nodes,
                 NIdMultMap* multNodes, NodeSet* orphanStations,
                 const AttrKeySet& keepAttrs, const FlatRels& flatRels,
                 const OsmReadOpts& opts) const;

  void readWriteNds(xml::File* i, util::xml::XmlWriter* o,
                    const RelMap& nodeRels, const OsmFilter& filter,
                    const OsmIdSet& bBoxNodes, NIdMap* nodes,
                    const AttrKeySet& keepAttrs, const FlatRels& f) const;

  void readWriteWays(xml::File* i, util::xml::XmlWriter* o, OsmIdList* ways,
                     const AttrKeySet& keepAttrs) const;

  void readWriteRels(xml::File* i, util::xml::XmlWriter* o, OsmIdList* ways,
                     NIdMap* nodes, const OsmFilter& filter,
                     const AttrKeySet& keepAttrs);

  void readEdges(xml::File* xml, Graph* g, const RelLst& rels,
                 const RelMap& wayRels, const OsmFilter& filter,
                 const OsmIdSet& bBoxNodes, NIdMap* nodes,
                 NIdMultMap* multNodes, const OsmIdSet& noHupNodes,
                 const AttrKeySet& keepAttrs, const Restrictions& rest,
                 Restrictor* restor, const FlatRels& flatRels,
                 EdgTracks* etracks, const OsmReadOpts& opts);

  void readEdges(xml::File* xml, const RelMap& wayRels, const OsmFilter& filter,
                 const OsmIdSet& bBoxNodes, const AttrKeySet& keepAttrs,
                 OsmIdList* ret, NIdMap* nodes, const FlatRels& flatRels);

  OsmWay nextWay(xml::File* xml, const RelMap& wayRels, const OsmFilter& filter,
                 const OsmIdSet& bBoxNodes, const AttrKeySet& keepAttrs,
                 const FlatRels& flatRels) const;

  bool keepWay(const OsmWay& w, const RelMap& wayRels, const OsmFilter& filter,
               const OsmIdSet& bBoxNodes, const FlatRels& fl) const;

  OsmWay nextWayWithId(xml::File* xml, osmid wid,
                       const AttrKeySet& keepAttrs) const;

  OsmNode nextNode(xml::File* xml, NIdMap* nodes, NIdMultMap* multNodes,
                   const RelMap& nodeRels, const OsmFilter& filter,
                   const OsmIdSet& bBoxNodes, const AttrKeySet& keepAttrs,
                   const FlatRels& flatRels) const;

  bool keepNode(const OsmNode& n, const NIdMap& nodes,
                const NIdMultMap& multNodes, const RelMap& nodeRels,
                const OsmIdSet& bBoxNodes, const OsmFilter& filter,
                const FlatRels& fl) const;

  OsmRel nextRel(xml::File* xml, const OsmFilter& filter,
                 const AttrKeySet& keepAttrs) const;

  Nullable<StatInfo> getStatInfo(Node* node, osmid nid, const FPoint& pos,
                                 const AttrMap& m, StAttrGroups* groups,
                                 const RelMap& nodeRels, const RelLst& rels,
                                 const OsmReadOpts& ops) const;

  void writeGeoms(Graph* g) const;
  void deleteOrphNds(Graph* g) const;
  void deleteOrphEdgs(Graph* g) const;
  double dist(const Node* a, const Node* b) const;
  double webMercDist(const Node* a, const Node* b) const;
  double webMercDistFactor(const FPoint& a) const;

  NodeGrid buildNodeIdx(Graph* g, size_t size,
                        const util::geo::Box<float>& webMercBox,
                        bool which) const;

  EdgeGrid buildEdgeIdx(Graph* g, size_t size,
                        const util::geo::Box<float>& webMercBox) const;

  void fixGaps(Graph* g, NodeGrid* ng) const;
  void collapseEdges(Graph* g) const;
  void writeODirEdgs(Graph* g, Restrictor* restor) const;
  void writeSelfEdgs(Graph* g) const;
  void writeEdgeTracks(const EdgTracks& tracks) const;
  void simplifyGeoms(Graph* g) const;
  uint32_t writeComps(Graph* g) const;
  bool edgesSim(const Edge* a, const Edge* b) const;
  const EdgePL& mergeEdgePL(Edge* a, Edge* b) const;
  void getEdgCands(const FPoint& s, EdgeCandPQ* ret, EdgeGrid* eg,
                   double d) const;

  std::set<Node*> getMatchingNds(const NodePL& s, NodeGrid* ng, double d) const;

  Node* getMatchingNd(const NodePL& s, NodeGrid* ng, double d) const;

  NodeSet snapStation(Graph* g, NodePL* s, EdgeGrid* eg, NodeGrid* sng,
                      const OsmReadOpts& opts, Restrictor* restor, bool surHeur,
                      double maxD) const;

  // Checks if from the edge e, a station similar to si can be reach with less
  // than maxD distance and less or equal to "maxFullTurns" full turns. If
  // such a station exists, it is returned. If not, 0 is returned.
  Node* eqStatReach(const Edge* e, const StatInfo* si, const FPoint& p,
                    double maxD, int maxFullTurns, double maxAng) const;

  Node* depthSearch(const Edge* e, const StatInfo* si,
                    const util::geo::FPoint& p, double maxD, int maxFullTurns,
                    double minAngle, const SearchFunc& sfunc) const;

  bool isBlocked(const Edge* e, const StatInfo* si, const FPoint& p,
                 double maxD, int maxFullTurns, double minAngle) const;

  StatGroup* groupStats(const NodeSet& s) const;

  std::vector<TransitEdgeLine*> getLines(const std::vector<size_t>& edgeRels,
                                         const RelLst& rels,
                                         const OsmReadOpts& ops);

  NodePL plFromGtfs(const Stop* s, const OsmReadOpts& ops) const;

  void getKeptAttrKeys(const OsmReadOpts& opts, AttrKeySet sets[3]) const;

  void skipUntil(xml::File* xml, const std::string& s) const;

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

  std::map<TransitEdgeLine, TransitEdgeLine*> _lines;
  std::map<size_t, TransitEdgeLine*> _relLines;
};
}  // namespace osm
}  // namespace pfaedle
#endif  // PFAEDLE_OSM_OSMBUILDER_H_

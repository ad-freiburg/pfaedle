// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_OSM_RESTRICTOR_H_
#define PFAEDLE_OSM_RESTRICTOR_H_

#include <unordered_map>
#include <map>
#include <vector>
#include <utility>
#include "pfaedle/osm/Osm.h"
#include "pfaedle/trgraph/Graph.h"

namespace pfaedle {
namespace osm {

typedef std::pair<const trgraph::Edge*, const trgraph::Edge*> RulePair;
typedef std::vector<RulePair> RuleVec;
typedef std::pair<const trgraph::Node*, size_t> DanglPath;
// very seldom, there are more than a handful of rules for a node. Use a
// vector here, should have lesser overhead and be faster for such small
// numbers
typedef std::unordered_map<const trgraph::Node*, RuleVec> Rules;
typedef std::pair<const trgraph::Node*, osmid> NodeOsmIdP;

/*
 * Stores restrictions between edges
 */
class Restrictor {
 public:
  Restrictor() {}

  void relax(osmid wid, const trgraph::Node* n, const trgraph::Edge* e);
  void add(const trgraph::Edge* from, osmid to, const trgraph::Node* via,
           bool pos);
  bool may(const trgraph::Edge* from, const trgraph::Edge* to,
           const trgraph::Node* via) const;
  void replaceEdge(const trgraph::Edge* old, const trgraph::Edge* newA,
                   const trgraph::Edge* newB);
  void duplicateEdge(const trgraph::Edge* old, const trgraph::Node* via,
                     const trgraph::Edge* newE);
  void duplicateEdge(const trgraph::Edge* old, const trgraph::Edge* newE);

 private:
  Rules _pos;
  Rules _neg;

  std::map<NodeOsmIdP, const trgraph::Edge*> _rlx;

  std::map<NodeOsmIdP, std::vector<DanglPath>> _posDangling;
  std::map<NodeOsmIdP, std::vector<DanglPath>> _negDangling;

  void replaceEdge(const trgraph::Edge* old, const trgraph::Node* via,
                   const trgraph::Edge* newE);
};
}  // namespace osm
}  // namespace pfaedle

#endif  // PFAEDLE_OSM_RESTRICTOR_H_

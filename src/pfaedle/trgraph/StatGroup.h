// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_TRGRAPH_STATGROUP_H_
#define PFAEDLE_TRGRAPH_STATGROUP_H_

#include <string>
#include <unordered_map>
#include <set>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "pfaedle/router/Router.h"
#include "pfaedle/trgraph/Graph.h"
#include "pfaedle/trgraph/Normalizer.h"

namespace pfaedle {
namespace trgraph {

using ad::cppgtfs::gtfs::Stop;

/*
 * A group of stations that belong together semantically (for example, multiple
 * stop points of a larger bus station)
 */
class StatGroup {
 public:
  StatGroup();
  StatGroup(const StatGroup& a) = delete;

  // Add a stop s to this station group
  void addStop(const Stop* s);

  // Add a node n to this station group
  void addNode(trgraph::Node* n);

  // Return all nodes contained in this group
  const std::set<trgraph::Node*>& getNodes() const;
  std::set<trgraph::Node*>& getNodes();

  // Return all stops contained in this group
  const std::set<const Stop*>& getStops() const;

  // Remove a node from this group
  void remNode(trgraph::Node* n);

  // All nodes in other will be in this group, their SI's updated, and the
  // "other" group deleted.
  void merge(StatGroup* other);

  // Return node candidates for stop s from this group
  const router::NodeCandGroup& getNodeCands(const Stop* s) const;

  // Write the penalties for all stops contained in this group so far.
  void writePens(const trgraph::Normalizer& platformNorm, double trackPen,
                 double distPenFac, double nonOsmPen);

 private:
  std::set<trgraph::Node*> _nodes;
  std::set<const Stop*> _stops;

  // for each stop in this group, a penalty for each of the nodes here, based on
  // its distance and optionally the track number
  std::unordered_map<const Stop*, router::NodeCandGroup> _stopNodePens;

  double getPen(const Stop* s, trgraph::Node* n,
                const trgraph::Normalizer& norm, double trackPen,
                double distPenFac, double nonOsmPen) const;
};
}  // namespace trgraph
}  // namespace pfaedle

#endif  // PFAEDLE_TRGRAPH_STATGROUP_H_

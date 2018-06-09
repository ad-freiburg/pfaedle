// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_CONFIG_MOTCONFIG_H_
#define PFAEDLE_CONFIG_MOTCONFIG_H_

#include "pfaedle/osm/OsmBuilder.h"
#include "pfaedle/router/Router.h"

namespace pfaedle {
namespace config {


struct MotConfig {
  router::MOTs mots;
  osm::OsmReadOpts osmBuildOpts;
  router::RoutingOpts routingOpts;
};

inline bool operator==(const MotConfig& a, const MotConfig& b) {
  return a.osmBuildOpts == b.osmBuildOpts && a.routingOpts == b.routingOpts;
}

}  // namespace config
}  // namespace pfaedle

#endif  // PFAEDLE_CONFIG_MOTCONFIG_H_

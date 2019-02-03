// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_CONFIG_MOTCONFIG_H_
#define PFAEDLE_CONFIG_MOTCONFIG_H_

#include <map>
#include <string>
#include "pfaedle/osm/OsmBuilder.h"
#include "pfaedle/router/Router.h"

namespace pfaedle {
namespace config {

struct MotConfig {
  router::MOTs mots;
  osm::OsmReadOpts osmBuildOpts;
  router::RoutingOpts routingOpts;
  std::map<std::string, std::string> unproced;
};

inline bool operator==(const MotConfig& a, const MotConfig& b) {
  bool unprocedEq = a.unproced.size() == b.unproced.size();
  for (const auto& kv : a.unproced) {
    if (!b.unproced.count(kv.first) ||
        b.unproced.find(kv.first)->second != kv.second) {
      unprocedEq = false;
      break;
    }
  }
  return a.osmBuildOpts == b.osmBuildOpts && a.routingOpts == b.routingOpts &&
         unprocedEq;
}

}  // namespace config
}  // namespace pfaedle

#endif  // PFAEDLE_CONFIG_MOTCONFIG_H_

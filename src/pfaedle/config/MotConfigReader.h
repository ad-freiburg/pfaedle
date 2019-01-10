// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_CONFIG_MOTCONFIGREADER_H_
#define PFAEDLE_CONFIG_MOTCONFIGREADER_H_

#include "pfaedle/_config.h"

#ifndef HOME_VAR
#define HOME_VAR "HOME"
#endif
#ifndef XDG_DATA_HOME_SUFFIX
#define XDG_DATA_HOME_SUFFIX "/.local/share"
#endif
#ifndef XDG_CONFIG_HOME_SUFFIX
#define XDG_CONFIG_HOME_SUFFIX "/.config"
#endif
#ifndef XDG_CACHE_HOME_SUFFIX
#define XDG_CACHE_HOME_SUFFIX "/.cache"
#endif
#ifndef XDG_DATA_DIRS_DEFAULT
#define XDG_DATA_DIRS_DEFAULT "/usr/local/share"
#endif
#ifndef XDG_CONFIG_DIRS_DEFAULT
#define XDG_CONFIG_DIRS_DEFAULT "/etc"
#endif
#ifndef CFG_FILE_NAME
#define CFG_FILE_NAME "pfaedle.cfg"
#endif

#include <unordered_map>
#include <vector>
#include <set>
#include <string>
#include "ad/cppgtfs/gtfs/Route.h"
#include "configparser/ConfigFileParser.h"
#include "pfaedle/config/MotConfig.h"
#include "pfaedle/osm/OsmBuilder.h"

namespace pfaedle {
namespace config {

using ad::cppgtfs::gtfs::Route;

class MotConfigReader {
 public:
  MotConfigReader();
  void parse(const std::vector<std::string>& paths);

  const std::vector<MotConfig>& getConfigs() const;

 private:
  std::vector<MotConfig> _cfgs;

  osm::KeyVal getKv(const std::string& kv) const;
  osm::FilterRule getFRule(const std::string& kv) const;

  trgraph::ReplRules getNormRules(const std::vector<std::string>& arr) const;
  osm::DeepAttrRule getDeepAttrRule(const std::string& rule) const;
  uint64_t getFlags(const std::set<string>& flags) const;
};
}  // namespace config
}  // namespace pfaedle

#endif  // PFAEDLE_CONFIG_MOTCONFIGREADER_H_

// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_CONFIG_PFAEDLECONFIG_H_
#define PFAEDLE_CONFIG_PFAEDLECONFIG_H_

#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "ad/cppgtfs/gtfs/Route.h"
#include "util/geo/Geo.h"

namespace pfaedle {
namespace config {

using ad::cppgtfs::gtfs::Route;

struct Config {
  Config()
      : dbgOutputPath("."),
        solveMethod("global"),
        outputPath("gtfs-out"),
        dropShapes(false),
        useHMM(false),
        writeGraph(false),
        buildTransitGraph(false),
        useCaching(false),
        writeOverpass(false),
        writeOsmfilter(false),
        inPlace(false),
        writeColors(false),
        noFastHops(false),
        noAStar(false),
        noTrie(false),
        noHopCache(false),
        writeStats(false),
        parseAdditionalGTFSFields(false),
        gridSize(2000 / util::geo::M_PER_DEG),
        gaussianNoise(0) {}
  std::string dbgOutputPath;
  std::string solveMethod;
  std::string shapeTripId;
  std::string outputPath;
  std::string writeOsm;
  std::string osmPath;
  std::string motCfgParam;
  std::vector<std::string> feedPaths;
  std::vector<std::string> configPaths;
  std::set<Route::TYPE> mots;
  bool dropShapes;
  bool useHMM;
  bool writeGraph;
  bool buildTransitGraph;
  bool useCaching;
  bool writeOverpass;
  bool writeOsmfilter;
  bool inPlace;
  bool writeColors;
  bool noFastHops;
  bool noAStar;
  bool noTrie;
  bool noHopCache;
  bool writeStats;
  bool parseAdditionalGTFSFields;
  double gridSize;
  double gaussianNoise;

  std::string toString() {
    std::stringstream ss;
    ss << "trip-id: " << shapeTripId << "\n"
       << "output-path: " << outputPath << "\n"
       << "write-osm-path: " << writeOsm << "\n"
       << "read-osm-path: " << osmPath << "\n"
       << "debug-output-path: " << dbgOutputPath << "\n"
       << "drop-shapes: " << dropShapes << "\n"
       << "use-hmm: " << useHMM << "\n"
       << "write-graph: " << writeGraph << "\n"
       << "grid-size: " << gridSize << "\n"
       << "use-cache: " << useCaching << "\n"
       << "write-overpass: " << writeOverpass << "\n"
       << "write-osmfilter: " << writeOsmfilter << "\n"
       << "inplace: " << inPlace << "\n"
       << "write-colors: " << writeColors << "\n"
       << "no-fast-hops: " << noFastHops << "\n"
       << "no-a-star: " << noAStar << "\n"
       << "no-trie: " << noTrie << "\n"
       << "no-hop-cache: " << noHopCache << "\n"
       << "parse-additional-gtfs-fields: " << parseAdditionalGTFSFields << "\n"
       << "write-stats: " << writeStats << "\n"
       << "feed-paths: ";

    for (const auto& p : feedPaths) {
      ss << p << " ";
    }

    ss << "\nconfig-paths: ";

    for (const auto& p : configPaths) {
      ss << p << " ";
    }

    ss << "\nmots: ";

    for (const auto& mot : mots) {
      ss << mot << " ";
    }

    ss << "\n";

    return ss.str();
  }
};

}  // namespace config
}  // namespace pfaedle

#endif  // PFAEDLE_CONFIG_PFAEDLECONFIG_H_

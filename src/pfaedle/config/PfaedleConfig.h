// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_CONFIG_PFAEDLECONFIG_H_
#define PFAEDLE_CONFIG_PFAEDLECONFIG_H_

#include <sstream>
#include <string>
#include <vector>
#include <set>
#include "ad/cppgtfs/gtfs/Route.h"

namespace pfaedle {
namespace config {

using ad::cppgtfs::gtfs::Route;

struct Config {
  Config()
      : dbgOutputPath("geo"),
        solveMethod("global"),
        evalPath("."),
        outputPath("gtfs-out"),
        dropShapes(false),
        useHMM(false),
        writeGraph(false),
        writeCombGraph(false),
        evaluate(false),
        buildTransitGraph(false),
        useCaching(false),
        writeOverpass(false),
        inPlace(false),
        gridSize(2000) {}
  std::string dbgOutputPath;
  std::string solveMethod;
  std::string evalPath;
  std::string shapeTripId;
  std::string outputPath;
  std::string writeOsm;
  std::string osmPath;
  std::string evalDfBins;
  std::vector<std::string> feedPaths;
  std::vector<std::string> configPaths;
  std::set<Route::TYPE> mots;
  bool dropShapes;
  bool useHMM;
  bool writeGraph;
  bool writeCombGraph;
  bool evaluate;
  bool buildTransitGraph;
  bool useCaching;
  bool writeOverpass;
  bool inPlace;
  double gridSize;

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
       << "write-cgraph: " << writeCombGraph << "\n"
       << "grid-size: " << gridSize << "\n"
       << "use-cache: " << useCaching << "\n"
       << "write-overpass: " << writeOverpass << "\n"
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

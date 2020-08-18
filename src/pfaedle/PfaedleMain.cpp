// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>
#include <climits>
#include <pwd.h>
#include <csignal>
#include <cstdio>
#include <sys/stat.h>
#include <unistd.h>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "ad/cppgtfs/Parser.h"
#include "ad/cppgtfs/Writer.h"
#include "pfaedle/config/ConfigReader.h"
#include "pfaedle/config/MotConfig.h"
#include "pfaedle/config/MotConfigReader.h"
#include "pfaedle/eval/Collector.h"
#include "pfaedle/gtfs/Feed.h"
#include "pfaedle/gtfs/Writer.h"
#include "pfaedle/netgraph/Graph.h"
#include "pfaedle/osm/OsmIdSet.h"
#include "pfaedle/router/ShapeBuilder.h"
#include "pfaedle/trgraph/Graph.h"
#include "pfaedle/trgraph/StatGroup.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/geo/output/GeoJsonOutput.h"
#include "util/json/Writer.h"
#include "util/log/Log.h"
#include "util/Misc.h"

#ifndef CFG_HOME_SUFFIX
#define CFG_HOME_SUFFIX "/.config"
#endif
#ifndef CFG_DIR
#define CFG_DIR "/etc"
#endif
#ifndef CFG_FILE_NAME
#define CFG_FILE_NAME "pfaedle.cfg"
#endif

using pfaedle::router::MOTs;
using pfaedle::osm::BBoxIdx;
using pfaedle::osm::OsmBuilder;
using pfaedle::config::MotConfig;
using pfaedle::config::Config;
using pfaedle::router::ShapeBuilder;
using configparser::ParseFileExc;
using pfaedle::config::MotConfigReader;
using pfaedle::config::ConfigReader;
using pfaedle::eval::Collector;

enum class RetCode {
  SUCCESS = 0,
  NO_INPUT_FEED = 1,
  MULT_FEEDS_NOT_ALWD = 2,
  TRIP_NOT_FOUND = 3,
  GTFS_PARSE_ERR = 4,
  NO_OSM_INPUT = 5,
  MOT_CFG_PARSE_ERR = 6,
  OSM_PARSE_ERR = 7,
  GTFS_WRITE_ERR = 8,
  NO_MOT_CFG = 9
};

std::string getFileNameMotStr(const MOTs& mots);
std::vector<std::string> getCfgPaths(const Config& cfg);

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, nullptr);

  // initialize randomness
  srand(time(NULL) + rand());  // NOLINT

  Config cfg;
  MotConfigReader motCfgReader;

  ConfigReader::read(&cfg, argc, argv);

  std::vector<pfaedle::gtfs::Feed> gtfs(cfg.feedPaths.size());
  // feed containing the shapes in memory for evaluation
  ad::cppgtfs::gtfs::Feed evalFeed;

  std::vector<std::string> cfgPaths = getCfgPaths(cfg);

  try {
    motCfgReader.parse(cfgPaths);
  } catch (const configparser::ParseExc& ex) {
    LOG(ERROR) << "Could not parse MOT configurations, reason was:";
    std::cerr << ex.what() << std::endl;
    exit(static_cast<int>(RetCode::MOT_CFG_PARSE_ERR));
  }

  if (cfg.osmPath.empty() && !cfg.writeOverpass) {
    std::cerr << "No OSM input file specified (-x), see --help." << std::endl;
    exit(static_cast<int>(RetCode::NO_OSM_INPUT));
  }

  if (motCfgReader.getConfigs().empty()) {
    LOG(ERROR) << "No MOT configurations specified and no implicit "
                  "configurations found, see --help.";
    exit(static_cast<int>(RetCode::NO_MOT_CFG));
  }

  if (cfg.feedPaths.size() == 1) {
    if (cfg.inPlace) cfg.outputPath = cfg.feedPaths[0];
    if (!cfg.writeOverpass)
      LOG(INFO) << "Reading " << cfg.feedPaths[0] << " ...";
    try {
      ad::cppgtfs::Parser p;
      p.parse(&gtfs[0], cfg.feedPaths[0]);
      if (cfg.evaluate) {
        // read the shapes and store them in memory
        p.parseShapes(&evalFeed, cfg.feedPaths[0]);
      }
    } catch (const ad::cppgtfs::ParserException& ex) {
      LOG(ERROR) << "Could not parse input GTFS feed, reason was:";
      std::cerr << ex.what() << std::endl;
      exit(static_cast<int>(RetCode::GTFS_PARSE_ERR));
    }
    if (!cfg.writeOverpass) LOG(INFO) << "Done.";
  } else if (!cfg.writeOsm.empty() || cfg.writeOverpass) {
    for (size_t i = 0; i < cfg.feedPaths.size(); i++) {
      if (!cfg.writeOverpass)
        LOG(INFO) << "Reading " << cfg.feedPaths[i] << " ...";
      ad::cppgtfs::Parser p;
      try {
        p.parse(&gtfs[i], cfg.feedPaths[i]);
      } catch (const ad::cppgtfs::ParserException& ex) {
        LOG(ERROR) << "Could not parse input GTFS feed, reason was:";
        std::cerr << ex.what() << std::endl;
        exit(static_cast<int>(RetCode::GTFS_PARSE_ERR));
      }
      if (!cfg.writeOverpass) LOG(INFO) << "Done.";
    }
  } else if (cfg.feedPaths.size() > 1) {
    std::cerr << "Multiple feeds only allowed in filter mode." << std::endl;
    exit(static_cast<int>(RetCode::MULT_FEEDS_NOT_ALWD));
  }

  LOG(DEBUG) << "Read " << motCfgReader.getConfigs().size()
             << " unique MOT configs.";
  MOTs cmdCfgMots = cfg.mots;
  pfaedle::gtfs::Trip* singleTrip = nullptr;

  if (!cfg.shapeTripId.empty()) {
    if (cfg.feedPaths.empty()) {
      std::cout << "No input feed specified, see --help" << std::endl;
      exit(static_cast<int>(RetCode::NO_INPUT_FEED));
    }
    singleTrip = gtfs[0].getTrips().get(cfg.shapeTripId);
    if (!singleTrip) {
      LOG(ERROR) << "Trip #" << cfg.shapeTripId << " not found.";
      exit(static_cast<int>(RetCode::TRIP_NOT_FOUND));
    }
  }

  if (!cfg.writeOsm.empty()) {
    LOG(INFO) << "Writing filtered XML to " << cfg.writeOsm << " ...";
    BBoxIdx box(BOX_PADDING);
    for (size_t i = 0; i < cfg.feedPaths.size(); i++) {
      ShapeBuilder::getGtfsBox(&gtfs[i], cmdCfgMots, cfg.shapeTripId, true,
                               &box);
    }
    OsmBuilder osmBuilder;
    std::vector<pfaedle::osm::OsmReadOpts> opts;
    for (const auto& o : motCfgReader.getConfigs()) {
      if (std::find_first_of(o.mots.begin(), o.mots.end(), cmdCfgMots.begin(),
                             cmdCfgMots.end()) != o.mots.end()) {
        opts.push_back(o.osmBuildOpts);
      }
    }
    try {
      osmBuilder.filterWrite(cfg.osmPath, cfg.writeOsm, opts, box);
    } catch (const pfxml::parse_exc& ex) {
      LOG(ERROR) << "Could not parse OSM data, reason was:";
      std::cerr << ex.what() << std::endl;
      exit(static_cast<int>(RetCode::OSM_PARSE_ERR));
    }
    exit(static_cast<int>(RetCode::SUCCESS));
  } else if (cfg.writeOverpass) {
    BBoxIdx box(BOX_PADDING);
    for (size_t i = 0; i < cfg.feedPaths.size(); i++) {
      ShapeBuilder::getGtfsBox(&gtfs[i], cmdCfgMots, cfg.shapeTripId, true,
                               &box);
    }
    OsmBuilder osmBuilder;
    std::vector<pfaedle::osm::OsmReadOpts> opts;
    for (const auto& o : motCfgReader.getConfigs()) {
      if (std::find_first_of(o.mots.begin(), o.mots.end(), cmdCfgMots.begin(),
                             cmdCfgMots.end()) != o.mots.end()) {
        opts.push_back(o.osmBuildOpts);
      }
    }
    osmBuilder.overpassQryWrite(&std::cout, opts, box);
    exit(static_cast<int>(RetCode::SUCCESS));
  } else if (cfg.feedPaths.empty()) {
    std::cout << "No input feed specified, see --help" << std::endl;
    exit(static_cast<int>(RetCode::NO_INPUT_FEED));
  }

  auto dfBinStrings = util::split(std::string(cfg.evalDfBins), ',');
  std::vector<double> dfBins(dfBinStrings.size());
  for (const auto& st : dfBinStrings) {
      dfBins.push_back(atof(st.c_str()));
  }
  Collector ecoll(cfg.evalPath, dfBins);

  for (const auto& motCfg : motCfgReader.getConfigs()) {
    std::string filePost;
    auto usedMots = pfaedle::router::motISect(motCfg.mots, cmdCfgMots);
    if (usedMots.empty())
        continue;
    if (singleTrip && !usedMots.count(singleTrip->getRoute()->getType()))
      continue;
    if (motCfgReader.getConfigs().size() > 1)
      filePost = getFileNameMotStr(usedMots);

    std::string motStr = pfaedle::router::getMotStr(usedMots);
    LOG(INFO) << "Calculating shapes for mots " << motStr;

    try {
      pfaedle::router::FeedStops fStops =
          pfaedle::router::writeMotStops(&gtfs[0], usedMots, cfg.shapeTripId);

      pfaedle::osm::Restrictor restr;
      pfaedle::trgraph::Graph graph;
      pfaedle::osm::OsmBuilder osmBuilder;

      pfaedle::osm::BBoxIdx box(BOX_PADDING);
      ShapeBuilder::getGtfsBox(&gtfs[0], cmdCfgMots, cfg.shapeTripId,
                               cfg.dropShapes, &box);

      if (!fStops.empty()) {
          osmBuilder.read(cfg.osmPath, motCfg.osmBuildOpts, &graph, box,
                          cfg.gridSize, &fStops, &restr);
      }

      // TODO(patrick): move this somewhere else
      for (auto& feedStop : fStops) {
        if (feedStop.second) {
          feedStop.second->pl().getSI()->getGroup()->writePens(
              motCfg.osmBuildOpts.trackNormzer,
              motCfg.routingOpts.platformUnmatchedPen,
              motCfg.routingOpts.stationDistPenFactor,
              motCfg.routingOpts.nonOsmPen);
        }
      }

      ShapeBuilder shapeBuilder(&gtfs[0], &evalFeed, cmdCfgMots, motCfg, &ecoll,
                                &graph, &fStops, &restr, cfg);

      if (cfg.writeGraph) {
        LOG(INFO) << "Outputting graph.json...";
        util::geo::output::GeoGraphJsonOutput out;
        mkdir(cfg.dbgOutputPath.c_str(), 775);
        std::ofstream fstr(cfg.dbgOutputPath + "/graph.json");
        out.printLatLng(*shapeBuilder.getGraph(), fstr);
        fstr.close();
      }

      if (singleTrip) {
        LOG(INFO) << "Outputting path.json...";
        mkdir(cfg.dbgOutputPath.c_str(), 775);
        std::ofstream pstr(cfg.dbgOutputPath + "/path.json");
        util::geo::output::GeoJsonOutput o(pstr);

        auto l = shapeBuilder.shapeL(singleTrip);

        // reproject to WGS84 to match RFC 7946
        o.printLatLng(l, {});

        o.flush();
        pstr.close();

        exit(static_cast<int>(RetCode::SUCCESS));
      }

      pfaedle::netgraph::Graph ng;
      shapeBuilder.shape(&ng);

      if (cfg.buildTransitGraph) {
        util::geo::output::GeoGraphJsonOutput out;
        LOG(INFO) << "Outputting trgraph" + filePost + ".json...";
        mkdir(cfg.dbgOutputPath.c_str(), 775);
        std::ofstream fstr(cfg.dbgOutputPath + "/trgraph" + filePost + ".json");
        out.printLatLng(ng, fstr);
        fstr.close();
      }
    } catch (const pfxml::parse_exc& ex) {
      LOG(ERROR) << "Could not parse OSM data, reason was:";
      std::cerr << ex.what() << std::endl;
      exit(static_cast<int>(RetCode::OSM_PARSE_ERR));
    }
  }

  if (cfg.evaluate)
      ecoll.printStats(&std::cout);

  if (!cfg.feedPaths.empty()) {
    try {
      mkdir(cfg.outputPath.c_str(), 775);
      LOG(INFO) << "Writing output GTFS to " << cfg.outputPath << " ...";
      pfaedle::gtfs::Writer w;
      w.write(&gtfs[0], cfg.outputPath);
    } catch (const ad::cppgtfs::WriterException& ex) {
      LOG(ERROR) << "Could not write final GTFS feed, reason was:";
      std::cerr << ex.what() << std::endl;
      exit(static_cast<int>(RetCode::GTFS_WRITE_ERR));
    }
  }

  return static_cast<int>(RetCode::SUCCESS);
}

// _____________________________________________________________________________
std::string getFileNameMotStr(const MOTs& mots) {
  std::string motStr;
  for (const auto& mot : mots) {
    motStr += "-" + ad::cppgtfs::gtfs::flat::Route::getTypeString(mot);
  }

  return motStr;
}

// _____________________________________________________________________________
std::vector<std::string> getCfgPaths(const Config& cfg) {
  if (!cfg.configPaths.empty()) {
      return cfg.configPaths;
  }
  std::vector<std::string> ret;


  // install prefix global configuration path, if available
  {
    auto path = std::string(INSTALL_PREFIX) +
                std::string(CFG_DIR) + "/" + "pfaedle" + "/" +
                CFG_FILE_NAME;
    std::ifstream is(path);

    LOG(DEBUG) << "Testing for config file at " << path;
    if (is.good()) {
      ret.push_back(path);
      LOG(DEBUG) << "Found implicit config file " << path;
    }
  }

  // local user configuration path, if available
  {
    auto path = util::getHomeDir() + CFG_HOME_SUFFIX + "/" +
                "pfaedle" + "/" + CFG_FILE_NAME;
    std::ifstream is(path);

    LOG(DEBUG) << "Testing for config file at " << path;
    if (is.good()) {
      ret.push_back(path);
      LOG(DEBUG) << "Found implicit config file " << path;
    }
  }

  // free this here, as we use homedir in the block above

  // CWD
  {
    char cwd[PATH_MAX];
    if (getcwd(cwd, sizeof(cwd))) {
      auto path = std::string(cwd) + "/" + CFG_FILE_NAME;
      std::ifstream is(path);

      LOG(DEBUG) << "Testing for config file at " << path;
      if (is.good()) {
        ret.push_back(path);
        LOG(DEBUG) << "Found implicit config file " << path;
      }
    }
  }

  return ret;
}

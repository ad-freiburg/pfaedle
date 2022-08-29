// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <limits.h>
#include <pwd.h>
#include <signal.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <clocale>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "ad/cppgtfs/Parser.h"
#include "ad/cppgtfs/Writer.h"
#include "pfaedle/config/ConfigReader.h"
#include "pfaedle/config/MotConfig.h"
#include "pfaedle/config/MotConfigReader.h"
#include "pfaedle/gtfs/Feed.h"
#include "pfaedle/gtfs/Writer.h"
#include "pfaedle/netgraph/Graph.h"
#include "pfaedle/osm/OsmIdSet.h"
#include "pfaedle/router/ShapeBuilder.h"
#include "pfaedle/router/Stats.h"
#include "pfaedle/statsimi-classifier/StatsimiClassifier.h"
#include "pfaedle/trgraph/Graph.h"
#include "util/Misc.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/geo/output/GeoJsonOutput.h"
#include "util/json/Writer.h"
#include "util/log/Log.h"

#ifndef CFG_HOME_SUFFIX
#define CFG_HOME_SUFFIX "/.config"
#endif
#ifndef CFG_DIR
#define CFG_DIR "/etc"
#endif
#ifndef CFG_FILE_NAME
#define CFG_FILE_NAME "pfaedle.cfg"
#endif

using configparser::ParseFileExc;
using pfaedle::config::Config;
using pfaedle::config::ConfigReader;
using pfaedle::config::MotConfig;
using pfaedle::config::MotConfigReader;
using pfaedle::osm::BBoxIdx;
using pfaedle::osm::OsmBuilder;
using pfaedle::router::DistDiffTransWeight;
using pfaedle::router::DistDiffTransWeightNoHeur;
using pfaedle::router::ExpoTransWeight;
using pfaedle::router::ExpoTransWeightNoHeur;
using pfaedle::router::MOTs;
using pfaedle::router::NormDistrTransWeight;
using pfaedle::router::NormDistrTransWeightNoHeur;
using pfaedle::router::Router;
using pfaedle::router::RouterImpl;
using pfaedle::router::ShapeBuilder;
using pfaedle::router::Stats;
using pfaedle::statsimiclassifier::BTSClassifier;
using pfaedle::statsimiclassifier::EDClassifier;
using pfaedle::statsimiclassifier::JaccardClassifier;
using pfaedle::statsimiclassifier::JaccardGeodistClassifier;
using pfaedle::statsimiclassifier::PEDClassifier;
using pfaedle::statsimiclassifier::StatsimiClassifier;

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
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());  // NOLINT

  // use utf8 locale
  std::setlocale(LC_ALL, "en_US.utf8");

  T_START(total);

  Config cfg;
  MotConfigReader motCfgReader;

  ConfigReader cr;
  cr.read(&cfg, argc, argv);

  std::vector<pfaedle::gtfs::Feed> gtfs(cfg.feedPaths.size());

  std::vector<std::string> cfgPaths = getCfgPaths(cfg);

  try {
    motCfgReader.parse(cfgPaths, cfg.motCfgParam);
  } catch (const configparser::ParseExc& ex) {
    LOG(ERROR) << "Could not parse MOT configurations, reason was:";
    std::cerr << ex.what() << std::endl;
    exit(static_cast<int>(RetCode::MOT_CFG_PARSE_ERR));
  }

  if (cfg.osmPath.empty() && !cfg.writeOverpass && !cfg.writeOsmfilter) {
    std::cerr << "No OSM input file specified (-x), see --help." << std::endl;
    exit(static_cast<int>(RetCode::NO_OSM_INPUT));
  }

  if (motCfgReader.getConfigs().size() == 0) {
    LOG(ERROR) << "No MOT configurations specified and no implicit "
                  "configurations found, see --help.";
    exit(static_cast<int>(RetCode::NO_MOT_CFG));
  }

  T_START(gtfsBuild);

  if (cfg.feedPaths.size() == 1) {
    if (cfg.inPlace) cfg.outputPath = cfg.feedPaths[0];
    if (!cfg.writeOverpass && !cfg.writeOsmfilter)
      LOG(INFO) << "Reading GTFS feed " << cfg.feedPaths[0] << " ...";
    try {
      ad::cppgtfs::Parser p;
      p.parse(&gtfs[0], cfg.feedPaths[0]);
    } catch (const ad::cppgtfs::ParserException& ex) {
      LOG(ERROR) << "Could not parse input GTFS feed, reason was:";
      std::cerr << ex.what() << std::endl;
      exit(static_cast<int>(RetCode::GTFS_PARSE_ERR));
    }
  } else if (cfg.writeOsm.size() || cfg.writeOverpass) {
    for (size_t i = 0; i < cfg.feedPaths.size(); i++) {
      if (!cfg.writeOverpass && !cfg.writeOsmfilter)
        LOG(INFO) << "Reading GTFS feed " << cfg.feedPaths[i] << " ...";
      ad::cppgtfs::Parser p;
      try {
        p.parse(&gtfs[i], cfg.feedPaths[i]);
      } catch (const ad::cppgtfs::ParserException& ex) {
        LOG(ERROR) << "Could not parse input GTFS feed, reason was:";
        std::cerr << ex.what() << std::endl;
        exit(static_cast<int>(RetCode::GTFS_PARSE_ERR));
      }
    }
  } else if (cfg.feedPaths.size() > 1) {
    std::cerr << "Multiple feeds only allowed in filter mode." << std::endl;
    exit(static_cast<int>(RetCode::MULT_FEEDS_NOT_ALWD));
  }

  auto tGtfsBuild = T_STOP(gtfsBuild);

  LOG(DEBUG) << "Read " << motCfgReader.getConfigs().size()
             << " unique MOT configs.";
  MOTs cmdCfgMots = cfg.mots;
  pfaedle::gtfs::Trip* singleTrip = 0;

  if (cfg.shapeTripId.size()) {
    if (!cfg.feedPaths.size()) {
      std::cout << "No input feed specified, see --help" << std::endl;
      exit(static_cast<int>(RetCode::NO_INPUT_FEED));
    }
    singleTrip = gtfs[0].getTrips().get(cfg.shapeTripId);
    if (!singleTrip) {
      LOG(ERROR) << "Trip #" << cfg.shapeTripId << " not found.";
      exit(static_cast<int>(RetCode::TRIP_NOT_FOUND));
    }
  }

  double maxSpeed = 0;
  for (const auto& c : motCfgReader.getConfigs()) {
    if (c.osmBuildOpts.maxSpeed > maxSpeed) {
      maxSpeed = c.osmBuildOpts.maxSpeed;
    }
  }

  if (cfg.writeOsm.size()) {
    LOG(INFO) << "Writing filtered XML to " << cfg.writeOsm << " ...";
    BBoxIdx box(BOX_PADDING);

    for (size_t i = 0; i < cfg.feedPaths.size(); i++) {
      ShapeBuilder::getGtfsBox(&gtfs[i], cmdCfgMots, cfg.shapeTripId, true,
                               &box, maxSpeed, 0);
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
                               &box, maxSpeed, 0);
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
  } else if (cfg.writeOsmfilter) {
    BBoxIdx box(BOX_PADDING);
    OsmBuilder osmBuilder;
    std::vector<pfaedle::osm::OsmReadOpts> opts;
    for (const auto& o : motCfgReader.getConfigs()) {
      if (std::find_first_of(o.mots.begin(), o.mots.end(), cmdCfgMots.begin(),
                             cmdCfgMots.end()) != o.mots.end()) {
        opts.push_back(o.osmBuildOpts);
      }
    }
    osmBuilder.osmfilterRuleWrite(&std::cout, opts, box);
    exit(static_cast<int>(RetCode::SUCCESS));
  } else if (!cfg.feedPaths.size()) {
    std::cout << "No input feed specified, see --help" << std::endl;
    exit(static_cast<int>(RetCode::NO_INPUT_FEED));
  }

  Stats stats;
  double tOsmBuild = 0;
  std::map<std::string, std::pair<size_t, size_t>> graphDimensions;
  std::vector<double> hopDists;

  for (const auto& motCfg : motCfgReader.getConfigs()) {
    std::string filePost;
    auto usedMots = pfaedle::router::motISect(motCfg.mots, cmdCfgMots);
    if (!usedMots.size()) continue;
    if (singleTrip && !usedMots.count(singleTrip->getRoute()->getType()))
      continue;
    if (motCfgReader.getConfigs().size() > 1)
      filePost = getFileNameMotStr(usedMots);

    std::string motStr = pfaedle::router::getMotStr(usedMots);
    LOG(INFO) << "Matching shapes for mots " << motStr;

    try {
      pfaedle::router::FeedStops fStops =
          pfaedle::router::writeMotStops(&gtfs[0], usedMots, cfg.shapeTripId);

      pfaedle::osm::Restrictor restr;
      pfaedle::trgraph::Graph graph;
      pfaedle::osm::OsmBuilder osmBuilder;

      pfaedle::osm::BBoxIdx box(BOX_PADDING);
      ShapeBuilder::getGtfsBox(&gtfs[0], usedMots, cfg.shapeTripId,
                               cfg.dropShapes, &box,
                               motCfg.osmBuildOpts.maxSpeed, &hopDists);

      T_START(osmBuild);

      if (fStops.size())
        osmBuilder.read(cfg.osmPath, motCfg.osmBuildOpts, &graph, box,
                        cfg.gridSize, &restr);

      tOsmBuild += T_STOP(osmBuild);
      graphDimensions[filePost].first = graph.getNds().size();

      for (const auto& nd : graph.getNds()) {
        graphDimensions[filePost].second += nd->getAdjListOut().size();
      }

      StatsimiClassifier* statsimiClassifier;

      if (motCfg.routingOpts.statsimiMethod == "bts") {
        statsimiClassifier = new BTSClassifier();
      } else if (motCfg.routingOpts.statsimiMethod == "jaccard") {
        statsimiClassifier = new JaccardClassifier();
      } else if (motCfg.routingOpts.statsimiMethod == "jaccard-geodist") {
        statsimiClassifier = new JaccardGeodistClassifier();
      } else if (motCfg.routingOpts.statsimiMethod == "ed") {
        statsimiClassifier = new EDClassifier();
      } else if (motCfg.routingOpts.statsimiMethod == "ped") {
        statsimiClassifier = new PEDClassifier();
      } else {
        LOG(ERROR) << "Unknown station similarity classifier "
                   << motCfg.routingOpts.statsimiMethod;
        exit(1);
      }

      Router* router = 0;

      if (motCfg.routingOpts.transPenMethod == "exp") {
        if (cfg.noAStar)
          router = new RouterImpl<ExpoTransWeightNoHeur>();
        else
          router = new RouterImpl<ExpoTransWeight>();
      } else if (motCfg.routingOpts.transPenMethod == "distdiff") {
        if (cfg.noAStar)
          router = new RouterImpl<DistDiffTransWeightNoHeur>();
        else
          router = new RouterImpl<DistDiffTransWeight>();
      } else if (motCfg.routingOpts.transPenMethod == "timenorm") {
        if (cfg.noAStar)
          router = new RouterImpl<NormDistrTransWeightNoHeur>();
        else
          router = new RouterImpl<NormDistrTransWeight>();
      } else {
        LOG(ERROR) << "Unknown routing method "
                   << motCfg.routingOpts.transPenMethod;
        exit(1);
      }

      ShapeBuilder shapeBuilder(&gtfs[0], usedMots, motCfg, &graph, &fStops,
                                &restr, statsimiClassifier, router, cfg);

      pfaedle::netgraph::Graph ng;

      if (singleTrip) {
        mkdir(cfg.dbgOutputPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::ofstream pstr(cfg.dbgOutputPath + "/path.json");
        util::geo::output::GeoJsonOutput o(pstr);

        auto l = shapeBuilder.shapeL(singleTrip);
        stats += l.second;

        LOG(INFO) << "Outputting path.json...";
        // reproject to WGS84 to match RFC 7946
        o.print(l.first, {});

        o.flush();
        pstr.close();
      } else {
        stats += shapeBuilder.shapeify(&ng);
      }

      if (router) delete router;
      if (statsimiClassifier) delete statsimiClassifier;

      if (cfg.writeGraph) {
        LOG(INFO) << "Outputting graph.json...";
        util::geo::output::GeoGraphJsonOutput out;
        mkdir(cfg.dbgOutputPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::ofstream fstr(cfg.dbgOutputPath + "/graph.json");
        out.print(*shapeBuilder.getGraph(), fstr);
        fstr.close();
      }

      if (singleTrip) exit(static_cast<int>(RetCode::SUCCESS));

      if (cfg.buildTransitGraph) {
        util::geo::output::GeoGraphJsonOutput out;
        LOG(INFO) << "Outputting trgraph-" + filePost + ".json...";
        mkdir(cfg.dbgOutputPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        std::ofstream fstr(cfg.dbgOutputPath + "/trgraph-" + filePost +
                           ".json");
        out.print(ng, fstr);
        fstr.close();
      }
    } catch (const pfxml::parse_exc& ex) {
      LOG(ERROR) << "Could not parse OSM data, reason was:";
      std::cerr << ex.what() << std::endl;
      exit(static_cast<int>(RetCode::OSM_PARSE_ERR));
    }
  }

  // outputting stats
  if (cfg.writeStats) {
    util::json::Dict graphSizes;

    double numNodesTot = 0;
    double numEdgesTot = 0;

    for (const auto& gd : graphDimensions) {
      util::json::Dict a;
      a["num_nodes"] = gd.second.first;
      a["num_edges"] = gd.second.second;
      numNodesTot += gd.second.first;
      numEdgesTot += gd.second.second;
      graphSizes[gd.first] = a;
    }

    double hopDistSum = 0;
    for (auto d : hopDists) hopDistSum += d;

    util::json::Dict jsonStats = {
        {"statistics",
         util::json::Dict{
             {"gtfs_num_stations", gtfs[0].getStops().size()},
             {"gtfs_num_trips", gtfs[0].getTrips().size()},
             {"gtfs_avg_hop_dist", hopDistSum / (hopDists.size() * 1.0)},
             {"graph_dimension", graphSizes},
             {"num_nodes_tot", numNodesTot},
             {"num_edges_tot", numEdgesTot},
             {"num_tries", stats.numTries},
             {"num_trie_leafs", stats.numTrieLeafs},
             {"dijkstra_iters", stats.dijkstraIters},
             {"time_solve", stats.solveTime},
             {"time_read_osm", tOsmBuild},
             {"time_read_gtfs", tGtfsBuild},
             {"time_tot", T_STOP(total)},
             {"peak-memory", util::readableSize(util::getPeakRSS())},
             {"peak-memory-bytes", util::getPeakRSS()}}}};

    std::ofstream ofs;
    ofs.open(cfg.dbgOutputPath + "/stats.json");
    util::json::Writer wr(&ofs, 10, true);
    wr.val(jsonStats);
    wr.closeAll();
  }

  if (cfg.feedPaths.size()) {
    try {
      mkdir(cfg.outputPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      LOG(INFO) << "Writing output GTFS to " << cfg.outputPath << " ...";
      pfaedle::gtfs::Writer w;
      w.write(&gtfs[0], cfg.outputPath);
    } catch (const ad::cppgtfs::WriterException& ex) {
      LOG(ERROR) << "Could not write output GTFS feed, reason was:";
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
    if (motStr.size()) motStr += "-";
    motStr += ad::cppgtfs::gtfs::flat::Route::getTypeString(mot);
  }

  return motStr;
}

// _____________________________________________________________________________
std::vector<std::string> getCfgPaths(const Config& cfg) {
  if (cfg.configPaths.size()) return cfg.configPaths;
  std::vector<std::string> ret;

  // install prefix global configuration path, if available
  {
    auto path = std::string(INSTALL_PREFIX) + std::string(CFG_DIR) + "/" +
                "pfaedle" + "/" + CFG_FILE_NAME;
    std::ifstream is(path);

    LOG(DEBUG) << "Testing for config file at " << path;
    if (is.good()) {
      ret.push_back(path);
      LOG(DEBUG) << "Found implicit config file " << path;
    }
  }

  // local user configuration path, if available
  {
    auto path = util::getHomeDir() + CFG_HOME_SUFFIX + "/" + "pfaedle" + "/" +
                CFG_FILE_NAME;
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

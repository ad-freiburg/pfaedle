// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <signal.h>
#include <stdio.h>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "ad/cppgtfs/Parser.h"
#include "ad/cppgtfs/Writer.h"
#include "ad/cppgtfs/gtfs/Feed.h"
#include "pfaedle/config/ConfigReader.h"
#include "pfaedle/config/MotConfig.h"
#include "pfaedle/config/MotConfigReader.h"
#include "pfaedle/eval/Collector.h"
#include "pfaedle/netgraph/Graph.h"
#include "pfaedle/osm/OsmIdSet.h"
#include "pfaedle/router/ShapeBuilder.h"
#include "pfaedle/trgraph/Graph.h"
#include "util/geo/output/GeoGraphJsonOutput.h"
#include "util/geo/output/GeoJsonOutput.h"
#include "util/json/JsonWriter.h"
#include "util/log/Log.h"

using std::string;
using pfaedle::router::MOTs;
using pfaedle::osm::BBoxIdx;
using pfaedle::osm::OsmBuilder;
using pfaedle::config::MotConfig;
using pfaedle::config::Config;
using pfaedle::router::ShapeBuilder;
using pfaedle::config::MotConfigReader;
using pfaedle::config::ConfigReader;
using pfaedle::eval::Collector;

std::string getMotStr(const MOTs& mots);
MOTs getContMots(const MotConfig& motCfg, const MOTs& mots);

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());  // NOLINT

  Config cfg;
  MotConfigReader motCfgReader;

  ConfigReader cr;
  cr.read(&cfg, argc, argv);

  ad::cppgtfs::gtfs::Feed gtfs;

  motCfgReader.parse(cfg.configPaths);

  if (cfg.feedPaths.size() == 1) {
    LOG(INFO) << "Reading " << cfg.feedPaths[0] << " ...";
    ad::cppgtfs::Parser p;
    p.parse(&gtfs, cfg.feedPaths[0]);
    LOG(INFO) << "Done.";
  } else if (cfg.feedPaths.size() > 1) {
    LOG(ERROR) << "Maximal one input feed allowed.";
    exit(1);
  }

  LOG(DEBUG) << "Read " << motCfgReader.getConfigs().size()
             << " unique MOT configs.";
  MOTs cmdCfgMots = cfg.mots;
  ad::cppgtfs::gtfs::Trip* singleTrip = 0;

  if (cfg.shapeTripId.size()) {
    singleTrip = gtfs.getTrips().get(cfg.shapeTripId);
    if (!singleTrip) {
      LOG(ERROR) << "Trip #" << cfg.shapeTripId << " not found.";
      exit(1);
    }
  }

  if (cfg.writeOsm.size()) {
    LOG(INFO) << "Writing filtered XML to " << cfg.writeOsm << " ...";
    BBoxIdx box(2500);
    if (cfg.feedPaths.size()) {
      box = ShapeBuilder::getPaddedGtfsBox(&gtfs, 2500, cmdCfgMots,
                                           cfg.shapeTripId);
    }
    OsmBuilder osmBuilder;
    std::vector<pfaedle::osm::OsmReadOpts> opts;
    for (const auto& o : motCfgReader.getConfigs()) {
      if (std::find_first_of(o.mots.begin(), o.mots.end(), cmdCfgMots.begin(),
                             cmdCfgMots.end()) != o.mots.end()) {
        opts.push_back(o.osmBuildOpts);
      }
    }
    osmBuilder.filterWrite(cfg.osmPath, cfg.writeOsm, opts, box);
    exit(0);
  }

  std::vector<double> dfBins;
  auto dfBinStrings = util::split(std::string(cfg.evalDfBins), ',');
  for (auto st : dfBinStrings) dfBins.push_back(atof(st.c_str()));
  Collector ecoll(cfg.evalPath, dfBins);

  for (const auto& motCfg : motCfgReader.getConfigs()) {
    auto usedMots = getContMots(motCfg, cmdCfgMots);
    if (!usedMots.size()) continue;

    std::string motStr = getMotStr(usedMots);
    LOG(INFO) << "Calculating shapes for mots " << motStr;

    ShapeBuilder shapeBuilder(&gtfs, cmdCfgMots, motCfg, &ecoll, cfg);

    if (cfg.writeGraph) {
      LOG(INFO) << "Outputting graph.json...";
      util::geo::output::GeoGraphJsonOutput out;
      std::ofstream fstr(cfg.dbgOutputPath + "/graph.json");
      out.print(*shapeBuilder.getGraph(), fstr);
      fstr.close();
    }

    if (singleTrip) {
      LOG(INFO) << "Outputting path.json...";
      std::ofstream pstr(cfg.dbgOutputPath + "/path.json");
      util::geo::output::GeoJsonOutput o(pstr);

      auto l = shapeBuilder.shapeL(singleTrip);

      if (singleTrip->getShape()) {
        auto orig = Collector::getWebMercLine(singleTrip->getShape(), -1, -1);
        o.print(orig, {{"ver", "old"}});
      }

      o.print(l, {{"ver", "new"}});
      o.flush();
      pstr.close();

      exit(0);
    }

    pfaedle::netgraph::Graph ng;
    shapeBuilder.shape(&ng);

    if (cfg.buildTransitGraph) {
      LOG(INFO) << "Outputting trgraph.json...";
      util::geo::output::GeoGraphJsonOutput out;
      std::ofstream fstr(cfg.dbgOutputPath + "/trgraph.json");
      out.print(ng, fstr);
      fstr.close();
    }
  }

  if (cfg.evaluate) ecoll.printStats(&std::cout);

  if (cfg.feedPaths.size()) {
    LOG(INFO) << "Writing output GTFS to " << cfg.outputPath << " ...";
    ad::cppgtfs::Writer w;
    w.write(&gtfs, cfg.outputPath);
  }

  return (0);
}

// _____________________________________________________________________________
std::string getMotStr(const MOTs& mots) {
  bool first = false;
  std::string motStr;
  for (const auto& mot : mots) {
    if (first) motStr += ", ";
    motStr += "<" + Route::getTypeString(mot) + ">";
    first = true;
  }

  return motStr;
}

// _____________________________________________________________________________
MOTs getContMots(const MotConfig& motCfg, const MOTs& mots) {
  MOTs ret;
  for (const auto& mot : mots) {
    if (motCfg.mots.count(mot)) {
      ret.insert(mot);
    }
  }

  return ret;
}

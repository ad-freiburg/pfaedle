// Copyright 2020, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <limits.h>
#include <stdlib.h>
#include <atomic>
#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include "ad/cppgtfs/Parser.h"
#include "shapevl/Collector.h"
#include "util/Misc.h"
#include "util/json/Writer.h"
#include "util/log/Log.h"

std::atomic<int> count(0);

// _____________________________________________________________________________
void printHelp(int argc, char** argv) {
  UNUSED(argc);
  std::cout << "Usage: " << argv[0]
            << " [-f <reportpath>] -g <gtfs> [-s] <test feeds>"
            << "\n";
  std::cout
      << "\nAllowed arguments:\n    -g <gtfs>     Ground truth GTFS file\n";
  std::cout << "     -s           Only output summary\n";
  std::cout << "     --json       Output JSON\n";
  std::cout << "     -f <folder>  Output full reports (per feed) to <folder>\n";
  std::cout
      << "     -m           MOTs to match (GTFS MOT or string, default: all)\n";
}

// _____________________________________________________________________________
void eval(const std::vector<std::string>* paths,
          std::vector<pfaedle::eval::Collector>* colls,
          const std::set<Route::TYPE>* mots,
          const ad::cppgtfs::gtfs::Feed* evalFeed) {
  while (1) {
    int myFeed = count-- - 1;
    if (myFeed < 0) return;
    std::string path = (*paths)[myFeed];
    LOG(DEBUG) << "Reading eval feed " << path << " "
               << " ...";
    ad::cppgtfs::gtfs::Feed feed;
    ad::cppgtfs::Parser p;

    try {
      p.parse(&feed, path);
    } catch (const ad::cppgtfs::ParserException& ex) {
      LOG(ERROR) << "Could not parse GTFS feed " << path << ", reason was:";
      std::cerr << ex.what() << std::endl;
      exit(1);
    }

    LOG(DEBUG) << "Evaluating " << path << "...";
    for (const auto& oldTrip : evalFeed->getTrips()) {
      if (!mots->count(oldTrip.second->getRoute()->getType())) continue;
      auto newTrip = feed.getTrips().get(oldTrip.first);
      if (!newTrip) {
        LOG(ERROR) << "Trip #" << oldTrip.first << " not present in " << path
                   << ", skipping...";
        continue;
      }
      (*colls)[myFeed].add(oldTrip.second, oldTrip.second->getShape(), newTrip,
                           newTrip->getShape());
    }
  }
}

// _____________________________________________________________________________
int main(int argc, char** argv) {
  // disable output buffering for standard output
  setbuf(stdout, NULL);

  // initialize randomness
  srand(time(NULL) + rand());  // NOLINT

  std::string groundTruthFeedPath, motStr;
  motStr = "all";
  ad::cppgtfs::gtfs::Feed groundTruthFeed;
  std::string fullReportPath = "";
  std::vector<std::string> evlFeedPaths;
  std::set<std::string> evlFeedPathsUniq;
  std::vector<pfaedle::eval::Collector> evalColls;
  std::vector<std::ofstream> reportStreams;
  bool summarize = false;
  bool json = false;

  for (int i = 1; i < argc; i++) {
    std::string cur = argv[i];
    if (cur == "-h" || cur == "--help") {
      printHelp(argc, argv);
      exit(0);
    } else if (cur == "-g") {
      if (++i >= argc) {
        LOG(ERROR) << "Missing argument for ground truth (-g).";
        exit(1);
      }
      groundTruthFeedPath = argv[i];
    } else if (cur == "-s") {
      summarize = true;
    } else if (cur == "--json") {
      json = true;
    } else if (cur == "-f") {
      if (++i >= argc) {
        LOG(ERROR) << "Missing argument for full reports (-f).";
        exit(1);
      }
      fullReportPath = argv[i];
    } else if (cur == "-m") {
      if (++i >= argc) {
        LOG(ERROR) << "Missing argument for mot (-m).";
        exit(1);
      }
      motStr = argv[i];
    } else {
      char fullPath[PATH_MAX + 1];
      if (!realpath(cur.c_str(), fullPath)) {
        LOG(ERROR) << "Error while reading " << fullPath;
        exit(1);
      }
      evlFeedPathsUniq.insert(fullPath);
    }
  }

  for (const auto& feedPath : evlFeedPathsUniq) {
    evlFeedPaths.push_back(feedPath);
    if (fullReportPath.size()) {
      reportStreams.emplace_back();
      reportStreams.back().exceptions(std::ios::failbit | std::ios::badbit);
      reportStreams.back().open(fullReportPath + "/" +
                                util::split(feedPath, '/').back() +
                                ".fullreport.tsv");
      evalColls.push_back({&reportStreams.back()});
    } else {
      evalColls.push_back({0});
    }
    count++;
  }

  if (groundTruthFeedPath.size() == 0) {
    LOG(ERROR) << "No ground truth feed path given (-g).";
    exit(1);
  }

  std::set<Route::TYPE> mots =
      ad::cppgtfs::gtfs::flat::Route::getTypesFromString(util::trim(motStr));

  std::vector<ad::cppgtfs::gtfs::Feed> evlFeeds(evlFeedPaths.size());

  try {
    LOG(DEBUG) << "Reading ground truth feed" << groundTruthFeedPath << " ...";
    ad::cppgtfs::Parser p;
    p.parse(&groundTruthFeed, groundTruthFeedPath);
  } catch (const ad::cppgtfs::ParserException& ex) {
    LOG(ERROR) << "Could not parse input GTFS feed, reason was:";
    std::cerr << ex.what() << std::endl;
    exit(1);
  }

  size_t THREADS = std::thread::hardware_concurrency();

  std::vector<std::thread> thrds(THREADS);
  for (auto& thr : thrds)
    thr =
        std::thread(&eval, &evlFeedPaths, &evalColls, &mots, &groundTruthFeed);

  for (auto& thr : thrds) thr.join();

  if (json) {
    util::json::Dict stats = {};

    for (size_t i = 0; i < evalColls.size(); i++) {
      stats[evlFeedPaths[i]] = evalColls[i].getJSONStats();
    }

    util::json::Dict jsonStats;

    if (evalColls.size() == 1) {
      jsonStats = {
          {"statistics", stats[evlFeedPaths[0]]
      }};
    } else {
      jsonStats = {
          {"statistics", stats
      }};
    }

    util::json::Writer wr(&std::cout, 10, true);
    wr.val(jsonStats);
    wr.closeAll();
  } else {
    for (size_t i = 0; i < evalColls.size(); i++) {
      if (summarize) {
        std::cout << evlFeedPaths[i] << ": ";
        evalColls[i].printShortStats(&std::cout);
        std::cout << std::endl;
      } else {
        std::cout << " == Evaluation results for " << evlFeedPaths[i]
                  << " ===" << std::endl;
        evalColls[i].printStats(&std::cout);
      }
    }
  }
}

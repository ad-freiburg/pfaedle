// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <float.h>
#include <getopt.h>
#include <exception>
#include <iostream>
#include <string>
#include "pfaedle/Def.h"
#include "pfaedle/_config.h"
#include "pfaedle/config/ConfigReader.h"
#include "util/String.h"
#include "util/log/Log.h"

using pfaedle::config::ConfigReader;

using std::string;
using std::exception;
using std::vector;

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wstring-plus-int"
static auto YEAR = __DATE__ + 7;
#pragma clang diagnostic pop
static auto COPY =
        "University of Freiburg - Chair of Algorithms and Data Structures";
static const char* AUTHORS = "Patrick Brosi <brosi@informatik.uni-freiburg.de>";

// _____________________________________________________________________________
void ConfigReader::help(const char* bin) {
  std::cout << std::setfill(' ') << std::left << "pfaedle GTFS map matcher "
            << VERSION_FULL << "\n(built " << __DATE__ << " " << __TIME__
            << " with geometry precision <" << PFAEDLE_PRECISION_STR << ">)\n\n"
            << "(C) " << YEAR << " " << COPY << "\n"
            << "Authors: " << AUTHORS << "\n\n"
            << "Usage: " << bin
            << " -x <OSM FILE> <GTFS FEED>\n\n"
            << "Allowed options:\n\n"
            << "General:\n"
            << std::setw(35) << "  -v [ --version ]"
            << "print version\n"
            << std::setw(35) << "  -h [ --help ]"
            << "show this help message\n"
            << std::setw(35) << "  -D [ --drop-shapes ]"
            << "drop shapes already present in the feed and\n"
            << std::setw(35) << " "
            << "  recalculate them\n"
            << "\nInput:\n"
            << std::setw(35) << "  -c [ --config ] arg"
            << "pfaedle config file\n"
            << std::setw(35) << "  -i [ --input ] arg"
            << "gtfs feed(s), may also be given as positional\n"
            << std::setw(35) << " "
            << "  parameter (see usage)\n"
            << std::setw(35) << "  -x [ --osm-file ] arg"
            << "OSM xml input file\n"
            << std::setw(35) << "  -m [ --mots ] arg (=all)"
            << "MOTs to calculate shapes for, comma sep.,\n"
            << std::setw(35) << " "
            << "  either as string "
               "{all, tram | streetcar,\n"
            << std::setw(35) << " "
            << "  subway | metro, rail | train, bus,\n"
            << std::setw(35) << " "
            << "  ferry | boat | ship, cablecar, gondola,\n"
            << std::setw(35) << " "
            << "  funicular, coach} or as GTFS mot codes\n"
            << "\nOutput:\n"
            << std::setw(35) << "  -o [ --output ] arg (=gtfs-out)"
            << "GTFS output path\n"
            << std::setw(35) << "  -X [ --osm-out ] arg"
            << "if specified, a filtered OSM file will be\n"
            << std::setw(35) << " "
            << "  written to <arg>\n"
            << std::setw(35) << "  --inplace"
            << "overwrite input GTFS feed with output feed\n"
            << "\nDebug Output:\n"
            << std::setw(35) << "  -d [ --dbg-path ] arg (=.)"
            << "output path for debug files\n"
            << std::setw(35) << "  --write-trgraph"
            << "write transit graph as GeoJSON to\n"
            << std::setw(35) << " "
            << "  <dbg-path>/trgraph.json\n"
            << std::setw(35) << "  --write-graph"
            << "write routing graph as GeoJSON to\n"
            << std::setw(35) << " "
            << "  <dbg-path>/graph.json\n"
            << std::setw(35) << "  --write-cgraph"
            << "if -T is set, write combination graph as\n"
            << std::setw(35) << " "
            << "  GeoJSON to "
               "<dbg-path>/combgraph.json\n"
            << std::setw(35) << "  --method arg (=global)"
            << "matching method to use, either 'global'\n"
            << std::setw(35) << " "
            << "  (based on HMM), 'greedy' or "
               "'greedy2'\n"
            << std::setw(35) << "  --eval"
            << "evaluate existing shapes against matched\n"
            << std::setw(35) << " "
            << "  shapes and print results\n"
            << std::setw(35) << "  --eval-path arg (=.)"
            << "path for eval file output\n"
            << std::setw(35) << "  --eval-df-bins arg (= )"
            << "bins to use for d_f histogram, comma sep.\n"
            << std::setw(35) << " "
            << "  (e.g. 10,20,30,40)\n"
            << "\nMisc:\n"
            << std::setw(35) << "  -T [ --trip-id ] arg"
            << "Do routing only for trip <arg>, write result \n"
            << std::setw(35) << " "
            << "  to <dbg-path>/path.json\n"
            << std::setw(35) << "  --overpass"
            << "Output overpass query for matching OSM data\n"
            << std::setw(35) << "  --grid-size arg (=2000)"
            << "Grid cell size\n"
            << std::setw(35) << "  --use-route-cache"
            << "(experimental) cache intermediate routing\n"
            << std::setw(35) << " "
            << "  results\n";
}

// _____________________________________________________________________________
void ConfigReader::read(Config* cfg, int argc, char** argv) {
  std::string motStr = "all";
  bool printOpts = false;

  struct option ops[] = {{"output", required_argument, 0, 'o'},
                         {"input", required_argument, 0, 'i'},
                         {"config", required_argument, 0, 'c'},
                         {"osm-file", required_argument, 0, 'x'},
                         {"drop-shapes", required_argument, 0, 'D'},
                         {"mots", required_argument, NULL, 'm'},
                         {"grid-size", required_argument, 0, 'g'},
                         {"overpass", no_argument, 0, 'a'},
                         {"osm-out", required_argument, 0, 'X'},
                         {"trip-id", required_argument, 0, 'T'},
                         {"write-graph", no_argument, 0, 1},
                         {"write-cgraph", no_argument, 0, 2},
                         {"write-trgraph", no_argument, 0, 4},
                         {"method", required_argument, 0, 5},
                         {"eval", no_argument, 0, 3},
                         {"eval-path", required_argument, 0, 6},
                         {"eval-df-bins", required_argument, 0, 7},
                         {"dbg-path", required_argument, 0, 'd'},
                         {"version", no_argument, 0, 'v'},
                         {"help", no_argument, 0, 'h'},
                         {"inplace", no_argument, 0, 9},
                         {"use-route-cache", no_argument, 0, 8},
                         {0, 0, 0, 0}};

  char c;
  while ((c = getopt_long(argc, argv, ":o:hvi:c:x:Dm:g:X:T:d:p", ops, nullptr)) != -1) {
    switch (c) {
      case 1:
        cfg->writeGraph = true;
        break;
      case 2:
        cfg->writeCombGraph = true;
        break;
      case 3:
        cfg->evaluate = true;
        break;
      case 4:
        cfg->buildTransitGraph = true;
        break;
      case 5:
        cfg->solveMethod = optarg;
        break;
      case 6:
        cfg->evalPath = optarg;
        break;
      case 7:
        cfg->evalDfBins = optarg;
        break;
      case 8:
        cfg->useCaching = true;
        break;
      case 'o':
        cfg->outputPath = optarg;
        break;
      case 'i':
        cfg->feedPaths.push_back(optarg);
        break;
      case 'c':
        cfg->configPaths.push_back(optarg);
        break;
      case 'x':
        cfg->osmPath = optarg;
        break;
      case 'D':
        cfg->dropShapes = true;
        break;
      case 'm':
        motStr = optarg;
        break;
      case 'g':
        cfg->gridSize = atof(optarg);
        break;
      case 'X':
        cfg->writeOsm = optarg;
        break;
      case 'T':
        cfg->shapeTripId = optarg;
        break;
      case 'd':
        cfg->dbgOutputPath = optarg;
        break;
      case 'a':
        cfg->writeOverpass = true;
        break;
      case 9:
        cfg->inPlace = true;
        break;
      case 'v':
        std::cout << "pfaedle " << VERSION_FULL << " (built " << __DATE__ << " "
                  << __TIME__ << " with geometry precision <"
                  << PFAEDLE_PRECISION_STR << ">)\n"
                  << "(C) " << YEAR << " " << COPY << "\n"
                  << "Authors: " << AUTHORS << "\nGNU General Public "
                                               "License v3.0\n";
        exit(0);
      case 'p':
        printOpts = true;
        break;
      case 'h':
        help(argv[0]);
        exit(0);
      case ':':
        std::cerr << argv[optind - 1];
        std::cerr << " requires an argument" << std::endl;
        exit(1);
      case '?':
        std::cerr << argv[optind - 1];
        std::cerr << " option unknown" << std::endl;
        exit(1);
        break;
      default:
        std::cerr << "Error while parsing arguments" << std::endl;
        exit(1);
        break;
    }
  }

  for (int i = optind; i < argc; i++) cfg->feedPaths.push_back(argv[i]);

  auto v = util::split(motStr, ',');
  for (const auto& motStr : v) {
    const auto& mots =
        ad::cppgtfs::gtfs::flat::Route::getTypesFromString(util::trim(motStr));
    cfg->mots.insert(mots.begin(), mots.end());
  }

  if (printOpts)
    std::cout << "\nConfigured options:\n\n" << cfg->toString() << std::endl;
}

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
#include "pfaedle/config/PfaedleConfig.h"
#include "util/String.h"
#include "util/geo/Geo.h"
#include "util/log/Log.h"

using pfaedle::config::ConfigReader;

using std::exception;
using std::string;
using std::vector;

static const char* YEAR = &__DATE__[7];
static const char* COPY =
    "University of Freiburg - Chair of Algorithms and Data Structures";
static const char* AUTHORS = "Patrick Brosi <brosi@informatik.uni-freiburg.de>";

// _____________________________________________________________________________
void ConfigReader::help(const char* bin) {
  std::cout << std::setfill(' ') << std::left << "pfaedle GTFS map matcher "
            << VERSION_FULL << "\n(built " << __DATE__ << " " << __TIME__
            << " with geometry precision <" << PFDL_PREC_STR << ">)\n\n"
            << "(C) " << YEAR << " " << COPY << "\n"
            << "Authors: " << AUTHORS << "\n\n"
            << "Usage: " << bin << " -x <OSM FILE> <GTFS FEED>\n\n"
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
            << std::setw(35) << "  --write-colors"
            << "write matched route line colors, where missing\n"
            << "\nInput:\n"
            << std::setw(35) << "  -c [ --config ] arg"
            << "pfaedle config file\n"
            << std::setw(35) << "  -i [ --input ] arg"
            << "gtfs feed(s), may also be given as positional\n"
            << std::setw(35) << "  -F [ --keep-additional-gtfs-fields ] arg"
            << "keep additional non-standard feeds in GTFS input\n"
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
            << "  funicular, coach, mono-rail | monorail,\n"
            << std::setw(35) << " "
            << "  trolley | trolleybus | trolley-bus} or\n"
            << std::setw(35) << " "
            << "  as GTFS mot codes\n"
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
            << "\nMisc:\n"
            << std::setw(35) << "  -T [ --trip-id ] arg"
            << "Do routing only for trip <arg>, write result \n"
            << std::setw(35) << " "
            << "  to <dbg-path>/path.json\n"
            << std::setw(35) << "  --overpass"
            << "Output overpass query for matching OSM data\n"
            << std::setw(35) << "  --osmfilter"
            << "Output osmfilter filter rules for matching OSM data\n"
            << std::setw(35) << "  -g [ --grid-size ] arg (=2000)"
            << "Approx. grid cell size in meters\n"
            << std::setw(35) << "  -b [ --box-padding ] arg (=20000)"
            << "Padding of bounding box used to crop input OSM data in meters\n"
            << std::setw(35) << "  --no-fast-hops"
            << "Disable fast hops technique\n"
            << std::setw(35) << "  --no-a-star"
            << "Disable A* heuristic \n"
            << std::setw(35) << "  --no-trie"
            << "Disable trip tries \n"
            << std::setw(35) << "  --no-hop-cache"
            << "Disable hop cache \n"
            << std::setw(35) << "  --stats"
            << "write stats to stats.json\n"
            << std::setw(35) << "  -W [ --warn ]"
            << "enable verbose warning messages\n"
            << std::setw(35) << "  -P"
            << "additional parameter string (in cfg file format)\n";
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
                         {"box-padding", required_argument, 0, 'b'},
                         {"overpass", no_argument, 0, 'a'},
                         {"osmfilter", no_argument, 0, 'f'},
                         {"osm-out", required_argument, 0, 'X'},
                         {"trip-id", required_argument, 0, 'T'},
                         {"write-graph", no_argument, 0, 1},
                         {"write-trgraph", no_argument, 0, 4},
                         {"dbg-path", required_argument, 0, 'd'},
                         {"version", no_argument, 0, 'v'},
                         {"help", no_argument, 0, 'h'},
                         {"inplace", no_argument, 0, 9},
                         {"no-fast-hops", no_argument, 0, 10},
                         {"no-a-star", no_argument, 0, 11},
                         {"no-trie", no_argument, 0, 12},
                         {"write-colors", no_argument, 0, 13},
                         {"stats", no_argument, 0, 14},
                         {"no-hop-cache", no_argument, 0, 15},
                         {"gaussian-noise", required_argument, 0, 16},
                         {"warn", no_argument, 0, 'W'},
                         {"keep-additional-gtfs-fields", no_argument, 0, 'F'},
                         {0, 0, 0, 0}};

  int c;
  while ((c = getopt_long(argc, argv, ":o:hvi:c:x:Dm:g:X:T:d:pP:FWb:", ops, 0)) !=
         -1) {
    switch (c) {
      case 1:
        cfg->writeGraph = true;
        break;
      case 4:
        cfg->buildTransitGraph = true;
        break;
      case 10:
        cfg->noFastHops = true;
        break;
      case 11:
        cfg->noAStar = true;
        break;
      case 12:
        cfg->noTrie = true;
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
        cfg->gridSize = atof(optarg) / util::geo::M_PER_DEG;
        break;
      case 'b':
        cfg->boxPadding = atof(optarg);
        break;
      case 'X':
        cfg->writeOsm = optarg;
        break;
      case 'T':
        cfg->shapeTripId = optarg;
        break;
      case 'P':
        cfg->motCfgParam += std::string("\n") + optarg;
        break;
      case 'd':
        cfg->dbgOutputPath = optarg;
        break;
      case 'a':
        cfg->writeOverpass = true;
        break;
      case 'f':
        cfg->writeOsmfilter = true;
        break;
      case 9:
        cfg->inPlace = true;
        break;
      case 13:
        cfg->writeColors = true;
        break;
      case 14:
        cfg->writeStats = true;
        break;
      case 15:
        cfg->noHopCache = true;
        break;
      case 16:
        cfg->gaussianNoise = atof(optarg);
        break;
      case 'W':
        cfg->verbosity = 1;
        break;
      case 'F':
        cfg->parseAdditionalGTFSFields = true;
        break;
      case 'v':
        std::cout << "pfaedle " << VERSION_FULL << std::endl;
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

// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_CONFIG_CONFIGREADER_H_
#define PFAEDLE_CONFIG_CONFIGREADER_H_

#include <vector>
#include "pfaedle/config/PfaedleConfig.h"

namespace pfaedle {
namespace config {

class ConfigReader {
 public:
  static void read(Config* targetConfig, int argc, char** argv);
  static void help(const char* bin);
};
}
}
#endif  // PFAEDLE_CONFIG_CONFIGREADER_H_

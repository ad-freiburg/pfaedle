// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <algorithm>
#include <iostream>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include "pfaedle/trgraph/Normalizer.h"

using pfaedle::trgraph::Normalizer;

// _____________________________________________________________________________
Normalizer::Normalizer(const ReplRules& rules) : _rulesOrig(rules) {
  buildRules(rules);
}

// _____________________________________________________________________________
std::string Normalizer::operator()(std::string sn) const {
  auto i = _cache.find(sn);
  if (i != _cache.end()) return i->second;

  std::string ret = sn;
  for (auto rule : _rules) {
    std::string tmp;
    std::regex_replace(std::back_inserter(tmp), ret.begin(), ret.end(),
                       rule.first, rule.second,
                       std::regex_constants::format_sed);
    std::swap(ret, tmp);
  }

  std::transform(ret.begin(), ret.end(), ret.begin(), ::tolower);

  _cache[sn] = ret;

  return ret;
}

// _____________________________________________________________________________
bool Normalizer::operator==(const Normalizer& b) const {
  return _rulesOrig == b._rulesOrig;
}

// _____________________________________________________________________________
void Normalizer::buildRules(const ReplRules& rules) {
  for (auto rule : rules) {
    try {
      _rules.push_back(ReplRuleComp(
          std::regex(rule.first, std::regex::ECMAScript | std::regex::icase |
                                     std::regex::optimize),
          rule.second));
    } catch (const std::regex_error& e) {
      std::stringstream ss;
      ss << "'" << rule.first << "'"
         << ": " << e.what();
      throw std::runtime_error(ss.str());
    }
  }
}

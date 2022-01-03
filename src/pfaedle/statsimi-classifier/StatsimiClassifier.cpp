// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <codecvt>
#include <exception>
#include <locale>
#include "pfaedle/Def.h"
#include "pfaedle/statsimi-classifier/StatsimiClassifier.h"
#include "util/geo/Geo.h"

using pfaedle::statsimiclassifier::JaccardClassifier;

// _____________________________________________________________________________
bool JaccardClassifier::similar(const std::string& nameA, const POINT& posA,
                                const std::string& nameB,
                                const POINT& posB) const {
  UNUSED(posA);
  UNUSED(posB);
  return similar(nameA, nameB);
}

// _____________________________________________________________________________
bool JaccardClassifier::similar(const std::string& nameA,
                                const std::string& nameB) const {
  // hard similarity
  if (nameA == nameB) return true;

  return util::jaccardSimi(nameA, nameB) > 0.45;  // 0.45 from paper
}

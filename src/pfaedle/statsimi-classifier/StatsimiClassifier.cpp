// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <codecvt>
#include <exception>
#include <locale>
#include "pfaedle/Def.h"
#include "pfaedle/statsimi-classifier/StatsimiClassifier.h"
#include "util/geo/Geo.h"

using pfaedle::statsimiclassifier::BTSClassifier;
using pfaedle::statsimiclassifier::EDClassifier;
using pfaedle::statsimiclassifier::JaccardClassifier;
using pfaedle::statsimiclassifier::JaccardGeodistClassifier;
using pfaedle::statsimiclassifier::PEDClassifier;

// _____________________________________________________________________________
bool JaccardGeodistClassifier::similar(const std::string& nameA,
                                       const POINT& posA,
                                       const std::string& nameB,
                                       const POINT& posB) const {
  const double THRES_M =
      0.00815467271246994481;    // ln 2/85 from statsimi evaluation
  const double THRES_JACC = .5;  // from statsimi evaluation

  const double m = exp(-THRES_M * util::geo::haversine(posA, posB));
  double jacc = util::jaccardSimi(nameA, nameB);

  if (jacc > THRES_JACC)
    jacc = .5 + (jacc - THRES_JACC) / (2.0 * (1.0 - THRES_JACC));
  else
    jacc = jacc / (2.0 * THRES_JACC);

  return ((m + jacc) / 2.0) > 0.5;
}

// _____________________________________________________________________________
bool JaccardGeodistClassifier::similar(const std::string& nameA,
                                       const std::string& nameB) const {
  return util::jaccardSimi(nameA, nameB) > 0.45;  // 0.45 from statsimi paper
}

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
  return util::jaccardSimi(nameA, nameB) > 0.45;  // 0.45 from statsimi paper
}

// _____________________________________________________________________________
bool BTSClassifier::similar(const std::string& nameA, const POINT& posA,
                            const std::string& nameB, const POINT& posB) const {
  UNUSED(posA);
  UNUSED(posB);
  return similar(nameA, nameB);
}

// _____________________________________________________________________________
bool BTSClassifier::similar(const std::string& nameA,
                            const std::string& nameB) const {
  return util::btsSimi(nameA, nameB) > 0.85;  // 0.85 from statsimi paper
}

// _____________________________________________________________________________
bool EDClassifier::similar(const std::string& nameA, const POINT& posA,
                           const std::string& nameB, const POINT& posB) const {
  UNUSED(posA);
  UNUSED(posB);
  return similar(nameA, nameB);
}

// _____________________________________________________________________________
bool EDClassifier::similar(const std::string& nameA,
                           const std::string& nameB) const {
  double edSimi = 1.0 - ((util::editDist(nameA, nameB) * 1.0) /
                         fmax(nameA.size(), nameB.size()));
  return edSimi > 0.85;  // 0.85 from statsimi paper
}

// _____________________________________________________________________________
bool PEDClassifier::similar(const std::string& nameA, const POINT& posA,
                            const std::string& nameB, const POINT& posB) const {
  UNUSED(posA);
  UNUSED(posB);
  return similar(nameA, nameB);
}

// _____________________________________________________________________________
bool PEDClassifier::similar(const std::string& nameA,
                            const std::string& nameB) const {
  double a = (util::prefixEditDist(nameA, nameB) * 1.0) / (nameA.size() * 1.0);
  double b = (util::prefixEditDist(nameB, nameA) * 1.0) / (nameB.size() * 1.0);
  double pedSimi = 1.0 - fmin(a, b);
  return pedSimi > 0.875;  // 0.875 average of values from statsimi paper
}

// Copyright 2020, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_STATSIMI_CLASSIFIER_STATSIMICLASSIFIER_H_
#define PFAEDLE_STATSIMI_CLASSIFIER_STATSIMICLASSIFIER_H_

#include <string>
#include "pfaedle/Def.h"
#include "util/geo/Geo.h"

namespace pfaedle {
namespace statsimiclassifier {

class StatsimiClassifier {
 public:
  virtual bool similar(const std::string& nameA, const POINT& posA,
                       const std::string& nameB, const POINT& posB) const = 0;

  virtual bool similar(const std::string& nameA,
                       const std::string& nameB) const = 0;
};

class JaccardClassifier : public StatsimiClassifier {
 public:
  virtual bool similar(const std::string& nameA, const POINT& posA,
                       const std::string& nameB, const POINT& posB) const;
  virtual bool similar(const std::string& nameA,
                       const std::string& nameB) const;
};

}  // namespace statsimiclassifier
}  // namespace pfaedle

#endif  // PFAEDLE_STATSIMI_CLASSIFIER_STATSIMICLASSIFIER_H_

// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_GTFS_SERVICE_H_
#define PFAEDLE_GTFS_SERVICE_H_

#include <string>
#include "ad/cppgtfs/gtfs/Service.h"
#include "util/Misc.h"

namespace pfaedle {
namespace gtfs {

class Service {
 public:
  typedef std::string Ref;
  static std::string getId(Ref r) { return r; }

  explicit Service(const string& id) : _id(id) {}
  Service(const string& id, uint8_t serviceDays,
          ad::cppgtfs::gtfs::ServiceDate start,
          ad::cppgtfs::gtfs::ServiceDate end)
      : _id(id) {
    UNUSED(serviceDays);
    UNUSED(start);
    UNUSED(end);
  }

  const std::string& getId() const { return _id; }
  void addException(const ad::cppgtfs::gtfs::ServiceDate& d,
                    ad::cppgtfs::gtfs::Service::EXCEPTION_TYPE t) {
    UNUSED(d);
    UNUSED(t);
  }

 private:
  std::string _id;
};
}  // namespace gtfs
}  // namespace pfaedle

#endif  // PFAEDLE_GTFS_SERVICE_H_

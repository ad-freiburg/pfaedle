// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_GTFS_ROUTE_H_
#define PFAEDLE_GTFS_ROUTE_H_

#include <stdint.h>
#include <algorithm>
#include <iomanip>
#include <set>
#include <sstream>
#include <string>
#include "ad/cppgtfs/gtfs/Agency.h"
#include "ad/cppgtfs/gtfs/Route.h"
#include "util/Misc.h"

using std::exception;
using std::string;

namespace pfaedle {
namespace gtfs {

class Route {
 public:
  typedef Route* Ref;
  static std::string getId(Ref r) { return r->getId(); }

  Route() {}

  Route(const string& id, ad::cppgtfs::gtfs::Agency* agency,
        const string& short_name, const string& long_name, const string& desc,
        ad::cppgtfs::gtfs::flat::Route::TYPE type, const string& url,
        uint32_t color, uint32_t text_color)
      : _id(id), _short_name(short_name), _long_name(long_name), _type(type) {
    UNUSED(agency);
    UNUSED(desc);
    UNUSED(url);
    UNUSED(color);
    UNUSED(text_color);
  }

  const std::string& getId() const { return _id; }

  const std::string& getShortName() const { return _short_name; }

  const std::string& getLongName() const { return _long_name; }

  ad::cppgtfs::gtfs::flat::Route::TYPE getType() const { return _type; }

 private:
  string _id;
  string _short_name;
  string _long_name;
  ad::cppgtfs::gtfs::flat::Route::TYPE _type;
};

}  // namespace gtfs
}  // namespace pfaedle

#endif  // PFAEDLE_GTFS_ROUTE_H_

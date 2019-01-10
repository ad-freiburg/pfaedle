// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_GTFS_WRITER_H_
#define PFAEDLE_GTFS_WRITER_H_

#include <string>
#include "Feed.h"

namespace pfaedle {
namespace gtfs {

class Writer {
 public:
  Writer() {}

  bool write(Feed* sourceFeed, const std::string& path) const;

 private:
  bool writeFeedInfo(Feed* f, std::ostream* os) const;
  bool writeAgency(Feed* f, std::ostream* os) const;
  bool writeStops(Feed* f, std::ostream* os) const;
  bool writeRoutes(Feed* f, std::ostream* os) const;
  bool writeCalendar(Feed* f, std::ostream* os) const;
  bool writeCalendarDates(Feed* f, std::ostream* os) const;
  bool writeFrequencies(Feed* f, std::ostream* os) const;
  bool writeTransfers(Feed* f, std::ostream* os) const;
  bool writeFares(Feed* f, std::ostream* os) const;
  bool writeFareRules(Feed* f, std::ostream* os) const;
  bool writeShapes(Feed* f, std::ostream* os) const;
  bool writeTrips(Feed* f, std::ostream* os) const;
  bool writeStopTimes(Feed* f, std::ostream* os) const;

  static void cannotWrite(const std::string& file, const std::string& file2);
};

}  // namespace gtfs
}  // namespace pfaedle

#endif  // PFAEDLE_GTFS_WRITER_H_

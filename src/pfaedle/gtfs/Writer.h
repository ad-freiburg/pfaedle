// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_GTFS_WRITER_H_
#define PFAEDLE_GTFS_WRITER_H_

#include <memory>
#include <string>
#ifdef LIBZIP_FOUND
#include <zip.h>
#endif

#include "Feed.h"
#include "ad/cppgtfs/Parser.h"
#include "ad/cppgtfs/Writer.h"

namespace pfaedle {
namespace gtfs {

class Writer {
 public:
  Writer() {}

  void write(Feed* sourceFeed, const std::string& path) const;

 private:
  void writeFeedInfo(Feed* f, std::ostream* os) const;
  void writeAgency(Feed* f, std::ostream* os) const;
  void writeStops(Feed* f, std::ostream* os) const;
  void writeRoutes(Feed* f, std::ostream* os) const;
  void writeCalendar(Feed* f, std::ostream* os) const;
  void writeCalendarDates(Feed* f, std::ostream* os) const;
  void writeFrequencies(Feed* f, std::ostream* os) const;
  void writeTransfers(Feed* f, std::ostream* os) const;
  void writeFares(Feed* f, std::ostream* os) const;
  void writeFareRules(Feed* f, std::ostream* os) const;
  void writeShapes(Feed* f, std::ostream* os) const;
  bool writeTrips(Feed* f, std::ostream* os) const;
  void writeStopTimes(Feed* f, std::ostream* os) const;
  void writeLevels(Feed* f, std::ostream* os) const;
  void writePathways(Feed* f, std::ostream* os) const;
  void writeAttribution(Feed* f, std::ostream* os) const;

  static void cannotWrite(const std::string& file, const std::string& file2);
  static void cannotWrite(const std::string& file);

#ifdef LIBZIP_FOUND
  static void moveIntoZip(zip* za, const std::string& sourcePath,
                          const std::string& targetPath);
#endif

  mutable std::ifstream _ifs;
};

}  // namespace gtfs
}  // namespace pfaedle

#endif  // PFAEDLE_GTFS_WRITER_H_

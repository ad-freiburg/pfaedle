// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_TRGRAPH_STATINFO_H_
#define PFAEDLE_TRGRAPH_STATINFO_H_

#include <string>
#include <unordered_map>
#include <vector>

namespace pfaedle {
namespace trgraph {

/*
 * Meta information (name, alternative names, track, ...) of a single stop
 */
class StatInfo {
 public:
  StatInfo();
  StatInfo(const StatInfo& si);
  StatInfo(const std::string& name, const std::string& track);

  // Return this stops names.
  const std::string& getName() const;

  // Return this stops track or empty string, if none.
  const std::string& getTrack() const;

  // Add an alternative name for this station.
  void addAltName(const std::string& name);

  // Return all alternative names for this station.
  const std::vector<std::string>& getAltNames() const;

  // Set the track of this stop.
  void setTrack(const std::string& tr);

  // Return the similarity between this stop and other
  double simi(const StatInfo* other) const;

#ifdef PFAEDLE_STATION_IDS
  const std::string& getId() const { return _id; }
  void setId(const std::string& id) { _id = id; }
#endif

 private:
  std::string _name;
  std::vector<std::string> _altNames;
  std::string _track;

#ifdef PFAEDLE_STATION_IDS
  // debug feature to store station ids from both OSM
  // and GTFS
  std::string _id;
#endif
};
}  // namespace trgraph
}  // namespace pfaedle

#endif  // PFAEDLE_TRGRAPH_STATINFO_H_

// Copyright 2017, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_TRGRAPH_STATINFO_H_
#define PFAEDLE_TRGRAPH_STATINFO_H_

#include <string>
#include <vector>
#include <unordered_map>

namespace pfaedle {
namespace trgraph {

// forward declaration
class StatGroup;

/*
 * Meta information (name, alternative names, track, group...) of a single stop
 */
class StatInfo {
 public:
  StatInfo();
  StatInfo(const StatInfo& si);
  StatInfo(const std::string& name, const std::string& track, bool _fromOsm);
  ~StatInfo();

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

  // Set this stations group.
  void setGroup(StatGroup* g);

  // Return this stations group.
  StatGroup* getGroup() const;

  // True if this stop was from osm
  bool isFromOsm() const;

  // Set this stop as coming from osm
  void setIsFromOsm(bool is);

 private:
  std::string _name;
  std::vector<std::string> _altNames;
  std::string _track;
  bool _fromOsm;
  StatGroup* _group;

  static std::unordered_map<const StatGroup*, size_t> _groups;
  static void unRefGroup(StatGroup* g);
};
}  // namespace trgraph
}  // namespace pfaedle

#endif  // PFAEDLE_TRGRAPH_STATINFO_H_

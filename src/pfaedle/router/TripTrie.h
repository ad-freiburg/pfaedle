// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_TRIPTRIE_H_
#define PFAEDLE_ROUTER_TRIPTRIE_H_

#include <vector>
#include <map>
#include <string>
#include "ad/cppgtfs/gtfs/Feed.h"
#include "pfaedle/gtfs/Feed.h"
#include "pfaedle/gtfs/StopTime.h"
#include "pfaedle/router/RoutingAttrs.h"

namespace pfaedle {
namespace router {

struct TripTrieNd {
  const ad::cppgtfs::gtfs::Stop* reprStop;
  std::string stopName;  // the stop name at this node
  std::string platform;  // the platform of node
  POINT pos;             // the position of this node
  double lat, lng;
  int time;
  bool arr;
  int accTime;
  size_t trips;
  size_t parent;
  std::vector<size_t> childs;
  RoutingAttrs rAttrs;
};

class TripTrie {
 public:
  // init node 0, this is the first decision node
  TripTrie() : _nds(1) {}
  bool addTrip(pfaedle::gtfs::Trip* trip, const RoutingAttrs& rAttrs,
               bool timeEx, bool degen);

  const std::vector<TripTrieNd>& getNds() const;
  const TripTrieNd& getNd(size_t nid) const;

  void toDot(std::ostream& os, const std::string& rootName, size_t gid) const;
  const std::map<size_t, std::vector<pfaedle::gtfs::Trip*>>& getNdTrips() const;

 private:
  std::vector<TripTrieNd> _nds;
  std::map<pfaedle::gtfs::Trip*, size_t> _tripNds;
  std::map<size_t, std::vector<pfaedle::gtfs::Trip*>> _ndTrips;

  bool add(pfaedle::gtfs::Trip* trip, const RoutingAttrs& rAttrs, bool timeEx);
  size_t get(pfaedle::gtfs::Trip* trip, bool timeEx);

  size_t getMatchChild(size_t parentNid, const std::string& stopName,
                       const std::string& platform, POINT pos, int time,
                       bool timeEx) const;
  size_t insert(const ad::cppgtfs::gtfs::Stop* stop, const RoutingAttrs& rAttrs,
                const POINT& pos, int time, bool arr, size_t parent);
};

}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_TRIPTRIE_H_

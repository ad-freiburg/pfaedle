// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <map>
#include <string>
#include <vector>

#include "TripTrie.h"
#include "ad/cppgtfs/gtfs/Feed.h"
#include "pfaedle/gtfs/Feed.h"
#include "pfaedle/gtfs/StopTime.h"

using pfaedle::gtfs::Trip;
using pfaedle::router::TripTrie;

// _____________________________________________________________________________
template <typename TRIP>
bool TripTrie<TRIP>::addTrip(TRIP* trip, const RoutingAttrs& rAttrs,
                             bool timeEx, bool degen) {
  if (!degen) return add(trip, rAttrs, timeEx);

  // check if trip is already fully and uniquely contained, if not, fail
  size_t existing = get(trip, timeEx);
  if (existing && _nds[existing].childs.size() == 0) {
    _tripNds[trip] = existing;
    _ndTrips[existing].push_back(trip);
    return true;
  } else {
    return false;
  }
}

// _____________________________________________________________________________
template <typename TRIP>
bool TripTrie<TRIP>::add(TRIP* trip, const RoutingAttrs& rAttrs, bool timeEx) {
  if (trip->getStopTimes().size() == 0) return false;

  int startSecs = 0;

  if (!trip->getStopTimes().front().getDepartureTime().empty()) {
    startSecs = trip->getStopTimes().front().getDepartureTime().seconds();
  }

  size_t curNdId = 0;
  for (size_t stId = 0; stId < trip->getStopTimes().size(); stId++) {
    const auto st = trip->getStopTimes()[stId];

    std::string name = st.getStop()->getName();
    std::string platform = st.getStop()->getPlatformCode();
    POINT pos = util::geo::latLngToWebMerc<PFDL_PREC>(st.getStop()->getLat(),
                                                      st.getStop()->getLng());

    if (stId > 0) {
      int arrTime = startSecs;

      if (!st.getArrivalTime().empty()) {
        arrTime = st.getArrivalTime().seconds() - startSecs;
      }

      size_t arrChild =
          getMatchChild(curNdId, name, platform, pos, arrTime, timeEx);

      if (arrChild) {
        curNdId = arrChild;

        _nds[arrChild].accTime += arrTime;
        _nds[arrChild].trips += 1;

        _nds[arrChild].rAttrs.merge(rAttrs);
      } else {
        curNdId = insert(st.getStop(), rAttrs, pos, arrTime, true, curNdId);
      }
    }

    if (stId < trip->getStopTimes().size() - 1) {
      int depTime = startSecs;

      if (!st.getDepartureTime().empty()) {
        depTime = st.getDepartureTime().seconds() - startSecs;
      }

      size_t depChild =
          getMatchChild(curNdId, name, platform, pos, depTime, timeEx);

      if (depChild) {
        curNdId = depChild;

        _nds[depChild].accTime += depTime;
        _nds[depChild].trips += 1;

        _nds[depChild].rAttrs.merge(rAttrs);
      } else {
        if (stId == 0 && _tripNds.size() > 0) return false;
        curNdId = insert(st.getStop(), rAttrs, pos, depTime, false, curNdId);
      }
    }
  }

  // curNdId is now the last matching node, insert the trip here
  _tripNds[trip] = curNdId;
  _ndTrips[curNdId].push_back(trip);

  return true;
}

// _____________________________________________________________________________
template <typename TRIP>
size_t TripTrie<TRIP>::get(TRIP* trip, bool timeEx) {
  if (trip->getStopTimes().size() == 0) return false;

  int startSecs = trip->getStopTimes().front().getDepartureTime().seconds();

  size_t curNdId = 0;
  for (size_t stId = 0; stId < trip->getStopTimes().size(); stId++) {
    const auto st = trip->getStopTimes()[stId];

    std::string name = st.getStop()->getName();
    std::string platform = st.getStop()->getPlatformCode();
    POINT pos = util::geo::latLngToWebMerc<PFDL_PREC>(st.getStop()->getLat(),
                                                      st.getStop()->getLng());

    if (stId > 0) {
      int arrTime = startSecs;

      if (!st.getArrivalTime().empty()) {
        arrTime = st.getArrivalTime().seconds() - startSecs;
      }

      size_t arrChild =
          getMatchChild(curNdId, name, platform, pos, arrTime, timeEx);

      if (arrChild) {
        curNdId = arrChild;
      } else {
        return 0;
      }
    }

    if (stId < trip->getStopTimes().size() - 1) {
      int depTime = startSecs;

      if (!st.getDepartureTime().empty()) {
        depTime = st.getDepartureTime().seconds() - startSecs;
      }

      size_t depChild =
          getMatchChild(curNdId, name, platform, pos, depTime, timeEx);

      if (depChild) {
        curNdId = depChild;
      } else {
        return 0;
      }
    }
  }

  return curNdId;
}

// _____________________________________________________________________________
template <typename TRIP>
size_t TripTrie<TRIP>::insert(const ad::cppgtfs::gtfs::Stop* stop,
                              const RoutingAttrs& rAttrs, const POINT& pos,
                              int time, bool arr, size_t parent) {
  _nds.emplace_back(TripTrieNd{stop,
                               stop->getName(),
                               stop->getPlatformCode(),
                               pos,
                               stop->getLat(),
                               stop->getLng(),
                               time,
                               arr,
                               time,
                               1,
                               parent,
                               {},
                               rAttrs});
  _nds[parent].childs.push_back(_nds.size() - 1);
  return _nds.size() - 1;
}

// _____________________________________________________________________________
template <typename TRIP>
const std::vector<pfaedle::router::TripTrieNd>& TripTrie<TRIP>::getNds() const {
  return _nds;
}

// _____________________________________________________________________________
template <typename TRIP>
size_t TripTrie<TRIP>::getMatchChild(size_t parentNid,
                                     const std::string& stopName,
                                     const std::string& platform, POINT pos,
                                     int time, bool timeEx) const {
  for (size_t child : _nds[parentNid].childs) {
    if (_nds[child].stopName == stopName && _nds[child].platform == platform &&
        util::geo::dist(_nds[child].pos, pos) < 1 &&
        (!timeEx || _nds[child].time == time)) {
      return child;
    }
  }

  return 0;
}

// _____________________________________________________________________________
template <typename TRIP>
void TripTrie<TRIP>::toDot(std::ostream& os, const std::string& rootName,
                           size_t gid) const {
  os << "digraph triptrie" << gid << " {";

  for (size_t nid = 0; nid < _nds.size(); nid++) {
    std::string color = "white";
    if (_ndTrips.count(nid)) color = "red";
    if (nid == 0) {
      os << "\"" << gid << ":0\" [label=\"" << rootName << "\"];\n";
    } else {
      os << "\"" << gid << ":" << nid
         << "\" [shape=\"box\" style=\"filled\" fillcolor=\"" << color
         << "\" label=\"#" << nid << ", " << _nds[nid].stopName << "@"
         << util::geo::getWKT(_nds[nid].pos) << " t=" << _nds[nid].time
         << "\"];\n";
    }
  }

  for (size_t nid = 0; nid < _nds.size(); nid++) {
    for (size_t child : _nds[nid].childs) {
      os << "\"" << gid << ":" << nid << "\" -> \"" << gid << ":" << child
         << "\";\n";
    }
  }

  os << "}";
}

// _____________________________________________________________________________
template <typename TRIP>
const std::map<size_t, std::vector<TRIP*>>& TripTrie<TRIP>::getNdTrips() const {
  return _ndTrips;
}

// _____________________________________________________________________________
template <typename TRIP>
const pfaedle::router::TripTrieNd& TripTrie<TRIP>::getNd(size_t nid) const {
  return _nds[nid];
}

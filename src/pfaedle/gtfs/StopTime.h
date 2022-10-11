// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_GTFS_STOPTIME_H_
#define PFAEDLE_GTFS_STOPTIME_H_

#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include "ad/cppgtfs/gtfs/Stop.h"
#include "ad/cppgtfs/gtfs/StopTime.h"
#include "util/Misc.h"

using std::exception;
using std::string;

namespace pfaedle {
namespace gtfs {

template <typename StopT>
class StopTime {
 public:
  StopTime(const ad::cppgtfs::gtfs::Time& at, const ad::cppgtfs::gtfs::Time& dt,
           typename StopT::Ref s, uint32_t seq, const std::string& hs,
           ad::cppgtfs::gtfs::flat::StopTime::PU_DO_TYPE put,
           ad::cppgtfs::gtfs::flat::StopTime::PU_DO_TYPE dot, float distTrav,
           bool isTp, uint8_t continuousDropOff,
           uint8_t continuousPickup)
      : _s(s), _sequence(seq), _dist(distTrav), _at(at), _dt(dt), _isTp(isTp) {
    UNUSED(hs);
    UNUSED(put);
    UNUSED(dot);
    UNUSED(distTrav);
    UNUSED(continuousDropOff);
    UNUSED(continuousPickup);
  }

  const typename StopT::Ref getStop() const { return _s; }
  typename StopT::Ref getStop() { return _s; }
  void setShapeDistanceTravelled(double d) { _dist = d; }

  ad::cppgtfs::gtfs::Time getArrivalTime() const {
    return _at;
  }
  ad::cppgtfs::gtfs::Time getDepartureTime() const {
    return _dt;
  }

  float getShapeDistanceTravelled() const { return _dist; }

  uint16_t getSeq() const { return _sequence; }
  bool isTp() const { return _isTp; }

 private:
  typename StopT::Ref _s;
  uint32_t _sequence;
  float _dist;
  ad::cppgtfs::gtfs::Time _at, _dt;
  bool _isTp;
};

template <typename StopTimeT>
struct StopTimeCompare {
  bool operator()(const StopTimeT& lh, const StopTimeT& rh) const {
    return lh.getSeq() < rh.getSeq();
  }
};

}  // namespace gtfs
}  // namespace pfaedle

#endif  // PFAEDLE_GTFS_STOPTIME_H_

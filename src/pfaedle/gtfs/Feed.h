// Copyright 2016, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_GTFS_FEED_H_
#define PFAEDLE_GTFS_FEED_H_

#include <string>
#include "Route.h"
#include "Service.h"
#include "ShapeContainer.h"
#include "StopTime.h"
#include "ad/cppgtfs/gtfs/ContContainer.h"
#include "ad/cppgtfs/gtfs/Feed.h"
#include "ad/cppgtfs/gtfs/NullContainer.h"
#include "ad/cppgtfs/gtfs/Stop.h"
#include "ad/cppgtfs/gtfs/StopTime.h"
#include "ad/cppgtfs/gtfs/Trip.h"

namespace pfaedle {
namespace gtfs {

typedef ad::cppgtfs::gtfs::FeedB<
    ad::cppgtfs::gtfs::Agency, Route, ad::cppgtfs::gtfs::Stop, Service,
    StopTime, Shape, ad::cppgtfs::gtfs::Fare, ad::cppgtfs::gtfs::Container,
    ad::cppgtfs::gtfs::ContContainer, ad::cppgtfs::gtfs::NullContainer,
    ad::cppgtfs::gtfs::ContContainer, ad::cppgtfs::gtfs::ContContainer,
    ShapeContainer, ad::cppgtfs::gtfs::NullContainer>
    Feed;
typedef ad::cppgtfs::gtfs::TripB<StopTime<ad::cppgtfs::gtfs::Stop>, Service,
                                 Route, Shape>
    Trip;

}  // namespace gtfs
}  // namespace pfaedle

#endif  // PFAEDLE_GTFS_FEED_H_

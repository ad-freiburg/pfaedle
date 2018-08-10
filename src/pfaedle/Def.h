// Copyright 2018
// Author: Patrick Brosi

#ifndef PFAEDLE_DEF_H_
#define PFAEDLE_DEF_H_

#include "util/geo/Geo.h"
#include "util/geo/PolyLine.h"

#define __str_a(s) __str_b(s)
#define __str_b(s) #s
#define __str_c(s) s ## 1
#define __str_d(s) __str_c(s)

#if !defined(PFAEDLE_PRECISION) || (__str_d(PFAEDLE_PRECISION) == 1)
#undef PFAEDLE_PRECISION
#define PFAEDLE_PRECISION double
#endif

#define PFAEDLE_PRECISION_STR __str_a(PFAEDLE_PRECISION)

// version number from cmake version module
#define POINT util::geo::Point<PFAEDLE_PRECISION>
#define LINE util::geo::Line<PFAEDLE_PRECISION>
#define BOX util::geo::Box<PFAEDLE_PRECISION>
#define POLYLINE util::geo::PolyLine<PFAEDLE_PRECISION>

#endif  // PFAEDLE_DEF_H_

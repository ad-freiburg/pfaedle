// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_OSM_BBOXIDX_H_
#define PFAEDLE_OSM_BBOXIDX_H_

#include <vector>
#include "pfaedle/Def.h"
#include "util/geo/Geo.h"

namespace pfaedle {
namespace osm {

using util::geo::Box;
using util::geo::Point;

struct BBoxIdxNd {
  BBoxIdxNd() : box(util::geo::minbox<double>()) {}
  explicit BBoxIdxNd(const Box<double>& box) : box(box) {}
  Box<double> box;
  std::vector<BBoxIdxNd> childs;
};

/*
 * Poor man's R-tree
 */
class BBoxIdx {
 public:
  explicit BBoxIdx(double padding);

  // Add a bounding box to this index
  void add(Box<double> box);

  // Check if a point is contained in this index
  bool contains(const Point<double>& box) const;

  // Return the full total bounding box of this index
  BOX getFullWebMercBox() const;

  // Return the size of this index
  size_t size() const;

 private:
  double _padding;
  size_t _size;

  BBoxIdxNd _root;

  void addToTree(const Box<double>& box, BBoxIdxNd* nd, size_t lvl);
  bool treeHas(const Point<double>& p, const BBoxIdxNd& nd) const;

  static const size_t MAX_LVL = 5;
  static constexpr double MIN_COM_AREA = 0.0;
};
}  // namespace osm
}  // namespace pfaedle

#endif  // PFAEDLE_OSM_BBOXIDX_H_

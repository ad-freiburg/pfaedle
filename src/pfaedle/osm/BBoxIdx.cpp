// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "pfaedle/osm/BBoxIdx.h"

using pfaedle::osm::BBoxIdx;

// _____________________________________________________________________________
BBoxIdx::BBoxIdx(double padding) : _padding(padding), _size(0) {}

// _____________________________________________________________________________
void BBoxIdx::add(Box<double> box) {
  // division by 83.000m is only correct here around a latitude deg of 25,
  // but should be a good heuristic. 1 deg is around 63km at latitude deg of 44,
  // and 110 at deg=0, since we usually dont do map matching in the arctic,
  // its okay to use 83km here.
  box = util::geo::pad(box, _padding / 83000);
  addToTree(box, &_root, 0);
  _size++;
}

// _____________________________________________________________________________
size_t BBoxIdx::size() const { return _size; }

// _____________________________________________________________________________
bool BBoxIdx::contains(const Point<double>& p) const {
  return treeHas(p, _root);
}

// _____________________________________________________________________________
BOX BBoxIdx::getFullWebMercBox() const {
  return BOX(
      util::geo::latLngToWebMerc<PFAEDLE_PRECISION>(
          _root.box.getLowerLeft().getY(), _root.box.getLowerLeft().getX()),
      util::geo::latLngToWebMerc<PFAEDLE_PRECISION>(
          _root.box.getUpperRight().getY(), _root.box.getUpperRight().getX()));
}

// _____________________________________________________________________________
BOX BBoxIdx::getFullBox() const { return _root.box; }

// _____________________________________________________________________________
std::vector<util::geo::Box<double>> BBoxIdx::getLeafs() const {
  std::vector<util::geo::Box<double>> ret;
  getLeafsRec(_root, &ret);
  return ret;
}

// _____________________________________________________________________________
void BBoxIdx::getLeafsRec(const BBoxIdxNd& nd,
                          std::vector<util::geo::Box<double>>* ret) const {
  if (!nd.childs.size()) {
    ret->push_back(nd.box);
    return;
  }

  for (const auto& child : nd.childs) {
    getLeafsRec(child, ret);
  }

  return;
}

// _____________________________________________________________________________
bool BBoxIdx::treeHas(const Point<double>& p, const BBoxIdxNd& nd) const {
  if (!nd.childs.size()) return util::geo::contains(p, nd.box);
  for (const auto& child : nd.childs) {
    if (util::geo::contains(p, child.box)) return treeHas(p, child);
  }

  return false;
}

// _____________________________________________________________________________
void BBoxIdx::addToTree(const Box<double>& box, BBoxIdxNd* nd, size_t lvl) {
  double bestCommonArea = 0;
  ssize_t bestChild = -1;

  // 1. update the bbox of this node
  nd->box = util::geo::extendBox(box, nd->box);

  if (lvl == MAX_LVL) return;

  // 2. find best candidate
  for (size_t i = 0; i < nd->childs.size(); i++) {
    double cur = util::geo::commonArea(box, nd->childs[i].box);
    if (cur > MIN_COM_AREA && cur > bestCommonArea) {
      bestChild = i;
      bestCommonArea = cur;
    }
  }

  if (bestChild < 0) {
    // 3. add a new node with the inserted bbox
    nd->childs.push_back(BBoxIdxNd(box));
    addToTree(box, &nd->childs.back(), lvl + 1);
  } else {
    // 3. add to best node
    addToTree(box, &nd->childs[bestChild], lvl + 1);
  }

  // TODO(patrick): some tree balancing by mergin overlapping bboxes in
  // non-leafs
}

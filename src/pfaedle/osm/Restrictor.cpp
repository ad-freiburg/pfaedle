// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include "pfaedle/osm/Restrictor.h"
#include "util/log/Log.h"

using pfaedle::osm::Restrictor;

// _____________________________________________________________________________
void Restrictor::relax(osmid wid, const trgraph::Node* n,
                       const trgraph::Edge* e) {
  // the wid is not unique here, because the OSM ways are split into
  // multiple edges. They are only unique as a pair with a "from-via" point.
  _rlx[NodeOsmIdP(n, wid)] = e;
  auto i = _posDangling.find(NodeOsmIdP(n, wid));
  if (i != _posDangling.end()) {
    for (const auto& path : i->second) {
      _pos[path.first][path.second].second = e;
      assert(path.first->hasEdge(e));
    }
  }

  auto j = _negDangling.find(NodeOsmIdP(n, wid));
  if (j != _negDangling.end()) {
    for (const auto& path : j->second) {
      _neg[path.first][path.second].second = e;
      assert(path.first->hasEdge(e));
    }
  }
}

// _____________________________________________________________________________
void Restrictor::add(const trgraph::Edge* from, osmid to,
                     const trgraph::Node* via, bool pos) {
  const trgraph::Edge* toE = 0;
  if (_rlx.count(NodeOsmIdP(via, to)))
    toE = _rlx.find(NodeOsmIdP(via, to))->second;
  if (pos) {
    _pos[via].push_back(RulePair(from, toE));
    if (!toE)
      _posDangling[NodeOsmIdP(via, to)].push_back(
          DanglPath(via, _pos[via].size() - 1));
  } else {
    _neg[via].push_back(RulePair(from, toE));
    if (!toE)
      _negDangling[NodeOsmIdP(via, to)].push_back(
          DanglPath(via, _neg[via].size() - 1));
  }
}

// _____________________________________________________________________________
bool Restrictor::may(const trgraph::Edge* from, const trgraph::Edge* to,
                     const trgraph::Node* via) const {
  auto posI = _pos.find(via);
  auto negI = _neg.find(via);

  if (posI != _pos.end()) {
    for (const auto& r : posI->second) {
      if (r.first == from && r.second && r.second != to)
        return false;
      else if (r.first == from && r.second == to)
        return true;
    }
  }
  if (negI != _neg.end()) {
    for (const auto& r : negI->second) {
      if (r.first == from && r.second == to) return false;
    }
  }
  return true;
}

// _____________________________________________________________________________
void Restrictor::replaceEdge(const trgraph::Edge* old,
                             const trgraph::Edge* newA,
                             const trgraph::Edge* newB) {
  const trgraph::Edge* newFrom;
  const trgraph::Edge* newTo;
  if (old->getFrom() == newA->getFrom() || old->getFrom() == newA->getTo()) {
    newFrom = newA;
    newTo = newB;
  } else {
    newFrom = newB;
    newTo = newA;
  }
  replaceEdge(old, old->getFrom(), newFrom);
  replaceEdge(old, old->getTo(), newTo);
}

// _____________________________________________________________________________
void Restrictor::duplicateEdge(const trgraph::Edge* old,
                               const trgraph::Edge* newE) {
  duplicateEdge(old, old->getFrom(), newE);
  duplicateEdge(old, old->getTo(), newE);
}

// _____________________________________________________________________________
void Restrictor::duplicateEdge(const trgraph::Edge* old,
                               const trgraph::Node* via,
                               const trgraph::Edge* newE) {
  auto posI = _pos.find(via);
  auto negI = _neg.find(via);

  assert(old->getFrom() == newE->getTo() && old->getTo() == newE->getFrom());

  if (posI != _pos.end()) {
    for (auto& r : posI->second) {
      if (r.first == old) {
        if (r.first->getTo() != via) {
          r.first = newE;
        }
      }
      if (r.second == old) {
        if (r.second->getFrom() != via) {
          r.second = newE;
        }
      }
    }
  }
  if (negI != _neg.end()) {
    for (auto& r : negI->second) {
      if (r.first == old) {
        if (r.first->getTo() != via) {
          r.first = newE;
        }
      }
      if (r.second == old) {
        if (r.second->getFrom() != via) {
          r.second = newE;
        }
      }
    }
  }
}

// _____________________________________________________________________________
void Restrictor::replaceEdge(const trgraph::Edge* old, const trgraph::Node* via,
                             const trgraph::Edge* newE) {
  auto posI = _pos.find(via);
  auto negI = _neg.find(via);

  if (posI != _pos.end()) {
    for (auto& r : posI->second) {
      if (r.first == old) r.first = newE;
      if (r.second == old) r.second = newE;
    }
  }
  if (negI != _neg.end()) {
    for (auto& r : negI->second) {
      if (r.first == old) r.first = newE;
      if (r.second == old) r.second = newE;
    }
  }
}

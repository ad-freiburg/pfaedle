// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_TRGRAPH_EDGEPL_H_
#define PFAEDLE_TRGRAPH_EDGEPL_H_

#include <map>
#include <set>
#include <string>
#include <vector>
#include "pfaedle/Def.h"
#include "pfaedle/router/Comp.h"
#include "util/geo/Geo.h"
#include "util/geo/GeoGraph.h"

using util::geograph::GeoEdgePL;



namespace pfaedle {
namespace trgraph {

/*
 * A line occuring on an edge
 */
struct TransitEdgeLine {
  std::string fromStr;
  std::string toStr;
  std::string shortName;
};

inline bool operator==(const TransitEdgeLine& a, const TransitEdgeLine& b) {
  return a.fromStr == b.fromStr && a.toStr == b.toStr &&
         a.shortName == b.shortName;
}

inline bool operator<(const TransitEdgeLine& a, const TransitEdgeLine& b) {
  return a.fromStr < b.fromStr ||
         (a.fromStr == b.fromStr && a.toStr < b.toStr) ||
         (a.fromStr == b.fromStr && a.toStr == b.toStr &&
          a.shortName < b.shortName);
}

/*
 * An edge payload class for the transit graph.
 */
class EdgePL {
 public:
  EdgePL();
  ~EdgePL();
  EdgePL(const EdgePL& pl);
  EdgePL(const EdgePL& pl, bool geoFlat);

  // Return the geometry of this edge.
  const LINE* getGeom() const;
  LINE* getGeom();

  // Extends this edge payload's geometry by Point p
  void addPoint(const POINT& p);

  // Fill obj with k/v pairs describing the parameters of this payload.
  util::json::Dict getAttrs() const;

  // Return the length in meters stored for this edge payload
  double getLength() const;

  // Set the length in meters for this edge payload
  void setLength(double d);

  // Set this edge as a one way node, either in the default direction of
  // the edge (no arg), or the direction specified in dir
  void setOneWay();
  void setOneWay(uint8_t dir);

  // Mark this payload' edge as having some restrictions
  void setRestricted();

  // Mark this payload' edge as being secondary to an inversed partner
  void setRev();

  // True if this edge is secondary to an inversed partner
  bool isRev() const;

  // True if this edge is restricted
  bool isRestricted() const;

  // Set the level of this edge.
  void setLvl(uint8_t lvl);

  // Return the level of this edge.
  uint8_t lvl() const;

  // Return the one-way code stored for this edge.
  uint8_t oneWay() const;

  // Add a TransitedgeLine to this payload's edge
  void addLine(const TransitEdgeLine* l);

  // Add multiple TransitedgeLine objects to this payload's edge
  void addLines(const std::vector<TransitEdgeLine*>& l);

  // Return the TransitEdgeLines stored for this payload
  const std::vector<const TransitEdgeLine*>& getLines() const;

  // Returns the last hop of the payload - this is the (n-2)th point in
  // the payload geometry of length n > 1
  const POINT& backHop() const;

  // Returns the first hop of the payload - this is the 2nd point in
  // the payload geometry of length n > 1
  const POINT& frontHop() const;

  // Obtain an exact copy of this edge, but in reverse.
  EdgePL revCopy() const;

 private:
  float _length;
  uint8_t _oneWay : 2;
  bool _hasRestr : 1;
  bool _rev : 1;
  uint8_t _lvl : 3;

  LINE* _l;

  std::vector<const TransitEdgeLine*> _lines;

  static void unRefTLine(const TransitEdgeLine* l);

  static std::map<LINE*, size_t> _flines;
  static std::map<const TransitEdgeLine*, size_t> _tlines;
};
}  // namespace trgraph
}  // namespace pfaedle

#endif  // PFAEDLE_TRGRAPH_EDGEPL_H_

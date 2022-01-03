// Copyright 2020, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_ROUTER_HOPCACHE_H_
#define PFAEDLE_ROUTER_HOPCACHE_H_

#include <map>
#include <set>
#include <utility>
#include "pfaedle/trgraph/Graph.h"
#include "util/Misc.h"

namespace pfaedle {
namespace router {

class HopCache {
 public:
  void setMin(const trgraph::Edge* a, const trgraph::Edge* b, uint32_t val);

  void setMin(const trgraph::Edge* a, const std::set<trgraph::Edge*>& b,
              uint32_t val);

  void setMin(const std::set<trgraph::Edge*>& a, const trgraph::Edge* b,
              uint32_t val);

  void setEx(const trgraph::Edge* a, const trgraph::Edge* b, uint32_t val);

  std::pair<uint32_t, bool> get(const trgraph::Edge* a,
                                const trgraph::Edge* b) const;

 private:
  util::SparseMatrix<const trgraph::Edge*, int64_t, 0> _cache;
};

}  // namespace router
}  // namespace pfaedle

#endif  // PFAEDLE_ROUTER_HOPCACHE_H_

// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_OSM_OSMIDSET_H_
#define PFAEDLE_OSM_OSMIDSET_H_

#include <fcntl.h>
#include <unistd.h>
#include <bitset>
#include <set>
#include <string>
#include <vector>
#include "pfaedle/osm/Osm.h"

namespace pfaedle {
namespace osm {

// buffer sizes _must_ be multiples of 8
static const size_t BUFFER_S = 8 * 64 * 1024;
static const size_t SORT_BUFFER_S = 8 * 64 * 1024;
static const size_t OBUFFER_S = 8 * 1024 * 1024;

#define BLOOMF_BITS 400000000

/*
 * A disk-based set for OSM ids. Read-access for checking the presence is
 * reduced by a bloom filter
 */
class OsmIdSet {
 public:
  OsmIdSet();
  ~OsmIdSet();

  // Add an OSM id
  void add(osmid id);

  // Check if an OSM id is contained
  bool has(osmid id) const;

  // Count the number of lookups and file lookups for debugging
  static size_t LOOKUPS;
  static size_t FLOOKUPS;

 private:
  std::set<osmid> _set;
  mutable bool _closed;
  mutable int _file;
  unsigned char* _buffer;
  unsigned char* _outBuffer;
  mutable bool _sorted;
  osmid _last;
  osmid _smallest;
  osmid _biggest;

  size_t _obufpos;
  mutable size_t _curBlock;
  mutable ssize_t _curBlockSize;

  // bloom filter
  std::bitset<BLOOMF_BITS>* _bitset;

  mutable std::vector<osmid> _blockEnds;

  mutable size_t _fsize;

  uint32_t knuth(uint32_t in) const;
  uint32_t jenkins(uint32_t in) const;
  uint32_t hash(uint32_t in, int i) const;
  void diskAdd(osmid id);
  void close() const;
  void sort() const;
  bool diskHas(osmid id) const;
  std::string getFName() const;
  size_t getBlock(osmid id) const;
  int openTmpFile() const;
  size_t cwrite(int f, const void* buf, size_t n) const;
  size_t cread(int f, void* buf, size_t n) const;

  static int qsortCmp(const void* a, const void* b) {
    if (*static_cast<const uint64_t*>(a) < *static_cast<const uint64_t*>(b))
      return -1;
    if (*static_cast<const uint64_t*>(a) > *static_cast<const uint64_t*>(b))
      return 1;
    return 0;
  }
};
}  // namespace osm
}  // namespace pfaedle

#endif  // PFAEDLE_OSM_OSMIDSET_H_

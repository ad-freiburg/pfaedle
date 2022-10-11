// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdio.h>
#include <unistd.h>
#include <algorithm>
#include <cassert>
#include <climits>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include "pfaedle/Def.h"
#include "pfaedle/osm/OsmIdSet.h"
#include "util/3rdparty/MurmurHash3.h"

using pfaedle::osm::OsmIdSet;

size_t OsmIdSet::LOOKUPS = 0;
size_t OsmIdSet::FLOOKUPS = 0;

// _____________________________________________________________________________
OsmIdSet::OsmIdSet()
    : _closed(false),
      _sorted(true),
      _last(0),
      _smallest(-1),
      _biggest(0),
      _hasInv(false),
      _obufpos(0),
      _curBlock(-1),
      _fsize(0) {
  _bitset = new std::bitset<BLOOMF_BITS>();
  _bitsetNotIn = new std::bitset<BLOOMF_BITS>();
  _file = openTmpFile();

  _buffer = new unsigned char[BUFFER_S];
  _outBuffer = new unsigned char[BUFFER_S];
}

// _____________________________________________________________________________
OsmIdSet::~OsmIdSet() {
  delete _bitset;
  delete _bitsetNotIn;
  delete[] _buffer;
  if (!_closed) delete[] _outBuffer;
}

// _____________________________________________________________________________
void OsmIdSet::nadd(osmid id) {
  if (_closed) throw std::exception();

  _hasInv = true;

  uint32_t h1, h2;
  MurmurHash3_x86_32(&id, 8, 469954432, &h1);
  h2 = jenkins(id);

  for (int i = 0; i < 5; i++) {
    uint32_t h = (h1 + i * h2) % BLOOMF_BITS;
    (*_bitsetNotIn)[h] = 1;
  }
}

// _____________________________________________________________________________
void OsmIdSet::add(osmid id) {
  if (_closed) throw std::exception();

  diskAdd(id);

  if (_last > id) _sorted = false;
  _last = id;
  if (id < _smallest) _smallest = id;
  if (id > _biggest) _biggest = id;

  uint32_t h1, h2;
  MurmurHash3_x86_32(&id, 8, 469954432, &h1);
  h2 = jenkins(id);

  for (int i = 0; i < 5; i++) {
    uint32_t h = (h1 + i * h2) % BLOOMF_BITS;
    (*_bitset)[h] = 1;
  }
}

// _____________________________________________________________________________
void OsmIdSet::diskAdd(osmid id) {
  memcpy(_outBuffer + _obufpos, &id, 8);

  _obufpos += 8;

  if (_obufpos % BUFFER_S == 0) {
    // this is the last value in this block
    _blockEnds.push_back(id);
  }

  if (_obufpos >= BUFFER_S) {
    ssize_t w = cwrite(_file, _outBuffer, BUFFER_S);
    _fsize += w;
    _obufpos = 0;
  }
}

// _____________________________________________________________________________
size_t OsmIdSet::getBlock(osmid id) const {
  auto it = std::upper_bound(_blockEnds.begin(), _blockEnds.end(), id);
  return (it - _blockEnds.begin());
}

// _____________________________________________________________________________
bool OsmIdSet::diskHas(osmid id) const {
  assert(_sorted);

  auto a = std::lower_bound(_blockEnds.begin(), _blockEnds.end(), id);
  if (a != _blockEnds.end() && *a == id) {
    return true;
  }

  size_t block = getBlock(id);

  if (block != _curBlock) {
    lseek(_file, block * BUFFER_S, SEEK_SET);

    ssize_t n = cread(_file, _buffer, BUFFER_S);
    _curBlockSize = n;
    FLOOKUPS++;
    _curBlock = block;
  }

  if (_curBlockSize <= 7) return false;
  if (*(reinterpret_cast<uint64_t*>(_buffer)) > id) return false;

  ssize_t l = 0;
  ssize_t r = _curBlockSize - 8;

  while (l <= r) {
    unsigned char* p = _buffer + (l + ((r - l) / 16) * 8);
    osmid cur = *(reinterpret_cast<uint64_t*>(p));
    if (cur == id) return true;
    if (cur < id)
      l = (p - _buffer) + 8;
    else
      r = (p - _buffer) - 8;
  }

  return false;
}

// _____________________________________________________________________________
bool OsmIdSet::has(osmid id) const {
  LOOKUPS++;
  if (!_closed) close();

  // trivial cases
  if (id < _smallest || id > _biggest) {
    return false;
  }

  uint32_t h1, h2;
  MurmurHash3_x86_32(&id, 8, 469954432, &h1);
  h2 = jenkins(id);

  for (int i = 0; i < 5; i++) {
    uint32_t h = (h1 + i * h2) % BLOOMF_BITS;
    if ((*_bitset)[h] == 0) {
      return false;
    }
    if (_hasInv && (*_bitsetNotIn)[h] == 0) {
      return true;
    }
  }

  bool has = diskHas(id);
  return has;
}

// _____________________________________________________________________________
void OsmIdSet::close() const {
  ssize_t w = cwrite(_file, _outBuffer, _obufpos);
  _fsize += w;
  _blockEnds.push_back(_biggest);
  delete[] _outBuffer;
  _closed = true;

  // if order was not sorted, sort now
  if (!_sorted) sort();
}

// _____________________________________________________________________________
void OsmIdSet::sort() const {
  // sort file via an external merge sort

  _blockEnds.clear();
  size_t parts = _fsize / SORT_BUFFER_S + 1;
  size_t partsBufSize = ((SORT_BUFFER_S / 8) / parts + 1) * 8;

  unsigned char* buf = new unsigned char[SORT_BUFFER_S];
  unsigned char** partbufs = new unsigned char*[parts];
  size_t* partpos = new size_t[parts];
  size_t* partsize = new size_t[parts];

  // sort the 'parts' number of file parts independently
  for (size_t i = 0; i < parts; i++) {
    partbufs[i] = new unsigned char[partsBufSize];
    partpos[i] = 0;
    partsize[i] = 0;
    lseek(_file, SORT_BUFFER_S * i, SEEK_SET);
    ssize_t n = read(_file, buf, SORT_BUFFER_S);
    if (n < 0) continue;
    qsort(buf, n / 8, 8, qsortCmp);
    lseek(_file, SORT_BUFFER_S * i, SEEK_SET);
    cwrite(_file, buf, n);

    memcpy(partbufs[i], buf, std::min<size_t>(n, partsBufSize));
    partsize[i] = n;
  }

  // now the individial parts are sorted
  int newFile = openTmpFile();

  for (size_t i = 0; i < _fsize; i += 8) {
    uint64_t smallest = UINT64_MAX;
    ssize_t smallestP = -1;

    // look for smallest element (not optimal, but running time is not
    // really critical here)
    for (size_t j = 0; j < parts; j++) {
      if (partpos[j] == partsize[j]) continue;  // bucket already empty
      if (*reinterpret_cast<uint64_t*>(
              &partbufs[j][partpos[j] % partsBufSize]) <= smallest) {
        smallestP = j;
        smallest = *reinterpret_cast<uint64_t*>(
            &partbufs[j][partpos[j] % partsBufSize]);
      }
    }

    assert(smallestP > -1);

    memcpy(buf + (i % SORT_BUFFER_S), &smallest, 8);

    if ((i + 8) % BUFFER_S == 0) _blockEnds.push_back(smallest);

    if ((i % SORT_BUFFER_S) == SORT_BUFFER_S - 8 || i == _fsize - 8) {
      // write to output file
      cwrite(newFile, buf, i % SORT_BUFFER_S + 8);
    }

    partpos[smallestP] += 8;

    if (partpos[smallestP] % partsBufSize == 0) {
      lseek(_file, SORT_BUFFER_S * smallestP + partpos[smallestP], SEEK_SET);
      cread(_file, partbufs[smallestP], partsBufSize);
    }
  }

  // cleanup
  delete[] buf;
  for (size_t j = 0; j < parts; j++) delete[] partbufs[j];
  delete[] partbufs;
  delete[] partpos;
  delete[] partsize;

  _file = newFile;
  _sorted = true;
}

// _____________________________________________________________________________
size_t OsmIdSet::cwrite(int f, const void* buf, size_t n) const {
  ssize_t w = write(f, buf, n);
  if (w < 0) {
    throw std::runtime_error("Could not write to tmp file.\n");
  }

  return w;
}

// _____________________________________________________________________________
size_t OsmIdSet::cread(int f, void* buf, size_t n) const {
  ssize_t w = read(f, buf, n);
  if (w < 0) {
    throw std::runtime_error("Could not read from tmp file.\n");
  }

  return w;
}

// _____________________________________________________________________________
uint32_t OsmIdSet::knuth(uint32_t in) const {
  const uint32_t a = 2654435769;
  return (in * a) >> 2;
}

// _____________________________________________________________________________
uint32_t OsmIdSet::jenkins(uint32_t in) const {
  in = (in + 0x7ed55d16) + (in << 12);
  in = (in ^ 0xc761c23c) ^ (in >> 19);
  in = (in + 0x165667b1) + (in << 5);
  in = (in + 0xd3a2646c) ^ (in << 9);
  in = (in + 0xfd7046c5) + (in << 3);
  in = (in ^ 0xb55a4f09) ^ (in >> 16);
  return in >> 2;
}

// _____________________________________________________________________________
int OsmIdSet::openTmpFile() const {
  const std::string& fname = util::getTmpFName("<tmp>", ".pfaedle-tmp", "");
  int file = open(fname.c_str(), O_RDWR | O_CREAT, 0666);

  // immediately unlink
  unlink(fname.c_str());

  if (file < 0) {
    std::cerr << "Could not open temporary file " << fname << std::endl;
    exit(1);
  }

#ifdef __unix__
  posix_fadvise(file, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif
  return file;
}

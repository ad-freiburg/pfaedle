// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_GTFS_SHAPECONTAINER_H_
#define PFAEDLE_GTFS_SHAPECONTAINER_H_

#include <unistd.h>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <string>
#include "ad/cppgtfs/gtfs/Shape.h"
#include "ad/cppgtfs/gtfs/flat/Shape.h"
#include "pfaedle/Def.h"
#include "util/Misc.h"

namespace pfaedle {
namespace gtfs {

struct Shape {
  explicit Shape(const std::string& id) : id(id) {}
  typedef std::string Ref;
  static std::string getId(Ref r) { return r; }

  template <typename T>
  bool addPoint(T p) {
    UNUSED(p);
    return true;
  }

  const std::string& getId() const { return id; }

  std::string id;
};

template <typename T>
class ShapeContainer {
 public:
  ShapeContainer();
  ~ShapeContainer();
  T* add(const T& obj);
  bool remove(const std::string& id);
  const T* get(const std::string& id) const;
  T* get(const std::string& id);
  const std::string getRef(const std::string& id) const;
  std::string getRef(const std::string& id);
  size_t size() const;
  void finalize() {}
  bool has(const std::string& id) const;

  std::string add(const ad::cppgtfs::gtfs::Shape& s);
  void open();
  bool nextStoragePt(ad::cppgtfs::gtfs::flat::ShapePoint* ret);

 private:
  std::set<std::string> _ids;
  std::fstream _storage;
  size_t _ptr;
  size_t _max;
  std::string _curId;
  std::stringstream _writeBuffer;
};

#include "ShapeContainer.tpp"

}  // namespace gtfs
}  // namespace pfaedle

#endif  // PFAEDLE_GTFS_SHAPECONTAINER_H_

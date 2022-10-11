// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <iostream>
#include <string>

// ____________________________________________________________________________
template <typename T>
ShapeContainer<T>::ShapeContainer() : _lastBuff(0) {
  std::string f = util::getTmpFName("<tmp>", ".pfaedle-tmp", "");
  _storage.open(f, std::fstream::in | std::fstream::out | std::fstream::trunc);

  // immediately unlink
  unlink(f.c_str());

  if (!_storage.good()) {
    std::cerr << "Could not open temporary file " << f << std::endl;
    exit(1);
  }
}

// ____________________________________________________________________________
template <typename T>
ShapeContainer<T>::~ShapeContainer() {
  _storage << _writeBuffer.rdbuf();
  _storage.flush();
  _storage.close();
}

// ____________________________________________________________________________
template <typename T>
T* ShapeContainer<T>::add(const T& ent) {
  _ids.insert(ent.getId());
  return reinterpret_cast<T*>(1);
}

// ____________________________________________________________________________
template <typename T>
bool ShapeContainer<T>::remove(const std::string& id) {
  _ids.erase(id);
  return true;
}

// ____________________________________________________________________________
template <typename T>
T* ShapeContainer<T>::get(const std::string& id) {
  UNUSED(id);
  return reinterpret_cast<T*>(0);
}

// ____________________________________________________________________________
template <typename T>
const T* ShapeContainer<T>::get(const std::string& id) const {
  UNUSED(id);
  return reinterpret_cast<T*>(0);
}

// ____________________________________________________________________________
template <typename T>
bool ShapeContainer<T>::has(const std::string& id) const {
  return _ids.count(id);
}

// ____________________________________________________________________________
template <typename T>
size_t ShapeContainer<T>::size() const {
  return _ids.size();
}

// ____________________________________________________________________________
template <typename T>
std::string ShapeContainer<T>::add(const ad::cppgtfs::gtfs::Shape& s) {
  if (has(s.getId())) return s.getId();
  size_t size = s.getPoints().size();
  _ids.insert(s.getId());

  _writeBuffer << s.getId();
  _writeBuffer.put(' ');
  _writeBuffer.write(reinterpret_cast<const char*>(&size), sizeof(size));

  for (const auto& p : s.getPoints()) {
    _writeBuffer.write(reinterpret_cast<const char*>(&p.lat), sizeof(p.lat));
    _writeBuffer.write(reinterpret_cast<const char*>(&p.lng), sizeof(p.lng));
    _writeBuffer.write(reinterpret_cast<const char*>(&p.travelDist),
                       sizeof(p.travelDist));
  }

  if (_writeBuffer.tellp() - _lastBuff > 1000 * 5000) {
    _lastBuff = _writeBuffer.tellp();
    _storage << _writeBuffer.rdbuf();
    _writeBuffer.clear();
  }

  return s.getId();
}

// ____________________________________________________________________________
template <typename T>
void ShapeContainer<T>::open() {
  _storage << _writeBuffer.rdbuf();
  _storage.flush();
  _writeBuffer.clear();

  _ptr = 0;
  _max = 0;
  _storage.clear();
  _storage.seekg(0, std::ios::beg);
}

// ____________________________________________________________________________
template <typename T>
bool ShapeContainer<T>::nextStoragePt(
    ad::cppgtfs::gtfs::flat::ShapePoint* ret) {
  while (_storage.good() && !_storage.fail()) {
    if (!_ptr) {
      _storage >> _curId;
      _storage.ignore();

      _storage.read(reinterpret_cast<char*>(&_max), sizeof(_max));
    }

    if (!_storage.good() || _storage.fail()) return false;

    _storage.read(reinterpret_cast<char*>(&ret->lat), sizeof(ret->lat));
    _storage.read(reinterpret_cast<char*>(&ret->lng), sizeof(ret->lng));
    _storage.read(reinterpret_cast<char*>(&ret->travelDist),
                  sizeof(ret->travelDist));

    ret->seq = _ptr + 1;
    ret->id = _curId;

    if (_ptr + 1 == _max)
      _ptr = 0;
    else
      _ptr++;

    if (!_storage.good() || _storage.fail()) return false;

    if (has(ret->id)) return true;
  }

  return false;
}

// ____________________________________________________________________________
template <typename T>
const std::string ShapeContainer<T>::getRef(const std::string& id) const {
  if (!has(id)) return "";
  return id;
}

// ____________________________________________________________________________
template <typename T>
std::string ShapeContainer<T>::getRef(const std::string& id) {
  if (!has(id)) return "";
  return id;
}

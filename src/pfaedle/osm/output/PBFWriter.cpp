// Copyright 2025, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <arpa/inet.h>
#include <fcntl.h>
#include <unistd.h>

#include "pfaedle/_config.h"
#include "pfaedle/osm/output/PBFWriter.h"
#ifdef ZLIB_FOUND
#include <zlib.h>
#endif

using pfaedle::osm::output::PBFWriter;
using pfaedle::osm::source::OsmSource;
using pfaedle::osm::source::PBFSource;

static const size_t BUFFER_S = 33 * 1024 * 1024;
static const double RES = .000000001;

// _____________________________________________________________________________
PBFWriter::PBFWriter(const std::string& path,
                     const util::geo::Box<double>& latLngBox,
                     const OsmSource* source)
    : _source(source), _path(path) {
  _file = open(_path.c_str(), O_WRONLY | O_TRUNC | O_CREAT, 0666);
  if (_file < 0) throw std::runtime_error(_path + ": " + "could not open file");

  _blockBuffer = new unsigned char[BUFFER_S + 1];
  _compressBuffer = new unsigned char[BUFFER_S + 1];
  _writeBuffer = new unsigned char[BUFFER_S + 1];

  PBFSource::OSMHeader header;
  header.bbox = latLngBox;
  header.writingProgram = std::string("pfaedle/") + VERSION_FULL;

  writeBlob("OSMHeader", writeOSMHeader(header));
}

// _____________________________________________________________________________
PBFWriter::~PBFWriter() {
  writeDenseNodes();
  writeWays();
  writeRels();
  delete[] _blockBuffer;
  delete[] _compressBuffer;
  delete[] _writeBuffer;

  close(_file);
}

// _____________________________________________________________________________
void PBFWriter::writeBlob(const std::string& type,
                          const PBFSource::Blob& blob) {
  size_t compressedLen = BUFFER_S;
#ifdef ZLIB_FOUND
  auto status = compress(_compressBuffer, &compressedLen,
                         reinterpret_cast<const unsigned char*>(blob.content),
                         blob.datasize);
  if (status != Z_OK)
    throw std::runtime_error(_path + ": " + "could not compress blob");
#else
  compressedLen = blob.datasize;
#endif

  // first 4 bytes will be BlobHeader length
  auto c = _writeBuffer + 4;

  // type
  writeTypeAndId({PBFSource::VarType::S, 1}, c);
  writeString(type, c);

  // datasize
  writeTypeAndId({PBFSource::VarType::V, 3}, c);
#ifdef ZLIB_FOUND
  writeVarUInt(compressedLen + 2 + varUIntNumBytes(compressedLen) +
                   varUIntNumBytes(blob.datasize),
               c);
#else
  writeVarUInt(compressedLen + 1 + varUIntNumBytes(compressedLen), c);
#endif

  // write blob header length
  uint32_t lenNBO = htonl(c - _writeBuffer - 4);
  memcpy(_writeBuffer, &lenNBO, 4);

  // write blob
  size_t written;

#ifdef ZLIB_FOUND
  // raw uncompressed size
  writeTypeAndId({PBFSource::VarType::V, 2}, c);
  writeVarUInt(blob.datasize, c);

  // zlib compressed data
  writeTypeAndId({PBFSource::VarType::S, 3}, c);
  writeVarUInt(compressedLen, c);

  // flush write buffer, which now contains BlobHeader and Blob w/o data
  written = util::writeAll(_file, _writeBuffer, c - _writeBuffer);
  if (written != static_cast<size_t>(c - _writeBuffer)) {
    throw std::runtime_error(_path + ": " + "could not write entire buffer");
  }

  // flush compressedBuffer, which contains the zlib data, compressLen
  // is guaranteed to be within _compressBuffer's bounds
  written = util::writeAll(_file, _compressBuffer, compressedLen);
  if (written != compressedLen) {
    throw std::runtime_error(_path + ": " + "could not write entire buffer");
  }
#else
  // if zlib not available, write raw
  writeTypeAndId({PBFSource::VarType::S, 1}, c);
  writeVarUInt(compressedLen, c);

  // flush write buffer, which now contains BlobHeader and Blob w/o data
  written = util::writeAll(_file, _writeBuffer, c - _writeBuffer);
  if (written != static_cast<size_t>(c - _writeBuffer)) {
    throw std::runtime_error(_path + ": " + "could not write entire buffer");
  }

  // write raw data
  written = util::writeAll(_file,
                           reinterpret_cast<const unsigned char*>(blob.content),
                           compressedLen);
  if (written != compressedLen) {
    throw std::runtime_error(_path + ": " + "could not write entire buffer");
  }
#endif
}

// _____________________________________________________________________________
void PBFWriter::writeDenseNodes() {
  if (_denseNdIds.size() == 0) return;

  PBFSource::Blob ret;
  ret.content = reinterpret_cast<const char*>(_blockBuffer);
  auto c = _blockBuffer;

  // string table
  writeTypeAndId({PBFSource::VarType::S, 1}, c);
  if (_strings.size() == 0) {
    checkBlockBufferBounds(c + 1);
    writeVarUInt(0, c);
  } else {
    checkBlockBufferBounds(c + varUIntNumBytes(_stringTable.size() + 2) + 1 +
                           1 + _stringTable.size());
    writeVarUInt(_stringTable.size() + 2, c);
    writeTypeAndId({PBFSource::VarType::S, 1}, c);
    writeVarUInt(0, c);
    writeBuf(_stringTable.data(), _stringTable.size(), c);
  }

  // primitive group
  checkBlockBufferBounds(c + 1);
  writeTypeAndId({PBFSource::VarType::S, 2}, c);
  size_t groupSize = 0;
  groupSize += _denseNdIds.size() + varUIntNumBytes(_denseNdIds.size()) + 1;
  groupSize += _denseNdLats.size() + varUIntNumBytes(_denseNdLats.size()) + 1;
  groupSize += _denseNdLngs.size() + varUIntNumBytes(_denseNdLngs.size()) + 1;
  groupSize += _denseNdAttrs.size() + varUIntNumBytes(_denseNdAttrs.size()) + 1;
  size_t s = groupSize + 1 + varUIntNumBytes(groupSize);
  checkBlockBufferBounds(c + varUIntNumBytes(s));
  writeVarUInt(s, c);

  // dense nodes
  checkBlockBufferBounds(c + 1 + varUIntNumBytes(groupSize));
  writeTypeAndId({PBFSource::VarType::S, 2}, c);
  writeVarUInt(groupSize, c);

  // ids
  checkBlockBufferBounds(c + 1 + varUIntNumBytes(_denseNdIds.size()) +
                         _denseNdIds.size());
  writeTypeAndId({PBFSource::VarType::S, 1}, c);
  writeVarUInt(_denseNdIds.size(), c);
  writeBuf(_denseNdIds.data(), _denseNdIds.size(), c);

  // lats
  checkBlockBufferBounds(c + 1 + varUIntNumBytes(_denseNdLats.size()) +
                         _denseNdLats.size());
  writeTypeAndId({PBFSource::VarType::S, 8}, c);
  writeVarUInt(_denseNdLats.size(), c);
  writeBuf(_denseNdLats.data(), _denseNdLats.size(), c);

  // lngs
  checkBlockBufferBounds(c + 1 + varUIntNumBytes(_denseNdLngs.size()) +
                         _denseNdLngs.size());
  writeTypeAndId({PBFSource::VarType::S, 9}, c);
  writeVarUInt(_denseNdLngs.size(), c);
  writeBuf(_denseNdLngs.data(), _denseNdLngs.size(), c);

  // attrs
  checkBlockBufferBounds(c + 1 + varUIntNumBytes(_denseNdAttrs.size()) +
                         _denseNdAttrs.size());
  writeTypeAndId({PBFSource::VarType::S, 10}, c);
  writeVarUInt(_denseNdAttrs.size(), c);
  writeBuf(_denseNdAttrs.data(), _denseNdAttrs.size(), c);

  ret.datasize = c - _blockBuffer;

  writeBlob("OSMData", ret);

  _denseNdIds.resize(0);
  _denseNdLats.resize(0);
  _denseNdLngs.resize(0);
  _denseNdAttrs.resize(0);
  _stringTable.resize(0);
  _strings.clear();

  _lastNid = 0;
  _lastLat = 0;
  _lastLng = 0;
}

// _____________________________________________________________________________
void PBFWriter::writeWays() {
  if (_ways.size() == 0) return;

  PBFSource::Blob ret;
  ret.content = reinterpret_cast<const char*>(_blockBuffer);
  auto c = _blockBuffer;

  // string table
  checkBlockBufferBounds(c + 1);
  writeTypeAndId({PBFSource::VarType::S, 1}, c);
  if (_strings.size() == 0) {
    checkBlockBufferBounds(c + 1);
    writeVarUInt(0, c);
  } else {
    checkBlockBufferBounds(c + varUIntNumBytes(_stringTable.size() + 2) + 1 +
                           1 + _stringTable.size());
    writeVarUInt(_stringTable.size() + 2, c);
    writeTypeAndId({PBFSource::VarType::S, 1}, c);
    writeVarUInt(0, c);
    writeBuf(_stringTable.data(), _stringTable.size(), c);
  }

  // primitive group
  checkBlockBufferBounds(c + 1 + varUIntNumBytes(_ways.size()) + _ways.size());
  writeTypeAndId({PBFSource::VarType::S, 2}, c);
  writeVarUInt(_ways.size(), c);

  // ways
  writeBuf(_ways.data(), _ways.size(), c);

  ret.datasize = c - _blockBuffer;

  writeBlob("OSMData", ret);

  _stringTable.resize(0);
  _strings.clear();
  _ways.resize(0);
}

// _____________________________________________________________________________
void PBFWriter::writeRels() {
  if (_rels.size() == 0) return;

  PBFSource::Blob ret;
  ret.content = reinterpret_cast<const char*>(_blockBuffer);
  auto c = _blockBuffer;

  // string table
  checkBlockBufferBounds(c + 1);
  writeTypeAndId({PBFSource::VarType::S, 1}, c);
  if (_strings.size() == 0) {
    checkBlockBufferBounds(c + 1);
    writeVarUInt(0, c);
  } else {
    checkBlockBufferBounds(c + varUIntNumBytes(_stringTable.size() + 2) + 1 +
                           1 + _stringTable.size());
    writeVarUInt(_stringTable.size() + 2, c);
    writeTypeAndId({PBFSource::VarType::S, 1}, c);
    writeVarUInt(0, c);
    writeBuf(_stringTable.data(), _stringTable.size(), c);
  }

  // primitive group
  checkBlockBufferBounds(c + 1 + varUIntNumBytes(_rels.size()) + _rels.size());
  writeTypeAndId({PBFSource::VarType::S, 2}, c);
  writeVarUInt(_rels.size(), c);

  // rels
  writeBuf(_rels.data(), _rels.size(), c);

  ret.datasize = c - _blockBuffer;

  writeBlob("OSMData", ret);

  _stringTable.resize(0);
  _strings.clear();
  _rels.resize(0);
}

// _____________________________________________________________________________
const PBFSource::Blob PBFWriter::writeOSMHeader(
    const PBFSource::OSMHeader& header) {
  PBFSource::Blob ret;
  ret.content = reinterpret_cast<const char*>(_blockBuffer);
  auto c = _blockBuffer;

  // header bbox
  writeTypeAndId({PBFSource::VarType::S, 1}, c);
  size_t s = varIntNumBytes(header.bbox.getLowerLeft().getX() / RES) +
             varIntNumBytes(header.bbox.getLowerLeft().getY() / RES) +
             varIntNumBytes(header.bbox.getUpperRight().getX() / RES) +
             varIntNumBytes(header.bbox.getUpperRight().getY() / RES) + 4;
  writeVarUInt(s, c);
  writeTypeAndId({PBFSource::VarType::V, 1}, c);
  writeVarInt(header.bbox.getLowerLeft().getX() / RES, c);
  writeTypeAndId({PBFSource::VarType::V, 4}, c);
  writeVarInt(header.bbox.getLowerLeft().getY() / RES, c);
  writeTypeAndId({PBFSource::VarType::V, 2}, c);
  writeVarInt(header.bbox.getUpperRight().getX() / RES, c);
  writeTypeAndId({PBFSource::VarType::V, 3}, c);
  writeVarInt(header.bbox.getUpperRight().getY() / RES, c);

  // writing program
  writeTypeAndId({PBFSource::VarType::S, 16}, c);
  writeString(header.writingProgram, c);

  // features
  writeTypeAndId({PBFSource::VarType::S, 4}, c);
  writeString("OsmSchema-V0.6", c);
  writeTypeAndId({PBFSource::VarType::S, 4}, c);
  writeString("DenseNodes", c);

  ret.datasize = c - _blockBuffer;

  return ret;
}

// _____________________________________________________________________________
size_t PBFWriter::varUIntNumBytes(uint64_t val) const {
  if (val == 0) return 1;
  // 7 bit payload per byte
  return ceil((floor(log2(val)) + 1) / 7);
}

// _____________________________________________________________________________
size_t PBFWriter::varIntNumBytes(int64_t val) const {
  // 7 bit payload per byte
  return varUIntNumBytes((static_cast<uint64_t>(val) << 1) ^
                         static_cast<uint64_t>(val >> 63));
}

// _____________________________________________________________________________
void PBFWriter::writeVarInt(int64_t val, unsigned char*& c) const {
  writeVarUInt(
      (static_cast<uint64_t>(val) << 1) ^ static_cast<uint64_t>(val >> 63), c);
}

// _____________________________________________________________________________
void PBFWriter::writeVarUInt(uint64_t val, unsigned char*& c) const {
  while (val >= 128) {  // as long as we have a value that doesnt fit in 7 bits
    // take one byte, ignore highest bit and set to 1 (continue)
    *c = static_cast<uint8_t>(val) | 128;
    c++;
    // shift by the written 7 bits
    val = val >> 7;
  }

  // remaining last byte, with unset continuation bit
  *c = static_cast<uint8_t>(val);
  c++;
}

// _____________________________________________________________________________
void PBFWriter::writeString(const std::string& str, unsigned char*& c) const {
  writeVarUInt(str.size(), c);
  memcpy(c, str.c_str(), str.size());
  c += str.size();
}

// _____________________________________________________________________________
void PBFWriter::writeBuf(const unsigned char* src, size_t s,
                         unsigned char*& c) const {
  memcpy(c, src, s);
  c += s;
}

// _____________________________________________________________________________
void PBFWriter::writeTypeAndId(
    const std::pair<PBFSource::VarType, uint8_t>& typeId,
    unsigned char*& c) const {
  uint64_t byte = 0;
  // set type
  byte |= typeId.first;
  // set id
  byte |= (static_cast<uint64_t>(typeId.second) << 3);

  writeVarUInt(byte, c);
}

// _____________________________________________________________________________
void PBFWriter::writeWay(const OsmWay& w) {
  // make sure dense nodes are flushed
  writeDenseNodes();

  _curWay.resize(0);
  _curKeys.resize(0);
  _curVals.resize(0);
  _curRefs.resize(0);

  int64_t lastRef = 0;

  for (const auto& kv : w.attrs) {
    auto kStrId = getStringId(_source->decode(kv.first));
    auto vStrId = getStringId(_source->decode(kv.second));

    auto c = need(_curKeys, varUIntNumBytes(kStrId));
    writeVarUInt(kStrId, c);

    c = need(_curVals, varUIntNumBytes(vStrId));
    writeVarUInt(vStrId, c);
  }

  for (auto nid : w.nodes) {
    auto dNid = nid - lastRef;
    lastRef = nid;
    auto c = need(_curRefs, varIntNumBytes(dNid));
    writeVarInt(dNid, c);
  }

  auto c =
      need(_curWay, varUIntNumBytes(w.id) + 1 +
                        varUIntNumBytes(_curKeys.size()) + 1 + _curKeys.size() +
                        varUIntNumBytes(_curVals.size()) + 1 + _curVals.size() +
                        varUIntNumBytes(_curRefs.size()) + 1 + _curRefs.size());
  writeTypeAndId({PBFSource::VarType::V, 1}, c);
  writeVarUInt(w.id, c);

  writeTypeAndId({PBFSource::VarType::S, 2}, c);
  writeVarUInt(_curKeys.size(), c);
  writeBuf(_curKeys.data(), _curKeys.size(), c);

  writeTypeAndId({PBFSource::VarType::S, 3}, c);
  writeVarUInt(_curVals.size(), c);
  writeBuf(_curVals.data(), _curVals.size(), c);

  writeTypeAndId({PBFSource::VarType::S, 8}, c);
  writeVarUInt(_curRefs.size(), c);
  writeBuf(_curRefs.data(), _curRefs.size(), c);

  // write way to ways
  c = need(_ways, 1 + varUIntNumBytes(_curWay.size()) + _curWay.size());
  writeTypeAndId({PBFSource::VarType::S, 3}, c);
  writeVarUInt(_curWay.size(), c);
  writeBuf(_curWay.data(), _curWay.size(), c);

  // check size, write ways if buffer is more than 1/3 full
  if (_stringTable.size() + _ways.size() > BUFFER_S / 3) writeWays();
}

// _____________________________________________________________________________
void PBFWriter::writeNode(const OsmNode& nd) {
  size_t dNid = nd.id - _lastNid;
  auto c = need(_denseNdIds, varIntNumBytes(dNid));
  writeVarInt(dNid, c);
  _lastNid = nd.id;

  int64_t dLat = ((nd.lat / RES) / 100) - _lastLat;
  c = need(_denseNdLats, varIntNumBytes(dLat));
  writeVarInt(dLat, c);
  _lastLat = dLat + _lastLat;

  int64_t dLng = (nd.lng / RES) / 100 - _lastLng;
  c = need(_denseNdLngs, varIntNumBytes(dLng));
  writeVarInt(dLng, c);
  _lastLng = dLng + _lastLng;

  for (const auto& kv : nd.attrs) {
    size_t keyStrId = getStringId(_source->decode(kv.first));
    size_t valStrId = getStringId(_source->decode(kv.second));
    c = need(_denseNdAttrs,
             varUIntNumBytes(keyStrId) + varUIntNumBytes(valStrId));
    writeVarUInt(keyStrId, c);
    writeVarUInt(valStrId, c);
  }

  _denseNdAttrs.push_back(0);
  c = _denseNdAttrs.data() + _denseNdAttrs.size() - 1;
  writeVarUInt(0, c);

  // check size, write nodes if buffer is more than 1/3 full
  if (_stringTable.size() + _denseNdIds.size() + _denseNdLats.size() +
          _denseNdLngs.size() + _denseNdAttrs.size() >
      BUFFER_S / 3) {
    writeDenseNodes();
  }
}

// _____________________________________________________________________________
void PBFWriter::writeRel(const OsmRel& rel, const OsmIdList& nodes,
                         const OsmIdList& ways,
                         std::vector<const char*> nodeRoles,
                         std::vector<const char*> wayRoles) {
  // make sure dense nodes are flushed
  writeDenseNodes();
  // make sure ways are flushed
  writeWays();

  _curRel.resize(0);
  _curKeys.resize(0);
  _curVals.resize(0);
  _curRefs.resize(0);
  _curRoles.resize(0);
  _curTypes.resize(0);

  int64_t lastRef = 0;

  for (const auto& kv : rel.attrs) {
    auto kStrId = getStringId(_source->decode(kv.first));
    auto vStrId = getStringId(_source->decode(kv.second));

    auto c = need(_curKeys, varUIntNumBytes(kStrId));
    writeVarUInt(kStrId, c);

    c = need(_curVals, varUIntNumBytes(vStrId));
    writeVarUInt(vStrId, c);
  }

  size_t i = 0;
  for (auto nid : nodes) {
    auto dNid = nid - lastRef;
    lastRef = nid;
    auto c = need(_curRefs, varIntNumBytes(dNid));
    writeVarInt(dNid, c);

    auto strId = getStringId(_source->decode(nodeRoles[i]));
    c = need(_curRoles, varUIntNumBytes(strId));
    writeVarUInt(strId, c);

    c = need(_curTypes, 1);
    writeVarUInt(0, c);

    i++;
  }

  i = 0;
  for (auto wid : ways) {
    auto dWid = wid - lastRef;
    lastRef = wid;
    auto c = need(_curRefs, varIntNumBytes(dWid));
    writeVarInt(dWid, c);

    auto strId = getStringId(_source->decode(wayRoles[i]));
    c = need(_curRoles, varUIntNumBytes(strId));
    writeVarUInt(strId, c);

    c = need(_curTypes, 1);
    writeVarUInt(1, c);

    i++;
  }

  auto c = need(_curRel,
                varUIntNumBytes(rel.id) + 1 + varUIntNumBytes(_curKeys.size()) +
                    1 + _curKeys.size() + varUIntNumBytes(_curVals.size()) + 1 +
                    _curVals.size() + varUIntNumBytes(_curRefs.size()) + 1 +
                    _curRefs.size() + varUIntNumBytes(_curRoles.size()) + 1 +
                    _curRoles.size() + varUIntNumBytes(_curTypes.size()) + 1 +
                    _curTypes.size());
  writeTypeAndId({PBFSource::VarType::V, 1}, c);
  writeVarUInt(rel.id, c);

  writeTypeAndId({PBFSource::VarType::S, 2}, c);
  writeVarUInt(_curKeys.size(), c);
  writeBuf(_curKeys.data(), _curKeys.size(), c);

  writeTypeAndId({PBFSource::VarType::S, 3}, c);
  writeVarUInt(_curVals.size(), c);
  writeBuf(_curVals.data(), _curVals.size(), c);

  writeTypeAndId({PBFSource::VarType::S, 9}, c);
  writeVarUInt(_curRefs.size(), c);
  writeBuf(_curRefs.data(), _curRefs.size(), c);

  writeTypeAndId({PBFSource::VarType::S, 8}, c);
  writeVarUInt(_curRoles.size(), c);
  writeBuf(_curRoles.data(), _curRoles.size(), c);

  writeTypeAndId({PBFSource::VarType::S, 10}, c);
  writeVarUInt(_curTypes.size(), c);
  writeBuf(_curTypes.data(), _curTypes.size(), c);

  // write rel to rels
  c = need(_rels, 1 + varUIntNumBytes(_curRel.size()) + _curRel.size());
  writeTypeAndId({PBFSource::VarType::S, 4}, c);
  writeVarUInt(_curRel.size(), c);
  writeBuf(_curRel.data(), _curRel.size(), c);

  // check size, write rels if buffer is more than 1/3 full
  if (_stringTable.size() + _rels.size() > BUFFER_S / 3) writeRels();
}

// _____________________________________________________________________________
unsigned char* PBFWriter::need(std::vector<unsigned char>& v, size_t n) const {
  v.resize(v.size() + n);
  return v.data() + v.size() - n;
}

// _____________________________________________________________________________
size_t PBFWriter::getStringId(const std::string& str) {
  auto kStr = _strings.find(str);
  if (kStr == _strings.end()) {
    auto c = need(_stringTable, varUIntNumBytes(str.size()) + str.size() + 1);
    writeTypeAndId({PBFSource::VarType::S, 1}, c);
    writeString(str, c);
    kStr = _strings
               .insert(std::pair<std::string, size_t>{str, _strings.size() + 1})
               .first;
  }

  return kStr->second;
}

// _____________________________________________________________________________
void PBFWriter::checkBlockBufferBounds(const unsigned char* c) const {
  if (c > _blockBuffer + BUFFER_S) {
    throw std::runtime_error(_path + ": " + "write out of buffer size");
  }
}

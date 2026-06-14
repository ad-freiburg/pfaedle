// Copyright 2024, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <arpa/inet.h>
#include <fcntl.h>
#include <unistd.h>

#include "pfaedle/osm/source/PBFSource.h"
#include "util/Misc.h"
#include "util/protobuf/Protobuf.h"
#ifdef ZLIB_FOUND
#include <zlib.h>
#endif

using pfaedle::osm::source::OsmSourceAttr;
using pfaedle::osm::source::OsmSourceNode;
using pfaedle::osm::source::OsmSourceRelation;
using pfaedle::osm::source::OsmSourceRelationMember;
using pfaedle::osm::source::OsmSourceWay;
using pfaedle::osm::source::PBFSource;
using util::readAll;
using namespace util::protobuf;

static const size_t BUFFER_S = 33 * 1024 * 1024;
static const double RES = .000000001;

// _____________________________________________________________________________
PBFSource::PBFSource(const std::string& path) : _path(path) {
  _file = open(_path.c_str(), O_RDONLY);
  if (_file < 0) throw std::runtime_error(_path + ": " + "could not open file");

  // single block of buffer memory, for easier bounds checking
  _buf = new unsigned char[2 * BUFFER_S];
  _blockbuf = _buf + BUFFER_S;

#ifdef __unix__
  posix_fadvise(_file, 0, 0, POSIX_FADV_SEQUENTIAL);
#endif

  // used for initialization here
  reset();

  // parse header
  getNextBlock();
}

// _____________________________________________________________________________
bool PBFSource::getNextBlock() {
  // block begins by big-endian 32bit length of BlobHeader
  if (readAll(_file, _buf, 4) != 4) return false;
  _c = _buf;

  uint32_t blobHeaderLength;
  memcpy(&blobHeaderLength, _c, 4);
  blobHeaderLength = ntohl(blobHeaderLength);

  if (blobHeaderLength > BUFFER_S)
    throw std::runtime_error(_path + ": " + "Too large block header length");

  if (readAll(_file, _buf, blobHeaderLength) != blobHeaderLength) return false;
  auto header = parseBlobHeader(blobHeaderLength);

  if (header.datasize > BUFFER_S)
    throw std::runtime_error(_path + ": " + "Too large blob");

  if (readAll(_file, _buf, header.datasize) != header.datasize) return false;
  _c = _buf;

  if (header.type == "OSMHeader") {
    _header = parseOSMHeader(parseBlob(header.datasize));
    return true;
  }

  if (header.type == "OSMData") {
    _curBlock = parseOSMData(parseBlob(header.datasize));
    return true;
  }

  return false;
}

// _____________________________________________________________________________
util::geo::Box<double> PBFSource::parseHeaderBBox(unsigned char*& c) {
  auto len = util::protobuf::parseVarUInt(c);
  auto start = c;
  checkBufferBounds(c + len);

  double lly = -90;
  double llx = -180;
  double ury = 90;
  double urx = 180;

  while ((size_t)(c - start) < len) {
    auto typeId = util::protobuf::nextTypeAndId(c);

    if (typeId.second == 1)
      llx = RES * util::protobuf::parseInt(typeId, c);
    else if (typeId.second == 4)
      lly = RES * util::protobuf::parseInt(typeId, c);
    else if (typeId.second == 2)
      urx = RES * util::protobuf::parseInt(typeId, c);
    else if (typeId.second == 3)
      ury = RES * util::protobuf::parseInt(typeId, c);
    else
      skipType(typeId.first, c);
  }

  return {{llx, lly}, {urx, ury}};
}

// _____________________________________________________________________________
PBFSource::PrimitiveBlock PBFSource::parseOSMData(const Blob& blob) {
  auto c = reinterpret_cast<unsigned char*>(const_cast<char*>(blob.content));
  auto start = c;
  checkBufferBounds(c + blob.datasize);

  PBFSource::PrimitiveBlock block;

  while ((size_t)(c - start) < blob.datasize) {
    auto typeId = util::protobuf::nextTypeAndId(c);

    if (typeId.second == 1) {
      block.stringTable = parseStringTable(c);
    } else if (typeId.second == 2) {
      block.primitiveGroups.push(c);
      c += util::protobuf::parseVarUInt(c);
    } else if (typeId.second == 17) {
      block.granularity = util::protobuf::parseUInt(typeId, c);
    } else if (typeId.second == 19) {
      block.latOffset = util::protobuf::parseUInt(typeId, c);
    } else if (typeId.second == 20) {
      block.lonOffset = util::protobuf::parseUInt(typeId, c);
    } else if (typeId.second == 18) {
      block.dateGranularity = util::protobuf::parseUInt(typeId, c);
    } else {
      skipType(typeId.first, c);
    }
  }

  return block;
}

// _____________________________________________________________________________
std::vector<std::string> PBFSource::parseStringTable(unsigned char*& c) {
  size_t len = util::protobuf::parseVarUInt(c);
  auto start = c;
  checkBufferBounds(c + len);

  std::vector<std::string> table;
  table.reserve(len / 10);

  while ((size_t)(c - start) < len) {
    auto typeId = util::protobuf::nextTypeAndId(c);
    if (typeId.second == 1)
      table.push_back(util::protobuf::parseString(c));
    else
      skipType(typeId.first, c);
  }

  return table;
}

// _____________________________________________________________________________
PBFSource::OSMHeader PBFSource::parseOSMHeader(const Blob& blob) {
  auto c = reinterpret_cast<unsigned char*>(const_cast<char*>(blob.content));
  auto start = c;
  checkBufferBounds(c + blob.datasize);

  PBFSource::OSMHeader header;

  header.bbox = {{-180, -90}, {180, 90}};

  while ((size_t)(c - start) < blob.datasize) {
    auto typeId = util::protobuf::nextTypeAndId(c);

    if (typeId.second == 1)
      header.bbox = parseHeaderBBox(c);
    else if (typeId.second == 4)
      header.requiredFeatures.push_back(util::protobuf::parseString(c));
    else if (typeId.second == 5)
      header.optionalFeatures.push_back(util::protobuf::parseString(c));
    else if (typeId.second == 16)
      header.writingProgram = util::protobuf::parseString(c);
    else
      skipType(typeId.first, c);
  }

  return header;
}

// _____________________________________________________________________________
PBFSource::Blob PBFSource::parseBlob(size_t len) {
  auto start = _c;
  checkBufferBounds(_c + len);

  PBFSource::Blob ret;

  while ((size_t)(_c - start) < len) {
    auto typeId = nextTypeAndId();

    if (typeId.second == 1) {
      // raw, no compression
      ret.datasize = parseVarUInt();
      ret.content = reinterpret_cast<const char*>(_c);
      _c += ret.datasize;
    } else if (typeId.second == 4) {
      throw std::runtime_error(_path + ": " + "LZMA compression not supported");
    } else if (typeId.second == 5) {
      throw std::runtime_error(_path + ": " + "BZ2 compression not supported");
    } else if (typeId.second == 6) {
      throw std::runtime_error(_path + ": " + "LZ4 compression not supported");
    } else if (typeId.second == 7) {
      throw std::runtime_error(_path + ": " + "ZSTD compression not supported");
    } else if (typeId.second == 3) {
#ifdef ZLIB_FOUND
      // ZLIB compression
      if (typeId.first != VarType::S) {
        throw std::runtime_error(_path + ": " + "expected byte array value");
      }
      size_t size = parseVarUInt();
      size_t uncompressedSize = BUFFER_S;
      auto status = uncompress(_blockbuf, &uncompressedSize, _c, size);
      if (status != Z_OK) {
        throw std::runtime_error(_path + ": " + "could not uncompress buffer");
      }

      ret.content = reinterpret_cast<const char*>(_blockbuf);
      ret.datasize = uncompressedSize;
      _c += size;
#else
      throw std::runtime_error(_path + ": " +
                               "ZLIB compression not supported, pfaedle was "
                               "compiled without zlib support");
#endif
    } else if (typeId.second == 2) {
      ret.datasize = parseVarUInt();
    } else {
      skipType(typeId.first);
    }
  }

  return ret;
}

// _____________________________________________________________________________
PBFSource::BlobHeader PBFSource::parseBlobHeader(size_t len) {
  auto start = _c;

  BlobHeader ret;

  while ((size_t)(_c - start) < len) {
    auto typeId = nextTypeAndId();

    if (typeId.second == 1) {
      ret.type = parseString();
    } else if (typeId.second == 3) {
      ret.datasize = parseVarUInt();
    } else {
      skipType(typeId.first);
    }
  }

  return ret;
}

// _____________________________________________________________________________
uint32_t PBFSource::parseFixedUInt32() {
  return util::protobuf::parseFixedUInt32(_c);
}

// _____________________________________________________________________________
uint64_t PBFSource::parseFixedUInt64() {
  return util::protobuf::parseFixedUInt64(_c);
}

// _____________________________________________________________________________
int32_t PBFSource::parseFixedInt32() {
  return util::protobuf::parseFixedInt32(_c);
}

// _____________________________________________________________________________
int64_t PBFSource::parseFixedInt64() {
  return util::protobuf::parseFixedInt64(_c);
}

// _____________________________________________________________________________
std::string PBFSource::parseString() { return util::protobuf::parseString(_c); }

// _____________________________________________________________________________
uint64_t PBFSource::parseUInt(std::pair<VarType, uint8_t> typeId) {
  return util::protobuf::parseUInt(typeId, _c);
}

// _____________________________________________________________________________
int64_t PBFSource::parseInt(std::pair<VarType, uint8_t> typeId) {
  return util::protobuf::parseInt(typeId, _c);
}

// _____________________________________________________________________________
int64_t PBFSource::parseVarInt() { return util::protobuf::parseVarInt(_c); }

// _____________________________________________________________________________
uint64_t PBFSource::parseVarUInt() { return util::protobuf::parseVarUInt(_c); }

// _____________________________________________________________________________
std::pair<VarType, uint8_t> PBFSource::nextTypeAndId() {
  return util::protobuf::nextTypeAndId(_c);
}

// _____________________________________________________________________________
void PBFSource::skipType(VarType type) { skipType(type, _c); }

// _____________________________________________________________________________
void PBFSource::skipType(VarType type, unsigned char*& c) {
  if (type == VarType::V)
    util::protobuf::parseVarUInt(c);
  else if (type == VarType::D)
    util::protobuf::parseFixedUInt64(c);
  else if (type == VarType::S)
    util::protobuf::parseString(c);
  else if (type == VarType::I)
    util::protobuf::parseFixedUInt32(c);
  else
    throw std::runtime_error(_path + ": " + "parse error");
}

// _____________________________________________________________________________
PBFSource::~PBFSource() { delete[] _buf; }

// _____________________________________________________________________________
bool PBFSource::cont() { return true; }
// _____________________________________________________________________________
const OsmSourceNode* PBFSource::nextNode() {
  // if we have a list of parsed dense nodes, first return them...
  if (_curBlock.denseNodePtr < _curBlock.curDenseNodes.size()) {
    const auto& a = _curBlock.curDenseNodes[_curBlock.denseNodePtr++];
    _curNode.id = a.id;
    _curNode.lat =
        RES * (_curBlock.latOffset + (a.lat * _curBlock.granularity));
    _curNode.lon =
        RES * (_curBlock.lonOffset + (a.lon * _curBlock.granularity));
    _curAttr = 0;
    _curAttrsPtr = &a.tags;
    return &_curNode;
  }

  while (checkGroup()) {
    auto typeId = util::protobuf::nextTypeAndId(_curBlock.c);
    if (typeId.second == 1) {
      parseNode(_curBlock.c);
      _curAttr = 0;
      _curAttrsPtr = &_curAttrs;
      return &_curNode;
    } else if (typeId.second == 2) {
      _curBlock.curDenseNodes = parseDenseNodes(_curBlock.c);
      _curBlock.denseNodePtr = 0;
      return nextNode();
    } else if (typeId.second == 4) {
      skip(util::protobuf::parseVarUInt(_curBlock.c));
    } else {
      skip(util::protobuf::parseVarUInt(_curBlock.c));
    }
  }

  return 0;
}

// _____________________________________________________________________________
void PBFSource::seekNodes() { reset(); }

// _____________________________________________________________________________
void PBFSource::seekWays() { reset(); }

// _____________________________________________________________________________
void PBFSource::seekRels() { reset(); }

// _____________________________________________________________________________
void PBFSource::reset() {
  _curBlock = {};
  lseek(_file, 0, SEEK_SET);

  // get header block
  getNextBlock();

  checkGroup();
}

// _____________________________________________________________________________
bool PBFSource::checkGroup() {
  while (true) {
    // skip to first non-empty block
    while (_curBlock.primitiveGroups.size() == 0 && getNextBlock()) {
    }

    // could not find primitive group, we are at the end...
    if (_curBlock.primitiveGroups.size() == 0) return false;

    if (_curBlock.c == 0) {
      _curBlock.c = _curBlock.primitiveGroups.front();
      _curBlock.curGroupLen = util::protobuf::parseVarUInt(_curBlock.c);
      _curBlock.primitiveGroups.front() = _curBlock.c;
    }

    if ((size_t)(_curBlock.c - _curBlock.primitiveGroups.front()) >=
        _curBlock.curGroupLen) {
      _curBlock.primitiveGroups.pop();
      _curBlock.c = 0;
      continue;
    }

    return true;
  }
}

// _____________________________________________________________________________
void PBFSource::parseWay(unsigned char*& c) {
  // normal node
  size_t len = util::protobuf::parseVarUInt(c);
  auto start = c;
  checkBufferBounds(c + len);

  size_t curKey = 0;
  size_t curValue = 0;
  size_t curMemberID = 0;

  _curAttrs = {};
  _curWayMembers = {};
  _curWayMember = 0;

  while ((size_t)(c - start) < len) {
    auto typeId = util::protobuf::nextTypeAndId(c);

    if (typeId.second == 1) {  // ID
      _curWay.id = util::protobuf::parseUInt(typeId, c);
    } else if (typeId.second == 8) {  // member IDs
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      while ((size_t)(c - start) < len) {
        auto memberId = util::protobuf::parseVarInt(c);
        if (curMemberID > 0)
          memberId = _curWayMembers[curMemberID - 1] + memberId;
        if (curMemberID >= _curWayMembers.size()) _curWayMembers.push_back(0);
        _curWayMembers[curMemberID++] = memberId;
      }
    } else if (typeId.second == 2) {  // attr keys
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      while ((size_t)(c - start) < len) {
        if (curKey >= _curAttrs.size()) _curAttrs.push_back({0, 0});
        _curAttrs[curKey++].first = util::protobuf::parseVarUInt(c);
      }
    } else if (typeId.second == 3) {  // attr values
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      while ((size_t)(c - start) < len) {
        if (curValue >= _curAttrs.size()) _curAttrs.push_back({0, 0});
        _curAttrs[curValue++].second = util::protobuf::parseVarUInt(c);
      }
    } else {
      c += util::protobuf::parseVarUInt(c);
    }
  }
}

// _____________________________________________________________________________
void PBFSource::parseRelation(unsigned char*& c) {
  // normal node
  size_t len = util::protobuf::parseVarUInt(c);
  auto start = c;
  checkBufferBounds(c + len);

  size_t curKey = 0;
  size_t curValue = 0;
  size_t curRole = 0;
  size_t curMemberID = 0;
  size_t curType = 0;

  _curAttrs = {};
  _curRelMembers = {};
  _curRelMember = 0;

  while ((size_t)(c - start) < len) {
    auto typeId = util::protobuf::nextTypeAndId(c);

    if (typeId.second == 1) {  // ID
      _curRelation.id = util::protobuf::parseVarUInt(c);
    } else if (typeId.second == 8) {  // roles
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      while ((size_t)(c - start) < len) {
        if (curRole >= _curRelMembers.size())
          _curRelMembers.push_back({0, 0, 0});
        _curRelMembers[curRole++].role =
            _curBlock.stringTable[util::protobuf::parseVarUInt(c)].c_str();
      }
    } else if (typeId.second == 9) {  // member IDs
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      while ((size_t)(c - start) < len) {
        auto memberId = util::protobuf::parseVarInt(c);
        if (curMemberID > 0)
          memberId = _curRelMembers[curMemberID - 1].id + memberId;
        if (curMemberID >= _curRelMembers.size())
          _curRelMembers.push_back({0, 0, 0});
        _curRelMembers[curMemberID++].id = memberId;
      }
    } else if (typeId.second == 10) {  // member types
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      while ((size_t)(c - start) < len) {
        if (curType >= _curRelMembers.size())
          _curRelMembers.push_back({0, 0, 0});
        _curRelMembers[curType++].type = util::protobuf::parseVarUInt(c);
      }
    } else if (typeId.second == 2) {  // attr keys
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      while ((size_t)(c - start) < len) {
        if (curKey >= _curAttrs.size()) _curAttrs.push_back({0, 0});
        _curAttrs[curKey++].first = util::protobuf::parseVarUInt(c);
      }
    } else if (typeId.second == 3) {  // attr values
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      while ((size_t)(c - start) < len) {
        if (curValue >= _curAttrs.size()) _curAttrs.push_back({0, 0});
        _curAttrs[curValue++].second = util::protobuf::parseVarUInt(c);
      }
    } else {
      c += util::protobuf::parseVarUInt(c);
    }
  }
}

// _____________________________________________________________________________
void PBFSource::parseNode(unsigned char*& c) {
  // normal node
  size_t len = util::protobuf::parseVarUInt(c);
  auto start = c;
  checkBufferBounds(c + len);

  size_t curKey = 0;
  size_t curValue = 0;

  _curAttrs = {};

  while ((size_t)(c - start) < len) {
    auto typeId = util::protobuf::nextTypeAndId(c);

    if (typeId.second == 1) {  // ID
      _curNode.id = util::protobuf::parseInt(typeId, c);
    } else if (typeId.second == 8) {  // lat
      _curNode.lat =
          RES * (_curBlock.latOffset +
                 ((double)util::protobuf::parseInt(typeId, c) * _curBlock.granularity));
    } else if (typeId.second == 9) {  // lon
      _curNode.lon =
          RES * (_curBlock.lonOffset +
                 ((double)util::protobuf::parseInt(typeId, c) * _curBlock.granularity));
    } else if (typeId.second == 2) {  // attr keys
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      while ((size_t)(c - start) < len) {
        if (curKey >= _curAttrs.size()) _curAttrs.push_back({0, 0});
        _curAttrs[curKey++].first = util::protobuf::parseVarUInt(c);
      }
    } else if (typeId.second == 3) {  // attr values
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      while ((size_t)(c - start) < len) {
        if (curValue >= _curAttrs.size()) _curAttrs.push_back({0, 0});
        _curAttrs[curValue++].second = util::protobuf::parseVarUInt(c);
      }
    } else {
      c += util::protobuf::parseVarUInt(c);
    }
  }
}

// _____________________________________________________________________________
std::vector<PBFSource::Node> PBFSource::parseDenseNodes(unsigned char*& c) {
  std::vector<PBFSource::Node> ret;

  size_t len = util::protobuf::parseVarUInt(c);
  auto start = c;
  checkBufferBounds(c + len);

  size_t curIdPos = 0;
  size_t curLatPos = 0;
  size_t curLonPos = 0;
  size_t curTagPos = 0;

  while ((size_t)(c - start) < len) {
    auto typeId = util::protobuf::nextTypeAndId(c);

    if (typeId.second == 1) {  // IDs
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      // delta encoded
      while ((size_t)(c - start) < len) {
        uint64_t nid = util::protobuf::parseVarInt(c);
        if (curIdPos > 0) nid = ret[curIdPos - 1].id + nid;

        if (curIdPos >= ret.size())
          ret.push_back({0, 0, nid, {}});
        else
          ret[curIdPos].id = nid;

        curIdPos++;
      }
    } else if (typeId.second == 5) {  // denseinfo
      c += util::protobuf::parseVarUInt(c);
    } else if (typeId.second == 8) {  // lat
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      // delta encoded
      while ((size_t)(c - start) < len) {
        double lat = util::protobuf::parseVarInt(c);
        if (curLatPos > 0) lat = ret[curLatPos - 1].lat + lat;

        if (curLatPos >= ret.size())
          ret.push_back({lat, 0, 0, {}});
        else
          ret[curLatPos].lat = lat;

        curLatPos++;
      }
    } else if (typeId.second == 9) {  // lon
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      // delta encoded
      while ((size_t)(c - start) < len) {
        double lon = util::protobuf::parseVarInt(c);
        if (curLonPos > 0) lon = ret[curLonPos - 1].lon + lon;

        if (curLonPos >= ret.size())
          ret.push_back({0, lon, 0, {}});
        else
          ret[curLonPos].lon = lon;

        curLonPos++;
      }
    } else if (typeId.second == 10) {  // attrs
      size_t len = util::protobuf::parseVarUInt(c);
      auto start = c;
      checkBufferBounds(c + len);

      // delta encoded
      while ((size_t)(c - start) < len) {
        auto key = util::protobuf::parseVarUInt(c);
        if (key == 0) {
          curTagPos++;
          continue;
        }

        if (curTagPos >= ret.size()) ret.push_back({0, 0, 0, {}});
        ret[curTagPos].tags.push_back({key, util::protobuf::parseVarUInt(c)});
      }
    }
  }

  return ret;
}

// _____________________________________________________________________________
void PBFSource::skip(size_t n) { _curBlock.c += n; }

// _____________________________________________________________________________
const OsmSourceWay* PBFSource::nextWay() {
  while (checkGroup()) {
    auto typeId = util::protobuf::nextTypeAndId(_curBlock.c);
    if (typeId.second == 3) {
      parseWay(_curBlock.c);

      _curAttr = 0;
      _curAttrsPtr = &_curAttrs;
      return &_curWay;
    } else {
      skip(util::protobuf::parseVarUInt(_curBlock.c));
    }
  }

  return 0;
}

// _____________________________________________________________________________
const OsmSourceRelationMember* PBFSource::nextMember() {
  _curRelMember++;
  if (_curRelMember > _curRelMembers.size()) return 0;
  return &_curRelMembers[_curRelMember - 1];
}

// _____________________________________________________________________________
uint64_t PBFSource::nextMemberNode() {
  _curWayMember++;
  if (_curWayMember > _curWayMembers.size()) return 0;
  return _curWayMembers[_curWayMember - 1];
}

// _____________________________________________________________________________
const OsmSourceRelation* PBFSource::nextRel() {
  while (checkGroup()) {
    auto typeId = util::protobuf::nextTypeAndId(_curBlock.c);
    if (typeId.second == 4) {
      parseRelation(_curBlock.c);

      _curAttr = 0;
      _curAttrsPtr = &_curAttrs;
      return &_curRelation;
    } else {
      skip(util::protobuf::parseVarUInt(_curBlock.c));
    }
  }

  return 0;
}

// _____________________________________________________________________________
const OsmSourceAttr PBFSource::nextAttr() {
  _curAttr++;
  if (_curAttr > _curAttrsPtr->size()) return {0, 0};
  return {_curBlock.stringTable[(*_curAttrsPtr)[_curAttr - 1].first].c_str(),
          _curBlock.stringTable[(*_curAttrsPtr)[_curAttr - 1].second].c_str()};
}

// _____________________________________________________________________________
util::geo::Box<double> PBFSource::getBounds() { return _header.bbox; }

// _____________________________________________________________________________
std::string PBFSource::decode(const char* str) const { return str; }

// _____________________________________________________________________________
std::string PBFSource::decode(const std::string& str) const { return str; }

// _____________________________________________________________________________
void PBFSource::checkBufferBounds(const unsigned char* c) const {
  if (c > _buf + 2 * BUFFER_S) {
    throw std::runtime_error(
        _path + ": " +
        "read out of buffer size, most likely caused by corrupt input file");
  }
}

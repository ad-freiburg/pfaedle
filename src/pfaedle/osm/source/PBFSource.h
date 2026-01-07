// Copyright 2024, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_OSM_SOURCE_PBFSOURCE_H_
#define PFAEDLE_OSM_SOURCE_PBFSOURCE_H_

#include <queue>

#include "pfaedle/osm/source/OsmSource.h"
#include "util/geo/Geo.h"

namespace pfaedle {
namespace osm {
namespace source {

class PBFSource : public OsmSource {
 public:
  enum VarType : uint8_t {
    V = 0,
    D = 1,
    S = 2,
    I = 5,
  };

  struct Node {
    double lat, lon;
    uint64_t id;
    std::vector<std::pair<size_t, size_t>> tags;
  };

  struct OSMHeader {
    util::geo::Box<double> bbox;
    std::vector<std::string> requiredFeatures;
    std::vector<std::string> optionalFeatures;
    std::string writingProgram;
  };

  struct PrimitiveBlock {
    std::vector<std::string> stringTable;
    std::queue<unsigned char*> primitiveGroups;
    uint32_t granularity = 100;
    uint64_t latOffset = 0;
    uint64_t lonOffset = 0;
    uint32_t dateGranularity = 1000;

    unsigned char* c = 0;
    size_t curGroupLen = 0;

    std::vector<Node> curDenseNodes;
    size_t denseNodePtr = 0;
  };

  struct BlobHeader {
    std::string type;
    uint32_t datasize;
  };

  struct Blob {
    const char* content;
    uint32_t datasize;
  };
  PBFSource(const std::string& path);
  virtual ~PBFSource();
  virtual const OsmSourceNode* nextNode();
  virtual const OsmSourceAttr nextAttr();
  virtual const OsmSourceWay* nextWay();
  virtual uint64_t nextMemberNode();
  virtual const OsmSourceRelationMember* nextMember();
  virtual const OsmSourceRelation* nextRel();
  virtual bool cont();

  virtual void seekNodes();
  virtual void seekWays();
  virtual void seekRels();

  virtual util::geo::Box<double> getBounds();

  virtual std::string decode(const char* str) const;
  virtual std::string decode(const std::string& str) const;

 private:
  std::string _path;
  int _file;

  unsigned char* _buf;
  unsigned char* _c;

  unsigned char* _blockbuf;

  PrimitiveBlock _curBlock;
  OSMHeader _header;

  bool getNextBlock();
  OSMHeader parseOSMHeader(const Blob& blob);
  PrimitiveBlock parseOSMData(const Blob& blob);
  BlobHeader parseBlobHeader(size_t len);
  std::vector<std::string> parseStringTable(unsigned char*& c);
  Blob parseBlob(size_t len);
  std::pair<VarType, uint8_t> nextTypeAndId();
  std::pair<VarType, uint8_t> nextTypeAndId(unsigned char*& c);
  std::string parseString();
  std::string parseString(unsigned char*& c);
  uint32_t parseFixedUInt32();
  uint32_t parseFixedUInt32(unsigned char*& c);
  uint64_t parseFixedUInt64();
  uint64_t parseFixedUInt64(unsigned char*& c);
  int32_t parseFixedInt32();
  int32_t parseFixedInt32(unsigned char*& c);
  int64_t parseFixedInt64();
  int64_t parseFixedInt64(unsigned char*& c);
  int64_t parseVarInt();
  int64_t parseVarInt(unsigned char*& c);
  uint64_t parseVarUInt();
  uint64_t parseVarUInt(unsigned char*& c);
  uint64_t parseUInt(std::pair<VarType, uint8_t> typeId);
  int64_t parseInt(std::pair<VarType, uint8_t> typeId);
  uint64_t parseUInt(std::pair<VarType, uint8_t> typeId, unsigned char*& c);
  int64_t parseInt(std::pair<VarType, uint8_t> typeId, unsigned char*& c);
  void skipType(VarType type, unsigned char*& c);
  void skipType(VarType type);
  void skip(size_t n);
  void reset();
  virtual bool checkGroup();
  util::geo::Box<double> parseHeaderBBox(unsigned char*& c);

  std::vector<Node> parseDenseNodes(unsigned char*& c);
  void parseNode(unsigned char*& c);
  void parseRelation(unsigned char*& c);
  void parseWay(unsigned char*& c);

  void checkBufferBounds(const unsigned char* c) const;

  OsmSourceNode _curNode;
  OsmSourceRelation _curRelation;
  OsmSourceWay _curWay;
  size_t _curAttr;
  const std::vector<std::pair<size_t, size_t>>* _curAttrsPtr;
  std::vector<std::pair<size_t, size_t>> _curAttrs;
  std::vector<OsmSourceRelationMember> _curRelMembers;
  size_t _curRelMember;
  std::vector<uint64_t> _curWayMembers;
  size_t _curWayMember;
};

}  // namespace source
}  // namespace osm
}  // namespace pfaedle

#endif

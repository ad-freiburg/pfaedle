// Copyright 2025, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#ifndef PFAEDLE_OSM_WRITER_PBFWRITER_H_
#define PFAEDLE_OSM_WRITER_PBFWRITER_H_

#include "pfaedle/osm/output/OsmWriter.h"
#include "pfaedle/osm/source/OsmSource.h"
#include "pfaedle/osm/source/PBFSource.h"

namespace pfaedle {
namespace osm {
namespace output {

class PBFWriter : public OsmWriter {
 public:
  PBFWriter(const std::string& path, const util::geo::Box<double>& latLngBox,
            const source::OsmSource* source);
  ~PBFWriter();
  virtual void writeNode(const OsmNode& node);
  virtual void writeWay(const OsmWay& node);
  virtual void writeRel(const OsmRel& nd, const OsmIdList& nodes,
                        const OsmIdList& ways,
                        std::vector<const char*> nodeRoles,
                        std::vector<const char*> wayRoles);

 private:
  const source::OsmSource* _source;
  std::string _path;
  int _file;

  unsigned char* _blockBuffer;
  unsigned char* _compressBuffer;
  unsigned char* _writeBuffer;

  // individual buffers
  std::vector<unsigned char> _denseNdIds, _denseNdLats, _denseNdLngs,
      _denseNdAttrs, _ways, _rels, _curWay, _curRel, _curKeys, _curVals,
      _curRefs, _curRoles, _curTypes, _stringTable;

  size_t _lastNid = 0;
  int64_t _lastLat = 0;
  int64_t _lastLng = 0;

  std::map<std::string, size_t> _strings;

  void writeDenseNodes();
  void writeWays();
  void writeRels();
  void writeBlob(const std::string& type, const source::PBFSource::Blob& blob);
  const source::PBFSource::Blob writeOSMHeader(
      const source::PBFSource::OSMHeader& header);
  void writeTypeAndId(
      const std::pair<source::PBFSource::VarType, uint8_t>& typeId,
      unsigned char*& c) const;
  void writeVarUInt(uint64_t val, unsigned char*& c) const;
  size_t varUIntNumBytes(uint64_t val) const;
  void writeVarInt(int64_t val, unsigned char*& c) const;
  size_t varIntNumBytes(int64_t val) const;
  void writeString(const std::string& str, unsigned char*& c) const;
  size_t getStringId(const std::string& str);
  unsigned char* need(std::vector<unsigned char>& v, size_t n) const;
  void writeBuf(const unsigned char* src, size_t s, unsigned char*& c) const;
  void checkBlockBufferBounds(const unsigned char* c) const;
};

}  // namespace output
}  // namespace osm
}  // namespace pfaedle

#endif

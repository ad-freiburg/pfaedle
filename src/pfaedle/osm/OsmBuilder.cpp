// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <osmium/io/any_input.hpp>
#include <osmium/geom/haversine.hpp>
#include <osmium/visitor.hpp>

#include <algorithm>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <stack>
#include <string>
#include <utility>
#include <vector>

#include "pfaedle/Def.h"
#include "pfaedle/osm/BBoxIdx.h"
#include "pfaedle/osm/Osm.h"
#include "pfaedle/osm/OsmBuilder.h"
#include "pfaedle/osm/OsmFilter.h"
#include "pfaedle/osm/Restrictor.h"
#include "pfaedle/trgraph/StatGroup.h"
#include "util/Misc.h"
#include "util/Nullable.h"
#include "util/log/Log.h"
#include "xml/pfxml.h"

using util::geo::webMercMeterDist;
using util::geo::Box;
using util::Nullable;
using pfaedle::trgraph::Normalizer;
using pfaedle::trgraph::Graph;
using pfaedle::trgraph::Node;
using pfaedle::trgraph::NodePL;
using pfaedle::trgraph::Edge;
using pfaedle::trgraph::EdgePL;
using pfaedle::trgraph::TransitEdgeLine;
using pfaedle::trgraph::StatInfo;
using pfaedle::trgraph::StatGroup;
using pfaedle::trgraph::Component;
using pfaedle::osm::OsmBuilder;
using pfaedle::osm::OsmWay;
using pfaedle::osm::OsmRel;
using pfaedle::osm::OsmNode;
using pfaedle::osm::EdgeGrid;
using pfaedle::osm::NodeGrid;
using pfaedle::osm::EqSearch;
using pfaedle::osm::BlockSearch;
using ad::cppgtfs::gtfs::Stop;

class NodeHandler : public osmium::handler::Handler {
        const pfaedle::osm::OsmFilter& _filter;
        const pfaedle::osm::BBoxIdx& _bbox;
        pfaedle::osm::OsmIdSet& _bboxNodes;
        pfaedle::osm::OsmIdSet& _noHupNodes;

    public:
        NodeHandler(const pfaedle::osm::OsmFilter& filter,
                    const pfaedle::osm::BBoxIdx& bbox,
                    pfaedle::osm::OsmIdSet& bboxNodes,
                    pfaedle::osm::OsmIdSet& noHupNodes) :
            _filter(filter),
            _bbox(bbox),
            _bboxNodes(bboxNodes),
            _noHupNodes(noHupNodes)
        {}

        void node(const osmium::Node& node) {
            bool ignored = false;

            for (const auto& tag : node.tags()) {
                if (_filter.nohup(tag.key(), tag.value())) {
                    _noHupNodes.add(node.id());
                    ignored = true;
                    break;
                }
            }

            if (!ignored) {
                Point<double> point{node.location().lon(), node.location().lat()};
                if (_bbox.contains(point)) {
                    _bboxNodes.add(node.id());
                }
            }
        }
};

class RelationHandler: public osmium::handler::Handler {
        const pfaedle::osm::OsmFilter& _filter;
        const pfaedle::osm::BBoxIdx& _bbox;
        const pfaedle::osm::AttrKeySet& _keepAttrs;
        pfaedle::osm::RelLst& _rels;
        pfaedle::osm::RelMap& _nodeRels;
        pfaedle::osm::RelMap& _wayRels;
        pfaedle::osm::Restrictions& _restrictions;

    public:
        RelationHandler(const pfaedle::osm::OsmFilter& filter,
                        const pfaedle::osm::BBoxIdx& bbox,
                        const pfaedle::osm::AttrKeySet& keepAttrs,
                        pfaedle::osm::RelLst& rels,
                        pfaedle::osm::RelMap& nodeRels,
                        pfaedle::osm::RelMap& wayRels,
                        pfaedle::osm::Restrictions& restrictions) :
            _filter(filter),
            _bbox(bbox),
            _keepAttrs(keepAttrs),
            _rels(rels),
            _nodeRels(nodeRels),
            _wayRels(wayRels),
            _restrictions(restrictions)
        {}

        void relation(const osmium::Relation &relation) {
            OsmRel rel;
            rel.id = relation.id();

            if (rel.id == 0) {
                return;
            }

            for (const auto& tag : relation.tags()) {
                if (_keepAttrs.count(tag.key())) {
                    rel.attrs[tag.key()] = tag.value();
                }
            }

            for (const auto& member : relation.members()) {
                using osmium::item_type;
                auto& obj = member.get_object();
                if (member.type() == item_type::node) {
                    rel.nodes.push_back(obj.id());
                    rel.nodeRoles.emplace_back(member.role());
                } else if (member.type() == item_type::way) {
                    rel.ways.push_back(obj.id());
                    rel.wayRoles.emplace_back(member.role());
                }
            }
            for (auto id : rel.nodes) {
                _nodeRels[id].push_back(_rels.rels.size() - 1);
            }

            for (auto id : rel.ways) {
                _wayRels[id].push_back(_rels.rels.size() - 1);
            }

            uint64_t keepFlags = _filter.keep(rel.attrs, pfaedle::osm::OsmFilter::REL);
            uint64_t dropFlags = _filter.drop(rel.attrs, pfaedle::osm::OsmFilter::REL);

            if (rel.id && !rel.attrs.empty() && keepFlags && !dropFlags) {
                rel.keepFlags = keepFlags;
                rel.dropFlags = dropFlags;
            }

            _rels.rels.push_back(rel.attrs);
            if (rel.keepFlags & pfaedle::osm::REL_NO_DOWN) {
                _rels.flat.insert(_rels.rels.size() - 1);
            }

            // more processing
            {
                if (!rel.attrs.count("type")) return;
                if (rel.attrs.find("type")->second != "restriction") return;

                bool pos = _filter.posRestr(rel.attrs);
                bool neg = _filter.negRestr(rel.attrs);

                if (!pos && !neg) return;

                uint64_t from = 0;
                uint64_t to = 0;
                uint64_t via = 0;

                for (size_t i = 0; i < rel.ways.size(); i++) {
                    if (rel.wayRoles[i] == "from") {
                        if (from) return;  // only one from member supported
                        from = rel.ways[i];
                    }
                    if (rel.wayRoles[i] == "to") {
                        if (to) return;  // only one to member supported
                        to = rel.ways[i];
                    }
                }

                for (size_t i = 0; i < rel.nodes.size(); i++) {
                    if (rel.nodeRoles[i] == "via") {
                        via = rel.nodes[i];
                        break;
                    }
                }

                if (from && to && via) {
                    if (pos) {
                        _restrictions.pos[via].emplace_back(from, to);
                    } else if (neg) {
                        _restrictions.neg[via].emplace_back(from, to);
                    }
                }
            }
        }
};

class WayHandler: public osmium::handler::Handler {
        pfaedle::trgraph::Graph& _graph;
        const pfaedle::osm::RelLst& _rels;
        const pfaedle::osm::RelMap& _wayRels;
        const pfaedle::osm::OsmFilter& _filter;
        const pfaedle::osm::OsmIdSet& _bBoxNodes;
        pfaedle::osm::NIdMap& _nodes;
        pfaedle::osm::NIdMultMap& _multiNodes;
        const pfaedle::osm::OsmIdSet& _noHupNodes;
        const pfaedle::osm::AttrKeySet& _keepAttrs;
        const pfaedle::osm::Restrictions& _restrictions;
        pfaedle::osm::Restrictor& _restrictor;
        const pfaedle::osm::FlatRels& _fl;
        pfaedle::osm::EdgTracks& _eTracks;
        const pfaedle::osm::OsmReadOpts& _opts;


        std::map<TransitEdgeLine, TransitEdgeLine *> _lines;
        std::map<size_t, TransitEdgeLine *> _relLines;

    public:
        WayHandler(pfaedle::trgraph::Graph &g,
                   const pfaedle::osm::RelLst &rels,
                   const pfaedle::osm::RelMap &wayRels,
                   const pfaedle::osm::OsmFilter &filter,
                   const pfaedle::osm::OsmIdSet &bBoxNodes,
                   pfaedle::osm::NIdMap &nodes,
                   pfaedle::osm::NIdMultMap &multiNodes,
                   const pfaedle::osm::OsmIdSet &noHupNodes,
                   const pfaedle::osm::AttrKeySet &keepAttrs,
                   const pfaedle::osm::Restrictions &rawRests,
                   pfaedle::osm::Restrictor &restor,
                   const pfaedle::osm::FlatRels &fl,
                   pfaedle::osm::EdgTracks &eTracks,
                   const pfaedle::osm::OsmReadOpts &opts) :
                _graph(g),
                _rels(rels),
                _wayRels(wayRels),
                _filter(filter),
                _bBoxNodes(bBoxNodes),
                _nodes(nodes),
                _multiNodes(multiNodes),
                _noHupNodes(noHupNodes),
                _keepAttrs(keepAttrs),
                _restrictions(rawRests),
                _restrictor(restor),
                _fl(fl),
                _eTracks(eTracks),
                _opts(opts)
        {}

        void way(const osmium::Way &way) {
            OsmWay w;
            w.id = way.id();
            for (const auto &node : way.nodes()) {
                w.nodes.emplace_back(node.ref());
            }
            for (const auto &tag : way.tags()) {
                if (_keepAttrs.count(tag.key())) {
                    w.attrs[tag.key()] = tag.value();
                }
            }

            bool valid = false;
            if (w.id && w.nodes.size() > 1 &&
                (relKeep(w.id, _wayRels, _fl) || _filter.keep(w.attrs, pfaedle::osm::OsmFilter::WAY)) &&
                !_filter.drop(w.attrs, pfaedle::osm::OsmFilter::WAY)) {
                for (auto nid : w.nodes) {
                    if (_bBoxNodes.has(nid)) {
                        valid = true;
                        break;
                    }
                }
            }

            if (valid) {
                Node *last = nullptr;
                std::vector<TransitEdgeLine *> lines;
                if (_wayRels.count(w.id)) {
                    lines = getLines(_wayRels.find(w.id)->second, _rels, _opts);
                }
                std::string track =
                        getAttrByFirstMatch(_opts.edgePlatformRules, w.id, w.attrs, _wayRels, _rels, _opts.trackNormzer);

                uint64_t lastnid = 0;

                for (auto nid : w.nodes) {
                    Node *n = nullptr;

                    if (_noHupNodes.has(nid)) {
                        n = _graph.addNd();
                        _multiNodes[nid].insert(n);
                    } else if (!_nodes.count(nid)) {
                        if (!_bBoxNodes.has(nid)) {
                            continue;
                        }
                        n = _graph.addNd();
                        _nodes[nid] = n;
                    } else {
                        n = _nodes[nid];
                    }

                    if (last) {
                        auto e = _graph.addEdg(last, n, EdgePL());
                        if (!e)
                            continue;

                        processRestrictor(nid, w.id, _restrictions, e, n, &_restrictor);
                        processRestrictor(lastnid, w.id, _restrictions, e, last, &_restrictor);

                        e->pl().addLines(lines);
                        e->pl().setLvl(_filter.level(w.attrs));
                        if (!track.empty()) {
                            _eTracks[e] = track;
                        }

                        if (_filter.oneway(w.attrs))
                            e->pl().setOneWay(1);
                        if (_filter.onewayrev(w.attrs))
                            e->pl().setOneWay(2);
                    }
                    lastnid = nid;
                    last = n;
                }
            }
        }

        static bool relKeep(uint64_t id, const pfaedle::osm::RelMap &rels, const pfaedle::osm::FlatRels &fl) {
            auto it = rels.find(id);

            if (it == rels.end())
                return false;

            return std::any_of(it->second.begin(), it->second.end(), [fl](decltype(*it->second.begin()) relId){
                return !fl.count(relId);
            });
        }

        std::vector<TransitEdgeLine *> getLines(const std::vector<size_t> &edgeRels,
                                                const pfaedle::osm::RelLst &rels,
                                                const pfaedle::osm::OsmReadOpts &ops)
        {
            std::vector<TransitEdgeLine *> ret;
            for (size_t relId : edgeRels) {
                TransitEdgeLine* transitEdgeLinePointer = nullptr;

                if (_relLines.count(relId)) {
                    transitEdgeLinePointer = _relLines[relId];
                } else {
                    TransitEdgeLine transitEdgeLine;

                    bool found = false;
                    for (const auto &r : ops.relLinerules.sNameRule) {
                        for (const auto &relAttr : rels.rels[relId]) {
                            if (relAttr.first == r) {
                                transitEdgeLine.shortName = ops.lineNormzer(pfxml::file::decode(relAttr.second));
                                if (!transitEdgeLine.shortName.empty()) found = true;
                            }
                        }
                        if (found) break;
                    }

                    found = false;
                    for (const auto &r : ops.relLinerules.fromNameRule) {
                        for (const auto &relAttr : rels.rels[relId]) {
                            if (relAttr.first == r) {
                                transitEdgeLine.fromStr = ops.statNormzer(pfxml::file::decode(relAttr.second));
                                if (!transitEdgeLine.fromStr.empty()) found = true;
                            }
                        }
                        if (found) break;
                    }

                    found = false;
                    for (const auto &r : ops.relLinerules.toNameRule) {
                        for (const auto &relAttr : rels.rels[relId]) {
                            if (relAttr.first == r) {
                                transitEdgeLine.toStr = ops.statNormzer(pfxml::file::decode(relAttr.second));
                                if (!transitEdgeLine.toStr.empty()) found = true;
                            }
                        }
                        if (found) break;
                    }

                    if (transitEdgeLine.shortName.empty() && transitEdgeLine.fromStr.empty() &&
                        transitEdgeLine.toStr.empty())
                        continue;

                    if (_lines.count(transitEdgeLine)) {
                        transitEdgeLinePointer = _lines[transitEdgeLine];
                        _relLines[relId] = transitEdgeLinePointer;
                    } else {
                        transitEdgeLinePointer = new TransitEdgeLine(transitEdgeLine);
                        _lines[transitEdgeLine] = transitEdgeLinePointer;
                        _relLines[relId] = transitEdgeLinePointer;
                    }
                }
                ret.push_back(transitEdgeLinePointer);
            }

            return ret;
        }

        static void processRestrictor(uint64_t nid, uint64_t wid,
                                 const pfaedle::osm::Restrictions &rawRests,
                                 Edge *e,
                                 Node *n,
                                 pfaedle::osm::Restrictor *restor) {
            if (rawRests.pos.count(nid)) {
                for (const auto &kv : rawRests.pos.find(nid)->second) {
                    if (kv.eFrom == wid) {
                        e->pl().setRestricted();
                        restor->add(e, kv.eTo, n, true);
                    } else if (kv.eTo == wid) {
                        e->pl().setRestricted();
                        restor->relax(wid, n, e);
                    }
                }
            }

            if (rawRests.neg.count(nid)) {
                for (const auto &kv : rawRests.neg.find(nid)->second) {
                    if (kv.eFrom == wid) {
                        e->pl().setRestricted();
                        restor->add(e, kv.eTo, n, false);
                    } else if (kv.eTo == wid) {
                        e->pl().setRestricted();
                        restor->relax(wid, n, e);
                    }
                }
            }
        }

        static std::string getAttr(const pfaedle::osm::DeepAttrRule &s, uint64_t id,
                                   const pfaedle::osm::AttrMap &attrs,
                                   const pfaedle::osm::RelMap &entRels,
                                   const pfaedle::osm::RelLst &rels) {
            if (s.relRule.kv.first.empty()) {
                if (attrs.find(s.attr) != attrs.end()) {
                    return attrs.find(s.attr)->second;
                }
            } else {
                if (entRels.count(id)) {
                    for (const auto &relId : entRels.find(id)->second) {
                        if (pfaedle::osm::OsmFilter::contained(rels.rels[relId], s.relRule.kv)) {
                            if (rels.rels[relId].count(s.attr)) {
                                return rels.rels[relId].find(s.attr)->second;
                            }
                        }
                    }
                }
            }
            return "";
        }

        static std::string getAttrByFirstMatch(const pfaedle::osm::DeepAttrLst &rule, uint64_t id,
                                               const pfaedle::osm::AttrMap &attrs,
                                               const pfaedle::osm::RelMap &entRels,
                                               const pfaedle::osm::RelLst &rels,
                                               const Normalizer &norm) {
            std::string ret;
            for (const auto &s : rule) {
                ret = norm(pfxml::file::decode(getAttr(s, id, attrs, entRels, rels)));
                if (!ret.empty()) return ret;
            }

            return ret;
        }
};

class NodeHandler2: public osmium::handler::Handler {
        pfaedle::trgraph::Graph& _graph;
        const pfaedle::osm::RelLst& rels;
        const pfaedle::osm::RelMap& nodeRels;
        const pfaedle::osm::OsmFilter& filter;
        const pfaedle::osm::OsmIdSet& bBoxNodes;
        pfaedle::osm::NIdMap& nodes;
        pfaedle::osm::NIdMultMap& multNodes;
        pfaedle::osm::NodeSet& orphanStations;
        const pfaedle::osm::AttrKeySet& keepAttrs;
        const pfaedle::osm::FlatRels& fl;
        const pfaedle::osm::OsmReadOpts& opts;

    public:
        NodeHandler2(Graph &g,
                     const pfaedle::osm::RelLst &rels,
                     const pfaedle::osm::RelMap &nodeRels,
                     const pfaedle::osm::OsmFilter &filter,
                     const pfaedle::osm::OsmIdSet &bBoxNodes,
                     pfaedle::osm::NIdMap &nodes,
                     pfaedle::osm::NIdMultMap &multNodes,
                     pfaedle::osm::NodeSet &orphanStations,
                     const pfaedle::osm::AttrKeySet &keepAttrs,
                     const pfaedle::osm::FlatRels &fl,
                     const pfaedle::osm::OsmReadOpts &opts) :
            _graph(g),
            rels(rels),
            nodeRels(nodeRels),
            filter(filter),
            bBoxNodes(bBoxNodes),
            nodes(nodes),
            multNodes(multNodes),
            orphanStations(orphanStations),
            keepAttrs(keepAttrs),
            fl(fl),
            opts(opts)
        {}

        void node(const osmium::Node &node) {
            OsmNode nd;
            nd.id = node.id();
            nd.lat = node.location().lat();
            nd.lng = node.location().lon();

            for (const auto &tag : node.tags()) {
                if (keepAttrs.count(tag.key()))
                    nd.attrs[tag.key()] = tag.value();
            }

            bool valid = false;
            if (nd.id &&
                    (nodes.count(nd.id) || multNodes.count(nd.id) ||
                     relKeep(nd.id, nodeRels, fl) || filter.keep(nd.attrs, pfaedle::osm::OsmFilter::NODE)) &&
                    (nodes.count(nd.id) || bBoxNodes.has(nd.id)) &&
                    (nodes.count(nd.id) || multNodes.count(nd.id) ||
                     !filter.drop(nd.attrs, pfaedle::osm::OsmFilter::NODE))) {
                valid = true;
            }

            pfaedle::osm::StAttrGroups attrGroups;
            if (valid) {
                Node *n = nullptr;
                auto pos = util::geo::latLngToWebMerc<PFAEDLE_PRECISION>(nd.lat, nd.lng);

                if (nodes.count(nd.id)) {
                    n = nodes[nd.id];
                    n->pl().setGeom(pos);

                    if (filter.station(nd.attrs)) {
                        auto si = getStatInfo(n, nd.id, pos, nd.attrs, &attrGroups, nodeRels, rels, opts);

                        if (!si.isNull())
                            n->pl().setSI(si);
                    } else if (filter.blocker(nd.attrs)) {
                        n->pl().setBlocker();
                    }
                } else if (multNodes.count(nd.id)) {
                    for (auto *n : multNodes[nd.id]) {
                        n->pl().setGeom(pos);
                        if (filter.station(nd.attrs)) {
                            auto si = getStatInfo(n, nd.id, pos, nd.attrs, &attrGroups, nodeRels, rels, opts);

                            if (!si.isNull())
                                n->pl().setSI(si);
                        } else if (filter.blocker(nd.attrs)) {
                            n->pl().setBlocker();
                        }
                    }
                } else {
                    // these are nodes without any connected edges
                    if (filter.station(nd.attrs)) {
                        auto tmp = _graph.addNd(NodePL(pos));
                        auto si = getStatInfo(tmp, nd.id, pos, nd.attrs, &attrGroups, nodeRels, rels, opts);

                        if (!si.isNull())
                            tmp->pl().setSI(si);

                        if (tmp->pl().getSI()) {
                            tmp->pl().getSI()->setIsFromOsm(false);
                            orphanStations.insert(tmp);
                        }
                    }
                }
            }
        }


        static bool relKeep(uint64_t id, const pfaedle::osm::RelMap &rels, const pfaedle::osm::FlatRels &fl) {
            auto it = rels.find(id);

            if (it == rels.end())
                return false;

            for (uint64_t relId : it->second) {
                // as soon as any of this entities relations is not flat, return true
                if (!fl.count(relId))
                    return true;
            }

            return false;
        }

        static Nullable<StatInfo> getStatInfo(Node *node, uint64_t nid,
                                              const POINT &pos,
                                              const pfaedle::osm::AttrMap &m,
                                              pfaedle::osm::StAttrGroups *groups,
                                              const pfaedle::osm::RelMap &nodeRels,
                                              const pfaedle::osm::RelLst &rels,
                                              const pfaedle::osm::OsmReadOpts &ops) {
            std::vector<std::string> names =
                    getAttrMatchRanked(ops.statAttrRules.nameRule, nid, m, nodeRels, rels, ops.statNormzer);
            std::string platform =
                    getAttrByFirstMatch(ops.statAttrRules.platformRule, nid, m, nodeRels, rels, ops.trackNormzer);

            if (names.empty()) {
                return Nullable<StatInfo>();
            }

            auto ret = StatInfo(names[0], platform, true);

#ifdef PFAEDLE_STATION_IDS
            ret.setId(getAttrByFirstMatch(ops.statAttrRules.idRule, nid, m, nodeRels, rels, ops.idNormzer));
#endif

            for (size_t i = 1; i < names.size(); i++) {
                ret.addAltName(names[i]);
            }
            bool groupFound = false;

            for (const auto &rule : ops.statGroupNAttrRules) {
                if (groupFound) break;
                std::string ruleVal = getAttr(rule.attr, nid, m, nodeRels, rels);
                if (!ruleVal.empty()) {
                    // check if a matching group exists
                    for (auto *group : (*groups)[rule.attr.attr][ruleVal]) {
                        if (groupFound) break;
                        for (const auto *member : group->getNodes()) {
                            if (webMercMeterDist(*member->pl().getGeom(), pos) <= rule.maxDist) {
                                // ok, group is matching
                                groupFound = true;
                                if (node) group->addNode(node);
                                ret.setGroup(group);
                                break;
                            }
                        }
                    }
                }
            }

            if (!groupFound) {
                for (const auto &rule : ops.statGroupNAttrRules) {
                    std::string ruleVal = getAttr(rule.attr, nid, m, nodeRels, rels);
                    if (!ruleVal.empty()) {
                        // add new group
                        auto *group = new StatGroup();
                        if (node) group->addNode(node);
                        ret.setGroup(group);
                        (*groups)[rule.attr.attr][ruleVal].push_back(group);
                        break;
                    }
                }
            }

            return ret;
        }

        static std::string getAttrByFirstMatch(const pfaedle::osm::DeepAttrLst &rule,
                                               uint64_t id,
                                               const pfaedle::osm::AttrMap &attrs,
                                               const pfaedle::osm::RelMap &entRels,
                                               const pfaedle::osm::RelLst &rels,
                                               const Normalizer &norm) {
            std::string ret;
            for (const auto &s : rule) {
                ret = norm(pfxml::file::decode(getAttr(s, id, attrs, entRels, rels)));
                if (!ret.empty()) {
                    break;
                }
            }

            return ret;
        }

        static std::vector<std::string> getAttrMatchRanked(const pfaedle::osm::DeepAttrLst &rule,
                                                           uint64_t id,
                                                           const pfaedle::osm::AttrMap &attrs,
                                                           const pfaedle::osm::RelMap &entRels,
                                                           const pfaedle::osm::RelLst &rels,
                                                           const Normalizer &norm) {
            std::vector<std::string> ret;
            for (const auto &s : rule) {
                std::string tmp = norm(pfxml::file::decode(getAttr(s, id, attrs, entRels, rels)));
                if (!tmp.empty()) {
                    ret.push_back(tmp);
                }
            }

            return ret;
        }

        static std::string getAttr(const pfaedle::osm::DeepAttrRule &s,
                                   uint64_t id,
                                   const pfaedle::osm::AttrMap &attrs,
                                   const pfaedle::osm::RelMap &entRels,
                                   const pfaedle::osm::RelLst &rels) {
            if (s.relRule.kv.first.empty()) {
                if (attrs.find(s.attr) != attrs.end()) {
                    return attrs.find(s.attr)->second;
                }
            } else {
                if (entRels.count(id)) {
                    for (const auto &relId : entRels.find(id)->second) {
                        if (pfaedle::osm::OsmFilter::contained(rels.rels[relId], s.relRule.kv)) {
                            if (rels.rels[relId].count(s.attr)) {
                                return rels.rels[relId].find(s.attr)->second;
                            }
                        }
                    }
                }
            }
            return "";
        }
};

// _____________________________________________________________________________
bool EqSearch::operator()(const Node* cand, const StatInfo* si) const {
    if (orphanSnap && cand->pl().getSI() &&
            (!cand->pl().getSI()->getGroup() || cand->pl().getSI()->getGroup()->getStops().empty())) {
        return true;
    }
    return cand->pl().getSI() && cand->pl().getSI()->simi(si) > minSimi;
}

// _____________________________________________________________________________
OsmBuilder::OsmBuilder() = default;

// _____________________________________________________________________________
void OsmBuilder::read(const std::string& path, const OsmReadOpts& opts,
                      Graph* g, const BBoxIdx& bbox, size_t gridSize,
                      router::FeedStops* fs, Restrictor* res) {
    if (!bbox.size()) return;

    LOG(INFO) << "Reading OSM file " << path << " ... ";

    NodeSet orphanStations;
    EdgTracks eTracks;
    {
        OsmIdSet bboxNodes, noHupNodes;

        NIdMap nodes;
        NIdMultMap multNodes;
        RelLst intmRels;
        RelMap nodeRels, wayRels;

        Restrictions rawRests;

        AttrKeySet attrKeys[3] = {};
        getKeptAttrKeys(opts, attrKeys);

        OsmFilter filter(opts);

        osmium::io::Reader reader {path,
                    osmium::osm_entity_bits::node |
                    osmium::osm_entity_bits::way |
                    osmium::osm_entity_bits::relation};
        osmium::io::Reader reader_nodes{path,
                    osmium::osm_entity_bits::node};

        NodeHandler nodeHandler(filter, bbox, bboxNodes, noHupNodes);
        RelationHandler relationHandler(filter, bbox, attrKeys[2], intmRels, nodeRels, wayRels, rawRests);
        WayHandler wayHandler(*g, intmRels, wayRels, filter, bboxNodes, nodes, multNodes, noHupNodes,
                              attrKeys[1], rawRests, *res, intmRels.flat, eTracks, opts);

        osmium::apply(reader, nodeHandler, relationHandler, wayHandler);
        NodeHandler2 nodeHandler2(*g, intmRels, nodeRels, filter, bboxNodes, nodes, multNodes, orphanStations,
                                  attrKeys[0], intmRels.flat, opts);
        osmium::apply(reader_nodes, nodeHandler2);

        // we do four passes of the file here to be as memory creedy as possible:
        // - the first pass collects all node IDs which are
        //    * inside the given bounding box
        //    * (TODO: maybe more filtering?)
        //   these nodes are stored on the HD via OsmIdSet (which implements a
        //   simple bloom filter / base 256 encoded id store
        // - the second pass collects filtered relations
        // - the third pass collects filtered ways which contain one of the nodes
        //   from pass 1
        // - the forth pass collects filtered nodes which were
        //    * collected as node ids in pass 1
        //    * match the filter criteria
        //    * have been used in a way in pass 3
    }

    LOG(VDEBUG) << "OSM ID set lookups: " << osm::OsmIdSet::LOOKUPS
                << ", file lookups: " << osm::OsmIdSet::FLOOKUPS;

    LOG(VDEBUG) << "Applying edge track numbers...";
    writeEdgeTracks(eTracks);
    eTracks.clear();

    {
        LOG(VDEBUG) << "Fixing gaps...";
        NodeGrid ng = buildNodeIdx(g, gridSize, bbox.getFullWebMercBox(), false);
        fixGaps(g, &ng);
    }

    LOG(VDEBUG) << "Writing edge geoms...";
    writeGeoms(g);

    LOG(VDEBUG) << "Snapping stations...";
    snapStats(opts, g, bbox, gridSize, fs, res, orphanStations);

    LOG(VDEBUG) << "Deleting orphan nodes...";
    deleteOrphNds(g);

    LOG(VDEBUG) << "Deleting orphan edges...";
    deleteOrphEdgs(g, opts);

    LOG(VDEBUG) << "Collapsing edges...";
    collapseEdges(g);

    LOG(VDEBUG) << "Deleting orphan nodes...";
    deleteOrphNds(g);

    LOG(VDEBUG) << "Deleting orphan edges...";
    deleteOrphEdgs(g, opts);

    LOG(VDEBUG) << "Writing graph components...";
    // the restrictor is needed here to prevent connections in the graph
    // which are not possible in reality
    uint32_t comps = writeComps(g);

    LOG(VDEBUG) << "Simplifying geometries...";
    simplifyGeoms(g);

    LOG(VDEBUG) << "Writing other-direction edges...";
    writeODirEdgs(g, res);

    LOG(VDEBUG) << "Write dummy node self-edges...";
    writeSelfEdgs(g);

    size_t numEdges = 0;

    for (auto* n : *g->getNds()) {
        numEdges += n->getAdjListOut().size();
    }

    LOG(DEBUG) << "Graph has " << g->getNds()->size() << " nodes, " << numEdges
               << " edges and " << comps << " connected component(s)";
}

// _____________________________________________________________________________
void OsmBuilder::overpassQryWrite(std::ostream* out,
                                  const std::vector<OsmReadOpts>& opts,
                                  const BBoxIdx& latLngBox) const {
    OsmIdSet bboxNodes, noHupNodes;
    MultAttrMap emptyF;

    RelLst rels;
    OsmIdList ways;
    RelMap nodeRels, wayRels;

    // TODO(patrick): not needed here!
    Restrictions rests;

    NIdMap nodes;

    // always empty
    NIdMultMap multNodes;
    util::xml::XmlWriter wr(out, true, 4);

    *out << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    wr.openComment();
    wr.writeText(" - written by pfaedle -");
    wr.closeTag();
    wr.openTag("osm-script");

    OsmFilter filter;

    for (const OsmReadOpts& o : opts) {
        filter = filter.merge(OsmFilter(o.keepFilter, o.dropFilter));
    }

    wr.openTag("union");
    size_t c = 0;
    for (auto box : latLngBox.getLeafs()) {
        if (box.getLowerLeft().getX() > box.getUpperRight().getX()) continue;
        c++;
        wr.openComment();
        wr.writeText(std::string("Bounding box #") + std::to_string(c) + " (" +
                     std::to_string(box.getLowerLeft().getY()) + ", " +
                     std::to_string(box.getLowerLeft().getX()) + ", " +
                     std::to_string(box.getUpperRight().getY()) + ", " +
                     std::to_string(box.getUpperRight().getX()) + ")");
        wr.closeTag();
        for (auto t : std::vector<std::string>{"way", "node", "relation"}) {
            for (auto r : filter.getKeepRules()) {
                for (auto val : r.second) {
                    if (t == "way" && (val.second & OsmFilter::WAY)) continue;
                    if (t == "relation" && (val.second & OsmFilter::REL)) continue;
                    if (t == "node" && (val.second & OsmFilter::NODE)) continue;

                    wr.openTag("query", {{"type", t}});
                    if (val.first == "*")
                        wr.openTag("has-kv", {{"k", r.first}});
                    else
                        wr.openTag("has-kv", {{"k", r.first}, {"v", val.first}});
                    wr.closeTag();
                    wr.openTag("bbox-query",
                    {{"s", std::to_string(box.getLowerLeft().getY())},
                     {"w", std::to_string(box.getLowerLeft().getX())},
                     {"n", std::to_string(box.getUpperRight().getY())},
                     {"e", std::to_string(box.getUpperRight().getX())}});
                    wr.closeTag();
                    wr.closeTag();
                }
            }
        }
    }

    wr.closeTag();

    wr.openTag("union");
    wr.openTag("item");
    wr.closeTag();
    wr.openTag("recurse", {{"type", "down"}});
    wr.closeTag();
    wr.closeTag();
    wr.openTag("print");

    wr.closeTags();
}

// _____________________________________________________________________________
void OsmBuilder::filterWrite(const std::string& in, const std::string& out,
                             const std::vector<OsmReadOpts>& opts,
                             const BBoxIdx& latLngBox) {
    OsmIdSet bboxNodes, noHupNodes;
    MultAttrMap emptyF;

    RelLst rels;
    OsmIdList ways;
    RelMap nodeRels, wayRels;

    // TODO(patrick): not needed here!
    Restrictions rests;

    NIdMap nodes;

    // always empty
    NIdMultMap multNodes;

    pfxml::file xml(in);
    std::ofstream outstr;
    outstr.open(out);

    util::xml::XmlWriter wr(&outstr, true, 4);

    outstr << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    wr.openTag("osm");
    wr.openTag(
                "bounds",
    {{"minlat", std::to_string(latLngBox.getFullBox().getLowerLeft().getY())},
     {"minlon", std::to_string(latLngBox.getFullBox().getLowerLeft().getX())},
     {"maxlat",
      std::to_string(latLngBox.getFullBox().getUpperRight().getY())},
     {"maxlon",
      std::to_string(latLngBox.getFullBox().getUpperRight().getX())}});
    wr.closeTag();

    OsmFilter filter;
    AttrKeySet attrKeys[3] = {};

    for (const OsmReadOpts& o : opts) {
        getKeptAttrKeys(o, attrKeys);
        filter = filter.merge(OsmFilter(o.keepFilter, o.dropFilter));
    }

    skipUntil(&xml, "node");
    pfxml::parser_state nodeBeg = xml.state();
    pfxml::parser_state edgesBeg =
            readBBoxNds(&xml, &bboxNodes, &noHupNodes, filter, latLngBox);

    skipUntil(&xml, "relation");
    readRels(&xml, &rels, &nodeRels, &wayRels, filter, attrKeys[2], &rests);

    xml.set_state(edgesBeg);
    readEdges(&xml, wayRels, filter, bboxNodes, attrKeys[1], &ways, &nodes,
            rels.flat);

    xml.set_state(nodeBeg);

    readWriteNds(&xml, &wr, nodeRels, filter, bboxNodes, &nodes, attrKeys[0],
            rels.flat);
    readWriteWays(&xml, &wr, &ways, attrKeys[1]);

    std::sort(ways.begin(), ways.end());
    skipUntil(&xml, "relation");
    readWriteRels(&xml, &wr, &ways, &nodes, filter, attrKeys[2]);

    wr.closeTags();
}

// _____________________________________________________________________________
void OsmBuilder::readWriteRels(pfxml::file* i, util::xml::XmlWriter* o,
                               OsmIdList* ways, NIdMap* nodes,
                               const OsmFilter& filter,
                               const AttrKeySet& keepAttrs) {
    OsmRel rel;
    while ((rel = nextRel(i, filter, keepAttrs)).id) {
        OsmIdList realNodes, realWays;
        std::vector<const char *> realNodeRoles, realWayRoles;

        for (size_t j = 0; j < rel.ways.size(); j++) {
            osmid wid = rel.ways[j];
            const auto& i = std::lower_bound(ways->begin(), ways->end(), wid);
            if (i != ways->end() && *i == wid) {
                realWays.push_back(wid);
                realWayRoles.push_back(rel.wayRoles[j].c_str());
            }
        }

        for (size_t j = 0; j < rel.nodes.size(); j++) {
            osmid nid = rel.nodes[j];
            if (nodes->count(nid)) {
                realNodes.push_back(nid);
                realNodeRoles.push_back(rel.nodeRoles[j].c_str());
            }
        }

        if (!realNodes.empty() || !realWays.empty()) {
            o->openTag("relation", "id", std::to_string(rel.id));

            for (size_t j = 0; j < realNodes.size(); j++) {
                osmid nid = realNodes[j];
                std::map<std::string, std::string> attrs;
                attrs["type"] = "node";
                if (strlen(realNodeRoles[j])) attrs["role"] = realNodeRoles[j];
                attrs["ref"] = std::to_string(nid);
                o->openTag("member", attrs);
                o->closeTag();
            }

            for (size_t j = 0; j < realWays.size(); j++) {
                osmid wid = realWays[j];
                std::map<std::string, std::string> attrs;
                attrs["type"] = "way";
                if (strlen(realWayRoles[j])) attrs["role"] = realWayRoles[j];
                attrs["ref"] = std::to_string(wid);
                o->openTag("member", attrs);
                o->closeTag();
            }

            for (const auto& kv : rel.attrs) {
                std::map<std::string, std::string> attrs = {
                    {"k", kv.first}, {"v", pfxml::file::decode(kv.second)}};
                o->openTag("tag", attrs);
                o->closeTag();
            }

            o->closeTag();
        }
    }
}

// _____________________________________________________________________________
void OsmBuilder::readWriteWays(pfxml::file* i, util::xml::XmlWriter* o,
                               OsmIdList* ways,
                               const AttrKeySet& keepAttrs) const {
    OsmWay w;
    NIdMultMap empty;
    for (auto wid : *ways) {
        w = nextWayWithId(i, wid, keepAttrs);
        assert(w.id);
        o->openTag("way", "id", std::to_string(wid));
        for (osmid nid : w.nodes) {
            o->openTag("nd", "ref", std::to_string(nid));
            o->closeTag();
        }
        for (const auto& kv : w.attrs) {
            std::map<std::string, std::string> attrs;
            attrs["k"] = kv.first;
            attrs["v"] = pfxml::file::decode(kv.second);
            o->openTag("tag", attrs);
            o->closeTag();
        }
        o->closeTag();
    }
}

// _____________________________________________________________________________
NodePL OsmBuilder::plFromGtfs(const Stop* s, const OsmReadOpts& ops) {
    NodePL ret(
                util::geo::latLngToWebMerc<PFAEDLE_PRECISION>(s->getLat(), s->getLng()),
                StatInfo(ops.statNormzer(s->getName()), ops.trackNormzer(s->getPlatformCode()), false));

#ifdef PFAEDLE_STATION_IDS
    // debug feature, store station id from GTFS
    ret.getSI()->setId(s->getId());
#endif

    if (s->getParentStation()) {
        ret.getSI()->addAltName(ops.statNormzer(s->getParentStation()->getName()));
    }

    return ret;
}

// _____________________________________________________________________________
pfxml::parser_state OsmBuilder::readBBoxNds(pfxml::file* xml, OsmIdSet* nodes,
                                            OsmIdSet* nohupNodes,
                                            const OsmFilter& filter,
                                            const BBoxIdx& bbox) const {
    bool inNodeBlock = false;
    uint64_t curId = 0;

    do {
        const pfxml::tag& cur = xml->get();

        if (inNodeBlock && xml->level() == 3 && curId &&
                strcmp(cur.name, "tag") == 0) {
            if (filter.nohup(cur.attrs.find("k")->second,
                             cur.attrs.find("v")->second)) {
                nohupNodes->add(curId);
            }
        }

        if (xml->level() != 2) continue;
        if (!inNodeBlock && strcmp(cur.name, "node") == 0) inNodeBlock = true;

        if (inNodeBlock) {
            // block ended
            if (strcmp(cur.name, "node") != 0)
                return xml->state();
            double y = util::atof(cur.attrs.find("lat")->second, 7);
            double x = util::atof(cur.attrs.find("lon")->second, 7);

            if (bbox.contains(Point<double>(x, y))) {
                curId = util::atoul(cur.attrs.find("id")->second);
                nodes->add(curId);
            }
        }
    } while (xml->next());

    return xml->state();
}

// _____________________________________________________________________________
OsmWay OsmBuilder::nextWayWithId(pfxml::file* xml, osmid wid,
                                 const AttrKeySet& keepAttrs) const {
    OsmWay w;

    do {
        const pfxml::tag& cur = xml->get();
        if (xml->level() == 2 || xml->level() == 0) {
            if (w.id || strcmp(cur.name, "way")) return w;

            osmid id = util::atoul(cur.attrs.find("id")->second);
            if (id == wid) w.id = id;
        }

        if (w.id && xml->level() == 3) {
            if (strcmp(cur.name, "nd") == 0) {
                w.nodes.push_back(util::atoul(cur.attrs.find("ref")->second));
            } else if (strcmp(cur.name, "tag") == 0) {
                if (keepAttrs.count(cur.attrs.find("k")->second))
                    w.attrs[cur.attrs.find("k")->second] = cur.attrs.find("v")->second;
            }
        }
    } while (xml->next());

    if (w.id) return w;

    return OsmWay();
}

// _____________________________________________________________________________
void OsmBuilder::skipUntil(pfxml::file* xml, const std::string& s) const {
    while (xml->next() && strcmp(xml->get().name, s.c_str()) != 0) {
    }
}

// _____________________________________________________________________________
bool OsmBuilder::relKeep(osmid id, const RelMap& rels,
                         const FlatRels& fl) const {
    auto it = rels.find(id);

    if (it == rels.end()) return false;

    for (osmid relId : it->second) {
        // as soon as any of this entities relations is not flat, return true
        if (!fl.count(relId)) return true;
    }

    return false;
}

// _____________________________________________________________________________
OsmWay OsmBuilder::nextWay(pfxml::file* xml, const RelMap& wayRels,
                           const OsmFilter& filter, const OsmIdSet& bBoxNodes,
                           const AttrKeySet& keepAttrs,
                           const FlatRels& fl) const {
    OsmWay w;

    do {
        const pfxml::tag& cur = xml->get();
        if (xml->level() == 2 || xml->level() == 0) {
            if (keepWay(w, wayRels, filter, bBoxNodes, fl))
                return w;

            if (strcmp(cur.name, "way") != 0)
                return OsmWay();

            w.id = util::atoul(cur.attrs.find("id")->second);
            w.nodes.clear();
            w.attrs.clear();
        }

        if (w.id && xml->level() == 3) {
            if (strcmp(cur.name, "nd") == 0) {
                osmid nid = util::atoul(cur.attrs.find("ref")->second);
                w.nodes.push_back(nid);
            } else if (strcmp(cur.name, "tag") == 0) {
                if (keepAttrs.count(cur.attrs.find("k")->second)) {
                    w.attrs[cur.attrs.find("k")->second] = cur.attrs.find("v")->second;
                }
            }
        }
    } while (xml->next());

    if (keepWay(w, wayRels, filter, bBoxNodes, fl))
        return w;

    return OsmWay();
}

// _____________________________________________________________________________
bool OsmBuilder::keepWay(const OsmWay& w, const RelMap& wayRels,
                         const OsmFilter& filter, const OsmIdSet& bBoxNodes,
                         const FlatRels& fl) const {
    if (w.id && w.nodes.size() > 1 &&
            (relKeep(w.id, wayRels, fl) || filter.keep(w.attrs, OsmFilter::WAY)) &&
            !filter.drop(w.attrs, OsmFilter::WAY)) {
        for (osmid nid : w.nodes) {
            if (bBoxNodes.has(nid)) {
                return true;
            }
        }
    }

    return false;
}

// _____________________________________________________________________________
void OsmBuilder::readEdges(pfxml::file* xml, const RelMap& wayRels,
                           const OsmFilter& filter, const OsmIdSet& bBoxNodes,
                           const AttrKeySet& keepAttrs, OsmIdList* ret,
                           NIdMap* nodes, const FlatRels& flat) {
    OsmWay w;
    while ((w = nextWay(xml, wayRels, filter, bBoxNodes, keepAttrs, flat)).id) {
        ret->push_back(w.id);
        for (auto n : w.nodes) {
            (*nodes)[n] = nullptr;
        }
    }
}

// _____________________________________________________________________________
void OsmBuilder::readEdges(pfxml::file* xml, Graph* g, const RelLst& rels,
                           const RelMap& wayRels, const OsmFilter& filter,
                           const OsmIdSet& bBoxNodes, NIdMap* nodes,
                           NIdMultMap* multiNodes, const OsmIdSet& noHupNodes,
                           const AttrKeySet& keepAttrs,
                           const Restrictions& rawRests, Restrictor* restor,
                           const FlatRels& fl, EdgTracks* eTracks,
                           const OsmReadOpts& opts) {
    OsmWay w;
    while ((w = nextWay(xml, wayRels, filter, bBoxNodes, keepAttrs, fl)).id) {
        Node* last = nullptr;
        std::vector<TransitEdgeLine*> lines;
        if (wayRels.count(w.id)) {
            lines = getLines(wayRels.find(w.id)->second, rels, opts);
        }
        std::string track =
                getAttrByFirstMatch(opts.edgePlatformRules, w.id, w.attrs, wayRels,
                                    rels, opts.trackNormzer);

        osmid lastnid = 0;
        for (osmid nid : w.nodes) {
            Node* n = nullptr;
            if (noHupNodes.has(nid)) {
                n = g->addNd();
                (*multiNodes)[nid].insert(n);
            } else if (!nodes->count(nid)) {
                if (!bBoxNodes.has(nid)) continue;
                n = g->addNd();
                (*nodes)[nid] = n;
            } else {
                n = (*nodes)[nid];
            }
            if (last) {
                auto e = g->addEdg(last, n, EdgePL());
                if (!e) continue;

                processRestr(nid, w.id, rawRests, e, n, restor);
                processRestr(lastnid, w.id, rawRests, e, last, restor);

                e->pl().addLines(lines);
                e->pl().setLvl(filter.level(w.attrs));
                if (!track.empty()) (*eTracks)[e] = track;

                if (filter.oneway(w.attrs)) e->pl().setOneWay(1);
                if (filter.onewayrev(w.attrs)) e->pl().setOneWay(2);
            }
            lastnid = nid;
            last = n;
        }
    }
}

// _____________________________________________________________________________
void OsmBuilder::processRestr(osmid nid, osmid wid,
                              const Restrictions& rawRests, Edge* e, Node* n,
                              Restrictor* restor) const {
    if (rawRests.pos.count(nid)) {
        for (const auto& kv : rawRests.pos.find(nid)->second) {
            if (kv.eFrom == wid) {
                e->pl().setRestricted();
                restor->add(e, kv.eTo, n, true);
            } else if (kv.eTo == wid) {
                e->pl().setRestricted();
                restor->relax(wid, n, e);
            }
        }
    }

    if (rawRests.neg.count(nid)) {
        for (const auto& kv : rawRests.neg.find(nid)->second) {
            if (kv.eFrom == wid) {
                e->pl().setRestricted();
                restor->add(e, kv.eTo, n, false);
            } else if (kv.eTo == wid) {
                e->pl().setRestricted();
                restor->relax(wid, n, e);
            }
        }
    }
}

// _____________________________________________________________________________
OsmNode OsmBuilder::nextNode(pfxml::file* xml, NIdMap* nodes,
                             NIdMultMap* multNodes, const RelMap& nodeRels,
                             const OsmFilter& filter, const OsmIdSet& bBoxNodes,
                             const AttrKeySet& keepAttrs,
                             const FlatRels& fl) const {
    OsmNode n;

    do {
        const pfxml::tag& cur = xml->get();
        if (xml->level() == 2 || xml->level() == 0) {
            if (keepNode(n, *nodes, *multNodes, nodeRels, bBoxNodes, filter, fl))
                return n;
            // block ended
            if (strcmp(cur.name, "node")) return OsmNode();

            n.attrs.clear();
            n.lat = util::atof(cur.attrs.find("lat")->second, 7);
            n.lng = util::atof(cur.attrs.find("lon")->second, 7);
            n.id = util::atoul(cur.attrs.find("id")->second);
        }

        if (xml->level() == 3 && n.id && strcmp(cur.name, "tag") == 0) {
            if (keepAttrs.count(cur.attrs.find("k")->second))
                n.attrs[cur.attrs.find("k")->second] = cur.attrs.find("v")->second;
        }
    } while (xml->next());

    if (keepNode(n, *nodes, *multNodes, nodeRels, bBoxNodes, filter, fl))
        return n;
    return OsmNode();
}

// _____________________________________________________________________________
bool OsmBuilder::keepNode(const OsmNode& n, const NIdMap& nodes,
                          const NIdMultMap& multNodes, const RelMap& nodeRels,
                          const OsmIdSet& bBoxNodes, const OsmFilter& filter,
                          const FlatRels& fl) const {
    if (n.id &&
            (nodes.count(n.id) || multNodes.count(n.id) ||
             relKeep(n.id, nodeRels, fl) || filter.keep(n.attrs, OsmFilter::NODE)) &&
            (nodes.count(n.id) || bBoxNodes.has(n.id)) &&
            (nodes.count(n.id) || multNodes.count(n.id) ||
             !filter.drop(n.attrs, OsmFilter::NODE))) {
        return true;
    }

    return false;
}

// _____________________________________________________________________________
void OsmBuilder::readWriteNds(pfxml::file* i, util::xml::XmlWriter* o,
                              const RelMap& nRels, const OsmFilter& filter,
                              const OsmIdSet& bBoxNds, NIdMap* nds,
                              const AttrKeySet& keepAttrs,
                              const FlatRels& f) const {
    OsmNode nd;
    NIdMultMap empt;
    while (
           (nd = nextNode(i, nds, &empt, nRels, filter, bBoxNds, keepAttrs, f)).id) {
        (*nds)[nd.id] = 0;
        o->openTag("node", {{"id", std::to_string(nd.id)},
                            {"lat", std::to_string(nd.lat)},
                            {"lon", std::to_string(nd.lng)}});
        for (const auto& kv : nd.attrs) {
            o->openTag("tag", {{"k", kv.first}, {"v", pfxml::file::decode(kv.second)}});
            o->closeTag();
        }
        o->closeTag();
    }
}

// _____________________________________________________________________________
void OsmBuilder::readNodes(pfxml::file* xml, Graph* g, const RelLst& rels,
                           const RelMap& nodeRels, const OsmFilter& filter,
                           const OsmIdSet& bBoxNodes, NIdMap* nodes,
                           NIdMultMap* multNodes, NodeSet* orphanStations,
                           const AttrKeySet& keepAttrs, const FlatRels& fl,
                           const OsmReadOpts& opts) const {
    StAttrGroups attrGroups;

    OsmNode nd;
    while ((nd = nextNode(xml, nodes, multNodes, nodeRels, filter, bBoxNodes,
                          keepAttrs, fl))
           .id) {
        Node* n = 0;
        auto pos = util::geo::latLngToWebMerc<PFAEDLE_PRECISION>(nd.lat, nd.lng);
        if (nodes->count(nd.id)) {
            n = (*nodes)[nd.id];
            n->pl().setGeom(pos);
            if (filter.station(nd.attrs)) {
                auto si = getStatInfo(n, nd.id, pos, nd.attrs, &attrGroups, nodeRels,
                                      rels, opts);
                if (!si.isNull()) n->pl().setSI(si);
            } else if (filter.blocker(nd.attrs)) {
                n->pl().setBlocker();
            }
        } else if ((*multNodes).count(nd.id)) {
            for (auto* n : (*multNodes)[nd.id]) {
                n->pl().setGeom(pos);
                if (filter.station(nd.attrs)) {
                    auto si = getStatInfo(n, nd.id, pos, nd.attrs, &attrGroups, nodeRels,
                                          rels, opts);
                    if (!si.isNull()) n->pl().setSI(si);
                } else if (filter.blocker(nd.attrs)) {
                    n->pl().setBlocker();
                }
            }
        } else {
            // these are nodes without any connected edges
            if (filter.station(nd.attrs)) {
                auto tmp = g->addNd(NodePL(pos));
                auto si = getStatInfo(tmp, nd.id, pos, nd.attrs, &attrGroups, nodeRels,
                                      rels, opts);
                if (!si.isNull()) tmp->pl().setSI(si);
                if (tmp->pl().getSI()) {
                    tmp->pl().getSI()->setIsFromOsm(false);
                    orphanStations->insert(tmp);
                }
            }
        }
    }
}

// _____________________________________________________________________________
OsmRel OsmBuilder::nextRel(pfxml::file* xml, const OsmFilter& filter,
                           const AttrKeySet& keepAttrs) const {
    OsmRel rel;

    do {
        const pfxml::tag& cur = xml->get();
        if (xml->level() == 2 || xml->level() == 0) {
            uint64_t keepFlags = 0;
            uint64_t dropFlags = 0;
            if (rel.id && !rel.attrs.empty() &&
                    (keepFlags = filter.keep(rel.attrs, OsmFilter::REL)) &&
                    !(dropFlags = filter.drop(rel.attrs, OsmFilter::REL))) {
                rel.keepFlags = keepFlags;
                rel.dropFlags = dropFlags;
                return rel;
            }

            // block ended
            if (strcmp(cur.name, "relation") != 0) {
                return OsmRel();
            }

            rel.attrs.clear();
            rel.nodes.clear();
            rel.ways.clear();
            rel.nodeRoles.clear();
            rel.wayRoles.clear();
            rel.keepFlags = 0;
            rel.dropFlags = 0;
            rel.id = util::atoul(cur.attrs.find("id")->second);
        }

        if (xml->level() == 3 && rel.id) {
            if (strcmp(cur.name, "member") == 0) {
                if (strcmp(cur.attrs.find("type")->second, "node") == 0) {
                    osmid id = util::atoul(cur.attrs.find("ref")->second);
                    // TODO(patrick): no need to push IDs that have been filtered out by
                    // the bounding box!!!!
                    rel.nodes.push_back(id);
                    if (cur.attrs.count("role")) {
                        rel.nodeRoles.emplace_back(cur.attrs.find("role")->second);
                    } else {
                        rel.nodeRoles.emplace_back("");
                    }
                }

                if (strcmp(cur.attrs.find("type")->second, "way") == 0) {
                    osmid id = util::atoul(cur.attrs.find("ref")->second);
                    rel.ways.push_back(id);
                    if (cur.attrs.count("role")) {
                        rel.wayRoles.emplace_back(cur.attrs.find("role")->second);
                    } else {
                        rel.wayRoles.emplace_back("");
                    }
                }

            } else if (strcmp(cur.name, "tag") == 0) {
                if (keepAttrs.count(cur.attrs.find("k")->second))
                    rel.attrs[cur.attrs.find("k")->second] = cur.attrs.find("v")->second;
            }
        }
    } while (xml->next());

    // dont forget last relation
    uint64_t keepFlags = 0;
    uint64_t dropFlags = 0;
    if (rel.id && !rel.attrs.empty() &&
            (keepFlags = filter.keep(rel.attrs, OsmFilter::REL)) &&
            !(dropFlags = filter.drop(rel.attrs, OsmFilter::REL))) {
        rel.keepFlags = keepFlags;
        rel.dropFlags = dropFlags;
        return rel;
    }

    return OsmRel();
}

// _____________________________________________________________________________
void OsmBuilder::readRels(pfxml::file* xml, RelLst* rels, RelMap* nodeRels,
                          RelMap* wayRels, const OsmFilter& filter,
                          const AttrKeySet& keepAttrs,
                          Restrictions* rests) const {
    OsmRel rel;
    while ((rel = nextRel(xml, filter, keepAttrs)).id) {
        rels->rels.push_back(rel.attrs);
        if (rel.keepFlags & osm::REL_NO_DOWN) {
            rels->flat.insert(rels->rels.size() - 1);
        }
        for (osmid id : rel.nodes) (*nodeRels)[id].push_back(rels->rels.size() - 1);
        for (osmid id : rel.ways) (*wayRels)[id].push_back(rels->rels.size() - 1);

        // TODO(patrick): this is not needed for the filtering - remove it here!
        readRestr(rel, rests, filter);
    }
}

// _____________________________________________________________________________
void OsmBuilder::readRestr(const OsmRel& rel, Restrictions* rests,
                           const OsmFilter& filter) const {
    if (!rel.attrs.count("type")) return;
    if (rel.attrs.find("type")->second != "restriction") return;

    bool pos = filter.posRestr(rel.attrs);
    bool neg = filter.negRestr(rel.attrs);

    if (!pos && !neg) return;

    osmid from = 0;
    osmid to = 0;
    osmid via = 0;

    for (size_t i = 0; i < rel.ways.size(); i++) {
        if (rel.wayRoles[i] == "from") {
            if (from) return;  // only one from member supported
            from = rel.ways[i];
        }
        if (rel.wayRoles[i] == "to") {
            if (to) return;  // only one to member supported
            to = rel.ways[i];
        }
    }

    for (size_t i = 0; i < rel.nodes.size(); i++) {
        if (rel.nodeRoles[i] == "via") {
            via = rel.nodes[i];
            break;
        }
    }

    if (from && to && via) {
        if (pos)
            rests->pos[via].push_back(Restriction{from, to});
        else if (neg)
            rests->neg[via].push_back(Restriction{from, to});
    }
}

// _____________________________________________________________________________
std::string OsmBuilder::getAttrByFirstMatch(const DeepAttrLst& rule, osmid id,
                                            const AttrMap& attrs,
                                            const RelMap& entRels,
                                            const RelLst& rels,
                                            const Normalizer& norm) const {
    std::string ret;
    for (const auto& s : rule) {
        ret = norm(pfxml::file::decode(getAttr(s, id, attrs, entRels, rels)));
        if (!ret.empty()) return ret;
    }

    return ret;
}

// _____________________________________________________________________________
std::vector<std::string> OsmBuilder::getAttrMatchRanked(
        const DeepAttrLst& rule, osmid id, const AttrMap& attrs,
        const RelMap& entRels, const RelLst& rels, const Normalizer& norm) const {
    std::vector<std::string> ret;
    for (const auto& s : rule) {
        std::string tmp =
                norm(pfxml::file::decode(getAttr(s, id, attrs, entRels, rels)));
        if (!tmp.empty()) ret.push_back(tmp);
    }

    return ret;
}

// _____________________________________________________________________________
std::string OsmBuilder::getAttr(const DeepAttrRule& s, osmid id,
                                const AttrMap& attrs, const RelMap& entRels,
                                const RelLst& rels) const {
    if (s.relRule.kv.first.empty()) {
        if (attrs.find(s.attr) != attrs.end()) {
            return attrs.find(s.attr)->second;
        }
    } else {
        if (entRels.count(id)) {
            for (const auto& relId : entRels.find(id)->second) {
                if (OsmFilter::contained(rels.rels[relId], s.relRule.kv)) {
                    if (rels.rels[relId].count(s.attr)) {
                        return rels.rels[relId].find(s.attr)->second;
                    }
                }
            }
        }
    }
    return "";
}

// _____________________________________________________________________________
Nullable<StatInfo> OsmBuilder::getStatInfo(Node* node, osmid nid,
                                           const POINT& pos, const AttrMap& m,
                                           StAttrGroups* groups,
                                           const RelMap& nodeRels,
                                           const RelLst& rels,
                                           const OsmReadOpts& ops) const {
    std::string platform;
    std::vector<std::string> names;

    names = getAttrMatchRanked(ops.statAttrRules.nameRule, nid, m, nodeRels, rels,
                               ops.statNormzer);
    platform = getAttrByFirstMatch(ops.statAttrRules.platformRule, nid, m,
                                   nodeRels, rels, ops.trackNormzer);

    if (!names.size()) return Nullable<StatInfo>();

    auto ret = StatInfo(names[0], platform, true);

#ifdef PFAEDLE_STATION_IDS
    ret.setId(getAttrByFirstMatch(ops.statAttrRules.idRule, nid, m, nodeRels,
                                  rels, ops.idNormzer));
#endif

    for (size_t i = 1; i < names.size(); i++) ret.addAltName(names[i]);

    bool groupFound = false;

    for (const auto& rule : ops.statGroupNAttrRules) {
        if (groupFound) break;
        std::string ruleVal = getAttr(rule.attr, nid, m, nodeRels, rels);
        if (!ruleVal.empty()) {
            // check if a matching group exists
            for (auto* group : (*groups)[rule.attr.attr][ruleVal]) {
                if (groupFound) break;
                for (const auto* member : group->getNodes()) {
                    if (webMercMeterDist(*member->pl().getGeom(), pos) <= rule.maxDist) {
                        // ok, group is matching
                        groupFound = true;
                        if (node) group->addNode(node);
                        ret.setGroup(group);
                        break;
                    }
                }
            }
        }
    }

    if (!groupFound) {
        for (const auto& rule : ops.statGroupNAttrRules) {
            std::string ruleVal = getAttr(rule.attr, nid, m, nodeRels, rels);
            if (!ruleVal.empty()) {
                // add new group
                StatGroup* g = new StatGroup();
                if (node) g->addNode(node);
                ret.setGroup(g);
                (*groups)[rule.attr.attr][ruleVal].push_back(g);
                break;
            }
        }
    }

    return ret;
}

// _____________________________________________________________________________
double OsmBuilder::dist(const Node* a, const Node* b) {
    return webMercMeterDist(*a->pl().getGeom(), *b->pl().getGeom());
}

// _____________________________________________________________________________
double OsmBuilder::webMercDist(const Node* a, const Node* b) {
    return webMercMeterDist(*a->pl().getGeom(), *b->pl().getGeom());
}

// _____________________________________________________________________________
void OsmBuilder::writeGeoms(Graph* g) {
    for (auto* n : *g->getNds()) {
        for (auto* e : n->getAdjListOut()) {
            e->pl().addPoint(*e->getFrom()->pl().getGeom());
            e->pl().setLength(dist(e->getFrom(), e->getTo()));
            e->pl().addPoint(*e->getTo()->pl().getGeom());
        }
    }
}

// _____________________________________________________________________________
void OsmBuilder::fixGaps(Graph* g, NodeGrid* ng) {
    double METER = 1;
    for (auto* n : *g->getNds()) {
        if (n->getInDeg() + n->getOutDeg() == 1) {
            // get all nodes in distance
            std::set<Node*> ret;
            double distor = util::geo::webMercDistFactor(*n->pl().getGeom());
            ng->get(util::geo::pad(util::geo::getBoundingBox(*n->pl().getGeom()),
                                   METER / distor),
                    &ret);
            for (auto* nb : ret) {
                if (nb != n && (nb->getInDeg() + nb->getOutDeg()) == 1 &&
                        webMercDist(nb, n) <= METER / distor) {
                    // special case: both node are non-stations, move
                    // the end point nb to n and delete nb
                    if (!nb->pl().getSI() && !n->pl().getSI()) {
                        Node* otherN;
                        if (nb->getOutDeg())
                            otherN = (*nb->getAdjListOut().begin())->getOtherNd(nb);
                        else
                            otherN = (*nb->getAdjListIn().begin())->getOtherNd(nb);
                        LINE l;
                        l.push_back(*otherN->pl().getGeom());
                        l.push_back(*n->pl().getGeom());

                        Edge* e;
                        if (nb->getOutDeg())
                            e = g->addEdg(otherN, n, (*nb->getAdjListOut().begin())->pl());
                        else
                            e = g->addEdg(otherN, n, (*nb->getAdjListIn().begin())->pl());
                        if (e) {
                            *e->pl().getGeom() = l;
                            g->delNd(nb);
                            ng->remove(nb);
                        }
                    } else {
                        // if one of the nodes is a station, just add an edge between them
                        if (nb->getOutDeg())
                            g->addEdg(n, nb, (*nb->getAdjListOut().begin())->pl());
                        else
                            g->addEdg(n, nb, (*nb->getAdjListIn().begin())->pl());
                    }
                }
            }
        }
    }
}

// _____________________________________________________________________________
EdgeGrid OsmBuilder::buildEdgeIdx(Graph* g, size_t size,
                                  const BOX& webMercBox) {
    EdgeGrid ret(size, size, webMercBox, false);
    for (auto* n : *g->getNds()) {
        for (auto* e : n->getAdjListOut()) {
            assert(e->pl().getGeom());
            ret.add(*e->pl().getGeom(), e);
        }
    }
    return ret;
}

// _____________________________________________________________________________
NodeGrid OsmBuilder::buildNodeIdx(Graph* g, size_t size, const BOX& webMercBox,
                                  bool which) {
    NodeGrid ret(size, size, webMercBox, false);
    for (auto* n : *g->getNds()) {
        if (!which && n->getInDeg() + n->getOutDeg() == 1)
            ret.add(*n->pl().getGeom(), n);
        else if (which && n->pl().getSI())
            ret.add(*n->pl().getGeom(), n);
    }
    return ret;
}

// _____________________________________________________________________________
Node* OsmBuilder::depthSearch(const Edge* e, const StatInfo* si, const POINT& p,
                              double maxD, int maxFullTurns, double minAngle,
                              const SearchFunc& sfunc) {
    // shortcuts
    double dFrom = webMercMeterDist(*e->getFrom()->pl().getGeom(), p);
    double dTo = webMercMeterDist(*e->getTo()->pl().getGeom(), p);
    if (dFrom > maxD && dTo > maxD) return 0;

    if (dFrom <= maxD && sfunc(e->getFrom(), si)) return e->getFrom();
    if (dTo <= maxD && sfunc(e->getTo(), si)) return e->getTo();

    NodeCandPQ pq;
    NodeSet closed;
    pq.push(NodeCand{dFrom, e->getFrom(), e, 0});
    if (e->getFrom() != e->getTo()) pq.push(NodeCand{dTo, e->getTo(), e, 0});

    while (!pq.empty()) {
        auto cur = pq.top();
        pq.pop();
        if (closed.count(cur.node)) continue;
        closed.insert(cur.node);

        for (size_t i = 0; i < cur.node->getInDeg() + cur.node->getOutDeg(); i++) {
            trgraph::Node* cand;
            trgraph::Edge* edg;

            if (i < cur.node->getInDeg()) {
                edg = cur.node->getAdjListIn()[i];
                cand = edg->getFrom();
            } else {
                edg = cur.node->getAdjListOut()[i - cur.node->getInDeg()];
                cand = edg->getTo();
            }

            if (cand == cur.node) continue;  // dont follow self edges

            int fullTurn = 0;

            if (cur.fromEdge &&
                    cur.node->getInDeg() + cur.node->getOutDeg() >
                    2) {  // only intersection angles
                const POINT& toP = *cand->pl().getGeom();
                const POINT& fromP =
                        *cur.fromEdge->getOtherNd(cur.node)->pl().getGeom();
                const POINT& nodeP = *cur.node->pl().getGeom();

                if (util::geo::innerProd(nodeP, fromP, toP) < minAngle) fullTurn = 1;
            }

            if ((maxFullTurns < 0 || cur.fullTurns + fullTurn <= maxFullTurns) &&
                    cur.dist + edg->pl().getLength() < maxD && !closed.count(cand)) {
                if (sfunc(cand, si)) {
                    return cand;
                } else {
                    pq.push(NodeCand{cur.dist + edg->pl().getLength(), cand, edg,
                                     cur.fullTurns + fullTurn});
                }
            }
        }
    }

    return 0;
}

// _____________________________________________________________________________
bool OsmBuilder::isBlocked(const Edge* e, const StatInfo* si, const POINT& p,
                           double maxD, int maxFullTurns, double minAngle) {
    return depthSearch(e, si, p, maxD, maxFullTurns, minAngle, BlockSearch());
}

// _____________________________________________________________________________
Node* OsmBuilder::eqStatReach(const Edge* e, const StatInfo* si, const POINT& p,
                              double maxD, int maxFullTurns, double minAngle,
                              bool orphanSnap) {
    return depthSearch(e, si, p, maxD, maxFullTurns, minAngle,
                       EqSearch(orphanSnap));
}

// _____________________________________________________________________________
void OsmBuilder::getEdgCands(const POINT& geom, EdgeCandPQ* ret, EdgeGrid* eg,
                             double d) {
    double distor = util::geo::webMercDistFactor(geom);
    std::set<Edge*> neighs;
    BOX box = util::geo::pad(util::geo::getBoundingBox(geom), d / distor);
    eg->get(box, &neighs);

    for (auto* e : neighs) {
        double dist = util::geo::distToSegment(*e->getFrom()->pl().getGeom(),
                                               *e->getTo()->pl().getGeom(), geom);

        if (dist * distor <= d) {
            ret->push(EdgeCand(-dist, e));
        }
    }
}

// _____________________________________________________________________________
std::set<Node*> OsmBuilder::getMatchingNds(const NodePL& s, NodeGrid* ng,
                                           double d) {
    std::set<Node*> ret;
    double distor = util::geo::webMercDistFactor(*s.getGeom());
    std::set<Node*> neighs;
    BOX box = util::geo::pad(util::geo::getBoundingBox(*s.getGeom()), d / distor);
    ng->get(box, &neighs);

    for (auto* n : neighs) {
        if (n->pl().getSI() && n->pl().getSI()->simi(s.getSI()) > 0.5) {
            double dist = webMercMeterDist(*n->pl().getGeom(), *s.getGeom());
            if (dist < d) ret.insert(n);
        }
    }

    return ret;
}

// _____________________________________________________________________________
Node* OsmBuilder::getMatchingNd(const NodePL& s, NodeGrid* ng, double d) {
    double distor = util::geo::webMercDistFactor(*s.getGeom());
    std::set<Node*> neighs;
    BOX box = util::geo::pad(util::geo::getBoundingBox(*s.getGeom()), d / distor);
    ng->get(box, &neighs);

    Node* ret = 0;
    double bestD = std::numeric_limits<double>::max();

    for (auto* n : neighs) {
        if (n->pl().getSI() && n->pl().getSI()->simi(s.getSI()) > 0.5) {
            double dist = webMercMeterDist(*n->pl().getGeom(), *s.getGeom());
            if (dist < d && dist < bestD) {
                bestD = dist;
                ret = n;
            }
        }
    }

    return ret;
}

// _____________________________________________________________________________
std::set<Node*> OsmBuilder::snapStation(Graph* g, NodePL* s, EdgeGrid* eg,
                                        NodeGrid* sng, const OsmReadOpts& opts,
                                        Restrictor* restor, bool surrHeur,
                                        bool orphSnap, double d) {
    assert(s->getSI());
    std::set<Node*> ret;

    EdgeCandPQ pq;

    getEdgCands(*s->getGeom(), &pq, eg, d);

    if (pq.empty() && surrHeur) {
        // no station found in the first round, try again with the nearest
        // surrounding station with matching name
        const Node* best = getMatchingNd(*s, sng, opts.maxSnapFallbackHeurDistance);
        if (best) {
            getEdgCands(*best->pl().getGeom(), &pq, eg, d);
        } else {
            // if still no luck, get edge cands in fallback snap distance
            getEdgCands(*s->getGeom(), &pq, eg, opts.maxSnapFallbackHeurDistance);
        }
    }

    while (!pq.empty()) {
        auto* e = pq.top().second;
        pq.pop();
        auto geom =
                util::geo::projectOn(*e->getFrom()->pl().getGeom(), *s->getGeom(),
                                     *e->getTo()->pl().getGeom());

        Node* eq = 0;
        if (!(eq = eqStatReach(e, s->getSI(), geom, 2 * d, 0,
                               opts.maxAngleSnapReach, orphSnap))) {
            if (e->pl().lvl() > opts.maxSnapLevel) continue;
            if (isBlocked(e, s->getSI(), geom, opts.maxBlockDistance, 0,
                          opts.maxAngleSnapReach)) {
                continue;
            }

            // if the projected position is near (< 2 meters) the end point of this
            // way and the endpoint is not already a station, place the station there.
            if (!e->getFrom()->pl().getSI() &&
                    webMercMeterDist(geom, *e->getFrom()->pl().getGeom()) < 2) {
                e->getFrom()->pl().setSI(*s->getSI());
                if (s->getSI()->getGroup())
                    s->getSI()->getGroup()->addNode(e->getFrom());
                ret.insert(e->getFrom());
            } else if (!e->getTo()->pl().getSI() &&
                       webMercMeterDist(geom, *e->getTo()->pl().getGeom()) < 2) {
                e->getTo()->pl().setSI(*s->getSI());
                if (s->getSI()->getGroup()) s->getSI()->getGroup()->addNode(e->getTo());
                ret.insert(e->getTo());
            } else {
                s->setGeom(geom);
                Node* n = g->addNd(*s);
                if (n->pl().getSI()->getGroup())
                    n->pl().getSI()->getGroup()->addNode(n);
                sng->add(geom, n);

                auto ne = g->addEdg(e->getFrom(), n, e->pl());
                ne->pl().setLength(webMercDist(n, e->getFrom()));
                LINE l;
                l.push_back(*e->getFrom()->pl().getGeom());
                l.push_back(*n->pl().getGeom());
                *ne->pl().getGeom() = l;
                eg->add(l, ne);

                auto nf = g->addEdg(n, e->getTo(), e->pl());
                nf->pl().setLength(webMercDist(n, e->getTo()));
                LINE ll;
                ll.push_back(*n->pl().getGeom());
                ll.push_back(*e->getTo()->pl().getGeom());
                *nf->pl().getGeom() = ll;
                eg->add(l, nf);

                // replace edge in restrictor
                restor->replaceEdge(e, ne, nf);

                g->delEdg(e->getFrom(), e->getTo());
                eg->remove(e);
                ret.insert(n);
            }
        } else {
            // if the snapped station is very near to the original OSM station
            // write additional info from this snap station to the equivalent stat
            if (webMercMeterDist(*s->getGeom(), *eq->pl().getGeom()) <
                    opts.maxOsmStationDistance) {
                if (eq->pl().getSI()->getTrack().empty())
                    eq->pl().getSI()->setTrack(s->getSI()->getTrack());
            }
            ret.insert(eq);
        }
    }

    return ret;
}

// _____________________________________________________________________________
StatGroup* OsmBuilder::groupStats(const NodeSet& s) {
    if (!s.size()) return 0;
    // reference group
    StatGroup* ret = new StatGroup();
    bool used = false;

    for (auto* n : s) {
        if (!n->pl().getSI()) continue;
        used = true;
        if (n->pl().getSI()->getGroup()) {
            // this node is already in a group - merge this group with this one
            ret->merge(n->pl().getSI()->getGroup());
        } else {
            ret->addNode(n);
            n->pl().getSI()->setGroup(ret);
        }
    }

    if (!used) {
        delete ret;
        return 0;
    }

    return ret;
}

// _____________________________________________________________________________
std::vector<TransitEdgeLine*> OsmBuilder::getLines(
        const std::vector<size_t>& edgeRels, const RelLst& rels,
        const OsmReadOpts& ops) {
    std::vector<TransitEdgeLine*> ret;
    for (size_t relId : edgeRels) {
        TransitEdgeLine* elp = 0;

        if (_relLines.count(relId)) {
            elp = _relLines[relId];
        } else {
            TransitEdgeLine el;

            bool found = false;
            for (const auto& r : ops.relLinerules.sNameRule) {
                for (const auto& relAttr : rels.rels[relId]) {
                    if (relAttr.first == r) {
                        el.shortName = ops.lineNormzer(pfxml::file::decode(relAttr.second));
                        if (!el.shortName.empty()) found = true;
                    }
                }
                if (found) break;
            }

            found = false;
            for (const auto& r : ops.relLinerules.fromNameRule) {
                for (const auto& relAttr : rels.rels[relId]) {
                    if (relAttr.first == r) {
                        el.fromStr = ops.statNormzer(pfxml::file::decode(relAttr.second));
                        if (!el.fromStr.empty()) found = true;
                    }
                }
                if (found) break;
            }

            found = false;
            for (const auto& r : ops.relLinerules.toNameRule) {
                for (const auto& relAttr : rels.rels[relId]) {
                    if (relAttr.first == r) {
                        el.toStr = ops.statNormzer(pfxml::file::decode(relAttr.second));
                        if (!el.toStr.empty()) found = true;
                    }
                }
                if (found) break;
            }

            if (!el.shortName.size() && !el.fromStr.size() && !el.toStr.size())
                continue;

            if (_lines.count(el)) {
                elp = _lines[el];
                _relLines[relId] = elp;
            } else {
                elp = new TransitEdgeLine(el);
                _lines[el] = elp;
                _relLines[relId] = elp;
            }
        }
        ret.push_back(elp);
    }

    return ret;
}

// _____________________________________________________________________________
void OsmBuilder::getKeptAttrKeys(const OsmReadOpts& opts,
                                 AttrKeySet sets[3]) const {
    for (const auto& i : opts.statGroupNAttrRules) {
        if (i.attr.relRule.kv.first.empty()) {
            sets[0].insert(i.attr.attr);
        } else {
            sets[2].insert(i.attr.relRule.kv.first);
            sets[2].insert(i.attr.attr);
        }
    }

    for (const auto& i : opts.keepFilter) {
        for (size_t j = 0; j < 3; j++) sets[j].insert(i.first);
    }

    for (const auto& i : opts.dropFilter) {
        for (size_t j = 0; j < 3; j++) sets[j].insert(i.first);
    }

    for (const auto& i : opts.noHupFilter) {
        sets[0].insert(i.first);
    }

    for (const auto& i : opts.oneWayFilter) {
        sets[1].insert(i.first);
    }

    for (const auto& i : opts.oneWayFilterRev) {
        sets[1].insert(i.first);
    }

    for (const auto& i : opts.twoWayFilter) {
        sets[1].insert(i.first);
    }

    for (const auto& i : opts.stationFilter) {
        sets[0].insert(i.first);
        sets[2].insert(i.first);
    }

    for (const auto& i : opts.stationBlockerFilter) {
        sets[0].insert(i.first);
    }

    for (uint8_t j = 0; j < 7; j++) {
        for (const auto& kv : *(opts.levelFilters + j)) {
            sets[1].insert(kv.first);
        }
    }

    // restriction system
    for (const auto& i : opts.restrPosRestr) {
        sets[2].insert(i.first);
    }
    for (const auto& i : opts.restrNegRestr) {
        sets[2].insert(i.first);
    }
    for (const auto& i : opts.noRestrFilter) {
        sets[2].insert(i.first);
    }

    sets[2].insert("from");
    sets[2].insert("via");
    sets[2].insert("to");

    sets[2].insert(opts.relLinerules.toNameRule.begin(),
                   opts.relLinerules.toNameRule.end());
    sets[2].insert(opts.relLinerules.fromNameRule.begin(),
                   opts.relLinerules.fromNameRule.end());
    sets[2].insert(opts.relLinerules.sNameRule.begin(),
                   opts.relLinerules.sNameRule.end());

    for (const auto& i : opts.statAttrRules.nameRule) {
        if (i.relRule.kv.first.empty()) {
            sets[0].insert(i.attr);
        } else {
            sets[2].insert(i.relRule.kv.first);
            sets[2].insert(i.attr);
        }
    }

    for (const auto& i : opts.edgePlatformRules) {
        if (i.relRule.kv.first.empty()) {
            sets[1].insert(i.attr);
        } else {
            sets[2].insert(i.relRule.kv.first);
            sets[2].insert(i.attr);
        }
    }

    for (const auto& i : opts.statAttrRules.platformRule) {
        if (i.relRule.kv.first.empty()) {
            sets[0].insert(i.attr);
        } else {
            sets[2].insert(i.relRule.kv.first);
            sets[2].insert(i.attr);
        }
    }

    for (const auto& i : opts.statAttrRules.idRule) {
        if (i.relRule.kv.first.empty()) {
            sets[0].insert(i.attr);
        } else {
            sets[2].insert(i.relRule.kv.first);
            sets[2].insert(i.attr);
        }
    }
}

// _____________________________________________________________________________
void OsmBuilder::deleteOrphEdgs(Graph* g, const OsmReadOpts& opts) {
    size_t ROUNDS = 3;
    for (size_t c = 0; c < ROUNDS; c++) {
        for (auto i = g->getNds()->begin(); i != g->getNds()->end();) {
            if ((*i)->getInDeg() + (*i)->getOutDeg() != 1 || (*i)->pl().getSI()) {
                ++i;
                continue;
            }

            // check if the removal of this edge would transform a steep angle
            // full turn at an intersection into a node 2 eligible for contraction
            // if so, dont delete
            if (keepFullTurn(*i, opts.fullTurnAngle)) {
                ++i;
                continue;
            }

            i = g->delNd(*i);
            continue;
            i++;
        }
    }
}

// _____________________________________________________________________________
void OsmBuilder::deleteOrphNds(Graph* g) {
    for (auto i = g->getNds()->begin(); i != g->getNds()->end();) {
        if ((*i)->getInDeg() + (*i)->getOutDeg() == 0 &&
                !((*i)->pl().getSI() && (*i)->pl().getSI()->getGroup())) {
            i = g->delNd(i);
            // TODO(patrick): maybe delete from node grid?
        } else {
            i++;
        }
    }
}

// _____________________________________________________________________________
bool OsmBuilder::edgesSim(const Edge* a, const Edge* b) {
    if (static_cast<bool>(a->pl().oneWay()) ^ static_cast<bool>(b->pl().oneWay()))
        return false;
    if (a->pl().lvl() != b->pl().lvl()) return false;
    if (a->pl().getLines().size() != b->pl().getLines().size()) return false;
    if (a->pl().getLines() != b->pl().getLines()) return false;
    if (a->pl().oneWay() && b->pl().oneWay()) {
        if (a->getFrom() != b->getTo() && a->getTo() != b->getFrom()) return false;
    }
    if (a->pl().isRestricted() || b->pl().isRestricted()) return false;

    return true;
}

// _____________________________________________________________________________
const EdgePL& OsmBuilder::mergeEdgePL(Edge* a, Edge* b) {
    const Node* n = 0;
    if (a->getFrom() == b->getFrom())
        n = a->getFrom();
    else if (a->getFrom() == b->getTo())
        n = a->getFrom();
    else
        n = a->getTo();

    if (a->getTo() == n && b->getTo() == n) {
        // --> n <--
        a->pl().getGeom()->insert(a->pl().getGeom()->end(),
                                  b->pl().getGeom()->rbegin(),
                                  b->pl().getGeom()->rend());
    } else if (a->getTo() == n && b->getFrom() == n) {
        // --> n -->
        a->pl().getGeom()->insert(a->pl().getGeom()->end(),
                                  b->pl().getGeom()->begin(),
                                  b->pl().getGeom()->end());
    } else if (a->getFrom() == n && b->getTo() == n) {
        // <-- n <--
        std::reverse(a->pl().getGeom()->begin(), a->pl().getGeom()->end());
        a->pl().getGeom()->insert(a->pl().getGeom()->end(),
                                  b->pl().getGeom()->rbegin(),
                                  b->pl().getGeom()->rend());
    } else {
        // <-- n -->
        std::reverse(a->pl().getGeom()->begin(), a->pl().getGeom()->end());
        a->pl().getGeom()->insert(a->pl().getGeom()->end(),
                                  b->pl().getGeom()->begin(),
                                  b->pl().getGeom()->end());
    }

    a->pl().setLength(a->pl().getLength() + b->pl().getLength());

    return a->pl();
}

// _____________________________________________________________________________
void OsmBuilder::collapseEdges(Graph* g) {
    for (auto* n : *g->getNds()) {
        if (n->getOutDeg() + n->getInDeg() != 2 || n->pl().getSI()) continue;

        Edge* ea;
        Edge* eb;
        if (n->getOutDeg() == 2) {
            ea = *n->getAdjListOut().begin();
            eb = *n->getAdjListOut().rbegin();
        } else if (n->getInDeg() == 2) {
            ea = *n->getAdjListIn().begin();
            eb = *n->getAdjListIn().rbegin();
        } else {
            ea = *n->getAdjListOut().begin();
            eb = *n->getAdjListIn().begin();
        }

        // important, we don't have a multigraph! if the same edge
        // will already exist, leave this node
        if (g->getEdg(ea->getOtherNd(n), eb->getOtherNd(n))) continue;
        if (g->getEdg(eb->getOtherNd(n), ea->getOtherNd(n))) continue;

        if (edgesSim(ea, eb)) {
            if (ea->pl().oneWay() && ea->getOtherNd(n) != ea->getFrom()) {
                g->addEdg(eb->getOtherNd(n), ea->getOtherNd(n), mergeEdgePL(eb, ea));
            } else {
                g->addEdg(ea->getOtherNd(n), eb->getOtherNd(n), mergeEdgePL(ea, eb));
            }
            g->delEdg(ea->getFrom(), ea->getTo());
            g->delEdg(eb->getFrom(), eb->getTo());
        }
    }
}

// _____________________________________________________________________________
void OsmBuilder::simplifyGeoms(Graph* g) {
    for (auto* n : *g->getNds()) {
        for (auto* e : n->getAdjListOut()) {
            (*e->pl().getGeom()) = util::geo::simplify(*e->pl().getGeom(), 0.5);
        }
    }
}

// _____________________________________________________________________________
uint32_t OsmBuilder::writeComps(Graph* g) {
    Component* comp = new Component{7};
    uint32_t numC = 0;

    for (auto* n : *g->getNds()) {
        if (!n->pl().getComp()) {
            std::stack<std::pair<Node*, Edge*>> q;
            q.push(std::pair<Node*, Edge*>(n, 0));
            while (!q.empty()) {
                std::pair<Node*, Edge*> cur = q.top();
                q.pop();

                cur.first->pl().setComp(comp);
                for (auto* e : cur.first->getAdjListOut()) {
                    if (e->pl().lvl() < comp->minEdgeLvl)
                        comp->minEdgeLvl = e->pl().lvl();
                    if (!e->getOtherNd(cur.first)->pl().getComp())
                        q.push(std::pair<Node*, Edge*>(e->getOtherNd(cur.first), e));
                }
                for (auto* e : cur.first->getAdjListIn()) {
                    if (e->pl().lvl() < comp->minEdgeLvl)
                        comp->minEdgeLvl = e->pl().lvl();
                    if (!e->getOtherNd(cur.first)->pl().getComp())
                        q.push(std::pair<Node*, Edge*>(e->getOtherNd(cur.first), e));
                }
            }

            numC++;
            comp = new Component{7};
        }
    }

    // the last comp was not used
    delete comp;

    return numC;
}

// _____________________________________________________________________________
void OsmBuilder::writeEdgeTracks(const EdgTracks& tracks) {
    for (const auto& tr : tracks) {
        if (tr.first->getTo()->pl().getSI() &&
                tr.first->getTo()->pl().getSI()->getTrack().empty()) {
            tr.first->getTo()->pl().getSI()->setTrack(tr.second);
        }
        if (tr.first->getFrom()->pl().getSI() &&
                tr.first->getFrom()->pl().getSI()->getTrack().empty()) {
            tr.first->getFrom()->pl().getSI()->setTrack(tr.second);
        }
    }
}

// _____________________________________________________________________________
void OsmBuilder::writeODirEdgs(Graph* g, Restrictor* restor) {
    for (auto* n : *g->getNds()) {
        for (auto* e : n->getAdjListOut()) {
            if (g->getEdg(e->getTo(), e->getFrom())) continue;
            auto newE = g->addEdg(e->getTo(), e->getFrom(), e->pl().revCopy());
            if (e->pl().isRestricted()) restor->duplicateEdge(e, newE);
        }
    }
}

// _____________________________________________________________________________
void OsmBuilder::writeSelfEdgs(Graph* g) {
    for (auto* n : *g->getNds()) {
        if (n->pl().getSI() && n->getAdjListOut().size() == 0) {
            g->addEdg(n, n);
        }
    }
}

// _____________________________________________________________________________
bool OsmBuilder::keepFullTurn(const trgraph::Node* n, double ang) {
    if (n->getInDeg() + n->getOutDeg() != 1) return false;

    const trgraph::Edge* e = 0;
    if (n->getOutDeg())
        e = n->getAdjListOut().front();
    else
        e = n->getAdjListIn().front();

    auto other = e->getOtherNd(n);

    if (other->getInDeg() + other->getOutDeg() == 3) {
        const trgraph::Edge* a = 0;
        const trgraph::Edge* b = 0;
        for (auto f : other->getAdjListIn()) {
            if (f != e && !a)
                a = f;
            else if (f != e && !b)
                b = f;
        }

        for (auto f : other->getAdjListOut()) {
            if (f != e && !a)
                a = f;
            else if (f != e && !b)
                b = f;
        }

        auto ap = a->pl().backHop();
        auto bp = b->pl().backHop();
        if (a->getTo() != other) ap = a->pl().frontHop();
        if (b->getTo() != other) bp = b->pl().frontHop();

        return router::angSmaller(ap, *other->pl().getGeom(), bp, ang);
    }

    return false;
}

// _____________________________________________________________________________
void OsmBuilder::snapStats(const OsmReadOpts& opts, Graph* g,
                           const BBoxIdx& bbox, size_t gridSize,
                           router::FeedStops* fs, Restrictor* res,
                           const NodeSet& orphanStations) {
    NodeGrid sng = buildNodeIdx(g, gridSize, bbox.getFullWebMercBox(), true);
    EdgeGrid eg = buildEdgeIdx(g, gridSize, bbox.getFullWebMercBox());

    LOG(DEBUG) << "Grid size of " << sng.getXWidth() << "x" << sng.getYHeight();

    for (double d : opts.maxSnapDistances) {
        for (auto s : orphanStations) {
            POINT geom = *s->pl().getGeom();
            NodePL pl = s->pl();
            pl.getSI()->setIsFromOsm(false);
            const auto& r =
                    snapStation(g, &pl, &eg, &sng, opts, res, false, false, d);
            groupStats(r);
            for (auto n : r) {
                // if the snapped station is very near to the original OSM
                // station, set is-from-osm to true
                if (webMercMeterDist(geom, *n->pl().getGeom()) <
                        opts.maxOsmStationDistance) {
                    if (n->pl().getSI()) n->pl().getSI()->setIsFromOsm(true);
                }
            }
        }
    }

    if (!fs) return;

    std::vector<const Stop*> notSnapped;

    for (auto& s : *fs) {
        bool snapped = false;
        auto pl = plFromGtfs(s.first, opts);
        for (size_t i = 0; i < opts.maxSnapDistances.size(); i++) {
            double d = opts.maxSnapDistances[i];

            StatGroup* group = groupStats(
                        snapStation(g, &pl, &eg, &sng, opts, res,
                                    i == opts.maxSnapDistances.size() - 1, false, d));

            if (group) {
                group->addStop(s.first);
                (*fs)[s.first] = *group->getNodes().begin();
                snapped = true;
            }
        }
        if (!snapped) {
            LOG(VDEBUG) << "Could not snap station "
                        << "(" << pl.getSI()->getName() << ")"
                        << " (" << s.first->getLat() << "," << s.first->getLng()
                        << ") in normal run, trying again later in orphan mode.";
            if (!bbox.contains(*pl.getGeom())) {
                LOG(VDEBUG) << "Note: '" << pl.getSI()->getName()
                            << "' does not lie within the bounds for this graph and "
                               "may be a stray station";
            }
            notSnapped.push_back(s.first);
        }
    }

    if (notSnapped.size())
        LOG(VDEBUG) << notSnapped.size() << " stations could not be snapped in "
                                            "normal run, trying again in orphan "
                                            "mode.";

    // try again, but aggressively snap to orphan OSM stations which have
    // not been assigned to any GTFS stop yet
    for (auto& s : notSnapped) {
        bool snapped = false;
        auto pl = plFromGtfs(s, opts);
        for (size_t i = 0; i < opts.maxSnapDistances.size(); i++) {
            double d = opts.maxSnapDistances[i];

            StatGroup* group = groupStats(
                        snapStation(g, &pl, &eg, &sng, opts, res,
                                    i == opts.maxSnapDistances.size() - 1, true, d));

            if (group) {
                group->addStop(s);
                // add the added station name as an alt name to ensure future
                // similarity
                for (auto n : group->getNodes()) {
                    if (n->pl().getSI())
                        n->pl().getSI()->addAltName(pl.getSI()->getName());
                }
                (*fs)[s] = *group->getNodes().begin();
                snapped = true;
            }
        }
        if (!snapped) {
            // finally give up

            // add a group with only this stop in it
            StatGroup* dummyGroup = new StatGroup();
            Node* dummyNode = g->addNd(pl);

            dummyNode->pl().getSI()->setGroup(dummyGroup);
            dummyGroup->addNode(dummyNode);
            dummyGroup->addStop(s);
            (*fs)[s] = dummyNode;
            if (!bbox.contains(*pl.getGeom())) {
                LOG(VDEBUG) << "Could not snap station "
                            << "(" << pl.getSI()->getName() << ")"
                            << " (" << s->getLat() << "," << s->getLng() << ")";
                LOG(VDEBUG) << "Note: '" << pl.getSI()->getName()
                            << "' does not lie within the bounds for this graph and "
                               "may be a stray station";
            } else {
                // only warn if it is contained in the BBOX for this graph
                LOG(WARN) << "Could not snap station "
                          << "(" << pl.getSI()->getName() << ")"
                          << " (" << s->getLat() << "," << s->getLng() << ")";
            }
        }
    }
}

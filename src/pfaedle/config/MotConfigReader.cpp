// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <set>
#include <string>
#include "pfaedle/config/MotConfigReader.h"
#include "util/Misc.h"
#include "util/String.h"
#include "util/log/Log.h"

using pfaedle::config::MotConfigReader;
using pfaedle::config::MotConfig;
using pfaedle::osm::FilterRule;
using pfaedle::osm::KeyVal;
using configparser::ConfigFileParser;
using configparser::ParseExc;
using pfaedle::osm::DeepAttrRule;
using pfaedle::trgraph::ReplRules;
using ad::cppgtfs::gtfs::Route;

// _____________________________________________________________________________
MotConfigReader::MotConfigReader() {}

// _____________________________________________________________________________
void MotConfigReader::parse(const std::vector<std::string>& paths) {
  ConfigFileParser p;

  // parse explicitely given paths
  for (const auto& s : paths) {
    LOG(DEBUG) << "Reading config file " << s;
    p.parse(s);
  }

  for (const auto& sec : p.getSecs()) {
    MotConfig curCfg;
    std::string secStr = sec.first;
    if (secStr.empty()) continue;
    std::set<std::string> procedKeys;

    if (p.hasKey(secStr, "osm_filter_keep")) {
      procedKeys.insert("osm_filter_keep");
      for (const auto& kvs : p.getStrArr(sec.first, "osm_filter_keep", ' ')) {
        auto fRule = getFRule(kvs);
        curCfg.osmBuildOpts.keepFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    for (uint8_t i = 0; i < 8; i++) {
      std::string name = std::string("osm_filter_lvl") + std::to_string(i);
      if (p.hasKey(secStr, name)) {
        procedKeys.insert(name);
        for (const auto& kvs : p.getStrArr(sec.first, name, ' ')) {
          auto fRule = getFRule(kvs);
          curCfg.osmBuildOpts.levelFilters[i][fRule.kv.first].insert(
              osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
        }
      }
    }

    if (p.hasKey(secStr, "osm_filter_drop")) {
      procedKeys.insert("osm_filter_drop");
      for (const auto& kvs : p.getStrArr(sec.first, "osm_filter_drop", ' ')) {
        auto fRule = getFRule(kvs);
        curCfg.osmBuildOpts.dropFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_max_snap_level")) {
      procedKeys.insert("osm_max_snap_level");
      curCfg.osmBuildOpts.maxSnapLevel =
          p.getInt(sec.first, "osm_max_snap_level");
    } else {
      curCfg.osmBuildOpts.maxSnapLevel = 7;
    }

    if (p.hasKey(secStr, "osm_filter_nohup")) {
      procedKeys.insert("osm_filter_nohup");
      for (const auto& kvs : p.getStrArr(sec.first, "osm_filter_nohup", ' ')) {
        auto fRule = getFRule(kvs);
        curCfg.osmBuildOpts.noHupFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_oneway")) {
      procedKeys.insert("osm_filter_oneway");
      for (const auto& kvs : p.getStrArr(sec.first, "osm_filter_oneway", ' ')) {
        auto fRule = getFRule(kvs);
        curCfg.osmBuildOpts.oneWayFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_oneway_reverse")) {
      procedKeys.insert("osm_filter_oneway_reverse");
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_filter_oneway_reverse", ' ')) {
        auto fRule = getFRule(kvs);
        curCfg.osmBuildOpts.oneWayFilterRev[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_undirected")) {
      procedKeys.insert("osm_filter_undirected");
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_filter_undirected", ' ')) {
        auto fRule = getFRule(kvs);
        curCfg.osmBuildOpts.twoWayFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_station")) {
      procedKeys.insert("osm_filter_station");
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_filter_station", ' ')) {
        auto fRule = getFRule(kvs);
        curCfg.osmBuildOpts.stationFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_station_blocker")) {
      procedKeys.insert("osm_filter_station_blocker");
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_filter_station_blocker", ' ')) {
        auto fRule = getFRule(kvs);
        curCfg.osmBuildOpts.stationBlockerFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_node_positive_restriction")) {
      procedKeys.insert("osm_node_positive_restriction");
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_node_positive_restriction", ' ')) {
        auto fRule = getFRule(kvs);
        curCfg.osmBuildOpts.restrPosRestr[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_node_negative_restriction")) {
      procedKeys.insert("osm_node_negative_restriction");
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_node_negative_restriction", ' ')) {
        auto fRule = getFRule(kvs);
        curCfg.osmBuildOpts.restrNegRestr[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_no_restriction")) {
      procedKeys.insert("osm_filter_no_restriction");
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_filter_no_restriction", ' ')) {
        auto fRule = getFRule(kvs);
        curCfg.osmBuildOpts.noRestrFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_station_name_attrs")) {
      procedKeys.insert("osm_station_name_attrs");
      for (const std::string& r :
           p.getStrArr(sec.first, "osm_station_name_attrs", ' ')) {
        curCfg.osmBuildOpts.statAttrRules.nameRule.push_back(
            getDeepAttrRule(r));
      }
    }

    if (p.hasKey(secStr, "osm_track_number_tags")) {
      procedKeys.insert("osm_track_number_tags");
      for (const std::string& r :
           p.getStrArr(sec.first, "osm_track_number_tags", ' ')) {
        curCfg.osmBuildOpts.statAttrRules.platformRule.push_back(
            getDeepAttrRule(r));
      }
    }

    if (p.hasKey(secStr, "osm_station_id_attrs")) {
      procedKeys.insert("osm_station_id_attrs");
      for (const std::string& r :
           p.getStrArr(sec.first, "osm_station_id_attrs", ' ')) {
        curCfg.osmBuildOpts.statAttrRules.idRule.push_back(getDeepAttrRule(r));
      }
    }

    if (p.hasKey(secStr, "osm_edge_track_number_tags")) {
      procedKeys.insert("osm_edge_track_number_tags");
      for (const std::string& r :
           p.getStrArr(sec.first, "osm_edge_track_number_tags", ' ')) {
        curCfg.osmBuildOpts.edgePlatformRules.push_back(getDeepAttrRule(r));
      }
    }

    if (p.hasKey(secStr, "osm_station_group_attrs")) {
      procedKeys.insert("osm_station_group_attrs");
      auto arr = p.getStrArr(secStr, "osm_station_group_attrs", ' ');

      for (const auto& ruleStr : arr) {
        auto deep = getDeepAttrRule(ruleStr);
        // TODO(patrick): getKv is misused here as a a=b parser
        auto attrD = getKv(deep.attr);
        deep.attr = attrD.first;
        double dist = atof(attrD.second.c_str());
        curCfg.osmBuildOpts.statGroupNAttrRules.push_back({deep, dist});
      }
    }

    if (p.hasKey(secStr, "osm_line_relation_tags")) {
      procedKeys.insert("osm_line_relation_tags");
      auto arr = p.getStrArr(secStr, "osm_line_relation_tags", ' ');

      for (const auto& ruleStr : arr) {
        auto rule = getKv(ruleStr);
        auto tags = util::split(rule.second, ',');
        if (rule.first == "from_name")
          curCfg.osmBuildOpts.relLinerules.fromNameRule = tags;
        else if (rule.first == "to_name")
          curCfg.osmBuildOpts.relLinerules.toNameRule = tags;
        else if (rule.first == "line_name")
          curCfg.osmBuildOpts.relLinerules.sNameRule = tags;
      }
    }

    if (p.hasKey(secStr, "osm_max_snap_distance")) {
      procedKeys.insert("osm_max_snap_distance");
      curCfg.osmBuildOpts.maxSnapDistances =
          p.getDoubleArr(secStr, "osm_max_snap_distance", ',');
    } else {
      curCfg.osmBuildOpts.maxSnapDistances.push_back(50);
    }

    if (p.hasKey(secStr, "osm_max_snap_fallback_distance")) {
      procedKeys.insert("osm_max_snap_fallback_distance");
      curCfg.osmBuildOpts.maxSnapFallbackHeurDistance =
          p.getDouble(secStr, "osm_max_snap_fallback_distance");
    } else {
      curCfg.osmBuildOpts.maxSnapFallbackHeurDistance =
          *std::max_element(curCfg.osmBuildOpts.maxSnapDistances.begin(),
                            curCfg.osmBuildOpts.maxSnapDistances.end()) *
          2;
    }

    if (p.hasKey(secStr, "osm_max_osm_station_distance")) {
      procedKeys.insert("osm_max_osm_station_distance");
      curCfg.osmBuildOpts.maxOsmStationDistance =
          p.getDouble(secStr, "osm_max_osm_station_distance");
    } else {
      curCfg.osmBuildOpts.maxOsmStationDistance = 5;
    }

    if (p.hasKey(secStr, "osm_max_node_block_distance")) {
      procedKeys.insert("osm_max_node_block_distance");
      curCfg.osmBuildOpts.maxBlockDistance =
          p.getDouble(secStr, "osm_max_node_block_distance");
    } else {
      curCfg.osmBuildOpts.maxBlockDistance =
          *std::max_element(curCfg.osmBuildOpts.maxSnapDistances.begin(),
                            curCfg.osmBuildOpts.maxSnapDistances.end()) /
          8;
    }

    for (uint8_t i = 0; i < 8; i++) {
      std::string name =
          std::string("routing_lvl") + std::to_string(i) + "_fac";
      if (p.hasKey(secStr, name)) {
        procedKeys.insert(name);
        double v = p.getDouble(sec.first, name);
        curCfg.routingOpts.levelPunish[i] = v;
      } else {
        curCfg.routingOpts.levelPunish[i] = 1;
      }
    }

    if (p.hasKey(secStr, "routing_full_turn_punish")) {
      procedKeys.insert("routing_full_turn_punish");
      curCfg.routingOpts.fullTurnPunishFac =
          p.getDouble(secStr, "routing_full_turn_punish");
    }

    if (p.hasKey(secStr, "routing_no_self_hops")) {
      procedKeys.insert("routing_no_self_hops");
      curCfg.routingOpts.noSelfHops = p.getBool(secStr, "routing_no_self_hops");
    }

    if (p.hasKey(secStr, "routing_full_turn_angle")) {
      procedKeys.insert("routing_full_turn_angle");
      double ang = p.getDouble(secStr, "routing_full_turn_angle");
      curCfg.routingOpts.fullTurnAngle = ang;
      curCfg.osmBuildOpts.fullTurnAngle = ang;
    } else {
      curCfg.routingOpts.fullTurnAngle = 5;
      curCfg.osmBuildOpts.fullTurnAngle = 5;
    }

    if (p.hasKey(secStr, "routing_snap_full_turn_angle")) {
      procedKeys.insert("routing_snap_full_turn_angle");
      double ang = p.getDouble(secStr, "routing_snap_full_turn_angle");
      curCfg.osmBuildOpts.maxAngleSnapReach = ang;
    } else {
      curCfg.osmBuildOpts.maxAngleSnapReach = curCfg.routingOpts.fullTurnAngle;
    }

    if (p.hasKey(secStr, "routing_pass_thru_station_punish")) {
      procedKeys.insert("routing_pass_thru_station_punish");
      curCfg.routingOpts.passThruStationsPunish =
          p.getDouble(secStr, "routing_pass_thru_station_punish");
    }

    if (p.hasKey(secStr, "routing_one_way_meter_punish_fac")) {
      procedKeys.insert("routing_one_way_meter_punish_fac");
      curCfg.routingOpts.oneWayPunishFac =
          p.getDouble(secStr, "routing_one_way_meter_punish_fac");
    }

    if (p.hasKey(secStr, "routing_one_way_edge_punish")) {
      procedKeys.insert("routing_one_way_edge_punish");
      curCfg.routingOpts.oneWayEdgePunish =
          p.getDouble(secStr, "routing_one_way_edge_punish");
    }

    if (p.hasKey(secStr, "routing_line_unmatched_punish_fac")) {
      procedKeys.insert("routing_line_unmatched_punish_fac");
      curCfg.routingOpts.lineUnmatchedPunishFact =
          p.getDouble(secStr, "routing_line_unmatched_punish_fac");
    }

    if (p.hasKey(secStr, "routing_platform_unmatched_punish")) {
      procedKeys.insert("routing_platform_unmatched_punish");
      curCfg.routingOpts.platformUnmatchedPen =
          p.getDouble(secStr, "routing_platform_unmatched_punish");
    }

    if (p.hasKey(secStr, "routing_non_osm_station_punish")) {
      procedKeys.insert("routing_non_osm_station_punish");
      curCfg.routingOpts.nonOsmPen =
          p.getDouble(secStr, "routing_non_osm_station_punish");
    } else {
      curCfg.routingOpts.nonOsmPen = 0;
    }

    if (p.hasKey(secStr, "routing_station_distance_punish_fac")) {
      procedKeys.insert("routing_station_distance_punish_fac");
      curCfg.routingOpts.stationDistPenFactor =
          p.getDouble(secStr, "routing_station_distance_punish_fac");
    } else {
      curCfg.routingOpts.stationDistPenFactor = 1;
    }

    if (p.hasKey(secStr, "station_normalize_chain")) {
      procedKeys.insert("station_normalize_chain");
      try {
        auto arr = p.getStrArr(secStr, "station_normalize_chain", ';');
        curCfg.osmBuildOpts.statNormzer =
            trgraph::Normalizer(getNormRules(arr));
      } catch (const std::exception& e) {
        throw ParseExc(p.getVal(secStr, "station_normalize_chain").line,
                       p.getVal(secStr, "station_normalize_chain").pos,
                       "<valid regular expression>",
                       std::string("<regex error: ") + e.what() + ">",
                       p.getVal(secStr, "station_normalize_chain").file);
      }
    }

    if (p.hasKey(secStr, "track_normalize_chain")) {
      procedKeys.insert("track_normalize_chain");
      try {
        auto arr = p.getStrArr(secStr, "track_normalize_chain", ';');
        curCfg.osmBuildOpts.trackNormzer =
            trgraph::Normalizer(getNormRules(arr));
      } catch (const std::exception& e) {
        throw ParseExc(p.getVal(secStr, "track_normalize_chain").line,
                       p.getVal(secStr, "track_normalize_chain").pos,
                       "<valid regular expression>",
                       std::string("<regex error: ") + e.what() + ">",
                       p.getVal(secStr, "track_normalize_chain").file);
      }
    }

    if (p.hasKey(secStr, "line_normalize_chain")) {
      procedKeys.insert("line_normalize_chain");
      try {
        auto arr = p.getStrArr(secStr, "line_normalize_chain", ';');
        curCfg.osmBuildOpts.lineNormzer =
            trgraph::Normalizer(getNormRules(arr));
      } catch (const std::exception& e) {
        throw ParseExc(p.getVal(secStr, "line_normalize_chain").line,
                       p.getVal(secStr, "line_normalize_chain").pos,
                       "<valid regular expression>",
                       std::string("<regex error: ") + e.what() + ">",
                       p.getVal(secStr, "line_normalize_chain").file);
      }
    }

    if (p.hasKey(secStr, "station_id_normalize_chain")) {
      procedKeys.insert("station_id_normalize_chain");
      try {
        auto arr = p.getStrArr(secStr, "station_id_normalize_chain", ';');
        curCfg.osmBuildOpts.idNormzer = trgraph::Normalizer(getNormRules(arr));
      } catch (const std::exception& e) {
        throw ParseExc(p.getVal(secStr, "station_id_normalize_chain").line,
                       p.getVal(secStr, "station_id_normalize_chain").pos,
                       "<valid regular expression>",
                       std::string("<regex error: ") + e.what() + ">",
                       p.getVal(secStr, "station_normalize_chain").file);
      }
    }

    for (const auto& kv : p.getKeyVals(secStr)) {
      if (!procedKeys.count(kv.first))
        curCfg.unproced[kv.first] = kv.second.val;
    }

    bool found = false;

    for (auto& cfg : _cfgs) {
      if (cfg == curCfg) {
        for (auto mot :
             ad::cppgtfs::gtfs::flat::Route::getTypesFromString(secStr)) {
          cfg.mots.insert(mot);
        }
        found = true;
        break;
      }
    }

    if (!found) {
      curCfg.mots = ad::cppgtfs::gtfs::flat::Route::getTypesFromString(secStr);
      _cfgs.push_back(curCfg);
    }
  }
}

// _____________________________________________________________________________
ReplRules MotConfigReader::getNormRules(
    const std::vector<std::string>& arr) const {
  trgraph::ReplRules ret;

  for (auto a : arr) {
    size_t p = a.find(" -> ");
    if (p == std::string::npos) continue;

    trgraph::ReplRule r;
    r.first = a.substr(0, p);
    r.second = a.substr(p + 4, std::string::npos);

    if (r.first.size() > 1 && r.first.front() == '\'' &&
        r.first.back() == '\'') {
      r.first = r.first.substr(1, r.first.size() - 2);
    }

    if (r.second.size() > 1 && r.second.front() == '\'' &&
        r.second.back() == '\'') {
      r.second = r.second.substr(1, r.second.size() - 2);
    }

    ret.push_back(r);
  }

  return ret;
}

// _____________________________________________________________________________
uint64_t MotConfigReader::getFlags(const std::set<string>& flags) const {
  uint64_t ret = osm::USE;

  for (const auto& flag : flags) {
    if (flag == "rel_flat") {
      ret |= osm::REL_NO_DOWN;
      continue;
    }
    if (flag == "no_match_nds") {
      ret |= osm::NO_NODES;
      continue;
    }
    if (flag == "no_match_rels") {
      ret |= osm::NO_RELATIONS;
      continue;
    }
    if (flag == "no_match_ways") {
      ret |= osm::NO_WAYS;
      continue;
    }
    if (flag == "mult_val_match") {
      ret |= osm::MULT_VAL_MATCH;
      continue;
    }
  }

  return ret;
}

// _____________________________________________________________________________
FilterRule MotConfigReader::getFRule(const std::string& r) const {
  osm::FilterRule ret;

  auto parts = util::split(util::trim(r), '|');

  ret.kv = getKv(parts[0]);
  ret.flags = std::set<std::string>(parts.begin() + 1, parts.end());

  return ret;
}

// _____________________________________________________________________________
KeyVal MotConfigReader::getKv(const std::string& kv) const {
  osm::KeyVal ret;
  size_t p = kv.find('=', 0);
  ret.first = kv.substr(0, p);

  if (p != std::string::npos) {
    ret.second = kv.substr(p + 1, std::string::npos);
  }

  return ret;
}

// _____________________________________________________________________________
const std::vector<MotConfig>& MotConfigReader::getConfigs() const {
  return _cfgs;
}

// _____________________________________________________________________________
DeepAttrRule MotConfigReader::getDeepAttrRule(const std::string& rule) const {
  if (rule[0] == '[' && rule.find(']') != std::string::npos) {
    auto kv = getFRule(rule.substr(1, rule.find(']') - 1));
    std::string attr = rule.substr(rule.find(']') + 1);
    return osm::DeepAttrRule{attr, kv};
  } else {
    return osm::DeepAttrRule{rule, osm::FilterRule()};
  }
}

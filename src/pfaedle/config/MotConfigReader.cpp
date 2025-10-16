// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <limits>
#include <set>
#include <string>
#include "pfaedle/config/MotConfigReader.h"
#include "pfaedle/osm/OsmReadOpts.h"
#include "util/Misc.h"
#include "util/String.h"
#include "util/log/Log.h"

using ad::cppgtfs::gtfs::Route;
using configparser::ConfigFileParser;
using configparser::ParseExc;
using pfaedle::config::MotConfig;
using pfaedle::config::MotConfigReader;
using pfaedle::osm::DeepAttrRule;
using pfaedle::osm::FilterRule;
using pfaedle::osm::KeyVal;
using pfaedle::trgraph::ReplRules;
using util::WARN;
using util::INFO;
using util::ERROR;
using util::DEBUG;
using util::VDEBUG;

double DEF_TRANS_PEN = 0.0083;

// _____________________________________________________________________________
MotConfigReader::MotConfigReader() {}

// _____________________________________________________________________________
void MotConfigReader::parse(const std::vector<std::string>& paths,
                            const std::string& literal) {
  ConfigFileParser p;

  // parse explicitely given paths
  for (const auto& s : paths) {
    LOG(DEBUG) << "Reading config file " << s;
    p.parse(s);
  }

  if (literal.size()) p.parseStr(literal);

  for (const auto& sec : p.getSecs()) {
    MotConfig cfg;

    cfg.transWeight = "expo";

    std::string secStr = sec.first;
    if (secStr.empty()) continue;

    if (p.hasKey(secStr, "routing_emission_method")) {
      cfg.routingOpts.emPenMethod =
          p.getStr(secStr, "routing_emission_method");
    } else {
      cfg.routingOpts.emPenMethod = "exp";
    }

    if (p.hasKey(secStr, "routing_transition_method")) {
      cfg.routingOpts.transPenMethod =
          p.getStr(secStr, "routing_transition_method");
    } else {
      cfg.routingOpts.transPenMethod = "exp";
    }

    if (p.hasKey(secStr, "station_similarity_classification_method")) {
      cfg.routingOpts.statsimiMethod =
          p.getStr(secStr, "station_similarity_classification_method");
    } else {
      cfg.routingOpts.statsimiMethod = "jaccard-geodist";
    }

    if (p.hasKey(secStr, "routing_use_stations")) {
      cfg.routingOpts.useStations = p.getBool(secStr, "routing_use_stations");
    } else {
      cfg.routingOpts.useStations = true;
    }

    if (p.hasKey(secStr, "osm_filter_keep")) {
      for (const auto& kvs : p.getStrArr(sec.first, "osm_filter_keep", ' ')) {
        auto fRule = getFRule(kvs);
        cfg.osmBuildOpts.keepFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    for (uint8_t i = 0; i < 8; i++) {
      std::string name = std::string("osm_filter_lvl") + std::to_string(i);
      if (p.hasKey(secStr, name)) {
        for (const auto& kvs : p.getStrArr(sec.first, name, ' ')) {
          auto fRule = getFRule(kvs);
          cfg.osmBuildOpts.levelFilters[i][fRule.kv.first].insert(
              osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
        }
      }
    }

    if (p.hasKey(secStr, "osm_filter_drop")) {
      for (const auto& kvs : p.getStrArr(sec.first, "osm_filter_drop", ' ')) {
        auto fRule = getFRule(kvs);
        cfg.osmBuildOpts.dropFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_max_snap_level")) {
      cfg.osmBuildOpts.maxSnapLevel = p.getInt(sec.first, "osm_max_snap_level");
    } else {
      cfg.osmBuildOpts.maxSnapLevel = 7;
    }

    if (p.hasKey(secStr, "osm_filter_nohup")) {
      for (const auto& kvs : p.getStrArr(sec.first, "osm_filter_nohup", ' ')) {
        auto fRule = getFRule(kvs);
        cfg.osmBuildOpts.noHupFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_oneway")) {
      for (const auto& kvs : p.getStrArr(sec.first, "osm_filter_oneway", ' ')) {
        auto fRule = getFRule(kvs);
        cfg.osmBuildOpts.oneWayFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_oneway_reverse")) {
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_filter_oneway_reverse", ' ')) {
        auto fRule = getFRule(kvs);
        cfg.osmBuildOpts.oneWayFilterRev[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_undirected")) {
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_filter_undirected", ' ')) {
        auto fRule = getFRule(kvs);
        cfg.osmBuildOpts.twoWayFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_station")) {
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_filter_station", ' ')) {
        auto fRule = getFRule(kvs);
        cfg.osmBuildOpts.stationFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_station_blocker")) {
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_filter_station_blocker", ' ')) {
        auto fRule = getFRule(kvs);
        cfg.osmBuildOpts.stationBlockerFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_turning_circle")) {
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_filter_turning_circle", ' ')) {
        auto fRule = getFRule(kvs);
        cfg.osmBuildOpts.turnCycleFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_node_positive_restriction")) {
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_node_positive_restriction", ' ')) {
        auto fRule = getFRule(kvs);
        cfg.osmBuildOpts.restrPosRestr[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_node_negative_restriction")) {
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_node_negative_restriction", ' ')) {
        auto fRule = getFRule(kvs);
        cfg.osmBuildOpts.restrNegRestr[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_filter_no_restriction")) {
      for (const auto& kvs :
           p.getStrArr(sec.first, "osm_filter_no_restriction", ' ')) {
        auto fRule = getFRule(kvs);
        cfg.osmBuildOpts.noRestrFilter[fRule.kv.first].insert(
            osm::AttrFlagPair(fRule.kv.second, getFlags(fRule.flags)));
      }
    }

    if (p.hasKey(secStr, "osm_station_name_attrs")) {
      for (const std::string& r :
           p.getStrArr(sec.first, "osm_station_name_attrs", ' ')) {
        cfg.osmBuildOpts.statAttrRules.nameRule.push_back(getDeepAttrRule(r));
      }
    }

    if (p.hasKey(secStr, "osm_track_number_tags")) {
      for (const std::string& r :
           p.getStrArr(sec.first, "osm_track_number_tags", ' ')) {
        cfg.osmBuildOpts.statAttrRules.platformRule.push_back(
            getDeepAttrRule(r));
      }
    }

    if (p.hasKey(secStr, "osm_station_id_attrs")) {
      for (const std::string& r :
           p.getStrArr(sec.first, "osm_station_id_attrs", ' ')) {
        cfg.osmBuildOpts.statAttrRules.idRule.push_back(getDeepAttrRule(r));
      }
    }

    if (p.hasKey(secStr, "osm_edge_track_number_tags")) {
      for (const std::string& r :
           p.getStrArr(sec.first, "osm_edge_track_number_tags", ' ')) {
        cfg.osmBuildOpts.edgePlatformRules.push_back(getDeepAttrRule(r));
      }
    }

    if (p.hasKey(secStr, "osm_station_group_attrs")) {
      LOG(WARN) << "Option osm_station_group_attrs has been removed.";
    }

    // default value, to enable color writing on old configs
    cfg.osmBuildOpts.relLinerules.colorRule = {"colour", "color"};

    if (p.hasKey(secStr, "osm_line_relation_tags")) {
      auto arr = p.getStrArr(secStr, "osm_line_relation_tags", ' ');

      for (const auto& ruleStr : arr) {
        auto rule = getKv(ruleStr);
        auto tags = util::split(rule.second, ',');
        if (rule.first == "from_name")
          cfg.osmBuildOpts.relLinerules.fromNameRule = tags;
        else if (rule.first == "to_name")
          cfg.osmBuildOpts.relLinerules.toNameRule = tags;
        else if (rule.first == "line_name")
          cfg.osmBuildOpts.relLinerules.sNameRule = tags;
        else if (rule.first == "line_color")
          cfg.osmBuildOpts.relLinerules.colorRule = tags;
      }
    }

    cfg.osmBuildOpts.maxSnapDistance = 50;
    if (p.hasKey(secStr, "osm_max_snap_distance")) {
      auto v = p.getDoubleArr(secStr, "osm_max_snap_distance", ',');
      if (v.size()) cfg.osmBuildOpts.maxSnapDistance = v.back();
    }

    cfg.osmBuildOpts.maxStationCandDistance =
        cfg.osmBuildOpts.maxSnapDistance * 2;
    if (p.hasKey(secStr, "osm_max_station_cand_distance")) {
      auto v = p.getDouble(secStr, "osm_max_station_cand_distance");
      cfg.osmBuildOpts.maxStationCandDistance = v;
    }

    if (p.hasKey(secStr, "osm_max_snap_fallback_distance")) {
      LOG(WARN) << "Option osm_max_snap_fallback_distance has been removed.";
    }

    if (p.hasKey(secStr, "osm_max_osm_station_distance")) {
      double ref = p.getDouble(secStr, "osm_max_osm_station_distance");
      cfg.osmBuildOpts.maxOsmStationDistances.push_back(ref);
    } else {
      cfg.osmBuildOpts.maxOsmStationDistances.push_back(15);
    }

    if (p.hasKey(secStr, "osm_max_node_block_distance")) {
      cfg.osmBuildOpts.maxBlockDistance =
          p.getDouble(secStr, "osm_max_node_block_distance");
    } else {
      cfg.osmBuildOpts.maxBlockDistance =
          *std::max_element(cfg.osmBuildOpts.maxOsmStationDistances.begin(),
                            cfg.osmBuildOpts.maxOsmStationDistances.end()) /
          8;
    }

    double DEF_SPEED = 85;
    for (uint8_t i = 0; i < 8; i++) {
      std::string name =
          std::string("routing_lvl") + std::to_string(i) + "_fac";
      if (p.hasKey(secStr, name)) {
        double f = p.getPosDouble(sec.first, name);
        LOG(WARN) << "Option " << name << " is deprecated, use osm_lvl"
                  << std::to_string(i) << "_avg_speed instead.";
        double v = DEF_SPEED / f;
        LOG(DEBUG) << " (using osm_lvl" << std::to_string(i) << "_avg_speed of "
                   << v << " instead)";
        cfg.osmBuildOpts.levelDefSpeed[i] = v * 0.2777;  // store in m/s
      }
    }

    for (uint8_t i = 0; i < 8; i++) {
      std::string name =
          std::string("osm_lvl") + std::to_string(i) + "_avg_speed";
      if (p.hasKey(secStr, name)) {
        double v = p.getPosDouble(sec.first, name);
        cfg.osmBuildOpts.levelDefSpeed[i] = v * 0.2777;  // store in m/s
      }
    }

    if (p.hasKey(secStr, "routing_one_way_meter_punish_fac")) {
      LOG(WARN) << "Option routing_one_way_meter_punish_fac is deprecated, use "
                   "osm_one_way_speed_penalty_fac instead.";
      cfg.osmBuildOpts.oneWaySpeedPen =
          1 + p.getPosDouble(secStr, "routing_one_way_meter_punish_fac");
      LOG(DEBUG) << " (using osm_one_way_speed_penalty_fac of "
                 << cfg.osmBuildOpts.oneWaySpeedPen << " instead)";
    } else {
      cfg.osmBuildOpts.oneWaySpeedPen = 1;
    }

    if (p.hasKey(secStr, "osm_one_way_speed_penalty_fac")) {
      cfg.osmBuildOpts.oneWaySpeedPen =
          p.getPosDouble(secStr, "osm_one_way_speed_penalty_fac");
    } else {
      // def already set above
    }

    if (p.hasKey(secStr, "osm_one_way_entry_cost")) {
      cfg.osmBuildOpts.oneWayEntryCost =
          p.getPosDouble(secStr, "osm_one_way_entry_cost");

    } else {
      cfg.osmBuildOpts.oneWayEntryCost = 0;
    }

    // take the same cost for taking restricted turns to keep
    // configuration simple
    double val = cfg.osmBuildOpts.oneWayEntryCost * 10.0;
    if (val > std::numeric_limits<uint32_t>::max()) {
      val = std::numeric_limits<uint32_t>::max();
    }

    cfg.routingOpts.turnRestrCost = val;

    if (p.hasKey(secStr, "routing_full_turn_punish")) {
      double val = p.getPosDouble(secStr, "routing_full_turn_punish");

      LOG(WARN) << "Option routing_full_turn_punish is deprecated, use "
                   "routing_full_turn_penalty instead.";

      val /= cfg.osmBuildOpts.levelDefSpeed[0];

      LOG(DEBUG) << " (using routing_full_turn_penalty of " << val
                 << " instead)";

      val *= 10.0;

      if (val > std::numeric_limits<uint32_t>::max()) {
        val = std::numeric_limits<uint32_t>::max();
      }

      cfg.routingOpts.fullTurnPunishFac = val;
    }

    if (p.hasKey(secStr, "routing_full_turn_penalty")) {
      double val = p.getPosDouble(secStr, "routing_full_turn_penalty") * 10.0;

      if (val > std::numeric_limits<uint32_t>::max()) {
        val = std::numeric_limits<uint32_t>::max();
      }

      cfg.routingOpts.fullTurnPunishFac = val;
    }

    if (p.hasKey(secStr, "routing_no_self_hops")) {
      cfg.routingOpts.noSelfHops = p.getBool(secStr, "routing_no_self_hops");
    }

    if (p.hasKey(secStr, "routing_full_turn_angle")) {
      double ang = p.getPosDouble(secStr, "routing_full_turn_angle");
      cfg.routingOpts.fullTurnAngle = ang;
      cfg.osmBuildOpts.fullTurnAngle = ang;
    } else {
      cfg.routingOpts.fullTurnAngle = 5;
      cfg.osmBuildOpts.fullTurnAngle = 5;
    }

    if (p.hasKey(secStr, "routing_snap_full_turn_angle")) {
      double ang = p.getPosDouble(secStr, "routing_snap_full_turn_angle");
      cfg.osmBuildOpts.maxAngleSnapReach = ang;
    } else {
      cfg.osmBuildOpts.maxAngleSnapReach = cfg.routingOpts.fullTurnAngle;
    }

    if (p.hasKey(secStr, "routing_pass_thru_station_punish")) {
      LOG(WARN) << "Option routing_pass_thru_station_punish has been removed.";
    }

    cfg.routingOpts.turnRestrCost *= 10.0;

    if (p.hasKey(secStr, "routing_no_lines_punish_fac")) {
      LOG(WARN) << "Option routing_no_lines_punish_fac is deprecated, use "
                   "routing_no_lines_penalty_fac instead.";

      cfg.routingOpts.noLinesPunishFact =
          1 + p.getPosDouble(secStr, "routing_no_lines_punish_fac");

      LOG(DEBUG) << " (using routing_no_lines_penalty_fac of "
                 << cfg.routingOpts.noLinesPunishFact << " instead)";
    } else {
      cfg.routingOpts.noLinesPunishFact = 1;
    }

    if (p.hasKey(secStr, "routing_no_lines_penalty_fac")) {
      cfg.routingOpts.noLinesPunishFact =
          p.getPosDouble(secStr, "routing_no_lines_penalty_fac");
    } else {
      // default already set above
    }

    // store this at two places, as we are writing the punishment into the graph
    cfg.osmBuildOpts.noLinesPunishFact = cfg.routingOpts.noLinesPunishFact;

    if (p.hasKey(secStr, "routing_line_unmatched_punish_fac")) {
      LOG(WARN)
          << "Option routing_line_unmatched_punish_fac is deprecated, use "
             "routing_line_unmatched_time_penalty_fac, "
             "routing_line_station_from_unmatched_time_penalty, and "
             "routing_line_station_to_unmatched_time_penalty instead.";

      cfg.routingOpts.lineUnmatchedPunishFact =
          1 + p.getPosDouble(secStr, "routing_line_unmatched_punish_fac") / 3;

      cfg.routingOpts.lineNameFromUnmatchedPunishFact =
          1 + p.getPosDouble(secStr, "routing_line_unmatched_punish_fac") / 3;

      cfg.routingOpts.lineNameToUnmatchedPunishFact =
          1 + p.getPosDouble(secStr, "routing_line_unmatched_punish_fac") / 3;

      LOG(DEBUG) << " (using routing_line_unmatched_punish_fac of "
                 << cfg.routingOpts.lineUnmatchedPunishFact << " instead)";
      LOG(DEBUG)
          << " (using routing_line_station_from_unmatched_time_penalty of "
          << cfg.routingOpts.lineNameFromUnmatchedPunishFact << " instead)";
      LOG(DEBUG) << " (using routing_line_station_to_unmatched_time_penalty of "
                 << cfg.routingOpts.lineNameToUnmatchedPunishFact
                 << " instead)";
    }

    if (p.hasKey(secStr, "routing_line_unmatched_time_penalty_fac")) {
      cfg.routingOpts.lineUnmatchedPunishFact =
          p.getPosDouble(secStr, "routing_line_unmatched_time_penalty_fac");
    }

    if (p.hasKey(secStr, "routing_line_station_from_unmatched_time_penalty")) {
      cfg.routingOpts.lineNameFromUnmatchedPunishFact = p.getPosDouble(
          secStr, "routing_line_station_from_unmatched_time_penalty");
    }

    if (p.hasKey(secStr, "routing_line_station_to_unmatched_time_penalty")) {
      cfg.routingOpts.lineNameToUnmatchedPunishFact = p.getPosDouble(
          secStr, "routing_line_station_to_unmatched_time_penalty");
    }

    if (p.hasKey(secStr, "routing_platform_unmatched_punish")) {
      LOG(WARN)
          << "Option routing_platform_unmatched_punish is deprecated, use "
             "routing_platform_unmatched_penalty instead.";
      cfg.routingOpts.platformUnmatchedPen =
          p.getPosDouble(secStr, "routing_platform_unmatched_punish");

      cfg.routingOpts.platformUnmatchedPen =
          cfg.routingOpts.platformUnmatchedPen *
          (DEF_TRANS_PEN / cfg.osmBuildOpts.levelDefSpeed[0]);

      LOG(DEBUG) << " (using routing_platform_unmatched_penalty of "
                 << cfg.routingOpts.platformUnmatchedPen << " instead)";
    } else {
      cfg.routingOpts.platformUnmatchedPen = 0;
    }

    if (p.hasKey(secStr, "routing_platform_unmatched_penalty")) {
      cfg.routingOpts.platformUnmatchedPen =
          p.getPosDouble(secStr, "routing_platform_unmatched_penalty");
    } else {
      // default already set above
    }

    if (p.hasKey(secStr, "routing_transition_penalty_fac")) {
      cfg.routingOpts.transitionPen =
          p.getPosDouble(secStr, "routing_transition_penalty_fac");
    } else {
      cfg.routingOpts.transitionPen = DEF_TRANS_PEN;
    }

    if (p.hasKey(secStr, "routing_station_distance_punish_fac")) {
      cfg.routingOpts.stationDistPenFactor =
          p.getPosDouble(secStr, "routing_station_distance_punish_fac");
      LOG(WARN) << "Option routing_station_distance_punish_fac is deprecated, "
                   "use routing_station_move_penalty_fac instead.";
      cfg.routingOpts.stationDistPenFactor =
          cfg.routingOpts.stationDistPenFactor *
          (DEF_TRANS_PEN / cfg.osmBuildOpts.levelDefSpeed[0]);
      LOG(DEBUG) << " (using routing_station_move_penalty_fac of "
                 << cfg.routingOpts.stationDistPenFactor << " instead)";
    } else {
      cfg.routingOpts.stationDistPenFactor =
          cfg.routingOpts.stationDistPenFactor *
          (DEF_TRANS_PEN / cfg.osmBuildOpts.levelDefSpeed[0]);
    }

    if (p.hasKey(secStr, "routing_station_move_penalty_fac")) {
      cfg.routingOpts.stationDistPenFactor =
          p.getPosDouble(secStr, "routing_station_move_penalty_fac");
    } else {
      // the default value was already set above
    }

    if (p.hasKey(secStr, "routing_non_osm_station_punish")) {
      cfg.routingOpts.nonStationPen =
          p.getPosDouble(secStr, "routing_non_osm_station_punish");
      LOG(WARN) << "Option routing_non_osm_station_punish is deprecated, use "
                   "routing_non_station_penalty instead.";
      cfg.routingOpts.nonStationPen =
          cfg.routingOpts.nonStationPen *
          (DEF_TRANS_PEN / cfg.osmBuildOpts.levelDefSpeed[0]);
      LOG(DEBUG) << " (using routing_non_station_penalty of "
                 << cfg.routingOpts.nonStationPen << " instead)";
    } else {
      cfg.routingOpts.nonStationPen = 0;
    }

    if (p.hasKey(secStr, "routing_non_station_penalty")) {
      cfg.routingOpts.nonStationPen =
          p.getPosDouble(secStr, "routing_non_station_penalty");
    } else {
      // default was already set above
    }

    if (p.hasKey(secStr, "routing_station_unmatched_penalty")) {
      cfg.routingOpts.stationUnmatchedPen =
          p.getPosDouble(secStr, "routing_station_unmatched_penalty");
    } else {
      cfg.routingOpts.stationUnmatchedPen = cfg.routingOpts.nonStationPen / 2;
    }

    if (p.hasKey(secStr, "station_normalize_chain")) {
      try {
        auto arr = p.getStrArr(secStr, "station_normalize_chain", ';');
        cfg.osmBuildOpts.statNormzer = trgraph::Normalizer(getNormRules(arr));
      } catch (const std::exception& e) {
        throw ParseExc(p.getVal(secStr, "station_normalize_chain").line,
                       p.getVal(secStr, "station_normalize_chain").pos,
                       "<valid regular expression>",
                       std::string("<regex error: ") + e.what() + ">",
                       p.getVal(secStr, "station_normalize_chain").file);
      }
    }

    if (p.hasKey(secStr, "track_normalize_chain")) {
      try {
        auto arr = p.getStrArr(secStr, "track_normalize_chain", ';');
        cfg.osmBuildOpts.trackNormzer = trgraph::Normalizer(getNormRules(arr));
      } catch (const std::exception& e) {
        throw ParseExc(p.getVal(secStr, "track_normalize_chain").line,
                       p.getVal(secStr, "track_normalize_chain").pos,
                       "<valid regular expression>",
                       std::string("<regex error: ") + e.what() + ">",
                       p.getVal(secStr, "track_normalize_chain").file);
      }
    }

    if (p.hasKey(secStr, "line_normalize_chain")) {
      try {
        auto arr = p.getStrArr(secStr, "line_normalize_chain", ';');
        cfg.osmBuildOpts.lineNormzer = trgraph::Normalizer(getNormRules(arr));
      } catch (const std::exception& e) {
        throw ParseExc(p.getVal(secStr, "line_normalize_chain").line,
                       p.getVal(secStr, "line_normalize_chain").pos,
                       "<valid regular expression>",
                       std::string("<regex error: ") + e.what() + ">",
                       p.getVal(secStr, "line_normalize_chain").file);
      }
    }

    if (p.hasKey(secStr, "station_id_normalize_chain")) {
      try {
        auto arr = p.getStrArr(secStr, "station_id_normalize_chain", ';');
        cfg.osmBuildOpts.idNormzer = trgraph::Normalizer(getNormRules(arr));
      } catch (const std::exception& e) {
        throw ParseExc(p.getVal(secStr, "station_id_normalize_chain").line,
                       p.getVal(secStr, "station_id_normalize_chain").pos,
                       "<valid regular expression>",
                       std::string("<regex error: ") + e.what() + ">",
                       p.getVal(secStr, "station_id_normalize_chain").file);
      }
    }

    // determine the maximum possible speed for this config, this is later
    // used to filter out station which are so far out of reach we don't
    // have to consider them for the bounding box calculation
    cfg.osmBuildOpts.maxSpeed = 0;
    cfg.osmBuildOpts.maxSpeedCorFac = 1;
    for (size_t i = 0; i < 8; i++) {
      if (cfg.osmBuildOpts.levelDefSpeed[i] > cfg.osmBuildOpts.maxSpeed)
        cfg.osmBuildOpts.maxSpeed = cfg.osmBuildOpts.levelDefSpeed[i];
    }

    if (cfg.routingOpts.lineUnmatchedPunishFact < 1)
      cfg.osmBuildOpts.maxSpeedCorFac *=
          cfg.routingOpts.lineUnmatchedPunishFact;
    if (cfg.routingOpts.lineNameFromUnmatchedPunishFact < 1)
      cfg.osmBuildOpts.maxSpeedCorFac *=
          cfg.routingOpts.lineNameFromUnmatchedPunishFact;
    if (cfg.routingOpts.lineNameToUnmatchedPunishFact < 1)
      cfg.osmBuildOpts.maxSpeedCorFac *=
          cfg.routingOpts.lineNameToUnmatchedPunishFact;

    if (cfg.routingOpts.noLinesPunishFact < 1)
      cfg.osmBuildOpts.maxSpeedCorFac *= cfg.routingOpts.noLinesPunishFact;

    if (cfg.osmBuildOpts.oneWaySpeedPen < 1)
      cfg.osmBuildOpts.maxSpeedCorFac *= cfg.osmBuildOpts.oneWaySpeedPen;

    cfg.osmBuildOpts.maxSpeed /= cfg.osmBuildOpts.maxSpeedCorFac;

    bool found = false;

    for (auto& exCfg : _cfgs) {
      if (cfg == exCfg) {
        for (auto mot :
             ad::cppgtfs::gtfs::flat::Route::getTypesFromString(secStr)) {
          exCfg.mots.insert(mot);
        }
        found = true;
        break;
      }
    }

    if (!found) {
      cfg.mots = ad::cppgtfs::gtfs::flat::Route::getTypesFromString(secStr);
      _cfgs.push_back(cfg);
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

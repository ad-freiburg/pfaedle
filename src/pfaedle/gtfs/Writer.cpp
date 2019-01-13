// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <utility>
#include "ad/cppgtfs/Parser.h"
#include "ad/cppgtfs/Writer.h"
#include "ad/cppgtfs/gtfs/flat/Agency.h"
#include "ad/util/CsvWriter.h"
#include "pfaedle/gtfs/Writer.h"

using ad::util::CsvWriter;
using ad::cppgtfs::Parser;
using pfaedle::gtfs::Writer;

// ____________________________________________________________________________
bool Writer::write(gtfs::Feed* sourceFeed, const std::string& path) const {
  std::ofstream fs;
  std::ifstream is;
  std::string gtfsPath(path);
  std::string curFile;
  std::string curFileTg;

  curFile = getTmpFName("agency.txt");
  curFileTg = gtfsPath + "/agency.txt";
  fs.open(curFile.c_str());
  if (!fs.good()) cannotWrite(curFile, curFileTg);
  writeAgency(sourceFeed, &fs);
  fs.close();
  std::rename(curFile.c_str(), curFileTg.c_str());

  curFile = getTmpFName("stops.txt");
  curFileTg = gtfsPath + "/stops.txt";
  fs.open(curFile.c_str());
  if (!fs.good()) cannotWrite(curFile, curFileTg);
  writeStops(sourceFeed, &fs);
  fs.close();
  std::rename(curFile.c_str(), curFileTg.c_str());

  curFile = getTmpFName("routes.txt");
  curFileTg = gtfsPath + "/routes.txt";
  fs.open(curFile.c_str());
  if (!fs.good()) cannotWrite(curFile, curFileTg);
  writeRoutes(sourceFeed, &fs);
  fs.close();
  std::rename(curFile.c_str(), curFileTg.c_str());

  is.open((sourceFeed->getPath() + "/calendar.txt").c_str());
  if (is.good()) {
    is.close();
    curFile = getTmpFName("calendar.txt");
    curFileTg = gtfsPath + "/calendar.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeCalendar(sourceFeed, &fs);
    fs.close();
    std::rename(curFile.c_str(), curFileTg.c_str());
  }

  is.open((sourceFeed->getPath() + "/calendar_dates.txt").c_str());
  if (is.good()) {
    is.close();
    curFile = getTmpFName("calendar_dates.txt");
    curFileTg = gtfsPath + "/calendar_dates.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeCalendarDates(sourceFeed, &fs);
    fs.close();
    std::rename(curFile.c_str(), curFileTg.c_str());
  }

  is.open((sourceFeed->getPath() + "/transfers.txt").c_str());
  if (is.good()) {
    is.close();
    curFile = getTmpFName("transfers.txt");
    curFileTg = gtfsPath + "/transfers.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeTransfers(sourceFeed, &fs);
    fs.close();
    std::rename(curFile.c_str(), curFileTg.c_str());
  }

  is.open((sourceFeed->getPath() + "/fare_attributes.txt").c_str());
  if (is.good()) {
    is.close();
    curFile = getTmpFName("fare_attributes.txt");
    curFileTg = gtfsPath + "/fare_attributes.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeFares(sourceFeed, &fs);
    fs.close();
    std::rename(curFile.c_str(), curFileTg.c_str());
  }

  is.open((sourceFeed->getPath() + "/fare_rules.txt").c_str());
  if (is.good()) {
    is.close();
    curFile = getTmpFName("fare_rules.txt");
    curFileTg = gtfsPath + "/fare_rules.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeFareRules(sourceFeed, &fs);
    fs.close();
    std::rename(curFile.c_str(), curFileTg.c_str());
  }

  is.close();
  curFile = getTmpFName("shapes.txt");
  curFileTg = gtfsPath + "/shapes.txt";
  fs.open(curFile.c_str());
  if (!fs.good()) cannotWrite(curFile, curFileTg);
  writeShapes(sourceFeed, &fs);
  fs.close();
  std::rename(curFile.c_str(), curFileTg.c_str());

  is.close();
  curFile = getTmpFName("trips.txt");
  curFileTg = gtfsPath + "/trips.txt";
  fs.open(curFile.c_str());
  if (!fs.good()) cannotWrite(curFile, curFileTg);
  bool hasFreqs = writeTrips(sourceFeed, &fs);
  fs.close();
  std::rename(curFile.c_str(), curFileTg.c_str());

  is.open((sourceFeed->getPath() + "/frequencies.txt").c_str());
  if (hasFreqs && is.good()) {
    is.close();
    curFile = getTmpFName("frequencies.txt");
    curFileTg = gtfsPath + "/frequencies.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeFrequencies(sourceFeed, &fs);
    fs.close();
    std::rename(curFile.c_str(), curFileTg.c_str());
  }

  is.close();
  curFile = getTmpFName("stop_times.txt");
  curFileTg = gtfsPath + "/stop_times.txt";
  fs.open(curFile.c_str());
  if (!fs.good()) cannotWrite(curFile, curFileTg);
  writeStopTimes(sourceFeed, &fs);
  fs.close();
  std::rename(curFile.c_str(), curFileTg.c_str());

  if (!sourceFeed->getPublisherUrl().empty() &&
      !sourceFeed->getPublisherName().empty()) {
    curFile = getTmpFName("feed_info.txt");
    curFileTg = gtfsPath + "/feed_info.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeFeedInfo(sourceFeed, &fs);
    fs.close();
    std::rename(curFile.c_str(), curFileTg.c_str());
  }

  return true;
}

// ____________________________________________________________________________
bool Writer::writeFeedInfo(gtfs::Feed* f, std::ostream* os) const {
  auto csvw = ad::cppgtfs::Writer::getFeedInfoCsvw(os);
  csvw.flushLine();
  csvw.writeString(f->getPublisherName());
  csvw.writeString(f->getPublisherUrl());
  csvw.writeString(f->getLang());
  if (!f->getStartDate().empty())
    csvw.writeInt(f->getStartDate().getYYYYMMDD());
  else
    csvw.skip();
  if (!f->getEndDate().empty())
    csvw.writeInt(f->getEndDate().getYYYYMMDD());
  else
    csvw.skip();
  csvw.writeString(f->getVersion());
  csvw.flushLine();

  return true;
}

// ____________________________________________________________________________
bool Writer::writeAgency(gtfs::Feed* sourceFeed, std::ostream* os) const {
  std::ifstream fs;
  fs.open((sourceFeed->getPath() + "/agency.txt").c_str());

  CsvParser csvp(&fs);
  Parser p;
  ad::cppgtfs::Writer w;

  CsvWriter csvw = ad::cppgtfs::Writer::getAgencyCsvw(os);
  csvw.flushLine();

  ad::cppgtfs::gtfs::flat::Agency fa;
  auto flds = Parser::getAgencyFlds(&csvp);

  while (p.nextAgency(&csvp, &fa, flds)) {
    w.writeAgency(fa, &csvw);
  }
  fs.close();

  return true;
}

// ____________________________________________________________________________
bool Writer::writeStops(gtfs::Feed* sourceFeed, std::ostream* os) const {
  std::ifstream fs;
  fs.open((sourceFeed->getPath() + "/stops.txt").c_str());

  CsvParser csvp(&fs);
  Parser p;
  ad::cppgtfs::Writer w;

  CsvWriter csvw = ad::cppgtfs::Writer::getStopsCsvw(os);
  csvw.flushLine();

  ad::cppgtfs::gtfs::flat::Stop s;
  auto flds = Parser::getStopFlds(&csvp);

  while (p.nextStop(&csvp, &s, flds)) {
    w.writeStop(s, &csvw);
  }
  fs.close();

  return true;
}

// ____________________________________________________________________________
bool Writer::writeRoutes(gtfs::Feed* sourceFeed, std::ostream* os) const {
  std::ifstream fs;
  fs.open((sourceFeed->getPath() + "/routes.txt").c_str());

  CsvParser csvp(&fs);
  Parser p;
  ad::cppgtfs::Writer w;

  CsvWriter csvw = ad::cppgtfs::Writer::getRoutesCsvw(os);
  csvw.flushLine();

  ad::cppgtfs::gtfs::flat::Route s;
  auto flds = Parser::getRouteFlds(&csvp);

  while (p.nextRoute(&csvp, &s, flds)) {
    w.writeRoute(s, &csvw);
  }
  fs.close();

  return true;
}

// ____________________________________________________________________________
bool Writer::writeCalendar(gtfs::Feed* sourceFeed, std::ostream* os) const {
  std::ifstream fs;
  fs.open((sourceFeed->getPath() + "/calendar.txt").c_str());

  CsvParser csvp(&fs);
  Parser p;
  ad::cppgtfs::Writer w;

  CsvWriter csvw = ad::cppgtfs::Writer::getCalendarCsvw(os);
  csvw.flushLine();

  ad::cppgtfs::gtfs::flat::Calendar c;
  auto flds = Parser::getCalendarFlds(&csvp);

  while (p.nextCalendar(&csvp, &c, flds)) {
    w.writeCalendar(c, &csvw);
  }
  fs.close();

  return true;
}

// ____________________________________________________________________________
bool Writer::writeCalendarDates(gtfs::Feed* sourceFeed,
                                std::ostream* os) const {
  std::ifstream fs;
  fs.open((sourceFeed->getPath() + "/calendar_dates.txt").c_str());

  CsvParser csvp(&fs);
  Parser p;
  ad::cppgtfs::Writer w;

  CsvWriter csvw = ad::cppgtfs::Writer::getCalendarDatesCsvw(os);
  csvw.flushLine();

  ad::cppgtfs::gtfs::flat::CalendarDate c;
  auto flds = Parser::getCalendarDateFlds(&csvp);

  while (p.nextCalendarDate(&csvp, &c, flds)) {
    w.writeCalendarDate(c, &csvw);
  }
  fs.close();

  return true;
}

// ____________________________________________________________________________
bool Writer::writeFrequencies(gtfs::Feed* sourceFeed, std::ostream* os) const {
  std::ifstream fs;
  fs.open((sourceFeed->getPath() + "/frequencies.txt").c_str());

  CsvParser csvp(&fs);
  Parser p;
  ad::cppgtfs::Writer w;

  CsvWriter csvw = ad::cppgtfs::Writer::getFrequencyCsvw(os);
  csvw.flushLine();

  ad::cppgtfs::gtfs::flat::Frequency f;
  auto flds = Parser::getFrequencyFlds(&csvp);

  while (p.nextFrequency(&csvp, &f, flds)) {
    w.writeFrequency(f, &csvw);
  }
  fs.close();

  return true;
}

// ____________________________________________________________________________
bool Writer::writeTransfers(gtfs::Feed* sourceFeed, std::ostream* os) const {
  std::ifstream fs;
  fs.open((sourceFeed->getPath() + "/transfers.txt").c_str());

  CsvParser csvp(&fs);
  Parser p;
  ad::cppgtfs::Writer w;

  CsvWriter csvw = ad::cppgtfs::Writer::getTransfersCsvw(os);
  csvw.flushLine();

  ad::cppgtfs::gtfs::flat::Transfer t;
  auto flds = Parser::getTransfersFlds(&csvp);

  while (p.nextTransfer(&csvp, &t, flds)) {
    w.writeTransfer(t, &csvw);
  }
  fs.close();

  return true;
}

// ____________________________________________________________________________
bool Writer::writeFares(gtfs::Feed* sourceFeed, std::ostream* os) const {
  std::ifstream fs;
  fs.open((sourceFeed->getPath() + "/fare_attributes.txt").c_str());

  CsvParser csvp(&fs);
  Parser p;
  ad::cppgtfs::Writer w;

  CsvWriter csvw = ad::cppgtfs::Writer::getFaresCsvw(os);
  csvw.flushLine();

  ad::cppgtfs::gtfs::flat::Fare f;
  auto flds = Parser::getFareFlds(&csvp);

  while (p.nextFare(&csvp, &f, flds)) {
    w.writeFare(f, &csvw);
  }
  fs.close();

  return true;
}

// ____________________________________________________________________________
bool Writer::writeFareRules(gtfs::Feed* sourceFeed, std::ostream* os) const {
  std::ifstream fs;
  fs.open((sourceFeed->getPath() + "/fare_rules.txt").c_str());

  CsvParser csvp(&fs);
  Parser p;
  ad::cppgtfs::Writer w;

  CsvWriter csvw = ad::cppgtfs::Writer::getFareRulesCsvw(os);
  csvw.flushLine();

  ad::cppgtfs::gtfs::flat::FareRule f;
  auto flds = Parser::getFareRuleFlds(&csvp);

  while (p.nextFareRule(&csvp, &f, flds)) {
    w.writeFareRule(f, &csvw);
  }
  fs.close();

  return true;
}

// ____________________________________________________________________________
bool Writer::writeShapes(gtfs::Feed* sourceFeed, std::ostream* os) const {
  std::ifstream fs;
  fs.open((sourceFeed->getPath() + "/shapes.txt").c_str());

  CsvWriter csvw = ad::cppgtfs::Writer::getShapesCsvw(os);
  csvw.flushLine();
  ad::cppgtfs::gtfs::flat::ShapePoint sp;
  ad::cppgtfs::Writer w;

  if (fs.good()) {
    CsvParser csvp(&fs);
    Parser p;

    auto flds = Parser::getShapeFlds(&csvp);

    std::string curShapeId;
    std::string curSkipShapeId;

    while (p.nextShapePoint(&csvp, &sp, flds)) {
      if (sp.id == curSkipShapeId) continue;
      if (sp.id != curShapeId) {
        if (sourceFeed->getShapes().has(sp.id)) {
          curShapeId = sp.id;
        } else {
          curSkipShapeId = sp.id;
          continue;
        }
      }

      w.writeShapePoint(sp, &csvw);
    }
  }

  sourceFeed->getShapes().open();
  while (sourceFeed->getShapes().nextStoragePt(&sp)) {
    w.writeShapePoint(sp, &csvw);
  }

  fs.close();

  return true;
}

// ____________________________________________________________________________
bool Writer::writeTrips(gtfs::Feed* sourceFeed, std::ostream* os) const {
  ad::cppgtfs::Writer w;
  bool hasFreqs = false;

  CsvWriter csvw = ad::cppgtfs::Writer::getTripsCsvw(os);
  csvw.flushLine();

  for (auto t : sourceFeed->getTrips()) {
    if (t.getFrequencies().size()) hasFreqs = true;
    w.writeTrip(t.getFlat(), &csvw);
  }

  return hasFreqs;
}

// ____________________________________________________________________________
bool Writer::writeStopTimes(gtfs::Feed* sourceFeed, std::ostream* os) const {
  std::ifstream fs;
  fs.open((sourceFeed->getPath() + "/stop_times.txt").c_str());

  CsvParser csvp(&fs);
  Parser p;
  ad::cppgtfs::Writer w;

  CsvWriter csvw = ad::cppgtfs::Writer::getStopTimesCsvw(os);
  csvw.flushLine();

  ad::cppgtfs::gtfs::flat::StopTime st;
  auto flds = Parser::getStopTimeFlds(&csvp);

  std::string curTripId;
  Trip* cur = 0;

  while (p.nextStopTime(&csvp, &st, flds)) {
    // we may have changed to distance field
    if (curTripId != st.trip) {
      cur = sourceFeed->getTrips().get(st.trip);
      curTripId = st.trip;
    }
    for (const auto& stN : cur->getStopTimes()) {
      if (stN.getSeq() == st.sequence)
        st.shapeDistTravelled = stN.getShapeDistanceTravelled();
    }

    w.writeStopTime(st, &csvw);
  }
  fs.close();

  return true;
}

// ___________________________________________________________________________
void Writer::cannotWrite(const std::string& file, const std::string& file2) {
  std::stringstream ss;
  ss << "(temporary file for " << file2 << ") Could not write to file";
  throw ad::cppgtfs::WriterException(ss.str(), file);
}

// _____________________________________________________________________________
std::string Writer::getTmpFName(const std::string& postf) {
  std::string f = ".pfaedle-tmp-" + postf;

  while (access(f.c_str(), F_OK) != -1) {
    std::stringstream ss;
    ss << ".pfaedle-tmp-";
    ss << postf << "-";
    ss << std::rand();
    f = ss.str().c_str();
  }

  return f;
}

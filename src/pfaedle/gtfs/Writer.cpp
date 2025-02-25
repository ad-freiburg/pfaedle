// Copyright 2018, University of Freiburg,
// Chair of Algorithms and Data Structures.
// Authors: Patrick Brosi <brosi@informatik.uni-freiburg.de>

#include <stdio.h>

#ifdef LIBZIP_FOUND
#include <zip.h>
#endif

#include <cstdio>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "ad/cppgtfs/Parser.h"
#include "ad/cppgtfs/Writer.h"
#include "ad/cppgtfs/gtfs/flat/Agency.h"
#include "ad/util/CsvWriter.h"
#include "pfaedle/gtfs/Writer.h"

using ad::cppgtfs::Parser;
using ad::util::CsvWriter;
#ifdef LIBZIP_FOUND
using ad::util::ZipCsvParser;
#endif
using pfaedle::gtfs::Writer;
using util::getTmpFName;

// ____________________________________________________________________________
void Writer::write(gtfs::Feed* sourceFeed, const std::string& path) const {
  bool toZip =
      (path.size() > 3 && 0 == path.compare(path.size() - 4, 4, ".zip"));

  std::ofstream fs;
  std::string gtfsPath(path);
  std::string curFile;
  std::string curFileTg;

  std::string tmpZip;
  std::string zipFileName;

  if (gtfsPath.size() == 0) gtfsPath = ".";

#ifdef LIBZIP_FOUND
  zip* za = 0;

  if (toZip) {
    const size_t slashIdx = path.rfind('/');
    if (slashIdx != std::string::npos) {
      zipFileName = path.substr(slashIdx + 1, -1);
      gtfsPath = path.substr(0, slashIdx);
    } else {
      zipFileName = path;
      gtfsPath = ".";
    }

    tmpZip = getTmpFName(gtfsPath, ".pfaedle-tmp", zipFileName);

    int zipErr = 0;
    za = zip_open(tmpZip.c_str(), ZIP_CREATE | ZIP_TRUNCATE, &zipErr);

    if (zipErr != 0) {
      char errBuf[100];
      zip_error_to_str(errBuf, sizeof(errBuf), zipErr, errno);
      cannotWrite(tmpZip, gtfsPath + "/" + zipFileName);
      std::stringstream ss;
      ss << "(temporary file for " << (gtfsPath + "/" + zipFileName)
         << ") Could not open ZIP file, reason was: " << errBuf;
      throw ad::cppgtfs::WriterException(ss.str(), tmpZip);
    }
#else
  if (toZip) {
    throw ad::cppgtfs::WriterException(
        "Could not output ZIP file, pfaedle was compiled without libzip", path);
#endif
  } else {
    mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  }

  try {
    Parser ip(sourceFeed->getPath());

    curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "agency.txt");
    curFileTg = gtfsPath + "/agency.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeAgency(sourceFeed, &fs);
    fs.close();

    if (toZip) {
#ifdef LIBZIP_FOUND
      moveIntoZip(za, curFile, "agency.txt");
#endif
    } else {
      if (std::rename(curFile.c_str(), curFileTg.c_str()))
        cannotWrite(curFileTg);
    }

    curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "stops.txt");
    curFileTg = gtfsPath + "/stops.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeStops(sourceFeed, &fs);
    fs.close();

    if (toZip) {
#ifdef LIBZIP_FOUND
      moveIntoZip(za, curFile, "stops.txt");
#endif
    } else {
      if (std::rename(curFile.c_str(), curFileTg.c_str()))
        cannotWrite(curFileTg);
    }

    curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "routes.txt");
    curFileTg = gtfsPath + "/routes.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeRoutes(sourceFeed, &fs);
    fs.close();

    if (toZip) {
#ifdef LIBZIP_FOUND
      moveIntoZip(za, curFile, "routes.txt");
#endif
    } else {
      if (std::rename(curFile.c_str(), curFileTg.c_str()))
        cannotWrite(curFileTg);
    }

    auto csvp = ip.getCsvParser("calendar.txt");
    if (csvp->isGood()) {
      curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "calendar.txt");
      curFileTg = gtfsPath + "/calendar.txt";
      fs.open(curFile.c_str());
      if (!fs.good()) cannotWrite(curFile, curFileTg);
      writeCalendar(sourceFeed, &fs);
      fs.close();
      if (toZip) {
#ifdef LIBZIP_FOUND
        moveIntoZip(za, curFile, "calendar.txt");
#endif
      } else {
        if (std::rename(curFile.c_str(), curFileTg.c_str()))
          cannotWrite(curFileTg);
      }
    }

    csvp = ip.getCsvParser("calendar_dates.txt");
    if (csvp->isGood()) {
      curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "calendar_dates.txt");
      curFileTg = gtfsPath + "/calendar_dates.txt";
      fs.open(curFile.c_str());
      if (!fs.good()) cannotWrite(curFile, curFileTg);
      writeCalendarDates(sourceFeed, &fs);
      fs.close();
      if (toZip) {
#ifdef LIBZIP_FOUND
        moveIntoZip(za, curFile, "calendar_dates.txt");
#endif
      } else {
        if (std::rename(curFile.c_str(), curFileTg.c_str()))
          cannotWrite(curFileTg);
      }
    }

    csvp = ip.getCsvParser("transfers.txt");
    if (csvp->isGood()) {
      curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "transfers.txt");
      curFileTg = gtfsPath + "/transfers.txt";
      fs.open(curFile.c_str());
      if (!fs.good()) cannotWrite(curFile, curFileTg);
      writeTransfers(sourceFeed, &fs);
      fs.close();
      if (toZip) {
#ifdef LIBZIP_FOUND
        moveIntoZip(za, curFile, "transfers.txt");
#endif
      } else {
        if (std::rename(curFile.c_str(), curFileTg.c_str()))
          cannotWrite(curFileTg);
      }
    }

    csvp = ip.getCsvParser("fare_attributes.txt");
    if (csvp->isGood()) {
      curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "fare_attributes.txt");
      curFileTg = gtfsPath + "/fare_attributes.txt";
      fs.open(curFile.c_str());
      if (!fs.good()) cannotWrite(curFile, curFileTg);
      writeFares(sourceFeed, &fs);
      fs.close();
      if (toZip) {
#ifdef LIBZIP_FOUND
        moveIntoZip(za, curFile, "fare_attributes.txt");
#endif
      } else {
        if (std::rename(curFile.c_str(), curFileTg.c_str()))
          cannotWrite(curFileTg);
      }
    }

    csvp = ip.getCsvParser("fare_rules.txt");
    if (csvp->isGood()) {
      curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "fare_rules.txt");
      curFileTg = gtfsPath + "/fare_rules.txt";
      fs.open(curFile.c_str());
      if (!fs.good()) cannotWrite(curFile, curFileTg);
      writeFareRules(sourceFeed, &fs);
      fs.close();
      if (toZip) {
#ifdef LIBZIP_FOUND
        moveIntoZip(za, curFile, "fare_rules.txt");
#endif
      } else {
        if (std::rename(curFile.c_str(), curFileTg.c_str()))
          cannotWrite(curFileTg);
      }
    }

    csvp = ip.getCsvParser("pathways.txt");
    if (csvp->isGood()) {
      curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "pathways.txt");
      curFileTg = gtfsPath + "/pathways.txt";
      fs.open(curFile.c_str());
      if (!fs.good()) cannotWrite(curFile, curFileTg);
      writePathways(sourceFeed, &fs);
      fs.close();

      if (toZip) {
#ifdef LIBZIP_FOUND
        moveIntoZip(za, curFile, "pathways.txt");
#endif
      } else {
        if (std::rename(curFile.c_str(), curFileTg.c_str()))
          cannotWrite(curFileTg);
      }
    }

    csvp = ip.getCsvParser("levels.txt");
    if (csvp->isGood()) {
      curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "levels.txt");
      curFileTg = gtfsPath + "/levels.txt";
      fs.open(curFile.c_str());
      if (!fs.good()) cannotWrite(curFile, curFileTg);
      writeLevels(sourceFeed, &fs);
      fs.close();

      if (toZip) {
#ifdef LIBZIP_FOUND
        moveIntoZip(za, curFile, "levels.txt");
#endif
      } else {
        if (std::rename(curFile.c_str(), curFileTg.c_str()))
          cannotWrite(curFileTg);
      }
    }

    curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "shapes.txt");
    curFileTg = gtfsPath + "/shapes.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeShapes(sourceFeed, &fs);
    fs.close();

    if (toZip) {
#ifdef LIBZIP_FOUND
      moveIntoZip(za, curFile, "shapes.txt");
#endif
    } else {
      if (std::rename(curFile.c_str(), curFileTg.c_str()))
        cannotWrite(curFileTg);
    }

    curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "trips.txt");
    curFileTg = gtfsPath + "/trips.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    bool hasFreqs = writeTrips(sourceFeed, &fs);
    fs.close();

    if (toZip) {
#ifdef LIBZIP_FOUND
      moveIntoZip(za, curFile, "trips.txt");
#endif
    } else {
      if (std::rename(curFile.c_str(), curFileTg.c_str()))
        cannotWrite(curFileTg);
    }

    csvp = ip.getCsvParser("frequencies.txt");
    if (hasFreqs && csvp->isGood()) {
      curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "frequencies.txt");
      curFileTg = gtfsPath + "/frequencies.txt";
      fs.open(curFile.c_str());
      if (!fs.good()) cannotWrite(curFile, curFileTg);
      writeFrequencies(sourceFeed, &fs);
      fs.close();

      if (toZip) {
#ifdef LIBZIP_FOUND
        moveIntoZip(za, curFile, "frequencies.txt");
#endif
      } else {
        if (std::rename(curFile.c_str(), curFileTg.c_str()))
          cannotWrite(curFileTg);
      }
    }

    curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "stop_times.txt");
    curFileTg = gtfsPath + "/stop_times.txt";
    fs.open(curFile.c_str());

    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeStopTimes(sourceFeed, &fs);
    fs.close();

    if (toZip) {
#ifdef LIBZIP_FOUND
      moveIntoZip(za, curFile, "stop_times.txt");
#endif
    } else {
      if (std::rename(curFile.c_str(), curFileTg.c_str()))
        cannotWrite(curFileTg);
    }

    if (!sourceFeed->getPublisherUrl().empty() &&
        !sourceFeed->getPublisherName().empty()) {
      curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "feed_info.txt");
      curFileTg = gtfsPath + "/feed_info.txt";
      fs.open(curFile.c_str());
      if (!fs.good()) cannotWrite(curFile, curFileTg);
      writeFeedInfo(sourceFeed, &fs);
      fs.close();

      if (toZip) {
#ifdef LIBZIP_FOUND
        moveIntoZip(za, curFile, "feed_info.txt");
#endif
      } else {
        if (std::rename(curFile.c_str(), curFileTg.c_str()))
          cannotWrite(curFileTg);
      }
    }

    curFile = getTmpFName(gtfsPath, ".pfaedle-tmp", "attributions.txt");
    curFileTg = gtfsPath + "/attributions.txt";
    fs.open(curFile.c_str());
    if (!fs.good()) cannotWrite(curFile, curFileTg);
    writeAttribution(sourceFeed, &fs);
    fs.close();

    if (toZip) {
#ifdef LIBZIP_FOUND
      moveIntoZip(za, curFile, "attributions.txt");
#endif
    } else {
      if (std::rename(curFile.c_str(), curFileTg.c_str()))
        cannotWrite(curFileTg);
    }
  } catch (...) {
#ifdef LIBZIP_FOUND
    zip_discard(za);
#endif
    throw;
  }

  if (toZip) {
#ifdef LIBZIP_FOUND
    std::string targetZipPath = gtfsPath + "/" + zipFileName;
    if (!za) cannotWrite(targetZipPath);
    zip_close(za);
    if (std::rename(tmpZip.c_str(), targetZipPath.c_str()))
      cannotWrite(targetZipPath);
#endif
  }
}

// ____________________________________________________________________________
void Writer::writeAttribution(gtfs::Feed*, std::ostream* os) const {
  auto csvw = ad::cppgtfs::Writer::getAttributionCsvw(os);

  csvw->flushLine();
  csvw->writeString("OpenStreetMap contributors");
  csvw->writeString("https://www.openstreetmap.org/copyright");
  csvw->writeInt(1);
  csvw->writeInt(0);
  csvw->writeInt(0);

  csvw->flushLine();
}

// ____________________________________________________________________________
void Writer::writeFeedInfo(gtfs::Feed* f, std::ostream* os) const {
  auto csvw = ad::cppgtfs::Writer::getFeedInfoCsvw(os);
  csvw->flushLine();
  csvw->writeString(f->getPublisherName());
  csvw->writeString(f->getPublisherUrl());
  csvw->writeString(f->getLang());
  if (!f->getStartDate().empty())
    csvw->writeInt(f->getStartDate().getYYYYMMDD());
  else
    csvw->skip();
  if (!f->getEndDate().empty())
    csvw->writeInt(f->getEndDate().getYYYYMMDD());
  else
    csvw->skip();
  csvw->writeString(f->getVersion());
  csvw->writeString(f->getContactEmail());
  csvw->writeString(f->getContactUrl());
  csvw->writeString(f->getDefaultLang());
  csvw->flushLine();
}

// ____________________________________________________________________________
void Writer::writePathways(gtfs::Feed* sourceFeed, std::ostream* os) const {
  Parser p(sourceFeed->getPath());
  auto csvp = p.getCsvParser("pathways.txt");
  ad::cppgtfs::Writer w;

  auto csvw = ad::cppgtfs::Writer::getPathwayCsvw(os);
  csvw->flushLine();

  ad::cppgtfs::gtfs::flat::Pathway fa;
  auto flds = Parser::getPathwayFlds(csvp.get());

  while (p.nextPathway(csvp.get(), &fa, flds)) {
    w.writePathway(fa, csvw.get());
  }
}

// ____________________________________________________________________________
void Writer::writeLevels(gtfs::Feed* sourceFeed, std::ostream* os) const {
  Parser p(sourceFeed->getPath());
  auto csvp = p.getCsvParser("levels.txt");
  ad::cppgtfs::Writer w;

  auto csvw = ad::cppgtfs::Writer::getLevelCsvw(os);
  csvw->flushLine();

  ad::cppgtfs::gtfs::flat::Level fa;
  auto flds = Parser::getLevelFlds(csvp.get());

  while (p.nextLevel(csvp.get(), &fa, flds)) {
    w.writeLevel(fa, csvw.get());
  }
}

// ____________________________________________________________________________
void Writer::writeAgency(gtfs::Feed* sourceFeed, std::ostream* os) const {
  Parser p(sourceFeed->getPath());
  auto csvp = p.getCsvParser("agency.txt");

  ad::cppgtfs::Writer w;

  auto csvw =
      ad::cppgtfs::Writer::getAgencyCsvw(os, sourceFeed->getAgencyAddFlds());
  csvw->flushLine();

  ad::cppgtfs::gtfs::flat::Agency fa;
  auto flds = Parser::getAgencyFlds(csvp.get());

  while (p.nextAgency(csvp.get(), &fa, flds)) {
    w.writeAgency(fa, csvw.get(), sourceFeed->getAgencyAddFlds());
  }
}

// ____________________________________________________________________________
void Writer::writeStops(gtfs::Feed* sourceFeed, std::ostream* os) const {
  Parser p(sourceFeed->getPath());
  auto csvp = p.getCsvParser("stops.txt");
  ad::cppgtfs::Writer w;

  auto csvw =
      ad::cppgtfs::Writer::getStopsCsvw(os, sourceFeed->getStopAddFlds());
  csvw->flushLine();

  ad::cppgtfs::gtfs::flat::Stop s;
  auto flds = Parser::getStopFlds(csvp.get());

  while (p.nextStop(csvp.get(), &s, flds)) {
    w.writeStop(s, csvw.get(), sourceFeed->getStopAddFlds());
  }
}

// ____________________________________________________________________________
void Writer::writeRoutes(gtfs::Feed* sourceFeed, std::ostream* os) const {
  ad::cppgtfs::Writer w;

  auto csvw =
      ad::cppgtfs::Writer::getRoutesCsvw(os, sourceFeed->getRouteAddFlds());
  csvw->flushLine();

  for (auto r : sourceFeed->getRoutes()) {
    w.writeRoute(r.second->getFlat(), csvw.get(),
                 sourceFeed->getRouteAddFlds());
  }
}

// ____________________________________________________________________________
void Writer::writeCalendar(gtfs::Feed* sourceFeed, std::ostream* os) const {
  Parser p(sourceFeed->getPath());
  auto csvp = p.getCsvParser("calendar.txt");
  ad::cppgtfs::Writer w;

  auto csvw = ad::cppgtfs::Writer::getCalendarCsvw(os);
  csvw->flushLine();

  ad::cppgtfs::gtfs::flat::Calendar c;
  auto flds = Parser::getCalendarFlds(csvp.get());

  while (p.nextCalendar(csvp.get(), &c, flds)) {
    w.writeCalendar(c, csvw.get());
  }
}

// ____________________________________________________________________________
void Writer::writeCalendarDates(gtfs::Feed* sourceFeed,
                                std::ostream* os) const {
  Parser p(sourceFeed->getPath());
  auto csvp = p.getCsvParser("calendar_dates.txt");
  ad::cppgtfs::Writer w;

  auto csvw = ad::cppgtfs::Writer::getCalendarDatesCsvw(os);
  csvw->flushLine();

  ad::cppgtfs::gtfs::flat::CalendarDate c;
  auto flds = Parser::getCalendarDateFlds(csvp.get());

  while (p.nextCalendarDate(csvp.get(), &c, flds)) {
    w.writeCalendarDate(c, csvw.get());
  }
}

// ____________________________________________________________________________
void Writer::writeFrequencies(gtfs::Feed* sourceFeed, std::ostream* os) const {
  Parser p(sourceFeed->getPath());
  auto csvp = p.getCsvParser("frequencies.txt");
  ad::cppgtfs::Writer w;

  auto csvw = ad::cppgtfs::Writer::getFrequencyCsvw(os);
  csvw->flushLine();

  ad::cppgtfs::gtfs::flat::Frequency f;
  auto flds = Parser::getFrequencyFlds(csvp.get());

  while (p.nextFrequency(csvp.get(), &f, flds)) {
    w.writeFrequency(f, csvw.get());
  }
}

// ____________________________________________________________________________
void Writer::writeTransfers(gtfs::Feed* sourceFeed, std::ostream* os) const {
  Parser p(sourceFeed->getPath());
  auto csvp = p.getCsvParser("transfers.txt");
  ad::cppgtfs::Writer w;

  auto csvw = ad::cppgtfs::Writer::getTransfersCsvw(os);
  csvw->flushLine();

  ad::cppgtfs::gtfs::flat::Transfer t;
  auto flds = Parser::getTransfersFlds(csvp.get());

  while (p.nextTransfer(csvp.get(), &t, flds)) {
    w.writeTransfer(t, csvw.get());
  }
}

// ____________________________________________________________________________
void Writer::writeFares(gtfs::Feed* sourceFeed, std::ostream* os) const {
  Parser p(sourceFeed->getPath());
  auto csvp = p.getCsvParser("fare_attributes.txt");
  ad::cppgtfs::Writer w;

  auto csvw = ad::cppgtfs::Writer::getFaresCsvw(os);
  csvw->flushLine();

  ad::cppgtfs::gtfs::flat::Fare f;
  auto flds = Parser::getFareFlds(csvp.get());

  while (p.nextFare(csvp.get(), &f, flds)) {
    w.writeFare(f, csvw.get());
  }
}

// ____________________________________________________________________________
void Writer::writeFareRules(gtfs::Feed* sourceFeed, std::ostream* os) const {
  Parser p(sourceFeed->getPath());
  auto csvp = p.getCsvParser("fare_rules.txt");
  ad::cppgtfs::Writer w;

  auto csvw = ad::cppgtfs::Writer::getFareRulesCsvw(os);
  csvw->flushLine();

  ad::cppgtfs::gtfs::flat::FareRule f;
  auto flds = Parser::getFareRuleFlds(csvp.get());

  while (p.nextFareRule(csvp.get(), &f, flds)) {
    w.writeFareRule(f, csvw.get());
  }
}

// ____________________________________________________________________________
void Writer::writeShapes(gtfs::Feed* sourceFeed, std::ostream* os) const {
  auto csvw = ad::cppgtfs::Writer::getShapesCsvw(os);
  csvw->flushLine();
  ad::cppgtfs::gtfs::flat::ShapePoint sp;
  ad::cppgtfs::Writer w;

  Parser p(sourceFeed->getPath());
  auto csvp = p.getCsvParser("shapes.txt");

  if (csvp->isGood()) {
    auto flds = Parser::getShapeFlds(csvp.get());

    std::string curShapeId;
    std::string curSkipShapeId;

    while (p.nextShapePoint(csvp.get(), &sp, flds)) {
      if (sp.id == curSkipShapeId) continue;
      if (sp.id != curShapeId) {
        if (sourceFeed->getShapes().has(sp.id)) {
          curShapeId = sp.id;
        } else {
          curSkipShapeId = sp.id;
          continue;
        }
      }

      w.writeShapePoint(sp, csvw.get());
    }
  }

  sourceFeed->getShapes().open();
  while (sourceFeed->getShapes().nextStoragePt(&sp)) {
    w.writeShapePoint(sp, csvw.get());
  }
}

// ____________________________________________________________________________
bool Writer::writeTrips(gtfs::Feed* sourceFeed, std::ostream* os) const {
  ad::cppgtfs::Writer w;
  bool hasFreqs = false;

  auto csvw =
      ad::cppgtfs::Writer::getTripsCsvw(os, sourceFeed->getTripAddFlds());
  csvw->flushLine();

  for (auto t : sourceFeed->getTrips()) {
    if (t.getFrequencies().size()) hasFreqs = true;
    w.writeTrip(t.getFlat(), csvw.get(), sourceFeed->getTripAddFlds());
  }

  return hasFreqs;
}

// ____________________________________________________________________________
void Writer::writeStopTimes(gtfs::Feed* sourceFeed, std::ostream* os) const {
  Parser p(sourceFeed->getPath());
  auto csvp = p.getCsvParser("stop_times.txt");
  ad::cppgtfs::Writer w;

  auto csvw = ad::cppgtfs::Writer::getStopTimesCsvw(os);
  csvw->flushLine();

  ad::cppgtfs::gtfs::flat::StopTime st;
  auto flds = Parser::getStopTimeFlds(csvp.get());

  std::string curTripId;
  Trip* cur = 0;

  while (p.nextStopTime(csvp.get(), &st, flds)) {
    // we may have changed to distance field
    if (curTripId != st.trip) {
      cur = sourceFeed->getTrips().get(st.trip);
      curTripId = st.trip;
    }
    for (const auto& stN : cur->getStopTimes()) {
      if (stN.getSeq() == st.sequence)
        st.shapeDistTravelled = stN.getShapeDistanceTravelled();
    }

    w.writeStopTime(st, csvw.get());
  }
}

// ___________________________________________________________________________
void Writer::cannotWrite(const std::string& file) {
  std::stringstream ss;
  ss << "Could not write to file";
  throw ad::cppgtfs::WriterException(ss.str(), file);
}

// ___________________________________________________________________________
void Writer::cannotWrite(const std::string& file, const std::string& file2) {
  std::stringstream ss;
  ss << "(temporary file for " << file2 << ") Could not write to file";
  throw ad::cppgtfs::WriterException(ss.str(), file);
}

// ___________________________________________________________________________
#ifdef LIBZIP_FOUND
void Writer::moveIntoZip(zip* za, const std::string& sourcePath,
                         const std::string& targetPath) {
  zip_source_t* s;
  FILE* fp = fopen(sourcePath.c_str(), "r");
  if (fp == 0) {
    std::stringstream ss;
    ss << "(temporary file for " << targetPath << ") Could not open file";
    throw ad::cppgtfs::WriterException(ss.str(), sourcePath);
  }

  // immediately unlink
  unlink(sourcePath.c_str());

  if ((s = zip_source_filep(za, fp, 0, -1)) == 0 ||
      zip_file_add(za, targetPath.c_str(), s, ZIP_FL_ENC_UTF_8) < 0) {
    zip_source_free(s);
    cannotWrite(targetPath);
  }
}
#endif

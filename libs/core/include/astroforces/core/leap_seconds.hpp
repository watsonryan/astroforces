/**
 * @file leap_seconds.hpp
 * @brief Leap-second table loading and UTC->TAI offset utilities.
 * @author Watosn
 */
#pragma once

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace astroforces::core::leap_seconds {

struct Entry {
  double utc_epoch_s{};
  double tai_minus_utc_s{};
};

using Table = std::vector<Entry>;

inline const Table& default_table() {
  static const Table table = {
      {63072000.0, 10.0},     // 1972-01-01
      {78796800.0, 11.0},     // 1972-07-01
      {94694400.0, 12.0},     // 1973-01-01
      {126230400.0, 13.0},    // 1974-01-01
      {157766400.0, 14.0},    // 1975-01-01
      {189302400.0, 15.0},    // 1976-01-01
      {220924800.0, 16.0},    // 1977-01-01
      {252460800.0, 17.0},    // 1978-01-01
      {283996800.0, 18.0},    // 1979-01-01
      {315532800.0, 19.0},    // 1980-01-01
      {362793600.0, 20.0},    // 1981-07-01
      {394329600.0, 21.0},    // 1982-07-01
      {425865600.0, 22.0},    // 1983-07-01
      {489024000.0, 23.0},    // 1985-07-01
      {567993600.0, 24.0},    // 1988-01-01
      {631152000.0, 25.0},    // 1990-01-01
      {662688000.0, 26.0},    // 1991-01-01
      {709948800.0, 27.0},    // 1992-07-01
      {741484800.0, 28.0},    // 1993-07-01
      {773020800.0, 29.0},    // 1994-07-01
      {820454400.0, 30.0},    // 1996-01-01
      {867715200.0, 31.0},    // 1997-07-01
      {915148800.0, 32.0},    // 1999-01-01
      {1136073600.0, 33.0},   // 2006-01-01
      {1230768000.0, 34.0},   // 2009-01-01
      {1341100800.0, 35.0},   // 2012-07-01
      {1435708800.0, 36.0},   // 2015-07-01
      {1483228800.0, 37.0},   // 2017-01-01
  };
  return table;
}

inline std::string trim(std::string s) {
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.front())) != 0) {
    s.erase(s.begin());
  }
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back())) != 0) {
    s.pop_back();
  }
  return s;
}

inline bool load_table_from_file(const std::string& path, Table* out) {
  if (out == nullptr) {
    return false;
  }
  std::ifstream in(path);
  if (!in.is_open()) {
    return false;
  }

  Table parsed{};
  std::string line{};
  while (std::getline(in, line)) {
    line = trim(line);
    if (line.empty() || line[0] == '#') {
      continue;
    }

    for (char& c : line) {
      if (c == ',') {
        c = ' ';
      }
    }

    std::istringstream iss(line);
    Entry e{};
    if (!(iss >> e.utc_epoch_s >> e.tai_minus_utc_s)) {
      return false;
    }
    parsed.push_back(e);
  }

  if (parsed.empty()) {
    return false;
  }

  std::sort(parsed.begin(), parsed.end(), [](const Entry& a, const Entry& b) { return a.utc_epoch_s < b.utc_epoch_s; });
  *out = std::move(parsed);
  return true;
}

inline const Table& active_table() {
  static const Table table = []() {
    const char* env = std::getenv("ASTROFORCES_LEAP_SECONDS_FILE");
    if (env != nullptr && env[0] != '\0') {
      Table loaded{};
      if (load_table_from_file(std::string(env), &loaded)) {
        return loaded;
      }
    }
    return default_table();
  }();
  return table;
}

inline double tai_minus_utc_seconds(const double utc_seconds, const Table& table) {
  double tai_utc = table.front().tai_minus_utc_s;
  for (const Entry& e : table) {
    if (utc_seconds >= e.utc_epoch_s) {
      tai_utc = e.tai_minus_utc_s;
    } else {
      break;
    }
  }
  return tai_utc;
}

inline double tai_minus_utc_seconds(const double utc_seconds) {
  return tai_minus_utc_seconds(utc_seconds, active_table());
}

}  // namespace astroforces::core::leap_seconds

/**
 * @file cip.hpp
 * @brief CIP (X, Y, s) loading and interpolation for rigorous GCRF/ITRF transforms.
 * @author Watosn
 */
#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

#include "astroforces/core/transforms.hpp"

namespace astroforces::core::cip {

enum class AngleUnit { Radians, Arcseconds };

struct Record {
  double mjd_utc{};
  CelestialIntermediatePole cip{};
};

class Series {
 public:
  static Series load_table(
      const std::filesystem::path& path,
      AngleUnit unit = AngleUnit::Radians,
      bool first_column_is_jd_utc = false) {
    Series series;
    std::ifstream in(path);
    if (!in.is_open()) {
      return series;
    }

    const double scale = (unit == AngleUnit::Radians) ? 1.0 : constants::kArcsecToRad;
    std::string line;
    while (std::getline(in, line)) {
      const auto rec = parse_line(line, scale, first_column_is_jd_utc);
      if (rec.has_value()) {
        series.records_.push_back(*rec);
      }
    }

    std::sort(series.records_.begin(), series.records_.end(), [](const Record& a, const Record& b) {
      return a.mjd_utc < b.mjd_utc;
    });
    series.records_.erase(std::unique(series.records_.begin(), series.records_.end(), [](const Record& a, const Record& b) {
      return a.mjd_utc == b.mjd_utc;
    }), series.records_.end());
    return series;
  }

  [[nodiscard]] bool empty() const noexcept { return records_.empty(); }
  [[nodiscard]] std::size_t size() const noexcept { return records_.size(); }
  [[nodiscard]] const std::vector<Record>& records() const noexcept { return records_; }

  [[nodiscard]] std::optional<CelestialIntermediatePole> at_mjd_utc(double mjd_utc) const noexcept {
    if (records_.empty()) {
      return std::nullopt;
    }
    if (records_.size() == 1) {
      return records_.front().cip;
    }
    if (mjd_utc <= records_.front().mjd_utc) {
      return records_.front().cip;
    }
    if (mjd_utc >= records_.back().mjd_utc) {
      return records_.back().cip;
    }

    const auto it = std::lower_bound(records_.begin(), records_.end(), mjd_utc, [](const Record& r, const double key) {
      return r.mjd_utc < key;
    });
    if (it == records_.end()) {
      return records_.back().cip;
    }
    if (it == records_.begin()) {
      return it->cip;
    }
    if (it->mjd_utc == mjd_utc) {
      return it->cip;
    }

    const Record& r1 = *it;
    const Record& r0 = *(it - 1);
    const double dt = r1.mjd_utc - r0.mjd_utc;
    if (dt <= 0.0) {
      return r0.cip;
    }
    const double alpha = (mjd_utc - r0.mjd_utc) / dt;

    CelestialIntermediatePole out{};
    out.x_rad = r0.cip.x_rad + alpha * (r1.cip.x_rad - r0.cip.x_rad);
    out.y_rad = r0.cip.y_rad + alpha * (r1.cip.y_rad - r0.cip.y_rad);
    out.s_rad = r0.cip.s_rad + alpha * (r1.cip.s_rad - r0.cip.s_rad);
    return out;
  }

  [[nodiscard]] std::optional<CelestialIntermediatePole> at_jd_utc(double jd_utc) const noexcept {
    return at_mjd_utc(jd_utc - 2400000.5);
  }

  [[nodiscard]] std::optional<CelestialIntermediatePole> at_utc_seconds(double utc_seconds) const noexcept {
    return at_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  }

 private:
  static std::optional<Record> parse_line(const std::string& line, double angle_scale, bool first_column_is_jd_utc) {
    std::string trimmed = line;
    if (trimmed.empty() || trimmed[0] == '#') {
      return std::nullopt;
    }
    for (char& c : trimmed) {
      if (c == ',') {
        c = ' ';
      }
    }
    std::istringstream iss(trimmed);
    std::vector<std::string> tokens;
    for (std::string tok; iss >> tok;) {
      tokens.push_back(tok);
    }
    if (tokens.size() < 4) {
      return std::nullopt;
    }

    try {
      const double t0 = std::stod(tokens[0]);
      const double x = std::stod(tokens[1]);
      const double y = std::stod(tokens[2]);
      const double s = std::stod(tokens[3]);
      Record rec{};
      rec.mjd_utc = first_column_is_jd_utc ? (t0 - 2400000.5) : t0;
      rec.cip.x_rad = x * angle_scale;
      rec.cip.y_rad = y * angle_scale;
      rec.cip.s_rad = s * angle_scale;
      if (!std::isfinite(rec.mjd_utc)) {
        return std::nullopt;
      }
      return rec;
    } catch (...) {
      return std::nullopt;
    }
  }

  std::vector<Record> records_{};
};

}  // namespace astroforces::core::cip

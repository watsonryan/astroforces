/**
 * @file eop.hpp
 * @brief IERS EOP loading and interpolation for frame transforms.
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

namespace astroforces::core::eop {

struct Record {
  double mjd_utc{};
  EarthOrientation eop{};
  bool predicted{};
};

class Series {
 public:
  static Series load_iers_finals(const std::filesystem::path& path) {
    Series series;
    std::ifstream in(path);
    if (!in.is_open()) {
      return series;
    }

    std::string line;
    while (std::getline(in, line)) {
      const auto rec = parse_iers_finals_line(line);
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

  [[nodiscard]] std::optional<EarthOrientation> at_mjd_utc(double mjd_utc) const noexcept {
    if (records_.empty()) {
      return std::nullopt;
    }
    if (records_.size() == 1) {
      return records_.front().eop;
    }
    if (mjd_utc <= records_.front().mjd_utc) {
      return records_.front().eop;
    }
    if (mjd_utc >= records_.back().mjd_utc) {
      return records_.back().eop;
    }

    const auto it = std::lower_bound(records_.begin(), records_.end(), mjd_utc, [](const Record& r, const double key) {
      return r.mjd_utc < key;
    });
    if (it == records_.end()) {
      return records_.back().eop;
    }
    if (it == records_.begin()) {
      return it->eop;
    }
    if (it->mjd_utc == mjd_utc) {
      return it->eop;
    }

    const Record& r1 = *it;
    const Record& r0 = *(it - 1);
    const double dt = r1.mjd_utc - r0.mjd_utc;
    if (dt <= 0.0) {
      return r0.eop;
    }
    const double alpha = (mjd_utc - r0.mjd_utc) / dt;

    EarthOrientation out{};
    out.xp_rad = r0.eop.xp_rad + alpha * (r1.eop.xp_rad - r0.eop.xp_rad);
    out.yp_rad = r0.eop.yp_rad + alpha * (r1.eop.yp_rad - r0.eop.yp_rad);
    out.dut1_s = r0.eop.dut1_s + alpha * (r1.eop.dut1_s - r0.eop.dut1_s);
    out.lod_s = r0.eop.lod_s + alpha * (r1.eop.lod_s - r0.eop.lod_s);
    out.dX_rad = r0.eop.dX_rad + alpha * (r1.eop.dX_rad - r0.eop.dX_rad);
    out.dY_rad = r0.eop.dY_rad + alpha * (r1.eop.dY_rad - r0.eop.dY_rad);
    return out;
  }

  [[nodiscard]] std::optional<EarthOrientation> at_jd_utc(double jd_utc) const noexcept {
    return at_mjd_utc(jd_utc - 2400000.5);
  }

  [[nodiscard]] std::optional<EarthOrientation> at_utc_seconds(double utc_seconds) const noexcept {
    return at_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  }

 private:
  static double parse_field(const std::string& line, const int start, const int width) {
    if (static_cast<int>(line.size()) < start + width) {
      return 0.0;
    }
    try {
      return std::stod(line.substr(static_cast<std::size_t>(start), static_cast<std::size_t>(width)));
    } catch (...) {
      return 0.0;
    }
  }

  static double mas_to_rad(const double mas) { return (mas * 1e-3) * constants::kArcsecToRad; }

  static std::optional<Record> parse_iers_finals_line(const std::string& line) {
    if (line.size() < 78) {
      return std::nullopt;
    }

    // Prefer token parsing for robustness across spacing variations in finals files.
    auto parse_from_tokens = [&]() -> std::optional<Record> {
      std::istringstream iss(line);
      std::vector<std::string> tokens;
      std::string tok;
      while (iss >> tok) {
        tokens.push_back(tok);
      }
      if (tokens.size() < 11) {
        return std::nullopt;
      }

      try {
        Record rec{};
        rec.mjd_utc = std::stod(tokens[3]);
        rec.eop.xp_rad = std::stod(tokens[5]) * constants::kArcsecToRad;
        rec.eop.yp_rad = std::stod(tokens[7]) * constants::kArcsecToRad;
        rec.eop.dut1_s = std::stod(tokens[10]);
        rec.eop.lod_s = (tokens.size() > 12) ? std::stod(tokens[12]) * 1e-3 : 0.0;
        // IERS finals2000A stores dX/dY in mas in later columns; parse by fixed
        // columns because token layouts vary with spacing/flags.
        rec.eop.dX_rad = (line.size() >= 106) ? mas_to_rad(parse_field(line, 97, 9)) : 0.0;
        rec.eop.dY_rad = (line.size() >= 125) ? mas_to_rad(parse_field(line, 116, 9)) : 0.0;
        rec.predicted = (tokens.size() > 9 && (tokens[4] == "P" || tokens[9] == "P"));
        if (rec.mjd_utc <= 0.0) {
          return std::nullopt;
        }
        return rec;
      } catch (...) {
        return std::nullopt;
      }
    };

    if (const auto token_rec = parse_from_tokens(); token_rec.has_value()) {
      return token_rec;
    }

    // Fallback: fixed-column parse for strict IERS layouts.
    const double mjd = parse_field(line, 7, 8);
    if (mjd <= 0.0) {
      return std::nullopt;
    }
    const double xp_arcsec = parse_field(line, 18, 9);
    const double yp_arcsec = parse_field(line, 37, 9);
    const double ut1_utc_s = parse_field(line, 58, 10);
    const double lod_s = (line.size() >= 93) ? parse_field(line, 79, 7) * 1e-3 : 0.0;

    Record rec{};
    rec.mjd_utc = mjd;
    rec.eop.xp_rad = xp_arcsec * constants::kArcsecToRad;
    rec.eop.yp_rad = yp_arcsec * constants::kArcsecToRad;
    rec.eop.dut1_s = ut1_utc_s;
    rec.eop.lod_s = lod_s;
    rec.eop.dX_rad = (line.size() >= 106) ? mas_to_rad(parse_field(line, 97, 9)) : 0.0;
    rec.eop.dY_rad = (line.size() >= 125) ? mas_to_rad(parse_field(line, 116, 9)) : 0.0;
    rec.predicted = (line.size() > 57 && (line[16] == 'P' || line[57] == 'P'));
    return rec;
  }

  std::vector<Record> records_{};
};

}  // namespace astroforces::core::eop

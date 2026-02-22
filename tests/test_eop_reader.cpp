/**
 * @file test_eop_reader.cpp
 * @brief EOP finals parser/interpolation tests.
 * @author Watosn
 */

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>

#include <spdlog/spdlog.h>

#include "astroforces/core/eop.hpp"

namespace {

bool approx(double a, double b, double tol) { return std::abs(a - b) <= tol; }

std::filesystem::path make_fixed_column_eop_line_file() {
  const auto path = std::filesystem::temp_directory_path() / "astroforces_eop_fixedcols_test.txt";
  std::ofstream out(path);
  std::string line(140, ' ');

  auto put_num = [&](const int start, const int width, const double value, const int prec) {
    std::ostringstream oss;
    oss << std::fixed << std::setw(width) << std::setprecision(prec) << value;
    const std::string s = oss.str();
    line.replace(static_cast<std::size_t>(start), static_cast<std::size_t>(width), s);
  };

  // Match parser field positions in astroforces::core::eop::Series.
  put_num(7, 8, 60680.0, 2);      // MJD
  put_num(18, 9, 0.123456, 6);    // xp arcsec
  put_num(37, 9, -0.234567, 6);   // yp arcsec
  put_num(58, 10, 0.1111111, 7);  // UT1-UTC [s]
  put_num(79, 7, 1.234, 3);       // LOD [ms]
  put_num(97, 9, 0.345678, 6);    // dX [mas]
  put_num(116, 9, -0.456789, 6);  // dY [mas]

  out << line << "\n";
  return path;
}

}  // namespace

int main() {
  using namespace astroforces::core;
  using namespace astroforces::core::eop;

#ifndef ASTROFORCES_SOURCE_DIR
  spdlog::error("ASTROFORCES_SOURCE_DIR missing");
  return 10;
#endif

  const std::filesystem::path sample = std::filesystem::path(ASTROFORCES_SOURCE_DIR) / "tests" / "data" / "eop_finals_sample.txt";
  const auto series = Series::load_iers_finals(sample);
  if (series.empty() || series.size() != 3) {
    spdlog::error("unexpected EOP series size");
    return 1;
  }

  const auto first = series.at_mjd_utc(60676.0);
  if (!first.has_value()) {
    spdlog::error("missing first EOP row");
    return 2;
  }
  if (!approx(first->xp_rad / astroforces::core::constants::kArcsecToRad, 0.123456, 1e-12)
      || !approx(first->yp_rad / astroforces::core::constants::kArcsecToRad, 0.234567, 1e-12)
      || !approx(first->dut1_s, 0.1, 1e-12)
      || !approx(first->lod_s, 0.0009, 1e-15)) {
    spdlog::error("first row values mismatch");
    return 3;
  }

  const auto mid = series.at_mjd_utc(60676.5);
  if (!mid.has_value()) {
    spdlog::error("missing interpolated row");
    return 4;
  }
  if (!approx(mid->xp_rad / astroforces::core::constants::kArcsecToRad, 0.173456, 1e-9)
      || !approx(mid->yp_rad / astroforces::core::constants::kArcsecToRad, 0.284567, 1e-9)
      || !approx(mid->dut1_s, 0.11, 1e-12)
      || !approx(mid->lod_s, 0.0010, 1e-15)) {
    spdlog::error("interpolation mismatch");
    return 5;
  }

  const auto clamped = series.at_mjd_utc(60690.0);
  if (!clamped.has_value() || !approx(clamped->xp_rad / astroforces::core::constants::kArcsecToRad, 0.323456, 1e-12)) {
    spdlog::error("clamp-to-end mismatch");
    return 6;
  }

  const auto fixed_cols = make_fixed_column_eop_line_file();
  const auto fixed_series = Series::load_iers_finals(fixed_cols);
  std::error_code ec;
  std::filesystem::remove(fixed_cols, ec);
  if (fixed_series.empty()) {
    spdlog::error("fixed-column EOP parse failed");
    return 7;
  }
  const auto fx = fixed_series.at_mjd_utc(60680.0);
  if (!fx.has_value()) {
    spdlog::error("missing fixed-column EOP row");
    return 8;
  }
  const double dX_mas = fx->dX_rad / astroforces::core::constants::kArcsecToRad * 1000.0;
  const double dY_mas = fx->dY_rad / astroforces::core::constants::kArcsecToRad * 1000.0;
  if (!approx(dX_mas, 0.345678, 1e-9) || !approx(dY_mas, -0.456789, 1e-9)) {
    spdlog::error("dX/dY parse mismatch in fixed-column EOP row");
    return 9;
  }

  return 0;
}

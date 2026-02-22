/**
 * @file test_cip_reader.cpp
 * @brief CIP X/Y/s parser/interpolation tests.
 * @author Watosn
 */

#include <cmath>
#include <filesystem>
#include <fstream>

#include <spdlog/spdlog.h>

#include "astroforces/core/cip.hpp"

namespace {

bool approx(double a, double b, double tol) { return std::abs(a - b) <= tol; }

}  // namespace

int main() {
  using namespace astroforces::core;
  namespace fs = std::filesystem;

  const fs::path p = fs::temp_directory_path() / "astroforces_test_cip.txt";
  {
    std::ofstream out(p);
    out << "# MJD_UTC X Y s (arcsec)\n";
    out << "60676.0  0.100  -0.200  0.010\n";
    out << "60677.0  0.120  -0.220  0.015\n";
    out << "60678.0  0.140  -0.240  0.020\n";
  }

  const auto series = cip::Series::load_table(p, cip::AngleUnit::Arcseconds, false);
  std::error_code ec;
  fs::remove(p, ec);

  if (series.empty() || series.size() != 3) {
    spdlog::error("unexpected CIP series size");
    return 1;
  }

  const auto first = series.at_mjd_utc(60676.0);
  if (!first.has_value()) {
    spdlog::error("missing first CIP row");
    return 2;
  }
  if (!approx(first->x_rad / constants::kArcsecToRad, 0.100, 1e-12)
      || !approx(first->y_rad / constants::kArcsecToRad, -0.200, 1e-12)
      || !approx(first->s_rad / constants::kArcsecToRad, 0.010, 1e-12)) {
    spdlog::error("first CIP row mismatch");
    return 3;
  }

  const auto mid = series.at_mjd_utc(60676.5);
  if (!mid.has_value()) {
    spdlog::error("missing interpolated CIP row");
    return 4;
  }
  if (!approx(mid->x_rad / constants::kArcsecToRad, 0.110, 1e-9)
      || !approx(mid->y_rad / constants::kArcsecToRad, -0.210, 1e-9)
      || !approx(mid->s_rad / constants::kArcsecToRad, 0.0125, 1e-9)) {
    spdlog::error("interpolated CIP mismatch");
    return 5;
  }

  return 0;
}

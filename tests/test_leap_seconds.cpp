/**
 * @file test_leap_seconds.cpp
 * @brief Leap-second provider parse/fallback tests.
 * @author Watson
 */

#include <cmath>
#include <filesystem>
#include <fstream>

#include <spdlog/spdlog.h>

#include "astroforces/core/leap_seconds.hpp"

namespace {

bool approx(double a, double b, double tol) { return std::abs(a - b) <= tol; }

}  // namespace

int main() {
  using namespace astroforces::core::leap_seconds;

  const auto& def = default_table();
  if (def.empty()) {
    spdlog::error("default leap-second table is empty");
    return 1;
  }

  const double tai_2016 = tai_minus_utc_seconds(1483228799.0, def);
  const double tai_2017 = tai_minus_utc_seconds(1483228800.0, def);
  if (!approx(tai_2016, 36.0, 1e-12) || !approx(tai_2017, 37.0, 1e-12)) {
    spdlog::error("default leap-second transition mismatch");
    return 2;
  }

  const auto path = std::filesystem::temp_directory_path() / "astroforces_test_leap_seconds.txt";
  {
    std::ofstream out(path);
    out << "# utc_epoch_s tai_minus_utc_s\n";
    out << "0 5\n";
    out << "100,6\n";
  }

  Table loaded{};
  if (!load_table_from_file(path.string(), &loaded)) {
    spdlog::error("failed to parse leap-second file");
    return 3;
  }

  std::error_code ec;
  std::filesystem::remove(path, ec);

  if (loaded.size() != 2 || !approx(tai_minus_utc_seconds(99.0, loaded), 5.0, 1e-12)
      || !approx(tai_minus_utc_seconds(101.0, loaded), 6.0, 1e-12)) {
    spdlog::error("parsed leap-second table values mismatch");
    return 4;
  }

  const auto& active = active_table();
  if (active.size() != def.size() || !approx(active.front().tai_minus_utc_s, def.front().tai_minus_utc_s, 1e-12)
      || !approx(active.back().tai_minus_utc_s, def.back().tai_minus_utc_s, 1e-12)) {
    spdlog::error("active table mismatch versus default table");
    return 5;
  }

  return 0;
}

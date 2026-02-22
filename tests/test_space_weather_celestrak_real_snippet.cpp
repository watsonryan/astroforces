/**
 * @file test_space_weather_celestrak_real_snippet.cpp
 * @brief Regression test from CelesTrak published format example row.
 * @author Watosn
 */

#include <cmath>
#include <filesystem>

#include <spdlog/spdlog.h>

#include "astroforces/atmo/types.hpp"
#include "astroforces/weather/celestrak_csv_provider.hpp"

namespace {

bool approx(double a, double b, double tol) { return std::abs(a - b) <= tol; }

}  // namespace

int main() {
  namespace fs = std::filesystem;
  const auto csv = fs::path(DRAGCPP_SOURCE_DIR) / "tests" / "data" / "celestrak_real_snippet.csv";
  const auto provider = astroforces::weather::CelesTrakCsvSpaceWeatherProvider::Create({.csv_file = csv});

  const astroforces::core::WeatherIndices w = provider->at(astroforces::core::Epoch{.utc_seconds = 946684800.0});  // 2000-01-01
  if (w.status != astroforces::core::Status::Ok) {
    spdlog::error("status");
    return 1;
  }
  if (!approx(w.f107, 130.1, 1e-12) || !approx(w.f107a, 131.0, 1e-12)) {
    spdlog::error("f107 mismatch");
    return 2;
  }
  if (!approx(w.kp, 3.0, 1e-12) || !approx(w.ap, 24.0, 1e-12)) {
    spdlog::error("daily geomag mismatch");
    return 3;
  }
  if (!approx(w.kp_3h_current, 4.3, 1e-12) || !approx(w.ap_3h_current, 56.0, 1e-12)) {
    spdlog::error("3h geomag mismatch");
    return 4;
  }
  if (!w.f107_observed || !w.geomagnetic_observed) {
    spdlog::error("quality flags mismatch");
    return 5;
  }
  return 0;
}

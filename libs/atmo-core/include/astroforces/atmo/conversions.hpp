/**
 * @file conversions.hpp
 * @brief Time and coordinate conversion helpers.
 * @author Watosn
 */
#pragma once

#include <cmath>
#include <ctime>
#include <utility>

#include "astroforces/atmo/constants.hpp"
#include "astroforces/atmo/types.hpp"

namespace astroforces::core {

inline GeodeticPoint spherical_geodetic_from_ecef(const Vec3& ecef_m) {
  constexpr double kPi = 3.1415926535897932384626433832795;
  const double r = norm(ecef_m);
  if (r <= 0.0) {
    return GeodeticPoint{};
  }
  const double lat = std::asin(ecef_m.z / r) * 180.0 / kPi;
  const double lon = std::atan2(ecef_m.y, ecef_m.x) * 180.0 / kPi;
  return GeodeticPoint{.lat_deg = lat, .lon_deg = lon, .alt_m = r - constants::kEarthRadiusWgs84M};
}

inline std::pair<int, double> utc_seconds_to_iyd_sec(double utc_seconds) {
  const std::time_t tt = static_cast<std::time_t>(utc_seconds);
  std::tm tm_utc{};
#if defined(_WIN32)
  gmtime_s(&tm_utc, &tt);
#else
  gmtime_r(&tt, &tm_utc);
#endif
  const int year2 = (tm_utc.tm_year + 1900) % 100;
  const int doy = tm_utc.tm_yday + 1;
  const int iyd = year2 * 1000 + doy;
  const double sec = static_cast<double>(tm_utc.tm_hour * 3600 + tm_utc.tm_min * 60 + tm_utc.tm_sec);
  return {iyd, sec};
}

inline double local_solar_time_hours(double utc_seconds, double lon_deg) {
  const auto iyd_sec = utc_seconds_to_iyd_sec(utc_seconds);
  const double ut_h = iyd_sec.second / 3600.0;
  double stl = std::fmod(ut_h + lon_deg / 15.0, 24.0);
  if (stl < 0.0) {
    stl += 24.0;
  }
  return stl;
}

inline double utc_seconds_to_julian_date_utc(double utc_seconds) {
  return utc_seconds / 86400.0 + 2440587.5;
}

}  // namespace astroforces::core

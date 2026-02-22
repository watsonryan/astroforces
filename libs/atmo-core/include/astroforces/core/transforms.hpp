/**
 * @file transforms.hpp
 * @brief Shared time and coordinate transform helpers.
 * @author Watosn
 */
#pragma once

#include <cmath>
#include <ctime>
#include <utility>

#include "astroforces/atmo/constants.hpp"
#include "astroforces/atmo/types.hpp"
#include "astroforces/core/math_utils.hpp"
#include "astroforces/core/sofa_utils.hpp"

namespace astroforces::core {

struct RotationWithDerivative {
  Mat3 r{};
  Mat3 dr{};
};

inline Vec3 cross(const Vec3& a, const Vec3& b) {
  return vec_cross(a, b);
}

inline GeodeticPoint spherical_geodetic_from_ecef(const Vec3& ecef_m) {
  const double r = norm(ecef_m);
  if (r <= 0.0) {
    return GeodeticPoint{};
  }
  const double lat = std::asin(ecef_m.z / r) * 180.0 / constants::kPi;
  const double lon = std::atan2(ecef_m.y, ecef_m.x) * 180.0 / constants::kPi;
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
  return utc_seconds / constants::kSecondsPerDay + 2440587.5;
}

inline double earth_rotation_angle_rad(double jd_ut1) {
  return sofa::era00(jd_ut1);
}

inline double gmst_rad_from_jd_utc(double jd_utc) {
  const double t = (jd_utc - constants::kJ2000Jd) / 36525.0;
  double gmst_deg = 280.46061837 + 360.98564736629 * (jd_utc - constants::kJ2000Jd) + 0.000387933 * t * t
                    - (t * t * t) / 38710000.0;
  gmst_deg = std::fmod(gmst_deg, 360.0);
  if (gmst_deg < 0.0) {
    gmst_deg += 360.0;
  }
  return gmst_deg * constants::kDegToRad;
}

inline Vec3 rotate_z(double angle_rad, const Vec3& v) {
  return mat_vec(sofa::rot_z(angle_rad), v);
}

// Approximate GMST-only transform (no precession/nutation/polar motion). Use
// gcrf_to_itrf_* with EOP/CIP inputs for production-quality navigation.
inline Vec3 eci_to_ecef_position(const Vec3& r_eci_m, double utc_seconds) {
  const double gmst = gmst_rad_from_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  return rotate_z(gmst, r_eci_m);
}

// Approximate GMST-only transform (no precession/nutation/polar motion). Use
// itrf_to_gcrf_* with EOP/CIP inputs for production-quality navigation.
inline Vec3 ecef_to_eci_position(const Vec3& r_ecef_m, double utc_seconds) {
  const double gmst = gmst_rad_from_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  return rotate_z(-gmst, r_ecef_m);
}

inline Vec3 eci_to_ecef_velocity(const Vec3& r_eci_m, const Vec3& v_eci_mps, double utc_seconds) {
  const double gmst = gmst_rad_from_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  const Vec3 omega_cross_r_eci{-constants::kEarthRotationRateRadPerSec * r_eci_m.y,
                                constants::kEarthRotationRateRadPerSec * r_eci_m.x,
                                0.0};
  return rotate_z(gmst, v_eci_mps - omega_cross_r_eci);
}

inline Vec3 ecef_to_eci_velocity(const Vec3& r_ecef_m, const Vec3& v_ecef_mps, double utc_seconds) {
  const double gmst = gmst_rad_from_jd_utc(utc_seconds_to_julian_date_utc(utc_seconds));
  const Vec3 r_eci_m = rotate_z(-gmst, r_ecef_m);
  const Vec3 omega_cross_r_eci{-constants::kEarthRotationRateRadPerSec * r_eci_m.y,
                                constants::kEarthRotationRateRadPerSec * r_eci_m.x,
                                0.0};
  return rotate_z(-gmst, v_ecef_mps) + omega_cross_r_eci;
}

inline double tio_locator_sp_rad(const double jd_tt) {
  return sofa::sp00(jd_tt);
}

inline Mat3 c2i_from_xys(const CelestialIntermediatePole& cip) {
  return sofa::c2ixys(cip.x_rad, cip.y_rad, cip.s_rad);
}

inline Mat3 polar_motion_matrix(const EarthOrientation& eop, const double sp_rad) {
  return sofa::pom00(eop.xp_rad, eop.yp_rad, sp_rad);
}

inline RotationWithDerivative gcrf_to_itrf_rotation_with_derivative(
    const double jd_utc,
    const double jd_tt,
    const CelestialIntermediatePole& cip,
    const EarthOrientation& eop) {
  const double jd_ut1 = jd_utc + eop.dut1_s / constants::kSecondsPerDay;
  CelestialIntermediatePole corrected = cip;
  corrected.x_rad += eop.dX_rad;
  corrected.y_rad += eop.dY_rad;

  const Mat3 rc2i = c2i_from_xys(corrected);
  const double era = earth_rotation_angle_rad(jd_ut1);
  const Mat3 r3era = sofa::rot_z(era);
  const Mat3 rpom = polar_motion_matrix(eop, tio_locator_sp_rad(jd_tt));

  const double omega = constants::kEarthRotationRateRadPerSec / (1.0 + eop.lod_s / constants::kSecondsPerDay);
  Mat3 s{};
  s(0, 1) = 1.0;
  s(1, 0) = -1.0;
  const Mat3 dr3 = mat_mul(s, r3era);
  Mat3 dr3_scaled{};
  for (int i = 0; i < 9; ++i) {
    dr3_scaled.v[static_cast<std::size_t>(i)] = omega * dr3.v[static_cast<std::size_t>(i)];
  }

  RotationWithDerivative out{};
  out.r = sofa::c2tcio(rc2i, era, rpom);
  out.dr = mat_mul(rpom, mat_mul(dr3_scaled, rc2i));
  return out;
}

inline Vec3 gcrf_to_itrf_position(
    const Vec3& r_gcrf,
    const double jd_utc,
    const double jd_tt,
    const CelestialIntermediatePole& cip,
    const EarthOrientation& eop) {
  return mat_vec(gcrf_to_itrf_rotation_with_derivative(jd_utc, jd_tt, cip, eop).r, r_gcrf);
}

inline Vec3 itrf_to_gcrf_position(
    const Vec3& r_itrf,
    const double jd_utc,
    const double jd_tt,
    const CelestialIntermediatePole& cip,
    const EarthOrientation& eop) {
  const Mat3 rt = mat_transpose(gcrf_to_itrf_rotation_with_derivative(jd_utc, jd_tt, cip, eop).r);
  return mat_vec(rt, r_itrf);
}

inline Vec3 gcrf_to_itrf_velocity(
    const Vec3& r_gcrf,
    const Vec3& v_gcrf,
    const double jd_utc,
    const double jd_tt,
    const CelestialIntermediatePole& cip,
    const EarthOrientation& eop) {
  const RotationWithDerivative rd = gcrf_to_itrf_rotation_with_derivative(jd_utc, jd_tt, cip, eop);
  return mat_vec(rd.r, v_gcrf) + mat_vec(rd.dr, r_gcrf);
}

inline Vec3 itrf_to_gcrf_velocity(
    const Vec3& r_itrf,
    const Vec3& v_itrf,
    const double jd_utc,
    const double jd_tt,
    const CelestialIntermediatePole& cip,
    const EarthOrientation& eop) {
  const RotationWithDerivative rd = gcrf_to_itrf_rotation_with_derivative(jd_utc, jd_tt, cip, eop);
  const Mat3 rt = mat_transpose(rd.r);
  const Vec3 r_gcrf = mat_vec(rt, r_itrf);
  const Vec3 correction = mat_vec(rd.dr, r_gcrf);
  return mat_vec(rt, v_itrf - correction);
}

}  // namespace astroforces::core

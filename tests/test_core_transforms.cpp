/**
 * @file test_core_transforms.cpp
 * @brief Core frame transform regression tests.
 * @author Watosn
 */

#include <cmath>

#include <spdlog/spdlog.h>

#include "astroforces/core/transforms.hpp"

namespace {

bool approx_abs(double a, double b, double tol) { return std::abs(a - b) <= tol; }

}  // namespace

int main() {
  using namespace astroforces::core;

  const double jd_utc = 2460737.5;  // 2025-03-16T00:00:00 UTC
  const double utc_seconds = (jd_utc - 2440587.5) * astroforces::core::constants::kSecondsPerDay;
  const double jd_tt = utc_seconds_to_julian_date_tt(utc_seconds);
  const double jd_tdb = utc_seconds_to_julian_date_tdb(utc_seconds);

  const double tdb_minus_tt_s = (jd_tdb - jd_tt) * astroforces::core::constants::kSecondsPerDay;
  if (!approx_abs(tdb_minus_tt_s, 0.00142, 5e-5)) {
    spdlog::error("utc->tdb conversion mismatch");
    return 4;
  }
  if (!(std::abs(tdb_minus_tt_s) <= 0.001679)) {
    spdlog::error("tdb-tt offset out of expected analytic bounds");
    return 5;
  }

  const CelestialIntermediatePole cip{
      .x_rad = 0.2 * astroforces::core::constants::kArcsecToRad,
      .y_rad = -0.15 * astroforces::core::constants::kArcsecToRad,
      .s_rad = -0.01 * astroforces::core::constants::kArcsecToRad,
  };
  const EarthOrientation eop{
      .xp_rad = 0.12 * astroforces::core::constants::kArcsecToRad,
      .yp_rad = 0.19 * astroforces::core::constants::kArcsecToRad,
      .dut1_s = 0.083,
      .lod_s = 0.0009,
      .dX_rad = 0.01 * astroforces::core::constants::kArcsecToRad,
      .dY_rad = -0.01 * astroforces::core::constants::kArcsecToRad,
  };

  const auto rd = gcrf_to_itrf_rotation_with_derivative(jd_utc, jd_tt, cip, eop);
  const auto should_be_i = mat_mul(rd.r, mat_transpose(rd.r));
  if (!approx_abs(should_be_i(0, 0), 1.0, 1e-14) || !approx_abs(should_be_i(1, 1), 1.0, 1e-14)
      || !approx_abs(should_be_i(2, 2), 1.0, 1e-14) || !approx_abs(should_be_i(0, 1), 0.0, 1e-14)
      || !approx_abs(should_be_i(0, 2), 0.0, 1e-14) || !approx_abs(should_be_i(1, 2), 0.0, 1e-14)) {
    spdlog::error("rotation matrix orthogonality failure");
    return 1;
  }

  const Vec3 r_gcrf{7000e3, -1200e3, 2500e3};
  const Vec3 v_gcrf{1200.0, 6800.0, -900.0};
  const Vec3 r_itrf = gcrf_to_itrf_position(r_gcrf, jd_utc, jd_tt, cip, eop);
  const Vec3 v_itrf = gcrf_to_itrf_velocity(r_gcrf, v_gcrf, jd_utc, jd_tt, cip, eop);

  const Vec3 r_back = itrf_to_gcrf_position(r_itrf, jd_utc, jd_tt, cip, eop);
  const Vec3 v_back = itrf_to_gcrf_velocity(r_itrf, v_itrf, jd_utc, jd_tt, cip, eop);

  if (!approx_abs(r_back.x, r_gcrf.x, 1e-6) || !approx_abs(r_back.y, r_gcrf.y, 1e-6)
      || !approx_abs(r_back.z, r_gcrf.z, 1e-6)) {
    spdlog::error("position round-trip failure");
    return 2;
  }
  if (!approx_abs(v_back.x, v_gcrf.x, 1e-9) || !approx_abs(v_back.y, v_gcrf.y, 1e-9)
      || !approx_abs(v_back.z, v_gcrf.z, 1e-9)) {
    spdlog::error("velocity round-trip failure");
    return 3;
  }

  return 0;
}

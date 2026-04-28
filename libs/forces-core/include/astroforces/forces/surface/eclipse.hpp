/**
 * @file eclipse.hpp
 * @brief Eclipse/visibility helpers for solar-dependent surface forces.
 * @author Watson
 */
#pragma once

#include <algorithm>
#include <cmath>
#include <numbers>

#include "astroforces/core/constants.hpp"
#include "astroforces/core/types.hpp"

namespace astroforces::forces {

inline double circle_visible_fraction(double r_sun_ang, double r_occ_ang, double center_sep_ang) {
  if (!(r_sun_ang > 0.0) || !std::isfinite(r_sun_ang) || !std::isfinite(r_occ_ang) || !std::isfinite(center_sep_ang)) {
    return 1.0;
  }
  if (r_occ_ang <= 0.0) {
    return 1.0;
  }

  const double d = std::max(0.0, center_sep_ang);
  const double pi = std::numbers::pi;

  double overlap = 0.0;
  if (d >= (r_sun_ang + r_occ_ang)) {
    overlap = 0.0;
  } else if (d <= std::abs(r_sun_ang - r_occ_ang)) {
    const double smaller = std::min(r_sun_ang, r_occ_ang);
    overlap = pi * smaller * smaller;
  } else {
    const double d1 = (r_sun_ang * r_sun_ang - r_occ_ang * r_occ_ang + d * d) / (2.0 * d);
    const double d2 = (r_occ_ang * r_occ_ang - r_sun_ang * r_sun_ang + d * d) / (2.0 * d);
    const double t1 = std::max(0.0, r_sun_ang * r_sun_ang - d1 * d1);
    const double t2 = std::max(0.0, r_occ_ang * r_occ_ang - d2 * d2);
    overlap = r_sun_ang * r_sun_ang * std::acos(std::clamp(d1 / r_sun_ang, -1.0, 1.0)) - d1 * std::sqrt(t1) +
              r_occ_ang * r_occ_ang * std::acos(std::clamp(d2 / r_occ_ang, -1.0, 1.0)) - d2 * std::sqrt(t2);
  }

  const double sun_area = pi * r_sun_ang * r_sun_ang;
  if (!(sun_area > 0.0)) {
    return 1.0;
  }
  return std::clamp((sun_area - overlap) / sun_area, 0.0, 1.0);
}

inline double sun_visibility_factor(const astroforces::core::Vec3& r_sc_eci_m,
                                    const astroforces::core::Vec3& r_sun_eci_m,
                                    const astroforces::core::Vec3* r_moon_eci_m = nullptr) {
  struct Occulter {
    astroforces::core::Vec3 pos{};
    double radius_m{};
  };

  const Occulter occluders[2] = {
      Occulter{astroforces::core::Vec3{0.0, 0.0, 0.0}, astroforces::core::constants::kEarthRadiusWgs84M},
      Occulter{r_moon_eci_m ? *r_moon_eci_m : astroforces::core::Vec3{0.0, 0.0, 0.0}, astroforces::core::constants::kMoonRadiusM},
  };

  double min_visibility = 1.0;
  for (int i = 0; i < 2; ++i) {
    if (i == 1 && r_moon_eci_m == nullptr) {
      continue;
    }
    const auto& occ = occluders[i];
    const auto sun_to_sat = r_sc_eci_m - r_sun_eci_m;
    const auto sun_to_occ = occ.pos - r_sun_eci_m;
    const auto occ_to_sat = r_sc_eci_m - occ.pos;

    const double d_sun_sat = astroforces::core::norm(sun_to_sat);
    const double d_sun_occ = astroforces::core::norm(sun_to_occ);
    if (!(d_sun_sat > 0.0) || !(d_sun_occ > 0.0) || d_sun_sat < d_sun_occ) {
      continue;
    }

    const auto sun_to_sat_hat = sun_to_sat / d_sun_sat;
    const double occ_proj = astroforces::core::dot(occ_to_sat, sun_to_sat_hat);
    if (!(occ_proj > 0.0)) {
      continue;
    }

    const astroforces::core::Vec3 center_sep{
        occ_to_sat.y * sun_to_sat_hat.z - occ_to_sat.z * sun_to_sat_hat.y,
        occ_to_sat.z * sun_to_sat_hat.x - occ_to_sat.x * sun_to_sat_hat.z,
        occ_to_sat.x * sun_to_sat_hat.y - occ_to_sat.y * sun_to_sat_hat.x,
    };
    const double r_sun_ang = astroforces::core::constants::kSunRadiusM / d_sun_sat;
    const double r_occ_ang = occ.radius_m / occ_proj;
    const double sep_ang = astroforces::core::norm(center_sep) / occ_proj;

    const double vis = circle_visible_fraction(r_sun_ang, r_occ_ang, sep_ang);
    min_visibility = std::min(min_visibility, vis);
    if (min_visibility <= 0.0) {
      return 0.0;
    }
  }
  return std::clamp(min_visibility, 0.0, 1.0);
}

}  // namespace astroforces::forces

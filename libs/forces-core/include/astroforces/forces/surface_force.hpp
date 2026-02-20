/**
 * @file surface_force.hpp
 * @brief Shared surface-force kernel used by drag and SRP models.
 * @author Watosn
 */
#pragma once

#include <cmath>

#include "astroforces/sc/spacecraft.hpp"

namespace astroforces::forces {

struct SurfaceForceResult {
  astroforces::atmo::Vec3 acceleration_mps2{};
  double area_m2{};
  double coeff{};
  astroforces::atmo::Status status{astroforces::atmo::Status::Ok};
};

inline astroforces::atmo::Vec3 unit_direction(const astroforces::atmo::Vec3& v, double* norm_out = nullptr) {
  const double n = astroforces::atmo::norm(v);
  if (norm_out) {
    *norm_out = n;
  }
  if (n <= 0.0 || !std::isfinite(n)) {
    return astroforces::atmo::Vec3{};
  }
  return v / n;
}

inline SurfaceForceResult evaluate_surface_force(const astroforces::sc::SpacecraftProperties& sc,
                                                 const astroforces::atmo::Vec3& flow_dir_frame_unit,
                                                 const astroforces::atmo::Vec3& flow_dir_body_unit,
                                                 double pressure_pa,
                                                 double reference_coeff,
                                                 astroforces::sc::SurfaceCoeffModel coeff_model,
                                                 double direction_sign) {
  if (sc.mass_kg <= 0.0 || !std::isfinite(pressure_pa)) {
    return SurfaceForceResult{.status = astroforces::atmo::Status::InvalidInput};
  }

  const auto proj =
      astroforces::sc::projected_area_and_coeff(sc, flow_dir_body_unit, reference_coeff, coeff_model);
  const double coeff = proj.coeff_effective;
  const double area = proj.area_m2;
  const double scale = direction_sign * pressure_pa * coeff * area / sc.mass_kg;

  return SurfaceForceResult{
      .acceleration_mps2 = scale * flow_dir_frame_unit,
      .area_m2 = area,
      .coeff = coeff,
      .status = astroforces::atmo::Status::Ok,
  };
}

}  // namespace astroforces::forces


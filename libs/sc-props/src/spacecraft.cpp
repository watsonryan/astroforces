/**
 * @file spacecraft.cpp
 * @brief Spacecraft geometry helper implementation.
 * @author Watosn
 */

#include "astroforces/sc/spacecraft.hpp"

#include <algorithm>

namespace astroforces::sc {

ForceProjection projected_area_and_coeff(const SpacecraftProperties& sc,
                                         const astroforces::core::Vec3& flow_dir_body,
                                         double reference_coeff,
                                         SurfaceCoeffModel coeff_model) {
  if (!sc.use_surface_model || sc.surfaces.empty()) {
    return ForceProjection{.area_m2 = sc.reference_area_m2, .coeff_effective = reference_coeff};
  }

  double area = 0.0;
  double weighted_coeff = 0.0;
  for (const auto& s : sc.surfaces) {
    const double c = -astroforces::core::dot(s.normal_body, flow_dir_body);
    const double proj = s.area_m2 * std::max(0.0, c);
    area += proj;
    double coeff_surface = reference_coeff;
    if (coeff_model == SurfaceCoeffModel::Drag) {
      const double cd_base = (s.cd > 0.0) ? s.cd : sc.cd;
      coeff_surface = cd_base;
      if (s.specularity > 0.0 || s.accommodation > 0.0) {
        const double incidence = std::clamp(c, 0.0, 1.0);
        const double spec = std::clamp(s.specularity, 0.0, 1.0);
        const double accom = std::clamp(s.accommodation, 0.0, 1.0);
        const double modifier = 1.0 + 0.25 * accom * (1.0 + incidence) - 0.15 * spec * (1.0 - incidence);
        coeff_surface = cd_base * modifier;
      }
    } else if (coeff_model == SurfaceCoeffModel::RadiationPressure) {
      coeff_surface = (s.cr > 0.0) ? s.cr : sc.cr;
    }
    weighted_coeff += coeff_surface * proj;
  }
  const double coeff_effective = (area > 0.0) ? (weighted_coeff / area) : reference_coeff;
  return ForceProjection{.area_m2 = area, .coeff_effective = coeff_effective};
}

AeroProjection projected_area_and_cd(const SpacecraftProperties& sc, const astroforces::core::Vec3& flow_dir_body) {
  const auto p = projected_area_and_coeff(sc, flow_dir_body, sc.cd, SurfaceCoeffModel::Drag);
  return AeroProjection{.area_m2 = p.area_m2, .cd_effective = p.coeff_effective};
}

ForceProjection projected_area_and_cr(const SpacecraftProperties& sc, const astroforces::core::Vec3& flow_dir_body) {
  return projected_area_and_coeff(sc, flow_dir_body, sc.cr, SurfaceCoeffModel::RadiationPressure);
}

double projected_area_m2(const SpacecraftProperties& sc, const astroforces::core::Vec3& flow_dir_body) {
  return projected_area_and_cd(sc, flow_dir_body).area_m2;
}

}  // namespace astroforces::sc

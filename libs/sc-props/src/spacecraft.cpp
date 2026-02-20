/**
 * @file spacecraft.cpp
 * @brief Spacecraft geometry helper implementation.
 * @author Watosn
 */

#include "dragcpp/sc/spacecraft.hpp"

#include <algorithm>

namespace dragcpp::sc {

AeroProjection projected_area_and_cd(const SpacecraftProperties& sc, const dragcpp::atmo::Vec3& flow_dir_body) {
  if (!sc.use_surface_model || sc.surfaces.empty()) {
    return AeroProjection{.area_m2 = sc.reference_area_m2, .cd_effective = sc.cd};
  }

  double area = 0.0;
  double weighted_cd = 0.0;
  for (const auto& s : sc.surfaces) {
    const double c = -dragcpp::atmo::dot(s.normal_body, flow_dir_body);
    const double proj = s.area_m2 * std::max(0.0, c);
    area += proj;
    const double cd_surface = (s.cd > 0.0) ? s.cd : sc.cd;
    weighted_cd += cd_surface * proj;
  }
  const double cd_effective = (area > 0.0) ? (weighted_cd / area) : sc.cd;
  return AeroProjection{.area_m2 = area, .cd_effective = cd_effective};
}

double projected_area_m2(const SpacecraftProperties& sc, const dragcpp::atmo::Vec3& flow_dir_body) {
  return projected_area_and_cd(sc, flow_dir_body).area_m2;
}

}  // namespace dragcpp::sc

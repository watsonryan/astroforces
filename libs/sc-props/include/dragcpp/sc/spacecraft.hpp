/**
 * @file spacecraft.hpp
 * @brief Spacecraft geometry and drag-relevant properties.
 * @author Watosn
 */
#pragma once

#include <vector>

#include "dragcpp/atmo/types.hpp"

namespace dragcpp::sc {

struct Surface {
  dragcpp::atmo::Vec3 normal_body{};
  double area_m2{};
  double cd{};
  double specularity{};
  double accommodation{};
};

struct SpacecraftProperties {
  double mass_kg{};
  double reference_area_m2{};
  double cd{};
  bool use_surface_model{};
  std::vector<Surface> surfaces{};
};

struct AeroProjection {
  double area_m2{};
  double cd_effective{};
};

[[nodiscard]] double projected_area_m2(const SpacecraftProperties& sc,
                                       const dragcpp::atmo::Vec3& flow_dir_body);
[[nodiscard]] AeroProjection projected_area_and_cd(const SpacecraftProperties& sc,
                                                   const dragcpp::atmo::Vec3& flow_dir_body);

}  // namespace dragcpp::sc

/**
 * @file spacecraft.hpp
 * @brief Spacecraft geometry and drag-relevant properties.
 * @author Watosn
 */
#pragma once

#include <vector>

#include "astroforces/atmo/types.hpp"

namespace astroforces::sc {

struct Surface {
  astroforces::core::Vec3 normal_body{};
  double area_m2{};
  double cd{};
  double cr{1.0};
  double specularity{};
  double accommodation{};
};

struct SpacecraftProperties {
  double mass_kg{};
  double reference_area_m2{};
  double cd{};
  double cr{1.0};
  bool use_surface_model{};
  std::vector<Surface> surfaces{};
};

struct AeroProjection {
  double area_m2{};
  double cd_effective{};
};

enum class SurfaceCoeffModel : unsigned char { Drag, RadiationPressure };

struct ForceProjection {
  double area_m2{};
  double coeff_effective{};
};

[[nodiscard]] double projected_area_m2(const SpacecraftProperties& sc,
                                       const astroforces::core::Vec3& flow_dir_body);
[[nodiscard]] AeroProjection projected_area_and_cd(const SpacecraftProperties& sc,
                                                   const astroforces::core::Vec3& flow_dir_body);
[[nodiscard]] ForceProjection projected_area_and_coeff(const SpacecraftProperties& sc,
                                                       const astroforces::core::Vec3& flow_dir_body,
                                                       double reference_coeff,
                                                       SurfaceCoeffModel coeff_model);
[[nodiscard]] ForceProjection projected_area_and_cr(const SpacecraftProperties& sc,
                                                    const astroforces::core::Vec3& flow_dir_body);

}  // namespace astroforces::sc

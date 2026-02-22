/**
 * @file spacecraft.hpp
 * @brief Spacecraft geometry and drag-relevant properties.
 * @author Watosn
 */
#pragma once

#include <vector>

#include "astroforces/core/types.hpp"

namespace astroforces::sc {

/**
 * @brief One spacecraft surface panel in BODY frame.
 */
struct Surface {
  astroforces::core::Vec3 normal_body{};
  double area_m2{};
  double cd{};
  double cr{1.0};
  double specularity{};
  double accommodation{};
};

/**
 * @brief Spacecraft mass and aerodynamic/radiation properties.
 */
struct SpacecraftProperties {
  double mass_kg{};
  double reference_area_m2{};
  double cd{};
  double cr{1.0};
  bool use_surface_model{};
  std::vector<Surface> surfaces{};
};

/**
 * @brief Projected aerodynamic area and effective drag coefficient.
 */
struct AeroProjection {
  double area_m2{};
  double cd_effective{};
};

/**
 * @brief Coefficient family used for surface-force aggregation.
 */
enum class SurfaceCoeffModel : unsigned char { Drag, RadiationPressure };

/**
 * @brief Projected area and effective scalar coefficient.
 */
struct ForceProjection {
  double area_m2{};
  double coeff_effective{};
};

/**
 * @brief Project surface model onto flow direction and return total area.
 */
[[nodiscard]] double projected_area_m2(const SpacecraftProperties& sc,
                                       const astroforces::core::Vec3& flow_dir_body);
/**
 * @brief Project surface model onto flow direction with effective drag coefficient.
 */
[[nodiscard]] AeroProjection projected_area_and_cd(const SpacecraftProperties& sc,
                                                   const astroforces::core::Vec3& flow_dir_body);
/**
 * @brief General projection helper for drag/SRP/ERP-style scalar coefficients.
 */
[[nodiscard]] ForceProjection projected_area_and_coeff(const SpacecraftProperties& sc,
                                                       const astroforces::core::Vec3& flow_dir_body,
                                                       double reference_coeff,
                                                       SurfaceCoeffModel coeff_model);
/**
 * @brief Convenience wrapper for radiation pressure coefficient projection.
 */
[[nodiscard]] ForceProjection projected_area_and_cr(const SpacecraftProperties& sc,
                                                    const astroforces::core::Vec3& flow_dir_body);

}  // namespace astroforces::sc

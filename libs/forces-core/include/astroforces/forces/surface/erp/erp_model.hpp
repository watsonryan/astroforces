/**
 * @file erp_model.hpp
 * @brief Earth radiation pressure acceleration model.
 * @author Watosn
 */
#pragma once

#include "astroforces/atmo/constants.hpp"
#include "astroforces/atmo/types.hpp"
#include "astroforces/sc/spacecraft.hpp"

namespace astroforces::erp {

struct ErpResult {
  astroforces::core::Vec3 acceleration_mps2{};
  double earth_radiation_pressure_pa{};
  double earth_distance_m{};
  double area_m2{};
  double cr{};
  astroforces::core::Status status{astroforces::core::Status::Ok};
};

class ErpAccelerationModel final {
 public:
  struct Config {
    double earth_radiation_pressure_ref_pa{astroforces::core::constants::kEarthRadiationPressureAtEarthRadiusPa};
    double earth_reference_radius_m{astroforces::core::constants::kEarthRadiusWgs84M};
  };

  ErpAccelerationModel() = default;
  explicit ErpAccelerationModel(const Config& config) : config_(config) {}

  [[nodiscard]] ErpResult evaluate(const astroforces::core::StateVector& state,
                                   const astroforces::sc::SpacecraftProperties& sc) const;

 private:
  Config config_{};
};

}  // namespace astroforces::erp

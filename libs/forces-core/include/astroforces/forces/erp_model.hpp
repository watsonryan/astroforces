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
  astroforces::atmo::Vec3 acceleration_mps2{};
  double earth_radiation_pressure_pa{};
  double earth_distance_m{};
  double area_m2{};
  double cr{};
  astroforces::atmo::Status status{astroforces::atmo::Status::Ok};
};

class ErpAccelerationModel final {
 public:
  struct Config {
    double earth_radiation_pressure_ref_pa{astroforces::atmo::constants::kEarthRadiationPressureAtEarthRadiusPa};
    double earth_reference_radius_m{astroforces::atmo::constants::kEarthRadiusWgs84M};
  };

  ErpAccelerationModel() = default;
  explicit ErpAccelerationModel(const Config& config) : config_(config) {}

  [[nodiscard]] ErpResult evaluate(const astroforces::atmo::StateVector& state,
                                   const astroforces::sc::SpacecraftProperties& sc) const;

 private:
  Config config_{};
};

}  // namespace astroforces::erp

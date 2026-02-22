/**
 * @file earth_radiation_model.hpp
 * @brief Earth radiation pressure acceleration model.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>

#include "astroforces/atmo/constants.hpp"
#include "astroforces/atmo/types.hpp"
#include "astroforces/sc/spacecraft.hpp"

namespace jpl::eph {
class Ephemeris;
class Workspace;
}  // namespace jpl::eph

namespace astroforces::forces {

struct EarthRadiationResult {
  astroforces::core::Vec3 acceleration_mps2{};
  double earth_radiation_pressure_pa{};
  double albedo_pressure_pa{};
  double ir_pressure_pa{};
  double albedo_phase_function{};
  double albedo_eclipse_factor{1.0};
  double earth_distance_m{};
  double area_m2{};
  double cr{};
  astroforces::core::Status status{astroforces::core::Status::Ok};
};

class EarthRadiationAccelerationModel final {
 public:
  struct Config {
    std::filesystem::path ephemeris_file{};
    double earth_reference_radius_m{astroforces::core::constants::kEarthRadiusWgs84M};
    double solar_flux_w_m2{astroforces::core::constants::kSolarIrradianceAt1AuWm2};
    double earth_albedo{astroforces::core::constants::kEarthBondAlbedo};
    double earth_ir_flux_w_m2{astroforces::core::constants::kEarthIrFluxWm2};
    double speed_of_light_mps{astroforces::core::constants::kSpeedOfLightMps};
    bool use_albedo{true};
    bool use_earth_ir{true};
    bool use_eclipse{true};
  };

  static std::unique_ptr<EarthRadiationAccelerationModel> Create(const Config& config);

  EarthRadiationAccelerationModel() = default;
  explicit EarthRadiationAccelerationModel(const Config& config);

  [[nodiscard]] EarthRadiationResult evaluate(const astroforces::core::StateVector& state,
                                   const astroforces::sc::SpacecraftProperties& sc) const;

 private:
  Config config_{};
  std::shared_ptr<jpl::eph::Ephemeris> ephemeris_{};
  mutable std::shared_ptr<jpl::eph::Workspace> workspace_{};
};

}  // namespace astroforces::forces

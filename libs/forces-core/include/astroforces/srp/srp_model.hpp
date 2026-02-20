/**
 * @file srp_model.hpp
 * @brief Solar radiation pressure acceleration model.
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
}

namespace astroforces::srp {

struct SrpResult {
  astroforces::atmo::Vec3 acceleration_mps2{};
  double solar_pressure_pa{};
  double sun_distance_m{};
  double area_m2{};
  double cr{};
  bool eclipsed{};
  astroforces::atmo::Status status{astroforces::atmo::Status::Ok};
};

class SrpAccelerationModel final {
 public:
  struct Config {
    std::filesystem::path ephemeris_file{};
    double solar_pressure_1au_pa{astroforces::atmo::constants::kSolarRadiationPressureAt1AuPa};
    double astronomical_unit_m{astroforces::atmo::constants::kAstronomicalUnitM};
    bool use_eclipse{false};
  };

  static std::unique_ptr<SrpAccelerationModel> Create(const Config& config);

  [[nodiscard]] SrpResult evaluate(const astroforces::atmo::StateVector& state,
                                   const astroforces::sc::SpacecraftProperties& sc) const;

 private:
  explicit SrpAccelerationModel(Config config) : config_(std::move(config)) {}

  Config config_{};
  std::shared_ptr<jpl::eph::Ephemeris> ephemeris_{};
  mutable std::shared_ptr<jpl::eph::Workspace> workspace_{};
};

}  // namespace astroforces::srp


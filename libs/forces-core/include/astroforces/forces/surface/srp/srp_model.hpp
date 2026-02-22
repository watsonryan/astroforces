/**
 * @file srp_model.hpp
 * @brief Solar radiation pressure acceleration model.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>

#include "astroforces/core/constants.hpp"
#include "astroforces/core/types.hpp"
#include "astroforces/sc/spacecraft.hpp"

namespace jpl::eph {
class Ephemeris;
class Workspace;
}

namespace astroforces::forces {

/**
 * @brief Output bundle for SRP evaluation.
 */
struct SrpResult {
  astroforces::core::Vec3 acceleration_mps2{};
  double solar_pressure_pa{};
  double sun_distance_m{};
  double eclipse_factor{1.0};
  double area_m2{};
  double cr{};
  bool eclipsed{};
  astroforces::core::Status status{astroforces::core::Status::Ok};
};

/**
 * @brief Solar radiation pressure acceleration model.
 */
class SrpAccelerationModel final {
 public:
  /**
   * @brief Configuration for SRP model construction.
   */
  struct Config {
    std::filesystem::path ephemeris_file{};
    double solar_pressure_1au_pa{astroforces::core::constants::kSolarRadiationPressureAt1AuPa};
    double astronomical_unit_m{astroforces::core::constants::kAstronomicalUnitM};
    bool use_eclipse{false};
  };

  /**
   * @brief Factory helper that loads ephemeris resources.
   */
  static std::unique_ptr<SrpAccelerationModel> Create(const Config& config);

  /**
   * @brief Evaluate SRP acceleration for one state and spacecraft.
   */
  [[nodiscard]] SrpResult evaluate(const astroforces::core::StateVector& state,
                                   const astroforces::sc::SpacecraftProperties& sc) const;

 private:
  explicit SrpAccelerationModel(Config config) : config_(std::move(config)) {}

  Config config_{};
  std::shared_ptr<jpl::eph::Ephemeris> ephemeris_{};
};

}  // namespace astroforces::forces

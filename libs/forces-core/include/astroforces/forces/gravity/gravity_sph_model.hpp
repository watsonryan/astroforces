/**
 * @file gravity_sph_model.hpp
 * @brief Full spherical harmonic gravity model (Cnm/Snm) with optional solid Earth tides.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>

#include "astroforces/atmo/constants.hpp"
#include "astroforces/atmo/types.hpp"

namespace jpl::eph {
class Ephemeris;
class Workspace;
}  // namespace jpl::eph

namespace astroforces::forces {

struct GravitySphResult {
  astroforces::core::Vec3 acceleration_mps2{};
  astroforces::core::Vec3 central_mps2{};
  astroforces::core::Vec3 sph_mps2{};
  astroforces::core::Vec3 solid_tide_sun_mps2{};
  astroforces::core::Vec3 solid_tide_moon_mps2{};
  astroforces::core::Status status{astroforces::core::Status::Ok};
};

class GravitySphAccelerationModel final {
 public:
  struct Config {
    std::filesystem::path gravity_model_file{};
    std::filesystem::path ephemeris_file{};
    int max_degree{360};
    bool use_central{true};
    bool use_sph{true};
    bool use_solid_earth_tides{true};
    bool use_sun_tide{true};
    bool use_moon_tide{true};
    bool convert_to_tide_free{true};
    bool use_simple_eci_to_ecef{true};
    double mu_earth_m3_s2{0.0};
    double earth_equatorial_radius_m{0.0};
    double mu_sun_m3_s2{astroforces::core::constants::kSunMuM3S2};
    double mu_moon_m3_s2{astroforces::core::constants::kMoonMuM3S2};
  };

  static std::unique_ptr<GravitySphAccelerationModel> Create(const Config& config);
  [[nodiscard]] GravitySphResult evaluate(const astroforces::core::StateVector& state) const;

 public:
  struct GravityFieldData;

 private:

  explicit GravitySphAccelerationModel(Config config) : config_(std::move(config)) {}

  Config config_{};
  std::shared_ptr<GravityFieldData> field_{};
  std::shared_ptr<jpl::eph::Ephemeris> ephemeris_{};
  mutable std::shared_ptr<jpl::eph::Workspace> workspace_{};
};

}  // namespace astroforces::forces

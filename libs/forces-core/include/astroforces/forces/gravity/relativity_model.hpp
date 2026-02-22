/**
 * @file relativity_model.hpp
 * @brief Post-Newtonian relativistic acceleration model (Eq. 25 decomposition).
 * @author Watosn
 */
#pragma once

#include <array>
#include <filesystem>
#include <memory>

#include "astroforces/atmo/constants.hpp"
#include "astroforces/atmo/types.hpp"

namespace jpl::eph {
class Ephemeris;
class Workspace;
}

namespace astroforces::forces {

struct RelativityResult {
  astroforces::core::Vec3 acceleration_mps2{};
  astroforces::core::Vec3 spherical_central_body_mps2{};
  astroforces::core::Vec3 geodesic_precession_mps2{};
  astroforces::core::Vec3 lense_thirring_mps2{};
  astroforces::core::Vec3 oblateness_j2_mps2{};
  astroforces::core::Vec3 rotational_energy_mps2{};
  astroforces::core::Status status{astroforces::core::Status::Ok};
};

class RelativityAccelerationModel final {
 public:
  struct Config {
    std::filesystem::path ephemeris_file{};
    bool use_spherical_central_body{true};
    bool use_geodesic_precession{true};
    bool use_lense_thirring{true};
    bool use_oblateness_j2{true};
    bool use_rotational_energy{true};

    double ppn_gamma{1.0};
    double ppn_beta{1.0};
    double lense_thirring_parameter{1.0};

    double mu_earth_m3_s2{astroforces::core::constants::kEarthMuM3S2};
    double mu_sun_m3_s2{astroforces::core::constants::kSunMuM3S2};
    double c_m_s{299792458.0};
    double earth_j2{1.08262668e-3};
    double earth_reference_radius_m{astroforces::core::constants::kEarthEquatorialRadiusM};
    double earth_angular_momentum_per_mass_m2_s{980000000.0};  // 980 km^2/s
    double earth_rotational_energy_per_mass_m2_s2{35500.0};    // 0.0355 km^2/s^2
    std::array<double, 3> earth_spin_unit{0.0, 0.0, 1.0};
  };

  static std::unique_ptr<RelativityAccelerationModel> Create(const Config& config);

  [[nodiscard]] RelativityResult evaluate(const astroforces::core::StateVector& state) const;

 private:
  explicit RelativityAccelerationModel(Config config) : config_(std::move(config)) {}

  Config config_{};
  std::shared_ptr<jpl::eph::Ephemeris> ephemeris_{};
  mutable std::shared_ptr<jpl::eph::Workspace> workspace_{};
};

}  // namespace astroforces::forces


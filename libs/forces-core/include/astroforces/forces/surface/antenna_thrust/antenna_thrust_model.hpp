/**
 * @file antenna_thrust_model.hpp
 * @brief Antenna thrust acceleration model.
 * @author Watosn
 */
#pragma once

#include <cstdint>

#include "astroforces/atmo/constants.hpp"
#include "astroforces/atmo/types.hpp"
#include "astroforces/sc/spacecraft.hpp"

namespace astroforces::forces {

enum class AntennaThrustDirectionMode : std::uint8_t { Velocity, Nadir, CustomEci, BodyFixed };

struct AntennaThrustResult {
  astroforces::core::Vec3 acceleration_mps2{};
  double thrust_n{};
  double effective_power_w{};
  double mass_kg{};
  astroforces::core::Vec3 direction_eci{};
  astroforces::core::Status status{astroforces::core::Status::Ok};
};

class AntennaThrustAccelerationModel final {
 public:
  struct Config {
    double transmit_power_w{20.0};
    double efficiency{1.0};
    double speed_of_light_mps{astroforces::core::constants::kSpeedOfLightMps};
    AntennaThrustDirectionMode direction_mode{AntennaThrustDirectionMode::Velocity};
    astroforces::core::Vec3 custom_direction_eci{1.0, 0.0, 0.0};
    astroforces::core::Vec3 body_axis{1.0, 0.0, 0.0};
  };

  AntennaThrustAccelerationModel() = default;
  explicit AntennaThrustAccelerationModel(const Config& config) : config_(config) {}

  [[nodiscard]] AntennaThrustResult evaluate(const astroforces::core::StateVector& state,
                                             const astroforces::sc::SpacecraftProperties& sc) const;

 private:
  [[nodiscard]] astroforces::core::Vec3 direction_unit_eci(const astroforces::core::StateVector& state,
                                                           astroforces::core::Status* status_out) const;

  Config config_{};
};

}  // namespace astroforces::forces

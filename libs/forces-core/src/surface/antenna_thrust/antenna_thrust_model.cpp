/**
 * @file antenna_thrust_model.cpp
 * @brief Antenna thrust acceleration model implementation.
 * @author Watosn
 */

#include "astroforces/forces/surface/antenna_thrust/antenna_thrust_model.hpp"

#include <cmath>

#include <Eigen/Dense>

namespace astroforces::forces {
namespace {

astroforces::core::Vec3 unit_direction(const astroforces::core::Vec3& v, astroforces::core::Status* status_out) {
  const double n = astroforces::core::norm(v);
  if (!(n > 0.0) || !std::isfinite(n)) {
    if (status_out) {
      *status_out = astroforces::core::Status::InvalidInput;
    }
    return astroforces::core::Vec3{};
  }
  if (status_out) {
    *status_out = astroforces::core::Status::Ok;
  }
  return v / n;
}

}  // namespace

astroforces::core::Vec3 AntennaThrustAccelerationModel::direction_unit_eci(const astroforces::core::StateVector& state,
                                                                            astroforces::core::Status* status_out) const {
  if (state.frame != astroforces::core::Frame::ECI) {
    if (status_out) {
      *status_out = astroforces::core::Status::InvalidInput;
    }
    return astroforces::core::Vec3{};
  }

  switch (static_cast<int>(config_.direction_mode)) {
    case static_cast<int>(AntennaThrustDirectionMode::Velocity):
      return unit_direction(state.velocity_mps, status_out);
    case static_cast<int>(AntennaThrustDirectionMode::Nadir):
      return unit_direction(-1.0 * state.position_m, status_out);
    case static_cast<int>(AntennaThrustDirectionMode::CustomEci):
      return unit_direction(config_.custom_direction_eci, status_out);
    case static_cast<int>(AntennaThrustDirectionMode::BodyFixed): {
      const Eigen::Vector3d axis_body(config_.body_axis.x, config_.body_axis.y, config_.body_axis.z);
      const Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> body_from_frame(state.body_from_frame_dcm.data());
      const Eigen::Vector3d axis_eci = body_from_frame.transpose() * axis_body;
      return unit_direction(astroforces::core::Vec3{axis_eci.x(), axis_eci.y(), axis_eci.z()}, status_out);
    }
    default:
      if (status_out) {
        *status_out = astroforces::core::Status::InvalidInput;
      }
      return astroforces::core::Vec3{};
  }
}

AntennaThrustResult AntennaThrustAccelerationModel::evaluate(const astroforces::core::StateVector& state,
                                                             const astroforces::sc::SpacecraftProperties& sc) const {
  if (!(sc.mass_kg > 0.0) || !(config_.speed_of_light_mps > 0.0) || !(config_.transmit_power_w >= 0.0) ||
      !(config_.efficiency >= 0.0 && config_.efficiency <= 1.0)) {
    return AntennaThrustResult{.status = astroforces::core::Status::InvalidInput};
  }

  astroforces::core::Status dir_status = astroforces::core::Status::Ok;
  const auto dir_eci = direction_unit_eci(state, &dir_status);
  if (dir_status != astroforces::core::Status::Ok) {
    return AntennaThrustResult{.status = dir_status};
  }

  const double effective_power_w = config_.transmit_power_w * config_.efficiency;
  const double thrust_n = effective_power_w / config_.speed_of_light_mps;
  const double accel_scale = thrust_n / sc.mass_kg;

  return AntennaThrustResult{
      .acceleration_mps2 = accel_scale * dir_eci,
      .thrust_n = thrust_n,
      .effective_power_w = effective_power_w,
      .mass_kg = sc.mass_kg,
      .direction_eci = dir_eci,
      .status = astroforces::core::Status::Ok,
  };
}

}  // namespace astroforces::forces

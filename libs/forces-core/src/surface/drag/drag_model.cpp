/**
 * @file drag_model.cpp
 * @brief Drag acceleration pipeline implementation.
 * @author Watosn
 */

#include "astroforces/forces/surface/drag/drag_model.hpp"

#include <cmath>
#include <optional>

#include "astroforces/core/transforms.hpp"
#include "astroforces/forces/surface/surface_force.hpp"

namespace astroforces::forces {
namespace {

astroforces::core::Vec3 body_from_frame_mul(const std::array<double, 9>& dcm, const astroforces::core::Vec3& v) {
  return astroforces::core::Vec3{
      dcm[0] * v.x + dcm[1] * v.y + dcm[2] * v.z,
      dcm[3] * v.x + dcm[4] * v.y + dcm[5] * v.z,
      dcm[6] * v.x + dcm[7] * v.y + dcm[8] * v.z,
  };
}

}  // namespace

DragResult DragAccelerationModel::evaluate(const astroforces::core::StateVector& state,
                                           const astroforces::sc::SpacecraftProperties& sc) const {
  if (sc.mass_kg <= 0.0) {
    return DragResult{.status = astroforces::core::Status::InvalidInput};
  }
  if (state.frame != astroforces::core::Frame::ECI && state.frame != astroforces::core::Frame::ECEF) {
    return DragResult{.status = astroforces::core::Status::InvalidInput};
  }

  const auto w = weather_.at(state.epoch);
  if (w.status != astroforces::core::Status::Ok) {
    return DragResult{.status = w.status};
  }
  std::optional<astroforces::core::ApproxEciEcefContext> frame_ctx{};
  const auto& get_frame_ctx = [&]() -> const astroforces::core::ApproxEciEcefContext& {
    if (!frame_ctx.has_value()) {
      frame_ctx = astroforces::core::build_approx_eci_ecef_context(state.epoch.utc_seconds);
    }
    return *frame_ctx;
  };

  astroforces::core::StateVector eval_state = state;
  if (state.frame == astroforces::core::Frame::ECI) {
    eval_state.frame = astroforces::core::Frame::ECEF;
    eval_state.position_m = astroforces::core::approx_eci_to_ecef_position(state.position_m, get_frame_ctx());
    eval_state.velocity_mps = astroforces::core::approx_eci_to_ecef_velocity(state.position_m, state.velocity_mps, get_frame_ctx());
  }

  const auto a = atmosphere_.evaluate(eval_state, w);
  if (a.status != astroforces::core::Status::Ok || a.density_kg_m3 < 0.0) {
    return DragResult{.status = a.status};
  }

  const auto wind = wind_.evaluate(eval_state, w);
  if (wind.status != astroforces::core::Status::Ok) {
    return DragResult{.status = wind.status};
  }

  astroforces::core::Vec3 wind_state_mps{};
  if (state.frame == astroforces::core::Frame::ECEF) {
    if (wind.frame == astroforces::core::Frame::ECEF) {
      wind_state_mps = wind.velocity_mps;
    } else if (wind.frame == astroforces::core::Frame::ECI) {
      const auto r_eci = astroforces::core::approx_ecef_to_eci_position(eval_state.position_m, get_frame_ctx());
      wind_state_mps = astroforces::core::approx_eci_to_ecef_velocity(r_eci, wind.velocity_mps, get_frame_ctx());
    } else {
      return DragResult{.status = astroforces::core::Status::InvalidInput};
    }
  } else {
    if (wind.frame == astroforces::core::Frame::ECI) {
      wind_state_mps = wind.velocity_mps;
    } else if (wind.frame == astroforces::core::Frame::ECEF) {
      wind_state_mps = astroforces::core::approx_ecef_to_eci_velocity(eval_state.position_m, wind.velocity_mps, get_frame_ctx());
    } else {
      return DragResult{.status = astroforces::core::Status::InvalidInput};
    }
  }

  const auto vrel = state.velocity_mps - wind_state_mps;
  double speed = 0.0;
  const auto flow_dir_frame = astroforces::forces::unit_direction(vrel, &speed);
  if (!std::isfinite(speed)) {
    return DragResult{.status = astroforces::core::Status::NumericalError};
  }

  astroforces::core::Vec3 flow_dir_body{};
  if (speed > 0.0) {
    flow_dir_body = body_from_frame_mul(state.body_from_frame_dcm, flow_dir_frame);
  }

  const double q_pa = 0.5 * a.density_kg_m3 * speed * speed;
  const auto sf = astroforces::forces::evaluate_surface_force(sc,
                                                              flow_dir_frame,
                                                              flow_dir_body,
                                                              q_pa,
                                                              sc.cd,
                                                              astroforces::sc::SurfaceCoeffModel::Drag,
                                                              -1.0);
  if (sf.status != astroforces::core::Status::Ok) {
    return DragResult{.status = sf.status};
  }

  return DragResult{.acceleration_mps2 = sf.acceleration_mps2,
                    .relative_velocity_mps = vrel,
                    .density_kg_m3 = a.density_kg_m3,
                    .temperature_k = a.temperature_k,
                    .relative_speed_mps = speed,
                    .dynamic_pressure_pa = q_pa,
                    .area_m2 = sf.area_m2,
                    .cd = sf.coeff,
                    .weather = w,
                    .status = astroforces::core::Status::Ok};
}

}  // namespace astroforces::forces

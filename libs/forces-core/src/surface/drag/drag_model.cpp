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

  std::optional<astroforces::core::ApproxEciEcefContext> approx_ctx{};
  const auto& get_approx_ctx = [&]() -> const astroforces::core::ApproxEciEcefContext& {
    if (!approx_ctx.has_value()) {
      approx_ctx = astroforces::core::build_approx_eci_ecef_context(state.epoch.utc_seconds);
    }
    return *approx_ctx;
  };

  std::optional<astroforces::core::GcrfItrfTransformContext> strict_ctx{};
  const auto get_strict_ctx = [&]() -> const astroforces::core::GcrfItrfTransformContext* {
    if (transform_mode_ != DragFrameTransformMode::StrictGcrfItrf) {
      return nullptr;
    }
    if (strict_ctx.has_value()) {
      return &(*strict_ctx);
    }
    if (!eop_series_.has_value() || !cip_series_.has_value()) {
      return nullptr;
    }
    const auto eop_sample = eop_series_->sample_at_utc_seconds(state.epoch.utc_seconds);
    const auto cip_sample = cip_series_->sample_at_utc_seconds(state.epoch.utc_seconds);
    if (!eop_sample.has_value() || !cip_sample.has_value()) {
      return nullptr;
    }
    const double jd_utc = astroforces::core::utc_seconds_to_julian_date_utc(state.epoch.utc_seconds);
    const double jd_tt = astroforces::core::utc_seconds_to_julian_date_tt(state.epoch.utc_seconds);
    strict_ctx = astroforces::core::gcrf_to_itrf_transform_context_exact(
        jd_utc,
        jd_tt,
        cip_sample->value,
        cip_sample->rate,
        eop_sample->value,
        eop_sample->rate);
    return &(*strict_ctx);
  };

  const auto eci_pos_to_ecef = [&](const astroforces::core::Vec3& r_eci_m, astroforces::core::Vec3* out) -> bool {
    if (out == nullptr) {
      return false;
    }
    if (transform_mode_ == DragFrameTransformMode::StrictGcrfItrf) {
      const auto* ctx = get_strict_ctx();
      if (ctx == nullptr) {
        return false;
      }
      *out = astroforces::core::gcrf_to_itrf_position(r_eci_m, *ctx);
      return true;
    }
    *out = astroforces::core::approx_eci_to_ecef_position(r_eci_m, get_approx_ctx());
    return true;
  };

  const auto eci_vel_to_ecef = [&](const astroforces::core::Vec3& r_eci_m,
                                   const astroforces::core::Vec3& v_eci_mps,
                                   astroforces::core::Vec3* out) -> bool {
    if (out == nullptr) {
      return false;
    }
    if (transform_mode_ == DragFrameTransformMode::StrictGcrfItrf) {
      const auto* ctx = get_strict_ctx();
      if (ctx == nullptr) {
        return false;
      }
      *out = astroforces::core::gcrf_to_itrf_velocity(r_eci_m, v_eci_mps, *ctx);
      return true;
    }
    *out = astroforces::core::approx_eci_to_ecef_velocity(r_eci_m, v_eci_mps, get_approx_ctx());
    return true;
  };

  const auto ecef_pos_to_eci = [&](const astroforces::core::Vec3& r_ecef_m, astroforces::core::Vec3* out) -> bool {
    if (out == nullptr) {
      return false;
    }
    if (transform_mode_ == DragFrameTransformMode::StrictGcrfItrf) {
      const auto* ctx = get_strict_ctx();
      if (ctx == nullptr) {
        return false;
      }
      *out = astroforces::core::itrf_to_gcrf_position(r_ecef_m, *ctx);
      return true;
    }
    *out = astroforces::core::approx_ecef_to_eci_position(r_ecef_m, get_approx_ctx());
    return true;
  };

  const auto ecef_vel_to_eci = [&](const astroforces::core::Vec3& r_ecef_m,
                                   const astroforces::core::Vec3& v_ecef_mps,
                                   astroforces::core::Vec3* out) -> bool {
    if (out == nullptr) {
      return false;
    }
    if (transform_mode_ == DragFrameTransformMode::StrictGcrfItrf) {
      const auto* ctx = get_strict_ctx();
      if (ctx == nullptr) {
        return false;
      }
      *out = astroforces::core::itrf_to_gcrf_velocity(r_ecef_m, v_ecef_mps, *ctx);
      return true;
    }
    *out = astroforces::core::approx_ecef_to_eci_velocity(r_ecef_m, v_ecef_mps, get_approx_ctx());
    return true;
  };

  astroforces::core::StateVector eval_state = state;
  if (state.frame == astroforces::core::Frame::ECI) {
    eval_state.frame = astroforces::core::Frame::ECEF;
    if (!eci_pos_to_ecef(state.position_m, &eval_state.position_m)) {
      return DragResult{.status = astroforces::core::Status::DataUnavailable};
    }
    if (!eci_vel_to_ecef(state.position_m, state.velocity_mps, &eval_state.velocity_mps)) {
      return DragResult{.status = astroforces::core::Status::DataUnavailable};
    }
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
      astroforces::core::Vec3 r_eci{};
      if (!ecef_pos_to_eci(eval_state.position_m, &r_eci)) {
        return DragResult{.status = astroforces::core::Status::DataUnavailable};
      }
      if (!eci_vel_to_ecef(r_eci, wind.velocity_mps, &wind_state_mps)) {
        return DragResult{.status = astroforces::core::Status::DataUnavailable};
      }
    } else {
      return DragResult{.status = astroforces::core::Status::InvalidInput};
    }
  } else {
    if (wind.frame == astroforces::core::Frame::ECI) {
      wind_state_mps = wind.velocity_mps;
    } else if (wind.frame == astroforces::core::Frame::ECEF) {
      if (!ecef_vel_to_eci(eval_state.position_m, wind.velocity_mps, &wind_state_mps)) {
        return DragResult{.status = astroforces::core::Status::DataUnavailable};
      }
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

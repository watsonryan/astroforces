/**
 * @file earth_radiation_model.cpp
 * @brief Earth radiation pressure acceleration model implementation.
 * @author Watosn
 */

#include "astroforces/forces/surface/earth_radiation/earth_radiation_model.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <unordered_map>

#include "astroforces/core/transforms.hpp"
#include "astroforces/forces/surface/eclipse.hpp"
#include "astroforces/forces/surface/surface_force.hpp"
#include "jpl_eph/jpl_eph.hpp"

namespace astroforces::forces {
namespace {

astroforces::core::Vec3 to_vec3(const std::array<double, 6>& pv) {
  return astroforces::core::Vec3{pv[0], pv[1], pv[2]};
}

astroforces::core::Status map_jpl_error(const jpl::eph::Status& s) {
  switch (s.code) {
    case jpl::eph::ErrorCode::kInvalidArgument:
      return astroforces::core::Status::InvalidInput;
    case jpl::eph::ErrorCode::kIo:
    case jpl::eph::ErrorCode::kCorruptFile:
    case jpl::eph::ErrorCode::kOutOfRange:
    case jpl::eph::ErrorCode::kUnsupported:
      return astroforces::core::Status::DataUnavailable;
    case jpl::eph::ErrorCode::kOk:
    default:
      return astroforces::core::Status::NumericalError;
  }
}

astroforces::core::Vec3 body_from_frame_mul(const std::array<double, 9>& dcm, const astroforces::core::Vec3& v) {
  return astroforces::core::Vec3{
      dcm[0] * v.x + dcm[1] * v.y + dcm[2] * v.z,
      dcm[3] * v.x + dcm[4] * v.y + dcm[5] * v.z,
      dcm[6] * v.x + dcm[7] * v.y + dcm[8] * v.z,
  };
}

double lambert_phase_function(double phase_angle_rad) {
  const double pi = std::numbers::pi;
  if (!std::isfinite(phase_angle_rad) || phase_angle_rad < 0.0 || phase_angle_rad > pi) {
    return 1.0;
  }
  const double phi = (std::sin(phase_angle_rad) + (pi - phase_angle_rad) * std::cos(phase_angle_rad)) / pi;
  return std::max(0.0, phi);
}

jpl::eph::Workspace& thread_local_workspace_for(const void* key) {
  thread_local std::unordered_map<const void*, jpl::eph::Workspace> workspaces{};
  return workspaces[key];
}

}  // namespace

std::unique_ptr<EarthRadiationAccelerationModel> EarthRadiationAccelerationModel::Create(const Config& config) {
  return std::unique_ptr<EarthRadiationAccelerationModel>(new EarthRadiationAccelerationModel(config));
}

EarthRadiationAccelerationModel::EarthRadiationAccelerationModel(const Config& config) : config_(config) {
  if (config_.ephemeris_file.empty()) {
    return;
  }
  auto opened = jpl::eph::Ephemeris::Open(config_.ephemeris_file.string());
  if (!opened.has_value()) {
    return;
  }
  ephemeris_ = opened.value();
}

EarthRadiationResult EarthRadiationAccelerationModel::evaluate(const astroforces::core::StateVector& state,
                                         const astroforces::sc::SpacecraftProperties& sc) const {
  if (sc.mass_kg <= 0.0 || config_.earth_reference_radius_m <= 0.0 || config_.solar_flux_w_m2 < 0.0 ||
      config_.earth_albedo < 0.0 || config_.earth_ir_flux_w_m2 < 0.0 || config_.speed_of_light_mps <= 0.0) {
    return EarthRadiationResult{.status = astroforces::core::Status::InvalidInput};
  }
  if (state.frame != astroforces::core::Frame::ECI && state.frame != astroforces::core::Frame::ECEF) {
    return EarthRadiationResult{.status = astroforces::core::Status::InvalidInput};
  }

  double earth_dist_m = 0.0;
  const auto flow_dir_frame = astroforces::forces::unit_direction(state.position_m, &earth_dist_m);
  if (!(earth_dist_m > 0.0)) {
    return EarthRadiationResult{.status = astroforces::core::Status::NumericalError};
  }

  const double ratio = config_.earth_reference_radius_m / earth_dist_m;
  const double distance_scale = ratio * ratio;

  double albedo_phase = 1.0;
  double albedo_eclipse_factor = 1.0;
  if (config_.use_albedo && ephemeris_ && state.frame == astroforces::core::Frame::ECI) {
    auto& workspace = thread_local_workspace_for(this);
    const double jd_tdb = astroforces::core::utc_seconds_to_julian_date_tdb(state.epoch.utc_seconds);
    const auto sun = ephemeris_->PlephSi(jd_tdb, jpl::eph::Body::Sun, jpl::eph::Body::Earth, false, workspace);
    if (!sun.has_value()) {
      return EarthRadiationResult{.status = map_jpl_error(sun.error())};
    }
    const auto r_sun = to_vec3(sun.value().pv);
    const double r_sun_n = astroforces::core::norm(r_sun);
    if (r_sun_n > 0.0) {
      const auto sun_hat = r_sun / r_sun_n;
      const auto sat_hat = flow_dir_frame;  // Earth -> spacecraft.
      const double cos_phase = std::clamp(astroforces::core::dot(sun_hat, sat_hat), -1.0, 1.0);
      const double phase = std::acos(cos_phase);
      albedo_phase = lambert_phase_function(phase);
    }
    if (config_.use_eclipse) {
      const auto moon = ephemeris_->PlephSi(jd_tdb, jpl::eph::Body::Moon, jpl::eph::Body::Earth, false, workspace);
      astroforces::core::Vec3 r_moon{};
      astroforces::core::Vec3* moon_ptr = nullptr;
      if (moon.has_value()) {
        r_moon = to_vec3(moon.value().pv);
        moon_ptr = &r_moon;
      }
      albedo_eclipse_factor = astroforces::forces::sun_visibility_factor(state.position_m, r_sun, moon_ptr);
    }
  }

  const double p_albedo = config_.use_albedo
                              ? (config_.earth_albedo * config_.solar_flux_w_m2 * distance_scale * albedo_phase /
                                 config_.speed_of_light_mps) *
                                    albedo_eclipse_factor
                              : 0.0;
  const double p_ir = config_.use_earth_ir
                          ? (config_.earth_ir_flux_w_m2 * distance_scale / config_.speed_of_light_mps)
                          : 0.0;
  const double pressure_pa = p_albedo + p_ir;
  if (!std::isfinite(pressure_pa) || pressure_pa < 0.0) {
    return EarthRadiationResult{.status = astroforces::core::Status::NumericalError};
  }

  const auto flow_dir_body = body_from_frame_mul(state.body_from_frame_dcm, flow_dir_frame);

  const auto sf = astroforces::forces::evaluate_surface_force(sc,
                                                              flow_dir_frame,
                                                              flow_dir_body,
                                                              pressure_pa,
                                                              sc.cr,
                                                              astroforces::sc::SurfaceCoeffModel::RadiationPressure,
                                                              +1.0);
  if (sf.status != astroforces::core::Status::Ok) {
    return EarthRadiationResult{.status = sf.status};
  }

  return EarthRadiationResult{
      .acceleration_mps2 = sf.acceleration_mps2,
      .earth_radiation_pressure_pa = pressure_pa,
      .albedo_pressure_pa = p_albedo,
      .ir_pressure_pa = p_ir,
      .albedo_phase_function = albedo_phase,
      .albedo_eclipse_factor = albedo_eclipse_factor,
      .earth_distance_m = earth_dist_m,
      .area_m2 = sf.area_m2,
      .cr = sf.coeff,
      .status = astroforces::core::Status::Ok,
  };
}

}  // namespace astroforces::forces

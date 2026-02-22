/**
 * @file srp_model.cpp
 * @brief Solar radiation pressure acceleration model implementation.
 * @author Watosn
 */

#include "astroforces/forces/surface/srp/srp_model.hpp"

#include <array>
#include <cmath>
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

jpl::eph::Workspace& thread_local_workspace_for(const void* key) {
  thread_local std::unordered_map<const void*, jpl::eph::Workspace> workspaces{};
  return workspaces[key];
}

}  // namespace

std::unique_ptr<SrpAccelerationModel> SrpAccelerationModel::Create(const Config& config) {
  auto out = std::unique_ptr<SrpAccelerationModel>(new SrpAccelerationModel(config));
  auto opened = jpl::eph::Ephemeris::Open(config.ephemeris_file.string());
  if (!opened.has_value()) {
    return out;
  }
  out->ephemeris_ = opened.value();
  return out;
}

SrpResult SrpAccelerationModel::evaluate(const astroforces::core::StateVector& state,
                                         const astroforces::sc::SpacecraftProperties& sc) const {
  if (!ephemeris_) {
    return SrpResult{.status = astroforces::core::Status::DataUnavailable};
  }
  if (state.frame != astroforces::core::Frame::ECI || sc.mass_kg <= 0.0) {
    return SrpResult{.status = astroforces::core::Status::InvalidInput};
  }

  const double jd_tdb = astroforces::core::utc_seconds_to_julian_date_tdb(state.epoch.utc_seconds);
  auto& workspace = thread_local_workspace_for(this);
  const auto sun = ephemeris_->PlephSi(jd_tdb, jpl::eph::Body::Sun, jpl::eph::Body::Earth, false, workspace);
  if (!sun.has_value()) {
    return SrpResult{.status = map_jpl_error(sun.error())};
  }
  const auto r_sun_eci_m = to_vec3(sun.value().pv);
  const auto sun_to_sc = state.position_m - r_sun_eci_m;  // Incoming photon direction.
  double sun_dist_m = 0.0;
  const auto flow_dir_frame = astroforces::forces::unit_direction(sun_to_sc, &sun_dist_m);
  if (!(sun_dist_m > 0.0)) {
    return SrpResult{.status = astroforces::core::Status::NumericalError};
  }

  astroforces::core::Vec3 r_moon_eci_m{};
  bool has_moon = false;
  if (config_.use_eclipse) {
    const auto moon = ephemeris_->PlephSi(jd_tdb, jpl::eph::Body::Moon, jpl::eph::Body::Earth, false, workspace);
    if (moon.has_value()) {
      r_moon_eci_m = to_vec3(moon.value().pv);
      has_moon = true;
    }
  }
  const double eclipse_factor =
      config_.use_eclipse ? astroforces::forces::sun_visibility_factor(state.position_m, r_sun_eci_m, has_moon ? &r_moon_eci_m : nullptr)
                          : 1.0;
  const bool eclipsed = eclipse_factor <= 0.0;
  const double inv_r2 = (config_.astronomical_unit_m / sun_dist_m) * (config_.astronomical_unit_m / sun_dist_m);
  const double pressure = config_.solar_pressure_1au_pa * inv_r2 * eclipse_factor;

  const auto flow_dir_body = body_from_frame_mul(state.body_from_frame_dcm, flow_dir_frame);

  const auto sf = astroforces::forces::evaluate_surface_force(sc,
                                                              flow_dir_frame,
                                                              flow_dir_body,
                                                              pressure,
                                                              sc.cr,
                                                              astroforces::sc::SurfaceCoeffModel::RadiationPressure,
                                                              -1.0);
  if (sf.status != astroforces::core::Status::Ok) {
    return SrpResult{.status = sf.status};
  }

  return SrpResult{
      .acceleration_mps2 = sf.acceleration_mps2,
      .solar_pressure_pa = pressure,
      .sun_distance_m = sun_dist_m,
      .eclipse_factor = eclipse_factor,
      .area_m2 = sf.area_m2,
      .cr = sf.coeff,
      .eclipsed = eclipsed,
      .status = astroforces::core::Status::Ok,
  };
}

}  // namespace astroforces::forces

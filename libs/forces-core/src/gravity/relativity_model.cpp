/**
 * @file relativity_model.cpp
 * @brief Post-Newtonian relativistic acceleration model implementation.
 * @author Watosn
 */

#include "astroforces/forces/gravity/relativity_model.hpp"

#include <array>
#include <cmath>
#include <unordered_map>

#include "astroforces/core/transforms.hpp"
#include "jpl_eph/jpl_eph.hpp"

namespace astroforces::forces {
namespace {

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

astroforces::core::Vec3 to_vec3(const std::array<double, 6>& pv) { return astroforces::core::Vec3{pv[0], pv[1], pv[2]}; }

bool finite_vec(const astroforces::core::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

jpl::eph::Workspace& thread_local_workspace_for(const void* key) {
  thread_local std::unordered_map<const void*, jpl::eph::Workspace> workspaces{};
  return workspaces[key];
}

}  // namespace

std::unique_ptr<RelativityAccelerationModel> RelativityAccelerationModel::Create(const Config& config) {
  auto out = std::unique_ptr<RelativityAccelerationModel>(new RelativityAccelerationModel(config));
  if (!config.use_geodesic_precession) {
    return out;
  }
  auto opened = jpl::eph::Ephemeris::Open(config.ephemeris_file.string());
  if (!opened.has_value()) {
    return out;
  }
  out->ephemeris_ = opened.value();
  return out;
}

RelativityResult RelativityAccelerationModel::evaluate(const astroforces::core::StateVector& state) const {
  if (state.frame != astroforces::core::Frame::ECI) {
    return RelativityResult{.status = astroforces::core::Status::InvalidInput};
  }
  if (!(config_.c_m_s > 0.0) || !(config_.mu_earth_m3_s2 > 0.0) || !(config_.earth_reference_radius_m > 0.0)) {
    return RelativityResult{.status = astroforces::core::Status::InvalidInput};
  }

  const auto& r = state.position_m;
  const auto& v = state.velocity_mps;
  const double r2 = astroforces::core::dot(r, r);
  const double rnorm = std::sqrt(r2);
  const double inv_r = 1.0 / rnorm;
  const double inv_r2 = inv_r * inv_r;
  const double inv_r3 = inv_r2 * inv_r;
  const double v2 = astroforces::core::dot(v, v);
  const double rv = astroforces::core::dot(r, v);
  if (!(rnorm > 0.0) || !std::isfinite(rnorm) || !std::isfinite(v2) || !std::isfinite(rv)) {
    return RelativityResult{.status = astroforces::core::Status::NumericalError};
  }

  const double c2 = config_.c_m_s * config_.c_m_s;
  RelativityResult out{};

  if (config_.use_spherical_central_body) {
    const double pref = config_.mu_earth_m3_s2 * inv_r3 / c2;
    const double scale_r = 2.0 * (config_.ppn_beta + config_.ppn_gamma) * config_.mu_earth_m3_s2 * inv_r - config_.ppn_gamma * v2;
    const double scale_v = 2.0 * (1.0 + config_.ppn_gamma) * rv;
    out.spherical_central_body_mps2 = pref * (scale_r * r + scale_v * v);
  }

  if (config_.use_geodesic_precession) {
    if (!ephemeris_) {
      return RelativityResult{.status = astroforces::core::Status::DataUnavailable};
    }
    auto& workspace = thread_local_workspace_for(this);
    const double jd_tdb = astroforces::core::utc_seconds_to_julian_date_tdb(state.epoch.utc_seconds);
    const auto sun_wrt_earth = ephemeris_->PlephSi(jd_tdb, jpl::eph::Body::Sun, jpl::eph::Body::Earth, true, workspace);
    if (!sun_wrt_earth.has_value()) {
      return RelativityResult{.status = map_jpl_error(sun_wrt_earth.error())};
    }
    const auto R_es = -1.0 * to_vec3(sun_wrt_earth.value().pv);
    const auto V_es =
        astroforces::core::Vec3{-sun_wrt_earth.value().pv[3], -sun_wrt_earth.value().pv[4], -sun_wrt_earth.value().pv[5]};
    const double R = astroforces::core::norm(R_es);
    if (!(R > 0.0)) {
      return RelativityResult{.status = astroforces::core::Status::NumericalError};
    }
    const double inv_r_sun = 1.0 / R;
    const double pref = -(1.0 + 2.0 * config_.ppn_gamma) * config_.mu_sun_m3_s2 * inv_r_sun * inv_r_sun * inv_r_sun / c2;
    out.geodesic_precession_mps2 = pref * astroforces::core::cross(astroforces::core::cross(V_es, R_es), v);
  }

  if (config_.use_lense_thirring) {
    const astroforces::core::Vec3 omega_hat{
        config_.earth_spin_unit[0], config_.earth_spin_unit[1], config_.earth_spin_unit[2]};
    const astroforces::core::Vec3 J = config_.earth_angular_momentum_per_mass_m2_s * omega_hat;
    const double pref = (1.0 + config_.ppn_gamma) * config_.lense_thirring_parameter * config_.mu_earth_m3_s2 * inv_r3 / c2;
    out.lense_thirring_mps2 = pref * ((3.0 * astroforces::core::dot(r, J) * inv_r2) * astroforces::core::cross(r, v)
                                      + astroforces::core::cross(v, J));
  }

  if (config_.use_oblateness_j2) {
    const double x = r.x;
    const double y = r.y;
    const double z = r.z;
    const double vz = v.z;
    const double z2_over_r2 = (z * z) * inv_r2;
    const double c = config_.mu_earth_m3_s2 * config_.earth_j2 * config_.earth_reference_radius_m * config_.earth_reference_radius_m *
                     inv_r2 * inv_r2 / c2;
    const double term_a = 2.5 * v2 - 4.0 * config_.mu_earth_m3_s2 * inv_r;
    const double term_b = 6.0 * rv * inv_r * vz;
    const double mu_over_r = config_.mu_earth_m3_s2 * inv_r;
    const double ax = c * x * inv_r *
                      (term_a * (1.0 - 5.0 * z2_over_r2) - term_b * (1.0 - 5.0 * z2_over_r2) + 2.0 * mu_over_r * z2_over_r2);
    const double ay = c * y * inv_r *
                      (term_a * (1.0 - 5.0 * z2_over_r2) - term_b * (1.0 - 5.0 * z2_over_r2) + 2.0 * mu_over_r * z2_over_r2);
    const double az = c * z * inv_r *
                      (term_a * (3.0 - 5.0 * z2_over_r2) + term_b * (1.0 - 5.0 * z2_over_r2) + 2.0 * mu_over_r * (3.0 - z2_over_r2));
    out.oblateness_j2_mps2 = astroforces::core::Vec3{ax, ay, az};
  }

  if (config_.use_rotational_energy) {
    const astroforces::core::Vec3 omega_hat{
        config_.earth_spin_unit[0], config_.earth_spin_unit[1], config_.earth_spin_unit[2]};
    const double z2_over_r2 = (r.z * r.z) * inv_r2;
    const double pref = -(1.0 + config_.ppn_gamma) * config_.earth_rotational_energy_per_mass_m2_s2 * config_.mu_earth_m3_s2 *
                        config_.earth_reference_radius_m * config_.earth_reference_radius_m * inv_r2 * inv_r2 / (4.0 * c2);
    const auto r_hat = r * inv_r;
    const astroforces::core::Vec3 term = (1.0 - 5.0 * z2_over_r2) * r_hat + 2.0 * astroforces::core::dot(r_hat, omega_hat) * omega_hat;
    out.rotational_energy_mps2 = pref * term;
  }

  out.acceleration_mps2 = out.spherical_central_body_mps2 + out.geodesic_precession_mps2 + out.lense_thirring_mps2 +
                          out.oblateness_j2_mps2 + out.rotational_energy_mps2;
  if (!finite_vec(out.acceleration_mps2)) {
    return RelativityResult{.status = astroforces::core::Status::NumericalError};
  }
  out.status = astroforces::core::Status::Ok;
  return out;
}

}  // namespace astroforces::forces

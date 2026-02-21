/**
 * @file relativity_model.cpp
 * @brief Post-Newtonian relativistic acceleration model implementation.
 * @author Watosn
 */

#include "astroforces/forces/relativity_model.hpp"

#include <array>
#include <cmath>

#include "astroforces/atmo/conversions.hpp"
#include "jpl_eph/jpl_eph.hpp"

namespace astroforces::forces {
namespace {

astroforces::atmo::Status map_jpl_error(const jpl::eph::Status& s) {
  switch (s.code) {
    case jpl::eph::ErrorCode::kInvalidArgument:
      return astroforces::atmo::Status::InvalidInput;
    case jpl::eph::ErrorCode::kIo:
    case jpl::eph::ErrorCode::kCorruptFile:
    case jpl::eph::ErrorCode::kOutOfRange:
    case jpl::eph::ErrorCode::kUnsupported:
      return astroforces::atmo::Status::DataUnavailable;
    case jpl::eph::ErrorCode::kOk:
    default:
      return astroforces::atmo::Status::NumericalError;
  }
}

astroforces::atmo::Vec3 to_vec3(const std::array<double, 6>& pv) { return astroforces::atmo::Vec3{pv[0], pv[1], pv[2]}; }

astroforces::atmo::Vec3 cross(const astroforces::atmo::Vec3& a, const astroforces::atmo::Vec3& b) {
  return astroforces::atmo::Vec3{
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x,
  };
}

bool finite_vec(const astroforces::atmo::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
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
  out->workspace_ = std::make_shared<jpl::eph::Workspace>();
  return out;
}

RelativityResult RelativityAccelerationModel::evaluate(const astroforces::atmo::StateVector& state) const {
  if (state.frame != astroforces::atmo::Frame::ECI) {
    return RelativityResult{.status = astroforces::atmo::Status::InvalidInput};
  }
  if (!(config_.c_m_s > 0.0) || !(config_.mu_earth_m3_s2 > 0.0) || !(config_.earth_reference_radius_m > 0.0)) {
    return RelativityResult{.status = astroforces::atmo::Status::InvalidInput};
  }

  const auto& r = state.position_m;
  const auto& v = state.velocity_mps;
  const double r2 = astroforces::atmo::dot(r, r);
  const double rnorm = std::sqrt(r2);
  const double v2 = astroforces::atmo::dot(v, v);
  const double rv = astroforces::atmo::dot(r, v);
  if (!(rnorm > 0.0) || !std::isfinite(rnorm) || !std::isfinite(v2) || !std::isfinite(rv)) {
    return RelativityResult{.status = astroforces::atmo::Status::NumericalError};
  }

  const double c2 = config_.c_m_s * config_.c_m_s;
  RelativityResult out{};

  if (config_.use_spherical_central_body) {
    const double pref = config_.mu_earth_m3_s2 / (c2 * rnorm * rnorm * rnorm);
    const double scale_r = 2.0 * (config_.ppn_beta + config_.ppn_gamma) * config_.mu_earth_m3_s2 / rnorm - config_.ppn_gamma * v2;
    const double scale_v = 2.0 * (1.0 + config_.ppn_gamma) * rv;
    out.spherical_central_body_mps2 = pref * (scale_r * r + scale_v * v);
  }

  if (config_.use_geodesic_precession) {
    if (!ephemeris_ || !workspace_) {
      return RelativityResult{.status = astroforces::atmo::Status::DataUnavailable};
    }
    const double jd_utc = astroforces::atmo::utc_seconds_to_julian_date_utc(state.epoch.utc_seconds);
    const auto sun_wrt_earth = ephemeris_->PlephSi(jd_utc, jpl::eph::Body::Sun, jpl::eph::Body::Earth, true, *workspace_);
    if (!sun_wrt_earth.has_value()) {
      return RelativityResult{.status = map_jpl_error(sun_wrt_earth.error())};
    }
    const auto R_es = -1.0 * to_vec3(sun_wrt_earth.value().pv);
    const auto V_es =
        astroforces::atmo::Vec3{-sun_wrt_earth.value().pv[3], -sun_wrt_earth.value().pv[4], -sun_wrt_earth.value().pv[5]};
    const double R = astroforces::atmo::norm(R_es);
    if (!(R > 0.0)) {
      return RelativityResult{.status = astroforces::atmo::Status::NumericalError};
    }
    const double pref = -(1.0 + 2.0 * config_.ppn_gamma) * config_.mu_sun_m3_s2 / (c2 * R * R * R);
    out.geodesic_precession_mps2 = pref * cross(cross(V_es, R_es), v);
  }

  if (config_.use_lense_thirring) {
    const astroforces::atmo::Vec3 omega_hat{
        config_.earth_spin_unit[0], config_.earth_spin_unit[1], config_.earth_spin_unit[2]};
    const astroforces::atmo::Vec3 J = config_.earth_angular_momentum_per_mass_m2_s * omega_hat;
    const double pref = (1.0 + config_.ppn_gamma) * config_.lense_thirring_parameter * config_.mu_earth_m3_s2 /
                        (c2 * rnorm * rnorm * rnorm);
    out.lense_thirring_mps2 =
        pref * ((3.0 * astroforces::atmo::dot(r, J) / r2) * cross(r, v) + cross(v, J));
  }

  if (config_.use_oblateness_j2) {
    const double x = r.x;
    const double y = r.y;
    const double z = r.z;
    const double vz = v.z;
    const double r4 = r2 * r2;
    const double z2_over_r2 = (z * z) / r2;
    const double c = config_.mu_earth_m3_s2 * config_.earth_j2 * config_.earth_reference_radius_m * config_.earth_reference_radius_m /
                     (c2 * r4);
    const double term_a = 2.5 * v2 - 4.0 * config_.mu_earth_m3_s2 / rnorm;
    const double term_b = 6.0 * rv / rnorm * vz;
    const double mu_over_r = config_.mu_earth_m3_s2 / rnorm;
    const double ax = c * x / rnorm *
                      (term_a * (1.0 - 5.0 * z2_over_r2) - term_b * (1.0 - 5.0 * z2_over_r2) + 2.0 * mu_over_r * z2_over_r2);
    const double ay = c * y / rnorm *
                      (term_a * (1.0 - 5.0 * z2_over_r2) - term_b * (1.0 - 5.0 * z2_over_r2) + 2.0 * mu_over_r * z2_over_r2);
    const double az = c * z / rnorm *
                      (term_a * (3.0 - 5.0 * z2_over_r2) + term_b * (1.0 - 5.0 * z2_over_r2) + 2.0 * mu_over_r * (3.0 - z2_over_r2));
    out.oblateness_j2_mps2 = astroforces::atmo::Vec3{ax, ay, az};
  }

  if (config_.use_rotational_energy) {
    const astroforces::atmo::Vec3 omega_hat{
        config_.earth_spin_unit[0], config_.earth_spin_unit[1], config_.earth_spin_unit[2]};
    const double z2_over_r2 = (r.z * r.z) / r2;
    const double pref = -(1.0 + config_.ppn_gamma) * config_.earth_rotational_energy_per_mass_m2_s2 * config_.mu_earth_m3_s2 *
                        config_.earth_reference_radius_m * config_.earth_reference_radius_m /
                        (4.0 * c2 * r2 * r2);
    const astroforces::atmo::Vec3 term = (1.0 - 5.0 * z2_over_r2) * (r / rnorm) + 2.0 * astroforces::atmo::dot(r / rnorm, omega_hat) * omega_hat;
    out.rotational_energy_mps2 = pref * term;
  }

  out.acceleration_mps2 = out.spherical_central_body_mps2 + out.geodesic_precession_mps2 + out.lense_thirring_mps2 +
                          out.oblateness_j2_mps2 + out.rotational_energy_mps2;
  if (!finite_vec(out.acceleration_mps2)) {
    return RelativityResult{.status = astroforces::atmo::Status::NumericalError};
  }
  out.status = astroforces::atmo::Status::Ok;
  return out;
}

}  // namespace astroforces::forces

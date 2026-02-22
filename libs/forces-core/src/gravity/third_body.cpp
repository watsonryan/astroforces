/**
 * @file third_body.cpp
 * @brief Third-body perturbation model implementation.
 * @author Watosn
 */

#include "astroforces/forces/gravity/third_body.hpp"

#include <cmath>

#include "astroforces/core/transforms.hpp"
#include "astroforces/core/constants.hpp"
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

astroforces::core::Vec3 third_body_direct_indirect(const astroforces::core::Vec3& r_sc_m,
                                               const astroforces::core::Vec3& r_tb_m,
                                               double mu_m3_s2) {
  const astroforces::core::Vec3 rho = r_tb_m - r_sc_m;
  const double rho2 = astroforces::core::dot(rho, rho);
  const double rtb2 = astroforces::core::dot(r_tb_m, r_tb_m);
  if (!(rho2 > 0.0) || !(rtb2 > 0.0)) {
    return astroforces::core::Vec3{};
  }
  const double inv_rho = 1.0 / std::sqrt(rho2);
  const double inv_rtb = 1.0 / std::sqrt(rtb2);
  const double inv_rho3 = inv_rho * inv_rho * inv_rho;
  const double inv_rtb3 = inv_rtb * inv_rtb * inv_rtb;
  return (mu_m3_s2 * inv_rho3) * rho - (mu_m3_s2 * inv_rtb3) * r_tb_m;
}

astroforces::core::Vec3 goce_eq79_indirect_j2(const astroforces::core::Vec3& r_tb_m, double mu_m3_s2) {
  const double r2 = astroforces::core::dot(r_tb_m, r_tb_m);
  if (!(r2 > 0.0)) {
    return astroforces::core::Vec3{};
  }
  const double inv_r = 1.0 / std::sqrt(r2);
  const double inv_r5 = inv_r * inv_r * inv_r * inv_r * inv_r;
  const double z2_over_r2 = (r_tb_m.z * r_tb_m.z) / r2;
  const double common_xy = 1.0 - 5.0 * z2_over_r2;
  const double common_z = 3.0 - 5.0 * z2_over_r2;
  const double ae2 = astroforces::core::constants::kEarthEquatorialRadiusM * astroforces::core::constants::kEarthEquatorialRadiusM;
  const double coeff = -mu_m3_s2 * astroforces::core::constants::kEarthC20FullyNormalized * ae2 * inv_r5;
  return astroforces::core::Vec3{
      coeff * r_tb_m.x * common_xy,
      coeff * r_tb_m.y * common_xy,
      coeff * r_tb_m.z * common_z,
  };
}

jpl::eph::Workspace& thread_local_workspace() {
  thread_local jpl::eph::Workspace workspace{};
  return workspace;
}

}  // namespace

ThirdBodyPerturbationModel::ThirdBodyPerturbationModel(Config config) : config_(std::move(config)) {}

std::unique_ptr<ThirdBodyPerturbationModel> ThirdBodyPerturbationModel::Create(const Config& config) {
  auto out = std::unique_ptr<ThirdBodyPerturbationModel>(new ThirdBodyPerturbationModel(config));
  auto opened = jpl::eph::Ephemeris::Open(config.ephemeris_file.string());
  if (!opened.has_value()) {
    return out;
  }
  out->ephemeris_ = opened.value();
  return out;
}

PerturbationContribution ThirdBodyPerturbationModel::evaluate(const PerturbationRequest& request) const {
  PerturbationContribution out{};
  out.name = config_.name;
  out.type = PerturbationType::ThirdBody;

  if (!ephemeris_) {
    out.status = astroforces::core::Status::DataUnavailable;
    return out;
  }
  auto& workspace = thread_local_workspace();
  if (request.state.frame != astroforces::core::Frame::ECI) {
    out.status = astroforces::core::Status::InvalidInput;
    return out;
  }

  const double jed_tdb = astroforces::core::utc_seconds_to_julian_date_tdb(request.state.epoch.utc_seconds);

  if (config_.use_sun) {
    const auto sun = ephemeris_->PlephSi(jed_tdb, jpl::eph::Body::Sun, jpl::eph::Body::Earth, false, workspace);
    if (!sun.has_value()) {
      out.status = map_jpl_error(sun.error());
      return out;
    }
    const auto r_sun = to_vec3(sun.value().pv);
    out.acceleration_mps2 = out.acceleration_mps2 + third_body_direct_indirect(request.state.position_m, r_sun, config_.mu_sun_m3_s2);
    if (config_.use_goce_eq79_indirect_j2) {
      // GOCE Standards GO-TN-HPF-GS-0111 Eq. (79): indirect J2 term for Sun/Moon.
      out.acceleration_mps2 = out.acceleration_mps2 + goce_eq79_indirect_j2(r_sun, config_.mu_sun_m3_s2);
    }
  }

  if (config_.use_moon) {
    const auto moon = ephemeris_->PlephSi(jed_tdb, jpl::eph::Body::Moon, jpl::eph::Body::Earth, false, workspace);
    if (!moon.has_value()) {
      out.status = map_jpl_error(moon.error());
      return out;
    }
    const auto r_moon = to_vec3(moon.value().pv);
    out.acceleration_mps2 = out.acceleration_mps2 + third_body_direct_indirect(request.state.position_m, r_moon, config_.mu_moon_m3_s2);
    if (config_.use_goce_eq79_indirect_j2) {
      // GOCE Standards GO-TN-HPF-GS-0111 Eq. (79): indirect J2 term for Sun/Moon.
      out.acceleration_mps2 = out.acceleration_mps2 + goce_eq79_indirect_j2(r_moon, config_.mu_moon_m3_s2);
    }
  }

  out.status = astroforces::core::Status::Ok;
  return out;
}

}  // namespace astroforces::forces

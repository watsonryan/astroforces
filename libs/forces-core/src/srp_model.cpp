/**
 * @file srp_model.cpp
 * @brief Solar radiation pressure acceleration model implementation.
 * @author Watosn
 */

#include "astroforces/srp/srp_model.hpp"

#include <array>
#include <cmath>

#include <Eigen/Dense>

#include "astroforces/atmo/conversions.hpp"
#include "astroforces/forces/surface_force.hpp"
#include "jpl_eph/jpl_eph.hpp"

namespace astroforces::srp {
namespace {

astroforces::atmo::Vec3 to_vec3(const std::array<double, 6>& pv) {
  return astroforces::atmo::Vec3{pv[0], pv[1], pv[2]};
}

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

bool in_cylindrical_umbra(const astroforces::atmo::Vec3& r_sc_eci_m,
                          const astroforces::atmo::Vec3& r_sun_eci_m) {
  const double r_sun = astroforces::atmo::norm(r_sun_eci_m);
  if (!(r_sun > 0.0)) {
    return false;
  }
  const auto s_hat = r_sun_eci_m / r_sun;  // Earth->Sun direction.
  const double proj = astroforces::atmo::dot(r_sc_eci_m, s_hat);
  if (proj >= 0.0) {
    return false;
  }
  const auto perp = r_sc_eci_m - proj * s_hat;
  return astroforces::atmo::norm(perp) < astroforces::atmo::constants::kEarthRadiusWgs84M;
}

}  // namespace

std::unique_ptr<SrpAccelerationModel> SrpAccelerationModel::Create(const Config& config) {
  auto out = std::unique_ptr<SrpAccelerationModel>(new SrpAccelerationModel(config));
  auto opened = jpl::eph::Ephemeris::Open(config.ephemeris_file.string());
  if (!opened.has_value()) {
    return out;
  }
  out->ephemeris_ = opened.value();
  out->workspace_ = std::make_shared<jpl::eph::Workspace>();
  return out;
}

SrpResult SrpAccelerationModel::evaluate(const astroforces::atmo::StateVector& state,
                                         const astroforces::sc::SpacecraftProperties& sc) const {
  if (!ephemeris_ || !workspace_) {
    return SrpResult{.status = astroforces::atmo::Status::DataUnavailable};
  }
  if (state.frame != astroforces::atmo::Frame::ECI || sc.mass_kg <= 0.0) {
    return SrpResult{.status = astroforces::atmo::Status::InvalidInput};
  }

  const double jd_utc = astroforces::atmo::utc_seconds_to_julian_date_utc(state.epoch.utc_seconds);
  const auto sun = ephemeris_->PlephSi(jd_utc, jpl::eph::Body::Sun, jpl::eph::Body::Earth, false, *workspace_);
  if (!sun.has_value()) {
    return SrpResult{.status = map_jpl_error(sun.error())};
  }
  const auto r_sun_eci_m = to_vec3(sun.value().pv);
  const auto sun_to_sc = state.position_m - r_sun_eci_m;  // Incoming photon direction.
  double sun_dist_m = 0.0;
  const auto flow_dir_frame = astroforces::forces::unit_direction(sun_to_sc, &sun_dist_m);
  if (!(sun_dist_m > 0.0)) {
    return SrpResult{.status = astroforces::atmo::Status::NumericalError};
  }

  const bool eclipsed = config_.use_eclipse && in_cylindrical_umbra(state.position_m, r_sun_eci_m);
  const double inv_r2 = (config_.astronomical_unit_m / sun_dist_m) * (config_.astronomical_unit_m / sun_dist_m);
  const double pressure = eclipsed ? 0.0 : (config_.solar_pressure_1au_pa * inv_r2);

  astroforces::atmo::Vec3 flow_dir_body{};
  {
    const Eigen::Vector3d flow_frame(flow_dir_frame.x, flow_dir_frame.y, flow_dir_frame.z);
    const Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> body_from_frame(state.body_from_frame_dcm.data());
    const Eigen::Vector3d flow_body = body_from_frame * flow_frame;
    flow_dir_body = astroforces::atmo::Vec3{flow_body.x(), flow_body.y(), flow_body.z()};
  }

  const auto sf = astroforces::forces::evaluate_surface_force(sc,
                                                              flow_dir_frame,
                                                              flow_dir_body,
                                                              pressure,
                                                              sc.cr,
                                                              astroforces::sc::SurfaceCoeffModel::RadiationPressure,
                                                              -1.0);
  if (sf.status != astroforces::atmo::Status::Ok) {
    return SrpResult{.status = sf.status};
  }

  return SrpResult{
      .acceleration_mps2 = sf.acceleration_mps2,
      .solar_pressure_pa = pressure,
      .sun_distance_m = sun_dist_m,
      .area_m2 = sf.area_m2,
      .cr = sf.coeff,
      .eclipsed = eclipsed,
      .status = astroforces::atmo::Status::Ok,
  };
}

}  // namespace astroforces::srp


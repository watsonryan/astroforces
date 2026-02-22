/**
 * @file erp_model.cpp
 * @brief Earth radiation pressure acceleration model implementation.
 * @author Watosn
 */

#include "astroforces/forces/surface/erp/erp_model.hpp"

#include <cmath>

#include <Eigen/Dense>

#include "astroforces/forces/surface/surface_force.hpp"

namespace astroforces::forces {

ErpResult ErpAccelerationModel::evaluate(const astroforces::core::StateVector& state,
                                         const astroforces::sc::SpacecraftProperties& sc) const {
  if (sc.mass_kg <= 0.0 || config_.earth_reference_radius_m <= 0.0 || config_.earth_radiation_pressure_ref_pa < 0.0) {
    return ErpResult{.status = astroforces::core::Status::InvalidInput};
  }
  if (state.frame != astroforces::core::Frame::ECI && state.frame != astroforces::core::Frame::ECEF) {
    return ErpResult{.status = astroforces::core::Status::InvalidInput};
  }

  double earth_dist_m = 0.0;
  const auto flow_dir_frame = astroforces::forces::unit_direction(state.position_m, &earth_dist_m);
  if (!(earth_dist_m > 0.0)) {
    return ErpResult{.status = astroforces::core::Status::NumericalError};
  }

  const double ratio = config_.earth_reference_radius_m / earth_dist_m;
  const double pressure_pa = config_.earth_radiation_pressure_ref_pa * ratio * ratio;
  if (!std::isfinite(pressure_pa)) {
    return ErpResult{.status = astroforces::core::Status::NumericalError};
  }

  astroforces::core::Vec3 flow_dir_body{};
  {
    const Eigen::Vector3d flow_frame(flow_dir_frame.x, flow_dir_frame.y, flow_dir_frame.z);
    const Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> body_from_frame(state.body_from_frame_dcm.data());
    const Eigen::Vector3d flow_body = body_from_frame * flow_frame;
    flow_dir_body = astroforces::core::Vec3{flow_body.x(), flow_body.y(), flow_body.z()};
  }

  const auto sf = astroforces::forces::evaluate_surface_force(sc,
                                                              flow_dir_frame,
                                                              flow_dir_body,
                                                              pressure_pa,
                                                              sc.cr,
                                                              astroforces::sc::SurfaceCoeffModel::RadiationPressure,
                                                              +1.0);
  if (sf.status != astroforces::core::Status::Ok) {
    return ErpResult{.status = sf.status};
  }

  return ErpResult{
      .acceleration_mps2 = sf.acceleration_mps2,
      .earth_radiation_pressure_pa = pressure_pa,
      .earth_distance_m = earth_dist_m,
      .area_m2 = sf.area_m2,
      .cr = sf.coeff,
      .status = astroforces::core::Status::Ok,
  };
}

}  // namespace astroforces::forces


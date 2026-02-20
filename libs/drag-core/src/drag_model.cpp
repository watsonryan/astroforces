/**
 * @file drag_model.cpp
 * @brief Drag acceleration pipeline implementation.
 * @author Watosn
 */

#include "dragcpp/drag/drag_model.hpp"

#include <cmath>

#include <Eigen/Dense>

namespace dragcpp::drag {

DragResult DragAccelerationModel::evaluate(const dragcpp::atmo::StateVector& state,
                                           const dragcpp::sc::SpacecraftProperties& sc) const {
  if (sc.mass_kg <= 0.0) {
    return DragResult{.status = dragcpp::atmo::Status::InvalidInput};
  }

  const auto w = weather_.at(state.epoch);
  if (w.status != dragcpp::atmo::Status::Ok) {
    return DragResult{.status = w.status};
  }

  const auto a = atmosphere_.evaluate(state, w);
  if (a.status != dragcpp::atmo::Status::Ok || a.density_kg_m3 < 0.0) {
    return DragResult{.status = a.status};
  }

  const auto wind = wind_.evaluate(state, w);
  if (wind.status != dragcpp::atmo::Status::Ok || wind.frame != state.frame) {
    return DragResult{.status = dragcpp::atmo::Status::InvalidInput};
  }

  const auto vrel = state.velocity_mps - wind.velocity_mps;
  const double speed = dragcpp::atmo::norm(vrel);
  if (!std::isfinite(speed)) {
    return DragResult{.status = dragcpp::atmo::Status::NumericalError};
  }

  dragcpp::atmo::Vec3 flow_dir_body{};
  if (speed > 0.0) {
    const Eigen::Vector3d flow_frame(vrel.x / speed, vrel.y / speed, vrel.z / speed);
    const Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor>> body_from_frame(state.body_from_frame_dcm.data());
    const Eigen::Vector3d flow_body = body_from_frame * flow_frame;
    flow_dir_body = dragcpp::atmo::Vec3{flow_body.x(), flow_body.y(), flow_body.z()};
  }
  const auto aero = dragcpp::sc::projected_area_and_cd(sc, flow_dir_body);
  const double area = aero.area_m2;
  const double cd = aero.cd_effective;

  const double coeff = -0.5 * a.density_kg_m3 * cd * area / sc.mass_kg;
  const auto accel = coeff * speed * vrel;
  const double q_pa = 0.5 * a.density_kg_m3 * speed * speed;

  return DragResult{.acceleration_mps2 = accel,
                    .relative_velocity_mps = vrel,
                    .density_kg_m3 = a.density_kg_m3,
                    .temperature_k = a.temperature_k,
                    .relative_speed_mps = speed,
                    .dynamic_pressure_pa = q_pa,
                    .area_m2 = area,
                    .cd = cd,
                    .weather = w,
                    .status = dragcpp::atmo::Status::Ok};
}

}  // namespace dragcpp::drag

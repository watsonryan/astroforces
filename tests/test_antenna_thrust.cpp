/**
 * @file test_antenna_thrust.cpp
 * @brief Antenna thrust model unit/integration checks.
 * @author Watosn
 */

#include <cmath>

#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/antenna_thrust/antenna_thrust_model.hpp"
#include "astroforces/forces/surface/antenna_thrust/antenna_thrust_perturbation.hpp"
#include "astroforces/sc/spacecraft.hpp"

namespace {

bool approx(double a, double b, double rel = 1e-12) {
  const double d = std::abs(a - b);
  const double n = std::max(std::abs(b), 1e-30);
  return d / n <= rel;
}

bool finite_vec(const astroforces::core::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

}  // namespace

int main() {
  astroforces::core::StateVector state{};
  state.frame = astroforces::core::Frame::ECI;
  state.epoch.utc_seconds = 1.0e9;
  state.position_m = astroforces::core::Vec3{7000e3, 0.0, 0.0};
  state.velocity_mps = astroforces::core::Vec3{0.0, 7500.0, 0.0};

  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = 600.0, .reference_area_m2 = 4.0, .cd = 2.2, .cr = 1.3, .use_surface_model = false, .surfaces = {}};

  const astroforces::forces::AntennaThrustAccelerationModel zero_power({
      .transmit_power_w = 0.0,
      .efficiency = 1.0,
      .direction_mode = astroforces::forces::AntennaThrustDirectionMode::Velocity,
  });
  const auto out0 = zero_power.evaluate(state, sc);
  if (out0.status != astroforces::core::Status::Ok) {
    spdlog::error("zero-power evaluation failed");
    return 1;
  }
  if (!(astroforces::core::norm(out0.acceleration_mps2) == 0.0)) {
    spdlog::error("zero-power acceleration should be zero");
    return 2;
  }

  const astroforces::forces::AntennaThrustAccelerationModel p20({
      .transmit_power_w = 20.0,
      .efficiency = 1.0,
      .direction_mode = astroforces::forces::AntennaThrustDirectionMode::Velocity,
  });
  const astroforces::forces::AntennaThrustAccelerationModel p40({
      .transmit_power_w = 40.0,
      .efficiency = 1.0,
      .direction_mode = astroforces::forces::AntennaThrustDirectionMode::Velocity,
  });
  const auto r20 = p20.evaluate(state, sc);
  const auto r40 = p40.evaluate(state, sc);
  if (r20.status != astroforces::core::Status::Ok || r40.status != astroforces::core::Status::Ok) {
    spdlog::error("power-scaling evaluation failed");
    return 3;
  }
  if (!approx(astroforces::core::norm(r40.acceleration_mps2), 2.0 * astroforces::core::norm(r20.acceleration_mps2), 1e-12)) {
    spdlog::error("power scaling mismatch");
    return 4;
  }
  if (!approx(r20.direction_eci.x, 0.0) || !approx(r20.direction_eci.y, 1.0) || !approx(r20.direction_eci.z, 0.0)) {
    spdlog::error("velocity direction mismatch");
    return 5;
  }

  const astroforces::forces::AntennaThrustAccelerationModel nadir({
      .transmit_power_w = 20.0,
      .efficiency = 1.0,
      .direction_mode = astroforces::forces::AntennaThrustDirectionMode::Nadir,
  });
  const auto rn = nadir.evaluate(state, sc);
  if (!approx(rn.direction_eci.x, -1.0) || !approx(rn.direction_eci.y, 0.0) || !approx(rn.direction_eci.z, 0.0)) {
    spdlog::error("nadir direction mismatch");
    return 6;
  }

  const astroforces::forces::AntennaThrustAccelerationModel custom({
      .transmit_power_w = 20.0,
      .efficiency = 1.0,
      .direction_mode = astroforces::forces::AntennaThrustDirectionMode::CustomEci,
      .custom_direction_eci = astroforces::core::Vec3{2.0, 0.0, 0.0},
  });
  const auto rc = custom.evaluate(state, sc);
  if (!approx(rc.direction_eci.x, 1.0) || !approx(rc.direction_eci.y, 0.0) || !approx(rc.direction_eci.z, 0.0)) {
    spdlog::error("custom direction mismatch");
    return 7;
  }

  state.body_from_frame_dcm = {0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
  const astroforces::forces::AntennaThrustAccelerationModel body({
      .transmit_power_w = 20.0,
      .efficiency = 1.0,
      .direction_mode = astroforces::forces::AntennaThrustDirectionMode::BodyFixed,
      .body_axis = astroforces::core::Vec3{1.0, 0.0, 0.0},
  });
  const auto rb = body.evaluate(state, sc);
  if (!finite_vec(rb.acceleration_mps2) || !approx(rb.direction_eci.x, 0.0) || !approx(rb.direction_eci.y, -1.0) ||
      !approx(rb.direction_eci.z, 0.0)) {
    spdlog::error("body-fixed direction mismatch");
    return 8;
  }

  const astroforces::forces::AntennaThrustPerturbationModel pert(p20, &sc, "antenna");
  const auto c = pert.evaluate(astroforces::forces::PerturbationRequest{.state = state, .spacecraft = &sc});
  if (c.status != astroforces::core::Status::Ok || c.type != astroforces::forces::PerturbationType::AntennaThrust) {
    spdlog::error("perturbation wrapper failed");
    return 9;
  }

  const astroforces::forces::AntennaThrustPerturbationModel pert_no_default(p20, nullptr, "antenna");
  const auto missing_sc = pert_no_default.evaluate(astroforces::forces::PerturbationRequest{.state = state, .spacecraft = nullptr});
  if (missing_sc.status != astroforces::core::Status::InvalidInput) {
    spdlog::error("missing spacecraft should fail");
    return 10;
  }

  return 0;
}


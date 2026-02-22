/**
 * @file test_earth_radiation.cpp
 * @brief Earth Radiation model unit/integration checks.
 * @author Watosn
 */

#include <cmath>

#include <spdlog/spdlog.h>

#include "astroforces/atmo/constants.hpp"
#include "astroforces/forces/surface/earth_radiation/earth_radiation_model.hpp"

namespace {

bool finite_vec(const astroforces::core::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

}  // namespace

int main() {
  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = 600.0, .reference_area_m2 = 4.0, .cd = 2.2, .cr = 1.3, .use_surface_model = false, .surfaces = {}};

  astroforces::core::StateVector near_state{};
  near_state.frame = astroforces::core::Frame::ECI;
  near_state.epoch.utc_seconds = 1.0e9;
  near_state.position_m = astroforces::core::Vec3{astroforces::core::constants::kEarthRadiusWgs84M + 400000.0, 0.0, 0.0};
  near_state.velocity_mps = astroforces::core::Vec3{0.0, 7670.0, 0.0};

  astroforces::forces::EarthRadiationAccelerationModel earth_radiation{};
  const auto near_out = earth_radiation.evaluate(near_state, sc);
  if (near_out.status != astroforces::core::Status::Ok) {
    spdlog::error("earth_radiation evaluate failed at near-state");
    return 1;
  }
  if (!finite_vec(near_out.acceleration_mps2) || !(astroforces::core::norm(near_out.acceleration_mps2) > 0.0)) {
    spdlog::error("earth_radiation acceleration invalid at near-state");
    return 2;
  }
  if (!(near_out.earth_radiation_pressure_pa > 0.0) || !(near_out.earth_distance_m > 0.0)) {
    spdlog::error("earth_radiation scalar outputs invalid at near-state");
    return 3;
  }
  if (std::abs(near_out.earth_radiation_pressure_pa - (near_out.albedo_pressure_pa + near_out.ir_pressure_pa)) > 1e-18) {
    spdlog::error("earth_radiation total pressure is not albedo+ir");
    return 7;
  }
  if (!(near_out.albedo_phase_function >= 0.0 && near_out.albedo_phase_function <= 1.0)) {
    spdlog::error("earth_radiation albedo phase is out of bounds: {}", near_out.albedo_phase_function);
    return 8;
  }
  if (!(near_out.albedo_eclipse_factor >= 0.0 && near_out.albedo_eclipse_factor <= 1.0)) {
    spdlog::error("earth_radiation albedo eclipse factor is out of bounds: {}", near_out.albedo_eclipse_factor);
    return 10;
  }

  auto far_state = near_state;
  far_state.position_m = 2.0 * near_state.position_m;
  const auto far_out = earth_radiation.evaluate(far_state, sc);
  if (far_out.status != astroforces::core::Status::Ok) {
    spdlog::error("earth_radiation evaluate failed at far-state");
    return 4;
  }
  const double near_mag = astroforces::core::norm(near_out.acceleration_mps2);
  const double far_mag = astroforces::core::norm(far_out.acceleration_mps2);
  const double ratio = far_mag / near_mag;
  if (std::abs(ratio - 0.25) > 1.0e-12) {
    spdlog::error("earth_radiation inverse-square scaling mismatch: ratio={}", ratio);
    return 5;
  }

  auto bad_state = near_state;
  bad_state.frame = astroforces::core::Frame::NED;
  const auto bad_out = earth_radiation.evaluate(bad_state, sc);
  if (bad_out.status != astroforces::core::Status::InvalidInput) {
    spdlog::error("earth_radiation invalid frame check failed");
    return 6;
  }

  const astroforces::forces::EarthRadiationAccelerationModel erp_ir_only({
      .use_albedo = false,
      .use_earth_ir = true,
  });
  const auto ir_only = erp_ir_only.evaluate(near_state, sc);
  if (ir_only.status != astroforces::core::Status::Ok || !(ir_only.ir_pressure_pa > 0.0) || ir_only.albedo_pressure_pa != 0.0) {
    spdlog::error("earth_radiation ir-only mode failed");
    return 9;
  }

  return 0;
}

/**
 * @file test_erp.cpp
 * @brief ERP model unit/integration checks.
 * @author Watosn
 */

#include <cmath>

#include <spdlog/spdlog.h>

#include "astroforces/atmo/constants.hpp"
#include "astroforces/forces/surface/erp/erp_model.hpp"

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

  astroforces::erp::ErpAccelerationModel erp{};
  const auto near_out = erp.evaluate(near_state, sc);
  if (near_out.status != astroforces::core::Status::Ok) {
    spdlog::error("erp evaluate failed at near-state");
    return 1;
  }
  if (!finite_vec(near_out.acceleration_mps2) || !(astroforces::core::norm(near_out.acceleration_mps2) > 0.0)) {
    spdlog::error("erp acceleration invalid at near-state");
    return 2;
  }
  if (!(near_out.earth_radiation_pressure_pa > 0.0) || !(near_out.earth_distance_m > 0.0)) {
    spdlog::error("erp scalar outputs invalid at near-state");
    return 3;
  }

  auto far_state = near_state;
  far_state.position_m = 2.0 * near_state.position_m;
  const auto far_out = erp.evaluate(far_state, sc);
  if (far_out.status != astroforces::core::Status::Ok) {
    spdlog::error("erp evaluate failed at far-state");
    return 4;
  }
  const double near_mag = astroforces::core::norm(near_out.acceleration_mps2);
  const double far_mag = astroforces::core::norm(far_out.acceleration_mps2);
  const double ratio = far_mag / near_mag;
  if (std::abs(ratio - 0.25) > 1.0e-12) {
    spdlog::error("erp inverse-square scaling mismatch: ratio={}", ratio);
    return 5;
  }

  auto bad_state = near_state;
  bad_state.frame = astroforces::core::Frame::NED;
  const auto bad_out = erp.evaluate(bad_state, sc);
  if (bad_out.status != astroforces::core::Status::InvalidInput) {
    spdlog::error("erp invalid frame check failed");
    return 6;
  }

  return 0;
}


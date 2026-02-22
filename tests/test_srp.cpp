/**
 * @file test_srp.cpp
 * @brief SRP model integration test.
 * @author Watosn
 */

#include <cmath>
#include <cstdlib>
#include <filesystem>

#include <spdlog/spdlog.h>

#include "astroforces/sc/spacecraft.hpp"
#include "astroforces/forces/surface/srp/srp_model.hpp"

namespace {

bool finite_vec(const astroforces::core::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

}  // namespace

int main() {
  namespace fs = std::filesystem;
  fs::path eph_path;
  if (const char* env = std::getenv("DRAGCPP_JPL_EPH_FILE")) {
    eph_path = env;
  } else {
    eph_path = fs::path(DRAGCPP_JPL_EPH_SOURCE_DIR) / "testdata" / "linux_p1550p2650.440";
  }
  if (!fs::exists(eph_path)) {
    spdlog::warn("srp test skipped: ephemeris not found: {}", eph_path.string());
    return 0;
  }

  astroforces::core::StateVector state{};
  state.frame = astroforces::core::Frame::ECI;
  state.epoch.utc_seconds = 1.0e9;
  state.position_m = astroforces::core::Vec3{6778137.0, 0.0, 0.0};
  state.velocity_mps = astroforces::core::Vec3{0.0, 7670.0, 0.0};

  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = 600.0, .reference_area_m2 = 4.0, .cd = 2.2, .cr = 1.3, .use_surface_model = false, .surfaces = {}};

  const auto srp = astroforces::srp::SrpAccelerationModel::Create({.ephemeris_file = eph_path, .use_eclipse = false});
  const auto out = srp->evaluate(state, sc);
  if (out.status != astroforces::core::Status::Ok) {
    spdlog::error("srp evaluate failed");
    return 1;
  }
  if (!finite_vec(out.acceleration_mps2) || !(astroforces::core::norm(out.acceleration_mps2) > 0.0)) {
    spdlog::error("srp acceleration invalid");
    return 2;
  }
  if (!(out.solar_pressure_pa > 0.0) || !(out.sun_distance_m > 0.0)) {
    spdlog::error("srp scalar outputs invalid");
    return 3;
  }

  const auto srp_e = astroforces::srp::SrpAccelerationModel::Create({.ephemeris_file = eph_path, .use_eclipse = true});
  const auto out_e = srp_e->evaluate(state, sc);
  if (out_e.status != astroforces::core::Status::Ok) {
    spdlog::error("srp eclipse evaluate failed");
    return 4;
  }

  return 0;
}

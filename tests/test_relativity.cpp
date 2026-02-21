/**
 * @file test_relativity.cpp
 * @brief Relativistic acceleration model tests.
 * @author Watosn
 */

#include <cmath>
#include <cstdlib>
#include <filesystem>

#include <spdlog/spdlog.h>

#include "astroforces/forces/relativity_model.hpp"

namespace {

bool finite_vec(const astroforces::atmo::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

double mag(const astroforces::atmo::Vec3& v) { return astroforces::atmo::norm(v); }

}  // namespace

int main() {
  astroforces::atmo::StateVector state{};
  state.frame = astroforces::atmo::Frame::ECI;
  state.epoch.utc_seconds = 1.0e9;
  state.position_m = astroforces::atmo::Vec3{6778137.0, 0.0, 0.0};
  state.velocity_mps = astroforces::atmo::Vec3{0.0, 7670.0, 0.0};

  const auto rel_no_geo = astroforces::forces::RelativityAccelerationModel::Create({.use_geodesic_precession = false});
  const auto out_no_geo = rel_no_geo->evaluate(state);
  if (out_no_geo.status != astroforces::atmo::Status::Ok) {
    spdlog::error("relativity no-geodesic evaluation failed");
    return 1;
  }
  if (!finite_vec(out_no_geo.acceleration_mps2) || !(mag(out_no_geo.acceleration_mps2) > 0.0)) {
    spdlog::error("relativity no-geodesic acceleration invalid");
    return 2;
  }

  const auto sum_terms = out_no_geo.spherical_central_body_mps2 + out_no_geo.lense_thirring_mps2 +
                         out_no_geo.oblateness_j2_mps2 + out_no_geo.rotational_energy_mps2;
  if (mag(out_no_geo.acceleration_mps2 - sum_terms) > 1e-18) {
    spdlog::error("relativity decomposition mismatch");
    return 3;
  }

  namespace fs = std::filesystem;
  fs::path eph_path;
  if (const char* env = std::getenv("DRAGCPP_JPL_EPH_FILE")) {
    eph_path = env;
  } else {
    eph_path = fs::path(DRAGCPP_JPL_EPH_SOURCE_DIR) / "testdata" / "linux_p1550p2650.440";
  }
  if (!fs::exists(eph_path)) {
    spdlog::warn("relativity geodesic test skipped: ephemeris not found: {}", eph_path.string());
    return 0;
  }
  const auto rel_with_geo =
      astroforces::forces::RelativityAccelerationModel::Create({.ephemeris_file = eph_path, .use_geodesic_precession = true});
  const auto out_with_geo = rel_with_geo->evaluate(state);
  if (out_with_geo.status != astroforces::atmo::Status::Ok) {
    spdlog::error("relativity geodesic evaluation failed");
    return 4;
  }
  if (!(mag(out_with_geo.geodesic_precession_mps2) > 0.0)) {
    spdlog::error("relativity geodesic term is zero");
    return 5;
  }
  return 0;
}


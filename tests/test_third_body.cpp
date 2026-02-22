/**
 * @file test_third_body.cpp
 * @brief Third-body perturbation integration tests.
 * @author Watosn
 */

#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <filesystem>

#include <spdlog/spdlog.h>

#include "astroforces/core/types.hpp"
#include "astroforces/forces/gravity/third_body.hpp"

namespace {

bool finite_vec(const astroforces::core::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

astroforces::core::Vec3 add(const astroforces::core::Vec3& a, const astroforces::core::Vec3& b) {
  return astroforces::core::Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

}  // namespace

int main() {
  namespace fs = std::filesystem;
  fs::path eph_path;
  if (const char* env = std::getenv("ASTROFORCES_JPL_EPH_FILE")) {
    eph_path = env;
  } else {
    eph_path = fs::path(ASTROFORCES_JPL_EPH_SOURCE_DIR) / "testdata" / "linux_p1550p2650.440";
    if (!fs::exists(eph_path)) {
      eph_path = fs::path(ASTROFORCES_SOURCE_DIR) / "data" / "required" / "linux_p1550p2650.440";
    }
  }
  if (!fs::exists(eph_path)) {
#if defined(ASTROFORCES_TEST_REQUIRE_EXTERNAL_DATA)
    spdlog::error("third-body test requires ephemeris but file not found: {}", eph_path.string());
    return 100;
#else
    spdlog::warn("third-body test skipped: ephemeris not found: {}", eph_path.string());
    return 0;
#endif
  }

  astroforces::core::StateVector state{};
  state.frame = astroforces::core::Frame::ECI;
  state.epoch.utc_seconds = 1.0e9;
  state.position_m = astroforces::core::Vec3{6778137.0, 0.0, 0.0};
  state.velocity_mps = astroforces::core::Vec3{0.0, 7670.0, 0.0};

  const auto both = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_path, .use_sun = true, .use_moon = true, .name = "third_both"});
  const auto sun = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_path, .use_sun = true, .use_moon = false, .name = "third_sun"});
  const auto moon = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_path, .use_sun = false, .use_moon = true, .name = "third_moon"});
  const auto both_no_eq79 = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_path,
       .use_sun = true,
       .use_moon = true,
       .use_goce_eq79_indirect_j2 = false,
       .name = "third_both_no_eq79"});

  const auto req = astroforces::forces::PerturbationRequest{.state = state, .spacecraft = nullptr};
  const auto c_both = both->evaluate(req);
  const auto c_sun = sun->evaluate(req);
  const auto c_moon = moon->evaluate(req);
  if (c_both.status != astroforces::core::Status::Ok || c_sun.status != astroforces::core::Status::Ok ||
      c_moon.status != astroforces::core::Status::Ok) {
    spdlog::error("third-body evaluation failed");
    return 1;
  }
  if (!finite_vec(c_both.acceleration_mps2) || !finite_vec(c_sun.acceleration_mps2) || !finite_vec(c_moon.acceleration_mps2)) {
    spdlog::error("third-body acceleration not finite");
    return 2;
  }
  if (!(astroforces::core::norm(c_both.acceleration_mps2) > 0.0) || !(astroforces::core::norm(c_sun.acceleration_mps2) > 0.0) ||
      !(astroforces::core::norm(c_moon.acceleration_mps2) > 0.0)) {
    spdlog::error("third-body acceleration unexpectedly zero");
    return 3;
  }

  const auto sum = add(c_sun.acceleration_mps2, c_moon.acceleration_mps2);
  const auto diff = astroforces::core::Vec3{c_both.acceleration_mps2.x - sum.x, c_both.acceleration_mps2.y - sum.y,
                                        c_both.acceleration_mps2.z - sum.z};
  if (!(astroforces::core::norm(diff) <= 1e-16 * std::max(1.0, astroforces::core::norm(c_both.acceleration_mps2)))) {
    spdlog::error("third-body superposition mismatch");
    return 4;
  }

  const auto c_both_no_eq79 = both_no_eq79->evaluate(req);
  if (c_both_no_eq79.status != astroforces::core::Status::Ok) {
    spdlog::error("third-body evaluation without eq79 failed");
    return 5;
  }
  const auto eq79_delta = astroforces::core::Vec3{
      c_both.acceleration_mps2.x - c_both_no_eq79.acceleration_mps2.x,
      c_both.acceleration_mps2.y - c_both_no_eq79.acceleration_mps2.y,
      c_both.acceleration_mps2.z - c_both_no_eq79.acceleration_mps2.z};
  if (!(astroforces::core::norm(eq79_delta) > 0.0)) {
    spdlog::error("goce eq79 indirect j2 term was not applied");
    return 6;
  }

  state.frame = astroforces::core::Frame::ECEF;
  const auto c_bad_frame = both->evaluate(astroforces::forces::PerturbationRequest{.state = state, .spacecraft = nullptr});
  if (c_bad_frame.status != astroforces::core::Status::InvalidInput) {
    spdlog::error("expected invalid frame error");
    return 7;
  }

  return 0;
}

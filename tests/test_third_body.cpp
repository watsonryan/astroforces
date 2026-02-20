/**
 * @file test_third_body.cpp
 * @brief Third-body perturbation integration tests.
 * @author Watosn
 */

#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <iostream>

#include "dragcpp/atmo/types.hpp"
#include "dragcpp/forces/third_body.hpp"

namespace {

bool finite_vec(const astroforces::atmo::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

astroforces::atmo::Vec3 add(const astroforces::atmo::Vec3& a, const astroforces::atmo::Vec3& b) {
  return astroforces::atmo::Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
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
    std::cout << "third-body test skipped: ephemeris not found: " << eph_path << "\n";
    return 0;
  }

  astroforces::atmo::StateVector state{};
  state.frame = astroforces::atmo::Frame::ECI;
  state.epoch.utc_seconds = 1.0e9;
  state.position_m = astroforces::atmo::Vec3{6778137.0, 0.0, 0.0};
  state.velocity_mps = astroforces::atmo::Vec3{0.0, 7670.0, 0.0};

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
  if (c_both.status != astroforces::atmo::Status::Ok || c_sun.status != astroforces::atmo::Status::Ok ||
      c_moon.status != astroforces::atmo::Status::Ok) {
    std::cerr << "third-body evaluation failed\n";
    return 1;
  }
  if (!finite_vec(c_both.acceleration_mps2) || !finite_vec(c_sun.acceleration_mps2) || !finite_vec(c_moon.acceleration_mps2)) {
    std::cerr << "third-body acceleration not finite\n";
    return 2;
  }
  if (!(astroforces::atmo::norm(c_both.acceleration_mps2) > 0.0) || !(astroforces::atmo::norm(c_sun.acceleration_mps2) > 0.0) ||
      !(astroforces::atmo::norm(c_moon.acceleration_mps2) > 0.0)) {
    std::cerr << "third-body acceleration unexpectedly zero\n";
    return 3;
  }

  const auto sum = add(c_sun.acceleration_mps2, c_moon.acceleration_mps2);
  const auto diff = astroforces::atmo::Vec3{c_both.acceleration_mps2.x - sum.x, c_both.acceleration_mps2.y - sum.y,
                                        c_both.acceleration_mps2.z - sum.z};
  if (!(astroforces::atmo::norm(diff) <= 1e-16 * std::max(1.0, astroforces::atmo::norm(c_both.acceleration_mps2)))) {
    std::cerr << "third-body superposition mismatch\n";
    return 4;
  }

  const auto c_both_no_eq79 = both_no_eq79->evaluate(req);
  if (c_both_no_eq79.status != astroforces::atmo::Status::Ok) {
    std::cerr << "third-body evaluation without eq79 failed\n";
    return 5;
  }
  const auto eq79_delta = astroforces::atmo::Vec3{
      c_both.acceleration_mps2.x - c_both_no_eq79.acceleration_mps2.x,
      c_both.acceleration_mps2.y - c_both_no_eq79.acceleration_mps2.y,
      c_both.acceleration_mps2.z - c_both_no_eq79.acceleration_mps2.z};
  if (!(astroforces::atmo::norm(eq79_delta) > 0.0)) {
    std::cerr << "goce eq79 indirect j2 term was not applied\n";
    return 6;
  }

  state.frame = astroforces::atmo::Frame::ECEF;
  const auto c_bad_frame = both->evaluate(astroforces::forces::PerturbationRequest{.state = state, .spacecraft = nullptr});
  if (c_bad_frame.status != astroforces::atmo::Status::InvalidInput) {
    std::cerr << "expected invalid frame error\n";
    return 7;
  }

  return 0;
}

/**
 * @file test_gravity_sph.cpp
 * @brief Full SPH gravity model tests.
 * @author Watosn
 */

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include <spdlog/spdlog.h>

#include "astroforces/forces/gravity/gravity_sph_model.hpp"

namespace {

bool finite_vec(const astroforces::core::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

std::filesystem::path write_test_gfc() {
  namespace fs = std::filesystem;
  const fs::path p = fs::temp_directory_path() / "astroforces_test_gm.gfc";
  std::ofstream out(p);
  out << "modelname TEST\n";
  out << "earth_gravity_constant 3.986004418e14\n";
  out << "radius 6378136.3\n";
  out << "max_degree 4\n";
  out << "tide_system tide_free\n";
  out << "end_of_head\n";
  out << "gfc 0 0 1.0 0.0 0.0 0.0\n";
  out << "gfc 2 0 -4.84165143790815e-4 0.0 0.0 0.0\n";
  out << "gfc 2 1 1.0e-10 -2.0e-10 0.0 0.0\n";
  out << "gfc 2 2 2.0e-6 -1.0e-6 0.0 0.0\n";
  out << "gfc 3 0 9.57161e-7 0.0 0.0 0.0\n";
  out << "gfc 3 1 2.0e-7 1.0e-7 0.0 0.0\n";
  out << "gfc 3 2 1.0e-7 -2.0e-7 0.0 0.0\n";
  out << "gfc 3 3 5.0e-8 4.0e-8 0.0 0.0\n";
  out.close();
  return p;
}

}  // namespace

int main() {
  namespace fs = std::filesystem;
  const auto gravity_file = write_test_gfc();

  fs::path eph_path;
  if (const char* env = std::getenv("DRAGCPP_JPL_EPH_FILE")) {
    eph_path = env;
  } else {
    eph_path = fs::path(DRAGCPP_JPL_EPH_SOURCE_DIR) / "testdata" / "linux_p1550p2650.440";
  }

  astroforces::core::StateVector state{};
  state.frame = astroforces::core::Frame::ECEF;
  state.epoch.utc_seconds = 1.0e9;
  state.position_m = astroforces::core::Vec3{6878137.0, 1200.0, -3400.0};
  state.velocity_mps = astroforces::core::Vec3{0.0, 7600.0, 0.0};

  const auto central_only = astroforces::forces::GravitySphAccelerationModel::Create(
      {.gravity_model_file = gravity_file, .max_degree = 4, .use_central = true, .use_sph = false, .use_solid_earth_tides = false});
  const auto full_sph = astroforces::forces::GravitySphAccelerationModel::Create(
      {.gravity_model_file = gravity_file, .max_degree = 4, .use_central = true, .use_sph = true, .use_solid_earth_tides = false});

  const auto out_c = central_only->evaluate(state);
  const auto out_s = full_sph->evaluate(state);
  if (out_c.status != astroforces::core::Status::Ok || out_s.status != astroforces::core::Status::Ok) {
    spdlog::error("gravity evaluation failed");
    return 1;
  }
  if (!finite_vec(out_c.acceleration_mps2) || !finite_vec(out_s.acceleration_mps2)) {
    spdlog::error("gravity acceleration not finite");
    return 2;
  }

  const auto sph_delta = astroforces::core::Vec3{
      out_s.acceleration_mps2.x - out_c.acceleration_mps2.x,
      out_s.acceleration_mps2.y - out_c.acceleration_mps2.y,
      out_s.acceleration_mps2.z - out_c.acceleration_mps2.z};
  if (!(astroforces::core::norm(sph_delta) > 0.0)) {
    spdlog::error("sph harmonics had no effect");
    return 3;
  }

  state.frame = astroforces::core::Frame::ECI;
  const auto out_eci = full_sph->evaluate(state);
  if (out_eci.status != astroforces::core::Status::Ok || !finite_vec(out_eci.acceleration_mps2)) {
    spdlog::error("eci evaluation failed");
    return 4;
  }

  if (fs::exists(eph_path)) {
    const auto no_tides = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .ephemeris_file = eph_path,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = false});
    const auto with_tides = astroforces::forces::GravitySphAccelerationModel::Create(
        {.gravity_model_file = gravity_file,
         .ephemeris_file = eph_path,
         .max_degree = 4,
         .use_central = true,
         .use_sph = true,
         .use_solid_earth_tides = true});
    const auto out_no_tides = no_tides->evaluate(state);
    const auto out_with_tides = with_tides->evaluate(state);
    if (out_no_tides.status != astroforces::core::Status::Ok || out_with_tides.status != astroforces::core::Status::Ok) {
      spdlog::error("tide model evaluation failed");
      return 5;
    }
    const auto tide_delta = astroforces::core::Vec3{
        out_with_tides.acceleration_mps2.x - out_no_tides.acceleration_mps2.x,
        out_with_tides.acceleration_mps2.y - out_no_tides.acceleration_mps2.y,
        out_with_tides.acceleration_mps2.z - out_no_tides.acceleration_mps2.z};
    if (!(astroforces::core::norm(tide_delta) > 0.0)) {
      spdlog::error("solid Earth tides had no effect");
      return 6;
    }
  }

  state.frame = astroforces::core::Frame::NED;
  const auto out_bad = full_sph->evaluate(state);
  if (out_bad.status != astroforces::core::Status::InvalidInput) {
    spdlog::error("expected invalid frame status");
    return 7;
  }

  return 0;
}

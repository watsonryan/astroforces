/**
 * @file gravity_sph_cli.cpp
 * @brief Single-state full SPH gravity CLI.
 * @author Watosn
 */

#include <cstdlib>
#include <filesystem>
#include <string>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/gravity/gravity_sph_model.hpp"

namespace {

double magnitude(const astroforces::core::Vec3& v) { return astroforces::core::norm(v); }

}  // namespace

int main(int argc, char** argv) {
  if (argc < 9 || argc > 13) {
    spdlog::error(
        "usage: gravity_sph_cli <x_m> <y_m> <z_m> <vx_mps> <vy_mps> <vz_mps> <epoch_utc_s> <gravity_gfc_file> [max_degree] [frame:eci|ecef] [jpl_ephemeris_file] [use_tides:0|1]");
    return 1;
  }

  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])};
  state.velocity_mps = astroforces::core::Vec3{std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
  state.epoch.utc_seconds = std::atof(argv[7]);

  const std::filesystem::path gravity_file = argv[8];
  if (!std::filesystem::exists(gravity_file)) {
    spdlog::error("gravity coefficient file not found: {}", gravity_file.string());
    return 2;
  }

  const int max_degree = (argc >= 10) ? std::atoi(argv[9]) : 360;
  const std::string frame = (argc >= 11) ? std::string(argv[10]) : "eci";
  if (frame == "ecef") {
    state.frame = astroforces::core::Frame::ECEF;
  } else if (frame == "eci") {
    state.frame = astroforces::core::Frame::ECI;
  } else {
    spdlog::error("invalid frame: {} (expected eci|ecef)", frame);
    return 3;
  }

  const std::filesystem::path eph_file = (argc >= 12) ? std::filesystem::path(argv[11]) : std::filesystem::path{};
  const bool use_tides = (argc >= 13) ? (std::atoi(argv[12]) != 0) : false;
  if (use_tides && (eph_file.empty() || !std::filesystem::exists(eph_file))) {
    spdlog::error("tides requested but ephemeris file not found: {}", eph_file.string());
    return 4;
  }

  const auto grav = astroforces::forces::GravitySphAccelerationModel::Create(
      {.gravity_model_file = gravity_file,
       .ephemeris_file = eph_file,
       .max_degree = max_degree,
       .use_central = true,
       .use_sph = true,
       .use_solid_earth_tides = use_tides});

  const auto out = grav->evaluate(state);
  if (out.status != astroforces::core::Status::Ok) {
    spdlog::error("gravity evaluation failed: status={}", static_cast<int>(out.status));
    return 5;
  }

  fmt::print("ax={} ay={} az={} amag={}\n", out.acceleration_mps2.x, out.acceleration_mps2.y, out.acceleration_mps2.z,
             magnitude(out.acceleration_mps2));
  fmt::print("central={} sph={} tide_sun={} tide_moon={}\n", magnitude(out.central_mps2), magnitude(out.sph_mps2),
             magnitude(out.solid_tide_sun_mps2), magnitude(out.solid_tide_moon_mps2));
  return 0;
}

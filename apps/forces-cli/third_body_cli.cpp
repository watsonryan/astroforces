/**
 * @file third_body_cli.cpp
 * @brief Single-state third-body perturbation CLI.
 * @author Watosn
 */

#include <cstdlib>
#include <filesystem>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/gravity/third_body.hpp"

namespace {

double magnitude(const astroforces::core::Vec3& v) { return astroforces::core::norm(v); }

}  // namespace

int main(int argc, char** argv) {
  if (argc < 9 || argc > 11) {
    spdlog::error(
        "usage: third_body_cli <x_eci_m> <y_eci_m> <z_eci_m> <vx_eci_mps> <vy_eci_mps> <vz_eci_mps> <epoch_utc_s> <jpl_ephemeris_file> [use_sun:0|1] [use_moon:0|1]");
    return 1;
  }

  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])};
  state.velocity_mps = astroforces::core::Vec3{std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
  state.epoch.utc_seconds = std::atof(argv[7]);
  state.frame = astroforces::core::Frame::ECI;

  const std::filesystem::path eph_file = argv[8];
  if (!std::filesystem::exists(eph_file)) {
    spdlog::error("ephemeris file not found: {}", eph_file.string());
    return 2;
  }

  const bool use_sun = (argc >= 10) ? (std::atoi(argv[9]) != 0) : true;
  const bool use_moon = (argc >= 11) ? (std::atoi(argv[10]) != 0) : true;

  const auto sun = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_file, .use_sun = use_sun, .use_moon = false, .name = "third_body_sun"});
  const auto moon = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_file, .use_sun = false, .use_moon = use_moon, .name = "third_body_moon"});
  const auto total = astroforces::forces::ThirdBodyPerturbationModel::Create(
      {.ephemeris_file = eph_file, .use_sun = use_sun, .use_moon = use_moon, .name = "third_body_total"});

  const auto req = astroforces::forces::PerturbationRequest{.state = state, .spacecraft = nullptr};
  const auto sun_c = sun->evaluate(req);
  const auto moon_c = moon->evaluate(req);
  const auto total_c = total->evaluate(req);

  if (sun_c.status != astroforces::core::Status::Ok || moon_c.status != astroforces::core::Status::Ok ||
      total_c.status != astroforces::core::Status::Ok) {
    spdlog::error("third-body evaluation failed: sun={}, moon={}, total={}", static_cast<int>(sun_c.status),
                  static_cast<int>(moon_c.status), static_cast<int>(total_c.status));
    return 3;
  }

  fmt::print("sun_ax={} sun_ay={} sun_az={} sun_mag={}\n", sun_c.acceleration_mps2.x, sun_c.acceleration_mps2.y,
             sun_c.acceleration_mps2.z, magnitude(sun_c.acceleration_mps2));
  fmt::print("moon_ax={} moon_ay={} moon_az={} moon_mag={}\n", moon_c.acceleration_mps2.x, moon_c.acceleration_mps2.y,
             moon_c.acceleration_mps2.z, magnitude(moon_c.acceleration_mps2));
  fmt::print("total_ax={} total_ay={} total_az={} total_mag={}\n", total_c.acceleration_mps2.x, total_c.acceleration_mps2.y,
             total_c.acceleration_mps2.z, magnitude(total_c.acceleration_mps2));

  return 0;
}


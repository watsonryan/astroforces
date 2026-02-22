/**
 * @file relativity_cli.cpp
 * @brief Single-state relativistic perturbation CLI.
 * @author Watosn
 */

#include <cstdlib>
#include <filesystem>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/gravity/relativity_model.hpp"

namespace {

double magnitude(const astroforces::core::Vec3& v) { return astroforces::core::norm(v); }

}  // namespace

int main(int argc, char** argv) {
  if (argc < 9 || argc > 10) {
    spdlog::error(
        "usage: relativity_cli <x_eci_m> <y_eci_m> <z_eci_m> <vx_eci_mps> <vy_eci_mps> <vz_eci_mps> <epoch_utc_s> <jpl_ephemeris_file> [use_geodesic:0|1]");
    return 1;
  }

  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])};
  state.velocity_mps = astroforces::core::Vec3{std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
  state.epoch.utc_seconds = std::atof(argv[7]);
  state.frame = astroforces::core::Frame::ECI;

  const std::filesystem::path eph_file = argv[8];
  const bool use_geodesic = (argc >= 10) ? (std::atoi(argv[9]) != 0) : true;
  if (use_geodesic && !std::filesystem::exists(eph_file)) {
    spdlog::error("ephemeris file not found: {}", eph_file.string());
    return 2;
  }

  const auto rel = astroforces::forces::RelativityAccelerationModel::Create(
      {.ephemeris_file = eph_file, .use_geodesic_precession = use_geodesic});
  const auto out = rel->evaluate(state);
  if (out.status != astroforces::core::Status::Ok) {
    spdlog::error("relativity evaluation failed: status={}", static_cast<int>(out.status));
    return 3;
  }

  fmt::print("ax={} ay={} az={} amag={}\n", out.acceleration_mps2.x, out.acceleration_mps2.y, out.acceleration_mps2.z,
             magnitude(out.acceleration_mps2));
  fmt::print("schw={} geo={} lt={} j2={} rot={}\n", magnitude(out.spherical_central_body_mps2), magnitude(out.geodesic_precession_mps2),
             magnitude(out.lense_thirring_mps2), magnitude(out.oblateness_j2_mps2), magnitude(out.rotational_energy_mps2));
  return 0;
}


/**
 * @file srp_cli.cpp
 * @brief Single-state SRP perturbation CLI.
 * @author Watosn
 */

#include <cstdlib>
#include <filesystem>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/srp/srp_model.hpp"

namespace {

double magnitude(const astroforces::core::Vec3& v) { return astroforces::core::norm(v); }

}  // namespace

int main(int argc, char** argv) {
  if (argc < 9 || argc > 13) {
    spdlog::error(
        "usage: srp_cli <x_eci_m> <y_eci_m> <z_eci_m> <vx_eci_mps> <vy_eci_mps> <vz_eci_mps> <epoch_utc_s> <jpl_ephemeris_file> [mass_kg] [area_m2] [cr] [use_eclipse:0|1]");
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

  const double mass_kg = (argc >= 10) ? std::atof(argv[9]) : 600.0;
  const double area_m2 = (argc >= 11) ? std::atof(argv[10]) : 4.0;
  const double cr = (argc >= 12) ? std::atof(argv[11]) : 1.3;
  const bool use_eclipse = (argc >= 13) ? (std::atoi(argv[12]) != 0) : false;

  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = mass_kg, .reference_area_m2 = area_m2, .cd = 2.2, .cr = cr, .use_surface_model = false, .surfaces = {}};

  auto srp = astroforces::srp::SrpAccelerationModel::Create(
      {.ephemeris_file = eph_file, .use_eclipse = use_eclipse});
  const auto out = srp->evaluate(state, sc);
  if (out.status != astroforces::core::Status::Ok) {
    spdlog::error("srp evaluation failed: status={}", static_cast<int>(out.status));
    return 3;
  }

  fmt::print("ax={} ay={} az={} amag={}\n", out.acceleration_mps2.x, out.acceleration_mps2.y, out.acceleration_mps2.z,
             magnitude(out.acceleration_mps2));
  fmt::print("p_pa={} r_sun_m={} area={} cr={} eclipsed={}\n", out.solar_pressure_pa, out.sun_distance_m, out.area_m2, out.cr,
             out.eclipsed ? 1 : 0);
  return 0;
}


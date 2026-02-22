/**
 * @file earth_radiation_cli.cpp
 * @brief Single-state Earth radiation pressure perturbation CLI.
 * @author Watosn
 */

#include <cstdlib>
#include <string>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/earth_radiation/earth_radiation_model.hpp"

namespace {

double magnitude(const astroforces::core::Vec3& v) { return astroforces::core::norm(v); }

}  // namespace

int main(int argc, char** argv) {
  if (argc < 8 || argc > 12) {
    spdlog::error(
        "usage: earth_radiation_cli <x_frame_m> <y_frame_m> <z_frame_m> <vx_frame_mps> <vy_frame_mps> <vz_frame_mps> <epoch_utc_s> [mass_kg] [area_m2] [cr] [jpl_ephemeris_file]");
    return 1;
  }

  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])};
  state.velocity_mps = astroforces::core::Vec3{std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
  state.epoch.utc_seconds = std::atof(argv[7]);
  state.frame = astroforces::core::Frame::ECI;

  const double mass_kg = (argc >= 9) ? std::atof(argv[8]) : 600.0;
  const double area_m2 = (argc >= 10) ? std::atof(argv[9]) : 4.0;
  const double cr = (argc >= 11) ? std::atof(argv[10]) : 1.3;
  const std::string eph_file = (argc >= 12) ? argv[11] : "";

  astroforces::sc::SpacecraftProperties sc{
      .mass_kg = mass_kg, .reference_area_m2 = area_m2, .cd = 2.2, .cr = cr, .use_surface_model = false, .surfaces = {}};

  auto earth_radiation = astroforces::forces::EarthRadiationAccelerationModel::Create({
      .ephemeris_file = eph_file,
  });
  if (!earth_radiation) {
    spdlog::error("failed to create Earth Radiation model");
    return 2;
  }
  const auto out = earth_radiation->evaluate(state, sc);
  if (out.status != astroforces::core::Status::Ok) {
    spdlog::error("earth_radiation evaluation failed: status={}", static_cast<int>(out.status));
    return 3;
  }

  fmt::print("ax={} ay={} az={} amag={}\n", out.acceleration_mps2.x, out.acceleration_mps2.y, out.acceleration_mps2.z,
             magnitude(out.acceleration_mps2));
  fmt::print("p_total_pa={} p_albedo_pa={} p_ir_pa={} albedo_phase={} albedo_eclipse_factor={} r_earth_m={} area={} cr={}\n",
             out.earth_radiation_pressure_pa, out.albedo_pressure_pa, out.ir_pressure_pa, out.albedo_phase_function,
             out.albedo_eclipse_factor, out.earth_distance_m, out.area_m2, out.cr);
  return 0;
}

/**
 * @file test_drag_core.cpp
 * @brief Core drag equation tests.
 * @author Watosn
 */

#include <cmath>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <random>
#include <sstream>
#include <string>

#include <spdlog/spdlog.h>

#include "astroforces/core/transforms.hpp"
#include "astroforces/core/eop.hpp"
#include "astroforces/core/cip.hpp"
#include "astroforces/forces/surface/drag/drag_model.hpp"
#include "astroforces/models/exponential_atmosphere.hpp"
#include "astroforces/sc/spacecraft.hpp"
#include "astroforces/weather/static_provider.hpp"

namespace {

bool approx(double a, double b, double rel) {
  const double d = std::abs(a - b);
  const double n = std::max(std::abs(b), 1e-30);
  return d / n <= rel;
}

std::filesystem::path make_fixed_column_eop_file() {
  const auto path = std::filesystem::temp_directory_path() / ("astroforces_drag_eop_test_" + std::to_string(
                                                                   std::chrono::high_resolution_clock::now().time_since_epoch().count())
                                                               + "_" + std::to_string(std::random_device{}()) + ".txt");
  std::ofstream out(path);

  auto make_line = [](double mjd,
                      double xp_arcsec,
                      double yp_arcsec,
                      double dut1_s,
                      double lod_ms,
                      double dX_mas,
                      double dY_mas) {
    std::string line(140, ' ');
    auto put_num = [&](const int start, const int width, const double value, const int prec) {
      std::ostringstream oss;
      oss << std::fixed << std::setw(width) << std::setprecision(prec) << value;
      const std::string s = oss.str();
      line.replace(static_cast<std::size_t>(start), static_cast<std::size_t>(width), s);
    };
    put_num(7, 8, mjd, 2);
    put_num(18, 9, xp_arcsec, 6);
    put_num(37, 9, yp_arcsec, 6);
    put_num(58, 10, dut1_s, 7);
    put_num(79, 7, lod_ms, 3);
    put_num(97, 9, dX_mas, 6);
    put_num(116, 9, dY_mas, 6);
    return line;
  };

  out << make_line(60676.0, 0.123456, -0.234567, 0.1000000, 0.900, 0.345678, -0.456789) << "\n";
  out << make_line(60677.0, 0.223456, -0.134567, 0.1200000, 1.100, 0.445678, -0.356789) << "\n";
  return path;
}

std::filesystem::path make_cip_file() {
  const auto path = std::filesystem::temp_directory_path() / ("astroforces_drag_cip_test_" + std::to_string(
                                                                   std::chrono::high_resolution_clock::now().time_since_epoch().count())
                                                               + "_" + std::to_string(std::random_device{}()) + ".txt");
  std::ofstream out(path);
  out << "# MJD_UTC X Y s (arcsec)\n";
  out << "60676.0 0.100 -0.200 0.010\n";
  out << "60677.0 0.120 -0.220 0.015\n";
  return path;
}

}  // namespace

int main() {
  using namespace astroforces;
  const core::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0, .status = core::Status::Ok};
  weather::StaticSpaceWeatherProvider weather(wx);
  models::ExponentialAtmosphereModel atmosphere(1.225, 0.0, 7000.0, 1000.0);
  models::ZeroWindModel wind;
  forces::DragAccelerationModel model(weather, atmosphere, wind);

  const double epoch_utc_s = (60676.5 + 2400000.5 - 2440587.5) * core::constants::kSecondsPerDay;

  core::StateVector state{};
  state.epoch.utc_seconds = epoch_utc_s;
  state.frame = core::Frame::ECEF;
  state.position_m = core::Vec3{6378137.0, 0.0, 0.0};
  state.velocity_mps = core::Vec3{7500.0, 0.0, 0.0};

  sc::SpacecraftProperties sc{.mass_kg = 1000.0, .reference_area_m2 = 10.0, .cd = 2.0, .use_surface_model = false};

  const auto out = model.evaluate(state, sc);
  if (out.status != core::Status::Ok) {
    spdlog::error("status failed");
    return 1;
  }

  const double expected_ax = -0.5 * 1.225 * 2.0 * 10.0 / 1000.0 * 7500.0 * 7500.0;
  if (!approx(out.acceleration_mps2.x, expected_ax, 1e-12)) {
    spdlog::error("ax mismatch");
    return 2;
  }
  if (!approx(out.acceleration_mps2.y, 0.0, 1e-12) || !approx(out.acceleration_mps2.z, 0.0, 1e-12)) {
    spdlog::error("vector mismatch");
    return 3;
  }
  if (!approx(out.relative_speed_mps, 7500.0, 1e-12) || !approx(out.dynamic_pressure_pa, 0.5 * 1.225 * 7500.0 * 7500.0, 1e-12)) {
    spdlog::error("derived drag scalars mismatch");
    return 7;
  }

  core::StateVector state_eci{};
  state_eci.epoch = state.epoch;
  state_eci.frame = core::Frame::ECI;
  state_eci.position_m = core::Vec3{state.position_m.x, state.position_m.y, state.position_m.z};
  state_eci.velocity_mps = core::Vec3{state.velocity_mps.x, state.velocity_mps.y, state.velocity_mps.z};
  state_eci.body_from_frame_dcm = state.body_from_frame_dcm;
  const auto out_eci = model.evaluate(state_eci, sc);
  if (out_eci.status != core::Status::Ok || !std::isfinite(out_eci.relative_speed_mps)) {
    spdlog::error("eci drag evaluation failed");
    return 9;
  }

  const auto eop_path = make_fixed_column_eop_file();
  const auto cip_path = make_cip_file();
  const auto eop_series = core::eop::Series::load_iers_finals(eop_path);
  const auto cip_series = core::cip::Series::load_table(cip_path, core::cip::AngleUnit::Arcseconds, false);
  std::error_code ec;
  std::filesystem::remove(eop_path, ec);
  std::filesystem::remove(cip_path, ec);

  if (eop_series.empty() || cip_series.empty()) {
    spdlog::error("failed to load strict transform support data");
    return 10;
  }

  forces::DragAccelerationModel strict_model(
      weather,
      atmosphere,
      wind,
      forces::DragFrameTransformMode::StrictGcrfItrf,
      &eop_series,
      &cip_series);

  const auto out_eci_strict = strict_model.evaluate(state_eci, sc);
  if (out_eci_strict.status != core::Status::Ok || !std::isfinite(out_eci_strict.relative_speed_mps)) {
    spdlog::error("strict eci drag evaluation failed");
    return 11;
  }

  sc::SpacecraftProperties macro_sc{
      .mass_kg = 1000.0,
      .reference_area_m2 = 99.0,
      .cd = 2.0,
      .use_surface_model = true,
      .surfaces = {
          sc::Surface{.normal_body = core::Vec3{-1.0, 0.0, 0.0}, .area_m2 = 2.0, .cd = 2.5},
          sc::Surface{.normal_body = core::Vec3{0.0, -1.0, 0.0}, .area_m2 = 4.0, .cd = 3.5},
      },
  };

  core::StateVector state_identity = state;
  state_identity.body_from_frame_dcm = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  const auto macro_identity = model.evaluate(state_identity, macro_sc);
  if (macro_identity.status != core::Status::Ok || !approx(macro_identity.area_m2, 2.0, 1e-12) || !approx(macro_identity.cd, 2.5, 1e-12)) {
    spdlog::error("macro identity area mismatch");
    return 4;
  }

  core::StateVector state_rot = state;
  state_rot.body_from_frame_dcm = {0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};  // +90deg yaw
  const auto macro_rot = model.evaluate(state_rot, macro_sc);
  if (macro_rot.status != core::Status::Ok || !approx(macro_rot.area_m2, 4.0, 1e-12) || !approx(macro_rot.cd, 3.5, 1e-12)) {
    spdlog::error("macro rotated area mismatch");
    return 5;
  }

  if (!approx(out.area_m2, 10.0, 1e-12)) {
    spdlog::error("cannonball area mismatch");
    return 6;
  }

  sc::SpacecraftProperties aero_sc{
      .mass_kg = 1000.0,
      .reference_area_m2 = 99.0,
      .cd = 2.0,
      .use_surface_model = true,
      .surfaces = {
          sc::Surface{.normal_body = core::Vec3{-1.0, 0.0, 0.0}, .area_m2 = 2.0, .cd = 2.0, .specularity = 0.0, .accommodation = 1.0},
      },
  };
  const auto aero_out = model.evaluate(state_identity, aero_sc);
  if (aero_out.status != core::Status::Ok || !(aero_out.cd > 2.0)) {
    spdlog::error("surface aero modifier mismatch");
    return 8;
  }

  return 0;
}

/**
 * @file test_adapters_integration.cpp
 * @brief Adapter-backed integration tests across NRLMSIS, DTM2020, and HWM14.
 * @author Watosn
 */

#include <cmath>
#include <algorithm>
#include <filesystem>

#include <spdlog/spdlog.h>

#include "astroforces/adapters/dtm2020_adapter.hpp"
#include "astroforces/adapters/hwm14_adapter.hpp"
#include "astroforces/adapters/nrlmsis21_adapter.hpp"
#include "astroforces/atmo/conversions.hpp"
#include "astroforces/forces/surface/drag/drag_model.hpp"
#include "astroforces/sc/spacecraft.hpp"
#include "astroforces/weather/static_provider.hpp"
#include "msis21/msis21.hpp"

namespace {

bool finite(const astroforces::core::Vec3& v) {
  return std::isfinite(v.x) && std::isfinite(v.y) && std::isfinite(v.z);
}

bool approx_rel(double a, double b, double rel = 1e-12) {
  const double d = std::abs(a - b);
  const double n = std::max(std::abs(b), 1e-30);
  return d / n <= rel;
}

}  // namespace

int main() {
  using namespace astroforces;
  namespace fs = std::filesystem;

  const auto nrl_parm = fs::path(DRAGCPP_NRLMSIS21_SOURCE_DIR) / "data" / "msis21.parm";
  const auto dtm_coeff = fs::path(DRAGCPP_DTM2020_SOURCE_DIR) / "testdata" / "operational_regression_coeff.dat";
  const auto hwm_data_dir = fs::path(DRAGCPP_HWM14_SOURCE_DIR) / "testdata";

  if (!fs::exists(nrl_parm) || !fs::exists(dtm_coeff) || !fs::exists(hwm_data_dir)) {
    spdlog::error("adapter test data paths missing");
    return 10;
  }

  const core::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0, .status = core::Status::Ok};

  core::StateVector state{};
  state.epoch.utc_seconds = 1.0e9;
  state.frame = core::Frame::ECEF;
  state.position_m = core::Vec3{6778137.0, 0.0, 0.0};
  state.velocity_mps = core::Vec3{0.0, 7670.0, 0.0};

  const auto nrl = adapters::Nrlmsis21AtmosphereAdapter::Create({.parm_file = nrl_parm});
  const auto nrl_out = nrl->evaluate(state, wx);
  if (nrl_out.status != core::Status::Ok || !(nrl_out.density_kg_m3 > 0.0) || !std::isfinite(nrl_out.temperature_k)) {
    spdlog::error("nrlmsis adapter evaluation failed");
    return 1;
  }

  core::WeatherIndices wx_hist = wx;
  wx_hist.ap = 18.0;
  wx_hist.ap_3h_current = 40.0;
  wx_hist.kp_3h_current = 5.0;
  wx_hist.ap_msis_history = {40.0, 36.0, 32.0, 28.0, 22.0, 16.0, 10.0};
  wx_hist.has_ap_msis_history = true;
  const auto nrl_hist_out = nrl->evaluate(state, wx_hist);
  if (nrl_hist_out.status != core::Status::Ok) {
    spdlog::error("nrlmsis adapter history evaluation failed");
    return 11;
  }

  msis21::Options msis_options{};
  auto msis_model = msis21::Model::load_from_file(nrl_parm, msis_options);
  const auto geo = core::spherical_geodetic_from_ecef(state.position_m);
  const auto iyd_sec = core::utc_seconds_to_iyd_sec(state.epoch.utc_seconds);
  msis21::Input msis_in{};
  msis_in.iyd = iyd_sec.first;
  msis_in.sec = iyd_sec.second;
  msis_in.alt_km = geo.alt_m * 1.0e-3;
  msis_in.glat_deg = geo.lat_deg;
  msis_in.glon_deg = geo.lon_deg;
  msis_in.stl_hr = core::local_solar_time_hours(state.epoch.utc_seconds, geo.lon_deg);
  msis_in.f107a = wx_hist.f107a;
  msis_in.f107 = wx_hist.f107;
  msis_in.ap = wx_hist.ap_3h_current;
  msis_in.ap_history = wx_hist.ap_msis_history;
  msis_in.has_ap_history = true;
  const auto msis_direct_hist = msis_model.evaluate(msis_in);
  if (msis_direct_hist.status != msis21::Status::Ok) {
    spdlog::error("direct nrlmsis history evaluate failed");
    return 12;
  }
  if (!approx_rel(nrl_hist_out.density_kg_m3, msis_direct_hist.out.rho * 1000.0) ||
      !approx_rel(nrl_hist_out.temperature_k, msis_direct_hist.out.t)) {
    spdlog::error("nrlmsis adapter parity mismatch for history mode");
    return 13;
  }

  msis_in.has_ap_history = false;
  const auto msis_direct_scalar = msis_model.evaluate(msis_in);
  if (msis_direct_scalar.status != msis21::Status::Ok) {
    spdlog::error("direct nrlmsis scalar evaluate failed");
    return 14;
  }

  const auto dtm = adapters::Dtm2020AtmosphereAdapter::Create({.coeff_file = dtm_coeff});
  const auto dtm_out = dtm->evaluate(state, wx);
  if (dtm_out.status != core::Status::Ok || !(dtm_out.density_kg_m3 > 0.0) || !std::isfinite(dtm_out.temperature_k)) {
    spdlog::error("dtm2020 adapter evaluation failed");
    return 2;
  }

  const auto hwm = adapters::Hwm14WindAdapter::Create({.data_dir = hwm_data_dir});
  const auto hwm_out = hwm->evaluate(state, wx);
  if (hwm_out.status != core::Status::Ok || hwm_out.frame != core::Frame::ECEF || !finite(hwm_out.velocity_mps)) {
    spdlog::error("hwm14 adapter evaluation failed");
    return 3;
  }

  weather::StaticSpaceWeatherProvider weather(wx);
  sc::SpacecraftProperties sc{.mass_kg = 600.0, .reference_area_m2 = 4.0, .cd = 2.25, .use_surface_model = false};
  drag::DragAccelerationModel drag_model(weather, *nrl, *hwm);
  const auto drag_out = drag_model.evaluate(state, sc);
  if (drag_out.status != core::Status::Ok || !(drag_out.density_kg_m3 > 0.0) || !finite(drag_out.acceleration_mps2)) {
    spdlog::error("adapter-backed drag pipeline failed");
    return 4;
  }

  return 0;
}

/**
 * @file drag_perf_benchmark.cpp
 * @brief Drag performance benchmark for single-point and batch throughput.
 * @author Watosn
 */

#include <chrono>
#include <cstdlib>
#include <vector>

#include <spdlog/spdlog.h>

#include "astroforces/forces/surface/drag/drag_model.hpp"
#include "astroforces/models/exponential_atmosphere.hpp"
#include "astroforces/sc/spacecraft.hpp"
#include "astroforces/weather/static_provider.hpp"

namespace {

int read_env_int(const char* name, int fallback) {
  const char* raw = std::getenv(name);
  if (!raw || *raw == '\0') {
    return fallback;
  }
  return std::atoi(raw);
}

}  // namespace

int main() {
  using namespace astroforces;
  using clock = std::chrono::steady_clock;

  const int samples = read_env_int("ASTRO_FORCES_PERF_SAMPLES", 40);
  const int iters = read_env_int("ASTRO_FORCES_PERF_ITERS", 5000);
  if (samples <= 0 || iters <= 0) {
    spdlog::error("invalid perf env configuration");
    return 1;
  }

  const core::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0, .status = core::Status::Ok};
  weather::StaticSpaceWeatherProvider weather(wx);
  models::ExponentialAtmosphereModel atmosphere(3.0e-11, 400e3, 65e3, 900.0);
  models::ZeroWindModel wind;
  drag::DragAccelerationModel model(weather, atmosphere, wind);
  sc::SpacecraftProperties sc{.mass_kg = 600.0, .reference_area_m2 = 4.0, .cd = 2.25, .use_surface_model = false};

  std::vector<core::StateVector> states;
  states.reserve(64);
  for (int i = 0; i < 64; ++i) {
    core::StateVector s{};
    s.epoch.utc_seconds = 1.0e9 + static_cast<double>(i) * 60.0;
    s.frame = core::Frame::ECEF;
    s.position_m = core::Vec3{6778137.0 + static_cast<double>(i), 10.0 * static_cast<double>(i), -5.0 * static_cast<double>(i)};
    s.velocity_mps = core::Vec3{20.0, 7670.0, -2.0};
    states.push_back(s);
  }

  double sink = 0.0;

  auto run_single = [&]() {
    auto t0 = clock::now();
    for (int i = 0; i < iters; ++i) {
      const auto r = model.evaluate(states[static_cast<std::size_t>(i % states.size())], sc);
      sink += r.acceleration_mps2.x;
    }
    auto t1 = clock::now();
    return std::chrono::duration<double, std::micro>(t1 - t0).count() / static_cast<double>(iters);
  };

  auto run_batch = [&]() {
    auto t0 = clock::now();
    for (int i = 0; i < iters; ++i) {
      for (const auto& s : states) {
        const auto r = model.evaluate(s, sc);
        sink += r.acceleration_mps2.y;
      }
    }
    auto t1 = clock::now();
    const double evals = static_cast<double>(iters) * static_cast<double>(states.size());
    return std::chrono::duration<double, std::micro>(t1 - t0).count() / evals;
  };

  double single_sum = 0.0;
  double batch_sum = 0.0;
  for (int i = 0; i < samples; ++i) {
    single_sum += run_single();
    batch_sum += run_batch();
  }

  const double single_mean_us = single_sum / static_cast<double>(samples);
  const double batch_mean_us = batch_sum / static_cast<double>(samples);
  const double single_hz = 1.0e6 / single_mean_us;
  const double batch_hz = 1.0e6 / batch_mean_us;

  spdlog::info("drag_perf_benchmark");
  spdlog::info("samples={} iters={} sink={}", samples, iters, sink);
  spdlog::info("single_eval_mean_us={} single_eval_hz={}", single_mean_us, single_hz);
  spdlog::info("batch_eval_mean_us={} batch_eval_hz={}", batch_mean_us, batch_hz);
  return 0;
}

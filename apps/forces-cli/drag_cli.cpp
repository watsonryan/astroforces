/**
 * @file drag_cli.cpp
 * @brief astroforces drag command-line entrypoint.
 * @author Watson
 */

#include <cstdlib>
#include <memory>
#include <string>

#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/adapters/dtm2020_adapter.hpp"
#include "astroforces/adapters/hwm14_adapter.hpp"
#include "astroforces/adapters/nrlmsis21_adapter.hpp"
#include "astroforces/forces/surface/drag/drag_model.hpp"
#include "astroforces/models/exponential_atmosphere.hpp"
#include "astroforces/sc/spacecraft.hpp"
#include "astroforces/weather/celestrak_csv_provider.hpp"
#include "astroforces/weather/static_provider.hpp"

int main(int argc, char** argv) {
  double x_m{};
  double y_m{};
  double z_m{};
  double vx_mps{};
  double vy_mps{};
  double vz_mps{};
  double epoch_utc_s{};
  std::string model_name{"basic"};
  std::string model_data{};
  std::string wind_name{"zero"};
  std::string wind_data{};
  std::string weather_csv{};

  CLI::App app{"astroforces drag command-line entrypoint"};
  app.add_option("x_m", x_m)->required();
  app.add_option("y_m", y_m)->required();
  app.add_option("z_m", z_m)->required();
  app.add_option("vx_mps", vx_mps)->required();
  app.add_option("vy_mps", vy_mps)->required();
  app.add_option("vz_mps", vz_mps)->required();
  app.add_option("epoch_utc_s", epoch_utc_s)->required();
  app.add_option("model", model_name)->check(CLI::IsMember({"basic", "nrlmsis", "dtm2020"}))->capture_default_str();
  app.add_option("model_data", model_data)->capture_default_str();
  app.add_option("wind", wind_name)->check(CLI::IsMember({"zero", "hwm14"}))->capture_default_str();
  app.add_option("wind_data", wind_data)->capture_default_str();
  app.add_option("weather_csv", weather_csv)->capture_default_str();
  CLI11_PARSE(app, argc, argv);

  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{x_m, y_m, z_m};
  state.velocity_mps = astroforces::core::Vec3{vx_mps, vy_mps, vz_mps};
  state.epoch.utc_seconds = epoch_utc_s;
  state.frame = astroforces::core::Frame::ECEF;

  std::unique_ptr<astroforces::core::ISpaceWeatherProvider> weather{};
  if (!weather_csv.empty()) {
    weather = astroforces::weather::CelesTrakCsvSpaceWeatherProvider::Create(
        astroforces::weather::CelesTrakCsvSpaceWeatherProvider::Config{.csv_file = weather_csv});
  } else {
    const astroforces::core::WeatherIndices wx{.f107 = 150.0, .f107a = 150.0, .ap = 4.0, .kp = 2.0,
                                            .status = astroforces::core::Status::Ok};
    weather = std::make_unique<astroforces::weather::StaticSpaceWeatherProvider>(wx);
  }

  std::unique_ptr<astroforces::core::IAtmosphereModel> atmosphere{};
  if (model_name == "nrlmsis") {
    atmosphere = astroforces::adapters::Nrlmsis21AtmosphereAdapter::Create(
        astroforces::adapters::Nrlmsis21AtmosphereAdapter::Config{.parm_file = model_data});
  } else if (model_name == "dtm2020") {
    atmosphere = astroforces::adapters::Dtm2020AtmosphereAdapter::Create(
        astroforces::adapters::Dtm2020AtmosphereAdapter::Config{.coeff_file = model_data});
  } else {
    atmosphere = std::make_unique<astroforces::models::ExponentialAtmosphereModel>(1.225, 0.0, 7000.0, 1000.0);
  }
  std::unique_ptr<astroforces::core::IWindModel> wind{};
  if (wind_name == "hwm14") {
    wind = astroforces::adapters::Hwm14WindAdapter::Create(
        astroforces::adapters::Hwm14WindAdapter::Config{.data_dir = wind_data});
  } else {
    wind = std::make_unique<astroforces::models::ZeroWindModel>();
  }

  astroforces::sc::SpacecraftProperties sc{.mass_kg = 1200.0,
                                       .reference_area_m2 = 12.0,
                                       .cd = 2.2,
                                       .use_surface_model = false,
                                       .surfaces = {}};

  astroforces::forces::DragAccelerationModel model(*weather, *atmosphere, *wind);
  const auto result = model.evaluate(state, sc);

  if (result.status != astroforces::core::Status::Ok) {
    spdlog::error("drag eval failed");
    return 2;
  }

  fmt::print("ax={} ay={} az={}\n", result.acceleration_mps2.x, result.acceleration_mps2.y, result.acceleration_mps2.z);
  fmt::print("rho={} temp_k={} vrel_mps={} q_pa={} area={} cd_eff={}\n", result.density_kg_m3, result.temperature_k,
             result.relative_speed_mps, result.dynamic_pressure_pa, result.area_m2, result.cd);
  fmt::print("wx_source={} wx_interp={} wx_extrap={} f107={} f107a={} ap={} kp={} ap3h={} kp3h={}\n",
             static_cast<int>(result.weather.source), result.weather.interpolated ? 1 : 0, result.weather.extrapolated ? 1 : 0,
             result.weather.f107, result.weather.f107a, result.weather.ap, result.weather.kp, result.weather.ap_3h_current,
             result.weather.kp_3h_current);
  return 0;
}

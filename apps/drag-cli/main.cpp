/**
 * @file main.cpp
 * @brief astrodynamics-forces-cpp drag command-line entrypoint.
 * @author Watosn
 */

#include <cstdlib>
#include <memory>
#include <string>

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
  if (argc < 8 || argc > 13) {
    spdlog::error("usage: drag_cli <x_m> <y_m> <z_m> <vx_mps> <vy_mps> <vz_mps> <epoch_utc_s> [model] [model_data] [wind] [wind_data] [weather_csv]");
    spdlog::error("models: basic | nrlmsis | dtm2020");
    spdlog::error("wind: zero | hwm14");
    spdlog::error("weather_csv: CelesTrak SW-Last5Years.csv (optional)");
    return 1;
  }

  astroforces::core::StateVector state{};
  state.position_m = astroforces::core::Vec3{std::atof(argv[1]), std::atof(argv[2]), std::atof(argv[3])};
  state.velocity_mps = astroforces::core::Vec3{std::atof(argv[4]), std::atof(argv[5]), std::atof(argv[6])};
  state.epoch.utc_seconds = std::atof(argv[7]);
  state.frame = astroforces::core::Frame::ECEF;

  const std::string weather_csv = (argc >= 13) ? argv[12] : "";
  const std::string model_name = (argc >= 9) ? argv[8] : "basic";
  const std::string model_data = (argc >= 10) ? argv[9] : "";
  const std::string wind_name = (argc >= 11) ? argv[10] : "zero";
const std::string wind_data = (argc >= 12) ? argv[11] : "";

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

  astroforces::drag::DragAccelerationModel model(*weather, *atmosphere, *wind);
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

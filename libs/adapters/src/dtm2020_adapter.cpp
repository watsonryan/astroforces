/**
 * @file dtm2020_adapter.cpp
 * @brief DTM2020 adapter implementation.
 * @author Watosn
 */

#include "astroforces/adapters/dtm2020_adapter.hpp"

#include <utility>

#include "astroforces/atmo/conversions.hpp"
#include "dtm2020/dtm2020_operational.hpp"

namespace astroforces::adapters {

class Dtm2020AtmosphereAdapter::Impl {
 public:
  explicit Impl(dtm2020::Dtm2020Operational model) : model_(std::move(model)) {}

  dtm2020::Dtm2020Operational model_;
};

std::unique_ptr<Dtm2020AtmosphereAdapter> Dtm2020AtmosphereAdapter::Create(const Config& config) {
  auto ptr = std::unique_ptr<Dtm2020AtmosphereAdapter>(new Dtm2020AtmosphereAdapter(config));
  auto loaded = dtm2020::Dtm2020Operational::LoadFromFile(config.coeff_file);
  if (!loaded.has_value()) {
    return ptr;
  }
  ptr->impl_ = std::make_shared<Impl>(loaded.value());
  return ptr;
}

astroforces::core::AtmosphereSample Dtm2020AtmosphereAdapter::evaluate(const astroforces::core::StateVector& state,
                                                                    const astroforces::core::WeatherIndices& weather) const {
  if (!impl_) {
    return astroforces::core::AtmosphereSample{.status = astroforces::core::Status::DataUnavailable};
  }
  if (weather.status != astroforces::core::Status::Ok) {
    return astroforces::core::AtmosphereSample{.status = weather.status};
  }

  const auto geo = astroforces::core::spherical_geodetic_from_ecef(state.position_m);
  const auto iyd_sec = astroforces::core::utc_seconds_to_iyd_sec(state.epoch.utc_seconds);
  const int doy = iyd_sec.first % 1000;

  dtm2020::OperationalInputs in{};
  in.altitude_km = geo.alt_m * 1e-3;
  in.latitude_deg = geo.lat_deg;
  in.longitude_deg = geo.lon_deg;
  in.local_time_h = astroforces::core::local_solar_time_hours(state.epoch.utc_seconds, geo.lon_deg);
  in.day_of_year = static_cast<double>(doy);
  in.f107 = weather.f107;
  in.f107m = weather.f107a;
  in.kp_delayed_3h = weather.kp_3h_current;
  in.kp_mean_24h = weather.kp;

  const auto out = impl_->model_.Evaluate(in);
  if (!out.has_value()) {
    return astroforces::core::AtmosphereSample{.status = astroforces::core::Status::NumericalError};
  }

  return astroforces::core::AtmosphereSample{
      .density_kg_m3 = out.value().density_g_cm3 * 1000.0,
      .temperature_k = out.value().temperature_k,
      .status = astroforces::core::Status::Ok};
}

}  // namespace astroforces::adapters

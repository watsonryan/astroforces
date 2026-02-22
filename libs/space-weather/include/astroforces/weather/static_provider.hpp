/**
 * @file static_provider.hpp
 * @brief Constant-value weather provider.
 * @author Watosn
 */
#pragma once

#include "astroforces/atmo/interfaces.hpp"

namespace astroforces::weather {

class StaticSpaceWeatherProvider final : public astroforces::core::ISpaceWeatherProvider {
 public:
  explicit StaticSpaceWeatherProvider(astroforces::core::WeatherIndices indices) : indices_(indices) {
    indices_.source = astroforces::core::WeatherSource::StaticProvider;
    indices_.ap_3h_current = indices_.ap;
    indices_.kp_3h_current = indices_.kp;
    indices_.ap_3h_utc.fill(indices_.ap);
    indices_.kp_3h_utc.fill(indices_.kp);
    indices_.ap_msis_history.fill(indices_.ap);
    indices_.has_ap_msis_history = true;
  }
  [[nodiscard]] astroforces::core::WeatherIndices at(const astroforces::core::Epoch& /*epoch*/) const override {
    return indices_;
  }

 private:
  astroforces::core::WeatherIndices indices_{};
};

}  // namespace astroforces::weather

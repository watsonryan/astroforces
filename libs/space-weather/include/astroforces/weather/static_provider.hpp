/**
 * @file static_provider.hpp
 * @brief Constant-value weather provider.
 * @author Watosn
 */
#pragma once

#include "astroforces/core/interfaces.hpp"

namespace astroforces::weather {

/**
 * @brief Constant weather provider for testing and deterministic runs.
 */
class StaticSpaceWeatherProvider final : public astroforces::core::ISpaceWeatherProvider {
 public:
  /**
   * @brief Construct provider from fixed weather indices.
   * @note History arrays and source metadata are populated automatically.
   */
  explicit StaticSpaceWeatherProvider(astroforces::core::WeatherIndices indices) : indices_(indices) {
    indices_.source = astroforces::core::WeatherSource::StaticProvider;
    indices_.ap_3h_current = indices_.ap;
    indices_.kp_3h_current = indices_.kp;
    indices_.ap_3h_utc.fill(indices_.ap);
    indices_.kp_3h_utc.fill(indices_.kp);
    indices_.ap_msis_history.fill(indices_.ap);
    indices_.has_ap_msis_history = true;
  }
  /**
   * @brief Return fixed weather values for any epoch.
   */
  [[nodiscard]] astroforces::core::WeatherIndices at(const astroforces::core::Epoch& /*epoch*/) const override {
    return indices_;
  }

 private:
  astroforces::core::WeatherIndices indices_{};
};

}  // namespace astroforces::weather

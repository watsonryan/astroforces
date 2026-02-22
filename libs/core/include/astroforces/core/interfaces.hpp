/**
 * @file interfaces.hpp
 * @brief Core model interfaces for atmosphere/wind/weather.
 * @author Watosn
 */
#pragma once

#include "astroforces/core/types.hpp"

namespace astroforces::core {

/**
 * @brief Interface for time-indexed space weather providers.
 */
class ISpaceWeatherProvider {
 public:
  virtual ~ISpaceWeatherProvider() = default;
  /**
   * @brief Return space weather indices at a given epoch.
   * @param epoch UTC epoch.
   * @return Weather indices with `status` set.
   */
  [[nodiscard]] virtual WeatherIndices at(const Epoch& epoch) const = 0;
};

/**
 * @brief Interface for thermospheric density/temperature models.
 */
class IAtmosphereModel {
 public:
  virtual ~IAtmosphereModel() = default;
  /**
   * @brief Evaluate atmospheric state for a spacecraft state and weather forcing.
   * @param state Input state vector.
   * @param weather Space weather forcing.
   * @return Atmosphere sample with `status` set.
   */
  [[nodiscard]] virtual AtmosphereSample evaluate(const StateVector& state,
                                                  const WeatherIndices& weather) const = 0;
};

/**
 * @brief Interface for neutral wind models.
 */
class IWindModel {
 public:
  virtual ~IWindModel() = default;
  /**
   * @brief Evaluate neutral wind for a spacecraft state and weather forcing.
   * @param state Input state vector.
   * @param weather Space weather forcing.
   * @return Wind sample with `status` set.
   */
  [[nodiscard]] virtual WindSample evaluate(const StateVector& state,
                                            const WeatherIndices& weather) const = 0;
};

}  // namespace astroforces::core

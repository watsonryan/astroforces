/**
 * @file nrlmsis21_adapter.hpp
 * @brief Adapter between astroforces atmosphere interface and NRLMSIS 2.1 model.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>

#include "astroforces/core/interfaces.hpp"

namespace msis21 {
class Model;
}

namespace astroforces::adapters {

/**
 * @brief Atmosphere adapter backed by NRLMSIS 2.1.
 */
class Nrlmsis21AtmosphereAdapter final : public astroforces::core::IAtmosphereModel {
 public:
  /**
   * @brief Configuration for NRLMSIS adapter construction.
   */
  struct Config {
    std::filesystem::path parm_file{};
  };

  /**
   * @brief Factory helper that loads model resources.
   */
  static std::unique_ptr<Nrlmsis21AtmosphereAdapter> Create(const Config& config);

  /**
   * @brief Evaluate density/temperature using NRLMSIS 2.1.
   */
  [[nodiscard]] astroforces::core::AtmosphereSample evaluate(const astroforces::core::StateVector& state,
                                                          const astroforces::core::WeatherIndices& weather) const override;

 private:
  class Impl;
  explicit Nrlmsis21AtmosphereAdapter(Config config) : config_(std::move(config)) {}

  Config config_{};
  std::shared_ptr<Impl> impl_{};
};

}  // namespace astroforces::adapters

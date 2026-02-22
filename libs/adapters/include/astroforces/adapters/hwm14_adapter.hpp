/**
 * @file hwm14_adapter.hpp
 * @brief Adapter between astroforces wind interface and HWM14 model.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>

#include "astroforces/core/interfaces.hpp"

namespace hwm14 {
class Model;
}

namespace astroforces::adapters {

/**
 * @brief Wind adapter backed by HWM14.
 */
class Hwm14WindAdapter final : public astroforces::core::IWindModel {
 public:
  /**
   * @brief Configuration for HWM14 adapter construction.
   */
  struct Config {
    std::filesystem::path data_dir{};
  };

  /**
   * @brief Factory helper that loads model resources.
   */
  static std::unique_ptr<Hwm14WindAdapter> Create(const Config& config);

  /**
   * @brief Evaluate neutral wind using HWM14.
   */
  [[nodiscard]] astroforces::core::WindSample evaluate(const astroforces::core::StateVector& state,
                                                    const astroforces::core::WeatherIndices& weather) const override;

 private:
  class Impl;
  explicit Hwm14WindAdapter(Config config) : config_(std::move(config)) {}

  Config config_{};
  std::shared_ptr<Impl> impl_{};
};

}  // namespace astroforces::adapters

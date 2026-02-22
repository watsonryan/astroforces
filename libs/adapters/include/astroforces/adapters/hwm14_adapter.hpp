/**
 * @file hwm14_adapter.hpp
 * @brief Adapter between dragcpp wind interface and HWM14 model.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>

#include "astroforces/atmo/interfaces.hpp"

namespace hwm14 {
class Model;
}

namespace astroforces::adapters {

class Hwm14WindAdapter final : public astroforces::core::IWindModel {
 public:
  struct Config {
    std::filesystem::path data_dir{};
  };

  static std::unique_ptr<Hwm14WindAdapter> Create(const Config& config);

  [[nodiscard]] astroforces::core::WindSample evaluate(const astroforces::core::StateVector& state,
                                                    const astroforces::core::WeatherIndices& weather) const override;

 private:
  class Impl;
  explicit Hwm14WindAdapter(Config config) : config_(std::move(config)) {}

  Config config_{};
  std::shared_ptr<Impl> impl_{};
};

}  // namespace astroforces::adapters

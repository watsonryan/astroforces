/**
 * @file dtm2020_adapter.hpp
 * @brief Adapter between dragcpp atmosphere interface and DTM2020 model.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>

#include "astroforces/atmo/interfaces.hpp"

namespace dtm2020 {
class Dtm2020Operational;
}

namespace astroforces::adapters {

class Dtm2020AtmosphereAdapter final : public astroforces::core::IAtmosphereModel {
 public:
  struct Config {
    std::filesystem::path coeff_file{};
  };

  static std::unique_ptr<Dtm2020AtmosphereAdapter> Create(const Config& config);

  [[nodiscard]] astroforces::core::AtmosphereSample evaluate(const astroforces::core::StateVector& state,
                                                          const astroforces::core::WeatherIndices& weather) const override;

 private:
  class Impl;
  explicit Dtm2020AtmosphereAdapter(Config config) : config_(std::move(config)) {}

  Config config_{};
  std::shared_ptr<Impl> impl_{};
};

}  // namespace astroforces::adapters

/**
 * @file celestrak_csv_provider.hpp
 * @brief Space weather provider backed by CelesTrak SW-Last5Years CSV.
 * @author Watosn
 */
#pragma once

#include <array>
#include <filesystem>
#include <memory>
#include <vector>

#include "astroforces/core/interfaces.hpp"

namespace astroforces::weather {

/**
 * @brief Space weather provider backed by CelesTrak daily CSV data.
 */
class CelesTrakCsvSpaceWeatherProvider final : public astroforces::core::ISpaceWeatherProvider {
 public:
  /**
   * @brief Normalized in-memory daily sample.
   */
  struct DailySample {
    double day_start_utc_s{};
    double f107_obs{};
    double f107_obs_center81{};
    double ap_avg{};
    double kp_avg{};
    std::array<double, 8> ap_3h_utc{};
    std::array<double, 8> kp_3h_utc{};
    bool f107_observed{true};
    bool geomagnetic_observed{true};
  };

  /**
   * @brief CSV provider configuration.
   */
  struct Config {
    std::filesystem::path csv_file{};
  };

  /**
   * @brief Factory helper that parses and validates CSV input.
   */
  static std::unique_ptr<CelesTrakCsvSpaceWeatherProvider> Create(const Config& config);

  /**
   * @brief Evaluate weather indices at an epoch via interpolation/selection.
   */
  [[nodiscard]] astroforces::core::WeatherIndices at(const astroforces::core::Epoch& epoch) const override;

 private:
  explicit CelesTrakCsvSpaceWeatherProvider(std::vector<DailySample> samples) : samples_(std::move(samples)) {}

  std::vector<DailySample> samples_{};
};

}  // namespace astroforces::weather

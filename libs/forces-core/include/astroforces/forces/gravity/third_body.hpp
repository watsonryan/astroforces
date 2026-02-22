/**
 * @file third_body.hpp
 * @brief Third-body perturbation model using JPL ephemerides.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>
#include <string>

#include "astroforces/forces/core/perturbation.hpp"

namespace jpl::eph {
class Ephemeris;
class Workspace;
}

namespace astroforces::forces {

/**
 * @brief Sun/Moon third-body perturbation model.
 */
class ThirdBodyPerturbationModel final : public IPerturbationModel {
 public:
  /**
   * @brief Configuration for third-body model construction.
   */
  struct Config {
    std::filesystem::path ephemeris_file{};
    bool use_sun{true};
    bool use_moon{true};
    bool use_goce_eq79_indirect_j2{true};
    double mu_sun_m3_s2{1.32712440018e20};
    double mu_moon_m3_s2{4.9048695e12};
    std::string name{"third_body"};
  };

  /**
   * @brief Factory helper that loads ephemeris resources.
   */
  static std::unique_ptr<ThirdBodyPerturbationModel> Create(const Config& config);

  /**
   * @brief Evaluate third-body perturbation contribution.
   */
  [[nodiscard]] PerturbationContribution evaluate(const PerturbationRequest& request) const override;

 private:
  explicit ThirdBodyPerturbationModel(Config config);

  Config config_{};
  std::shared_ptr<jpl::eph::Ephemeris> ephemeris_{};
};

}  // namespace astroforces::forces

/**
 * @file third_body.hpp
 * @brief Third-body perturbation model using JPL ephemerides.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>
#include <string>

#include "astroforces/forces/perturbation.hpp"

namespace jpl::eph {
class Ephemeris;
class Workspace;
}

namespace astroforces::forces {

class ThirdBodyPerturbationModel final : public IPerturbationModel {
 public:
  struct Config {
    std::filesystem::path ephemeris_file{};
    bool use_sun{true};
    bool use_moon{true};
    bool use_goce_eq79_indirect_j2{true};
    double mu_sun_m3_s2{1.32712440018e20};
    double mu_moon_m3_s2{4.9048695e12};
    std::string name{"third_body"};
  };

  static std::unique_ptr<ThirdBodyPerturbationModel> Create(const Config& config);

  [[nodiscard]] PerturbationContribution evaluate(const PerturbationRequest& request) const override;

 private:
  explicit ThirdBodyPerturbationModel(Config config);

  Config config_{};
  std::shared_ptr<jpl::eph::Ephemeris> ephemeris_{};
  mutable std::shared_ptr<jpl::eph::Workspace> workspace_{};
};

}  // namespace astroforces::forces

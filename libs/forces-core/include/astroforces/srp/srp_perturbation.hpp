/**
 * @file srp_perturbation.hpp
 * @brief SRP perturbation wrapper for generic force stack.
 * @author Watosn
 */
#pragma once

#include <memory>
#include <string>

#include "astroforces/forces/perturbation.hpp"
#include "astroforces/srp/srp_model.hpp"

namespace astroforces::srp {

class SrpPerturbationModel final : public astroforces::forces::IPerturbationModel {
 public:
  SrpPerturbationModel(std::unique_ptr<SrpAccelerationModel> srp,
                       const astroforces::sc::SpacecraftProperties* default_spacecraft = nullptr,
                       std::string name = "srp")
      : srp_(std::move(srp)), default_spacecraft_(default_spacecraft), name_(std::move(name)) {}

  [[nodiscard]] astroforces::forces::PerturbationContribution evaluate(
      const astroforces::forces::PerturbationRequest& request) const override;

 private:
  std::unique_ptr<SrpAccelerationModel> srp_{};
  const astroforces::sc::SpacecraftProperties* default_spacecraft_{nullptr};
  std::string name_{};
};

}  // namespace astroforces::srp


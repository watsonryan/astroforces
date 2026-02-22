/**
 * @file relativity_perturbation.hpp
 * @brief Relativity perturbation wrapper for generic force stack.
 * @author Watosn
 */
#pragma once

#include <memory>
#include <string>

#include "astroforces/forces/core/perturbation.hpp"
#include "astroforces/forces/gravity/relativity_model.hpp"

namespace astroforces::forces {

/**
 * @brief Perturbation-stack wrapper around RelativityAccelerationModel.
 */
class RelativityPerturbationModel final : public IPerturbationModel {
 public:
  /**
   * @brief Construct wrapper with owned relativity model.
   */
  explicit RelativityPerturbationModel(std::unique_ptr<RelativityAccelerationModel> relativity,
                                       std::string name = "relativity")
      : relativity_(std::move(relativity)), name_(std::move(name)) {}

  /**
   * @brief Evaluate relativity contribution.
   */
  [[nodiscard]] PerturbationContribution evaluate(const PerturbationRequest& request) const override;

 private:
  std::unique_ptr<RelativityAccelerationModel> relativity_{};
  std::string name_{};
};

}  // namespace astroforces::forces

/**
 * @file erp_perturbation.hpp
 * @brief ERP perturbation wrapper for generic force stack.
 * @author Watosn
 */
#pragma once

#include <string>

#include "astroforces/forces/surface/erp/erp_model.hpp"
#include "astroforces/forces/core/perturbation.hpp"

namespace astroforces::erp {

class ErpPerturbationModel final : public astroforces::forces::IPerturbationModel {
 public:
  ErpPerturbationModel(ErpAccelerationModel erp = ErpAccelerationModel{},
                       const astroforces::sc::SpacecraftProperties* default_spacecraft = nullptr,
                       std::string name = "erp")
      : erp_(erp), default_spacecraft_(default_spacecraft), name_(std::move(name)) {}

  [[nodiscard]] astroforces::forces::PerturbationContribution evaluate(
      const astroforces::forces::PerturbationRequest& request) const override;

 private:
  ErpAccelerationModel erp_{};
  const astroforces::sc::SpacecraftProperties* default_spacecraft_{nullptr};
  std::string name_{};
};

}  // namespace astroforces::erp


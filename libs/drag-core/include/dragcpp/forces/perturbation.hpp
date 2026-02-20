/**
 * @file perturbation.hpp
 * @brief Generic perturbation model interfaces and aggregation utilities.
 * @author Watosn
 */
#pragma once

#include <memory>
#include <string>
#include <vector>

#include "dragcpp/atmo/types.hpp"
#include "dragcpp/sc/spacecraft.hpp"

namespace dragcpp::forces {

enum class PerturbationType : unsigned char { Unknown, Drag, Gravity, SRP, ThirdBody };

struct PerturbationRequest {
  dragcpp::atmo::StateVector state{};
  const dragcpp::sc::SpacecraftProperties* spacecraft{nullptr};
};

struct PerturbationContribution {
  std::string name{};
  PerturbationType type{PerturbationType::Unknown};
  dragcpp::atmo::Vec3 acceleration_mps2{};
  dragcpp::atmo::Status status{dragcpp::atmo::Status::Ok};
};

struct PerturbationResult {
  dragcpp::atmo::Vec3 total_acceleration_mps2{};
  std::vector<PerturbationContribution> contributions{};
  dragcpp::atmo::Status status{dragcpp::atmo::Status::Ok};
};

class IPerturbationModel {
 public:
  virtual ~IPerturbationModel() = default;
  [[nodiscard]] virtual PerturbationContribution evaluate(const PerturbationRequest& request) const = 0;
};

class PerturbationStack final {
 public:
  void add(std::unique_ptr<IPerturbationModel> model);
  [[nodiscard]] PerturbationResult evaluate(const PerturbationRequest& request) const;

 private:
  std::vector<std::unique_ptr<IPerturbationModel>> models_{};
};

}  // namespace dragcpp::forces

/**
 * @file constituent_tide.hpp
 * @brief Constituent ocean/atmospheric tide SH correction model.
 * @author Watosn
 */
#pragma once

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

namespace astroforces::forces::tides {

class ConstituentTideModel final {
 public:
  static std::unique_ptr<ConstituentTideModel> load_from_file(const std::filesystem::path& path, int max_degree);

  [[nodiscard]] bool empty() const noexcept;
  [[nodiscard]] int max_degree() const noexcept;

  void add_delta_coefficients(double jd_utc, double gmst_rad, Eigen::MatrixXd& dC, Eigen::MatrixXd& dS) const;

 private:
  struct Wave {
    std::string name{};
    std::array<double, 6> doodson{};
    Eigen::MatrixXd C1{};
    Eigen::MatrixXd C2{};
    Eigen::MatrixXd S1{};
    Eigen::MatrixXd S2{};
  };

  explicit ConstituentTideModel(int max_degree) : max_degree_(max_degree) {}

  int max_degree_{0};
  std::vector<Wave> waves_{};
};

}  // namespace astroforces::forces::tides

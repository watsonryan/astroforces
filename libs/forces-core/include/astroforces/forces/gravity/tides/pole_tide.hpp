/**
 * @file pole_tide.hpp
 * @brief Pole tide correction helpers for gravity SH models.
 * @author Watosn
 */
#pragma once

#include <filesystem>
#include <memory>

#include <Eigen/Dense>

namespace astroforces::forces::tides {

class OceanPoleTideModel final {
 public:
  static std::unique_ptr<OceanPoleTideModel> load_from_file(const std::filesystem::path& path, int max_degree);
  [[nodiscard]] bool empty() const noexcept;
  void add_delta_coefficients(double m1_arcsec, double m2_arcsec, Eigen::MatrixXd& dC, Eigen::MatrixXd& dS) const;

 private:
  explicit OceanPoleTideModel(int max_degree) : max_degree_(max_degree) {}
  int max_degree_{0};
  Eigen::MatrixXd cnmp_{};
  Eigen::MatrixXd cnmm_{};
  Eigen::MatrixXd snmp_{};
  Eigen::MatrixXd snmm_{};
};

void add_pole_solid_tide_delta(double mjd_tt,
                               double xp_rad,
                               double yp_rad,
                               Eigen::MatrixXd& dC,
                               Eigen::MatrixXd& dS);

void add_pole_ocean_tide_delta(double mjd_tt,
                               double xp_rad,
                               double yp_rad,
                               Eigen::MatrixXd& dC,
                               Eigen::MatrixXd& dS);

}  // namespace astroforces::forces::tides

/**
 * @file orbital_elements_cli.cpp
 * @brief ECI/osculating-orbital-element conversion utility.
 * @author Watson
 */

#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <CLI/CLI.hpp>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include "astroforces/core/orbital_elements.hpp"

namespace {

bool parse_csv_values(const std::string& line, std::vector<double>& values) {
  std::stringstream ss(line);
  std::string tok;
  values.clear();
  while (std::getline(ss, tok, ',')) {
    if (tok.empty()) {
      return false;
    }
    char* end = nullptr;
    const double v = std::strtod(tok.c_str(), &end);
    if (end == tok.c_str() || *end != '\0') {
      return false;
    }
    values.push_back(v);
  }
  return !values.empty();
}

std::vector<std::string> split_csv_tokens(const std::string& line) {
  std::stringstream ss(line);
  std::string tok;
  std::vector<std::string> tokens;
  while (std::getline(ss, tok, ',')) {
    tokens.push_back(tok);
  }
  return tokens;
}

int find_column(const std::vector<std::string>& header, const std::string& name) {
  for (std::size_t i = 0; i < header.size(); ++i) {
    if (header[i] == name) {
      return static_cast<int>(i);
    }
  }
  return -1;
}

bool use_stdio_path(const std::filesystem::path& path) {
  return path.string() == "-";
}

std::istream* open_input_stream(const std::filesystem::path& input_csv, std::unique_ptr<std::ifstream>& owned_stream) {
  if (use_stdio_path(input_csv)) {
    return &std::cin;
  }
  owned_stream = std::make_unique<std::ifstream>(input_csv);
  if (!*owned_stream) {
    return nullptr;
  }
  return owned_stream.get();
}

std::ostream* open_output_stream(const std::filesystem::path& output_csv, std::unique_ptr<std::ofstream>& owned_stream) {
  if (use_stdio_path(output_csv)) {
    return &std::cout;
  }
  owned_stream = std::make_unique<std::ofstream>(output_csv);
  if (!*owned_stream) {
    return nullptr;
  }
  return owned_stream.get();
}

int run_eci_to_osc_batch(const std::filesystem::path& input_csv, const std::filesystem::path& output_csv, const double mu) {
  using namespace astroforces::core;

  std::unique_ptr<std::ifstream> input_file;
  std::istream* in = open_input_stream(input_csv, input_file);
  if (!in) {
    spdlog::error("failed to open input csv: {}", input_csv.string());
    return 10;
  }
  std::unique_ptr<std::ofstream> output_file;
  std::ostream* out = open_output_stream(output_csv, output_file);
  if (!out) {
    spdlog::error("failed to open output csv: {}", output_csv.string());
    return 11;
  }

  *out << "a_m,p_m,e,i_rad,raan_rad,argp_rad,nu_rad,E_rad,M_rad,circular,equatorial,status\n";

  std::string line;
  std::size_t line_no = 0;
  std::vector<double> values;
  int rx_col = 0;
  int ry_col = 1;
  int rz_col = 2;
  int vx_col = 3;
  int vy_col = 4;
  int vz_col = 5;
  while (std::getline(*in, line)) {
    ++line_no;
    if (line.empty()) {
      continue;
    }
    if (!parse_csv_values(line, values)) {
      if (line_no == 1 && line.find("rx_eci_m") != std::string::npos) {
        const auto header = split_csv_tokens(line);
        rx_col = find_column(header, "rx_eci_m");
        ry_col = find_column(header, "ry_eci_m");
        rz_col = find_column(header, "rz_eci_m");
        vx_col = find_column(header, "vx_eci_mps");
        vy_col = find_column(header, "vy_eci_mps");
        vz_col = find_column(header, "vz_eci_mps");
        continue;
      }
      spdlog::warn("skipping malformed row {}", line_no);
      continue;
    }
    const int max_col = std::max({rx_col, ry_col, rz_col, vx_col, vy_col, vz_col});
    if (max_col < 0 || values.size() <= static_cast<std::size_t>(max_col)) {
      if (line_no == 1 && line.find("rx_eci_m") != std::string::npos) {
        continue;
      }
      spdlog::warn("skipping row {} with {} columns, expected ECI state columns", line_no, values.size());
      continue;
    }

    const auto elements = state_eci_to_osculating_orbital_elements(
        Vec3{values[static_cast<std::size_t>(rx_col)],
             values[static_cast<std::size_t>(ry_col)],
             values[static_cast<std::size_t>(rz_col)]},
        Vec3{values[static_cast<std::size_t>(vx_col)],
             values[static_cast<std::size_t>(vy_col)],
             values[static_cast<std::size_t>(vz_col)]},
        mu);

    *out << fmt::format(
        "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{},{},{}\n",
        elements.semi_major_axis_m,
        elements.semi_latus_rectum_m,
        elements.eccentricity,
        elements.inclination_rad,
        elements.raan_rad,
        elements.argument_of_periapsis_rad,
        elements.true_anomaly_rad,
        elements.eccentric_anomaly_rad,
        elements.mean_anomaly_rad,
        elements.circular ? 1 : 0,
        elements.equatorial ? 1 : 0,
        static_cast<int>(elements.status));
  }

  spdlog::info("wrote orbital-elements batch output: {}", use_stdio_path(output_csv) ? "stdout" : output_csv.string());
  return 0;
}

int run_osc_to_eci_batch(const std::filesystem::path& input_csv, const std::filesystem::path& output_csv, const double mu) {
  using namespace astroforces::core;

  std::unique_ptr<std::ifstream> input_file;
  std::istream* in = open_input_stream(input_csv, input_file);
  if (!in) {
    spdlog::error("failed to open input csv: {}", input_csv.string());
    return 12;
  }
  std::unique_ptr<std::ofstream> output_file;
  std::ostream* out = open_output_stream(output_csv, output_file);
  if (!out) {
    spdlog::error("failed to open output csv: {}", output_csv.string());
    return 13;
  }

  *out << "rx_eci_m,ry_eci_m,rz_eci_m,vx_eci_mps,vy_eci_mps,vz_eci_mps,status\n";

  std::string line;
  std::size_t line_no = 0;
  std::vector<double> values;
  int a_col = 0;
  int e_col = 1;
  int i_col = 2;
  int raan_col = 3;
  int argp_col = 4;
  int nu_col = 5;
  int p_col = -1;
  while (std::getline(*in, line)) {
    ++line_no;
    if (line.empty()) {
      continue;
    }
    if (!parse_csv_values(line, values)) {
      if (line_no == 1 && line.find("a_m") != std::string::npos) {
        const auto header = split_csv_tokens(line);
        a_col = find_column(header, "a_m");
        e_col = find_column(header, "e");
        i_col = find_column(header, "i_rad");
        raan_col = find_column(header, "raan_rad");
        argp_col = find_column(header, "argp_rad");
        nu_col = find_column(header, "nu_rad");
        p_col = find_column(header, "p_m");
        continue;
      }
      spdlog::warn("skipping malformed row {}", line_no);
      continue;
    }
    const int max_col = std::max({a_col, e_col, i_col, raan_col, argp_col, nu_col, p_col});
    if (a_col < 0 || e_col < 0 || i_col < 0 || raan_col < 0 || argp_col < 0 || nu_col < 0
        || values.size() <= static_cast<std::size_t>(max_col)) {
      if (line_no == 1 && line.find("a_m") != std::string::npos) {
        continue;
      }
      spdlog::warn("skipping row {} with {} columns, expected osculating-element columns", line_no, values.size());
      continue;
    }

    ClassicalOrbitalElements elements{};
    elements.semi_major_axis_m = values[static_cast<std::size_t>(a_col)];
    elements.eccentricity = values[static_cast<std::size_t>(e_col)];
    elements.inclination_rad = values[static_cast<std::size_t>(i_col)];
    elements.raan_rad = values[static_cast<std::size_t>(raan_col)];
    elements.argument_of_periapsis_rad = values[static_cast<std::size_t>(argp_col)];
    elements.true_anomaly_rad = values[static_cast<std::size_t>(nu_col)];
    if (p_col >= 0) {
      elements.semi_latus_rectum_m = values[static_cast<std::size_t>(p_col)];
    } else if (values.size() >= 7U) {
      elements.semi_latus_rectum_m = values[6];
    }

    const auto state = osculating_orbital_elements_to_state_eci(elements, mu);
    *out << fmt::format(
        "{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e},{}\n",
        state.position_m.x,
        state.position_m.y,
        state.position_m.z,
        state.velocity_mps.x,
        state.velocity_mps.y,
        state.velocity_mps.z,
        static_cast<int>(state.status));
  }

  spdlog::info("wrote ECI-state batch output: {}", use_stdio_path(output_csv) ? "stdout" : output_csv.string());
  return 0;
}

int run_eci_to_osc_single(const std::vector<std::string>& args, const double mu) {
  using namespace astroforces::core;

  if (args.size() != 6U) {
    spdlog::error("eci_to_osc expects either 6 scalar values or <input_csv> <output_csv>");
    return 2;
  }

  const Vec3 position_m{
      std::stod(args[0]),
      std::stod(args[1]),
      std::stod(args[2]),
  };
  const Vec3 velocity_mps{
      std::stod(args[3]),
      std::stod(args[4]),
      std::stod(args[5]),
  };

  const auto elements = state_eci_to_osculating_orbital_elements(position_m, velocity_mps, mu);
  if (elements.status != Status::Ok) {
    spdlog::error("eci_to_osc failed: status={}", static_cast<int>(elements.status));
    return 3;
  }

  fmt::print("a_m={:.17e}\n", elements.semi_major_axis_m);
  fmt::print("p_m={:.17e}\n", elements.semi_latus_rectum_m);
  fmt::print("e={:.17e}\n", elements.eccentricity);
  fmt::print("i_rad={:.17e}\n", elements.inclination_rad);
  fmt::print("raan_rad={:.17e}\n", elements.raan_rad);
  fmt::print("argp_rad={:.17e}\n", elements.argument_of_periapsis_rad);
  fmt::print("nu_rad={:.17e}\n", elements.true_anomaly_rad);
  fmt::print("E_rad={:.17e}\n", elements.eccentric_anomaly_rad);
  fmt::print("M_rad={:.17e}\n", elements.mean_anomaly_rad);
  fmt::print("circular={} equatorial={}\n", elements.circular ? 1 : 0, elements.equatorial ? 1 : 0);
  return 0;
}

int run_osc_to_eci_single(const std::vector<std::string>& args, const double mu) {
  using namespace astroforces::core;

  if (args.size() != 6U && args.size() != 7U) {
    spdlog::error("osc_to_eci expects 6 or 7 scalar values, or <input_csv> <output_csv>");
    return 4;
  }

  ClassicalOrbitalElements elements{};
  elements.semi_major_axis_m = std::stod(args[0]);
  elements.eccentricity = std::stod(args[1]);
  elements.inclination_rad = std::stod(args[2]);
  elements.raan_rad = std::stod(args[3]);
  elements.argument_of_periapsis_rad = std::stod(args[4]);
  elements.true_anomaly_rad = std::stod(args[5]);
  if (args.size() == 7U) {
    elements.semi_latus_rectum_m = std::stod(args[6]);
  }

  const auto state = osculating_orbital_elements_to_state_eci(elements, mu);
  if (state.status != Status::Ok) {
    spdlog::error("osc_to_eci failed: status={}", static_cast<int>(state.status));
    return 5;
  }

  fmt::print("r_eci_m={:.17e},{:.17e},{:.17e}\n", state.position_m.x, state.position_m.y, state.position_m.z);
  fmt::print("v_eci_mps={:.17e},{:.17e},{:.17e}\n", state.velocity_mps.x, state.velocity_mps.y, state.velocity_mps.z);
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  using namespace astroforces::core;

  CLI::App app{"ECI/osculating-orbital-element conversion utility"};
  app.require_subcommand(1);
  app.footer(
      "Examples:\n"
      "  orbital_elements_cli eci_to_osc 7000000 0 0 0 7546.053290107542 0\n"
      "  orbital_elements_cli eci_to_osc states.csv kepler.csv\n"
      "  cat states.csv | orbital_elements_cli eci_to_osc - - > kepler.csv\n"
      "  orbital_elements_cli osc_to_eci 7000000 0 0 0 0 0\n"
      "  cat kepler.csv | orbital_elements_cli osc_to_eci - - > states_back.csv");

  std::vector<std::string> eci_args;
  double eci_mu = constants::kEarthMuM3S2;
  auto* eci_to_osc = app.add_subcommand("eci_to_osc", "Convert an ECI state to standard classical osculating elements");
  eci_to_osc->add_option("args", eci_args, "Either 6 scalar values or <input_csv> <output_csv>")->required();
  eci_to_osc->add_option("--mu", eci_mu, "Central-body gravitational parameter in m^3/s^2")
      ->capture_default_str();

  std::vector<std::string> osc_args;
  double osc_mu = constants::kEarthMuM3S2;
  auto* osc_to_eci = app.add_subcommand("osc_to_eci", "Convert standard classical osculating elements to an ECI state");
  osc_to_eci->add_option("args", osc_args, "Either 6-7 scalar values or <input_csv> <output_csv>")->required();
  osc_to_eci->add_option("--mu", osc_mu, "Central-body gravitational parameter in m^3/s^2")
      ->capture_default_str();

  CLI11_PARSE(app, argc, argv);

  if (*eci_to_osc) {
    if (eci_args.size() == 2U) {
      return run_eci_to_osc_batch(eci_args[0], eci_args[1], eci_mu);
    }
    return run_eci_to_osc_single(eci_args, eci_mu);
  }

  if (*osc_to_eci) {
    if (osc_args.size() == 2U) {
      return run_osc_to_eci_batch(osc_args[0], osc_args[1], osc_mu);
    }
    return run_osc_to_eci_single(osc_args, osc_mu);
  }

  return 0;
}

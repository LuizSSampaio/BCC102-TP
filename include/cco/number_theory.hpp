#pragma once

#include "cco/model.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace cco {

struct ModularClassStats {
  std::size_t residue = 0;
  std::size_t count = 0;
  double mean_length = 0.0;
  double mean_radius = 0.0;
  std::vector<std::size_t> segment_ids;
};

struct ModularAnalysisB0 {
  std::string filename;
  std::size_t modulus = 7;
  std::size_t total_segments = 0;
  std::vector<ModularClassStats> classes;
  bool pairwise_disjoint = false;
  bool sum_equals_total = false;
  bool covers_all_indices = false;
  bool partition_valid = false;
  std::size_t pigeonhole_lower_bound = 0;
  std::size_t max_class_size = 0;
  std::size_t max_class_residue = 0;
  bool pigeonhole_verified = false;
};

struct EuclidStep {
  std::int64_t a = 0;
  std::int64_t b = 0;
  std::int64_t quotient = 0;
  std::int64_t remainder = 0;
};

struct ExtendedEuclidResult {
  std::int64_t gcd = 0;
  std::int64_t x = 0;
  std::int64_t y = 0;
  std::vector<EuclidStep> steps;
};

struct CrtSolution {
  std::int64_t a1 = 0;
  std::int64_t m1 = 0;
  std::int64_t a2 = 0;
  std::int64_t m2 = 0;
  ExtendedEuclidResult euclid_m1_m2;
  bool solvable = false;
  std::int64_t combined_modulus = 0; // M = m1 * m2
  std::int64_t M1 = 0;               // M / m1 = m2
  std::int64_t M2 = 0;               // M / m2 = m1
  std::int64_t y1 = 0;               // M1^{-1} mod m1
  std::int64_t y2 = 0;               // M2^{-1} mod m2
  std::int64_t x0 = 0;               // unique solution in [0, M-1]
};

struct CrtQueryResult {
  std::string filename;
  std::size_t total_segments = 0;
  CrtSolution solution;
  std::vector<std::size_t> direct_ids;
  std::vector<std::size_t> equivalent_ids;
  bool methods_coincide = false;
  std::size_t selected_count = 0;
  double mean_length = 0.0;
  double mean_radius = 0.0;
};

// Part B0: residue classes modulo m (default m = 7), partition validation, and
// pigeonhole principle lower bound.
ModularAnalysisB0 analyze_modular_b0(const VascularModel& model, std::size_t modulus = 7);

// Euclidean algorithm + extended coefficients (a*x + b*y = gcd) with recorded steps.
ExtendedEuclidResult extended_euclid(std::int64_t a, std::int64_t b);

// Chinese Remainder Theorem for two congruences:
//   i === a1 (mod m1) and i === a2 (mod m2).
CrtSolution solve_crt(std::int64_t a1, std::int64_t m1, std::int64_t a2, std::int64_t m2);

// Part B1 (Track 1): runs direct selection (two congruences) and equivalent
// selection (single CRT congruence) and verifies both sets coincide.
CrtQueryResult run_crt_query(const VascularModel& model, std::int64_t a1 = 1, std::int64_t m1 = 4,
                             std::int64_t a2 = 2, std::int64_t m2 = 5);

} // namespace cco

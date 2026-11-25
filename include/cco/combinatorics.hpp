#pragma once

#include "cco/model.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace cco {

struct ModelCountingRow {
  std::string filename;
  std::size_t k = 0;
  std::size_t T = 0;
  std::size_t B = 0;
  std::size_t V = 0;
  std::size_t E = 0;
  bool rel_B_eq_T_minus_1 = false;
  bool rel_E_eq_2T_minus_1 = false;
  bool rel_V_eq_2T = false;
  bool rel_V_eq_E_plus_1 = false;
  bool rel_V_eq_1_plus_B_plus_T = false;
  bool handshaking_out_matches = false;
  bool handshaking_in_matches = false;
  std::uint64_t pairs_two_terminals = 0; // C(T, 2)
  std::uint64_t ordered_pairs_B_T = 0;   // B * T
};

struct InclusionExclusionResult {
  std::string filename;
  std::size_t total_segments = 0;
  double median_radius = 0.0;
  double median_length = 0.0;
  std::size_t size_A = 0;         // |{i : r_i >= r_med}|
  std::size_t size_B = 0;         // |{i : l_i >= l_med}|
  std::size_t size_A_inter_B = 0; // |A \cap B|
  std::size_t size_A_union_B = 0; // |A \cup B|
  std::size_t size_A_minus_B = 0;
  std::size_t size_B_minus_A = 0;
  std::size_t size_neither = 0;
  bool identity_holds = false; // |A \cup B| == |A| + |B| - |A \cap B|
};

struct RecurrenceRow {
  std::size_t k = 0;
  std::size_t n_term_expected = 0; // 16 * 2^k
  std::string filename;
  std::size_t observed_V = 0;
  std::size_t observed_E = 0;
  std::size_t observed_B = 0;
  std::size_t observed_T = 0;
  std::size_t recursive_E = 0;      // e_0 = 31, e_k = 2*e_{k-1} + 1
  std::size_t iteration_E = 0;      // 2^k * e_0 + (2^k - 1)
  std::size_t characteristic_E = 0; // 32 * 2^k - 1
  std::size_t recursive_B = 0;      // b_0 = 15, b_k = 2*b_{k-1} + 1
  std::size_t closed_B = 0;         // 16 * 2^k - 1
  std::size_t recursive_V = 0;      // v_0 = 32, v_k = 2*v_{k-1}
  std::size_t closed_V = 0;         // 32 * 2^k
  bool matches_all = false;
};

struct RecurrenceAnalysis {
  std::vector<RecurrenceRow> rows;
  bool all_models_match = false;
};

// Computes C(n, 2) = n*(n-1)/2.
std::uint64_t binomial_choose_2(std::size_t n);

// Computes the exact median of a non-empty vector of doubles.
double compute_median(std::vector<double> values);

// Part C1 (C1.1-C1.3): counting invariants, C(T,2), and B*T across all models.
std::vector<ModelCountingRow> analyze_counting_c1(const std::vector<VascularModel>& models);

// Part C1.4: Inclusion-Exclusion on 3D geometry (median radius and median length).
InclusionExclusionResult analyze_inclusion_exclusion(const VascularModel& model);

// Part C2 (C2.1-C2.4): recurrences for e_k, b_k, and v_k across the ordered sequence.
RecurrenceAnalysis analyze_recurrences_c2(const std::vector<VascularModel>& models);

} // namespace cco

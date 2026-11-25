#include "cco/combinatorics.hpp"

#include <algorithm>
#include <stdexcept>

namespace cco {

std::uint64_t binomial_choose_2(std::size_t n) {
  if (n < 2) {
    return 0;
  }
  const auto nu = static_cast<std::uint64_t>(n);
  return (nu * (nu - 1)) / 2;
}

double compute_median(std::vector<double> values) {
  if (values.empty()) {
    throw std::invalid_argument("Nao e possivel calcular a mediana de um vetor vazio.");
  }
  std::ranges::sort(values);
  const std::size_t n = values.size();
  if (n % 2 == 1) {
    return values[n / 2];
  }
  return (values[(n / 2) - 1] + values[n / 2]) / 2.0;
}

std::vector<ModelCountingRow> analyze_counting_c1(const std::vector<VascularModel>& models) {
  std::vector<ModelCountingRow> rows;
  rows.reserve(models.size());

  for (std::size_t idx = 0; idx < models.size(); ++idx) {
    const auto& m = models[idx];
    ModelCountingRow row;
    row.filename = m.filename;
    row.k = idx;
    row.T = m.structure.terminals.size();
    row.B = m.structure.bifurcations.size();
    row.V = m.points.size();
    row.E = m.segments.size();

    row.rel_B_eq_T_minus_1 = (row.T > 0) && (row.B == row.T - 1);
    row.rel_E_eq_2T_minus_1 = (row.T > 0) && (row.E == 2 * row.T - 1);
    row.rel_V_eq_2T = (row.V == 2 * row.T);
    row.rel_V_eq_E_plus_1 = (row.V == row.E + 1);
    row.rel_V_eq_1_plus_B_plus_T = (row.V == 1 + row.B + row.T);

    // Double counting (handshaking lemma for directed trees)
    std::size_t sum_out = 0;
    std::size_t sum_in = 0;
    for (const std::size_t d : m.structure.out_degree) {
      sum_out += d;
    }
    for (const std::size_t d : m.structure.in_degree) {
      sum_in += d;
    }
    row.handshaking_out_matches = (sum_out == row.E) && (1 + 2 * row.B == row.E);
    row.handshaking_in_matches = (sum_in == row.E) && (row.B + row.T == row.E);

    // C1.2: C(T, 2)
    row.pairs_two_terminals = binomial_choose_2(row.T);
    // C1.3: B * T
    row.ordered_pairs_B_T = static_cast<std::uint64_t>(row.B) * static_cast<std::uint64_t>(row.T);

    rows.push_back(row);
  }
  return rows;
}

InclusionExclusionResult analyze_inclusion_exclusion(const VascularModel& model) {
  InclusionExclusionResult res;
  res.filename = model.filename;
  res.total_segments = model.segments.size();
  if (model.segments.empty()) {
    return res;
  }

  std::vector<double> radii;
  std::vector<double> lengths;
  radii.reserve(model.segments.size());
  lengths.reserve(model.segments.size());
  for (const auto& seg : model.segments) {
    radii.push_back(seg.radius);
    lengths.push_back(seg.length);
  }

  res.median_radius = compute_median(radii);
  res.median_length = compute_median(lengths);

  for (const auto& seg : model.segments) {
    const bool in_A = (seg.radius >= res.median_radius);
    const bool in_B = (seg.length >= res.median_length);
    if (in_A) {
      res.size_A++;
    }
    if (in_B) {
      res.size_B++;
    }
    if (in_A && in_B) {
      res.size_A_inter_B++;
    }
    if (in_A || in_B) {
      res.size_A_union_B++;
    }
    if (in_A && !in_B) {
      res.size_A_minus_B++;
    }
    if (!in_A && in_B) {
      res.size_B_minus_A++;
    }
    if (!in_A && !in_B) {
      res.size_neither++;
    }
  }

  res.identity_holds = (res.size_A_union_B == (res.size_A + res.size_B - res.size_A_inter_B));
  return res;
}

RecurrenceAnalysis analyze_recurrences_c2(const std::vector<VascularModel>& models) {
  RecurrenceAnalysis analysis;
  analysis.all_models_match = !models.empty();
  analysis.rows.reserve(models.size());

  std::size_t prev_E = 0;
  std::size_t prev_B = 0;
  std::size_t prev_V = 0;

  for (std::size_t k = 0; k < models.size(); ++k) {
    const auto& m = models[k];
    const std::size_t pow2_k = (std::size_t{1} << k);

    RecurrenceRow row;
    row.k = k;
    row.n_term_expected = 16 * pow2_k;
    row.filename = m.filename;
    row.observed_V = m.points.size();
    row.observed_E = m.segments.size();
    row.observed_B = m.structure.bifurcations.size();
    row.observed_T = m.structure.terminals.size();

    // Primary recurrence: e_0 = 31, e_k = 2 * e_{k-1} + 1
    if (k == 0) {
      row.recursive_E = 31;
      row.recursive_B = 15;
      row.recursive_V = 32;
    } else {
      row.recursive_E = 2 * prev_E + 1;
      row.recursive_B = 2 * prev_B + 1;
      row.recursive_V = 2 * prev_V;
    }
    prev_E = row.recursive_E;
    prev_B = row.recursive_B;
    prev_V = row.recursive_V;

    // Closed form via iteration: 2^k * e_0 + (2^k - 1)
    row.iteration_E = pow2_k * 31 + (pow2_k - 1);
    // Closed form via characteristic equation: 32 * 2^k - 1
    row.characteristic_E = 32 * pow2_k - 1;

    // Secondary closed forms: b_k = 16 * 2^k - 1 = T_k - 1, v_k = 32 * 2^k = 2 * T_k
    row.closed_B = 16 * pow2_k - 1;
    row.closed_V = 32 * pow2_k;

    row.matches_all = (row.observed_T == row.n_term_expected) &&
                      (row.observed_E == row.recursive_E) && (row.observed_E == row.iteration_E) &&
                      (row.observed_E == row.characteristic_E) &&
                      (row.observed_B == row.recursive_B) && (row.observed_B == row.closed_B) &&
                      (row.observed_V == row.recursive_V) && (row.observed_V == row.closed_V);

    if (!row.matches_all) {
      analysis.all_models_match = false;
    }
    analysis.rows.push_back(row);
  }

  return analysis;
}

} // namespace cco

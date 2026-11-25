#include "cco/number_theory.hpp"

#include <stdexcept>

namespace cco {

namespace {

std::int64_t positive_mod(std::int64_t value, std::int64_t modulus) {
  const std::int64_t rem = value % modulus;
  return rem < 0 ? rem + modulus : rem;
}

} // namespace

ModularAnalysisB0 analyze_modular_b0(const VascularModel& model, std::size_t modulus) {
  if (modulus == 0) {
    throw std::invalid_argument("Modulo deve ser estritamente positivo.");
  }

  ModularAnalysisB0 result;
  result.filename = model.filename;
  result.modulus = modulus;
  result.total_segments = model.segments.size();
  result.classes.resize(modulus);

  std::vector<double> sum_lengths(modulus, 0.0);
  std::vector<double> sum_radii(modulus, 0.0);

  for (std::size_t j = 0; j < modulus; ++j) {
    result.classes[j].residue = j;
  }

  for (const auto& seg : model.segments) {
    const std::size_t r = seg.id % modulus;
    result.classes[r].segment_ids.push_back(seg.id);
    sum_lengths[r] += seg.length;
    sum_radii[r] += seg.radius;
  }

  for (std::size_t j = 0; j < modulus; ++j) {
    const std::size_t cnt = result.classes[j].segment_ids.size();
    result.classes[j].count = cnt;
    if (cnt > 0) {
      result.classes[j].mean_length = sum_lengths[j] / static_cast<double>(cnt);
      result.classes[j].mean_radius = sum_radii[j] / static_cast<double>(cnt);
    }
    if (cnt > result.max_class_size) {
      result.max_class_size = cnt;
      result.max_class_residue = j;
    }
  }

  // B0.2: Explicitly verify partition properties
  // (i) Pairwise disjointness C_j \cap C_k = \emptyset for j != k
  bool disjoint = true;
  std::vector<std::size_t> membership_count(result.total_segments, 0);
  for (std::size_t j = 0; j < modulus; ++j) {
    for (const std::size_t id : result.classes[j].segment_ids) {
      if (id >= result.total_segments) {
        disjoint = false;
      } else {
        membership_count[id]++;
        if (membership_count[id] > 1) {
          disjoint = false;
        }
      }
    }
  }

  // (ii) Sum of cardinalities equals |E| and every index appears once
  std::size_t sum_sizes = 0;
  for (const auto& cls : result.classes) {
    sum_sizes += cls.count;
  }
  bool covers_all = true;
  for (const std::size_t cnt : membership_count) {
    if (cnt != 1) {
      covers_all = false;
      break;
    }
  }

  result.pairwise_disjoint = disjoint;
  result.sum_equals_total = (sum_sizes == result.total_segments);
  result.covers_all_indices = covers_all;
  result.partition_valid =
      result.pairwise_disjoint && result.sum_equals_total && result.covers_all_indices;

  // B0.3: Pigeonhole Principle lower bound ceil(|E| / m)
  result.pigeonhole_lower_bound =
      result.total_segments == 0 ? 0 : (result.total_segments + modulus - 1) / modulus;
  result.pigeonhole_verified = (result.max_class_size >= result.pigeonhole_lower_bound);

  return result;
}

ExtendedEuclidResult extended_euclid(std::int64_t a, std::int64_t b) {
  ExtendedEuclidResult res;

  std::int64_t old_r = a;
  std::int64_t r = b;
  std::int64_t old_s = 1;
  std::int64_t s = 0;
  std::int64_t old_t = 0;
  std::int64_t t = 1;

  while (r != 0) {
    const std::int64_t q = old_r / r;
    const std::int64_t rem = old_r - q * r;
    res.steps.push_back(EuclidStep{old_r, r, q, rem});

    old_r = r;
    r = rem;

    const std::int64_t next_s = old_s - q * s;
    old_s = s;
    s = next_s;

    const std::int64_t next_t = old_t - q * t;
    old_t = t;
    t = next_t;
  }

  if (old_r < 0) {
    old_r = -old_r;
    old_s = -old_s;
    old_t = -old_t;
  }

  res.gcd = old_r;
  res.x = old_s;
  res.y = old_t;
  return res;
}

CrtSolution solve_crt(std::int64_t a1, std::int64_t m1, std::int64_t a2, std::int64_t m2) {
  if (m1 <= 0 || m2 <= 0) {
    throw std::invalid_argument("Os modulos do TCR devem ser positivos.");
  }

  CrtSolution sol;
  sol.a1 = positive_mod(a1, m1);
  sol.m1 = m1;
  sol.a2 = positive_mod(a2, m2);
  sol.m2 = m2;
  sol.euclid_m1_m2 = extended_euclid(m1, m2);

  if (sol.euclid_m1_m2.gcd != 1) {
    sol.solvable = false;
    return sol;
  }

  sol.solvable = true;
  sol.combined_modulus = m1 * m2;
  sol.M1 = m2; // M / m1
  sol.M2 = m1; // M / m2

  // M1 * y1 === 1 (mod m1) and M2 * y2 === 1 (mod m2)
  const ExtendedEuclidResult inv1 = extended_euclid(sol.M1, m1);
  const ExtendedEuclidResult inv2 = extended_euclid(sol.M2, m2);
  sol.y1 = positive_mod(inv1.x, m1);
  sol.y2 = positive_mod(inv2.x, m2);

  const std::int64_t raw_x = sol.a1 * sol.M1 * sol.y1 + sol.a2 * sol.M2 * sol.y2;
  sol.x0 = positive_mod(raw_x, sol.combined_modulus);
  return sol;
}

CrtQueryResult run_crt_query(const VascularModel& model, std::int64_t a1, std::int64_t m1,
                             std::int64_t a2, std::int64_t m2) {
  CrtQueryResult result;
  result.filename = model.filename;
  result.total_segments = model.segments.size();
  result.solution = solve_crt(a1, m1, a2, m2);

  if (!result.solution.solvable) {
    return result;
  }

  const auto m1_u = static_cast<std::size_t>(result.solution.m1);
  const auto m2_u = static_cast<std::size_t>(result.solution.m2);
  const auto a1_u = static_cast<std::size_t>(result.solution.a1);
  const auto a2_u = static_cast<std::size_t>(result.solution.a2);
  const auto M_u = static_cast<std::size_t>(result.solution.combined_modulus);
  const auto x0_u = static_cast<std::size_t>(result.solution.x0);

  double sum_len = 0.0;
  double sum_rad = 0.0;

  for (const auto& seg : model.segments) {
    // Method 1: Direct selection by the two simultaneous congruences
    if ((seg.id % m1_u == a1_u) && (seg.id % m2_u == a2_u)) {
      result.direct_ids.push_back(seg.id);
      sum_len += seg.length;
      sum_rad += seg.radius;
    }
    // Method 2: Selection by the single equivalent CRT congruence
    if (seg.id % M_u == x0_u) {
      result.equivalent_ids.push_back(seg.id);
    }
  }

  result.methods_coincide = (result.direct_ids == result.equivalent_ids);
  result.selected_count = result.direct_ids.size();
  if (result.selected_count > 0) {
    result.mean_length = sum_len / static_cast<double>(result.selected_count);
    result.mean_radius = sum_rad / static_cast<double>(result.selected_count);
  }
  return result;
}

} // namespace cco

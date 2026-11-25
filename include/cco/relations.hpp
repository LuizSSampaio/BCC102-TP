#pragma once

#include "cco/model.hpp"
#include "cco/number_theory.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

namespace cco {

using BoolMatrix = std::vector<std::vector<std::uint8_t>>;

struct DirectedPair {
  std::size_t u = 0;
  std::size_t v = 0;
};

struct TransitivityWitness {
  std::size_t u = 0;
  std::size_t v = 0;
  std::size_t w = 0;
};

struct RelationProperties {
  bool reflexive = false;
  std::optional<std::size_t> non_reflexive_example;
  bool irreflexive = false;
  bool symmetric = false;
  std::optional<DirectedPair> non_symmetric_example;
  bool antisymmetric = false;
  std::optional<DirectedPair> non_antisymmetric_example;
  bool asymmetric = false;
  bool transitive = false;
  std::optional<TransitivityWitness> non_transitive_example;
};

struct ParentFunctionAnalysis {
  bool is_well_defined_function = false;
  bool root_zero_has_no_parent = false;
  bool all_non_root_have_unique_parent = false;
  bool is_injective = false;
  std::optional<std::tuple<std::size_t, std::size_t, std::size_t>> non_injective_witness;
  bool is_surjective_onto_V = false;
  std::optional<std::size_t> non_surjective_witness;
  std::vector<std::optional<std::size_t>> parent_of;
};

struct ClosureAndOrderAnalysis {
  std::size_t num_vertices = 0;
  BoolMatrix matrix_R;
  BoolMatrix matrix_R_plus;
  BoolMatrix matrix_leq;
  std::size_t closure_iterations = 0;
  std::vector<std::size_t> pairs_added_per_iteration;
  std::size_t count_R = 0;
  std::size_t count_R_plus = 0;
  std::size_t count_indirect = 0;
  std::size_t count_leq = 0;
  std::vector<DirectedPair> pairs_R_plus;
  std::vector<TransitivityWitness> indirect_witnesses;
  bool leq_reflexive = false;
  bool leq_antisymmetric = false;
  bool leq_transitive = false;
  bool is_partial_order = false;
  bool is_total_order = false;
  bool root_zero_is_minimum = false;
  std::size_t incomparable_pair_count = 0;
  std::vector<DirectedPair> incomparable_examples;
};

struct ModularEquivalenceAnalysis {
  std::string filename;
  std::size_t modulus = 7;
  std::size_t set_size = 0;
  bool reflexive = false;
  bool symmetric = false;
  bool transitive = false;
  bool is_equivalence_relation = false;
  std::vector<std::vector<std::size_t>> equivalence_classes;
  bool matches_module_b_classes = false;
};

struct RelationAnalysisD {
  std::string filename;
  std::size_t num_vertices = 0;
  std::vector<DirectedPair> pairs_R;
  std::vector<std::size_t> domain_R;
  std::vector<std::size_t> image_R;
  std::vector<DirectedPair> pairs_R_inv;
  std::vector<std::size_t> domain_R_inv;
  std::vector<std::size_t> image_R_inv;
  RelationProperties props_R;
  RelationProperties props_R_inv;
  bool R_inv_is_functional = false;
  ParentFunctionAnalysis parent_function;
  std::optional<ClosureAndOrderAnalysis> closure_and_order;
};

// Module D (D1-D5): parent-child relation R, inverse R^-1, properties,
// parent function p: V\{0} -> V, and (when compute_dense_matrix is true)
// dense matrix A, iterative transitive closure R^+, and partial order <=.
RelationAnalysisD analyze_relations_d(const VascularModel& model, bool compute_dense_matrix);

// Module D (D6): modular equivalence relation i ~ j <=> i === j (mod m) on S = {0,...,|E|-1}.
ModularEquivalenceAnalysis analyze_modular_equivalence_d6(const VascularModel& model,
                                                          const ModularAnalysisB0& b0);

} // namespace cco

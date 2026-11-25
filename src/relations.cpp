#include "cco/relations.hpp"

#include <algorithm>
#include <set>
#include <utility>

namespace cco {

namespace {

RelationProperties analyze_properties_sparse(std::size_t n_vertices,
                                             const std::vector<DirectedPair>& pairs) {
  RelationProperties props;

  // Build adjacency lists for efficient property verification without dense O(n^2) storage
  std::vector<std::vector<std::size_t>> adj(n_vertices);
  std::set<std::pair<std::size_t, std::size_t>> edge_set;
  for (const auto& p : pairs) {
    if (p.u < n_vertices && p.v < n_vertices) {
      adj[p.u].push_back(p.v);
      edge_set.insert({p.u, p.v});
    }
  }

  // 1. Reflexivity & Irreflexivity
  props.reflexive = true;
  props.irreflexive = true;
  for (std::size_t v = 0; v < n_vertices; ++v) {
    const bool has_loop = edge_set.contains({v, v});
    if (!has_loop) {
      if (props.reflexive) {
        props.reflexive = false;
        props.non_reflexive_example = v;
      }
    } else {
      props.irreflexive = false;
    }
  }

  // 2. Symmetry, Antisymmetry, Asymmetry
  props.symmetric = true;
  props.antisymmetric = true;
  props.asymmetric = true;
  for (const auto& p : pairs) {
    const bool has_rev = edge_set.contains({p.v, p.u});
    if (!has_rev && props.symmetric) {
      props.symmetric = false;
      props.non_symmetric_example = p;
    }
    if (has_rev) {
      props.asymmetric = false;
      if (p.u != p.v && props.antisymmetric) {
        props.antisymmetric = false;
        props.non_antisymmetric_example = p;
      }
    }
  }

  // 3. Transitivity: for every (u,v) and (v,w), check if (u,w) exists
  props.transitive = true;
  for (const auto& p : pairs) {
    for (const std::size_t w : adj[p.v]) {
      if (!edge_set.contains({p.u, w})) {
        props.transitive = false;
        props.non_transitive_example = TransitivityWitness{p.u, p.v, w};
        break;
      }
    }
    if (!props.transitive) {
      break;
    }
  }

  return props;
}

ClosureAndOrderAnalysis compute_closure_and_order(std::size_t n_vertices,
                                                  const std::vector<DirectedPair>& pairs) {
  ClosureAndOrderAnalysis res;
  res.num_vertices = n_vertices;
  res.matrix_R.assign(n_vertices, std::vector<std::uint8_t>(n_vertices, 0));

  for (const auto& p : pairs) {
    if (p.u < n_vertices && p.v < n_vertices && res.matrix_R[p.u][p.v] == 0) {
      res.matrix_R[p.u][p.v] = 1;
      res.count_R++;
    }
  }

  // D5: Iterative fixed-point transitive closure (as instructed in Section 7):
  // Start from R and repeatedly add (u, w) whenever (u, v) and (v, w) are present,
  // until no new pairs are added.
  res.matrix_R_plus = res.matrix_R;
  while (true) {
    std::vector<DirectedPair> to_add;
    std::vector<TransitivityWitness> new_witnesses;
    for (std::size_t u = 0; u < n_vertices; ++u) {
      for (std::size_t v = 0; v < n_vertices; ++v) {
        if (res.matrix_R_plus[u][v] == 0) {
          continue;
        }
        for (std::size_t w = 0; w < n_vertices; ++w) {
          if (res.matrix_R_plus[v][w] == 1 && res.matrix_R_plus[u][w] == 0) {
            to_add.push_back(DirectedPair{u, w});
            if (new_witnesses.size() < 20) {
              new_witnesses.push_back(TransitivityWitness{u, v, w});
            }
          }
        }
      }
    }

    std::size_t actually_added = 0;
    for (const auto& p : to_add) {
      if (res.matrix_R_plus[p.u][p.v] == 0) {
        res.matrix_R_plus[p.u][p.v] = 1;
        actually_added++;
      }
    }
    if (actually_added == 0) {
      break;
    }
    res.closure_iterations++;
    res.pairs_added_per_iteration.push_back(actually_added);
    for (const auto& wit : new_witnesses) {
      if (res.indirect_witnesses.size() < 15) {
        res.indirect_witnesses.push_back(wit);
      }
    }
  }

  // Collect R^+ pairs and counts
  for (std::size_t u = 0; u < n_vertices; ++u) {
    for (std::size_t v = 0; v < n_vertices; ++v) {
      if (res.matrix_R_plus[u][v] == 1) {
        res.count_R_plus++;
        res.pairs_R_plus.push_back(DirectedPair{u, v});
        if (res.matrix_R[u][v] == 0) {
          res.count_indirect++;
        }
      }
    }
  }

  // Construct reflexive-transitive closure <= = R^+ U {(v, v) : v in V}
  res.matrix_leq = res.matrix_R_plus;
  for (std::size_t v = 0; v < n_vertices; ++v) {
    res.matrix_leq[v][v] = 1;
  }
  for (std::size_t u = 0; u < n_vertices; ++u) {
    for (std::size_t v = 0; v < n_vertices; ++v) {
      if (res.matrix_leq[u][v] == 1) {
        res.count_leq++;
      }
    }
  }

  // Verify partial order axioms directly on matrix_leq
  res.leq_reflexive = true;
  res.leq_antisymmetric = true;
  res.leq_transitive = true;

  for (std::size_t u = 0; u < n_vertices; ++u) {
    if (res.matrix_leq[u][u] == 0) {
      res.leq_reflexive = false;
    }
    for (std::size_t v = 0; v < n_vertices; ++v) {
      if (u != v && res.matrix_leq[u][v] == 1 && res.matrix_leq[v][u] == 1) {
        res.leq_antisymmetric = false;
      }
      if (res.matrix_leq[u][v] == 1) {
        for (std::size_t w = 0; w < n_vertices; ++w) {
          if (res.matrix_leq[v][w] == 1 && res.matrix_leq[u][w] == 0) {
            res.leq_transitive = false;
          }
        }
      }
    }
  }
  res.is_partial_order = res.leq_reflexive && res.leq_antisymmetric && res.leq_transitive;

  // Check if root 0 is the minimum element (0 <= v for all v in V)
  res.root_zero_is_minimum = (n_vertices > 0);
  for (std::size_t v = 0; v < n_vertices; ++v) {
    if (res.matrix_leq[0][v] == 0) {
      res.root_zero_is_minimum = false;
      break;
    }
  }

  // Count incomparable unordered pairs {u, v} with u < v
  for (std::size_t u = 0; u < n_vertices; ++u) {
    for (std::size_t v = u + 1; v < n_vertices; ++v) {
      if (res.matrix_leq[u][v] == 0 && res.matrix_leq[v][u] == 0) {
        res.incomparable_pair_count++;
        if (res.incomparable_examples.size() < 10) {
          res.incomparable_examples.push_back(DirectedPair{u, v});
        }
      }
    }
  }
  res.is_total_order = res.is_partial_order && (res.incomparable_pair_count == 0);

  return res;
}

} // namespace

RelationAnalysisD analyze_relations_d(const VascularModel& model, bool compute_dense_matrix) {
  RelationAnalysisD res;
  res.filename = model.filename;
  res.num_vertices = model.points.size();

  // D1: Parent-child relation R and inverse relation R^-1
  res.pairs_R.reserve(model.segments.size());
  res.pairs_R_inv.reserve(model.segments.size());

  std::vector<bool> in_dom(res.num_vertices, false);
  std::vector<bool> in_img(res.num_vertices, false);

  for (const auto& seg : model.segments) {
    res.pairs_R.push_back(DirectedPair{seg.u, seg.v});
    res.pairs_R_inv.push_back(DirectedPair{seg.v, seg.u});
    if (seg.u < res.num_vertices) {
      in_dom[seg.u] = true;
    }
    if (seg.v < res.num_vertices) {
      in_img[seg.v] = true;
    }
  }

  for (std::size_t v = 0; v < res.num_vertices; ++v) {
    if (in_dom[v]) {
      res.domain_R.push_back(v);
      res.image_R_inv.push_back(v);
    }
    if (in_img[v]) {
      res.image_R.push_back(v);
      res.domain_R_inv.push_back(v);
    }
  }

  // D3: Properties of R and R^-1
  res.props_R = analyze_properties_sparse(res.num_vertices, res.pairs_R);
  res.props_R_inv = analyze_properties_sparse(res.num_vertices, res.pairs_R_inv);

  // D4: Parent function p: V\{0} -> V
  res.parent_function.parent_of.assign(res.num_vertices, std::nullopt);
  std::vector<std::size_t> parent_counts(res.num_vertices, 0);
  for (const auto& p : res.pairs_R) {
    if (p.v < res.num_vertices) {
      parent_counts[p.v]++;
      res.parent_function.parent_of[p.v] = p.u;
    }
  }

  res.parent_function.root_zero_has_no_parent = (res.num_vertices > 0) && (parent_counts[0] == 0);
  res.parent_function.all_non_root_have_unique_parent = (res.num_vertices > 0);
  for (std::size_t v = 1; v < res.num_vertices; ++v) {
    if (parent_counts[v] != 1) {
      res.parent_function.all_non_root_have_unique_parent = false;
      break;
    }
  }
  res.parent_function.is_well_defined_function =
      res.parent_function.root_zero_has_no_parent &&
      res.parent_function.all_non_root_have_unique_parent;
  res.R_inv_is_functional = res.parent_function.is_well_defined_function;

  // Check injectivity of p on V\{0}
  std::vector<std::optional<std::size_t>> first_child_of(res.num_vertices, std::nullopt);
  res.parent_function.is_injective = true;
  for (std::size_t v = 1; v < res.num_vertices; ++v) {
    const auto& parent_opt = res.parent_function.parent_of[v];
    if (!parent_opt.has_value()) {
      continue;
    }
    const std::size_t parent = *parent_opt;
    const auto& existing_child = first_child_of[parent];
    if (existing_child.has_value()) {
      res.parent_function.is_injective = false;
      if (!res.parent_function.non_injective_witness.has_value()) {
        res.parent_function.non_injective_witness = std::make_tuple(*existing_child, v, parent);
      }
    } else {
      first_child_of[parent] = v;
    }
  }

  // Check surjectivity of p onto V
  res.parent_function.is_surjective_onto_V = true;
  for (std::size_t v = 0; v < res.num_vertices; ++v) {
    if (!first_child_of[v].has_value()) {
      res.parent_function.is_surjective_onto_V = false;
      res.parent_function.non_surjective_witness = v;
      break;
    }
  }

  // D2 & D5: Dense matrix A, iterative transitive closure R^+, and partial order <=
  if (compute_dense_matrix) {
    res.closure_and_order = compute_closure_and_order(res.num_vertices, res.pairs_R);
  }

  return res;
}

ModularEquivalenceAnalysis analyze_modular_equivalence_d6(const VascularModel& model,
                                                          const ModularAnalysisB0& b0) {
  ModularEquivalenceAnalysis res;
  res.filename = model.filename;
  res.modulus = b0.modulus;
  res.set_size = model.segments.size();

  const auto equiv = [m = res.modulus](std::size_t i, std::size_t j) { return (i % m) == (j % m); };

  // Exhaustive verification of reflexivity, symmetry, and transitivity on S = {0,...,|E|-1}
  res.reflexive = true;
  res.symmetric = true;
  res.transitive = true;

  for (std::size_t i = 0; i < res.set_size; ++i) {
    if (!equiv(i, i)) {
      res.reflexive = false;
      break;
    }
  }

  for (std::size_t i = 0; i < res.set_size && res.symmetric; ++i) {
    for (std::size_t j = 0; j < res.set_size; ++j) {
      if (equiv(i, j) && !equiv(j, i)) {
        res.symmetric = false;
        break;
      }
    }
  }

  for (std::size_t i = 0; i < res.set_size && res.transitive; ++i) {
    for (std::size_t j = 0; j < res.set_size && res.transitive; ++j) {
      if (!equiv(i, j)) {
        continue;
      }
      for (std::size_t k = 0; k < res.set_size; ++k) {
        if (equiv(j, k) && !equiv(i, k)) {
          res.transitive = false;
          break;
        }
      }
    }
  }

  res.is_equivalence_relation = res.reflexive && res.symmetric && res.transitive;

  // Build equivalence classes [0], ..., [m-1]
  res.equivalence_classes.resize(res.modulus);
  for (std::size_t r = 0; r < res.modulus; ++r) {
    for (std::size_t i = 0; i < res.set_size; ++i) {
      if (equiv(i, r)) {
        res.equivalence_classes[r].push_back(i);
      }
    }
  }

  // Verify exact match with C_0, ..., C_{m-1} from Module B
  res.matches_module_b_classes = (res.equivalence_classes.size() == b0.classes.size());
  if (res.matches_module_b_classes) {
    for (std::size_t r = 0; r < res.modulus; ++r) {
      if (res.equivalence_classes[r] != b0.classes[r].segment_ids) {
        res.matches_module_b_classes = false;
        break;
      }
    }
  }

  return res;
}

} // namespace cco

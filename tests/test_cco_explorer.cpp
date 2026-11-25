#include "cco/combinatorics.hpp"
#include "cco/number_theory.hpp"
#include "cco/relations.hpp"
#include "cco/vtk_parser.hpp"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace {

int g_failures = 0;

void expect(bool condition, const std::string& label) {
  if (condition) {
    std::cout << "[PASS] " << label << "\n";
  } else {
    std::cerr << "[FAIL] " << label << "\n";
    g_failures++;
  }
}

cco::BoolMatrix warshall_closure(const cco::BoolMatrix& adj) {
  const std::size_t n = adj.size();
  cco::BoolMatrix reach = adj;
  for (std::size_t k = 0; k < n; ++k) {
    for (std::size_t i = 0; i < n; ++i) {
      if (reach[i][k] == 0) {
        continue;
      }
      for (std::size_t j = 0; j < n; ++j) {
        if (reach[k][j] == 1) {
          reach[i][j] = 1;
        }
      }
    }
  }
  return reach;
}

const cco::VascularModel& require_model(const std::vector<cco::VascularModel>& models,
                                        std::string_view name) {
  for (const auto& m : models) {
    if (m.filename == name) {
      return m;
    }
  }
  throw std::runtime_error("Modelo obrigatorio ausente nos testes: " + std::string(name));
}

void test_parser_rejects_malformed_vtk() {
  const auto tmp_path = std::filesystem::temp_directory_path() / "cco_test_malformed.vtk";
  {
    std::ofstream out(tmp_path);
    out << "# vtk DataFile Version 3.0\n"
        << "bad test\n"
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS 2 float\n"
        << "0.0 0.0 0.0\n"
        << "1.0 0.0 0.0\n"
        << "LINES 1 3\n"
        << "2 0 5\n"
        << "CELL_DATA 1\n"
        << "scalars raio float\n"
        << "LOOKUP_TABLE default\n"
        << "0.5\n";
  }
  const cco::VascularModel bad = cco::parse_vtk_file(tmp_path);
  std::filesystem::remove(tmp_path);
  expect(!bad.validation.ok && !bad.validation.all_endpoints_valid,
         "Unitario: parser rejeita ponto inexistente em LINES");
}

void test_math_primitives() {
  // Extended Euclid & CRT
  const auto eg = cco::extended_euclid(4, 5);
  expect(eg.gcd == 1 && (4 * eg.x + 5 * eg.y == 1),
         "Unitario: Algoritmo de Euclides Estendido para (4, 5)");

  const auto crt = cco::solve_crt(1, 4, 2, 5);
  expect(crt.solvable && crt.combined_modulus == 20 && crt.x0 == 17,
         "Unitario: TCR resolve i===1(mod 4), i===2(mod 5) => i===17(mod 20)");

  const auto crt_unsolvable = cco::solve_crt(1, 4, 2, 6);
  expect(!crt_unsolvable.solvable, "Unitario: TCR detecta modulos nao coprimos (4, 6)");

  // Median (odd & even)
  const double med_odd = cco::compute_median({5.0, 1.0, 3.0});
  const double med_even = cco::compute_median({10.0, 2.0, 4.0, 8.0});
  expect(std::abs(med_odd - 3.0) < 1e-12 && std::abs(med_even - 6.0) < 1e-12,
         "Unitario: mediana exata para vetores impares e pares");

  // Binomial C(T, 2)
  expect(cco::binomial_choose_2(16) == 120 && cco::binomial_choose_2(128) == 8128,
         "Unitario: C(16,2)=120 e C(128,2)=8128");
}

} // namespace

int main(int argc, char* argv[]) {
  try {
    std::filesystem::path input_dir = "dados";
    for (int i = 1; i < argc; ++i) {
      const std::string arg = argv[i];
      if ((arg == "--entrada" || arg == "--input") && i + 1 < argc) {
        input_dir = argv[++i];
      }
    }

    std::cout << "Executando suite de testes CCO Explorer (entrada: " << input_dir.string()
              << ")\n";
    test_math_primitives();
    test_parser_rejects_malformed_vtk();

    const std::vector<cco::VascularModel> models = cco::load_all_vtk_models(input_dir);

    // T1: Process all six VTKs without manual alterations
    bool all_valid = (models.size() == 6);
    for (const auto& m : models) {
      if (!m.validation.ok) {
        all_valid = false;
      }
    }
    expect(all_valid, "T1: Processa todos os seis arquivos VTK com validacao OK");

    // T2: Points, segments, terminals, and bifurcations match Table 1
    bool all_match_table1 = (models.size() == 6);
    for (const auto& m : models) {
      const auto ref = cco::get_table1_reference(m.filename);
      if (!ref.has_value()) {
        all_match_table1 = false;
        continue;
      }
      const bool ok = (m.points.size() == ref->points) && (m.segments.size() == ref->segments) &&
                      (m.structure.terminals.size() == ref->n_term) &&
                      (m.structure.bifurcations.size() == ref->bifurcations) &&
                      (m.structure.entry_points.size() == 1 && m.structure.entry_points[0] == 0) &&
                      (m.structure.out_degree[0] == 1 && m.structure.in_degree[0] == 0);
      if (!ok) {
        all_match_table1 = false;
      }
    }
    expect(all_match_table1,
           "T2: |V|, |E|, T, B e ponto de entrada 0 coincidem com a Tabela 1 nos 6 modelos");

    const auto& m016 = require_model(models, "tree3D_Nterm016.vtk");
    const auto& m128 = require_model(models, "tree3D_Nterm128.vtk");
    const auto& m512 = require_model(models, "tree3D_Nterm512.vtk");

    // T3: Parent-child matrix of Nterm016 is 32 x 32 with exactly 31 entries equal to 1
    const auto rel16 = cco::analyze_relations_d(m016, true);
    expect(
        rel16.closure_and_order.has_value() && rel16.closure_and_order->matrix_R.size() == 32 &&
            rel16.closure_and_order->matrix_R[0].size() == 32 &&
            rel16.closure_and_order->count_R == 31,
        "T3: Matriz pai-filho de Nterm016 tem dimensoes 32x32 e exatamente 31 entradas iguais a 1");

    // T4: R^+ for Nterm016 produces verified indirect ancestry & matches Warshall; <= is partial
    // order
    bool t4_ok = false;
    if (rel16.closure_and_order.has_value()) {
      const auto& cl16 = *rel16.closure_and_order;
      const auto warshall_16 = warshall_closure(cl16.matrix_R);
      const bool indirect_0_12 = (cl16.matrix_R[0][6] == 1) && (cl16.matrix_R[6][12] == 1) &&
                                 (cl16.matrix_R[0][12] == 0) && (cl16.matrix_R_plus[0][12] == 1);
      t4_ok = (cl16.matrix_R_plus == warshall_16) && (cl16.count_indirect > 0) && indirect_0_12 &&
              cl16.is_partial_order && !cl16.is_total_order && cl16.root_zero_is_minimum &&
              !cl16.incomparable_examples.empty();
    }
    expect(t4_ok,
           "T4: R^+ de Nterm016 gera ancestralidade indireta correta, coincide com Warshall e <= e "
           "ordem parcial nao total");

    // Extra D3/D4 checks on Nterm016
    expect(!rel16.props_R.reflexive && rel16.props_R.irreflexive && !rel16.props_R.symmetric &&
               rel16.props_R.antisymmetric && !rel16.props_R.transitive &&
               rel16.parent_function.is_well_defined_function &&
               !rel16.parent_function.is_injective && !rel16.parent_function.is_surjective_onto_V,
           "Modulo D (D3-D4): propriedades de R e funcao pai p: V\\{0} -> V verificadas em "
           "Nterm016");

    // T5: Inclusion-exclusion closes numerically on Nterm128
    const auto ie128 = cco::analyze_inclusion_exclusion(m128);
    const bool ie_partition = (ie128.size_A_minus_B + ie128.size_B_minus_A + ie128.size_A_inter_B +
                               ie128.size_neither) == 255;
    expect(ie128.identity_holds && ie_partition && ie128.size_A >= 128 && ie128.size_B >= 128,
           "T5: Identidade de inclusao-exclusao fecha numericamente em Nterm128");

    // Extra C1.2 / C1.3 checks
    const auto counting = cco::analyze_counting_c1(models);
    expect(counting[0].pairs_two_terminals == 120 && counting[0].ordered_pairs_B_T == 240 &&
               counting[3].pairs_two_terminals == 8128 && counting[3].ordered_pairs_B_T == 16256,
           "Modulo C (C1.2-C1.3): combinacoes C(T,2) e pares ordenados B*T conferidos");

    // T6: Recurrences reproduce all six observed values
    const auto rec = cco::analyze_recurrences_c2(models);
    expect(rec.rows.size() == 6 && rec.all_models_match,
           "T6: Recorrencias deduzidas (e_k, b_k, v_k) reproduzem exatamente os 6 modelos");

    // T7: B0 sum of 7 frequencies equals |E| of Nterm128 & B1 Track 1 CRT verifiable
    const auto b0_128 = cco::analyze_modular_b0(m128, 7);
    const auto crt_128 = cco::run_crt_query(m128, 1, 4, 2, 5);
    const auto crt_512 = cco::run_crt_query(m512, 1, 4, 2, 5);
    // Count of i in [0, 254] with i === 17 (mod 20): 17, 37, ..., 237 => (237-17)/20 + 1 = 12
    // Count of i in [0, 1022] with i === 17 (mod 20): 17, 37, ..., 1017 => (1017-17)/20 + 1 = 51
    expect(b0_128.partition_valid && b0_128.pigeonhole_lower_bound == 37 &&
               b0_128.max_class_size == 37 && b0_128.pigeonhole_verified &&
               crt_128.methods_coincide && crt_128.selected_count == 12 &&
               crt_512.methods_coincide && crt_512.selected_count == 51,
           "T7: Particao B0 (soma=255, PCP>=37) e Trilha B1 TCR (12 em Nterm128, 51 em Nterm512) "
           "verificadas");

    // Extra D6 check: modular equivalence relation on Nterm128
    const auto eq128 = cco::analyze_modular_equivalence_d6(m128, b0_128);
    expect(eq128.is_equivalence_relation && eq128.matches_module_b_classes,
           "Modulo D (D6): relacao de equivalencia modular verificada e identica a C_0..C_6");

    if (g_failures == 0) {
      std::cout << "\nTODOS OS TESTES PASSARAM COM SUCESSO!\n";
      return 0;
    }
    std::cerr << "\nTOTAL DE FALHAS: " << g_failures << "\n";
    return 1;
  } catch (const std::exception& ex) {
    std::cerr << "[EXCECAO] " << ex.what() << "\n";
    return 1;
  }
}

#include "cco/exporter.hpp"
#include "cco/vtk_parser.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace cco {

namespace {

std::string bool_pt(bool v) {
  return v ? "SIM" : "NAO";
}

std::string format_double(double val, int precision = 6) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(precision) << val;
  return oss.str();
}

} // namespace

void write_text_file(const std::filesystem::path& filepath, const std::string& content) {
  if (filepath.has_parent_path()) {
    std::filesystem::create_directories(filepath.parent_path());
  }
  std::ofstream out(filepath);
  if (!out.is_open()) {
    throw std::runtime_error("Nao foi possivel criar o arquivo: " + filepath.string());
  }
  out << content;
}

std::string format_mathematical_summary(const VascularModel& model, const ModularAnalysisB0& b0,
                                        const CrtQueryResult& crt,
                                        const RelationAnalysisD& relations,
                                        const RecurrenceAnalysis* recurrences) {
  std::ostringstream oss;
  oss << "CCO EXPLORER -- RESUMO MATEMATICO\n";
  oss << "Arquivo: " << model.filename << "\n";
  oss << "Pontos (V):                 " << model.points.size() << "\n";
  oss << "Segmentos (E):              " << model.segments.size() << "\n";
  oss << "Terminais (T):              " << model.structure.terminals.size() << "\n";
  oss << "Bifurcacoes (B):            " << model.structure.bifurcations.size() << "\n";
  if (!model.structure.entry_points.empty()) {
    oss << "Ponto de entrada:           " << model.structure.entry_points[0] << "\n";
  } else {
    oss << "Ponto de entrada:           INDEFINIDO\n";
  }

  if (relations.closure_and_order.has_value()) {
    oss << "Matriz da relacao:          " << model.points.size() << " x " << model.points.size()
        << " (matriz densa gerada; " << relations.closure_and_order->count_R
        << " entradas iguais a 1; |R+|=" << relations.closure_and_order->count_R_plus << ")\n";
  } else {
    oss << "Matriz da relacao:          " << model.points.size() << " x " << model.points.size()
        << " (representacao esparsa por lista de pares)\n";
  }

  oss << "Comprimentos 3D:            min=" << format_double(model.structure.min_length)
      << ", medio=" << format_double(model.structure.mean_length)
      << ", max=" << format_double(model.structure.max_length)
      << ", soma=" << format_double(model.structure.total_length) << " (unidades originais VTK)\n";

  oss << "Raios por segmento:         min=" << format_double(model.structure.min_radius)
      << ", medio=" << format_double(model.structure.mean_radius)
      << ", max=" << format_double(model.structure.max_radius) << " (unidades originais VTK)\n";

  oss << "Relacao pai-filho:          " << relations.pairs_R.size()
      << " pares ordenados (|dom(R)|=" << relations.domain_R.size()
      << ", |im(R)|=" << relations.image_R.size() << "; funcao p: V\\{0}->V="
      << (relations.parent_function.is_well_defined_function ? "SIM" : "NAO") << ")\n";

  if (recurrences != nullptr && !recurrences->rows.empty()) {
    oss << "Recorrencia entre modelos:  e_k = 2*e_{k-1} + 1 (e_0=31) => e_k = 32*2^k - 1; "
        << "b_k = 16*2^k - 1; v_k = 32*2^k [validacao nos " << recurrences->rows.size()
        << " modelos: " << (recurrences->all_models_match ? "OK" : "DIVERGENCIA") << "]\n";
  } else {
    oss << "Recorrencia entre modelos:  e_k = 2*e_{k-1} + 1 (e_0=31) => e_k = 32*2^k - 1 "
        << "(|E| = 2*T - 1 = " << (2 * model.structure.terminals.size() - 1) << ")\n";
  }

  oss << "Consulta modular (mod " << b0.modulus << "):   frequencias=[";
  for (std::size_t j = 0; j < b0.classes.size(); ++j) {
    if (j > 0) {
      oss << ", ";
    }
    oss << "C" << j << ":" << b0.classes[j].count;
  }
  oss << "]; soma=" << b0.total_segments << "; PCP max>=ceil(|E|/" << b0.modulus
      << ")=" << b0.pigeonhole_lower_bound << " (obs=" << b0.max_class_size << ")\n";

  oss << "Trilha B1 (TCR):            i===" << crt.solution.a1 << " (mod " << crt.solution.m1
      << ") e i===" << crt.solution.a2 << " (mod " << crt.solution.m2
      << ") <=> i===" << crt.solution.x0 << " (mod " << crt.solution.combined_modulus
      << "); selecionados=" << crt.selected_count
      << ", comp_medio=" << format_double(crt.mean_length)
      << ", raio_medio=" << format_double(crt.mean_radius)
      << " (coincidencia=" << (crt.methods_coincide ? "OK" : "FALHA") << ")\n";

  if (model.validation.ok) {
    oss << "Validacao:                  OK";
    if (const auto ref = get_table1_reference(model.filename); ref.has_value()) {
      const bool matches_t1 = (model.points.size() == ref->points) &&
                              (model.segments.size() == ref->segments) &&
                              (model.structure.terminals.size() == ref->n_term) &&
                              (model.structure.bifurcations.size() == ref->bifurcations);
      oss << " (confere com Tabela 1: " << (matches_t1 ? "SIM" : "NAO") << ")";
    }
    oss << "\n";
  } else {
    oss << "Validacao:                  INCONSISTENCIAS (" << model.validation.errors.size()
        << ")\n";
    for (const auto& err : model.validation.errors) {
      oss << "  - " << err << "\n";
    }
  }

  return oss.str();
}

void export_consolidated_table_csv(const std::filesystem::path& filepath,
                                   const std::vector<VascularModel>& models) {
  std::ostringstream oss;
  oss << "arquivo,N_term,pontos_V,segmentos_E,bifurcacoes_B,ponto_entrada,"
         "comp_min,comp_medio,comp_max,comp_total,raio_min,raio_medio,raio_max,"
         "validacao_ok,confere_tabela1\n";
  oss << std::fixed << std::setprecision(8);
  for (const auto& m : models) {
    bool matches_t1 = false;
    if (const auto ref = get_table1_reference(m.filename); ref.has_value()) {
      matches_t1 = (m.points.size() == ref->points) && (m.segments.size() == ref->segments) &&
                   (m.structure.terminals.size() == ref->n_term) &&
                   (m.structure.bifurcations.size() == ref->bifurcations);
    }
    const std::string entry_str =
        m.structure.entry_points.empty() ? "-1" : std::to_string(m.structure.entry_points[0]);
    oss << m.filename << "," << m.structure.terminals.size() << "," << m.points.size() << ","
        << m.segments.size() << "," << m.structure.bifurcations.size() << "," << entry_str << ","
        << m.structure.min_length << "," << m.structure.mean_length << "," << m.structure.max_length
        << "," << m.structure.total_length << "," << m.structure.min_radius << ","
        << m.structure.mean_radius << "," << m.structure.max_radius << ","
        << (m.validation.ok ? "OK" : "ERRO") << "," << (matches_t1 ? "SIM" : "NAO") << "\n";
  }
  write_text_file(filepath, oss.str());
}

void export_segments_csv(const std::filesystem::path& filepath, const VascularModel& model) {
  std::ostringstream oss;
  oss << "id_segmento,u,v,x_u,y_u,z_u,x_v,y_v,z_v,raio,comprimento_3d,mx,my,mz\n";
  oss << std::fixed << std::setprecision(8);
  for (const auto& seg : model.segments) {
    const Point3D& pu = model.points[seg.u];
    const Point3D& pv = model.points[seg.v];
    oss << seg.id << "," << seg.u << "," << seg.v << "," << pu.x << "," << pu.y << "," << pu.z
        << "," << pv.x << "," << pv.y << "," << pv.z << "," << seg.radius << "," << seg.length
        << "," << seg.midpoint.x << "," << seg.midpoint.y << "," << seg.midpoint.z << "\n";
  }
  write_text_file(filepath, oss.str());
}

void export_modular_b0_csv(const std::filesystem::path& filepath, const ModularAnalysisB0& b0) {
  std::ostringstream oss;
  oss << "classe,residuo_j,cardinalidade_Cj,comprimento_medio,raio_medio\n";
  oss << std::fixed << std::setprecision(8);
  for (const auto& cls : b0.classes) {
    oss << "C_" << cls.residue << "," << cls.residue << "," << cls.count << "," << cls.mean_length
        << "," << cls.mean_radius << "\n";
  }
  write_text_file(filepath, oss.str());
}

void export_crt_b1_csv(const std::filesystem::path& filepath,
                       const std::vector<CrtQueryResult>& results) {
  std::ostringstream oss;
  oss << "arquivo,total_segmentos,congruencia_equivalente,selecionados_direto,"
         "selecionados_equivalente,coincidem,comprimento_medio,raio_medio\n";
  oss << std::fixed << std::setprecision(8);
  for (const auto& r : results) {
    oss << r.filename << "," << r.total_segments << ","
        << "i===" << r.solution.x0 << "(mod " << r.solution.combined_modulus << "),"
        << r.direct_ids.size() << "," << r.equivalent_ids.size() << ","
        << (r.methods_coincide ? "SIM" : "NAO") << "," << r.mean_length << "," << r.mean_radius
        << "\n";
  }
  write_text_file(filepath, oss.str());
}

void export_module_b_report_txt(const std::filesystem::path& filepath, const ModularAnalysisB0& b0,
                                const std::vector<CrtQueryResult>& crt_results) {
  std::ostringstream oss;
  oss << "====================================================================\n";
  oss << "MODULO B -- TEORIA DOS NUMEROS APLICADA A CONSULTA DOS SEGMENTOS\n";
  oss << "====================================================================\n\n";

  oss << "1. PARTE B0 -- NUCLEO OBRIGATORIO (" << b0.filename << ", m = " << b0.modulus << ")\n";
  oss << "Total de segmentos |E| = " << b0.total_segments << "\n\n";
  oss << std::left << std::setw(10) << "Classe" << std::setw(12) << "|C_j|" << std::setw(20)
      << "Comp. Medio" << std::setw(20) << "Raio Medio"
      << "\n";
  oss << std::string(62, '-') << "\n";
  oss << std::fixed << std::setprecision(8);
  for (const auto& cls : b0.classes) {
    oss << std::left << std::setw(10) << ("C_" + std::to_string(cls.residue)) << std::setw(12)
        << cls.count << std::setw(20) << cls.mean_length << std::setw(20) << cls.mean_radius
        << "\n";
  }
  oss << "\nValidacao da particao (B0.2):\n";
  oss << "  - Disjuncao par a par (C_j \\cap C_k = vazio): " << bool_pt(b0.pairwise_disjoint)
      << "\n";
  oss << "  - Soma das cardinalidades == |E|:              " << bool_pt(b0.sum_equals_total)
      << "\n";
  oss << "  - Cobertura exata de {0, ..., |E|-1}:          " << bool_pt(b0.covers_all_indices)
      << "\n";
  oss << "  - Particao valida:                             " << bool_pt(b0.partition_valid) << "\n";

  oss << "\nPrincipio da Casa dos Pombos (B0.3):\n";
  oss << "  - Limite inferior: ceil(" << b0.total_segments << " / " << b0.modulus
      << ") = " << b0.pigeonhole_lower_bound << "\n";
  oss << "  - Maior classe observada: C_" << b0.max_class_residue << " com " << b0.max_class_size
      << " segmentos\n";
  oss << "  - Verificacao (max >= limite inferior): " << bool_pt(b0.pigeonhole_verified) << "\n\n";

  oss << "2. PARTE B1 -- TRILHA ELETIVA ESCOLHIDA: TRILHA 1 (ALGORITMO DE EUCLIDES E TCR)\n";
  if (!crt_results.empty()) {
    const auto& sol = crt_results[0].solution;
    oss << "Sistema de congruencias simultaneas:\n";
    oss << "  i === " << sol.a1 << " (mod " << sol.m1 << ")\n";
    oss << "  i === " << sol.a2 << " (mod " << sol.m2 << ")\n\n";
    oss << "Passos do Algoritmo de Euclides para mdc(" << sol.m1 << ", " << sol.m2 << "):\n";
    for (const auto& step : sol.euclid_m1_m2.steps) {
      oss << "  " << step.a << " = " << step.quotient << " * " << step.b << " + " << step.remainder
          << "\n";
    }
    oss << "  => mdc(" << sol.m1 << ", " << sol.m2 << ") = " << sol.euclid_m1_m2.gcd
        << " (identidade de Bezout: " << sol.m1 << "*(" << sol.euclid_m1_m2.x << ") + " << sol.m2
        << "*(" << sol.euclid_m1_m2.y << ") = " << sol.euclid_m1_m2.gcd << ")\n\n";
    oss << "Resolucao pelo Teorema Chines do Resto:\n";
    oss << "  M  = " << sol.m1 << " * " << sol.m2 << " = " << sol.combined_modulus << "\n";
    oss << "  M1 = M / m1 = " << sol.M1 << ", inverso y1 = M1^{-1} mod " << sol.m1 << " = "
        << sol.y1 << "\n";
    oss << "  M2 = M / m2 = " << sol.M2 << ", inverso y2 = M2^{-1} mod " << sol.m2 << " = "
        << sol.y2 << "\n";
    oss << "  x0 = (" << sol.a1 << "*" << sol.M1 << "*" << sol.y1 << " + " << sol.a2 << "*"
        << sol.M2 << "*" << sol.y2 << ") mod " << sol.combined_modulus << " = " << sol.x0 << "\n";
    oss << "  Congruencia equivalente: i === " << sol.x0 << " (mod " << sol.combined_modulus
        << ")\n\n";
  }

  for (const auto& r : crt_results) {
    oss << "Resultados em " << r.filename << " (|E|=" << r.total_segments << "):\n";
    oss << "  - Selecionados (direto pelas 2 congruencias): " << r.direct_ids.size() << "\n";
    oss << "  - Selecionados (congruencia equivalente TCR): " << r.equivalent_ids.size() << "\n";
    oss << "  - Conjuntos coincidem exatamente:             " << bool_pt(r.methods_coincide)
        << "\n";
    oss << "  - Comprimento medio:                          " << r.mean_length << "\n";
    oss << "  - Raio medio:                                 " << r.mean_radius << "\n";
    oss << "  - IDs selecionados (" << r.selected_count << "): ";
    for (std::size_t idx = 0; idx < r.direct_ids.size(); ++idx) {
      if (idx > 0) {
        oss << ", ";
      }
      oss << r.direct_ids[idx];
    }
    oss << "\n\n";
  }

  oss << "3. CRITERIO DE PERTINENCIA E LIMITACOES:\n";
  oss << "Os identificadores i in {0, ..., |E|-1} sao rotulos computacionais derivados da ordem\n";
  oss << "de armazenamento das entradas LINES no arquivo VTK. Portanto, classes modulares e\n";
  oss << "consultas via TCR servem para particionamento, amostragem deterministica e verificacao\n";
  oss << "reprodutivel de registros, nao implicando propriedades fisiologicas especiais.\n";

  write_text_file(filepath, oss.str());
}

void export_counting_c1_csv(const std::filesystem::path& filepath,
                            const std::vector<ModelCountingRow>& rows) {
  std::ostringstream oss;
  oss << "k,arquivo,T,B,V,E,B_eq_T_minus_1,E_eq_2T_minus_1,V_eq_2T,V_eq_E_plus_1,"
         "V_eq_1_plus_B_plus_T,C_T_2,B_vezes_T\n";
  for (const auto& r : rows) {
    oss << r.k << "," << r.filename << "," << r.T << "," << r.B << "," << r.V << "," << r.E << ","
        << bool_pt(r.rel_B_eq_T_minus_1) << "," << bool_pt(r.rel_E_eq_2T_minus_1) << ","
        << bool_pt(r.rel_V_eq_2T) << "," << bool_pt(r.rel_V_eq_E_plus_1) << ","
        << bool_pt(r.rel_V_eq_1_plus_B_plus_T) << "," << r.pairs_two_terminals << ","
        << r.ordered_pairs_B_T << "\n";
  }
  write_text_file(filepath, oss.str());
}

void export_inclusion_exclusion_csv(const std::filesystem::path& filepath,
                                    const InclusionExclusionResult& ie) {
  std::ostringstream oss;
  oss << "arquivo,total_segmentos,r_med,l_med,card_A,card_B,card_A_inter_B,card_A_union_B,"
         "card_A_menos_B,card_B_menos_A,card_nenhum,identidade_valida\n";
  oss << std::fixed << std::setprecision(8);
  oss << ie.filename << "," << ie.total_segments << "," << ie.median_radius << ","
      << ie.median_length << "," << ie.size_A << "," << ie.size_B << "," << ie.size_A_inter_B << ","
      << ie.size_A_union_B << "," << ie.size_A_minus_B << "," << ie.size_B_minus_A << ","
      << ie.size_neither << "," << bool_pt(ie.identity_holds) << "\n";
  write_text_file(filepath, oss.str());
}

void export_recurrences_c2_csv(const std::filesystem::path& filepath,
                               const RecurrenceAnalysis& rec) {
  std::ostringstream oss;
  oss << "k,N_term,arquivo,V_obs,E_obs,B_obs,T_obs,E_recursivo,E_iteracao,E_caracteristica,"
         "B_recursivo,B_fechada,V_recursivo,V_fechada,confere_tudo\n";
  for (const auto& r : rec.rows) {
    oss << r.k << "," << r.n_term_expected << "," << r.filename << "," << r.observed_V << ","
        << r.observed_E << "," << r.observed_B << "," << r.observed_T << "," << r.recursive_E << ","
        << r.iteration_E << "," << r.characteristic_E << "," << r.recursive_B << "," << r.closed_B
        << "," << r.recursive_V << "," << r.closed_V << "," << bool_pt(r.matches_all) << "\n";
  }
  write_text_file(filepath, oss.str());
}

void export_module_c_report_txt(const std::filesystem::path& filepath,
                                const std::vector<ModelCountingRow>& counting,
                                const InclusionExclusionResult& ie, const RecurrenceAnalysis& rec) {
  std::ostringstream oss;
  oss << "====================================================================\n";
  oss << "MODULO C -- CONTAGEM E RELACOES DE RECORRENCIA\n";
  oss << "====================================================================\n\n";

  oss << "1. PARTE C1 -- TECNICAS DE CONTAGEM\n";
  oss << std::left << std::setw(22) << "Arquivo" << std::setw(6) << "T" << std::setw(6) << "B"
      << std::setw(6) << "|V|" << std::setw(6) << "|E|" << std::setw(12) << "C(T,2)"
      << std::setw(12) << "B * T"
      << "\n";
  oss << std::string(70, '-') << "\n";
  for (const auto& r : counting) {
    oss << std::left << std::setw(22) << r.filename << std::setw(6) << r.T << std::setw(6) << r.B
        << std::setw(6) << r.V << std::setw(6) << r.E << std::setw(12) << r.pairs_two_terminals
        << std::setw(12) << r.ordered_pairs_B_T << "\n";
  }
  oss << "\nRelacoes estruturais verificadas em todos os modelos:\n";
  oss << "  B = T - 1, |E| = 2T - 1, |V| = 2T, |V| = |E| + 1, |V| = 1 + B + T\n";
  oss << "  Dupla contagem de graus: soma(out_degree) = 1 + 2B = |E| e soma(in_degree) = B + T = "
         "|E|\n\n";

  oss << "Inclusao-Exclusao 3D em " << ie.filename << " (C1.4):\n";
  oss << std::fixed << std::setprecision(8);
  oss << "  r_med (mediana dos raios):        " << ie.median_radius << "\n";
  oss << "  l_med (mediana dos comprimentos): " << ie.median_length << "\n";
  oss << "  |A| = |{i : r_i >= r_med}| =      " << ie.size_A << "\n";
  oss << "  |B| = |{i : l_i >= l_med}| =      " << ie.size_B << "\n";
  oss << "  |A \\cap B| =                      " << ie.size_A_inter_B << "\n";
  oss << "  |A \\cup B| =                      " << ie.size_A_union_B << "\n";
  oss << "  |A| + |B| - |A \\cap B| =          " << (ie.size_A + ie.size_B - ie.size_A_inter_B)
      << "\n";
  oss << "  Identidade |A U B| = |A| + |B| - |A n B| fechada: " << bool_pt(ie.identity_holds)
      << "\n";
  oss << "  Complementares: |A\\B|=" << ie.size_A_minus_B << ", |B\\A|=" << ie.size_B_minus_A
      << ", |(A U B)^c|=" << ie.size_neither << "\n\n";

  oss << "2. PARTE C2 -- RECORRENCIAS E VALIDACAO\n";
  oss << "Recorrencia principal (segmentos): e_0 = 31, e_k = 2*e_{k-1} + 1 (k >= 1)\n";
  oss << "  Solucao por iteracao:            e_k = 2^k * e_0 + (2^k - 1) = 32 * 2^k - 1\n";
  oss << "  Solucao por equacao caract.:     e_k = 32 * 2^k - 1 = 2^{k+5} - 1\n";
  oss << "Segundas relacoes (C2.4):\n";
  oss << "  Bifurcacoes: b_0 = 15, b_k = 2*b_{k-1} + 1 => b_k = 16 * 2^k - 1 = T_k - 1\n";
  oss << "  Pontos:      v_0 = 32, v_k = 2*v_{k-1}     => v_k = 32 * 2^k     = 2 * T_k\n\n";
  oss << std::left << std::setw(4) << "k" << std::setw(8) << "Nterm" << std::setw(8) << "|V|"
      << std::setw(8) << "e_k" << std::setw(8) << "B" << std::setw(8) << "T" << std::setw(10)
      << "32*2^k-1" << std::setw(10) << "Confere"
      << "\n";
  oss << std::string(64, '-') << "\n";
  for (const auto& r : rec.rows) {
    oss << std::left << std::setw(4) << r.k << std::setw(8) << r.n_term_expected << std::setw(8)
        << r.observed_V << std::setw(8) << r.observed_E << std::setw(8) << r.observed_B
        << std::setw(8) << r.observed_T << std::setw(10) << r.characteristic_E << std::setw(10)
        << bool_pt(r.matches_all) << "\n";
  }
  oss << "\nValidacao global das recorrencias nos 6 modelos: " << bool_pt(rec.all_models_match)
      << "\n";

  write_text_file(filepath, oss.str());
}

void export_relation_pairs_csv(const std::filesystem::path& filepath,
                               const std::vector<DirectedPair>& pairs) {
  std::ostringstream oss;
  oss << "u,v\n";
  for (const auto& p : pairs) {
    oss << p.u << "," << p.v << "\n";
  }
  write_text_file(filepath, oss.str());
}

void export_bool_matrix_csv(const std::filesystem::path& filepath, const BoolMatrix& matrix) {
  std::ostringstream oss;
  const std::size_t n = matrix.size();
  oss << "u\\v";
  for (std::size_t v = 0; v < n; ++v) {
    oss << "," << v;
  }
  oss << "\n";
  for (std::size_t u = 0; u < n; ++u) {
    oss << u;
    for (std::size_t v = 0; v < n; ++v) {
      oss << "," << static_cast<int>(matrix[u][v]);
    }
    oss << "\n";
  }
  write_text_file(filepath, oss.str());
}

void export_module_d_report_txt(const std::filesystem::path& filepath,
                                const RelationAnalysisD& rel_16,
                                const ModularEquivalenceAnalysis& eq_128) {
  std::ostringstream oss;
  oss << "====================================================================\n";
  oss << "MODULO D -- RELACOES, REPRESENTACOES E PROPRIEDADES\n";
  oss << "====================================================================\n\n";

  oss << "1. RELACAO PAI-FILHO R EM " << rel_16.filename << " (D1-D2)\n";
  oss << "  |V| = " << rel_16.num_vertices << ", |R| = " << rel_16.pairs_R.size() << "\n";
  oss << "  |dom(R)| = " << rel_16.domain_R.size() << ", dom(R) = {";
  for (std::size_t i = 0; i < rel_16.domain_R.size(); ++i) {
    if (i > 0) {
      oss << ", ";
    }
    oss << rel_16.domain_R[i];
  }
  oss << "}\n";
  oss << "  |im(R)|  = " << rel_16.image_R.size() << " (igual a V \\ {0})\n";
  oss << "  Primeiros 8 pares de R: ";
  for (std::size_t i = 0; i < std::min<std::size_t>(8, rel_16.pairs_R.size()); ++i) {
    if (i > 0) {
      oss << ", ";
    }
    oss << "(" << rel_16.pairs_R[i].u << "," << rel_16.pairs_R[i].v << ")";
  }
  oss << "\n\n";

  oss << "2. PROPRIEDADES DE R E R^{-1} (D3)\n";
  oss << "  Relacao pai-filho R:\n";
  oss << "    - Reflexiva:     " << bool_pt(rel_16.props_R.reflexive);
  if (rel_16.props_R.non_reflexive_example.has_value()) {
    const std::size_t v = *rel_16.props_R.non_reflexive_example;
    oss << " (contraexemplo: (" << v << "," << v << ") nao pertence a R)";
  }
  oss << "\n";
  oss << "    - Irreflexiva:   " << bool_pt(rel_16.props_R.irreflexive) << "\n";
  oss << "    - Simetrica:     " << bool_pt(rel_16.props_R.symmetric);
  if (rel_16.props_R.non_symmetric_example.has_value()) {
    const auto p = *rel_16.props_R.non_symmetric_example;
    oss << " (contraexemplo: (" << p.u << "," << p.v << ") in R, mas (" << p.v << "," << p.u
        << ") not in R)";
  }
  oss << "\n";
  oss << "    - Antissimetrica:" << " " << bool_pt(rel_16.props_R.antisymmetric)
      << " (assimetrica: " << bool_pt(rel_16.props_R.asymmetric) << ")\n";
  oss << "    - Transitiva:    " << bool_pt(rel_16.props_R.transitive);
  if (rel_16.props_R.non_transitive_example.has_value()) {
    const auto w = *rel_16.props_R.non_transitive_example;
    oss << " (contraexemplo: (" << w.u << "," << w.v << ") in R e (" << w.v << "," << w.w
        << ") in R, mas (" << w.u << "," << w.w << ") not in R)";
  }
  oss << "\n";
  oss << "  Relacao inversa R^{-1} (filho-pai):\n";
  oss << "    - |R^{-1}| = " << rel_16.pairs_R_inv.size()
      << ", |dom(R^{-1})| = " << rel_16.domain_R_inv.size()
      << ", |im(R^{-1})| = " << rel_16.image_R_inv.size() << "\n";
  oss << "    - Reflexiva: " << bool_pt(rel_16.props_R_inv.reflexive)
      << ", Simetrica: " << bool_pt(rel_16.props_R_inv.symmetric)
      << ", Antissimetrica: " << bool_pt(rel_16.props_R_inv.antisymmetric)
      << ", Transitiva: " << bool_pt(rel_16.props_R_inv.transitive)
      << ", Funcional em V\\{0}: " << bool_pt(rel_16.R_inv_is_functional) << "\n\n";

  oss << "3. FUNCAO ASSOCIADA p: V\\{0} -> V (D4)\n";
  oss << "  - Ponto 0 sem pai (grau de entrada 0):                "
      << bool_pt(rel_16.parent_function.root_zero_has_no_parent) << "\n";
  oss << "  - Todo v in V\\{0} possui pai unico (grau entrada 1): "
      << bool_pt(rel_16.parent_function.all_non_root_have_unique_parent) << "\n";
  oss << "  - p e funcao bem definida:                            "
      << bool_pt(rel_16.parent_function.is_well_defined_function) << "\n";
  oss << "  - p e injetora:                                       "
      << bool_pt(rel_16.parent_function.is_injective);
  if (rel_16.parent_function.non_injective_witness.has_value()) {
    const auto [v1, v2, par] = *rel_16.parent_function.non_injective_witness;
    oss << " (contraexemplo: p(" << v1 << ") = p(" << v2 << ") = " << par << ")";
  }
  oss << "\n";
  oss << "  - p e sobrejetora em V:                               "
      << bool_pt(rel_16.parent_function.is_surjective_onto_V);
  if (rel_16.parent_function.non_surjective_witness.has_value()) {
    oss << " (contraexemplo: terminal " << *rel_16.parent_function.non_surjective_witness
        << " nao e imagem de nenhum ponto)";
  }
  oss << "\n\n";

  if (rel_16.closure_and_order.has_value()) {
    const auto& cl = *rel_16.closure_and_order;
    oss << "4. FECHO TRANSITIVO R^+ E ORDEM PARCIAL <= EM " << rel_16.filename << " (D5)\n";
    oss << "  - Dimensao da matriz booleana: " << cl.num_vertices << " x " << cl.num_vertices
        << "\n";
    oss << "  - |R| = " << cl.count_R << ", |R^+| = " << cl.count_R_plus
        << ", pares indiretos |R^+ \\ R| = " << cl.count_indirect << ", |<=| = " << cl.count_leq
        << "\n";
    oss << "  - Iteracoes do fecho iterativo ate ponto fixo: " << cl.closure_iterations << " (";
    for (std::size_t i = 0; i < cl.pairs_added_per_iteration.size(); ++i) {
      if (i > 0) {
        oss << ", ";
      }
      oss << "it" << (i + 1) << ":+" << cl.pairs_added_per_iteration[i];
    }
    oss << ")\n";
    oss << "  - Exemplos verificados de ancestralidade indireta (u -> mid -> w):\n";
    for (std::size_t i = 0; i < std::min<std::size_t>(6, cl.indirect_witnesses.size()); ++i) {
      const auto& w = cl.indirect_witnesses[i];
      oss << "      (" << w.u << "," << w.v << ") e (" << w.v << "," << w.w << ") => (" << w.u
          << "," << w.w << ") in R^+ \\ R\n";
    }
    oss << "  - Axiomas da relacao <= = R^+ U Delta_V:\n";
    oss << "      Reflexiva: " << bool_pt(cl.leq_reflexive)
        << ", Antissimetrica: " << bool_pt(cl.leq_antisymmetric)
        << ", Transitiva: " << bool_pt(cl.leq_transitive) << "\n";
    oss << "      E ordem parcial: " << bool_pt(cl.is_partial_order)
        << "; E ordem total: " << bool_pt(cl.is_total_order) << "\n";
    oss << "      Ponto 0 e elemento minimo (0 <= v para todo v): "
        << bool_pt(cl.root_zero_is_minimum) << "\n";
    oss << "  - Pares nao ordenados incomparaveis: " << cl.incomparable_pair_count << "\n";
    oss << "    Exemplos: ";
    for (std::size_t i = 0; i < std::min<std::size_t>(6, cl.incomparable_examples.size()); ++i) {
      if (i > 0) {
        oss << ", ";
      }
      oss << "{" << cl.incomparable_examples[i].u << ", " << cl.incomparable_examples[i].v << "}";
    }
    oss << "\n\n";
  }

  oss << "5. RELACAO DE EQUIVALENCIA MODULAR EM " << eq_128.filename << " (D6)\n";
  oss << "  Conjunto S = {0, ..., " << (eq_128.set_size > 0 ? eq_128.set_size - 1 : 0)
      << "}, i ~ j <=> i === j (mod " << eq_128.modulus << ")\n";
  oss << "  - Reflexiva: " << bool_pt(eq_128.reflexive)
      << ", Simetrica: " << bool_pt(eq_128.symmetric)
      << ", Transitiva: " << bool_pt(eq_128.transitive) << "\n";
  oss << "  - E relacao de equivalencia: " << bool_pt(eq_128.is_equivalence_relation) << "\n";
  oss << "  - Classes [0]..[6] coincidem exatamente com C_0..C_6 do Modulo B: "
      << bool_pt(eq_128.matches_module_b_classes) << "\n";

  write_text_file(filepath, oss.str());
}

} // namespace cco

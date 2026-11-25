#include "cco/runner.hpp"

#include "cco/combinatorics.hpp"
#include "cco/exporter.hpp"
#include "cco/number_theory.hpp"
#include "cco/relations.hpp"
#include "cco/vtk_parser.hpp"

#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace cco {

namespace {

const VascularModel* find_model_by_name(const std::vector<VascularModel>& models,
                                        const std::string& name) {
  for (const auto& m : models) {
    if (m.filename == name) {
      return &m;
    }
  }
  return nullptr;
}

std::string strip_extension(const std::string& filename) {
  const auto pos = filename.rfind('.');
  return pos == std::string::npos ? filename : filename.substr(0, pos);
}

void print_help() {
  std::cout
      << "CCO Explorer -- Analise de Redes Vasculares por Matematica Discreta (BCC102/UFOP)\n\n"
      << "Uso:\n"
      << "  cco_explorer [opcoes]\n\n"
      << "Opcoes:\n"
      << "  --entrada <dir>      Diretorio contendo os arquivos .vtk (padrao: dados)\n"
      << "  --modelo <arquivo>   Exibe o Resumo Matematico de um modelo VTK especifico\n"
      << "                       (ex: tree3D_Nterm128.vtk)\n"
      << "  --trilha <nome>      Trilha eletiva da Parte B1 (padrao e implementada: crt)\n"
      << "  --saida <dir>        Exporta todos os arquivos CSV e TXT para o diretorio\n"
      << "  --todos              Executa analise integrada de todos os modelos VTK no diretorio\n"
      << "  --ajuda, -h          Exibe esta mensagem de ajuda\n\n"
      << "Exemplos:\n"
      << "  cco_explorer --entrada dados/ --modelo tree3D_Nterm128.vtk --trilha crt\n"
      << "  cco_explorer --entrada dados/ --todos --trilha crt --saida resultados/\n";
}

} // namespace

int execute_cli(const CliOptions& options) {
  if (options.show_help) {
    print_help();
    return 0;
  }

  if (options.track != "crt") {
    std::cerr << "Erro: trilha '" << options.track
              << "' nao suportada. A trilha eletiva implementada e 'crt' "
                 "(Trilha 1: Algoritmo de Euclides e Teorema Chines do Resto).\n";
    return 2;
  }

  try {
    // Mode 1: Single model requested (--modelo) without --todos
    if (options.model_filename.has_value() && !options.run_all) {
      std::filesystem::path target_path(*options.model_filename);
      if (!std::filesystem::exists(target_path)) {
        target_path = options.input_dir / *options.model_filename;
      }
      if (!std::filesystem::exists(target_path)) {
        std::cerr << "Erro: modelo nao encontrado em '" << target_path.string() << "'.\n";
        return 1;
      }

      const VascularModel model = parse_vtk_file(target_path);
      const ModularAnalysisB0 b0 = analyze_modular_b0(model, 7);
      const CrtQueryResult crt = run_crt_query(model, 1, 4, 2, 5);
      const bool dense_matrix = (model.points.size() <= 32);
      const RelationAnalysisD rel = analyze_relations_d(model, dense_matrix);

      // If input_dir contains the full 6-model collection, also verify recurrences
      std::optional<RecurrenceAnalysis> rec_opt;
      if (std::filesystem::exists(options.input_dir) &&
          std::filesystem::is_directory(options.input_dir)) {
        const auto all_models = load_all_vtk_models(options.input_dir);
        if (all_models.size() == 6) {
          rec_opt = analyze_recurrences_c2(all_models);
        }
      }

      const std::string summary =
          format_mathematical_summary(model, b0, crt, rel, rec_opt ? &(*rec_opt) : nullptr);
      std::cout << summary;

      if (options.output_dir.has_value()) {
        const auto& out_dir = *options.output_dir;
        const std::string stem = strip_extension(model.filename);
        write_text_file(out_dir / ("resumo_matematico_" + stem + ".txt"), summary);
        export_segments_csv(out_dir / ("segmentos_" + stem + ".csv"), model);
        export_modular_b0_csv(out_dir / ("classes_mod7_" + stem + ".csv"), b0);
        export_relation_pairs_csv(out_dir / ("pares_R_" + stem + ".csv"), rel.pairs_R);
        if (rel.closure_and_order.has_value()) {
          export_bool_matrix_csv(out_dir / ("matriz_R_" + stem + ".csv"),
                                 rel.closure_and_order->matrix_R);
          export_bool_matrix_csv(out_dir / ("matriz_Rplus_" + stem + ".csv"),
                                 rel.closure_and_order->matrix_R_plus);
        }
      }
      return model.validation.ok ? 0 : 1;
    }

    // Mode 2: Full collection analysis (all .vtk files in input_dir)
    const std::vector<VascularModel> models = load_all_vtk_models(options.input_dir);
    if (models.empty()) {
      std::cerr << "Erro: nenhum arquivo .vtk encontrado em '" << options.input_dir.string()
                << "'.\n";
      return 1;
    }

    const RecurrenceAnalysis recurrences = analyze_recurrences_c2(models);
    const std::vector<ModelCountingRow> counting = analyze_counting_c1(models);

    const VascularModel* m016 = find_model_by_name(models, "tree3D_Nterm016.vtk");
    const VascularModel* m128 = find_model_by_name(models, "tree3D_Nterm128.vtk");
    const VascularModel* m512 = find_model_by_name(models, "tree3D_Nterm512.vtk");

    if (m016 == nullptr) {
      m016 = &models.front();
    }
    if (m128 == nullptr) {
      m128 = &models[models.size() / 2];
    }
    if (m512 == nullptr) {
      m512 = &models.back();
    }

    // Module B
    const ModularAnalysisB0 b0_128 = analyze_modular_b0(*m128, 7);
    const CrtQueryResult crt_128 = run_crt_query(*m128, 1, 4, 2, 5);
    const CrtQueryResult crt_512 = run_crt_query(*m512, 1, 4, 2, 5);
    const std::vector<CrtQueryResult> crt_results = {crt_128, crt_512};

    // Module C (Inclusion-Exclusion on Nterm128)
    const InclusionExclusionResult ie_128 = analyze_inclusion_exclusion(*m128);

    // Module D (dense matrix & closure only on Nterm016; modular equivalence on Nterm128)
    const RelationAnalysisD rel_16 = analyze_relations_d(*m016, true);
    const RelationAnalysisD rel_128 = analyze_relations_d(*m128, false);
    const ModularEquivalenceAnalysis eq_128 = analyze_modular_equivalence_d6(*m128, b0_128);

    std::ostringstream all_summaries;
    for (std::size_t i = 0; i < models.size(); ++i) {
      const auto& m = models[i];
      const ModularAnalysisB0 b0 = analyze_modular_b0(m, 7);
      const CrtQueryResult crt = run_crt_query(m, 1, 4, 2, 5);
      const bool dense = (m.points.size() <= 32);
      const RelationAnalysisD rel = analyze_relations_d(m, dense);
      if (i > 0) {
        all_summaries << "\n";
      }
      all_summaries << format_mathematical_summary(m, b0, crt, rel, &recurrences);
    }

    std::cout << all_summaries.str();

    if (options.output_dir.has_value()) {
      const auto& out_dir = *options.output_dir;
      std::filesystem::create_directories(out_dir);

      // Module A outputs
      export_consolidated_table_csv(out_dir / "resumo_consolidado.csv", models);
      for (const auto& m : models) {
        const std::string tag = strip_extension(m.filename);
        export_segments_csv(out_dir / "segmentos" / ("segmentos_" + tag + ".csv"), m);
      }

      // Module B outputs
      export_modular_b0_csv(out_dir / "modulo_b_classes_mod7_Nterm128.csv", b0_128);
      export_crt_b1_csv(out_dir / "modulo_b_trilha1_crt.csv", crt_results);
      export_module_b_report_txt(out_dir / "modulo_b_relatorio.txt", b0_128, crt_results);

      // Module C outputs
      export_counting_c1_csv(out_dir / "modulo_c_contagem.csv", counting);
      export_inclusion_exclusion_csv(out_dir / "modulo_c_inclusao_exclusao_Nterm128.csv", ie_128);
      export_recurrences_c2_csv(out_dir / "modulo_c_recorrencias.csv", recurrences);
      export_module_c_report_txt(out_dir / "modulo_c_relatorio.txt", counting, ie_128, recurrences);

      // Module D outputs
      export_relation_pairs_csv(out_dir / "modulo_d_pares_R_Nterm016.csv", rel_16.pairs_R);
      export_relation_pairs_csv(out_dir / "modulo_d_pares_R_Nterm128.csv", rel_128.pairs_R);
      if (rel_16.closure_and_order.has_value()) {
        export_bool_matrix_csv(out_dir / "modulo_d_matriz_R_Nterm016.csv",
                               rel_16.closure_and_order->matrix_R);
        export_bool_matrix_csv(out_dir / "modulo_d_matriz_Rplus_Nterm016.csv",
                               rel_16.closure_and_order->matrix_R_plus);
        export_relation_pairs_csv(out_dir / "modulo_d_pares_Rplus_Nterm016.csv",
                                  rel_16.closure_and_order->pairs_R_plus);
      }
      export_module_d_report_txt(out_dir / "modulo_d_relatorio.txt", rel_16, eq_128);

      // Module E outputs (Section 8.2 summaries)
      write_text_file(out_dir / "resumo_completo_todos_modelos.txt", all_summaries.str());
      const ModularAnalysisB0 b0_16 = analyze_modular_b0(*m016, 7);
      const CrtQueryResult crt_16 = run_crt_query(*m016, 1, 4, 2, 5);
      write_text_file(out_dir / "resumo_matematico_Nterm016.txt",
                      format_mathematical_summary(*m016, b0_16, crt_16, rel_16, &recurrences));
      write_text_file(out_dir / "resumo_matematico_Nterm128.txt",
                      format_mathematical_summary(*m128, b0_128, crt_128, rel_128, &recurrences));

      std::cout << "\nResultados exportados com sucesso para: " << out_dir.string() << "\n";
    }

    bool all_ok = true;
    for (const auto& m : models) {
      if (!m.validation.ok) {
        all_ok = false;
      }
    }
    return all_ok ? 0 : 1;
  } catch (const std::exception& ex) {
    std::cerr << "Erro durante a execucao: " << ex.what() << "\n";
    return 1;
  }
}

} // namespace cco

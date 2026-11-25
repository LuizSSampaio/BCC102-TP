#pragma once

#include "cco/combinatorics.hpp"
#include "cco/model.hpp"
#include "cco/number_theory.hpp"
#include "cco/relations.hpp"

#include <filesystem>
#include <string>
#include <vector>

namespace cco {

// Formats the Section 8.2 "CCO EXPLORER -- RESUMO MATEMATICO" for a single model.
std::string format_mathematical_summary(const VascularModel& model, const ModularAnalysisB0& b0,
                                        const CrtQueryResult& crt,
                                        const RelationAnalysisD& relations,
                                        const RecurrenceAnalysis* recurrences = nullptr);

// Module A exports
void export_consolidated_table_csv(const std::filesystem::path& filepath,
                                   const std::vector<VascularModel>& models);
void export_segments_csv(const std::filesystem::path& filepath, const VascularModel& model);

// Module B exports
void export_modular_b0_csv(const std::filesystem::path& filepath, const ModularAnalysisB0& b0);
void export_crt_b1_csv(const std::filesystem::path& filepath,
                       const std::vector<CrtQueryResult>& results);
void export_module_b_report_txt(const std::filesystem::path& filepath, const ModularAnalysisB0& b0,
                                const std::vector<CrtQueryResult>& crt_results);

// Module C exports
void export_counting_c1_csv(const std::filesystem::path& filepath,
                            const std::vector<ModelCountingRow>& rows);
void export_inclusion_exclusion_csv(const std::filesystem::path& filepath,
                                    const InclusionExclusionResult& ie);
void export_recurrences_c2_csv(const std::filesystem::path& filepath,
                               const RecurrenceAnalysis& rec);
void export_module_c_report_txt(const std::filesystem::path& filepath,
                                const std::vector<ModelCountingRow>& counting,
                                const InclusionExclusionResult& ie, const RecurrenceAnalysis& rec);

// Module D exports
void export_relation_pairs_csv(const std::filesystem::path& filepath,
                               const std::vector<DirectedPair>& pairs);
void export_bool_matrix_csv(const std::filesystem::path& filepath, const BoolMatrix& matrix);
void export_module_d_report_txt(const std::filesystem::path& filepath,
                                const RelationAnalysisD& rel_16,
                                const ModularEquivalenceAnalysis& eq_128);

// Helper to write any text file (creates parent directories automatically).
void write_text_file(const std::filesystem::path& filepath, const std::string& content);

} // namespace cco

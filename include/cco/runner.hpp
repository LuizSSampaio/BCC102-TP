#pragma once

#include <filesystem>
#include <optional>
#include <string>

namespace cco {

struct CliOptions {
  std::filesystem::path input_dir = "dados";
  std::optional<std::string> model_filename;
  std::string track = "crt";
  std::optional<std::filesystem::path> output_dir;
  bool run_all = false;
  bool show_help = false;
};

// Executes the CCO Explorer workflow according to CLI options.
// Returns 0 on success, non-zero on error.
int execute_cli(const CliOptions& options);

} // namespace cco

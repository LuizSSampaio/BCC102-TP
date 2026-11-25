#include "cco/runner.hpp"

#include <iostream>
#include <string>

int main(int argc, char* argv[]) {
  cco::CliOptions opts;

  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "--ajuda" || arg == "-h" || arg == "--help") {
      opts.show_help = true;
    } else if (arg == "--todos" || arg == "--all") {
      opts.run_all = true;
    } else if ((arg == "--entrada" || arg == "--input") && i + 1 < argc) {
      opts.input_dir = argv[++i];
    } else if ((arg == "--modelo" || arg == "--model") && i + 1 < argc) {
      opts.model_filename = argv[++i];
    } else if ((arg == "--trilha" || arg == "--track") && i + 1 < argc) {
      opts.track = argv[++i];
    } else if ((arg == "--saida" || arg == "--output") && i + 1 < argc) {
      opts.output_dir = argv[++i];
    } else {
      std::cerr << "Argumento desconhecido ou sem valor: " << arg << "\n";
      std::cerr << "Use --ajuda para listar as opcoes validas.\n";
      return 2;
    }
  }

  return cco::execute_cli(opts);
}

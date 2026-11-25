{ pkgs, lib, config, inputs, ... }:

{
  # https://devenv.sh/packages/
  packages = [
    pkgs.git
    pkgs.clang
    pkgs.clang-tools
    pkgs.cmake
  ];

  # https://devenv.sh/languages/
  languages.cplusplus = {
    enable = true;
    lsp.package = pkgs.clang;
  };

  languages.python.enable = true;

  languages.texlive.enable = true;

  # Set compiler environment variables to Clang
  env.CXX = "clang++";
  env.CC = "clang";

  # https://devenv.sh/scripts/
  scripts = {
    build = {
      description = "Build C++ targets (cco_explorer and cco_tests) using CMake";
      exec = ''
        cmake -B build -S .
        cmake --build build "$@"
        if [ -f build/compile_commands.json ] && [ ! -e compile_commands.json ]; then
          ln -sf build/compile_commands.json compile_commands.json
        fi
      '';
    };

    "run-tests" = {
      description = "Run C++ acceptance and unit tests via CTest";
      exec = ''
        [ -f build/CMakeCache.txt ] || cmake -B build -S .
        cmake --build build
        ctest --test-dir build --output-on-failure "$@"
      '';
    };

    run-all = {
      description = "Build and run CCO Explorer on all six VTK models, exporting outputs to resultados/";
      exec = ''
        [ -f build/CMakeCache.txt ] || cmake -B build -S .
        cmake --build build --target cco_explorer
        ./build/cco_explorer --entrada dados --todos --trilha crt --saida resultados "$@"
      '';
    };

    format = {
      description = "Format C++ source and header files using clang-format via CMake";
      exec = ''
        [ -f build/CMakeCache.txt ] || cmake -B build -S .
        cmake --build build --target format "$@"
      '';
    };

    "format-check" = {
      description = "Check formatting of C++ source and header files via CMake";
      exec = ''
        [ -f build/CMakeCache.txt ] || cmake -B build -S .
        cmake --build build --target format-check "$@"
      '';
    };

    tidy = {
      description = "Run clang-tidy on C++ source files via CMake";
      exec = ''
        [ -f build/CMakeCache.txt ] || cmake -B build -S .
        cmake --build build --target tidy "$@"
      '';
    };

    lint = {
      description = "Run linter (alias for tidy)";
      exec = ''
        [ -f build/CMakeCache.txt ] || cmake -B build -S .
        cmake --build build --target lint "$@"
      '';
    };

    clean = {
      description = "Remove build artifacts using CMake";
      exec = ''
        if [ "$1" = "--all" ]; then
          rm -rf build
        elif [ -f build/CMakeCache.txt ]; then
          cmake --build build --target clean "$@"
        else
          rm -rf build
        fi
      '';
    };
  };

  # See full reference at https://devenv.sh/reference/options/
}

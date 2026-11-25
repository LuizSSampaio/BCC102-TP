cmake_minimum_required(VERSION 3.20)
project(cco_explorer VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CCO_CORE_HEADERS
  include/cco/model.hpp
  include/cco/vtk_parser.hpp
  include/cco/number_theory.hpp
  include/cco/combinatorics.hpp
  include/cco/relations.hpp
  include/cco/exporter.hpp
  include/cco/runner.hpp
)

set(CCO_CORE_SOURCES
  src/vtk_parser.cpp
  src/number_theory.cpp
  src/combinatorics.cpp
  src/relations.cpp
  src/exporter.cpp
  src/runner.cpp
)

add_library(cco_core STATIC ${CCO_CORE_SOURCES} ${CCO_CORE_HEADERS})
target_include_directories(cco_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(cco_core PRIVATE
  -Wall
  -Wextra
  -Wpedantic
  -Wshadow
  -Wconversion
  -Wsign-conversion
)

add_executable(cco_explorer src/main.cpp)
target_link_libraries(cco_explorer PRIVATE cco_core)
target_compile_options(cco_explorer PRIVATE
  -Wall
  -Wextra
  -Wpedantic
  -Wshadow
  -Wconversion
  -Wsign-conversion
)

enable_testing()
add_executable(cco_tests tests/test_cco_explorer.cpp)
target_link_libraries(cco_tests PRIVATE cco_core)
target_compile_options(cco_tests PRIVATE
  -Wall
  -Wextra
  -Wpedantic
  -Wshadow
  -Wconversion
  -Wsign-conversion
)
add_test(
  NAME cco_acceptance_and_unit_tests
  COMMAND cco_tests --entrada ${CMAKE_CURRENT_SOURCE_DIR}/dados
)

# Tooling targets integrated with devenv.nix
set(CCO_ALL_CXX_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cco_explorer.cpp
)
foreach(f ${CCO_CORE_HEADERS} ${CCO_CORE_SOURCES})
  list(APPEND CCO_ALL_CXX_FILES ${CMAKE_CURRENT_SOURCE_DIR}/${f})
endforeach()

set(CCO_TIDY_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cco_explorer.cpp
)
foreach(f ${CCO_CORE_SOURCES})
  list(APPEND CCO_TIDY_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/${f})
endforeach()

find_program(CLANG_FORMAT_EXE NAMES clang-format)
if(CLANG_FORMAT_EXE)
  add_custom_target(format
    COMMAND ${CLANG_FORMAT_EXE} -i -style=file ${CCO_ALL_CXX_FILES}
    COMMENT "Formatando arquivos C++ com clang-format"
    VERBATIM
  )
  add_custom_target(format-check
    COMMAND ${CLANG_FORMAT_EXE} --dry-run --Werror -style=file ${CCO_ALL_CXX_FILES}
    COMMENT "Verificando formatacao C++ com clang-format"
    VERBATIM
  )
endif()

find_program(CLANG_TIDY_EXE NAMES clang-tidy)
if(CLANG_TIDY_EXE)
  add_custom_target(tidy
    COMMAND ${CLANG_TIDY_EXE} -p ${CMAKE_BINARY_DIR} ${CCO_TIDY_SOURCES}
    COMMENT "Executando analise estatica com clang-tidy"
    VERBATIM
  )
  add_custom_target(lint
    DEPENDS tidy
    COMMENT "Executando linter (alias de tidy)"
  )
endif()

#pragma once

#include "cco/model.hpp"

#include <cstddef>
#include <filesystem>
#include <optional>
#include <string>
#include <vector>

namespace cco {

struct Table1Reference {
  std::string filename;
  std::size_t n_term = 0;
  std::size_t points = 0;
  std::size_t segments = 0;
  std::size_t bifurcations = 0;
};

// Module A (A1-A4): reads an ASCII POLYDATA VTK file, computes 3D lengths and
// midpoints, validates integrity, and classifies points structurally.
VascularModel parse_vtk_file(const std::filesystem::path& filepath);

// Loads all .vtk files found in dir_path, sorted by terminal count T ascending.
std::vector<VascularModel> load_all_vtk_models(const std::filesystem::path& dir_path);

// Returns the expected Table 1 reference row if filename is one of the six
// provided CCO models.
std::optional<Table1Reference> get_table1_reference(const std::string& filename);

// Returns all six reference rows from Table 1 of the specification.
std::vector<Table1Reference> all_table1_references();

} // namespace cco

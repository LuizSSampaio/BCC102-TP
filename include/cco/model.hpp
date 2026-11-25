#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace cco {

struct Point3D {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
};

struct Segment {
  std::size_t id = 0;
  std::size_t u = 0;
  std::size_t v = 0;
  double radius = 0.0;
  double length = 0.0;
  Point3D midpoint{};
};

struct ValidationReport {
  bool ok = false;
  std::size_t declared_points = 0;
  std::size_t read_points = 0;
  std::size_t declared_lines = 0;
  std::size_t read_lines = 0;
  std::size_t declared_line_ints = 0;
  std::size_t expected_line_ints = 0;
  std::size_t declared_cell_data = 0;
  std::size_t read_radii = 0;
  bool header_valid = false;
  bool all_lines_have_two_vertices = false;
  bool all_endpoints_valid = false;
  bool no_self_loops = false;
  bool one_radius_per_segment = false;
  bool all_radii_positive_finite = false;
  bool topology_consistent = false;
  std::vector<std::string> errors;
};

struct StructuralAnalysis {
  std::vector<std::size_t> in_degree;
  std::vector<std::size_t> out_degree;
  std::vector<std::size_t> entry_points;
  std::vector<std::size_t> terminals;
  std::vector<std::size_t> bifurcations;
  std::vector<std::size_t> other_points;
  double min_length = 0.0;
  double max_length = 0.0;
  double mean_length = 0.0;
  double total_length = 0.0;
  double min_radius = 0.0;
  double max_radius = 0.0;
  double mean_radius = 0.0;
};

struct VascularModel {
  std::string filename;
  std::string filepath;
  std::vector<Point3D> points;
  std::vector<Segment> segments;
  ValidationReport validation;
  StructuralAnalysis structure;
};

} // namespace cco

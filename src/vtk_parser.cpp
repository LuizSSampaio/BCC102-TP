#include "cco/vtk_parser.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cco {

namespace {

std::string to_upper_ascii(std::string s) {
  for (char& ch : s) {
    ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
  }
  return s;
}

std::string trim(const std::string& s) {
  const auto first = s.find_first_not_of(" \t\r\n");
  if (first == std::string::npos) {
    return "";
  }
  const auto last = s.find_last_not_of(" \t\r\n");
  return s.substr(first, last - first + 1);
}

bool read_non_empty_line(std::istream& in, std::string& out_line) {
  while (std::getline(in, out_line)) {
    out_line = trim(out_line);
    if (!out_line.empty()) {
      return true;
    }
  }
  return false;
}

} // namespace

std::vector<Table1Reference> all_table1_references() {
  return {
      {"tree3D_Nterm016.vtk", 16, 32, 31, 15},     {"tree3D_Nterm032.vtk", 32, 64, 63, 31},
      {"tree3D_Nterm064.vtk", 64, 128, 127, 63},   {"tree3D_Nterm128.vtk", 128, 256, 255, 127},
      {"tree3D_Nterm256.vtk", 256, 512, 511, 255}, {"tree3D_Nterm512.vtk", 512, 1024, 1023, 511},
  };
}

std::optional<Table1Reference> get_table1_reference(const std::string& filename) {
  for (const auto& ref : all_table1_references()) {
    if (ref.filename == filename) {
      return ref;
    }
  }
  return std::nullopt;
}

VascularModel parse_vtk_file(const std::filesystem::path& filepath) {
  VascularModel model;
  model.filename = filepath.filename().string();
  model.filepath = filepath.string();

  std::ifstream file(filepath);
  if (!file.is_open()) {
    model.validation.errors.push_back("Nao foi possivel abrir o arquivo VTK: " + filepath.string());
    model.validation.ok = false;
    return model;
  }

  // A1: Read VTK ASCII / DATASET POLYDATA header
  std::string line1;
  std::string line2;
  std::string line3;
  std::string line4;

  const bool got_header = static_cast<bool>(std::getline(file, line1)) &&
                          static_cast<bool>(std::getline(file, line2)) &&
                          static_cast<bool>(std::getline(file, line3)) &&
                          static_cast<bool>(std::getline(file, line4));

  if (!got_header) {
    model.validation.errors.push_back("Cabecalho VTK incompleto.");
    model.validation.ok = false;
    return model;
  }

  line1 = trim(line1);
  line3 = to_upper_ascii(trim(line3));
  line4 = to_upper_ascii(trim(line4));

  bool header_ok = true;
  if (!line1.starts_with("# vtk DataFile")) {
    model.validation.errors.push_back("Linha 1 nao inicia com '# vtk DataFile'.");
    header_ok = false;
  }
  if (line3 != "ASCII") {
    model.validation.errors.push_back("Formato VTK nao e ASCII: " + line3);
    header_ok = false;
  }
  {
    std::istringstream iss(line4);
    std::string kw;
    std::string dtype;
    iss >> kw >> dtype;
    if (kw != "DATASET" || dtype != "POLYDATA") {
      model.validation.errors.push_back("Tipo de dataset nao e POLYDATA: " + line4);
      header_ok = false;
    }
  }
  model.validation.header_valid = header_ok;

  // Section parsing: POINTS, LINES, CELL_DATA
  bool found_points = false;
  bool found_lines = false;
  bool found_cell_data = false;

  std::vector<double> radii;
  bool all_lines_two = true;
  bool all_endpoints_valid = true;
  bool no_self_loops = true;
  bool all_radii_pos_finite = true;

  std::string line;
  while (read_non_empty_line(file, line)) {
    std::istringstream iss(line);
    std::string keyword;
    iss >> keyword;
    const std::string upper_kw = to_upper_ascii(keyword);

    if (upper_kw == "POINTS") {
      found_points = true;
      std::size_t n_points = 0;
      std::string data_type;
      if (!(iss >> n_points >> data_type)) {
        model.validation.errors.push_back("Cabecalho POINTS malformado: " + line);
        break;
      }
      model.validation.declared_points = n_points;
      model.points.reserve(n_points);

      for (std::size_t i = 0; i < n_points; ++i) {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        if (!(file >> x >> y >> z)) {
          model.validation.errors.push_back("Falha ao ler coordenadas do ponto " +
                                            std::to_string(i) + ".");
          break;
        }
        if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z)) {
          model.validation.errors.push_back("Coordenada nao finita no ponto " + std::to_string(i) +
                                            ".");
        }
        model.points.push_back(Point3D{x, y, z});
      }
      model.validation.read_points = model.points.size();
      std::getline(file, line); // consume trailing newline
    } else if (upper_kw == "LINES") {
      found_lines = true;
      std::size_t n_lines = 0;
      std::size_t total_ints = 0;
      if (!(iss >> n_lines >> total_ints)) {
        model.validation.errors.push_back("Cabecalho LINES malformado: " + line);
        break;
      }
      model.validation.declared_lines = n_lines;
      model.validation.declared_line_ints = total_ints;
      model.validation.expected_line_ints = 3 * n_lines;
      model.segments.reserve(n_lines);

      for (std::size_t i = 0; i < n_lines; ++i) {
        std::int64_t vertex_count = -1;
        std::int64_t u_raw = -1;
        std::int64_t v_raw = -1;
        if (!(file >> vertex_count)) {
          model.validation.errors.push_back("Falha ao ler segmento " + std::to_string(i) +
                                            " em LINES.");
          break;
        }
        if (vertex_count != 2) {
          all_lines_two = false;
          model.validation.errors.push_back("Segmento " + std::to_string(i) + " possui " +
                                            std::to_string(vertex_count) +
                                            " vertices (esperado: 2).");
          // Skip remaining indices of this malformed line
          for (std::int64_t skip = 0; skip < vertex_count; ++skip) {
            std::int64_t dummy = 0;
            file >> dummy;
          }
          continue;
        }
        if (!(file >> u_raw >> v_raw)) {
          model.validation.errors.push_back("Falha ao ler extremos do segmento " +
                                            std::to_string(i) + ".");
          break;
        }
        if (u_raw < 0 || v_raw < 0 || static_cast<std::size_t>(u_raw) >= model.points.size() ||
            static_cast<std::size_t>(v_raw) >= model.points.size()) {
          all_endpoints_valid = false;
          model.validation.errors.push_back(
              "Segmento " + std::to_string(i) + " referencia ponto inexistente: (" +
              std::to_string(u_raw) + ", " + std::to_string(v_raw) + ").");
        }
        if (u_raw == v_raw) {
          no_self_loops = false;
          model.validation.errors.push_back("Segmento " + std::to_string(i) +
                                            " e um laco proprio no ponto " + std::to_string(u_raw) +
                                            ".");
        }

        Segment seg;
        seg.id = i;
        seg.u = u_raw >= 0 ? static_cast<std::size_t>(u_raw) : 0;
        seg.v = v_raw >= 0 ? static_cast<std::size_t>(v_raw) : 0;
        model.segments.push_back(seg);
      }
      model.validation.read_lines = model.segments.size();
      std::getline(file, line); // consume trailing newline
    } else if (upper_kw == "CELL_DATA") {
      found_cell_data = true;
      std::size_t n_cells = 0;
      if (!(iss >> n_cells)) {
        model.validation.errors.push_back("Cabecalho CELL_DATA malformado: " + line);
        break;
      }
      model.validation.declared_cell_data = n_cells;

      std::string scalars_line;
      std::string lookup_line;
      if (!read_non_empty_line(file, scalars_line) || !read_non_empty_line(file, lookup_line)) {
        model.validation.errors.push_back("Metadados SCALARS/LOOKUP_TABLE ausentes em CELL_DATA.");
        break;
      }
      {
        std::istringstream s_iss(scalars_line);
        std::string s_kw;
        std::string attr_name;
        s_iss >> s_kw >> attr_name;
        if (to_upper_ascii(s_kw) != "SCALARS" || attr_name != "raio") {
          model.validation.errors.push_back("Esperado 'scalars raio ...', encontrado: " +
                                            scalars_line);
        }
      }
      {
        std::istringstream l_iss(lookup_line);
        std::string l_kw;
        l_iss >> l_kw;
        if (to_upper_ascii(l_kw) != "LOOKUP_TABLE") {
          model.validation.errors.push_back("Esperado 'LOOKUP_TABLE ...', encontrado: " +
                                            lookup_line);
        }
      }

      radii.reserve(n_cells);
      for (std::size_t i = 0; i < n_cells; ++i) {
        double r = 0.0;
        if (!(file >> r)) {
          model.validation.errors.push_back("Falha ao ler raio do segmento " + std::to_string(i) +
                                            ".");
          break;
        }
        if (!std::isfinite(r) || r <= 0.0) {
          all_radii_pos_finite = false;
          model.validation.errors.push_back("Raio invalido (<= 0 ou nao finito) no segmento " +
                                            std::to_string(i) + ".");
        }
        radii.push_back(r);
      }
      model.validation.read_radii = radii.size();

      // Check for unexpected extra numeric tokens right after radii
      double extra_val = 0.0;
      if (file >> extra_val) {
        model.validation.errors.push_back(
            "Existem valores extras de raio alem do declarado em CELL_DATA.");
      }
      break;
    }
  }

  if (!found_points) {
    model.validation.errors.push_back("Secao POINTS nao encontrada.");
  }
  if (!found_lines) {
    model.validation.errors.push_back("Secao LINES nao encontrada.");
  }
  if (!found_cell_data) {
    model.validation.errors.push_back("Secao CELL_DATA nao encontrada.");
  }

  if (model.validation.declared_points != model.validation.read_points) {
    model.validation.errors.push_back("Numero de pontos lidos diverge do declarado em POINTS.");
  }
  if (model.validation.declared_lines != model.validation.read_lines) {
    model.validation.errors.push_back("Numero de segmentos lidos diverge do declarado em LINES.");
  }
  if (model.validation.declared_line_ints != model.validation.expected_line_ints) {
    model.validation.errors.push_back("Total de inteiros declarado em LINES diverge de 3 * |E|.");
  }

  const bool one_radius =
      (model.validation.declared_cell_data == model.validation.declared_lines) &&
      (model.validation.read_radii == model.segments.size());
  if (!one_radius) {
    model.validation.errors.push_back(
        "Quantidade de valores de raio nao coincide com o numero de segmentos.");
  }

  model.validation.all_lines_have_two_vertices = all_lines_two;
  model.validation.all_endpoints_valid = all_endpoints_valid;
  model.validation.no_self_loops = no_self_loops;
  model.validation.one_radius_per_segment = one_radius;
  model.validation.all_radii_positive_finite = all_radii_pos_finite;

  // A2: Assign radii and compute 3D geometric length and midpoint for each segment
  const std::size_t n_seg = model.segments.size();
  double sum_len = 0.0;
  double sum_rad = 0.0;
  double min_len = std::numeric_limits<double>::infinity();
  double max_len = -std::numeric_limits<double>::infinity();
  double min_rad = std::numeric_limits<double>::infinity();
  double max_rad = -std::numeric_limits<double>::infinity();

  for (std::size_t i = 0; i < n_seg; ++i) {
    Segment& seg = model.segments[i];
    if (i < radii.size()) {
      seg.radius = radii[i];
    }
    if (seg.u < model.points.size() && seg.v < model.points.size()) {
      const Point3D& pu = model.points[seg.u];
      const Point3D& pv = model.points[seg.v];
      const double dx = pv.x - pu.x;
      const double dy = pv.y - pu.y;
      const double dz = pv.z - pu.z;
      seg.length = std::sqrt(dx * dx + dy * dy + dz * dz);
      seg.midpoint = Point3D{(pu.x + pv.x) / 2.0, (pu.y + pv.y) / 2.0, (pu.z + pv.z) / 2.0};
    }

    sum_len += seg.length;
    sum_rad += seg.radius;
    min_len = std::min(min_len, seg.length);
    max_len = std::max(max_len, seg.length);
    min_rad = std::min(min_rad, seg.radius);
    max_rad = std::max(max_rad, seg.radius);
  }

  if (n_seg > 0) {
    model.structure.min_length = min_len;
    model.structure.max_length = max_len;
    model.structure.mean_length = sum_len / static_cast<double>(n_seg);
    model.structure.total_length = sum_len;
    model.structure.min_radius = min_rad;
    model.structure.max_radius = max_rad;
    model.structure.mean_radius = sum_rad / static_cast<double>(n_seg);
  }

  // A4: Structural identification of points via simple loops and counters
  const std::size_t n_pts = model.points.size();
  model.structure.in_degree.assign(n_pts, 0);
  model.structure.out_degree.assign(n_pts, 0);

  if (all_endpoints_valid) {
    for (const auto& seg : model.segments) {
      model.structure.out_degree[seg.u]++;
      model.structure.in_degree[seg.v]++;
    }
  }

  bool topo_ok = true;
  for (std::size_t v = 0; v < n_pts; ++v) {
    const std::size_t in_d = model.structure.in_degree[v];
    const std::size_t out_d = model.structure.out_degree[v];

    if (in_d == 0 && out_d == 1) {
      model.structure.entry_points.push_back(v);
    } else if (in_d == 1 && out_d == 0) {
      model.structure.terminals.push_back(v);
    } else if (in_d == 1 && out_d == 2) {
      model.structure.bifurcations.push_back(v);
    } else {
      model.structure.other_points.push_back(v);
      topo_ok = false;
      model.validation.errors.push_back(
          "Ponto " + std::to_string(v) + " possui graus fora do padrao CCO (in=" +
          std::to_string(in_d) + ", out=" + std::to_string(out_d) + ").");
    }
  }

  if (model.structure.entry_points.size() != 1 || model.structure.entry_points[0] != 0) {
    topo_ok = false;
    model.validation.errors.push_back(
        "Topologia esperada requer ponto de entrada unico com identificador 0 (in=0, out=1).");
  }

  model.validation.topology_consistent = topo_ok;
  model.validation.ok = model.validation.errors.empty();
  return model;
}

std::vector<VascularModel> load_all_vtk_models(const std::filesystem::path& dir_path) {
  if (!std::filesystem::exists(dir_path) || !std::filesystem::is_directory(dir_path)) {
    throw std::runtime_error("Diretorio de entrada invalido: " + dir_path.string());
  }

  std::vector<std::filesystem::path> vtk_files;
  for (const auto& entry : std::filesystem::directory_iterator(dir_path)) {
    if (entry.is_regular_file() && entry.path().extension() == ".vtk") {
      vtk_files.push_back(entry.path());
    }
  }

  std::ranges::sort(vtk_files);

  std::vector<VascularModel> models;
  models.reserve(vtk_files.size());
  for (const auto& path : vtk_files) {
    models.push_back(parse_vtk_file(path));
  }

  // Order by terminal count T ascending (16, 32, 64, 128, 256, 512)
  std::ranges::sort(models, [](const VascularModel& a, const VascularModel& b) {
    if (a.structure.terminals.size() != b.structure.terminals.size()) {
      return a.structure.terminals.size() < b.structure.terminals.size();
    }
    return a.filename < b.filename;
  });

  return models;
}

} // namespace cco

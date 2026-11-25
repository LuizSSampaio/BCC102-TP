#!/usr/bin/env python3
"""
Script auxiliar de verificacao cruzada independente para o CCO Explorer (BCC102/UFOP).

Le os seis arquivos VTK diretamente em Python puro (biblioteca padrao), recalcula
de forma independente todas as metricas dos Modulos A, B, C e D e confronta os
resultados com os arquivos CSV gerados pelo executavel C++ `cco_explorer`.
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import sys
from pathlib import Path


TABLE1_EXPECTED = {
    "tree3D_Nterm016.vtk": (16, 32, 31, 15),
    "tree3D_Nterm032.vtk": (32, 64, 63, 31),
    "tree3D_Nterm064.vtk": (64, 128, 127, 63),
    "tree3D_Nterm128.vtk": (128, 256, 255, 127),
    "tree3D_Nterm256.vtk": (256, 512, 511, 255),
    "tree3D_Nterm512.vtk": (512, 1024, 1023, 511),
}


def parse_vtk(path: Path) -> dict:
    lines = [ln.strip() for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    assert lines[0].startswith("# vtk DataFile"), f"Cabecalho invalido em {path}"
    assert lines[2].upper() == "ASCII", f"Nao ASCII em {path}"
    assert lines[3].upper() == "DATASET POLYDATA", f"Nao POLYDATA em {path}"

    idx = 4
    pts_hdr = lines[idx].split()
    assert pts_hdr[0].upper() == "POINTS"
    n_pts = int(pts_hdr[1])
    idx += 1

    points = []
    for _ in range(n_pts):
        x, y, z = map(float, lines[idx].split())
        points.append((x, y, z))
        idx += 1

    lines_hdr = lines[idx].split()
    assert lines_hdr[0].upper() == "LINES"
    n_seg = int(lines_hdr[1])
    n_ints = int(lines_hdr[2])
    assert n_ints == 3 * n_seg
    idx += 1

    edges = []
    for _ in range(n_seg):
        c, u, v = map(int, lines[idx].split())
        assert c == 2
        assert 0 <= u < n_pts and 0 <= v < n_pts and u != v
        edges.append((u, v))
        idx += 1

    cd_hdr = lines[idx].split()
    assert cd_hdr[0].upper() == "CELL_DATA"
    n_cd = int(cd_hdr[1])
    assert n_cd == n_seg
    idx += 1
    assert lines[idx].lower().startswith("scalars raio")
    idx += 1
    assert lines[idx].upper().startswith("LOOKUP_TABLE")
    idx += 1

    radii = [float(lines[idx + i]) for i in range(n_cd)]
    assert idx + n_cd == len(lines), f"Sobras inesperadas em {path}"

    segments = []
    in_deg = [0] * n_pts
    out_deg = [0] * n_pts
    for i, ((u, v), r) in enumerate(zip(edges, radii)):
        pu, pv = points[u], points[v]
        length = math.sqrt((pv[0] - pu[0]) ** 2 + (pv[1] - pu[1]) ** 2 + (pv[2] - pu[2]) ** 2)
        mid = ((pu[0] + pv[0]) / 2.0, (pu[1] + pv[1]) / 2.0, (pu[2] + pv[2]) / 2.0)
        segments.append({"id": i, "u": u, "v": v, "r": r, "l": length, "mid": mid})
        out_deg[u] += 1
        in_deg[v] += 1

    entries = [v for v in range(n_pts) if in_deg[v] == 0 and out_deg[v] == 1]
    terminals = [v for v in range(n_pts) if in_deg[v] == 1 and out_deg[v] == 0]
    bifurcations = [v for v in range(n_pts) if in_deg[v] == 1 and out_deg[v] == 2]
    assert entries == [0]
    assert len(entries) + len(terminals) + len(bifurcations) == n_pts

    return {
        "filename": path.name,
        "points": points,
        "segments": segments,
        "T": len(terminals),
        "B": len(bifurcations),
        "V": n_pts,
        "E": n_seg,
    }


def approx_eq(a: float, b: float, tol: float = 1e-6) -> bool:
    return abs(a - b) <= tol


def main() -> int:
    parser = argparse.ArgumentParser(description="Verificacao cruzada dos resultados do CCO Explorer")
    parser.add_argument("--entrada", type=Path, default=Path("dados"))
    parser.add_argument("--resultados", type=Path, default=Path("resultados"))
    args = parser.parse_args()

    vtk_paths = sorted(args.entrada.glob("*.vtk"))
    models = sorted((parse_vtk(p) for p in vtk_paths), key=lambda m: m["T"])
    assert len(models) == 6, f"Esperados 6 modelos VTK, encontrados {len(models)}"
    by_name = {m["filename"]: m for m in models}

    # T1 & T2 + Modulo A CSV check
    for m in models:
        exp = TABLE1_EXPECTED[m["filename"]]
        assert (m["T"], m["V"], m["E"], m["B"]) == exp
    with (args.resultados / "resumo_consolidado.csv").open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 6
    for row, m in zip(rows, models):
        assert row["arquivo"] == m["filename"]
        assert int(row["N_term"]) == m["T"]
        assert int(row["pontos_V"]) == m["V"]
        assert int(row["segmentos_E"]) == m["E"]
        assert int(row["bifurcacoes_B"]) == m["B"]
        mean_l = statistics.fmean(s["l"] for s in m["segments"])
        mean_r = statistics.fmean(s["r"] for s in m["segments"])
        assert approx_eq(float(row["comp_medio"]), mean_l)
        assert approx_eq(float(row["raio_medio"]), mean_r)
    print("[OK] T1/T2 e Modulo A (resumo_consolidado.csv) conferidos com Python independente.")

    # Modulo B (B0 mod 7 em Nterm128 + B1 CRT em Nterm128 e Nterm512)
    m128 = by_name["tree3D_Nterm128.vtk"]
    m512 = by_name["tree3D_Nterm512.vtk"]
    with (args.resultados / "modulo_b_classes_mod7_Nterm128.csv").open(encoding="utf-8") as f:
        b0_rows = list(csv.DictReader(f))
    assert len(b0_rows) == 7
    total_b0 = 0
    for j, row in enumerate(b0_rows):
        cls = [s for s in m128["segments"] if s["id"] % 7 == j]
        assert int(row["residuo_j"]) == j
        assert int(row["cardinalidade_Cj"]) == len(cls)
        assert approx_eq(float(row["comprimento_medio"]), statistics.fmean(s["l"] for s in cls))
        assert approx_eq(float(row["raio_medio"]), statistics.fmean(s["r"] for s in cls))
        total_b0 += len(cls)
    assert total_b0 == m128["E"] == 255

    with (args.resultados / "modulo_b_trilha1_crt.csv").open(encoding="utf-8") as f:
        crt_rows = list(csv.DictReader(f))
    assert len(crt_rows) == 2
    for row, m in zip(crt_rows, [m128, m512]):
        sel_direct = [s for s in m["segments"] if s["id"] % 4 == 1 and s["id"] % 5 == 2]
        sel_equiv = [s for s in m["segments"] if s["id"] % 20 == 17]
        assert sel_direct == sel_equiv
        assert int(row["selecionados_direto"]) == len(sel_direct)
        assert int(row["selecionados_equivalente"]) == len(sel_equiv)
        assert row["coincidem"] == "SIM"
        assert approx_eq(float(row["comprimento_medio"]), statistics.fmean(s["l"] for s in sel_direct))
        assert approx_eq(float(row["raio_medio"]), statistics.fmean(s["r"] for s in sel_direct))
    print("[OK] T7 e Modulo B (B0 mod 7 e B1 Trilha 1 TCR) conferidos com Python independente.")

    # Modulo C (C1 contagem, C1.4 inclusao-exclusao, C2 recorrencias)
    r_med = statistics.median(s["r"] for s in m128["segments"])
    l_med = statistics.median(s["l"] for s in m128["segments"])
    set_A = {s["id"] for s in m128["segments"] if s["r"] >= r_med}
    set_B = {s["id"] for s in m128["segments"] if s["l"] >= l_med}
    with (args.resultados / "modulo_c_inclusao_exclusao_Nterm128.csv").open(encoding="utf-8") as f:
        ie_row = next(csv.DictReader(f))
    assert approx_eq(float(ie_row["r_med"]), r_med)
    assert approx_eq(float(ie_row["l_med"]), l_med)
    assert int(ie_row["card_A"]) == len(set_A)
    assert int(ie_row["card_B"]) == len(set_B)
    assert int(ie_row["card_A_inter_B"]) == len(set_A & set_B)
    assert int(ie_row["card_A_union_B"]) == len(set_A | set_B)
    assert len(set_A | set_B) == len(set_A) + len(set_B) - len(set_A & set_B)
    assert ie_row["identidade_valida"] == "SIM"

    with (args.resultados / "modulo_c_recorrencias.csv").open(encoding="utf-8") as f:
        rec_rows = list(csv.DictReader(f))
    assert len(rec_rows) == 6
    for k, (row, m) in enumerate(zip(rec_rows, models)):
        assert int(row["k"]) == k
        assert int(row["E_obs"]) == 32 * (2**k) - 1 == m["E"]
        assert int(row["B_obs"]) == 16 * (2**k) - 1 == m["B"]
        assert int(row["V_obs"]) == 32 * (2**k) == m["V"]
        assert row["confere_tudo"] == "SIM"
    print("[OK] T5/T6 e Modulo C (contagem, inclusao-exclusao, recorrencias) conferidos.")

    # Modulo D (matriz 32x32 e R^+ de Nterm016)
    m016 = by_name["tree3D_Nterm016.vtk"]
    n16 = m016["V"]
    mat_R = [[0] * n16 for _ in range(n16)]
    for s in m016["segments"]:
        mat_R[s["u"]][s["v"]] = 1

    with (args.resultados / "modulo_d_matriz_R_Nterm016.csv").open(encoding="utf-8") as f:
        csv_rows = list(csv.reader(f))
    assert len(csv_rows) == 33 and len(csv_rows[0]) == 33
    ones_R = 0
    for u in range(n16):
        parsed = [int(x) for x in csv_rows[u + 1][1:]]
        assert parsed == mat_R[u]
        ones_R += sum(parsed)
    assert ones_R == 31

    # Independent Warshall for R^+
    mat_Rp = [row[:] for row in mat_R]
    for k in range(n16):
        for i in range(n16):
            if mat_Rp[i][k]:
                for j in range(n16):
                    if mat_Rp[k][j]:
                        mat_Rp[i][j] = 1
    with (args.resultados / "modulo_d_matriz_Rplus_Nterm016.csv").open(encoding="utf-8") as f:
        csv_rp_rows = list(csv.reader(f))
    assert len(csv_rp_rows) == 33
    for u in range(n16):
        parsed = [int(x) for x in csv_rp_rows[u + 1][1:]]
        assert parsed == mat_Rp[u]
    print("[OK] T3/T4 e Modulo D (matriz 32x32 de R e fecho transitivo R^+) conferidos.")

    print("\nVERIFICACAO CRUZADA PYTHON <-> C++ CONCLUIDA COM 100% DE CONFORMIDADE!")
    return 0


if __name__ == "__main__":
    sys.exit(main())

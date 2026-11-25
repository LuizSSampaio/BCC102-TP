#!/usr/bin/env python3
"""
Gera figuras vetoriais em TikZ (LaTeX) diretamente a partir dos arquivos CSV
produzidos pelo executavel C++ `cco_explorer`.

Nao requer bibliotecas Python externas (usa apenas biblioteca padrao).
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def generate_recurrence_figure(res_dir: Path, out_path: Path) -> None:
    with (res_dir / "modulo_c_recorrencias.csv").open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    # Plot log2 scale on Y so k=0..5 (31..1024) are all clearly legible,
    # alongside exact numeric labels at every point.
    width = 11.2
    height = 5.4
    y_min_log = 3.5   # 2^3.5 ~ 11.3
    y_max_log = 10.6  # 2^10.6 ~ 1552

    def sx(k: float) -> float:
        return 0.7 + (k / 5.0) * (width - 1.4)

    def sy(val: float) -> float:
        lv = math.log2(max(val, 1.0))
        return 0.5 + ((lv - y_min_log) / (y_max_log - y_min_log)) * (height - 0.9)

    lines = [
        "% Figura gerada automaticamente por scripts/gerar_figuras_tikz.py",
        r"\begin{tikzpicture}[font=\footnotesize]",
        f"  % Eixos e grade horizontal (potencias de 2)",
        f"  \\draw[gray!25] (0.5,0.5) rectangle ({width:.2f},{height:.2f});",
    ]

    for p in range(4, 11):
        val = 2**p
        y = sy(val)
        lines.append(f"  \\draw[gray!25, dashed] (0.5,{y:.3f}) -- ({width:.2f},{y:.3f});")
        lines.append(f"  \\node[anchor=east, font=\\scriptsize] at (0.45,{y:.3f}) {{{val}}};")

    for r in rows:
        k = int(r["k"])
        nterm = int(r["N_term"])
        x = sx(k)
        lines.append(f"  \\draw[gray!20] ({x:.3f},0.5) -- ({x:.3f},{height:.2f});")
        lines.append(f"  \\node[anchor=north, font=\\scriptsize, align=center] at ({x:.3f},0.42) "
                     f"{{$k={k}$\\\\({nterm})}};")

    # Axis labels
    lines.append(f"  \\node[font=\\footnotesize, anchor=north] at ({width/2:.2f},-0.35) "
                 r"{Índice da coleção $k$ e quantidade de terminais $N_{\mathrm{term}} = 16 \cdot 2^k$};")
    lines.append(f"  \\node[font=\\footnotesize, rotate=90, anchor=south] at (-0.65,{height/2:.2f}) "
                 r"{Contagem de elementos (escala $\log_2$)};")

    # Series: |V|, |E| (= e_k), B
    series_spec = [
        ("V_obs", "teal!75!black", "solid", r"$|V| = v_k = 32\cdot 2^k$", "above", 0.14),
        ("E_obs", "blue!75!black", "dashed", r"$|E| = e_k = 32\cdot 2^k - 1$", "below right", -0.05),
        ("B_obs", "orange!85!black", "solid", r"$B = b_k = 16\cdot 2^k - 1$", "below", -0.14),
    ]

    for col, color, style, _label, pos, dy in series_spec:
        coords = " -- ".join(f"({sx(int(r['k'])):.3f},{sy(float(r[col])):.3f})" for r in rows)
        lines.append(f"  \\draw[{color}, thick, {style}] {coords};")
        for r in rows:
            x = sx(int(r["k"]))
            v = int(r[col])
            y = sy(v)
            lines.append(f"  \\fill[{color}] ({x:.3f},{y:.3f}) circle (2.0pt);")
            if col in ("V_obs", "B_obs"):
                lines.append(f"  \\node[{color}, font=\\tiny, {pos}] at ({x:.3f},{y + dy:.3f}) {{{v}}};")

    # Legend box
    lx, ly = 0.85, height - 0.35
    lines.append(f"  \\fill[white, draw=gray!50, rounded corners=2pt] ({lx:.2f},{ly-1.35:.2f}) "
                 f"rectangle ({lx+4.35:.2f},{ly+0.15:.2f});")
    for idx, (_col, color, style, label, _pos, _dy) in enumerate(series_spec):
        yy = ly - 0.18 - idx * 0.44
        lines.append(f"  \\draw[{color}, thick, {style}] ({lx+0.15:.2f},{yy:.2f}) -- ({lx+0.65:.2f},{yy:.2f});")
        lines.append(f"  \\fill[{color}] ({lx+0.40:.2f},{yy:.2f}) circle (1.8pt);")
        lines.append(f"  \\node[anchor=west, font=\\scriptsize] at ({lx+0.72:.2f},{yy:.2f}) {{{label}}};")

    lines.append(r"\end{tikzpicture}")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_modular_b0_figure(res_dir: Path, out_path: Path) -> None:
    with (res_dir / "modulo_b_classes_mod7_Nterm128.csv").open(encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    # Two side-by-side panels: (a) cardinalities |C_j| with PCP=37 line; (b) mean radius per class
    lines = [
        "% Figura gerada automaticamente por scripts/gerar_figuras_tikz.py",
        r"\begin{tikzpicture}[font=\footnotesize]",
        "  % Painel esquerdo: Cardinalidades |C_j|",
        r"  \begin{scope}[xshift=0cm]",
        r"    \draw[->] (0,0) -- (6.2,0);",
        r"    \draw[->] (0,0) -- (0,4.2);",
        r"    \node[font=\footnotesize\bfseries, anchor=south] at (3.1,4.3) {(a) Frequências $|C_j|$ ($m=7$, $|E|=255$)};",
        r"    \node[font=\scriptsize, rotate=90, anchor=south] at (-0.65,2.1) {Segmentos $|C_j|$};",
    ]

    max_y = 42.0
    scale_y = 3.8 / max_y
    for ytick in (0, 10, 20, 30, 37, 40):
        yy = ytick * scale_y
        lines.append(f"    \\draw[gray!25] (0,{yy:.3f}) -- (6.0,{yy:.3f});")
        lines.append(f"    \\node[anchor=east, font=\\tiny] at (-0.06,{yy:.3f}) {{{ytick}}};")

    for r in rows:
        j = int(r["residuo_j"])
        cnt = int(r["cardinalidade_Cj"])
        x0 = 0.35 + j * 0.80
        h = cnt * scale_y
        fill_col = "teal!65!white" if cnt == 37 else "teal!35!white"
        lines.append(f"    \\fill[{fill_col}, draw=teal!80!black] ({x0:.2f},0) rectangle ({x0+0.55:.2f},{h:.3f});")
        lines.append(f"    \\node[anchor=south, font=\\tiny] at ({x0+0.275:.2f},{h+0.03:.3f}) {{{cnt}}};")
        lines.append(f"    \\node[anchor=north, font=\\scriptsize] at ({x0+0.275:.2f},-0.05) {{$C_{j}$}};")

    pcp_y = 37 * scale_y
    lines.append(f"    \\draw[red!75!black, dashed, thick] (0,{pcp_y:.3f}) -- (6.0,{pcp_y:.3f});")
    lines.append(f"    \\node[red!75!black, anchor=south west, font=\\tiny] at (3.2,{pcp_y+0.12:.3f}) "
                 r"{PCP: $\lceil 255/7 \rceil = 37$};")
    lines.append(r"  \end{scope}")

    # Right panel: Mean radius per class (showing stability around global mean ~0.3534)
    lines.extend([
        "  % Painel direito: Raio medio por classe",
        r"  \begin{scope}[xshift=7.4cm]",
        r"    \draw[->] (0,0) -- (6.2,0);",
        r"    \draw[->] (0,0) -- (0,4.2);",
        r"    \node[font=\footnotesize\bfseries, anchor=south] at (3.1,4.3) {(b) Raio médio por classe (unid. VTK)};",
        r"    \node[font=\scriptsize, rotate=90, anchor=south] at (-0.75,2.1) {Raio médio $\bar{r}_j$};",
    ])

    r_min, r_max = 0.30, 0.39
    def sy_r(val: float) -> float:
        return ((val - r_min) / (r_max - r_min)) * 3.8

    for rtick in (0.30, 0.32, 0.34, 0.36, 0.38):
        yy = sy_r(rtick)
        lines.append(f"    \\draw[gray!25] (0,{yy:.3f}) -- (6.0,{yy:.3f});")
        lines.append(f"    \\node[anchor=east, font=\\tiny] at (-0.06,{yy:.3f}) {{{rtick:.2f}}};")

    radii_vals = [float(r["raio_medio"]) for r in rows]
    global_mean_r = sum(float(r["raio_medio"]) * int(r["cardinalidade_Cj"]) for r in rows) / 255.0

    for r in rows:
        j = int(r["residuo_j"])
        rm = float(r["raio_medio"])
        x0 = 0.35 + j * 0.80
        h = sy_r(rm)
        lines.append(f"    \\fill[blue!45!white, draw=blue!75!black] ({x0:.2f},0) rectangle ({x0+0.55:.2f},{h:.3f});")
        lines.append(f"    \\node[anchor=south, font=\\tiny] at ({x0+0.275:.2f},{h+0.03:.3f}) {{{rm:.3f}}};")
        lines.append(f"    \\node[anchor=north, font=\\scriptsize] at ({x0+0.275:.2f},-0.05) {{$C_{j}$}};")

    gm_y = sy_r(global_mean_r)
    lines.append(f"    \\draw[orange!85!black, dashed, thick] (0,{gm_y:.3f}) -- (6.0,{gm_y:.3f});")
    lines.append(f"    \\node[orange!85!black, fill=white, fill opacity=0.9, text opacity=1, "
                 f"inner sep=1.5pt, anchor=north east, font=\\tiny] at (5.95,3.88) "
                 f"{{Linha tracejada: média global = {global_mean_r:.4f}}};")
    lines.append(r"  \end{scope}")
    lines.append(r"\end{tikzpicture}")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_inclusion_exclusion_figure(res_dir: Path, out_path: Path) -> None:
    with (res_dir / "modulo_c_inclusao_exclusao_Nterm128.csv").open(encoding="utf-8") as f:
        ie = next(csv.DictReader(f))
    r_med = float(ie["r_med"])
    l_med = float(ie["l_med"])

    with (res_dir / "segmentos" / "segmentos_tree3D_Nterm128.csv").open(encoding="utf-8") as f:
        segs = list(csv.DictReader(f))

    width = 11.4
    height = 6.2
    r_min, r_max = 0.10, 1.30
    l_min, l_max = 0.0, 0.0155

    def sx(r: float) -> float:
        return ((r - r_min) / (r_max - r_min)) * width

    def sy(l: float) -> float:
        return ((l - l_min) / (l_max - l_min)) * height

    xr_med = sx(r_med)
    yl_med = sy(l_med)

    lines = [
        "% Figura gerada automaticamente por scripts/gerar_figuras_tikz.py",
        r"\begin{tikzpicture}[font=\footnotesize]",
        # Subtle quadrant shading
        f"  \\fill[teal!8] ({xr_med:.3f},{yl_med:.3f}) rectangle ({width:.2f},{height:.2f});",
        f"  \\fill[blue!6] ({xr_med:.3f},0) rectangle ({width:.2f},{yl_med:.3f});",
        f"  \\fill[orange!7] (0,{yl_med:.3f}) rectangle ({xr_med:.3f},{height:.2f});",
        f"  \\fill[gray!6] (0,0) rectangle ({xr_med:.3f},{yl_med:.3f});",
        f"  \\draw[gray!50] (0,0) rectangle ({width:.2f},{height:.2f});",
    ]

    for rt in (0.2, 0.4, 0.6, 0.8, 1.0, 1.2):
        xx = sx(rt)
        lines.append(f"  \\draw[gray!20] ({xx:.3f},0) -- ({xx:.3f},{height:.2f});")
        lines.append(f"  \\node[anchor=north, font=\\scriptsize] at ({xx:.3f},-0.08) {{{rt:.1f}}};")

    for lt in (0.0, 0.003, 0.006, 0.009, 0.012, 0.015):
        yy = sy(lt)
        lines.append(f"  \\draw[gray!20] (0,{yy:.3f}) -- ({width:.2f},{yy:.3f});")
        lines.append(f"  \\node[anchor=east, font=\\scriptsize] at (-0.08,{yy:.3f}) {{{lt:.3f}}};")

    # Median lines
    lines.append(f"  \\draw[red!75!black, dashed, thick] ({xr_med:.3f},0) -- ({xr_med:.3f},{height:.2f});")
    lines.append(f"  \\draw[red!75!black, dashed, thick] (0,{yl_med:.3f}) -- ({width:.2f},{yl_med:.3f});")
    lines.append(f"  \\node[red!75!black, anchor=south, font=\\scriptsize] at ({xr_med:.3f},{height+0.06:.3f}) "
                 f"{{$r_{{\\mathrm{{med}}}} = {r_med:.4f}$}};")
    lines.append(f"  \\node[red!75!black, fill=white, fill opacity=0.88, text opacity=1, inner sep=1.5pt, "
                 f"anchor=south east, font=\\scriptsize] at ({width-0.12:.3f},{yl_med+0.08:.3f}) "
                 f"{{$\\ell_{{\\mathrm{{med}}}} = {l_med:.5f}$}};")

    # Plot all 255 segments colored by quadrant
    for s in segs:
        r = float(s["raio"])
        l = float(s["comprimento_3d"])
        in_a = r >= r_med
        in_b = l >= l_med
        if in_a and in_b:
            col = "teal!80!black"
        elif in_a and not in_b:
            col = "blue!70!black"
        elif not in_a and in_b:
            col = "orange!85!black"
        else:
            col = "gray!70!black"
        lines.append(f"  \\fill[{col}] ({sx(r):.3f},{sy(l):.3f}) circle (1.35pt);")

    # Quadrant summary badges
    lines.append(f"  \\node[fill=white, fill opacity=0.88, text opacity=1, draw=teal!70!black, rounded corners=2pt, "
                 f"font=\\scriptsize, anchor=north east] at ({width-0.18:.2f},{height-0.18:.2f}) "
                 f"{{$A \\cap B$: {ie['card_A_inter_B']} seg.}};")
    lines.append(f"  \\node[fill=white, fill opacity=0.88, text opacity=1, draw=blue!70!black, rounded corners=2pt, "
                 f"font=\\scriptsize, anchor=south east] at ({width-0.18:.2f},0.18) "
                 f"{{$A \\setminus B$: {ie['card_A_menos_B']} seg.}};")
    lines.append(f"  \\node[fill=white, fill opacity=0.88, text opacity=1, draw=orange!80!black, rounded corners=2pt, "
                 f"font=\\scriptsize, anchor=north west] at (0.14,{height-0.18:.2f}) "
                 f"{{$B \\setminus A$: {ie['card_B_menos_A']} seg.}};")
    lines.append(f"  \\node[fill=white, fill opacity=0.88, text opacity=1, draw=gray!70!black, rounded corners=2pt, "
                 f"font=\\scriptsize, anchor=south west] at (0.14,0.18) "
                 f"{{$(A \\cup B)^c$: {ie['card_nenhum']} seg.}};")

    lines.append(f"  \\node[font=\\footnotesize, anchor=north] at ({width/2:.2f},-0.55) "
                 r"{Raio do segmento $r_i$ (unidades originais VTK)};")
    lines.append(f"  \\node[font=\\footnotesize, rotate=90, anchor=south] at (-1.05,{height/2:.2f}) "
                 r"{Comprimento 3D $\ell_i$ (unidades originais VTK)};")
    lines.append(r"\end{tikzpicture}")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def generate_matrices_d_figure(res_dir: Path, out_path: Path) -> None:
    def load_matrix(path: Path) -> list[list[int]]:
        with path.open(encoding="utf-8") as f:
            rows = list(csv.reader(f))[1:]
        return [[int(x) for x in r[1:]] for r in rows]

    mat_r = load_matrix(res_dir / "modulo_d_matriz_R_Nterm016.csv")
    mat_rp = load_matrix(res_dir / "modulo_d_matriz_Rplus_Nterm016.csv")
    n = len(mat_r)
    cell = 0.175  # 32 * 0.175 = 5.60 cm per matrix

    lines = [
        "% Figura gerada automaticamente por scripts/gerar_figuras_tikz.py",
        r"\begin{tikzpicture}[font=\scriptsize]",
    ]

    # Left panel: Matrix A of R (32x32)
    size = n * cell
    lines.append(r"  \begin{scope}[xshift=0cm]")
    lines.append(f"    \\fill[gray!6] (0,0) rectangle ({size:.3f},{size:.3f});")
    for k in range(0, n + 1, 4):
        pos = k * cell
        lines.append(f"    \\draw[gray!25, very thin] ({pos:.3f},0) -- ({pos:.3f},{size:.3f});")
        lines.append(f"    \\draw[gray!25, very thin] (0,{size - pos:.3f}) -- ({size:.3f},{size - pos:.3f});")
    for u in range(n):
        for v in range(n):
            if mat_r[u][v] == 1:
                x = v * cell
                y = (n - 1 - u) * cell
                lines.append(f"    \\fill[teal!80!black] ({x:.3f},{y:.3f}) rectangle ({x+cell:.3f},{y+cell:.3f});")
    lines.append(f"    \\draw[gray!70] (0,0) rectangle ({size:.3f},{size:.3f});")
    for t in (0, 8, 16, 24, 31):
        xt = (t + 0.5) * cell
        yt = (n - 1 - t + 0.5) * cell
        lines.append(f"    \\node[anchor=north, font=\\tiny] at ({xt:.3f},-0.05) {{{t}}};")
        lines.append(f"    \\node[anchor=east, font=\\tiny] at (-0.05,{yt:.3f}) {{{t}}};")
    lines.append(f"    \\node[font=\\footnotesize\\bfseries, anchor=south] at ({size/2:.3f},{size+0.12:.3f}) "
                 r"{(a) Matriz $A$ de $\mathcal{R}$ ($|\mathcal{R}| = 31$)};")
    lines.append(f"    \\node[font=\\scriptsize, anchor=north] at ({size/2:.3f},-0.42) {{Coluna $v$ (filho)}};")
    lines.append(f"    \\node[font=\\scriptsize, rotate=90, anchor=south] at (-0.48,{size/2:.3f}) {{Linha $u$ (pai)}};")
    lines.append(r"  \end{scope}")

    # Right panel: Matrix of R^+ (direct in teal, indirect in orange)
    lines.append(r"  \begin{scope}[xshift=7.1cm]")
    lines.append(f"    \\fill[gray!6] (0,0) rectangle ({size:.3f},{size:.3f});")
    for k in range(0, n + 1, 4):
        pos = k * cell
        lines.append(f"    \\draw[gray!25, very thin] ({pos:.3f},0) -- ({pos:.3f},{size:.3f});")
        lines.append(f"    \\draw[gray!25, very thin] (0,{size - pos:.3f}) -- ({size:.3f},{size - pos:.3f});")
    for u in range(n):
        for v in range(n):
            if mat_rp[u][v] == 1:
                x = v * cell
                y = (n - 1 - u) * cell
                col = "teal!80!black" if mat_r[u][v] == 1 else "orange!85!black"
                lines.append(f"    \\fill[{col}] ({x:.3f},{y:.3f}) rectangle ({x+cell:.3f},{y+cell:.3f});")
    lines.append(f"    \\draw[gray!70] (0,0) rectangle ({size:.3f},{size:.3f});")
    for t in (0, 8, 16, 24, 31):
        xt = (t + 0.5) * cell
        yt = (n - 1 - t + 0.5) * cell
        lines.append(f"    \\node[anchor=north, font=\\tiny] at ({xt:.3f},-0.05) {{{t}}};")
        lines.append(f"    \\node[anchor=east, font=\\tiny] at (-0.05,{yt:.3f}) {{{t}}};")
    lines.append(f"    \\node[font=\\footnotesize\\bfseries, anchor=south] at ({size/2:.3f},{size+0.12:.3f}) "
                 r"{(b) Fecho $\mathcal{R}^+$ ($31$ diretos + $148$ indiretos)};")
    lines.append(f"    \\node[font=\\scriptsize, anchor=north] at ({size/2:.3f},-0.42) {{Coluna $v$ (descendente)}};")
    lines.append(f"    \\node[font=\\scriptsize, rotate=90, anchor=south] at (-0.48,{size/2:.3f}) {{Linha $u$ (ancestral)}};")
    lines.append(r"  \end{scope}")

    # Shared compact legend below (well separated)
    lines.append(r"  \fill[teal!80!black] (0.9,-1.05) rectangle (1.25,-0.80);")
    lines.append(r"  \node[anchor=west, font=\scriptsize] at (1.35,-0.92) {Relação direta $(u,v)\in \mathcal{R}$ (31 pares)};")
    lines.append(r"  \fill[orange!85!black] (6.9,-1.05) rectangle (7.25,-0.80);")
    lines.append(r"  \node[anchor=west, font=\scriptsize] at (7.35,-0.92) {Ancestralidade indireta $\mathcal{R}^+\setminus \mathcal{R}$ (148 pares)};")

    lines.append(r"\end{tikzpicture}")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Gera figuras TikZ a partir dos CSVs do CCO Explorer")
    parser.add_argument("--resultados", type=Path, default=Path("resultados"))
    parser.add_argument("--saida", type=Path, default=Path("relatorio/figuras"))
    args = parser.parse_args()

    args.saida.mkdir(parents=True, exist_ok=True)
    generate_recurrence_figure(args.resultados, args.saida / "fig_recorrencia.tex")
    generate_modular_b0_figure(args.resultados, args.saida / "fig_modular_b0.tex")
    generate_inclusion_exclusion_figure(args.resultados, args.saida / "fig_inclusao_exclusao.tex")
    generate_matrices_d_figure(args.resultados, args.saida / "fig_matrizes_d.tex")
    print(f"Figuras TikZ geradas com sucesso em: {args.saida}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

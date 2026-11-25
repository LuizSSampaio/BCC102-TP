#!/usr/bin/env python3
"""
Empacota os entregaveis do Trabalho Pratico CCO Explorer (Secao 9 do enunciado)
em um arquivo .zip reprodutivel, sem duplicar os arquivos .vtk originais.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import zipfile


INCLUDE_FILES = [
    "CMakeLists.txt",
    ".clang-format",
    ".clang-tidy",
    ".envrc",
    ".gitignore",
    "devenv.nix",
    "devenv.yaml",
    "README.md",
]

INCLUDE_DIRS = [
    "include",
    "src",
    "tests",
    "scripts",
    "resultados",
]


def main() -> int:
    parser = argparse.ArgumentParser(description="Cria arquivo ZIP de entrega do CCO Explorer")
    parser.add_argument("--raiz", type=Path, default=Path("."))
    parser.add_argument("--saida", type=Path, default=Path("cco_explorer_entrega.zip"))
    args = parser.parse_args()

    root = args.raiz.resolve()
    out_zip = args.saida.resolve()

    added = 0
    with zipfile.ZipFile(out_zip, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for rel in INCLUDE_FILES:
            p = root / rel
            if p.is_file():
                zf.write(p, arcname=rel)
                added += 1

        for drel in INCLUDE_DIRS:
            dpath = root / drel
            if not dpath.is_dir():
                continue
            for f in sorted(dpath.rglob("*")):
                if f.is_file() and "__pycache__" not in f.parts and f.suffix != ".pyc":
                    zf.write(f, arcname=str(f.relative_to(root)))
                    added += 1

        # Relatorio PDF e fonte LaTeX (sem auxiliares)
        rel_dir = root / "relatorio"
        if rel_dir.is_dir():
            for f in sorted(rel_dir.rglob("*")):
                if f.is_file() and f.suffix in {".pdf", ".tex"}:
                    zf.write(f, arcname=str(f.relative_to(root)))
                    added += 1

    print(f"Pacote de entrega gerado em {out_zip} ({added} arquivos; VTK originais excluidos).")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

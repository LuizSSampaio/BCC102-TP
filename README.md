# CCO Explorer -- BCC102 (Matemática Discreta II | UFOP)

Implementação em **C++20** (com scripts auxiliares em **Python 3** e ambiente reprodutível via **devenv**) do Trabalho Prático **CCO Explorer: Aplicações de Matemática Discreta na Análise de Redes Vasculares** (UFOP -- Período 2026.2).

---

## 1. Ambiente, Linguagem e Versões

O projeto utiliza [`devenv`](https://devenv.sh/) (`devenv.nix` e `devenv.yaml`) para fixar versões equivalentes de compilador, ferramentas de build, formatador, linter, Python e TeX Live:

- **Linguagem principal:** C++20 (`CMake >= 3.20`, compilado com `clang++` / `g++` do `devenv`, sem bibliotecas externas de parsing VTK ou grafos).
- **Scripts auxiliares:** Python 3 (somente biblioteca padrão: verificação cruzada, geração de figuras vetoriais TikZ e empacotamento).
- **Relatório:** LaTeX (`pdflatex` via TeX Live no `devenv`).
- **Trilha eletiva escolhida (Parte B1):** **Trilha 1 -- Algoritmo de Euclides e Teorema Chinês do Resto (`crt`)**.

---

## 2. Estrutura do Projeto

```text
.
├── CMakeLists.txt               # Build CMake (cco_core, cco_explorer, cco_tests, format, tidy)
├── devenv.nix / devenv.yaml     # Ambiente reprodutível (C++, CMake, clang-tools, Python, TeX Live)
├── .clang-format / .clang-tidy  # Padronização de estilo e análise estática
├── include/cco/                 # Cabeçalhos dos módulos A, B, C, D e E
│   ├── model.hpp                # Estruturas Point3D, Segment, VascularModel, ValidationReport
│   ├── vtk_parser.hpp           # Módulo A: leitura ASCII POLYDATA, geometria 3D, validação e graus
│   ├── number_theory.hpp        # Módulo B: B0 (mod 7, PCP) e B1 Trilha 1 (Euclides + TCR)
│   ├── combinatorics.hpp        # Módulo C: C1 (contagem, inclusão-exclusão) e C2 (recorrências)
│   ├── relations.hpp            # Módulo D: R, matriz 32x32, propriedades, p, R+, ordem parcial, ~
│   ├── exporter.hpp             # Exportadores CSV/TXT e Resumo Matemático (Seção 8.2)
│   └── runner.hpp               # Orquestração da linha de comando
├── src/                         # Implementações C++20 e entrypoint CLI (main.cpp)
├── tests/
│   └── test_cco_explorer.cpp    # Suíte automatizada CTest (T1..T7 + testes unitários)
├── scripts/
│   ├── verificar_resultados.py  # Auditoria cruzada independente em Python puro
│   ├── gerar_figuras_tikz.py    # Geração de figuras TikZ a partir dos CSVs gerados pelo C++
│   └── empacotar_entrega.py     # Geração do .zip final (Seção 9, sem duplicar os .vtk)
├── resultados/                  # Saídas geradas (CSV e TXT) para os seis modelos
├── relatorio/
│   ├── relatorio.tex            # Fonte LaTeX do relatório científico-técnico
│   └── relatorio_cco_explorer.pdf
└── dados/                       # Modelos VTK originais (não duplicados no ZIP de entrega)
```

---

## 3. Instalação e Compilação

### Usando `devenv` (recomendado)

```bash
# Entra no ambiente ou executa diretamente via devenv shell
devenv shell -- build
```

### Compilação direta com CMake (fora do `devenv`)

```bash
cmake -B build -S .
cmake --build build
```

---

## 4. Comandos de Execução

O executável `build/cco_explorer` não depende de caminhos absolutos particulares e aceita qualquer diretório de modelos via `--entrada`.

### Modelo individual (Resumo Matemático -- Seção 8.2)

```bash
./build/cco_explorer --entrada dados/ --modelo tree3D_Nterm128.vtk --trilha crt
./build/cco_explorer --entrada dados/ --modelo tree3D_Nterm016.vtk --trilha crt
```

### Análise completa dos seis modelos + exportação CSV/TXT

```bash
./build/cco_explorer --entrada dados/ --todos --trilha crt --saida resultados/
# ou via script do devenv:
devenv shell -- run-all
```

### Ajuda da linha de comando

```bash
./build/cco_explorer --ajuda
```

---

## 5. Testes, Qualidade de Código e Reprodutibilidade

```bash
# 1. Suíte C++ de testes mínimos de aceitação (T1..T7) e unitários via CTest
devenv shell -- run-tests

# 2. Verificação cruzada independente em Python puro (VTK <-> CSVs do C++)
devenv shell -- python3 scripts/verificar_resultados.py --entrada dados --resultados resultados

# 3. Verificação de formatação e análise estática (clang-format e clang-tidy)
devenv shell -- format-check
devenv shell -- tidy

# 4. Regeneração das figuras TikZ e compilação do relatório PDF
devenv shell -- python3 scripts/gerar_figuras_tikz.py --resultados resultados --saida relatorio/figuras
devenv shell -- bash -c 'cd relatorio && pdflatex -interaction=nonstopmode -jobname=relatorio_cco_explorer relatorio.tex && pdflatex -interaction=nonstopmode -jobname=relatorio_cco_explorer relatorio.tex'

# 5. Geração do pacote .zip de entrega (sem duplicar arquivos .vtk)
devenv shell -- python3 scripts/empacotar_entrega.py --saida cco_explorer_entrega.zip
```

---

## 6. Módulos e Trilha Eletiva Implementados

- **Módulo A (Modelagem e validação):** leitura dos 6 arquivos VTK ASCII `POLYDATA`, cálculo de $\ell_i$ e $M_i$ nas unidades originais, validação de cabeçalho/índices/raios/topologia e comparação automática com a Tabela 1.
- **Módulo B0 (Núcleo modular):** classes $C_0, \dots, C_6$ módulo $7$ em `Nterm128` ($|E| = 255$), validação explícita da partição e verificação do Princípio da Casa dos Pombos ($\max |C_j| \ge \lceil 255/7 \rceil = 37$).
- **Módulo B1 -- Trilha 1 (`crt`):** Algoritmo de Euclides e Euclides Estendido para $\gcd(4,5)=1$, resolução pelo Teorema Chinês do Resto ($i \equiv 1 \pmod 4 \land i \equiv 2 \pmod 5 \iff i \equiv 17 \pmod{20}$) e confirmação de coincidência entre seleção direta e equivalente em `Nterm128` (12 segmentos) e `Nterm512` (51 segmentos).
- **Módulo C (Contagem e recorrências):** dupla contagem estrutural ($B = T - 1$, $|E| = 2T - 1$, $|V| = 2T$), $\binom{T}{2}$ e $B \times T$, Inclusão--Exclusão sobre medianas 3D em `Nterm128` ($128 + 128 - 42 = 214$), e recorrências $e_k = 2e_{k-1} + 1 \implies 32\cdot 2^k - 1$, $b_k = 16\cdot 2^k - 1$ e $v_k = 32\cdot 2^k$ validadas nos 6 modelos.
- **Módulo D (Relações):** relação pai--filho $\mathcal{R}$, inversa $\mathcal{R}^{-1}$, matriz densa $32 \times 32$ exclusivamente para `Nterm016`, propriedades com contraexemplos concretos, função pai $p : V \setminus \{0\} \to V$, fecho transitivo iterativo $\mathcal{R}^+$ ($179$ pares), ordem parcial $\preceq$ ($317$ pares incomparáveis) e relação de equivalência modular $\sim$ em `Nterm128`.
- **Módulo E (Integração):** CLI unificada, saídas CSV/TXT em `resultados/`, relatório científico-técnico em PDF (`relatorio/relatorio_cco_explorer.pdf`) e suíte completa de testes T1--T7.

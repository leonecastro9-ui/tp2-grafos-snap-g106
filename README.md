# TP2 — Analisando Grafos do SNAP · Grafo **G106 = `memetracker9`**

**Aluno:** Leone Vicentine de Castro Souza  
**Matrícula:** 224216278  
**Instituição:** Instituto de Computação — Universidade Federal da Bahia (UFBA)  
**Disciplina:** Algoritmos em Grafos / Teoria dos Grafos — Prof. Bruno Pereira dos Santos  
**Repositório:** <https://github.com/leonecastro9-ui/tp2-grafos-snap-g106>

---

## 0. Formato do relatório

O relatório segue o **template oficial do WebMedia** indicado no enunciado
(`webmedia.cls` v1.71, estilo `sigconf` — duas colunas, derivado do
`acmart` da ACM). A classe acompanha o repositório em `relatorio/`.
Para compilar é preciso ter o pacote `acmart` do TeX Live instalado
(`texlive-publishers`), de onde vêm as dependências `libertine`,
`newtxmath` e `ACM-Reference-Format`.

---

## 1. Determinação do grafo

Regra do enunciado:

```
f(M) = ((soma dos dígitos de M) × (últimos 2 dígitos de M + 1)) mod 129
```

Para `M = 224216278`:

| termo | valor |
|---|---|
| soma dos dígitos | 2+2+4+2+1+6+2+7+8 = **34** |
| últimos dois dígitos | **78** |
| cálculo | (34 × (78 + 1)) mod 129 = 34 × 79 = 2686 |
| 2686 mod 129 | 2686 − 20×129 = 2686 − 2580 = **106** |

**f(224216278) = 106 → grafo G106 = `memetracker9`**
(<https://snap.stanford.edu/data/memetracker9.html>)

Verificação reproduzível:

```bash
python3 -c "M='224216278'; print((sum(map(int,M))*(int(M[-2:])+1))%129)"
```

---

## 2. ⚠️ Aviso de recorte

O grafo original tem **96.608.034 documentos** e **418.237.269 hiperlinks**
(≈ 15 GB comprimidos, 9 arquivos mensais). Ele **não cabe em memória** como
grafo de documentos. Foram aplicados dois recortes, ambos previstos no
enunciado ("filtragem por intervalo temporal", "grafo temporal → grafo
estático", "redução do tamanho por questões de memória"):

1. **Recorte temporal** — usou-se apenas o mês de **agosto de 2008**
   (`quotes_2008-08.txt.gz`, 1,05 GB comprimido), o primeiro mês completo
   da coleta. São 7.340.805 documentos e 27.810.398 hiperlinks.
2. **Agregação por site** — cada URL foi reduzida ao seu *hostname*
   normalizado (minúsculas, sem `www.`, sem porta). O grafo deixa de ser
   "documento → documento" e passa a ser **"site → site"**, granularidade
   padrão na literatura de redes de hiperlinks.

Demais tratamentos: dirigido → não dirigido; multigrafo → grafo simples
(peso = nº de hiperlinks); remoção de auto-loops (links internos do próprio
site); extração da maior componente conexa para as métricas de distância.

Tudo está detalhado na Seção "Tratamento dos dados" do relatório.

---

## 3. Estrutura do repositório

```
.
├── README.md
├── relatorio/
│   ├── 224216278.tex        # fonte LaTeX (template oficial do WebMedia)
│   ├── 224216278.pdf        # relatório final (duas colunas, sigconf)
│   ├── webmedia.cls         # classe oficial do WebMedia (v1.71)
│   ├── numeros.tex          # macros com todos os números (geradas)
│   └── tab_*.tex            # tabelas (geradas)
├── scripts/
│   ├── lib_grafo.py             # E/S, CSR, utilidades
│   ├── lib_algoritmos.py        # implementações PRÓPRIAS dos algoritmos
│   ├── lib_dist.py              # diâmetro/raio exatos (BoundingDiameters)
│   ├── 01_construir_grafo.py    # bruto SNAP → arestas dirigidas (site→site)
│   ├── 02_preparar_grafo.py     # não dirigido, simples, componentes
│   ├── 03_analise_estrutural.py # PARTE I
│   ├── 04_algoritmos.py         # PARTE II (+ tempos, IC, escalabilidade)
│   ├── 05_small_world.py        # PARTE III.1 (implementação própria)
│   ├── 06_lei_potencia.py       # PARTE III.2
│   ├── 07_robustez.py           # PARTE III.3
│   ├── 08_modelos_classicos.py  # PARTE III (extra)
│   ├── 09_validacao.py          # validação das implementações vs. igraph
│   ├── 10_gerar_tex.py          # JSON de resultados → macros/tabelas LaTeX
│   ├── 11_validacao_recorte.py  # repete o pipeline em setembro/2008
│   ├── 12_betweenness.py        # Brandes amostrado + ataque por intermediação
│   └── 13_triangulos.py         # contagem exata de triângulos (própria)
├── dados/            # conjunto de dados TRATADO (ver seção 5)
├── resultados/       # JSON/CSV com todos os números do relatório
└── figuras/          # figuras em PDF usadas no relatório
```

---

## 4. Como executar

### Requisitos

```bash
python3 -m pip install numpy scipy matplotlib python-igraph openpyxl
```

Testado com Python 3.13, numpy 2.x, scipy 1.18, igraph 1.0, matplotlib 3.11,
em contêiner Linux x86-64 com 4 vCPU e 15 GB de RAM.

### Passo 1 — baixar o dado bruto (1,05 GB)

```bash
mkdir -p dados_brutos && cd dados_brutos
curl -O https://snap.stanford.edu/data/bigdata/memetracker9/quotes_2008-08.txt.gz
cd ..
```

### Passo 2 — pipeline completo

```bash
python3 scripts/01_construir_grafo.py dados_brutos/quotes_2008-08.txt.gz \
        dados/arestas_dirigidas.csv 2008-08-01 2008-09-01
python3 scripts/02_preparar_grafo.py
python3 scripts/03_analise_estrutural.py     # Parte I
python3 scripts/04_algoritmos.py             # Parte II
python3 scripts/05_small_world.py            # Parte III.1
python3 scripts/06_lei_potencia.py           # Parte III.2
python3 scripts/07_robustez.py               # Parte III.3
python3 scripts/08_modelos_classicos.py      # Parte III (extra)
python3 scripts/09_validacao.py              # validação das implementações
python3 scripts/13_triangulos.py             # triângulos (contagem própria)
python3 scripts/12_betweenness.py            # intermediação (~20 min)
python3 scripts/11_validacao_recorte.py      # validação do recorte (~15 min)
python3 scripts/10_gerar_tex.py              # gera numeros.tex e tab_*.tex
```

O passo `11_validacao_recorte.py` precisa também do arquivo de setembro:

```bash
cd dados_brutos
curl -O https://snap.stanford.edu/data/bigdata/memetracker9/quotes_2008-09.txt.gz
cd ..
```

### Passo 3 — compilar o relatório

```bash
cd relatorio && pdflatex 224216278.tex && pdflatex 224216278.tex
```

---

## 5. Conjunto de dados tratado (`dados/`)

| arquivo | conteúdo |
|---|---|
| `vertices.csv` | `id,hostname,componente` — 673.017 sites |
| `arestas.csv` | `u,v,peso` — 1.614.246 arestas não dirigidas simples |
| `arestas_cg.csv` | idem, restrito à componente gigante (1.590.365 arestas) |
| `arestas_dirigidas.csv` | `origem,destino,peso` — 1.621.841 arestas dirigidas (hostnames) |
| `arestas_dirigidas_id.csv` | idem, com ids inteiros (usado no algoritmo de Tarjan) |
| `componentes.csv` | distribuição de tamanhos das componentes |

`dados_brutos/` (1,05 GB) **não** é versionado; o Passo 1 o reconstrói.
Também não são versionados `arestas_dirigidas.csv` (71 MB) e
`arestas_dirigidas_id.csv` (27 MB), que são intermediários regenerados
pelos Passos 1 e 2 — o conjunto tratado final (`vertices.csv`,
`arestas.csv`, `arestas_cg.csv`, `componentes.csv`) está versionado.

---

## 6. Implementações próprias

Nenhuma biblioteca de grafos é usada nos algoritmos da Parte II nem no
cálculo de *small-world* da Parte III — apenas Python puro sobre uma
representação CSR (`lib_algoritmos.py`, `lib_dist.py`):

BFS · DFS · verificação de eulerianidade · Dijkstra · Bellman-Ford ·
Floyd-Warshall (puro e vetorizado) · Tarjan (componentes fortemente conexas,
vértices de articulação e pontes) · Prim · Kruskal ·
contagem exata de triângulos e coeficiente de clusterização ·
comprimento médio dos caminhos · gerador Erdős–Rényi ·
modelo de configuração por trocas duplas de arestas ·
BoundingDiameters (diâmetro e raio exatos) ·
Brandes amostrado (centralidade de intermediação).

O igraph é usado **apenas** como (a) referência de validação cruzada
(`09_validacao.py`) e (b) motor de BFS em C nas rotinas de amostragem
massiva da Parte I e da análise de robustez, onde se executam dezenas de
milhares de buscas. Isso está declarado no relatório.

`python3 scripts/09_validacao.py` confere todas as implementações próprias
contra o igraph em grafos aleatórios pequenos.

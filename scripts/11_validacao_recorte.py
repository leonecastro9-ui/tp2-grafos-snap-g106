#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VALIDACAO DO RECORTE TEMPORAL.

O relatorio usa apenas agosto de 2008.  A objecao obvia e: as conclusoes
seriam as mesmas em outro mes, ou sao um artefato da janela escolhida?

Este script repete o pipeline inteiro sobre SETEMBRO de 2008
(quotes_2008-09.txt.gz) e compara, lado a lado, todas as grandezas das
quais o relatorio tira conclusoes:

  * ordem, tamanho, grau medio/maximo, densidade;
  * numero de componentes e fracao na componente gigante;
  * clusterizacao media e transitividade global;
  * diametro e raio EXATOS (BoundingDiameters);
  * comprimento medio dos caminhos (amostrado, com IC 95%);
  * expoente gamma da cauda (MLE de Clauset com k_min por KS);
  * robustez: S apos remocao aleatoria de 5% (T repeticoes) e apos
    ataque aos 5% de maior grau.

Saida: resultados/11_validacao_recorte.json
"""
import os, sys, json, math, random, csv, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import igraph as ig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lib_grafo import BASE, DADOS, RES, FIG
from lib_dist import bounding_diameters

random.seed(31416); np.random.seed(31416)
Z = 1.959963985
AMOSTRA_L = 500
T_ROB = 10

BRUTO_SET = os.path.join(BASE, "dados_brutos", "quotes_2008-09.txt.gz")
ARESTAS_SET = os.path.join(DADOS, "arestas_dirigidas_set.csv")


def gamma_mle(x, kmin):
    xt = x[x >= kmin]
    if len(xt) < 10:
        return None, 0
    return 1.0 + len(xt) / np.log(xt / (kmin - 0.5)).sum(), len(xt)


def ks_dist(x, kmin, g):
    xt = np.sort(x[x >= kmin]); nt = len(xt)
    emp = np.arange(nt) / nt
    teo = 1.0 - (xt / (kmin - 0.5)) ** (1.0 - g)
    return float(np.max(np.abs(emp - teo)))


def ajusta_cauda(graus):
    x = graus[graus > 0].astype(float)
    kmins = np.unique(np.concatenate([
        np.arange(1, 51),
        np.unique(np.round(np.logspace(np.log10(51), np.log10(500), 40)))]))
    melhor = None
    for kmin in kmins:
        g, nt = gamma_mle(x, kmin)
        if g is None:
            continue
        D = ks_dist(x, kmin, g)
        if melhor is None or D < melhor[2]:
            melhor = (int(kmin), float(g), D, int(nt))
    kmin, g, D, nt = melhor
    return {"k_min": kmin, "gamma": g, "ks_D": D, "n_cauda": nt,
            "erro_padrao_gamma": (g - 1) / math.sqrt(nt)}


def caminho_medio(gc, k=AMOSTRA_L):
    fontes = random.sample(range(gc.vcount()), min(k, gc.vcount()))
    medias = []
    for s in fontes:
        d = np.array(gc.distances(source=[s])[0], dtype=float)
        d = d[d > 0]
        medias.append(d.mean())
    a = np.array(medias)
    e = Z * a.std(ddof=1) / math.sqrt(len(a))
    return {"L": float(a.mean()), "ic95": [float(a.mean() - e), float(a.mean() + e)],
            "n_fontes": len(a)}


def metricas(nome, arestas_csv):
    """Carrega um CSV de arestas dirigidas (hostnames) e mede tudo."""
    und = {}
    with open(arestas_csv, encoding="utf-8") as f:
        r = csv.reader(f); next(r)
        for a, b, p in r:
            k = (a, b) if a < b else (b, a)
            und[k] = und.get(k, 0) + int(p)
    nomes = sorted({x for k in und for x in k})
    idx = {s: i for i, s in enumerate(nomes)}
    n = len(nomes)
    arestas = [(idx[a], idx[b]) for a, b in und]
    G = ig.Graph(n=n, edges=arestas, directed=False)
    m = G.ecount()
    graus = np.array(G.degree(), dtype=np.int64)
    comp = G.connected_components()
    tam = np.array(comp.sizes())
    gc = comp.giant()

    print(f"   [{nome}] n={n:,} m={m:,} CG={gc.vcount():,}", flush=True)
    bd = bounding_diameters(gc, verbose=False)
    L = caminho_medio(gc)
    cauda = ajusta_cauda(graus)

    # robustez resumida
    ncg = gc.vcount(); r = int(0.05 * ncg)
    s_rand = []
    for _ in range(T_ROB):
        H = gc.copy(); H.delete_vertices(random.sample(range(ncg), r))
        s_rand.append(max(H.connected_components().sizes()) / ncg)
    s_rand = np.array(s_rand)
    gcg = np.array(gc.degree())
    H = gc.copy(); H.delete_vertices(np.argsort(-gcg)[:r].tolist())
    s_cent = max(H.connected_components().sizes()) / ncg

    return {
        "mes": nome,
        "n_vertices": n, "n_arestas": m,
        "grau_medio": float(graus.mean()), "grau_max": int(graus.max()),
        "grau_mediano": float(np.median(graus)),
        "densidade": 2.0 * m / (n * (n - 1)),
        "n_componentes": int(len(tam)),
        "fracao_componente_gigante": float(tam.max() / n),
        "clusterizacao_media": float(np.mean(
            G.transitivity_local_undirected(mode="zero"))),
        "transitividade_global": float(G.transitivity_undirected(mode="zero")),
        "diametro": bd["diametro"], "raio": bd["raio"],
        "diametro_exato": bd["exato"],
        "caminho_medio": L,
        "cauda": cauda,
        "robustez": {
            "r": r, "T": T_ROB,
            "S_aleatorio_media": float(s_rand.mean()),
            "S_aleatorio_desvio": float(s_rand.std(ddof=1)),
            "S_ataque_grau": float(s_cent),
            "razao": float(s_cent / s_rand.mean()),
        },
    }


def main():
    if not os.path.exists(ARESTAS_SET):
        print(">> construindo o grafo de setembro/2008 a partir do bruto", flush=True)
        subprocess.run([sys.executable, "-I",
                        os.path.join(BASE, "scripts", "01_construir_grafo.py"),
                        BRUTO_SET, ARESTAS_SET, "2008-09-01", "2008-10-01"],
                       check=True)

    print(">> medindo agosto/2008 (referencia do relatorio)", flush=True)
    ago = metricas("agosto/2008", os.path.join(DADOS, "arestas_dirigidas.csv"))
    print(">> medindo setembro/2008 (validacao)", flush=True)
    set_ = metricas("setembro/2008", ARESTAS_SET)

    def var(a, b):
        return 100.0 * (b - a) / a if a else float("nan")

    comp = {}
    for ch in ["n_vertices", "n_arestas", "grau_medio", "grau_max", "densidade",
               "n_componentes", "fracao_componente_gigante",
               "clusterizacao_media", "transitividade_global", "diametro", "raio"]:
        comp[ch] = {"agosto": ago[ch], "setembro": set_[ch],
                    "variacao_pct": var(ago[ch], set_[ch])}
    comp["caminho_medio"] = {"agosto": ago["caminho_medio"]["L"],
                             "setembro": set_["caminho_medio"]["L"],
                             "variacao_pct": var(ago["caminho_medio"]["L"],
                                                 set_["caminho_medio"]["L"])}
    comp["gamma"] = {"agosto": ago["cauda"]["gamma"],
                     "setembro": set_["cauda"]["gamma"],
                     "variacao_pct": var(ago["cauda"]["gamma"], set_["cauda"]["gamma"])}
    comp["k_min"] = {"agosto": ago["cauda"]["k_min"], "setembro": set_["cauda"]["k_min"]}
    comp["S_aleatorio"] = {"agosto": ago["robustez"]["S_aleatorio_media"],
                           "setembro": set_["robustez"]["S_aleatorio_media"],
                           "variacao_pct": var(ago["robustez"]["S_aleatorio_media"],
                                               set_["robustez"]["S_aleatorio_media"])}
    comp["S_ataque_grau"] = {"agosto": ago["robustez"]["S_ataque_grau"],
                             "setembro": set_["robustez"]["S_ataque_grau"],
                             "variacao_pct": var(ago["robustez"]["S_ataque_grau"],
                                                 set_["robustez"]["S_ataque_grau"])}

    saida = {"agosto": ago, "setembro": set_, "comparacao": comp}
    with open(os.path.join(RES, "11_validacao_recorte.json"), "w", encoding="utf-8") as f:
        json.dump(saida, f, indent=2, ensure_ascii=False)

    # figura: variacao percentual de cada grandeza
    rot = ["$n$", "$m$", r"$\langle k\rangle$", "densidade", "comp. gigante",
           "$C$", "$L$", "diâmetro", r"$\gamma$", "$S_{rand}$", "$S_{cent}$"]
    chs = ["n_vertices", "n_arestas", "grau_medio", "densidade",
           "fracao_componente_gigante", "clusterizacao_media", "caminho_medio",
           "diametro", "gamma", "S_aleatorio", "S_ataque_grau"]
    vals = [comp[c]["variacao_pct"] for c in chs]
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    LIM = 60.0                       # corte do eixo; barras maiores sao anotadas
    plot = [min(max(v, -LIM), LIM) for v in vals]
    cores = ["#1f4e79" if abs(v) < 15 else "#c00000" for v in vals]
    ax.bar(range(len(vals)), plot, color=cores)
    ax.axhline(0, color="k", lw=.8)
    ax.axhspan(-15, 15, color="#1f4e79", alpha=.08)
    for i, v in enumerate(vals):
        if abs(v) > LIM:             # barra cortada: anota o valor real
            ax.annotate(f"{v:+.0f}%", xy=(i, LIM), xytext=(i, LIM * 0.45),
                        ha="center", fontsize=7.5, color="white",
                        fontweight="bold", rotation=90)
    ax.set_ylim(-LIM * 1.05, LIM * 1.05)
    ax.set_xticks(range(len(vals)))
    ax.set_xticklabels(rot, fontsize=8, rotation=30)
    ax.set_ylabel("variacao de setembro\nem relacao a agosto (%)", fontsize=8)
    ax.set_title("Estabilidade das medidas entre dois meses "
                 "(eixo cortado em $\\pm$60%)", fontsize=10)
    ax.grid(alpha=.3, axis="y")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_recorte.pdf"))
    plt.close(fig)

    print(json.dumps(comp, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

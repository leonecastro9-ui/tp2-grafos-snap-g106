#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PARTE III (EXTRA) - O grafo se aproxima de algum modelo classico?

Compara G106 com tres modelos gerados com o MESMO n e (aproximadamente)
o MESMO m:
  * Erdos-Renyi  G(n, m)
  * Barabasi-Albert (crescimento + ligacao preferencial)
  * Watts-Strogatz (anel regular com religacao p)
Metricas comparadas: grau medio, grau maximo, coeficiente de clusterizacao
medio, comprimento medio dos caminhos (amostrado), assortatividade de grau
e expoente da cauda da distribuicao de graus.
"""
import os, sys, json, math, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import igraph as ig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lib_grafo import ler_arestas, compactar, DADOS, RES, FIG

random.seed(99); np.random.seed(99)
AMOSTRA = 300
Z = 1.959963985


def resumo(g, nome):
    comp = g.connected_components()
    gc = comp.giant()
    graus = np.array(g.degree(), dtype=float)
    k = min(AMOSTRA, gc.vcount())
    medias = []
    for s in random.sample(range(gc.vcount()), k):
        d = np.array(gc.distances(source=[s])[0], dtype=float)
        d = d[d > 0]
        medias.append(d.mean())
    a = np.array(medias)
    x = graus[graus > 0]
    kmin = max(2.0, np.percentile(x, 95))
    xt = x[x >= kmin]
    gam = 1 + len(xt) / np.log(xt / (kmin - 0.5)).sum() if len(xt) > 20 else float("nan")
    return {
        "modelo": nome, "n": g.vcount(), "m": g.ecount(),
        "grau_medio": float(graus.mean()), "grau_max": int(graus.max()),
        "clusterizacao_media": float(np.mean(g.transitivity_local_undirected(mode="zero"))),
        "transitividade_global": float(g.transitivity_undirected(mode="zero")),
        "caminho_medio": float(a.mean()),
        "caminho_medio_ic95": [float(a.mean() - Z * a.std(ddof=1) / math.sqrt(len(a))),
                               float(a.mean() + Z * a.std(ddof=1) / math.sqrt(len(a)))],
        "fracao_componente_gigante": gc.vcount() / g.vcount(),
        "assortatividade_grau": float(g.assortativity_degree()),
        "gamma_cauda_5pct": float(gam),
    }


def main():
    U, V, P = ler_arestas(os.path.join(DADOS, "arestas.csv"))
    u, v, rot = compactar(U, V)
    n = len(rot)
    G = ig.Graph(n=n, edges=list(zip(u.tolist(), v.tolist())), directed=False)
    m = G.ecount()
    print(f"n={n:,} m={m:,}", flush=True)

    R = {"real": resumo(G, "G106 (memetracker9, ago/2008)")}
    print(">> Erdos-Renyi", flush=True)
    R["er"] = resumo(ig.Graph.Erdos_Renyi(n=n, m=m), "Erdos-Renyi G(n,m)")
    print(">> Barabasi-Albert", flush=True)
    mba = max(1, round(m / n))
    R["ba"] = resumo(ig.Graph.Barabasi(n=n, m=mba), f"Barabasi-Albert (m={mba})")
    print(">> Watts-Strogatz", flush=True)
    nei = max(1, round(m / n))
    R["ws"] = resumo(ig.Graph.Watts_Strogatz(1, n, nei, 0.05), f"Watts-Strogatz (k={2*nei}, p=0.05)")

    with open(os.path.join(RES, "08_modelos.json"), "w", encoding="utf-8") as f:
        json.dump(R, f, indent=2, ensure_ascii=False)

    fig, ax = plt.subplots(figsize=(5.6, 4))
    cores = {"real": "#1f4e79", "er": "#c00000", "ba": "#375623", "ws": "#bf8f00"}
    for ch, lbl in [("real", "G106 (real)"), ("er", "Erdos-Renyi"),
                    ("ba", "Barabasi-Albert"), ("ws", "Watts-Strogatz")]:
        pass
    # distribuicoes de grau comparadas
    gs = {"G106 (real)": np.array(G.degree()),
          "Erdos-Renyi": np.array(ig.Graph.Erdos_Renyi(n=n, m=m).degree()),
          "Barabasi-Albert": np.array(ig.Graph.Barabasi(n=n, m=mba).degree()),
          "Watts-Strogatz": np.array(ig.Graph.Watts_Strogatz(1, n, nei, 0.05).degree())}
    for (lbl, arr), c in zip(gs.items(), ["#1f4e79", "#c00000", "#375623", "#bf8f00"]):
        a = arr[arr > 0]
        xs = np.sort(a); cc = 1.0 - np.arange(len(xs)) / len(xs)
        ax.loglog(xs, cc, "-", lw=1.6, color=c, label=lbl)
    ax.set_xlabel("grau $k$"); ax.set_ylabel("$P(K \\geq k)$")
    ax.set_title("CCDF dos graus: grafo real vs. modelos classicos")
    ax.legend(fontsize=8); ax.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_modelos.pdf")); plt.close(fig)
    print(json.dumps(R, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

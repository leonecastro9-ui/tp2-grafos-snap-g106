#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PARTE I - Analise estrutural obrigatoria.

Calcula: numero de vertices e arestas, graus (min/max/medio), distribuicao
de graus, densidade, componentes conexas e seus tamanhos, diametro, raio,
comprimento medio dos caminhos, coeficiente de clusterizacao medio e numero
de triangulos.  Gera tambem a visualizacao reduzida do grafo.
"""
import os, sys, json, time, math, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import igraph as ig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lib_grafo import ler_arestas, compactar, DADOS, RES, FIG, nomes_vertices
from lib_dist import bounding_diameters

random.seed(42); np.random.seed(42)
AMOSTRA_CAMINHOS = 1000          # fontes amostradas para o caminho medio
os.makedirs(FIG, exist_ok=True); os.makedirs(RES, exist_ok=True)


def monta(arquivo):
    U, V, P = ler_arestas(arquivo)
    u, v, lab = compactar(U, V)
    g = ig.Graph(n=len(lab), edges=list(zip(u.tolist(), v.tolist())), directed=False)
    g.es["peso"] = P.tolist()
    return g, lab


def main():
    R = {}
    print(">> carregando grafo completo")
    G, labG = monta(os.path.join(DADOS, "arestas.csv"))
    n, m = G.vcount(), G.ecount()
    graus = np.array(G.degree(), dtype=np.int64)

    R["n_vertices"] = n
    R["n_arestas"] = m
    R["grau_min"] = int(graus.min())
    R["grau_max"] = int(graus.max())
    R["grau_medio"] = float(graus.mean())
    R["grau_mediano"] = float(np.median(graus))
    R["grau_desvio"] = float(graus.std(ddof=1))
    R["densidade"] = 2.0 * m / (n * (n - 1))

    comp = G.connected_components()
    tam = np.array(comp.sizes())
    R["n_componentes"] = len(tam)
    R["maiores_componentes"] = sorted(tam.tolist(), reverse=True)[:10]
    R["componentes_tamanho_1"] = int((tam == 1).sum())
    R["componentes_tamanho_2"] = int((tam == 2).sum())
    R["componentes_menores_que_10"] = int((tam < 10).sum())

    # --- distribuicao de graus ---
    vals, cnt = np.unique(graus, return_counts=True)
    with open(os.path.join(RES, "distribuicao_graus.csv"), "w", encoding="utf-8") as f:
        f.write("grau,frequencia,probabilidade\n")
        for k, c in zip(vals, cnt):
            f.write(f"{k},{c},{c/n:.10g}\n")

    # --- hubs ---
    nomes = nomes_vertices()
    ordem = np.argsort(-graus)[:25]
    R["hubs"] = [{"site": nomes[int(labG[i])], "grau": int(graus[i])} for i in ordem]

    # --- clusterizacao e triangulos ---
    print(">> clusterizacao")
    t0 = time.time()
    R["transitividade_global"] = float(G.transitivity_undirected(mode="zero"))
    loc = np.array(G.transitivity_local_undirected(mode="zero"))
    R["clusterizacao_media"] = float(loc.mean())
    triplas_v = graus * (graus - 1) / 2.0
    R["n_triplas_conectadas"] = float(triplas_v.sum())
    R["n_triangulos"] = int(round(float((loc * triplas_v).sum()) / 3.0))
    R["tempo_clusterizacao_s"] = time.time() - t0

    # --- componente gigante ---
    print(">> componente gigante")
    CG, labCG = monta(os.path.join(DADOS, "arestas_cg.csv"))
    ncg, mcg = CG.vcount(), CG.ecount()
    gcg = np.array(CG.degree(), dtype=np.int64)
    R["cg"] = {
        "n_vertices": ncg, "n_arestas": mcg,
        "fracao_vertices": ncg / n, "fracao_arestas": mcg / m,
        "grau_min": int(gcg.min()), "grau_max": int(gcg.max()),
        "grau_medio": float(gcg.mean()),
        "densidade": 2.0 * mcg / (ncg * (ncg - 1)),
        "clusterizacao_media": float(np.array(
            CG.transitivity_local_undirected(mode="zero")).mean()),
        "transitividade_global": float(CG.transitivity_undirected(mode="zero")),
    }

    # --- diametro e raio exatos ---
    print(">> diametro e raio (BoundingDiameters)")
    t0 = time.time()
    bd = bounding_diameters(CG)
    bd["tempo_s"] = time.time() - t0
    R["cg"]["distancias"] = bd

    # --- comprimento medio dos caminhos (estimativa por amostragem) ---
    print(">> comprimento medio dos caminhos (amostragem)")
    t0 = time.time()
    fontes = random.sample(range(ncg), AMOSTRA_CAMINHOS)
    medias, maxs = [], []
    hist = np.zeros(64, dtype=np.int64)
    for i, s in enumerate(fontes):
        d = np.array(CG.distances(source=[s])[0], dtype=np.int64)
        d = d[d > 0]
        medias.append(d.mean()); maxs.append(int(d.max()))
        hist[:d.max() + 1] += np.bincount(d, minlength=d.max() + 1)
        if (i + 1) % 200 == 0:
            print(f"    {i+1}/{AMOSTRA_CAMINHOS}", flush=True)
    medias = np.array(medias)
    mu, sd = medias.mean(), medias.std(ddof=1)
    ic = 1.959963985 * sd / math.sqrt(len(medias))
    R["cg"]["caminho_medio"] = {
        "estimativa": float(mu), "desvio_padrao": float(sd),
        "ic95": [float(mu - ic), float(mu + ic)],
        "n_fontes": AMOSTRA_CAMINHOS,
        "excentricidade_max_amostrada": int(max(maxs)),
        "tempo_s": time.time() - t0,
    }
    hist = hist[:int(max(maxs)) + 1]
    with open(os.path.join(RES, "distribuicao_distancias.csv"), "w", encoding="utf-8") as f:
        f.write("distancia,pares_amostrados,fracao\n")
        tot = hist.sum()
        for d_, c in enumerate(hist):
            if d_ > 0:
                f.write(f"{d_},{c},{c/tot:.10g}\n")

    with open(os.path.join(RES, "03_estrutura.json"), "w", encoding="utf-8") as f:
        json.dump(R, f, indent=2, ensure_ascii=False)

    # ================= FIGURAS =================
    print(">> figuras")
    # 1. distribuicao de graus (log-log) + CCDF
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].loglog(vals[vals > 0], cnt[vals > 0] / n, "o", ms=3, alpha=.6, color="#1f4e79")
    ax[0].set_xlabel("grau $k$"); ax[0].set_ylabel("$P(k)$")
    ax[0].set_title("Distribuicao de graus (log-log)")
    ax[0].grid(alpha=.3, which="both")
    ks = np.sort(graus[graus > 0])
    ccdf = 1.0 - np.arange(len(ks)) / len(ks)
    ax[1].loglog(ks, ccdf, "-", color="#c00000")
    ax[1].set_xlabel("grau $k$"); ax[1].set_ylabel("$P(K \\geq k)$")
    ax[1].set_title("CCDF dos graus")
    ax[1].grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_graus.pdf")); plt.close(fig)

    # 2. distribuicao de distancias
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    dd = np.arange(1, len(hist))
    ax.bar(dd, hist[1:] / hist.sum(), color="#1f4e79")
    ax.set_xlabel("distancia $d$"); ax.set_ylabel("fracao dos pares")
    ax.set_title("Distribuicao de distancias (componente gigante)")
    ax.grid(alpha=.3, axis="y")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_distancias.pdf")); plt.close(fig)

    # 3. tamanho das componentes
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    tvals, tcnt = np.unique(tam, return_counts=True)
    ax.loglog(tvals, tcnt, "s", ms=4, color="#375623")
    ax.set_xlabel("tamanho da componente"); ax.set_ylabel("quantidade")
    ax.set_title("Tamanho das componentes conexas")
    ax.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_componentes.pdf")); plt.close(fig)

    # 4. visualizacao reduzida: subgrafo dos 400 vertices de maior grau
    top = np.argsort(-gcg)[:400]
    sub = CG.subgraph(top.tolist())
    sub.delete_vertices([v.index for v in sub.vs if v.degree() == 0])
    lay = sub.layout_fruchterman_reingold(niter=500)
    gs = np.array(sub.degree(), dtype=float)
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    coords = np.array(lay.coords)
    for e in sub.es:
        a, b = e.tuple
        ax.plot(coords[[a, b], 0], coords[[a, b], 1], "-", lw=.25, color="#999999", alpha=.45, zorder=1)
    ax.scatter(coords[:, 0], coords[:, 1], s=8 + 60 * gs / gs.max(),
               c=gs, cmap="viridis", zorder=2, edgecolors="none")
    ax.set_axis_off()
    ax.set_title(f"Subgrafo dos {sub.vcount()} sites de maior grau\n"
                 f"({sub.ecount()} arestas; cor/tamanho = grau no subgrafo)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_visualizacao.pdf"), dpi=300); plt.close(fig)

    print(json.dumps({k: v for k, v in R.items() if k != "hubs"}, indent=2, ensure_ascii=False)[:4000])


if __name__ == "__main__":
    main()

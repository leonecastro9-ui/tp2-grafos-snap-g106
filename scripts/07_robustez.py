#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PARTE III.3 - Robustez (Apendice C).

Politica adotada: todas as analises sao feitas sobre a COMPONENTE GIGANTE
do grafo tratado.  Define-se n = |V| e r = floor(0.05 n).

Cenarios:
  (0) grafo intacto (referencia);
  (1) remocao ALEATORIA de r vertices, com T repeticoes;
  (2) ATAQUE DIRECIONADO: remocao dos r vertices de maior grau
      (ranking fixo, conforme o apendice);
  (3) ATAQUE DIRECIONADO alternativo: remocao dos r vertices de maior
      PageRank.

Metricas apos a remocao:
  A - tamanho da maior componente (absoluto e relativo a n-r)
  B - numero de componentes
  C - distancias: caminho medio estimado por amostragem + diametro exato
      (BoundingDiameters) da maior componente
  D - fracao de vertices isolados

Observacao sobre betweenness: a centralidade de intermediacao exata custa
O(n*m) ~ 1e12 operacoes neste grafo, inviavel; adotou-se grau (recomendado
no enunciado) e PageRank como segunda centralidade.
"""
import os, sys, json, math, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import igraph as ig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lib_grafo import ler_arestas, compactar, DADOS, RES, FIG, nomes_vertices
from lib_dist import bounding_diameters

random.seed(2024); np.random.seed(2024)
T = 30                 # repeticoes da remocao aleatoria
AMOSTRA_L = 150        # fontes de BFS para estimar o caminho medio
Z = 1.959963985


def metricas(g, n_orig, calcular_distancias=True):
    comp = g.connected_components()
    tam = np.array(comp.sizes())
    maior = int(tam.max())
    graus = np.array(g.degree())
    res = {
        "A_maior_componente": maior,
        "A_fracao": maior / g.vcount(),
        "A_fracao_do_original": maior / n_orig,
        "B_n_componentes": int(len(tam)),
        "D_isolados": int((graus == 0).sum()),
        "D_fracao_isolados": float((graus == 0).mean()),
    }
    if calcular_distancias:
        gc = comp.giant()
        k = min(AMOSTRA_L, gc.vcount())
        fontes = random.sample(range(gc.vcount()), k)
        medias = []
        for s in fontes:
            d = np.array(gc.distances(source=[s])[0], dtype=float)
            d = d[d > 0]
            medias.append(d.mean())
        a = np.array(medias)
        res["C_caminho_medio"] = float(a.mean())
        res["C_caminho_medio_ic95"] = [float(a.mean() - Z * a.std(ddof=1) / math.sqrt(len(a))),
                                       float(a.mean() + Z * a.std(ddof=1) / math.sqrt(len(a)))]
        bd = bounding_diameters(gc, verbose=False)
        res["C_diametro"] = bd["diametro"]
        res["C_diametro_exato"] = bd["exato"]
    return res


def main():
    U, V, P = ler_arestas(os.path.join(DADOS, "arestas_cg.csv"))
    u, v, rot = compactar(U, V)
    n = len(rot)
    G = ig.Graph(n=n, edges=list(zip(u.tolist(), v.tolist())), directed=False)
    r = int(math.floor(0.05 * n))
    print(f"n={n:,}  m={G.ecount():,}  r={r:,}", flush=True)
    R = {"n": n, "m": G.ecount(), "r": r, "T": T, "politica": "componente gigante"}

    print(">> grafo intacto", flush=True)
    R["intacto"] = metricas(G, n)
    print(json.dumps(R["intacto"], indent=2), flush=True)

    # ---------------- remocao aleatoria ----------------
    print(">> remocao aleatoria (T repeticoes)", flush=True)
    acum = []
    t0 = time.time()
    for t in range(T):
        alvo = random.sample(range(n), r)
        H = G.copy(); H.delete_vertices(alvo)
        acum.append(metricas(H, n))
        print(f"   rep {t+1}/{T}  S={acum[-1]['A_fracao_do_original']:.4f} "
              f"comp={acum[-1]['B_n_componentes']}  ({time.time()-t0:.0f}s)", flush=True)
    R["aleatorio"] = {"repeticoes": acum}
    for ch in ["A_maior_componente", "A_fracao_do_original", "B_n_componentes",
               "C_caminho_medio", "C_diametro", "D_fracao_isolados"]:
        a = np.array([x[ch] for x in acum], dtype=float)
        e = Z * a.std(ddof=1) / math.sqrt(len(a))
        R["aleatorio"][ch] = {"media": float(a.mean()), "desvio": float(a.std(ddof=1)),
                              "ic95": [float(a.mean() - e), float(a.mean() + e)],
                              "min": float(a.min()), "max": float(a.max())}

    # ---------------- ataque por grau ----------------
    print(">> ataque direcionado por GRAU", flush=True)
    graus = np.array(G.degree())
    alvo_g = np.argsort(-graus)[:r].tolist()
    nomes = nomes_vertices()
    R["ataque_grau_top10"] = [{"site": nomes[int(rot[i])], "grau": int(graus[i])}
                              for i in alvo_g[:10]]
    R["ataque_grau_arestas_removidas"] = int(graus[alvo_g].sum())
    H = G.copy(); H.delete_vertices(alvo_g)
    R["ataque_grau"] = metricas(H, n)
    print(json.dumps(R["ataque_grau"], indent=2), flush=True)

    # ---------------- ataque por PageRank ----------------
    print(">> ataque direcionado por PAGERANK", flush=True)
    pr = np.array(G.pagerank(damping=0.85))
    alvo_p = np.argsort(-pr)[:r].tolist()
    R["ataque_pagerank_sobreposicao_com_grau"] = len(set(alvo_p) & set(alvo_g)) / r
    H = G.copy(); H.delete_vertices(alvo_p)
    R["ataque_pagerank"] = metricas(H, n)
    print(json.dumps(R["ataque_pagerank"], indent=2), flush=True)

    # ---------------- curva de degradacao ----------------
    print(">> curva de degradacao (0% a 20%)", flush=True)
    fracoes = [0.0, 0.01, 0.02, 0.05, 0.10, 0.15, 0.20]
    curva = {"fracoes": fracoes, "aleatorio": [], "grau": []}
    ordem_grau = np.argsort(-graus)
    for fr in fracoes:
        k = int(fr * n)
        if k == 0:
            s_al = s_gr = R["intacto"]["A_fracao_do_original"]
        else:
            H = G.copy(); H.delete_vertices(random.sample(range(n), k))
            s_al = max(H.connected_components().sizes()) / n
            H = G.copy(); H.delete_vertices(ordem_grau[:k].tolist())
            s_gr = max(H.connected_components().sizes()) / n
        curva["aleatorio"].append(s_al); curva["grau"].append(s_gr)
        print(f"   f={fr:.2f}  S_rand={s_al:.4f}  S_cent={s_gr:.4f}", flush=True)
    R["curva_degradacao"] = curva

    with open(os.path.join(RES, "07_robustez.json"), "w", encoding="utf-8") as f:
        json.dump(R, f, indent=2, ensure_ascii=False)

    # ---------------- figuras ----------------
    fig, axs = plt.subplots(2, 2, figsize=(7.6, 6.0))
    axs = axs.ravel()
    chaves = [("A_fracao_do_original", "A: maior componente / $n$", False),
              ("B_n_componentes", "B: n$^{o}$ de componentes", True),
              ("C_caminho_medio", "C: caminho medio", False),
              ("D_fracao_isolados", "D: fracao de isolados", False)]
    for ax, (ch, titulo, logy) in zip(axs, chaves):
        dados = np.array([x[ch] for x in acum], dtype=float)
        ax.boxplot([dados], widths=.45, patch_artist=True,
                   boxprops=dict(facecolor="#cfe2f3"),
                   medianprops=dict(color="#1f4e79"))
        mu, sd = dados.mean(), dados.std(ddof=1)
        ax.errorbar([1], [mu], yerr=[sd], fmt="o", color="#1f4e79", capsize=4,
                    label="aleatorio: media $\\pm$ dp")
        intac, atac = R["intacto"][ch], R["ataque_grau"][ch]
        lo, hi = min(dados.min(), intac), max(dados.max(), intac)
        span = max(hi - lo, 1e-9); lo -= 0.6 * span; hi += 0.6 * span
        ax.axhline(intac, color="#777777", ls=":", label="intacto")
        if logy:
            ax.set_yscale("log")
            ax.axhline(atac, color="#c00000", ls="--", label="ataque por grau")
        elif lo <= atac <= hi:
            ax.axhline(atac, color="#c00000", ls="--", label="ataque por grau")
        else:
            # valor do ataque fora da escala: anota com seta em vez de linha
            ax.set_ylim(lo, hi)
            ax.annotate(f"ataque por grau: {atac:,.4g}".replace(",", "\u2009"),
                        xy=(1, hi), xytext=(1, hi - 0.12 * (hi - lo)),
                        ha="center", fontsize=8, color="#c00000",
                        fontweight="bold")
            ax.annotate("", xy=(1.0, hi), xytext=(1.0, hi - 0.33 * (hi - lo)),
                        arrowprops=dict(arrowstyle="-|>", color="#c00000", lw=1.4))
        ax.set_title(titulo, fontsize=10); ax.set_xticks([])
        ax.tick_params(axis="y", labelsize=8); ax.grid(alpha=.3, axis="y")
    axs[0].legend(fontsize=7.5, loc="lower left")
    fig.suptitle(f"Robustez a remocao de 5% dos vertices "
                 f"(r={r:,}; T={T} repeticoes aleatorias)".replace(",", "\u2009"),
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig_robustez_box.pdf")); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    ax.plot([f * 100 for f in fracoes], curva["aleatorio"], "o-", color="#1f4e79",
            label="falhas aleatorias")
    ax.plot([f * 100 for f in fracoes], curva["grau"], "s-", color="#c00000",
            label="ataque direcionado (grau)")
    ax.set_xlabel("% de vertices removidos"); ax.set_ylabel("$S$ = maior componente / $n$")
    ax.set_title("Degradacao da componente gigante"); ax.legend(fontsize=8); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_robustez_curva.pdf")); plt.close(fig)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ATAQUE DIRIGIDO POR CENTRALIDADE DE INTERMEDIACAO (BETWEENNESS).

O Apendice C do enunciado recomenda "betweenness ou degree".  O calculo
EXATO da intermediacao e o algoritmo de Brandes, O(n*m), que neste grafo
daria ~1e12 operacoes -- inviavel.  Usa-se entao a versao AMOSTRADA do
proprio algoritmo de Brandes (Brandes & Pich, 2007): em vez de acumular
dependencias a partir dos n vertices, acumula-se a partir de K pivos
sorteados uniformemente e multiplica-se o resultado por n/K.  O estimador
e nao enviesado, e para o nosso fim -- ORDENAR os vertices e remover os
5% maiores -- a precisao do ranking no topo e o que importa.

Implementacao propria, em Python puro sobre o CSR (nenhuma biblioteca de
grafos), em duas fases por pivo:
  1. BFS que calcula distancia d[] e numero de caminhos minimos sigma[];
  2. acumulacao em ordem inversa de distancia, sem armazenar listas de
     predecessores (eles sao redescobertos varrendo os vizinhos w com
     d[w] == d[v] - 1), o que mantem a memoria em O(n).

Saida: resultados/12_betweenness.json  (+ figura comparando os ataques)
"""
import os, sys, json, math, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import igraph as ig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lib_grafo import ler_arestas, compactar, lista_adj, DADOS, RES, FIG, nomes_vertices
from lib_dist import bounding_diameters

random.seed(777); np.random.seed(777)
K_PIVOS = 400
AMOSTRA_L = 150
Z = 1.959963985


# ----------------------------------------------------------------------
def brandes_amostrado(n, indptr, indices, pivos):
    """
    Centralidade de intermediacao estimada pelo algoritmo de Brandes a
    partir de um conjunto de pivos.  Devolve o vetor nao normalizado
    (soma das dependencias), que ja basta para ordenar.
    """
    CB = [0.0] * n
    INF = -1
    for s in pivos:
        d = [INF] * n
        sigma = [0.0] * n
        delta = [0.0] * n
        d[s] = 0
        sigma[s] = 1.0
        fila = [s]
        cab = 0
        ordem = []
        while cab < len(fila):                       # fase 1: BFS
            v = fila[cab]; cab += 1
            ordem.append(v)
            dv = d[v]; sv = sigma[v]
            for i in range(indptr[v], indptr[v + 1]):
                w = indices[i]
                dw = d[w]
                if dw == INF:
                    d[w] = dv + 1
                    sigma[w] = sv
                    fila.append(w)
                elif dw == dv + 1:
                    sigma[w] += sv
        for v in reversed(ordem):                    # fase 2: acumulacao
            dv = d[v]
            if dv == 0:
                continue
            coef = (1.0 + delta[v]) / sigma[v]
            dvm = dv - 1
            for i in range(indptr[v], indptr[v + 1]):
                w = indices[i]
                if d[w] == dvm:
                    delta[w] += sigma[w] * coef
            CB[v] += delta[v]
    return np.array(CB)


def validar_contra_igraph():
    """Com TODOS os vertices como pivos, Brandes proprio == igraph exato."""
    falhas = []
    for t in range(5):
        n0 = random.randint(40, 120)
        g = ig.Graph.Erdos_Renyi(n=n0, m=random.randint(2 * n0, 3 * n0))
        g.simplify()
        g = g.connected_components().giant()
        n = g.vcount()
        e = np.array(g.get_edgelist())
        ip, ix = lista_adj(n, e[:, 0], e[:, 1])
        CB = brandes_amostrado(n, ip.tolist(), ix.tolist(), list(range(n)))
        ref = np.array(g.betweenness())
        if not np.allclose(CB / 2.0, ref, atol=1e-6):
            falhas.append(t)
    return {"grafos_testados": 5, "falhas": falhas,
            "status": "OK" if not falhas else "FALHOU"}


def metricas(g, n_orig, nome):
    comp = g.connected_components()
    tam = np.array(comp.sizes())
    graus = np.array(g.degree())
    gc = comp.giant()
    k = min(AMOSTRA_L, gc.vcount())
    medias = []
    for s in random.sample(range(gc.vcount()), k):
        d = np.array(gc.distances(source=[s])[0], dtype=float)
        d = d[d > 0]
        medias.append(d.mean())
    a = np.array(medias)
    bd = bounding_diameters(gc, verbose=False)
    return {
        "cenario": nome,
        "A_maior_componente": int(tam.max()),
        "A_fracao_do_original": float(tam.max() / n_orig),
        "B_n_componentes": int(len(tam)),
        "C_caminho_medio": float(a.mean()),
        "C_diametro": bd["diametro"],
        "D_fracao_isolados": float((graus == 0).mean()),
    }


def main():
    R = {}
    print(">> validando a implementacao de Brandes contra o igraph", flush=True)
    R["validacao"] = validar_contra_igraph()
    print("   ", R["validacao"], flush=True)

    U, V, P = ler_arestas(os.path.join(DADOS, "arestas_cg.csv"))
    u, v, rot = compactar(U, V)
    n = len(rot)
    indptr, indices = lista_adj(n, u, v)
    indptr, indices = indptr.tolist(), indices.tolist()
    G = ig.Graph(n=n, edges=list(zip(u.tolist(), v.tolist())), directed=False)
    r = int(math.floor(0.05 * n))
    R["n"] = n; R["m"] = G.ecount(); R["r"] = r; R["K_pivos"] = K_PIVOS

    print(f">> Brandes amostrado com K={K_PIVOS} pivos "
          f"(n={n:,}, m={G.ecount():,})", flush=True)
    pivos = random.sample(range(n), K_PIVOS)
    t0 = time.time()
    CB = brandes_amostrado(n, indptr, indices, pivos)
    dt = time.time() - t0
    CB_est = CB * (n / K_PIVOS) / 2.0          # escala e correcao de grafo nao dirigido
    R["tempo_brandes_s"] = dt
    R["tempo_estimado_exato_h"] = dt / K_PIVOS * n / 3600.0
    print(f"   {dt:.0f}s  (o calculo exato levaria "
          f"~{R['tempo_estimado_exato_h']:.0f} h)", flush=True)

    # vetores completos para diagnostico posterior
    np.save(os.path.join(RES, "12_cb_est.npy"), CB_est)
    np.save(os.path.join(RES, "12_pivos.npy"), np.array(pivos, dtype=np.int64))

    nomes = nomes_vertices()
    ordem_bt = np.argsort(-CB_est)
    R["top_betweenness"] = [{"site": nomes[int(rot[i])],
                             "betweenness_estimada": float(CB_est[i]),
                             "grau": int(G.degree(int(i)))}
                            for i in ordem_bt[:15]]

    graus = np.array(G.degree())
    alvo_bt = ordem_bt[:r].tolist()
    alvo_gr = np.argsort(-graus)[:r].tolist()
    R["sobreposicao_com_grau"] = len(set(alvo_bt) & set(alvo_gr)) / r
    R["correlacao_spearman_grau_betweenness"] = float(
        np.corrcoef(np.argsort(np.argsort(-CB_est)),
                    np.argsort(np.argsort(-graus)))[0, 1])

    print(">> ataque por betweenness", flush=True)
    H = G.copy(); H.delete_vertices(alvo_bt)
    R["ataque_betweenness"] = metricas(H, n, "ataque por betweenness")
    print(json.dumps(R["ataque_betweenness"], indent=2), flush=True)

    # referencia: o mesmo numero de vertices removidos por grau
    print(">> ataque por grau (referencia, mesmo r)", flush=True)
    H = G.copy(); H.delete_vertices(alvo_gr)
    R["ataque_grau"] = metricas(H, n, "ataque por grau")

    with open(os.path.join(RES, "12_betweenness.json"), "w", encoding="utf-8") as f:
        json.dump(R, f, indent=2, ensure_ascii=False)

    # ---- figura: grau x betweenness ----
    esc = (n / K_PIVOS) / 2.0
    fig, ax = plt.subplots(figsize=(5.6, 4.0))
    msk = CB_est > 0
    ax.loglog(graus[msk], CB_est[msk], ".", ms=1.5, alpha=.25, color="#1f4e79")
    # as duas faixas horizontais visiveis tem explicacao estrutural:
    ax.axhline(K_PIVOS * esc, ls="--", lw=1, color="#bf8f00")
    ax.axhline(n * esc, ls="--", lw=1, color="#c00000")
    ax.text(0.98, K_PIVOS * esc * 1.6, "uma folha pendente ($C_B$ bruto $=K$)",
            transform=ax.get_yaxis_transform(), ha="right", fontsize=7,
            color="#bf8f00")
    ax.text(0.98, n * esc * 1.6, "porteiro de um pivô ($C_B$ bruto $\\approx n$)",
            transform=ax.get_yaxis_transform(), ha="right", fontsize=7,
            color="#c00000")
    ax.set_xlabel("grau $k$")
    ax.set_ylabel("intermediação estimada $C_B$")
    ax.set_title(f"Grau vs. intermediação ($K={K_PIVOS}$ pivôs)", fontsize=10)
    ax.grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_betweenness.pdf"))
    plt.close(fig)

    print("concluido")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CONTAGEM EXATA DE TRIANGULOS E CLUSTERIZACAO -- IMPLEMENTACAO PROPRIA.

Na primeira versao do relatorio o numero de triangulos era derivado da
transitividade local calculada pelo igraph.  Como o trabalho afirma que
os algoritmos sao de implementacao propria, aqui o valor e recalculado do
zero e confrontado com o da biblioteca.

Algoritmo: para cada aresta (u,v), conta-se |N(u) cap N(v)| percorrendo a
MENOR das duas vizinhancas e testando pertinencia na maior.  O custo e
  sum_{(u,v) in E} min(grau(u), grau(v)),
que neste grafo vale ~4,0e7 -- tratavel em Python puro, ao contrario do
algoritmo "forward" ingenuo (~4,5e9 aqui).

Cada triangulo e descoberto uma vez por aresta, logo
  T = (1/3) * sum_{(u,v)} |N(u) cap N(v)|,
e cada vertice recebe +1 por aresta incidente que participa do triangulo,
logo o numero de triangulos que contem v e tri[v]/2.  Dai sai tambem
  C_local(v) = 2*t_v / (k_v (k_v - 1))   e   C = media de C_local,
  transitividade global = 3T / (numero de triplas conectadas).

Saida: resultados/13_triangulos.json
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import igraph as ig

from lib_grafo import ler_arestas, compactar, lista_adj, DADOS, RES, nomes_vertices


def contar(n, indptr, indices, u, v):
    """Devolve (T, tri_por_vertice) com T exato."""
    viz = [None] * n
    for x in range(n):
        viz[x] = set(indices[indptr[x]:indptr[x + 1]])
    tri = [0] * n
    total = 0
    graus = [indptr[i + 1] - indptr[i] for i in range(n)]
    for i in range(len(u)):
        a = u[i]; b = v[i]
        if graus[a] <= graus[b]:
            menor, Smaior = a, viz[b]
        else:
            menor, Smaior = b, viz[a]
        c = 0
        for w in viz[menor]:
            if w in Smaior:
                c += 1
        if c:
            total += c
            tri[a] += c
            tri[b] += c
    assert total % 3 == 0, "soma das intersecoes deveria ser multipla de 3"
    return total // 3, np.array(tri, dtype=np.int64)


def main():
    U, V, P = ler_arestas(os.path.join(DADOS, "arestas.csv"))
    u, v, rot = compactar(U, V)
    n = len(rot)
    indptr, indices = lista_adj(n, u, v)
    indptr, indices = indptr.tolist(), indices.tolist()
    uu, vv = u.tolist(), v.tolist()
    graus = np.diff(np.array(indptr))

    print(f">> contando triangulos  (n={n:,}, m={len(uu):,})", flush=True)
    t0 = time.time()
    T, tri = contar(n, indptr, indices, uu, vv)
    dt = time.time() - t0
    print(f"   T = {T:,}   ({dt:.0f}s)", flush=True)

    t_v = tri / 2.0                                   # triangulos que contem v
    denom = graus * (graus - 1) / 2.0
    com_denom = denom > 0
    C_local = np.zeros(n)
    C_local[com_denom] = t_v[com_denom] / denom[com_denom]
    C_media = float(C_local.mean())
    triplas = float(denom.sum())
    transitividade = 3.0 * T / triplas

    print(">> conferindo com o igraph", flush=True)
    G = ig.Graph(n=n, edges=list(zip(uu, vv)), directed=False)
    C_ig = float(np.mean(G.transitivity_local_undirected(mode="zero")))
    tr_ig = float(G.transitivity_undirected(mode="zero"))

    nomes = nomes_vertices()
    ordem = np.argsort(-t_v)[:10]

    R = {
        "n_vertices": n, "n_arestas": len(uu),
        "n_triangulos_exato": int(T),
        "custo_do_algoritmo_operacoes": int(
            np.minimum(graus[u], graus[v]).sum()),
        "tempo_s": dt,
        "clusterizacao_media_propria": C_media,
        "clusterizacao_media_igraph": C_ig,
        "diferenca_clusterizacao": abs(C_media - C_ig),
        "transitividade_propria": transitividade,
        "transitividade_igraph": tr_ig,
        "diferenca_transitividade": abs(transitividade - tr_ig),
        "n_triplas_conectadas": triplas,
        "vertices_em_algum_triangulo": int((t_v > 0).sum()),
        "fracao_vertices_em_triangulo": float((t_v > 0).mean()),
        "top_triangulos": [{"site": nomes[int(rot[i])],
                            "triangulos": int(t_v[i]),
                            "grau": int(graus[i])} for i in ordem],
        "confere_com_igraph": bool(abs(C_media - C_ig) < 1e-9
                                   and abs(transitividade - tr_ig) < 1e-9),
    }
    with open(os.path.join(RES, "13_triangulos.json"), "w", encoding="utf-8") as f:
        json.dump(R, f, indent=2, ensure_ascii=False)
    print(json.dumps({k: R[k] for k in R if k != "top_triangulos"},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

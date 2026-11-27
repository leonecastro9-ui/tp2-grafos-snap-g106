#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Validacao das implementacoes proprias (lib_algoritmos.py e lib_dist.py)
contra a biblioteca igraph, em grafos aleatorios pequenos onde a resposta
de referencia pode ser calculada exatamente.

Sao verificados: BFS/DFS (mesmo conjunto alcancado), Dijkstra, Bellman-Ford,
Floyd-Warshall (versoes pura e NumPy), Tarjan (SCC, articulacoes e pontes),
Prim, Kruskal e o calculo exato de diametro/raio (BoundingDiameters).
"""
import os, sys, json, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import igraph as ig
import lib_algoritmos as A
from lib_grafo import lista_adj, RES
from lib_dist import bounding_diameters

random.seed(11); np.random.seed(11)
N_TESTES = 15


def csr_dirigido(n, u, v):
    o = np.argsort(u, kind="stable")
    return (np.concatenate([[0], np.cumsum(np.bincount(u, minlength=n))]).tolist(),
            v[o].tolist())


def main():
    falhas = []
    testes = 0
    for t in range(N_TESTES):
        n0 = random.randint(40, 250)
        m0 = random.randint(2 * n0, 4 * n0)

        # ---- grafo dirigido: Tarjan SCC ----
        gd = ig.Graph.Erdos_Renyi(n=n0, m=m0, directed=True); gd.simplify()
        e = np.array(gd.get_edgelist())
        ip, ix = csr_dirigido(n0, e[:, 0], e[:, 1])
        _, tams = A.tarjan_scc(n0, ip, ix)
        testes += 1
        if sorted(tams) != sorted(gd.connected_components(mode="strong").sizes()):
            falhas.append(f"teste {t}: Tarjan SCC")

        # ---- grafo nao dirigido conexo ----
        gu = ig.Graph.Erdos_Renyi(n=n0, m=m0, directed=False); gu.simplify()
        gu = gu.connected_components().giant()
        n = gu.vcount()
        eu = np.array(gu.get_edgelist()); uu, vv = eu[:, 0], eu[:, 1]
        w = np.random.rand(len(uu)) + 0.1
        ipp, ixx, cst = lista_adj(n, uu, vv, w)
        ipp, ixx, cst = ipp.tolist(), ixx.tolist(), cst.tolist()
        ref = np.array(gu.distances(source=[0], weights=w.tolist())[0])

        testes += 1
        if sorted(A.bfs(n, ipp, ixx, 0)[1]) != sorted(A.dfs(n, ipp, ixx, 0)):
            falhas.append(f"teste {t}: BFS/DFS alcancam conjuntos diferentes")

        testes += 1
        if not np.allclose(np.array(A.dijkstra(n, ipp, ixx, cst, 0)), ref):
            falhas.append(f"teste {t}: Dijkstra")

        testes += 1
        db, _ = A.bellman_ford(n, [(int(uu[i]), int(vv[i]), float(w[i]))
                                   for i in range(len(uu))], 0)
        if not np.allclose(np.array(db), ref):
            falhas.append(f"teste {t}: Bellman-Ford")

        D = [[float("inf")] * n for _ in range(n)]
        for i in range(n):
            D[i][i] = 0.0
        for i in range(len(uu)):
            a, b = int(uu[i]), int(vv[i])
            D[a][b] = D[b][a] = min(D[a][b], float(w[i]))
        testes += 2
        if not np.allclose(np.array(A.floyd_warshall([r[:] for r in D])[0]), ref):
            falhas.append(f"teste {t}: Floyd-Warshall")
        if not np.allclose(A.floyd_warshall_np(np.array(D), np)[0], ref):
            falhas.append(f"teste {t}: Floyd-Warshall (NumPy)")

        testes += 1
        if A.tarjan_articulacoes(n, ipp, ixx) != (len(gu.articulation_points()),
                                                  len(gu.bridges())):
            falhas.append(f"teste {t}: Tarjan articulacoes/pontes")

        mst = sum(w[i] for i in gu.spanning_tree(weights=w.tolist(), return_tree=False))
        tp, ap = A.prim(n, ipp, ixx, cst, 0)
        tk, ak = A.kruskal(n, [(float(w[i]), int(uu[i]), int(vv[i]))
                               for i in range(len(uu))])
        testes += 2
        if abs(tp - mst) > 1e-9 or ap != n - 1:
            falhas.append(f"teste {t}: Prim")
        if abs(tk - mst) > 1e-9 or ak != n - 1:
            falhas.append(f"teste {t}: Kruskal")

        ecc = gu.eccentricity()
        bd = bounding_diameters(gu, verbose=False)
        testes += 1
        if bd["diametro"] != int(max(ecc)) or bd["raio"] != int(min(ecc)):
            falhas.append(f"teste {t}: BoundingDiameters")

        testes += 1
        eul = A.eulerianidade(gu.degree(), True)
        impares = sum(1 for d in gu.degree() if d % 2 == 1)
        if eul["vertices_grau_impar"] != impares:
            falhas.append(f"teste {t}: eulerianidade")

    saida = {"n_grafos_testados": N_TESTES, "n_verificacoes": testes,
             "falhas": falhas, "status": "OK" if not falhas else "FALHOU"}
    with open(os.path.join(RES, "09_validacao.json"), "w", encoding="utf-8") as f:
        json.dump(saida, f, indent=2, ensure_ascii=False)
    print(json.dumps(saida, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

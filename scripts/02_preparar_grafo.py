#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Etapa 2: tratamento final e materializacao dos grafos de trabalho.

Entrada : dados/arestas_dirigidas.csv  (origem,destino,peso)
Saidas  :
  dados/vertices.csv              id,hostname,componente
  dados/arestas.csv               u,v,peso           (nao dirigido, simples)
  dados/arestas_cg.csv            u,v,peso           (componente gigante)
  dados/arestas_dirigidas_id.csv  origem,destino,peso (dirigido, ids inteiros)
  dados/componentes.csv           tamanho,qtd
  resultados/02_preparacao.json

Tratamentos (ver relatorio):
  - dirigido -> nao dirigido (a aresta u-v existe se ha hiperlink em
    qualquer sentido);  peso = soma dos hiperlinks nos dois sentidos;
  - multigrafo -> grafo simples;
  - remocao de auto-loops (ja feita na etapa 1);
  - identificacao das componentes conexas e extracao da maior (CG).
"""
import csv, json, os, sys
from collections import defaultdict

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ENT = os.path.join(BASE, "dados", "arestas_dirigidas.csv")
DD = os.path.join(BASE, "dados")
RD = os.path.join(BASE, "resultados")
os.makedirs(RD, exist_ok=True)


def main():
    und = {}
    n_dir = 0
    reciprocas = 0
    with open(ENT, encoding="utf-8") as f:
        r = csv.reader(f); next(r)
        for a, b, p in r:
            n_dir += 1
            k = (a, b) if a < b else (b, a)
            if k in und:
                reciprocas += 1
                und[k] += int(p)
            else:
                und[k] = int(p)

    nomes = sorted({x for k in und for x in k})
    idx = {n: i for i, n in enumerate(nomes)}
    N = len(nomes)

    adj = [[] for _ in range(N)]
    arestas = []
    for (a, b), p in und.items():
        u, v = idx[a], idx[b]
        arestas.append((u, v, p))
        adj[u].append(v); adj[v].append(u)

    # --- componentes conexas (BFS iterativa propria) ---
    comp = [-1] * N
    tamanhos = []
    for s in range(N):
        if comp[s] != -1:
            continue
        c = len(tamanhos); comp[s] = c; fila = [s]; t = 0
        while fila:
            prox = []
            for u in fila:
                t += 1
                for v in adj[u]:
                    if comp[v] == -1:
                        comp[v] = c; prox.append(v)
            fila = prox
        tamanhos.append(t)

    cg = max(range(len(tamanhos)), key=lambda c: tamanhos[c])

    # --- gravacao ---
    with open(os.path.join(DD, "vertices.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["id", "hostname", "componente"])
        for nome in nomes:
            i = idx[nome]; w.writerow([i, nome, comp[i]])

    with open(os.path.join(DD, "arestas.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["u", "v", "peso"])
        w.writerows(arestas)

    with open(os.path.join(DD, "arestas_cg.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["u", "v", "peso"])
        for u, v, p in arestas:
            if comp[u] == cg:
                w.writerow([u, v, p])

    # arestas DIRIGIDAS com ids inteiros (usadas pelo algoritmo de Tarjan)
    with open(os.path.join(DD, "arestas_dirigidas_id.csv"), "w", newline="",
              encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["origem", "destino", "peso"])
        with open(ENT, encoding="utf-8") as fin:
            r = csv.reader(fin); next(r)
            for a, b, p in r:
                w.writerow([idx[a], idx[b], p])

    dist = defaultdict(int)
    for t in tamanhos:
        dist[t] += 1
    with open(os.path.join(DD, "componentes.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f); w.writerow(["tamanho", "quantidade"])
        for t in sorted(dist, reverse=True):
            w.writerow([t, dist[t]])

    info = {
        "arestas_dirigidas_distintas": n_dir,
        "pares_reciprocos_fundidos": reciprocas,
        "n_vertices": N,
        "n_arestas_nao_dirigidas": len(arestas),
        "n_componentes": len(tamanhos),
        "tamanho_componente_gigante": tamanhos[cg],
        "fracao_componente_gigante": tamanhos[cg] / N,
        "maiores_componentes": sorted(tamanhos, reverse=True)[:10],
        "arestas_componente_gigante": sum(1 for u, v, _ in arestas if comp[u] == cg),
    }
    with open(os.path.join(RD, "02_preparacao.json"), "w", encoding="utf-8") as f:
        json.dump(info, f, indent=2, ensure_ascii=False)
    print(json.dumps(info, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

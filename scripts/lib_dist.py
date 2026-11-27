#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Diametro e raio EXATOS pelo algoritmo BoundingDiameters
(Takes & Kosters, "Determining the diameter of small world networks", 2011),
que mantem limites inferior/superior de excentricidade por vertice e poda
o conjunto candidato, exigindo apenas algumas centenas de BFS em vez de n.
"""
import numpy as np


def bounding_diameters(g, max_bfs=20000, verbose=True):
    """g: igraph.Graph conexo e nao dirigido. Devolve dict com diametro/raio."""
    n = g.vcount()
    INF = 10 ** 9
    ecc_lo = np.zeros(n, dtype=np.int64)
    ecc_hi = np.full(n, INF, dtype=np.int64)
    ativo = np.ones(n, dtype=bool)
    dL, dU = 0, INF            # limites do diametro
    rL, rU = 0, INF            # limites do raio
    graus = np.array(g.degree(), dtype=np.int64)
    nbfs = 0
    alto = True                # alterna entre maior ecc_hi e menor ecc_lo

    while ativo.any() and nbfs < max_bfs and (dL < dU or rL < rU):
        cand = np.flatnonzero(ativo)
        if alto:
            chave = ecc_hi[cand]
            v = cand[np.lexsort((-graus[cand], -chave))[0]]
        else:
            chave = ecc_lo[cand]
            v = cand[np.lexsort((-graus[cand], chave))[0]]
        alto = not alto

        d = np.array(g.distances(source=[int(v)])[0], dtype=np.int64)
        nbfs += 1
        e = int(d.max())
        ativo[v] = False
        ecc_lo[v] = ecc_hi[v] = e

        dL = max(dL, e); dU = min(dU, 2 * e)
        rU = min(rU, e); rL = max(rL, (e + 1) // 2)

        idx = np.flatnonzero(ativo)
        dv = d[idx]
        ecc_lo[idx] = np.maximum(ecc_lo[idx], np.maximum(e - dv, dv))
        ecc_hi[idx] = np.minimum(ecc_hi[idx], e + dv)
        lo, hi = ecc_lo[idx], ecc_hi[idx]
        podar = (lo == hi) | ((hi <= dL) & (2 * lo >= dU) & (lo >= rU) & ((hi + 1) // 2 <= rL))
        ativo[idx[podar]] = False

        if verbose and nbfs % 25 == 0:
            print(f"    BFS {nbfs:5d} | diam in [{dL},{dU}] | raio in [{rL},{rU}] "
                  f"| candidatos {int(ativo.sum()):,}", flush=True)

    # Quando o conjunto candidato se esvazia, nenhum vertice restante pode
    # ter excentricidade maior que dL nem menor que rU: os valores sao exatos.
    exato = (not ativo.any()) or (dL == dU and rL == rU)
    return {"diametro": int(dL), "raio": int(rU),
            "diametro_lb": int(dL), "diametro_ub": int(dU),
            "raio_lb": int(rL), "raio_ub": int(rU),
            "bfs_executadas": int(nbfs),
            "candidatos_restantes": int(ativo.sum()),
            "exato": bool(exato)}

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Funcoes utilitarias compartilhadas pelos scripts do TP2."""
import csv, os, numpy as np

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DADOS = os.path.join(BASE, "dados")
RES = os.path.join(BASE, "resultados")
FIG = os.path.join(BASE, "figuras")


def ler_arestas(arquivo):
    """Le um CSV u,v,peso e devolve (u, v, peso) como arrays numpy."""
    U, V, P = [], [], []
    with open(arquivo, encoding="utf-8") as f:
        r = csv.reader(f); next(r)
        for a, b, p in r:
            U.append(int(a)); V.append(int(b)); P.append(int(p))
    return np.array(U, dtype=np.int64), np.array(V, dtype=np.int64), np.array(P, dtype=np.int64)


def compactar(U, V):
    """Renumera os vertices para 0..n-1. Devolve (u, v, rotulos_originais)."""
    nodes = np.unique(np.concatenate([U, V]))
    remap = np.full(nodes.max() + 1, -1, dtype=np.int64)
    remap[nodes] = np.arange(len(nodes))
    return remap[U], remap[V], nodes


def lista_adj(n, u, v, custo=None):
    """
    Lista de adjacencia em formato CSR.  Devolve (indptr, indices) ou
    (indptr, indices, custos) quando um vetor de custos por aresta e dado.
    Os tres vetores ficam alinhados: custos[i] e o custo da aresta que
    leva a indices[i].
    """
    ii = np.concatenate([u, v])
    jj = np.concatenate([v, u])
    ordem = np.argsort(ii, kind="stable")
    indices = jj[ordem]
    indptr = np.concatenate([[0], np.cumsum(np.bincount(ii, minlength=n))])
    if custo is None:
        return indptr, indices
    cc = np.concatenate([custo, custo])[ordem]
    return indptr, indices, cc


def nomes_vertices():
    """Dicionario id -> hostname."""
    d = {}
    with open(os.path.join(DADOS, "vertices.csv"), encoding="utf-8") as f:
        r = csv.reader(f); next(r)
        for i, h, c in r:
            d[int(i)] = h
    return d

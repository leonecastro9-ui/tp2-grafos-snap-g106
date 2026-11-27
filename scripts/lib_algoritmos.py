#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Implementacoes proprias (Python puro, iterativas) dos algoritmos da
disciplina.  Nenhuma biblioteca de grafos e usada aqui: o grafo e
representado em CSR (indptr, indices, pesos), o que mantem o custo de
acesso a vizinhanca em O(grau).
"""
import heapq
from math import inf


# ---------------------------------------------------------------- BFS
def bfs(n, indptr, indices, s):
    """Busca em largura. Complexidade O(n + m)."""
    dist = [-1] * n
    dist[s] = 0
    fila = [s]
    cabeca = 0
    ordem = []
    while cabeca < len(fila):
        u = fila[cabeca]; cabeca += 1
        ordem.append(u)
        du = dist[u] + 1
        for i in range(indptr[u], indptr[u + 1]):
            v = indices[i]
            if dist[v] < 0:
                dist[v] = du
                fila.append(v)
    return dist, ordem


# ---------------------------------------------------------------- DFS
def dfs(n, indptr, indices, s):
    """Busca em profundidade iterativa. Complexidade O(n + m)."""
    visitado = bytearray(n)
    ponteiro = list(indptr[:n])          # ponteiro[u] = proxima posicao de u no CSR
    pilha = [s]
    visitado[s] = 1
    ordem = [s]
    while pilha:
        u = pilha[-1]
        avancou = False
        i = ponteiro[u]
        fim = indptr[u + 1]
        while i < fim:
            v = indices[i]; i += 1
            if not visitado[v]:
                visitado[v] = 1
                ponteiro[u] = i
                pilha.append(v)
                ordem.append(v)
                avancou = True
                break
        if not avancou:
            ponteiro[u] = i
            pilha.pop()
    return ordem


# ------------------------------------------------------- Eulerianidade
def eulerianidade(graus, conexo_nos_nao_isolados):
    """
    Grafo nao dirigido tem circuito euleriano  <=> conexo (entre os
    vertices de grau > 0) e todos os graus pares.
    Tem trilha euleriana aberta <=> conexo e exatamente 2 vertices impares.
    Complexidade O(n).
    """
    impares = 0
    for g in graus:
        if g & 1:
            impares += 1
    return {
        "vertices_grau_impar": impares,
        "conexo": bool(conexo_nos_nao_isolados),
        "circuito_euleriano": bool(conexo_nos_nao_isolados and impares == 0),
        "trilha_euleriana": bool(conexo_nos_nao_isolados and impares == 2),
    }


# ------------------------------------------------------------ Dijkstra
def dijkstra(n, indptr, indices, custos, s):
    """Caminhos minimos com fila de prioridade. O((n + m) log n)."""
    dist = [inf] * n
    dist[s] = 0.0
    fechado = bytearray(n)
    h = [(0.0, s)]
    while h:
        d, u = heapq.heappop(h)
        if fechado[u]:
            continue
        fechado[u] = 1
        for i in range(indptr[u], indptr[u + 1]):
            v = indices[i]
            if fechado[v]:
                continue
            nd = d + custos[i]
            if nd < dist[v]:
                dist[v] = nd
                heapq.heappush(h, (nd, v))
    return dist


# -------------------------------------------------------- Bellman-Ford
def bellman_ford(n, arestas, s):
    """
    Caminhos minimos tolerando pesos negativos. O(n*m).
    'arestas' e uma lista de (u, v, custo); por ser nao dirigido cada
    aresta e relaxada nos dois sentidos.
    """
    dist = [inf] * n
    dist[s] = 0.0
    for _ in range(n - 1):
        mudou = False
        for u, v, c in arestas:
            du = dist[u]; dv = dist[v]
            if du + c < dv:
                dist[v] = du + c; mudou = True
            elif dv + c < du:
                dist[u] = dv + c; mudou = True
        if not mudou:
            break
    # deteccao de ciclo negativo
    for u, v, c in arestas:
        if dist[u] + c < dist[v] - 1e-12 or dist[v] + c < dist[u] - 1e-12:
            return dist, True
    return dist, False


# ------------------------------------------------------ Floyd-Warshall
def floyd_warshall(D):
    """
    Caminhos minimos entre todos os pares, O(n^3), sobre matriz de
    adjacencia D (lista de listas, inf para ausencia de aresta).
    """
    n = len(D)
    for k in range(n):
        Dk = D[k]
        for i in range(n):
            Di = D[i]
            dik = Di[k]
            if dik == inf:
                continue
            for j in range(n):
                nd = dik + Dk[j]
                if nd < Di[j]:
                    Di[j] = nd
    return D


def floyd_warshall_np(D, np):
    """Mesma recorrencia, vetorizada em NumPy (mesma ordem O(n^3))."""
    n = D.shape[0]
    for k in range(n):
        D = np.minimum(D, D[:, k, None] + D[None, k, :])
    return D


# --------------------------------------------------------------- Tarjan
def tarjan_scc(n, indptr, indices):
    """
    Componentes fortemente conexas de um grafo DIRIGIDO (Tarjan, versao
    iterativa).  Complexidade O(n + m).
    """
    indice = [-1] * n
    low = [0] * n
    na_pilha = bytearray(n)
    pilha = []
    comp = [-1] * n
    contador = 0
    ncomp = 0
    tamanhos = []
    for raiz in range(n):
        if indice[raiz] != -1:
            continue
        trabalho = [(raiz, indptr[raiz])]
        indice[raiz] = low[raiz] = contador; contador += 1
        pilha.append(raiz); na_pilha[raiz] = 1
        while trabalho:
            u, i = trabalho[-1]
            avancou = False
            while i < indptr[u + 1]:
                v = indices[i]; i += 1
                if indice[v] == -1:
                    trabalho[-1] = (u, i)
                    indice[v] = low[v] = contador; contador += 1
                    pilha.append(v); na_pilha[v] = 1
                    trabalho.append((v, indptr[v]))
                    avancou = True
                    break
                elif na_pilha[v] and indice[v] < low[u]:
                    low[u] = indice[v]
            if avancou:
                continue
            trabalho[-1] = (u, i)
            if low[u] == indice[u]:
                t = 0
                while True:
                    w = pilha.pop(); na_pilha[w] = 0
                    comp[w] = ncomp; t += 1
                    if w == u:
                        break
                tamanhos.append(t); ncomp += 1
            trabalho.pop()
            if trabalho:
                p, pi = trabalho[-1]
                if low[u] < low[p]:
                    low[p] = low[u]
    return comp, tamanhos


def tarjan_articulacoes(n, indptr, indices):
    """
    Vertices de articulacao e pontes de um grafo NAO DIRIGIDO
    (DFS com low-link, versao iterativa). O(n + m).
    """
    disc = [-1] * n
    low = [0] * n
    pai = [-1] * n
    art = bytearray(n)
    pontes = 0
    tempo = 0
    for raiz in range(n):
        if disc[raiz] != -1:
            continue
        filhos_raiz = 0
        pilha = [(raiz, indptr[raiz])]
        disc[raiz] = low[raiz] = tempo; tempo += 1
        while pilha:
            u, i = pilha[-1]
            avancou = False
            while i < indptr[u + 1]:
                v = indices[i]; i += 1
                if disc[v] == -1:
                    pilha[-1] = (u, i)
                    pai[v] = u
                    disc[v] = low[v] = tempo; tempo += 1
                    if u == raiz:
                        filhos_raiz += 1
                    pilha.append((v, indptr[v]))
                    avancou = True
                    break
                elif v != pai[u] and disc[v] < low[u]:
                    low[u] = disc[v]
            if avancou:
                continue
            pilha[-1] = (u, i)
            pilha.pop()
            if pilha:
                p = pilha[-1][0]
                if low[u] < low[p]:
                    low[p] = low[u]
                if low[u] > disc[p]:
                    pontes += 1
                if p != raiz and low[u] >= disc[p]:
                    art[p] = 1
        if filhos_raiz > 1:
            art[raiz] = 1
    return int(sum(art)), pontes


# ------------------------------------------------------------------ MST
def prim(n, indptr, indices, custos, s=0):
    """Arvore geradora minima por Prim com heap. O(m log n)."""
    emarvore = bytearray(n)
    total = 0.0
    nar = 0
    h = [(0.0, s)]
    while h:
        c, u = heapq.heappop(h)
        if emarvore[u]:
            continue
        emarvore[u] = 1
        total += c; nar += 1
        for i in range(indptr[u], indptr[u + 1]):
            v = indices[i]
            if not emarvore[v]:
                heapq.heappush(h, (custos[i], v))
    return total, nar - 1


def kruskal(n, arestas):
    """Arvore (floresta) geradora minima por Kruskal. O(m log m)."""
    pai = list(range(n))
    rank = [0] * n

    def find(x):
        r = x
        while pai[r] != r:
            r = pai[r]
        while pai[x] != r:
            pai[x], x = r, pai[x]
        return r

    total = 0.0
    nar = 0
    for c, u, v in sorted(arestas):
        ru, rv = find(u), find(v)
        if ru == rv:
            continue
        if rank[ru] < rank[rv]:
            ru, rv = rv, ru
        pai[rv] = ru
        if rank[ru] == rank[rv]:
            rank[ru] += 1
        total += c; nar += 1
    return total, nar

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PARTE III.1 - Propriedade de mundo pequeno (small-world).

Implementacao PROPRIA (sem bibliotecas de grafos) de:
  * coeficiente de clusterizacao local e medio  C(G);
  * comprimento medio dos caminhos  L(G), estimado por BFS a partir de
    uma amostra aleatoria de fontes, com intervalo de confianca de 95%;
  * dois modelos nulos aleatorios equivalentes:
      (a) Erdos-Renyi G(n, m)  - mesmo n e mesmo m;
      (b) modelo de configuracao - mesma sequencia de graus, obtido por
          trocas duplas de arestas (double edge swap) do proprio grafo.

Criterio (Apendice A):  L(G) ~ L_rand  e  C(G) >> C_rand.
Coeficiente small-world de Humphries & Gurney:  sigma = (C/C_r)/(L/L_r).
"""
import os, sys, json, math, random, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lib_grafo import ler_arestas, compactar, lista_adj, DADOS, RES, FIG

random.seed(123); np.random.seed(123)
AMOSTRA_L = 500        # fontes de BFS para estimar L
AMOSTRA_C = 50000      # vertices amostrados para estimar C
Z = 1.959963985


# ----------------------------------------------------------------------
def csr(n, u, v):
    indptr, indices = lista_adj(n, u, v)
    return indptr.tolist(), indices.tolist()


def bfs_media(n, indptr, indices, s):
    """BFS propria; devolve (soma das distancias, nº de alcancados, ecc)."""
    dist = [-1] * n
    dist[s] = 0
    fila = [s]; cab = 0; soma = 0; alc = 0; ecc = 0
    while cab < len(fila):
        u = fila[cab]; cab += 1
        du = dist[u] + 1
        for i in range(indptr[u], indptr[u + 1]):
            w = indices[i]
            if dist[w] < 0:
                dist[w] = du; soma += du; alc += 1
                if du > ecc:
                    ecc = du
                fila.append(w)
    return soma, alc, ecc


def comprimento_medio(n, indptr, indices, fontes):
    medias = []
    for s in fontes:
        soma, alc, _ = bfs_media(n, indptr, indices, s)
        if alc:
            medias.append(soma / alc)
    a = np.array(medias)
    mu, sd = a.mean(), a.std(ddof=1)
    e = Z * sd / math.sqrt(len(a))
    return {"L": float(mu), "desvio": float(sd), "ic95": [float(mu - e), float(mu + e)],
            "n_fontes": len(a)}


def clusterizacao(n, indptr, indices, amostra=None):
    """
    Coeficiente de clusterizacao local  C_v = 2*t_v / (k_v (k_v - 1))
    e sua media.  t_v e obtido contando, para cada par de vizinhos,
    quantos sao adjacentes entre si (interseccao de listas de adjacencia
    percorrendo sempre a menor delas).
    """
    vizinhos = [None] * n
    alvos = range(n) if amostra is None else amostra
    locais = []
    triangulos_total = 0
    for v in alvos:
        ini, fim = indptr[v], indptr[v + 1]
        k = fim - ini
        if k < 2:
            locais.append(0.0)
            continue
        Sv = vizinhos[v]
        if Sv is None:
            Sv = vizinhos[v] = set(indices[ini:fim])
        t = 0
        for i in range(ini, fim):
            u = indices[i]
            iu, fu = indptr[u], indptr[u + 1]
            if fu - iu <= k:
                for j in range(iu, fu):
                    if indices[j] in Sv:
                        t += 1
            else:
                Su = vizinhos[u]
                if Su is None:
                    Su = vizinhos[u] = set(indices[iu:fu])
                for w in Sv:
                    if w in Su:
                        t += 1
        t //= 2                      # cada triangulo foi contado 2x
        triangulos_total += t
        locais.append(2.0 * t / (k * (k - 1)))
    a = np.array(locais)
    mu, sd = a.mean(), a.std(ddof=1)
    e = Z * sd / math.sqrt(len(a))
    return {"C": float(mu), "desvio": float(sd),
            "ic95": [float(mu - e), float(mu + e)],
            "n_vertices_avaliados": len(a),
            "triangulos_nos_vertices_avaliados": int(triangulos_total)}


def erdos_renyi(n, m, rng):
    """G(n, m) proprio: m pares distintos sorteados sem auto-loops."""
    vistos = set()
    U = np.empty(m, dtype=np.int64); V = np.empty(m, dtype=np.int64)
    i = 0
    while i < m:
        falta = m - i
        a = rng.integers(0, n, size=falta * 2)
        b = rng.integers(0, n, size=falta * 2)
        for x, y in zip(a.tolist(), b.tolist()):
            if x == y:
                continue
            k = (x, y) if x < y else (y, x)
            if k in vistos:
                continue
            vistos.add(k); U[i] = k[0]; V[i] = k[1]; i += 1
            if i == m:
                break
    return U, V


def configuracao_por_trocas(u, v, rng, fator=10):
    """
    Modelo de configuracao: preserva exatamente a sequencia de graus,
    aplicando  fator*m  trocas duplas de arestas (u1-v1, u2-v2) ->
    (u1-v2, u2-v1), rejeitando as que criam auto-loop ou aresta repetida.
    """
    u = u.copy(); v = v.copy()
    m = len(u)
    conjunto = set()
    for i in range(m):
        a, b = int(u[i]), int(v[i])
        conjunto.add((a, b) if a < b else (b, a))
    tentativas = fator * m
    aceitas = 0
    i1s = rng.integers(0, m, size=tentativas)
    i2s = rng.integers(0, m, size=tentativas)
    for t in range(tentativas):
        i, j = int(i1s[t]), int(i2s[t])
        if i == j:
            continue
        a, b = int(u[i]), int(v[i])
        c, d = int(u[j]), int(v[j])
        if a == d or c == b or a == c or b == d:
            continue
        k1 = (a, d) if a < d else (d, a)
        k2 = (c, b) if c < b else (b, c)
        if k1 in conjunto or k2 in conjunto:
            continue
        o1 = (a, b) if a < b else (b, a)
        o2 = (c, d) if c < d else (d, c)
        conjunto.discard(o1); conjunto.discard(o2)
        conjunto.add(k1); conjunto.add(k2)
        v[i] = d; v[j] = b
        aceitas += 1
    return u, v, aceitas, tentativas


# ----------------------------------------------------------------------
def main():
    rng = np.random.default_rng(123)
    R = {}
    U, V, P = ler_arestas(os.path.join(DADOS, "arestas_cg.csv"))
    u, v, rot = compactar(U, V)
    n, m = len(rot), len(u)
    indptr, indices = csr(n, u, v)
    graus = np.array([indptr[i + 1] - indptr[i] for i in range(n)])
    R["grafo"] = {"n": n, "m": m, "grau_medio": float(graus.mean())}
    print(f"n={n:,} m={m:,}", flush=True)

    fontes = random.sample(range(n), AMOSTRA_L)
    amostraC = random.sample(range(n), AMOSTRA_C)

    print(">> C e L do grafo real", flush=True)
    t0 = time.time()
    R["real"] = {"C": clusterizacao(n, indptr, indices, amostraC),
                 "L": comprimento_medio(n, indptr, indices, fontes)}
    print(f"   C={R['real']['C']['C']:.5f}  L={R['real']['L']['L']:.4f} "
          f"({time.time()-t0:.0f}s)", flush=True)

    print(">> modelo nulo 1: Erdos-Renyi G(n,m)", flush=True)
    Ue, Ve = erdos_renyi(n, m, rng)
    ue, ve, rote = compactar(Ue, Ve)
    ne = len(rote)
    ip_e, ix_e = csr(ne, ue, ve)
    fe = random.sample(range(ne), AMOSTRA_L)
    ae = random.sample(range(ne), min(AMOSTRA_C, ne))
    R["er"] = {"n": ne, "m": len(ue),
               "C": clusterizacao(ne, ip_e, ix_e, ae),
               "L": comprimento_medio(ne, ip_e, ix_e, fe)}
    print(f"   C={R['er']['C']['C']:.6f}  L={R['er']['L']['L']:.4f}", flush=True)

    print(">> modelo nulo 2: configuracao (mesma sequencia de graus)", flush=True)
    uc, vc, aceitas, tent = configuracao_por_trocas(u, v, rng, fator=10)
    ip_c, ix_c = csr(n, uc, vc)
    R["configuracao"] = {
        "trocas_aceitas": aceitas, "trocas_tentadas": tent,
        "C": clusterizacao(n, ip_c, ix_c, amostraC),
        "L": comprimento_medio(n, ip_c, ix_c, fontes)}
    print(f"   C={R['configuracao']['C']['C']:.6f}  "
          f"L={R['configuracao']['L']['L']:.4f}", flush=True)

    # ---------- criterios ----------
    C, L = R["real"]["C"]["C"], R["real"]["L"]["L"]
    kmed = float(graus.mean())
    for nome in ("er", "configuracao"):
        Cr, Lr = R[nome]["C"]["C"], R[nome]["L"]["L"]
        R[nome]["razao_C"] = C / Cr if Cr > 0 else float("inf")
        R[nome]["razao_L"] = L / Lr
        R[nome]["sigma_small_world"] = (C / Cr) / (L / Lr) if Cr > 0 else float("inf")
    R["analitico"] = {
        "C_rand_esperado_<k>/n": kmed / n,
        "L_rand_esperado_ln(n)/ln(<k>)": math.log(n) / math.log(kmed),
    }
    R["conclusao"] = {
        "L_proximo_de_L_rand": abs(L - R["er"]["L"]["L"]) / R["er"]["L"]["L"] < 0.5,
        "C_muito_maior_que_C_rand": C > 10 * R["er"]["C"]["C"],
    }

    with open(os.path.join(RES, "05_small_world.json"), "w", encoding="utf-8") as f:
        json.dump(R, f, indent=2, ensure_ascii=False)

    # ---------- figura ----------
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.6))
    nomes = ["G106 (real)", "Erdos-Renyi", "Configuracao"]
    Cs = [C, R["er"]["C"]["C"], R["configuracao"]["C"]["C"]]
    Ls = [L, R["er"]["L"]["L"], R["configuracao"]["L"]["L"]]
    cores = ["#1f4e79", "#c00000", "#375623"]
    ax[0].bar(nomes, Cs, color=cores); ax[0].set_yscale("log")
    ax[0].set_ylabel("$C$ (clusterizacao media)"); ax[0].grid(alpha=.3, axis="y")
    for i, c in enumerate(Cs):
        ax[0].text(i, c, f"{c:.2e}", ha="center", va="bottom", fontsize=7)
    ax[1].bar(nomes, Ls, color=cores)
    ax[1].set_ylabel("$L$ (caminho medio)"); ax[1].grid(alpha=.3, axis="y")
    for i, l in enumerate(Ls):
        ax[1].text(i, l, f"{l:.2f}", ha="center", va="bottom", fontsize=7)
    for a in ax:
        a.tick_params(axis="x", labelsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_small_world.pdf")); plt.close(fig)
    print(json.dumps({k: R[k].get("sigma_small_world") for k in ("er", "configuracao")},
                     indent=2))


if __name__ == "__main__":
    main()

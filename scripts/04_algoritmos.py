#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PARTE II - Algoritmos da disciplina.

Executa implementacoes PROPRIAS de BFS, DFS, verificacao de eulerianidade,
Dijkstra, Bellman-Ford, Floyd-Warshall, Tarjan (componentes fortemente
conexas e vertices de articulacao), Prim e Kruskal.

Para cada algoritmo mede-se o tempo real de execucao em T repeticoes e
reporta-se media, desvio padrao e intervalo de confianca de 95%
(Normal/Z quando T >= 30; t-Student quando T < 30).

Alem da tabela principal, e feito um estudo de escalabilidade sobre
subgrafos induzidos de tamanho crescente, do qual se estima o expoente
empirico b em  t ~ x^b,  comparado ao expoente teorico.
"""
import os, sys, json, time, math, random
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lib_grafo import ler_arestas, compactar, lista_adj, DADOS, RES, FIG
import lib_algoritmos as A

random.seed(7); np.random.seed(7)
sys.setrecursionlimit(10000)

T_PADRAO = 30          # repeticoes -> IC pela Normal
T_CARO = 10            # repeticoes para algoritmos caros -> IC pela t-Student

# valores criticos
Z_975 = 1.959963985
T_975 = {5: 2.776, 10: 2.262, 15: 2.145, 20: 2.093, 25: 2.064, 29: 2.045}


def estatisticas(tempos):
    """Media, desvio padrao amostral e IC 95% (Z se n>=30, t caso contrario)."""
    n = len(tempos)
    x = np.array(tempos, dtype=float)
    mu = x.mean()
    sd = x.std(ddof=1) if n > 1 else 0.0
    if n >= 30:
        crit, dist = Z_975, "Normal (Z)"
    else:
        crit, dist = T_975.get(n, 2.262), "t-Student"
    e = crit * sd / math.sqrt(n)
    return {"n_repeticoes": n, "media_s": mu, "desvio_s": sd,
            "ic95": [mu - e, mu + e], "distribuicao_ic": dist,
            "min_s": float(x.min()), "max_s": float(x.max())}


def cronometrar(func, reps, args_por_rep=None):
    tempos = []
    saida = None
    for i in range(reps):
        a = args_por_rep[i] if args_por_rep else ()
        t0 = time.perf_counter()
        saida = func(*a)
        tempos.append(time.perf_counter() - t0)
    return estatisticas(tempos), saida


# --------------------------------------------------------------------
class Grafo:
    """Grafo nao dirigido em CSR, com custos = 1/peso."""

    def __init__(self, U, V, P):
        u, v, self.rot = compactar(U, V)
        self.n = len(self.rot)
        self.u, self.v = u, v
        self.m = len(u)
        self.peso = P
        custo_arr = 1.0 / P
        indptr, indices, custos = lista_adj(self.n, u, v, custo_arr)
        self.indptr = indptr.tolist()
        self.indices = indices.tolist()
        self.custos = custos.tolist()
        custo = custo_arr.tolist()
        uu, vv = u.tolist(), v.tolist()
        self.arestas_kruskal = [(custo[i], uu[i], vv[i]) for i in range(self.m)]
        self.arestas_bf = [(uu[i], vv[i], custo[i]) for i in range(self.m)]
        self.graus = [self.indptr[i + 1] - self.indptr[i] for i in range(self.n)]

    def subgrafo_bfs(self, k, semente=0):
        """Subgrafo induzido pelos k primeiros vertices visitados em BFS."""
        _, ordem = A.bfs(self.n, self.indptr, self.indices, semente)
        dentro = np.zeros(self.n, dtype=bool)
        dentro[np.array(ordem[:k], dtype=np.int64)] = True
        idx = np.flatnonzero(dentro[self.u] & dentro[self.v])
        return Grafo(self.u[idx], self.v[idx], self.peso[idx])

    def matriz(self):
        inf = math.inf
        D = [[inf] * self.n for _ in range(self.n)]
        for i in range(self.n):
            D[i][i] = 0.0
        for i in range(self.m):
            a, b = int(self.u[i]), int(self.v[i])
            c = 1.0 / float(self.peso[i])
            if c < D[a][b]:
                D[a][b] = D[b][a] = c
        return D


# --------------------------------------------------------------------
def main():
    R = {"config": {"T_padrao": T_PADRAO, "T_caro": T_CARO}}

    print(">> montando grafos", flush=True)
    U, V, P = ler_arestas(os.path.join(DADOS, "arestas_cg.csv"))
    G = Grafo(U, V, P)
    print(f"   componente gigante: n={G.n:,}  m={G.m:,}", flush=True)
    R["grafo"] = {"n": G.n, "m": G.m, "custo": "1/peso (peso = nº de hiperlinks)"}

    fontes = [random.randrange(G.n) for _ in range(T_PADRAO)]

    # ---------------- BFS ----------------
    print(">> BFS", flush=True)
    st, _ = cronometrar(lambda s: A.bfs(G.n, G.indptr, G.indices, s),
                        T_PADRAO, [(s,) for s in fontes])
    dist, _ = A.bfs(G.n, G.indptr, G.indices, fontes[0])
    st["complexidade_teorica"] = "O(n + m)"
    st["observacao"] = (f"a partir do vertice {fontes[0]} alcanca "
                        f"{sum(1 for d in dist if d >= 0):,} vertices; "
                        f"excentricidade = {max(dist)}")
    R["BFS"] = st

    # ---------------- DFS ----------------
    print(">> DFS", flush=True)
    st, ordem = cronometrar(lambda s: A.dfs(G.n, G.indptr, G.indices, s),
                            T_PADRAO, [(s,) for s in fontes])
    st["complexidade_teorica"] = "O(n + m)"
    st["observacao"] = f"percorre {len(ordem):,} vertices (versao iterativa)"
    R["DFS"] = st

    # ---------------- Eulerianidade ----------------
    print(">> Eulerianidade", flush=True)
    st, eul = cronometrar(lambda: A.eulerianidade(G.graus, True), T_PADRAO)
    st["complexidade_teorica"] = "O(n)"
    st["resultado"] = eul
    R["Eulerianidade"] = st

    # ---------------- Dijkstra ----------------
    print(">> Dijkstra", flush=True)
    st, _ = cronometrar(lambda s: A.dijkstra(G.n, G.indptr, G.indices, G.custos, s),
                        T_PADRAO, [(s,) for s in fontes])
    st["complexidade_teorica"] = "O((n + m) log n)"
    R["Dijkstra"] = st

    # ---------------- Prim ----------------
    print(">> Prim", flush=True)
    st, (custo_prim, nar_prim) = cronometrar(
        lambda s: A.prim(G.n, G.indptr, G.indices, G.custos, s),
        T_PADRAO, [(s,) for s in fontes])
    st["complexidade_teorica"] = "O(m log n)"
    st["resultado"] = {"custo_total": custo_prim, "arestas_na_arvore": nar_prim}
    R["Prim"] = st

    # ---------------- Kruskal ----------------
    print(">> Kruskal", flush=True)
    st, (custo_kr, nar_kr) = cronometrar(lambda: A.kruskal(G.n, G.arestas_kruskal), T_CARO)
    st["complexidade_teorica"] = "O(m log m)"
    st["resultado"] = {"custo_total": custo_kr, "arestas_na_arvore": nar_kr}
    R["Kruskal"] = st

    # ---------------- Tarjan: articulacoes/pontes ----------------
    print(">> Tarjan (articulacoes e pontes)", flush=True)
    st, (nart, npontes) = cronometrar(
        lambda: A.tarjan_articulacoes(G.n, G.indptr, G.indices), T_CARO)
    st["complexidade_teorica"] = "O(n + m)"
    st["resultado"] = {"vertices_de_articulacao": nart, "pontes": npontes,
                       "fracao_articulacoes": nart / G.n}
    R["Tarjan_articulacoes"] = st

    # ---------------- Tarjan: SCC no grafo dirigido ----------------
    print(">> Tarjan (componentes fortemente conexas)", flush=True)
    Ud, Vd, Pd = [], [], []
    import csv as _csv
    with open(os.path.join(DADOS, "arestas_dirigidas_id.csv"), encoding="utf-8") as f:
        r = _csv.reader(f); next(r)
        for a, b, p in r:
            Ud.append(int(a)); Vd.append(int(b)); Pd.append(int(p))
    Ud = np.array(Ud); Vd = np.array(Vd)
    nodes = np.unique(np.concatenate([Ud, Vd]))
    remap = np.full(nodes.max() + 1, -1, dtype=np.int64); remap[nodes] = np.arange(len(nodes))
    ud, vd = remap[Ud], remap[Vd]
    nd = len(nodes)
    ordem_ar = np.argsort(ud, kind="stable")
    grau_s = np.bincount(ud, minlength=nd)
    indptr_d = np.concatenate([[0], np.cumsum(grau_s)]).tolist()
    indices_d = vd[ordem_ar].tolist()
    st, (comp, tams) = cronometrar(lambda: A.tarjan_scc(nd, indptr_d, indices_d), T_CARO)
    st["complexidade_teorica"] = "O(n + m)"
    st["resultado"] = {"n_vertices_dirigido": nd, "m_arestas_dirigido": len(ud),
                       "n_scc": len(tams),
                       "maior_scc": int(max(tams)),
                       "fracao_maior_scc": max(tams) / nd,
                       "scc_triviais": int(sum(1 for t in tams if t == 1))}
    R["Tarjan_SCC"] = st

    # ---------------- Bellman-Ford (subgrafo) ----------------
    print(">> Bellman-Ford (subgrafo)", flush=True)
    SUB_BF = 3000
    g_bf = G.subgrafo_bfs(SUB_BF)
    st, (_, neg) = cronometrar(lambda: A.bellman_ford(g_bf.n, g_bf.arestas_bf, 0), T_CARO)
    st["complexidade_teorica"] = "O(n·m)"
    st["observacao"] = (f"executado sobre subgrafo induzido com n={g_bf.n:,} e "
                        f"m={g_bf.m:,}; inviavel na componente gigante inteira "
                        f"(n·m ≈ {G.n*G.m:.2e} relaxamentos)")
    st["resultado"] = {"ciclo_negativo": neg, "n": g_bf.n, "m": g_bf.m}
    R["Bellman_Ford"] = st

    # ---------------- Floyd-Warshall (subgrafo) ----------------
    print(">> Floyd-Warshall (subgrafo)", flush=True)
    SUB_FW = 400
    g_fw = G.subgrafo_bfs(SUB_FW)
    st, _ = cronometrar(lambda: A.floyd_warshall(g_fw.matriz()), T_CARO)
    st["complexidade_teorica"] = "O(n^3)"
    st["observacao"] = (f"executado sobre subgrafo induzido com n={g_fw.n:,} e "
                        f"m={g_fw.m:,}; na componente gigante seriam "
                        f"{G.n**3:.2e} operacoes e O(n^2) = "
                        f"{G.n**2*8/1e12:.1f} TB de memoria")
    R["Floyd_Warshall"] = st

    D = np.array(g_fw.matriz(), dtype=float)
    st, Dn = cronometrar(lambda: A.floyd_warshall_np(D.copy(), np), T_CARO)
    st["complexidade_teorica"] = "O(n^3) (vetorizado em NumPy)"
    finito = Dn[np.isfinite(Dn) & (Dn > 0)]
    st["resultado"] = {"n": g_fw.n, "distancia_media_ponderada": float(finito.mean()),
                       "maior_distancia_ponderada": float(finito.max())}
    R["Floyd_Warshall_NumPy"] = st

    # =============== estudo de escalabilidade ===============
    print(">> estudo de escalabilidade", flush=True)
    esc = {}
    tamanhos = [2000, 5000, 10000, 25000, 50000, 100000, 200000, 400000, G.n]
    subs = []
    for k in tamanhos:
        sg = G if k >= G.n else G.subgrafo_bfs(k)
        subs.append(sg)
        print(f"   subgrafo n={sg.n:,} m={sg.m:,}", flush=True)

    def medir(nome, f, reps=5, lista=None):
        linhas = []
        for sg in (lista or subs):
            ts = []
            for _ in range(reps):
                t0 = time.perf_counter(); f(sg); ts.append(time.perf_counter() - t0)
            linhas.append({"n": sg.n, "m": sg.m, "tempo_medio_s": float(np.mean(ts))})
            print(f"     {nome} n={sg.n:,} -> {np.mean(ts):.4f}s", flush=True)
        esc[nome] = linhas

    medir("BFS", lambda g: A.bfs(g.n, g.indptr, g.indices, 0))
    medir("DFS", lambda g: A.dfs(g.n, g.indptr, g.indices, 0))
    medir("Dijkstra", lambda g: A.dijkstra(g.n, g.indptr, g.indices, g.custos, 0), reps=3)
    medir("Prim", lambda g: A.prim(g.n, g.indptr, g.indices, g.custos, 0), reps=3)
    medir("Kruskal", lambda g: A.kruskal(g.n, g.arestas_kruskal), reps=3)
    medir("Tarjan_articulacoes", lambda g: A.tarjan_articulacoes(g.n, g.indptr, g.indices), reps=3)

    subs_bf = [G.subgrafo_bfs(k) for k in (250, 500, 1000, 2000, 4000)]
    medir("Bellman_Ford", lambda g: A.bellman_ford(g.n, g.arestas_bf, 0), reps=3, lista=subs_bf)
    subs_fw = [G.subgrafo_bfs(k) for k in (50, 100, 200, 300, 400)]
    medir("Floyd_Warshall", lambda g: A.floyd_warshall(g.matriz()), reps=3, lista=subs_fw)

    # expoentes empiricos
    def expoente(linhas, chave):
        x = np.array([chave(l) for l in linhas], dtype=float)
        y = np.array([l["tempo_medio_s"] for l in linhas], dtype=float)
        ok = (x > 0) & (y > 0)
        b, a = np.polyfit(np.log(x[ok]), np.log(y[ok]), 1)
        return float(b)

    R["escalabilidade"] = esc
    R["expoentes_empiricos"] = {
        "BFS_em_(n+m)": expoente(esc["BFS"], lambda l: l["n"] + l["m"]),
        "DFS_em_(n+m)": expoente(esc["DFS"], lambda l: l["n"] + l["m"]),
        "Dijkstra_em_(n+m)logn": expoente(esc["Dijkstra"], lambda l: (l["n"] + l["m"]) * math.log2(max(l["n"], 2))),
        "Prim_em_mlogn": expoente(esc["Prim"], lambda l: l["m"] * math.log2(max(l["n"], 2))),
        "Kruskal_em_mlogm": expoente(esc["Kruskal"], lambda l: max(l["m"], 2) * math.log2(max(l["m"], 2))),
        "Tarjan_em_(n+m)": expoente(esc["Tarjan_articulacoes"], lambda l: l["n"] + l["m"]),
        "Bellman_Ford_em_nm": expoente(esc["Bellman_Ford"], lambda l: l["n"] * l["m"]),
        "Floyd_Warshall_em_n3": expoente(esc["Floyd_Warshall"], lambda l: l["n"] ** 3),
    }

    with open(os.path.join(RES, "04_algoritmos.json"), "w", encoding="utf-8") as f:
        json.dump(R, f, indent=2, ensure_ascii=False, default=str)

    # ------- figura de escalabilidade -------
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    for nome in ["BFS", "DFS", "Dijkstra", "Prim", "Kruskal", "Tarjan_articulacoes"]:
        L = esc[nome]
        ax[0].loglog([l["n"] + l["m"] for l in L], [l["tempo_medio_s"] for l in L],
                     "o-", ms=4, label=nome.replace("_", " "))
    ax[0].set_xlabel("$n + m$"); ax[0].set_ylabel("tempo (s)")
    ax[0].set_title("Algoritmos de custo quase linear"); ax[0].grid(alpha=.3, which="both")
    ax[0].legend(fontsize=7)
    for nome, xk in [("Bellman_Ford", lambda l: l["n"] * l["m"]),
                     ("Floyd_Warshall", lambda l: l["n"] ** 3)]:
        L = esc[nome]
        ax[1].loglog([xk(l) for l in L], [l["tempo_medio_s"] for l in L],
                     "s-", ms=4, label=nome.replace("_", " "))
    ax[1].set_xlabel("$n\\cdot m$  /  $n^3$"); ax[1].set_ylabel("tempo (s)")
    ax[1].set_title("Algoritmos de custo superlinear"); ax[1].grid(alpha=.3, which="both")
    ax[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_escalabilidade.pdf")); plt.close(fig)

    print(json.dumps(R["expoentes_empiricos"], indent=2))


if __name__ == "__main__":
    main()

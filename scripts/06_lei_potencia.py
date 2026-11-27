#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PARTE III.2 - A distribuicao de graus segue lei de potencia?

Procedimento de Clauset, Shalizi & Newman (2009):
  1. para cada k_min candidato, estima-se gamma por maxima verossimilhanca
     (estimador discreto aproximado);
  2. escolhe-se o k_min que minimiza a distancia de Kolmogorov-Smirnov (D)
     entre a CCDF empirica e a teorica;
  3. avalia-se a qualidade do ajuste por um teste de bondade bootstrap,
     obtendo um valor-p (p >= 0.10 => a lei de potencia nao e rejeitada);
  4. compara-se com um ajuste log-normal e com uma exponencial por razao
     de verossimilhanca (teste de Vuong).
Tambem e feito o ajuste ingenuo por minimos quadrados no grafico log-log,
apenas para contraste com o estimador de maxima verossimilhanca.
"""
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lib_grafo import ler_arestas, compactar, lista_adj, DADOS, RES, FIG

rng = np.random.default_rng(2024)
N_BOOT = 200


def gamma_mle(x, kmin):
    """Estimador de MV (aproximacao continua com correcao de -0.5)."""
    xt = x[x >= kmin]
    if len(xt) < 10:
        return None, 0
    s = np.log(xt / (kmin - 0.5)).sum()
    return 1.0 + len(xt) / s, len(xt)


def ks_dist(x, kmin, g):
    xt = np.sort(x[x >= kmin])
    nt = len(xt)
    emp = np.arange(nt) / nt                     # CDF empirica
    teo = 1.0 - (xt / (kmin - 0.5)) ** (1.0 - g)  # CDF teorica
    return float(np.max(np.abs(emp - teo)))


def ajusta(x, kmins):
    melhor = None
    curva = []
    for kmin in kmins:
        g, nt = gamma_mle(x, kmin)
        if g is None:
            continue
        D = ks_dist(x, kmin, g)
        curva.append((int(kmin), float(g), D, int(nt)))
        if melhor is None or D < melhor[2]:
            melhor = (int(kmin), float(g), D, int(nt))
    return melhor, curva


def amostra_powerlaw(nt, kmin, g, rng):
    """Amostra discreta de lei de potencia por inversao (aprox. continua)."""
    r = rng.random(nt)
    return np.floor((kmin - 0.5) * (1 - r) ** (-1.0 / (g - 1.0)) + 0.5)


def main():
    U, V, P = ler_arestas(os.path.join(DADOS, "arestas.csv"))
    u, v, rot = compactar(U, V)
    n = len(rot)
    indptr, _ = lista_adj(n, u, v)
    graus = np.diff(indptr).astype(float)
    x = graus[graus > 0]
    R = {"n_vertices": int(n), "n_vertices_grau_positivo": int(len(x)),
         "grau_medio": float(graus.mean()), "grau_max": int(graus.max())}

    kmins = np.unique(np.concatenate([np.arange(1, 51),
                                      np.unique(np.round(np.logspace(np.log10(51), np.log10(500), 40)))]))
    (kmin, g, D, nt), curva = ajusta(x, kmins)
    R["mle"] = {"k_min": kmin, "gamma": g, "ks_D": D, "n_cauda": nt,
                "fracao_cauda": nt / len(x),
                "erro_padrao_gamma": (g - 1) / math.sqrt(nt)}
    print(f"MLE: k_min={kmin} gamma={g:.4f} D={D:.5f} n_cauda={nt:,}", flush=True)

    # ---- teste de bondade do ajuste (bootstrap semi-parametrico) ----
    corpo = x[x < kmin]
    p_frac = nt / len(x)
    maiores = 0
    for b in range(N_BOOT):
        nb = len(x)
        esc = rng.random(nb) < p_frac
        sint = np.empty(nb)
        n_cauda = int(esc.sum())
        sint[esc] = amostra_powerlaw(n_cauda, kmin, g, rng)
        if len(corpo) > 0:
            sint[~esc] = rng.choice(corpo, size=nb - n_cauda, replace=True)
        else:
            sint[~esc] = kmin
        m2, _ = ajusta(sint, kmins[kmins <= max(2, kmin * 3)])
        if m2 and m2[2] >= D:
            maiores += 1
        if (b + 1) % 50 == 0:
            print(f"   bootstrap {b+1}/{N_BOOT}", flush=True)
    R["bondade_ajuste"] = {"n_bootstrap": N_BOOT, "p_valor": maiores / N_BOOT,
                           "criterio": "p >= 0.10 => lei de potencia NAO e rejeitada"}
    print("p-valor:", maiores / N_BOOT, flush=True)

    # ---- comparacao com log-normal e exponencial (razao de verossimilhanca) ----
    # Ambas as alternativas sao ajustadas CONDICIONADAS a x >= k_min (ou seja,
    # truncadas e renormalizadas), para que o teste de Vuong seja justo.
    from scipy.optimize import minimize
    from scipy.special import erfc as _erfc
    xt = x[x >= kmin]
    nt_ = len(xt)
    lx = np.log(xt)

    def ll_powerlaw_vec(gg):
        return np.log(gg - 1) - np.log(kmin - 0.5) - gg * np.log(xt / (kmin - 0.5))

    def ll_lognormal_vec(par):
        mu_, s_ = par[0], abs(par[1]) + 1e-12
        z = (math.log(kmin) - mu_) / (s_ * math.sqrt(2))
        cauda = 0.5 * _erfc(z)                      # P(X >= kmin)
        return (-np.log(xt * s_ * math.sqrt(2 * math.pi))
                - (lx - mu_) ** 2 / (2 * s_ ** 2) - math.log(max(cauda, 1e-300)))

    # Observacao importante: sem restricao, o MLE da log-normal truncada
    # diverge para mu -> -inf e sigma -> +inf, limite em que ela CONVERGE
    # para a propria lei de potencia (os modelos deixam de ser distintos e o
    # teste de Vuong perde validade).  Reportamos o ajuste restrito a uma
    # regiao log-normal "legitima" e registramos tambem o valor irrestrito.
    res_livre = minimize(lambda par: -ll_lognormal_vec(par).sum(),
                         x0=[lx.mean(), lx.std(ddof=0)], method="Nelder-Mead")
    res = minimize(lambda par: -ll_lognormal_vec(par).sum(),
                   x0=[lx.mean(), lx.std(ddof=0)], method="L-BFGS-B",
                   bounds=[(-20.0, 20.0), (0.01, 5.0)])
    v_ln = ll_lognormal_vec(res.x)
    lam = 1.0 / (xt - kmin).mean() if (xt - kmin).mean() > 0 else 1.0
    v_ex = np.log(lam) - lam * (xt - kmin)
    v_pl = ll_powerlaw_vec(g)
    ll_pl, ll_ln, ll_ex = float(v_pl.sum()), float(v_ln.sum()), float(v_ex.sum())

    def vuong(v1, v2):
        d = v1 - v2
        Rr = float(d.sum())
        sg = float(d.std(ddof=0))
        if sg == 0:
            return Rr, 1.0
        Zv = Rr / (sg * math.sqrt(len(d)))
        return Rr, float(math.erfc(abs(Zv) / math.sqrt(2)))

    r_ln, p_ln = vuong(v_pl, v_ln)
    r_ex, p_ex = vuong(v_pl, v_ex)
    R["comparacao_modelos"] = {
        "logverossimilhanca_potencia": ll_pl,
        "logverossimilhanca_lognormal_truncada": ll_ln,
        "logverossimilhanca_exponencial_truncada": ll_ex,
        "lognormal_mu": float(res.x[0]), "lognormal_sigma": float(abs(res.x[1])),
        "lognormal_irrestrita_mu": float(res_livre.x[0]),
        "lognormal_irrestrita_sigma": float(abs(res_livre.x[1])),
        "lognormal_irrestrita_loglik": float(-res_livre.fun),
        "nota_degenerescencia": ("o MLE irrestrito da log-normal foge para "
                                 "mu -> -inf, sigma -> +inf, limite em que ela "
                                 "coincide com a lei de potencia; nesse regime "
                                 "os modelos sao indistinguiveis"),
        "exponencial_lambda": float(lam),
        "LR_potencia_vs_lognormal": r_ln, "p_vuong_lognormal": p_ln,
        "LR_potencia_vs_exponencial": r_ex, "p_vuong_exponencial": p_ex,
        "leitura": ("LR > 0 favorece a lei de potencia; LR < 0 favorece a "
                    "alternativa; p < 0.05 torna o sinal significativo"),
    }

    # ---- ajuste ingenuo por minimos quadrados no log-log ----
    vals, cnt = np.unique(x, return_counts=True)
    pk = cnt / len(x)
    msk = vals >= kmin
    b_ols, a_ols = np.polyfit(np.log(vals[msk]), np.log(pk[msk]), 1)
    ss = 1 - np.var(np.log(pk[msk]) - (a_ols + b_ols * np.log(vals[msk]))) / np.var(np.log(pk[msk]))
    R["ols_loglog"] = {"gamma": float(-b_ols), "R2": float(ss),
                       "observacao": "estimador enviesado; apresentado so para contraste"}

    R["curva_kmin"] = [{"k_min": c[0], "gamma": c[1], "D": c[2], "n_cauda": c[3]}
                       for c in curva]
    with open(os.path.join(RES, "06_lei_potencia.json"), "w", encoding="utf-8") as f:
        json.dump(R, f, indent=2, ensure_ascii=False)

    # ---- figuras ----
    xs = np.sort(x)
    ccdf = 1.0 - np.arange(len(xs)) / len(xs)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].loglog(vals, pk, "o", ms=3, alpha=.5, color="#1f4e79", label="dados")
    kk = np.logspace(np.log10(kmin), np.log10(vals.max()), 50)
    norm = pk[np.argmin(np.abs(vals - kmin))]
    ax[0].loglog(kk, norm * (kk / kmin) ** (-g), "-", lw=2, color="#c00000",
                 label=f"lei de potencia $\\gamma$={g:.2f}")
    ax[0].axvline(kmin, ls=":", color="k", lw=1)
    ax[0].set_xlabel("grau $k$"); ax[0].set_ylabel("$P(k)$")
    ax[0].set_title("Distribuicao de graus"); ax[0].legend(fontsize=8)
    ax[0].grid(alpha=.3, which="both")

    ax[1].loglog(xs, ccdf, "-", lw=1.5, color="#1f4e79", label="CCDF empirica")
    msk2 = xs >= kmin
    ccdf_kmin = ccdf[np.argmax(msk2)]
    ax[1].loglog(xs[msk2], ccdf_kmin * (xs[msk2] / kmin) ** (1 - g), "--", lw=2,
                 color="#c00000", label=f"ajuste MV ($k_{{min}}$={kmin})")
    ax[1].axvline(kmin, ls=":", color="k", lw=1)
    ax[1].set_xlabel("grau $k$"); ax[1].set_ylabel("$P(K \\geq k)$")
    ax[1].set_title("CCDF e ajuste de maxima verossimilhanca")
    ax[1].legend(fontsize=8); ax[1].grid(alpha=.3, which="both")
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_lei_potencia.pdf")); plt.close(fig)

    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    cc = np.array(curva, dtype=object)
    ax.plot([c[0] for c in curva], [c[2] for c in curva], "-", color="#1f4e79")
    ax.axvline(kmin, ls="--", color="#c00000", label=f"$k_{{min}}$ otimo = {kmin}")
    ax.set_xscale("log"); ax.set_xlabel("$k_{min}$"); ax.set_ylabel("distancia KS $D$")
    ax.set_title("Escolha de $k_{min}$"); ax.legend(fontsize=8); ax.grid(alpha=.3)
    fig.tight_layout(); fig.savefig(os.path.join(FIG, "fig_kmin.pdf")); plt.close(fig)
    print(json.dumps({k: R[k] for k in ("mle", "bondade_ajuste", "comparacao_modelos", "ols_loglog")},
                     indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

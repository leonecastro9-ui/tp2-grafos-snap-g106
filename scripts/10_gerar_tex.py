#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Converte os JSON de resultados em macros e tabelas LaTeX consumidas pelo
relatorio (relatorio/numeros.tex e relatorio/tab_*.tex).
"""
import os, sys, json, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_grafo import RES, BASE

REL = os.path.join(BASE, "relatorio")
os.makedirs(REL, exist_ok=True)
macros = []


def J(nome):
    with open(os.path.join(RES, nome), encoding="utf-8") as f:
        return json.load(f)


def mil(x):
    """Separador de milhar em padrao brasileiro (espaco fino)."""
    return f"{int(round(x)):,}".replace(",", "\\,")


def dec(x, c=2):
    return f"{x:.{c}f}".replace(".", ",")


def sci(x, c=2):
    s = f"{x:.{c}e}"
    m, e = s.split("e")
    return f"{m.replace('.', ',')} \\times 10^{{{int(e)}}}"


def M(nome, valor):
    macros.append(f"\\newcommand{{\\{nome}}}{{{valor}}}")


# ====================== dados brutos / tratamento ======================
prep = J("02_preparacao.json")
M("recorteArquivo", "\\texttt{quotes\\_2008-08.txt.gz}")
M("docsLidos", mil(7340810))
M("docsRecorte", mil(7340805))
M("hiperlinksRecorte", mil(27810398))
M("autoLoops", mil(18626762))
M("urlsInvalidas", mil(8804))
M("arestasDirigidas", mil(prep["arestas_dirigidas_distintas"]))
M("paresReciprocos", mil(prep["pares_reciprocos_fundidos"]))

# ====================== Parte I ======================
e = J("03_estrutura.json")
cg = e["cg"]
M("numVertices", mil(e["n_vertices"]))
M("numArestas", mil(e["n_arestas"]))
M("grauMin", str(e["grau_min"]))
M("grauMax", mil(e["grau_max"]))
M("grauMedio", dec(e["grau_medio"], 3))
M("grauMediano", dec(e["grau_mediano"], 0))
M("grauDesvio", dec(e["grau_desvio"], 1))
M("densidade", f"\\ensuremath{{{sci(e['densidade'])}}}")
M("numComponentes", mil(e["n_componentes"]))
M("compDois", mil(e["componentes_tamanho_2"]))
M("compMenoresDez", mil(e["componentes_menores_que_10"]))
M("segundaComp", mil(e["maiores_componentes"][1]))
M("terceiraComp", mil(e["maiores_componentes"][2]))
M("clusterMedia", dec(e["clusterizacao_media"], 4))
M("transitividade", dec(e["transitividade_global"], 5))
M("numTriangulos", mil(e["n_triangulos"]))
M("numTriplas", f"\\ensuremath{{{sci(e['n_triplas_conectadas'])}}}")
M("tamCG", mil(cg["n_vertices"]))
M("arestasCG", mil(cg["n_arestas"]))
M("fracCG", dec(100 * cg["fracao_vertices"], 2))
M("fracArestasCG", dec(100 * cg["fracao_arestas"], 2))
M("grauMedioCG", dec(cg["grau_medio"], 3))
M("densidadeCG", f"\\ensuremath{{{sci(cg['densidade'])}}}")
M("clusterMediaCG", dec(cg["clusterizacao_media"], 4))
M("transitividadeCG", dec(cg["transitividade_global"], 5))
d = cg["distancias"]
M("diametro", str(d["diametro"]))
M("raio", str(d["raio"]))
M("bfsDiametro", str(d["bfs_executadas"]))
M("tempoDiametro", dec(d["tempo_s"], 1))
cm = cg["caminho_medio"]
M("caminhoMedio", dec(cm["estimativa"], 3))
M("caminhoMedioICa", dec(cm["ic95"][0], 3))
M("caminhoMedioICb", dec(cm["ic95"][1], 3))
M("caminhoMedioDP", dec(cm["desvio_padrao"], 3))
M("caminhoMedioFontes", mil(cm["n_fontes"]))
import csv as _c
_rows = list(_c.DictReader(open(os.path.join(RES, "distribuicao_graus.csv"), encoding="utf-8")))
_n = sum(int(r["frequencia"]) for r in _rows)
M("fracGrauUm", dec(100 * int(_rows[0]["frequencia"]) / _n, 1))
M("fracGrauAteDois", dec(100 * sum(int(r["frequencia"]) for r in _rows
                                   if int(r["grau"]) <= 2) / _n, 1))
M("hubPrimeiro", e["hubs"][0]["site"].replace("_", "\\_"))
M("hubPrimeiroGrau", mil(e["hubs"][0]["grau"]))

with open(os.path.join(REL, "tab_hubs.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{rlr}\n\\hline\n"
            "\\# & \\textbf{Site (\\textit{hostname})} & \\textbf{Grau} \\\\ \\hline\n")
    for i, h in enumerate(e["hubs"][:15], 1):
        f.write(f"{i} & \\texttt{{{h['site'].replace('_', chr(92)+'_')}}} & {mil(h['grau'])} \\\\\n")
    f.write("\\hline\n\\end{tabular}\n")

# ====================== Parte II ======================
a = J("04_algoritmos.json")
ordem = [("BFS", "Busca em largura (BFS)", "$O(n+m)$"),
         ("DFS", "Busca em profundidade (DFS)", "$O(n+m)$"),
         ("Eulerianidade", "Verificação de eulerianidade", "$O(n)$"),
         ("Dijkstra", "Dijkstra (fila de prioridade)", "$O((n+m)\\log n)$"),
         ("Bellman_Ford", "Bellman--Ford$^{\\dagger}$", "$O(n\\cdot m)$"),
         ("Floyd_Warshall", "Floyd--Warshall$^{\\ddagger}$", "$O(n^3)$"),
         ("Floyd_Warshall_NumPy", "Floyd--Warshall (NumPy)$^{\\ddagger}$", "$O(n^3)$"),
         ("Tarjan_SCC", "Tarjan (comp.\\ fortemente conexas)", "$O(n+m)$"),
         ("Tarjan_articulacoes", "Tarjan (articulações e pontes)", "$O(n+m)$"),
         ("Prim", "Prim (AGM)", "$O(m\\log n)$"),
         ("Kruskal", "Kruskal (AGM)", "$O(m\\log m)$")]
with open(os.path.join(REL, "tab_tempos.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lccrrrc}\n\\hline\n"
            "\\textbf{Algoritmo} & \\textbf{Complexidade} & $T$ & "
            "\\textbf{Média (s)} & \\textbf{DP (s)} & "
            "\\textbf{IC 95\\% (s)} & \\textbf{Dist.} \\\\ \\hline\n")
    for ch, nome, cx in ordem:
        s = a[ch]
        dist = "Z" if s["distribuicao_ic"].startswith("Normal") else "$t$"
        f.write(f"{nome} & {cx} & {s['n_repeticoes']} & {dec(s['media_s'],4)} & "
                f"{dec(s['desvio_s'],4)} & "
                f"[{dec(s['ic95'][0],4)}; {dec(s['ic95'][1],4)}] & {dist} \\\\\n")
    f.write("\\hline\n\\end{tabular}\n")

exp = a["expoentes_empiricos"]
with open(os.path.join(REL, "tab_expoentes.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{llcc}\n\\hline\n"
            "\\textbf{Algoritmo} & \\textbf{Variável de custo $x$} & "
            "$b$ \\textbf{teórico} & $b$ \\textbf{observado} \\\\ \\hline\n")
    linhas = [("BFS", "$n+m$", "BFS_em_(n+m)"),
              ("DFS", "$n+m$", "DFS_em_(n+m)"),
              ("Dijkstra", "$(n+m)\\log n$", "Dijkstra_em_(n+m)logn"),
              ("Prim", "$m\\log n$", "Prim_em_mlogn"),
              ("Kruskal", "$m\\log m$", "Kruskal_em_mlogm"),
              ("Tarjan", "$n+m$", "Tarjan_em_(n+m)"),
              ("Bellman--Ford", "$n\\cdot m$", "Bellman_Ford_em_nm"),
              ("Floyd--Warshall", "$n^3$", "Floyd_Warshall_em_n3")]
    for nome, var, ch in linhas:
        f.write(f"{nome} & {var} & 1,00 & {dec(exp[ch],3)} \\\\\n")
    f.write("\\hline\n\\end{tabular}\n")

eul = a["Eulerianidade"]["resultado"]
M("expBF", dec(exp["Bellman_Ford_em_nm"], 3))
M("fwMemoria", dec(a["grafo"]["n"] ** 2 * 8 / 1e12, 1) + "~TB "
   "para a matriz de dist\u00e2ncias")
M("grauImpar", mil(eul["vertices_grau_impar"]))
M("fracGrauImpar", dec(100 * eul["vertices_grau_impar"] / e["cg"]["n_vertices"], 1))
scc = a["Tarjan_SCC"]["resultado"]
M("numSCC", mil(scc["n_scc"]))
M("maiorSCC", mil(scc["maior_scc"]))
M("fracSCC", dec(100 * scc["fracao_maior_scc"], 2))
M("sccTriviais", mil(scc["scc_triviais"]))
M("nDirigido", mil(scc["n_vertices_dirigido"]))
M("mDirigido", mil(scc["m_arestas_dirigido"]))
art = a["Tarjan_articulacoes"]["resultado"]
M("numArticulacoes", mil(art["vertices_de_articulacao"]))
M("fracArticulacoes", dec(100 * art["fracao_articulacoes"], 2))
M("numPontes", mil(art["pontes"]))
M("custoMST", dec(a["Prim"]["resultado"]["custo_total"], 1))
M("arestasMST", mil(a["Prim"]["resultado"]["arestas_na_arvore"]))
M("custoMSTkruskal", dec(a["Kruskal"]["resultado"]["custo_total"], 1))
M("subBF", mil(a["Bellman_Ford"]["resultado"]["n"]))
M("subBFm", mil(a["Bellman_Ford"]["resultado"]["m"]))
M("subFW", mil(a["Floyd_Warshall_NumPy"]["resultado"]["n"]))

# ====================== Parte III.1 ======================
sw = J("05_small_world.json")
M("swC", dec(sw["real"]["C"]["C"], 4))
M("swL", dec(sw["real"]["L"]["L"], 3))
M("swLICa", dec(sw["real"]["L"]["ic95"][0], 3))
M("swLICb", dec(sw["real"]["L"]["ic95"][1], 3))
M("swCer", f"\\ensuremath{{{sci(sw['er']['C']['C'])}}}")
M("swLer", dec(sw["er"]["L"]["L"], 3))
M("swCconf", dec(sw["configuracao"]["C"]["C"], 4))
M("swLconf", dec(sw["configuracao"]["L"]["L"], 3))
M("swRazaoCer", dec(sw["er"]["razao_C"], 1))
M("swRazaoLer", dec(sw["er"]["razao_L"], 3))
M("swSigmaEr", dec(sw["er"]["sigma_small_world"], 1))
M("swRazaoCconf", dec(sw["configuracao"]["razao_C"], 1))
M("swRazaoLconf", dec(sw["configuracao"]["razao_L"], 3))
M("swSigmaConf", dec(sw["configuracao"]["sigma_small_world"], 1))
M("swCanalitico", f"\\ensuremath{{{sci(sw['analitico']['C_rand_esperado_<k>/n'])}}}")
M("swLanalitico", dec(sw["analitico"]["L_rand_esperado_ln(n)/ln(<k>)"], 2))
M("swRazaoCanalitico", mil(sw["real"]["C"]["C"] / sw["analitico"]["C_rand_esperado_<k>/n"]))
M("swCerICa", f"\\ensuremath{{{sci(sw['er']['C']['ic95'][0])}}}")
M("swCerICb", f"\\ensuremath{{{sci(sw['er']['C']['ic95'][1])}}}")
M("swTrocas", mil(sw["configuracao"]["trocas_aceitas"]))
M("swAmostraC", mil(sw["real"]["C"]["n_vertices_avaliados"]))
M("LtresA", dec(cg["caminho_medio"]["estimativa"], 3))
M("LtresAICa", dec(cg["caminho_medio"]["ic95"][0], 3))
M("LtresAICb", dec(cg["caminho_medio"]["ic95"][1], 3))
M("LtresAn", mil(cg["caminho_medio"]["n_fontes"]))
M("LtresB", dec(sw["real"]["L"]["L"], 3))
M("LtresBICa", dec(sw["real"]["L"]["ic95"][0], 3))
M("LtresBICb", dec(sw["real"]["L"]["ic95"][1], 3))
M("LtresBn", mil(sw["real"]["L"]["n_fontes"]))
M("swAmostraL", mil(sw["real"]["L"]["n_fontes"]))

# ====================== Parte III.2 ======================
pl = J("06_lei_potencia.json")
M("plGamma", dec(pl["mle"]["gamma"], 3))
M("plGammaErro", dec(pl["mle"]["erro_padrao_gamma"], 3))
M("plKmin", str(pl["mle"]["k_min"]))
M("plNcauda", mil(pl["mle"]["n_cauda"]))
M("plFracCauda", dec(100 * pl["mle"]["fracao_cauda"], 2))
M("plKS", dec(pl["mle"]["ks_D"], 4))
M("plPvalor", dec(pl["bondade_ajuste"]["p_valor"], 2))
M("plBoot", str(pl["bondade_ajuste"]["n_bootstrap"]))
cm_ = pl["comparacao_modelos"]
M("plLRlognormal", dec(cm_["LR_potencia_vs_lognormal"], 1))
M("plLRexp", dec(cm_["LR_potencia_vs_exponencial"], 1))
M("plGammaOLS", dec(pl["ols_loglog"]["gamma"], 2))
M("plROLS", dec(pl["ols_loglog"]["R2"], 2))

# ====================== Parte III.3 ======================
rb = J("07_robustez.json")
M("robR", mil(rb["r"]))
M("robT", str(rb["T"]))
M("robSintacto", mil(rb["intacto"]["A_maior_componente"]))
M("robCompIntacto", mil(rb["intacto"]["B_n_componentes"]))
M("robLintacto", dec(rb["intacto"]["C_caminho_medio"], 3))
M("robDiamIntacto", str(rb["intacto"]["C_diametro"]))
al = rb["aleatorio"]
M("robSrand", dec(al["A_fracao_do_original"]["media"], 4))
M("robSrandDP", f"\\ensuremath{{{sci(al['A_fracao_do_original']['desvio'])}}}")
M("robSrandICa", dec(al["A_fracao_do_original"]["ic95"][0], 4))
M("robSrandICb", dec(al["A_fracao_do_original"]["ic95"][1], 4))
M("robSrandAbs", mil(al["A_maior_componente"]["media"]))
M("robCompRand", mil(al["B_n_componentes"]["media"]))
M("robLrand", dec(al["C_caminho_medio"]["media"], 3))
M("robDiamRand", dec(al["C_diametro"]["media"], 1))
M("robIsoRand", dec(100 * al["D_fracao_isolados"]["media"], 2))
at = rb["ataque_grau"]
M("robScent", dec(at["A_fracao_do_original"], 4))
M("robScentAbs", mil(at["A_maior_componente"]))
M("robCompCent", mil(at["B_n_componentes"]))
M("robLcent", dec(at["C_caminho_medio"], 3))
M("robDiamCent", str(at["C_diametro"]))
M("robIsoCent", dec(100 * at["D_fracao_isolados"], 2))
ap = rb["ataque_pagerank"]
M("robSpr", dec(ap["A_fracao_do_original"], 4))
M("robCompPr", mil(ap["B_n_componentes"]))
M("robLpr", dec(ap["C_caminho_medio"], 3))
M("robDiamPr", str(ap["C_diametro"]))
M("robIsoPr", dec(100 * ap["D_fracao_isolados"], 2))
M("robSobreposicao", dec(100 * rb["ataque_pagerank_sobreposicao_com_grau"], 1))
M("robArestasRemovidas", dec(100 * rb["ataque_grau_arestas_removidas"] / (2 * rb["m"]), 1))
M("robRazaoS", dec(at["A_fracao_do_original"] / al["A_fracao_do_original"]["media"], 3))

with open(os.path.join(REL, "tab_robustez.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lrrrr}\n\\hline\n"
            "\\textbf{Métrica} & \\textbf{Intacto} & "
            "\\textbf{Aleatório ($\\mu\\pm\\sigma$)} & \\textbf{Ataque: grau} & "
            "\\textbf{Ataque: PageRank} \\\\ \\hline\n")
    lin = [("A --- maior componente", "A_maior_componente", 0),
           ("A --- $S$ (fração de $n$)", "A_fracao_do_original", 4),
           ("B --- n\\textsuperscript{o} de componentes", "B_n_componentes", 0),
           ("C --- caminho médio", "C_caminho_medio", 3),
           ("C --- diâmetro", "C_diametro", 1),
           ("D --- fração de isolados", "D_fracao_isolados", 4)]
    for nome, ch, c in lin:
        i_ = rb["intacto"][ch]; a_ = al[ch]; g_ = at[ch]; p_ = ap[ch]
        fmt = (lambda x: mil(x)) if c == 0 else (lambda x: dec(x, c))
        f.write(f"{nome} & {fmt(i_)} & {fmt(a_['media'])} $\\pm$ "
                f"{dec(a_['desvio'], max(c,4) if c else 1)} & {fmt(g_)} & {fmt(p_)} \\\\\n")
    f.write("\\hline\n\\end{tabular}\n")

# ====================== Parte III (extra) ======================
md = J("08_modelos.json")
with open(os.path.join(REL, "tab_modelos.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lrrrr}\n\\hline\n"
            " & \\textbf{G106 (real)} & \\textbf{Erd\\H{o}s--Rényi} & "
            "\\textbf{Barabási--Albert} & \\textbf{Watts--Strogatz} \\\\ \\hline\n")
    def fmt_mod(x, modo):
        if modo == "int":
            return mil(x)
        if modo == "peq":                       # ordens de grandeza distintas
            return ("\\ensuremath{" + sci(x, 1) + "}" if x != 0 and abs(x) < 1e-3
                    else dec(x, 5))
        return dec(x, modo)
    lin = [("Grau médio $\\langle k\\rangle$", "grau_medio", 3),
           ("Grau máximo", "grau_max", "int"),
           ("Clusterização média $C$", "clusterizacao_media", "peq"),
           ("Transitividade global", "transitividade_global", "peq"),
           ("Caminho médio $L$", "caminho_medio", 3),
           ("Fração na comp.\\ gigante", "fracao_componente_gigante", 4),
           ("Assortatividade de grau", "assortatividade_grau", 4),
           ("$\\gamma$ (cauda 5\\%)", "gamma_cauda_5pct", 2)]
    for nome, ch, modo in lin:
        vals = [md[k][ch] for k in ("real", "er", "ba", "ws")]
        f.write(nome + " & " + " & ".join(fmt_mod(v, modo) for v in vals) + " \\\\\n")
    f.write("\\hline\n\\end{tabular}\n")
M("modRealAssort", dec(md["real"]["assortatividade_grau"], 3))
M("modBaC", f"\\ensuremath{{{sci(md['ba']['clusterizacao_media'], 1)}}}")
M("modBaL", dec(md["ba"]["caminho_medio"], 3))
M("modWsC", dec(md["ws"]["clusterizacao_media"], 3))
M("modWsL", dec(md["ws"]["caminho_medio"], 2))
M("modErC", f"\\ensuremath{{{sci(md['er']['clusterizacao_media'], 1)}}}")
M("modErL", dec(md["er"]["caminho_medio"], 3))




# ====================== triangulos (implementacao propria) ======================
tg = J("13_triangulos.json")
M("triExato", mil(tg["n_triangulos_exato"]))
M("triTempo", dec(tg["tempo_s"], 1))
M("triOps", f"\\ensuremath{{{sci(tg['custo_do_algoritmo_operacoes'], 1)}}}")
M("triConfere", "coincidem at\u00e9 o \u00faltimo d\u00edgito"
  if tg["confere_com_igraph"] else "DIVERGEM")
M("triVertices", mil(tg["vertices_em_algum_triangulo"]))
M("triFracVertices", dec(100 * tg["fracao_vertices_em_triangulo"], 1))

# ====================== betweenness ======================
bt = J("12_betweenness.json")
M("btK", str(bt["K_pivos"]))
M("btTempo", mil(bt["tempo_brandes_s"]))
M("btTempoExato", mil(bt["tempo_estimado_exato_h"]))
M("btSobreposicao", dec(100 * bt["sobreposicao_com_grau"], 1))
M("btSpearman", dec(bt["correlacao_spearman_grau_betweenness"], 3))
M("btValidacao", "confere" if bt["validacao"]["status"] == "OK" else "DIVERGEM")
_ab, _ag = bt["ataque_betweenness"], bt["ataque_grau"]
M("btS", dec(_ab["A_fracao_do_original"], 4))
M("btComp", mil(_ab["B_n_componentes"]))
M("btL", dec(_ab["C_caminho_medio"], 2))
M("btDiam", str(_ab["C_diametro"]))
M("btIso", dec(100 * _ab["D_fracao_isolados"], 2))
M("btGrauS", dec(_ag["A_fracao_do_original"], 4))
M("btGrauL", dec(_ag["C_caminho_medio"], 2))
M("btGrauDiam", str(_ag["C_diametro"]))
M("btGrauIso", dec(100 * _ag["D_fracao_isolados"], 2))
_dg = J("12_diagnostico_linhas.json")
M("btZeros", dec(_dg["zeros_pct"], 1))
M("btFolhas", mil(_dg["linha_folhas_CB_raw_igual_K"]))
M("btPorteiros", mil(_dg["linha_porteiros_n_vertices"]))
M("btPorteirosRazao", dec(_dg["linha_porteiros_CB_raw_sobre_n_mediana"], 2))
M("btPorteirosPct", dec(100 * _dg["linha_porteiros_n_vertices"] / bt["r"], 2))
M("btCorte", f"\\ensuremath{{{sci(_dg['corte_top5pct'], 1)}}}")
M("btFatorCorte", mil(5.0e8 / _dg["corte_top5pct"]))
M("btRazaoS", dec(_ab["A_fracao_do_original"] / _ag["A_fracao_do_original"], 2))

with open(os.path.join(REL, "tab_betweenness.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lrr}\n\\hline\n"
            "\\textbf{Métrica} & \\textbf{Ataque: grau} & "
            "\\textbf{Ataque: intermediação} \\\\ \\hline\n")
    _l = [("A --- maior componente", "A_maior_componente", "int"),
          ("A --- $S$ (fração de $n$)", "A_fracao_do_original", 4),
          ("B --- n\\textsuperscript{o} de componentes", "B_n_componentes", "int"),
          ("C --- caminho médio", "C_caminho_medio", 2),
          ("C --- diâmetro", "C_diametro", "int"),
          ("D --- fração de isolados", "D_fracao_isolados", 4)]
    for nome, ch, modo in _l:
        g_, b_ = _ag[ch], _ab[ch]
        fm = (lambda x: mil(x)) if modo == "int" else (lambda x: dec(x, modo))
        f.write(f"{nome} & {fm(g_)} & {fm(b_)} \\\\\n")
    f.write("\\hline\n\\end{tabular}\n")

with open(os.path.join(REL, "tab_topbt.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{rlr}\n\\hline\n"
            "\\# & \\textbf{Site} & \\textbf{Grau} \\\\ \\hline\n")
    for i, h in enumerate(bt["top_betweenness"][:10], 1):
        f.write(f"{i} & \\texttt{{{h['site'].replace('_', chr(92)+'_')}}} & "
                f"{mil(h['grau'])} \\\\\n")
    f.write("\\hline\n\\end{tabular}\n")

# ====================== validacao do recorte temporal ======================
vr = J("11_validacao_recorte.json")
cp = vr["comparacao"]
M("setN", mil(vr["setembro"]["n_vertices"]))
M("setM", mil(vr["setembro"]["n_arestas"]))
M("setCG", dec(100 * vr["setembro"]["fracao_componente_gigante"], 2))
M("setDiam", str(vr["setembro"]["diametro"]))
M("setRaio", str(vr["setembro"]["raio"]))
M("setL", dec(vr["setembro"]["caminho_medio"]["L"], 3))
M("setGamma", dec(vr["setembro"]["cauda"]["gamma"], 3))
M("setKmin", str(vr["setembro"]["cauda"]["k_min"]))
M("setSrand", dec(vr["setembro"]["robustez"]["S_aleatorio_media"], 4))
M("setScent", dec(vr["setembro"]["robustez"]["S_ataque_grau"], 4))
M("setVarN", dec(cp["n_vertices"]["variacao_pct"], 1))
M("setVarM", dec(cp["n_arestas"]["variacao_pct"], 1))
M("setVarCG", dec(abs(cp["fracao_componente_gigante"]["variacao_pct"]), 2))
M("setVarL", dec(cp["caminho_medio"]["variacao_pct"], 1))
M("setVarGamma", dec(abs(cp["gamma"]["variacao_pct"]), 1))
M("setVarSrand", dec(cp["S_aleatorio"]["variacao_pct"], 1))
M("setVarScent", mil(cp["S_ataque_grau"]["variacao_pct"]))
M("agoSrandDez", dec(vr["agosto"]["robustez"]["S_aleatorio_media"], 4))
M("agoLdez", dec(vr["agosto"]["caminho_medio"]["L"], 3))
_cc = vr["concentracao_de_arestas"]
M("concAgo", dec(100 * _cc["agosto"]["fracao_pontas_no_top5"], 1))
M("concSet", dec(100 * _cc["setembro"]["fracao_pontas_no_top5"], 1))
M("grauMedioCGset", dec(_cc["setembro"]["grau_medio_cg"], 3))

_lin = [("Vértices $n$", "n_vertices", "int"),
        ("Arestas $m$", "n_arestas", "int"),
        ("Grau médio $\\langle k\\rangle$", "grau_medio", 3),
        ("Grau máximo", "grau_max", "int"),
        ("Densidade", "densidade", "peq"),
        ("Componentes conexas", "n_componentes", "int"),
        ("Fração na comp.\\ gigante", "fracao_componente_gigante", 4),
        ("Clusterização média $C$", "clusterizacao_media", 4),
        ("Diâmetro (exato)", "diametro", "int"),
        ("Raio (exato)", "raio", "int"),
        ("Caminho médio $L$", "caminho_medio", 3),
        ("Expoente $\\gamma$ da cauda", "gamma", 3),
        ("$S$ após remoção aleatória", "S_aleatorio", 4),
        ("$S$ após ataque por grau", "S_ataque_grau", 4)]
with open(os.path.join(REL, "tab_recorte.tex"), "w", encoding="utf-8") as f:
    f.write("\\begin{tabular}{lrrr}\n\\hline\n"
            "\\textbf{Grandeza} & \\textbf{ago/2008} & \\textbf{set/2008} & "
            "\\textbf{variação} \\\\ \\hline\n")
    for nome, ch, modo in _lin:
        a_, s_ = cp[ch]["agosto"], cp[ch]["setembro"]
        v_ = cp[ch]["variacao_pct"]
        if modo == "int":
            fa, fs = mil(a_), mil(s_)
        elif modo == "peq":
            fa = "\\ensuremath{" + sci(a_, 1) + "}"
            fs = "\\ensuremath{" + sci(s_, 1) + "}"
        else:
            fa, fs = dec(a_, modo), dec(s_, modo)
        f.write(f"{nome} & {fa} & {fs} & {dec(v_, 1)}\\% \\\\\n")
    f.write("\\hline\n\\end{tabular}\n")

# ====================== validacao ======================
vl = J("09_validacao.json")
M("valGrafos", str(vl["n_grafos_testados"]))
M("valChecagens", str(vl["n_verificacoes"]))
M("valStatus", "todas as verificações passaram" if vl["status"] == "OK"
   else "ATENÇÃO: houve falhas")

with open(os.path.join(REL, "numeros.tex"), "w", encoding="utf-8") as f:
    f.write("% gerado automaticamente por scripts/10_gerar_tex.py\n")
    f.write("\n".join(macros) + "\n")
print(f"{len(macros)} macros gravadas em relatorio/numeros.tex")

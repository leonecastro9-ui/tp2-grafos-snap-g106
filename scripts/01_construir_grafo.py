#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TP2 - Analisando Grafos do SNAP
Grafo G106 = memetracker9 (96 million Memetracker memes)

Etapa 1: construcao do grafo a partir do arquivo bruto do SNAP.

O arquivo bruto e um grafo TEMPORAL de hiperlinks entre DOCUMENTOS
(posts de blog / artigos de midia).  Formato (registros separados por
linha em branco):

    P <url do documento>
    T <timestamp>
    Q <frase/meme>        (0..n)
    L <url apontada>      (0..n)

Tratamento aplicado (ver relatorio, secao "Tratamento dos dados"):
  1. Recorte temporal   : apenas um intervalo [DATA_INI, DATA_FIM).
  2. Agregacao por site : cada URL e reduzida ao seu HOSTNAME normalizado
                          (minusculas, sem "www.", sem porta).  O grafo
                          passa de "documento -> documento" para
                          "site -> site", que e a granularidade usada na
                          literatura para redes de hiperlinks.
  3. Grafo temporal -> estatico: uniao de todas as arestas do intervalo,
                          com peso = numero de hiperlinks observados.
  4. Remocao de auto-loops (links internos do proprio site).

Saida: CSV "origem,destino,peso" (arestas DIRIGIDAS agregadas).
"""

import gzip
import sys
import csv
from collections import defaultdict

ENTRADA = sys.argv[1] if len(sys.argv) > 1 else "dados_brutos/quotes_2008-08.txt.gz"
SAIDA = sys.argv[2] if len(sys.argv) > 2 else "dados/arestas_dirigidas.csv"
DATA_INI = sys.argv[3] if len(sys.argv) > 3 else "2008-08-01"
DATA_FIM = sys.argv[4] if len(sys.argv) > 4 else "2008-09-01"   # exclusivo


def hostname(url):
    """Extrai o hostname normalizado de uma URL. Retorna None se invalida."""
    i = url.find("://")
    if i >= 0:
        url = url[i + 3:]
    # corta no primeiro /, ? ou #
    for ch in ("/", "?", "#"):
        j = url.find(ch)
        if j >= 0:
            url = url[:j]
    # remove credenciais e porta
    k = url.rfind("@")
    if k >= 0:
        url = url[k + 1:]
    k = url.find(":")
    if k >= 0:
        url = url[:k]
    url = url.strip().lower()
    if url.startswith("www."):
        url = url[4:]
    if not url or "." not in url or " " in url:
        return None
    return url


def main():
    arestas = defaultdict(int)
    n_doc = 0          # documentos no recorte
    n_doc_total = 0    # documentos lidos
    n_link = 0         # hiperlinks no recorte
    n_loop = 0         # auto-loops descartados
    n_invalid = 0

    src = None
    dentro = False

    with gzip.open(ENTRADA, "rt", encoding="utf-8", errors="replace") as f:
        for linha in f:
            if len(linha) < 2:
                continue
            tag = linha[0]
            if tag == "P":
                n_doc_total += 1
                src = hostname(linha[2:].rstrip("\n"))
                dentro = False
            elif tag == "T":
                data = linha[2:12]
                dentro = (DATA_INI <= data < DATA_FIM) and (src is not None)
                if dentro:
                    n_doc += 1
            elif tag == "L" and dentro:
                dst = hostname(linha[2:].rstrip("\n"))
                if dst is None:
                    n_invalid += 1
                    continue
                n_link += 1
                if dst == src:
                    n_loop += 1
                    continue
                arestas[(src, dst)] += 1

            if n_doc_total % 2000000 == 0 and tag == "P":
                print(f"  ... {n_doc_total:,} documentos lidos, "
                      f"{len(arestas):,} arestas dirigidas distintas",
                      file=sys.stderr, flush=True)

    with open(SAIDA, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["origem", "destino", "peso"])
        for (a, b), p in arestas.items():
            w.writerow([a, b, p])

    sites = set()
    for a, b in arestas:
        sites.add(a); sites.add(b)

    print("=" * 60)
    print(f"Arquivo .............................. {ENTRADA}")
    print(f"Recorte temporal ..................... [{DATA_INI}, {DATA_FIM})")
    print(f"Documentos lidos (arquivo inteiro) ... {n_doc_total:,}")
    print(f"Documentos no recorte ................ {n_doc:,}")
    print(f"Hiperlinks no recorte ................ {n_link:,}")
    print(f"  auto-loops descartados ............. {n_loop:,}")
    print(f"  URLs invalidas descartadas ......... {n_invalid:,}")
    print(f"Arestas dirigidas distintas .......... {len(arestas):,}")
    print(f"Sites (hostnames) distintos .......... {len(sites):,}")
    print("=" * 60)


if __name__ == "__main__":
    main()

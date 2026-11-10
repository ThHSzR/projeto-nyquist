"""Análise no terminal de um caso pronto ou de uma malha informada."""
import argparse
import json

import numpy as np

from modelos import MODELOS
from nucleo import Malha, margens, nyquist, pade, routh


def main():
    p = argparse.ArgumentParser(description=__doc__)
    origem = p.add_mutually_exclusive_group(required=True)
    origem.add_argument("--caso", choices=MODELOS)
    origem.add_argument("--num", type=float, nargs="+", help="Coeficientes em ordem decrescente.")
    p.add_argument("--den", type=float, nargs="+")
    p.add_argument("--atraso", type=float, default=0)
    p.add_argument("--ordem-pade", type=int, default=2)
    args = p.parse_args()
    if args.num is not None and not args.den:
        p.error("--num requer --den")
    if args.caso and (args.den or args.atraso):
        p.error("Use --num e --den para informar uma planta personalizada.")
    m = MODELOS[args.caso] if args.caso else Malha(tuple(args.num), tuple(args.den), args.atraso)
    racional = pade(m, args.ordem_pade) if m.atraso else m
    mg, mf = margens(m)
    pchar = racional.caracteristico()
    tabela, n, eventos = routh(pchar)
    print(json.dumps({
        "Nyquist": nyquist(m),
        "margem_ganho_primeira": mg[0] if mg else None,
        "margem_fase_primeira": mf[0] if mf else None,
        "faixa_margens_rad_s": [1e-4, 100],
        "unidades_margens": ["rad/s", "dB ou graus"],
        "aproximacao_racional_Pade": args.ordem_pade if m.atraso else None,
        "polinomio_caracteristico_racional": pchar.tolist(),
        "Routh": {"RHP": n, "eventos": eventos, "tabela_normalizada": tabela.tolist()},
        "polos_racionais": [[float(z.real), float(z.imag)] for z in np.roots(pchar)],
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

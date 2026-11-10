"""Gera um gráfico e resume o envelope e o RTT dos ensaios da bancada."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("csv", type=Path, nargs="+")
    p.add_argument("--saida", type=Path, default=Path("comparacao_hil.png"))
    args = p.parse_args()
    fig, ax = plt.subplots(figsize=(9, 5), layout="constrained")
    for caminho in args.csv:
        dados = np.genfromtxt(caminho, delimiter=",", names=True, ndmin=1)
        t, y, rtt = dados["t_s"], dados["y"], dados["rtt_ms"]
        ax.plot(t, y, label=caminho.stem)
        anterior = (t >= 30) & (t < 40)
        recente = (t >= 50) & (t <= 60)
        msg = f"{caminho.name} | n={len(t)} | RTT p50={np.median(rtt):.3f} ms"
        msg += f" | RTT p95={np.percentile(rtt, 95):.3f} ms"
        if anterior.any() and recente.any():
            a = max(abs(y[anterior] - 1))
            b = max(abs(y[recente] - 1))
            if a > 0:
                msg += f" | razão dos envelopes={b/a:.3f}"
        print(msg)
    ax.axhline(1, color="#64748b", ls="--", lw=.8)
    ax.set(xlabel="Tempo (s)", ylabel="Saída y(t)")
    ax.grid(alpha=.25)
    ax.legend()
    fig.savefig(args.saida, dpi=160)


if __name__ == "__main__":
    main()

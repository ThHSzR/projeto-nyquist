"""Reproduz a validação, o benchmark, a bancada e todos os gráficos."""
import argparse
import csv
from math import comb, pi, sqrt
from pathlib import Path
import platform
import statistics
import subprocess
import tempfile
from time import process_time_ns

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import linalg, optimize, signal

from bancada import canal, simular, TS
from modelos import D, MODELOS, TAU_CRITICO
from nucleo import Malha, lugar_raizes, margens, nyquist, pade, routh

PASTA = Path(__file__).resolve().parent
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.22})
AZUL, VERMELHO, VERDE = "#174f87", "#a92d3d", "#28744e"


def csv_gravar(pasta, nome, chaves, linhas):
    with (pasta / nome).open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(chaves)
        w.writerows(linhas)


def salvar(pasta, nome, figura):
    figura.savefig(pasta / f"{nome}.png", dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(figura)


def planta_polos(malha):
    return linalg.eigvals(linalg.companion(malha.caracteristico()))


def principal(pasta):
    fig, axs = plt.subplots(2, 2, figsize=(9, 6.6), layout="constrained")
    bode, ab = plt.subplots(2, 4, figsize=(14, 5.8), layout="constrained")
    locus, al = plt.subplots(2, 2, figsize=(9, 6.6), layout="constrained")
    degrau, ad = plt.subplots(2, 2, figsize=(9, 6.6), layout="constrained")
    freq = np.geomspace(1e-3, 100, 3600)
    linhas = []
    for idx, (nome, m) in enumerate(MODELOS.items()):
        diag = nyquist(m)
        mg, mf = margens(m)
        racional = pade(m, 2) if m.atraso else m
        w, ref = signal.freqresp((m.num, m.den), w=freq)
        erro_freq = np.max(np.abs(m.avaliar(1j*w) - ref*np.exp(-1j*w*m.atraso)))
        polos = planta_polos(racional)
        tab, trocas, eventos = routh(racional.caracteristico())
        assert trocas == np.count_nonzero(polos.real > 1e-7) == diag["Z"]
        assert erro_freq < 1e-10
        t, y = signal.step((racional.num, racional.caracteristico()), T=np.linspace(0, 25, 1501))
        linhas.append((nome, diag["P"], diag["N"], diag["Z"], diag["estado"],
                       mg[0][0], mg[0][1], mf[0][0], mf[0][1], trocas, max(polos.real),
                       erro_freq, diag["amostras"], diag["distancia"]))
        csv_gravar(pasta, f"routh_{nome}.csv", ("potencia", *[f"c{i+1}" for i in range(tab.shape[1])]),
                   ((len(tab)-1-i, *row) for i, row in enumerate(tab)))
        csv_gravar(pasta, f"polos_{nome}.csv", ("real", "imag"), ((p.real,p.imag) for p in polos))
        csv_gravar(pasta, f"degrau_{nome}.csv", ("t_s", "y"), zip(t, y))

        ax = axs.flat[idx]
        z = m.avaliar(1j * freq)
        for vetor, cor, label in [(z, AZUL, "ω > 0"), (z.conj(), "#6396b8", "ω < 0")]:
            ax.plot(vetor.real, vetor.imag, color=cor, lw=1.5, label=label)
        ax.plot(-1, 0, "x", color=VERMELHO, ms=9, mew=2, label="−1 + j0")
        ax.axhline(0, color="#64748b", lw=.8)
        ax.axvline(0, color="#64748b", lw=.8)
        ax.set(title=f"Sistema {nome} | P={diag['P']}, N={diag['N']}, Z={diag['Z']}",
               xlabel="Re L(jω)", ylabel="Im L(jω)")
        if m.atraso:
            ax.set(xlim=(-1.65, .35), ylim=(-1.5, 1.5))
        ax.legend(fontsize=8, loc="best")

        ab[0, idx].semilogx(freq, 20*np.log10(abs(z)), color=AZUL)
        ab[1, idx].semilogx(freq, np.rad2deg(np.unwrap(np.angle(z))), color=AZUL)
        ab[0, idx].set(title=f"Sistema {nome}", ylabel="Magnitude (dB)")
        ab[1, idx].set(xlabel="ω (rad/s)", ylabel="Fase (°)")
        ab[0, idx].axhline(0, color=VERMELHO, lw=.8, ls="--")
        ab[0, idx].axvline(mf[0][0], color=VERDE, lw=.8, ls=":")
        faixa_fase = (freq <= 7) if m.atraso else (freq <= 100)
        ab[1, idx].lines[0].set_data(freq[faixa_fase], np.rad2deg(np.unwrap(np.angle(z)))[faixa_fase])
        ab[1, idx].axhline(-180, color=VERMELHO, lw=.8, ls="--")
        ab[1, idx].axvline(mf[0][0], color=VERDE, lw=.8, ls=":")
        ab[1, idx].relim()
        ab[1, idx].autoscale_view()

        ax = al.flat[idx]
        final = {"A": 7, "B": 2.5, "C": 2, "D": 2.5}[nome]
        r = lugar_raizes(racional, np.linspace(0, final, 501))
        for ramo in r.T:
            ax.plot(ramo.real, ramo.imag, color=AZUL, lw=1.3)
        ax.plot(r[0].real, r[0].imag, "x", color="#151b2e", label="Polos abertos")
        zeros = np.roots(racional.num)
        if len(zeros):
            ax.plot(zeros.real, zeros.imag, "o", mfc="none", mec="#151b2e", label="Zeros abertos")
        ax.plot(polos.real, polos.imag, "s", color=VERMELHO, ms=4, label="Polos em K=1")
        ax.axvline(0, color="#64748b", lw=.8)
        ax.set(title=f"Sistema {nome}" + (" | Padé [2/2]" if m.atraso else ""),
               xlabel="Re s (1/s)", ylabel="Im s (1/s)")
        ax.legend(fontsize=7, loc="best")

        ax = ad.flat[idx]
        ax.plot(t, y, color=AZUL)
        ax.axhline(0, color="#64748b", lw=.8)
        ax.set(title=f"Sistema {nome}" + (" | Padé [2/2]" if m.atraso else ""),
               xlabel="Tempo (s)", ylabel="Saída y(t)")
    for nome, figura in [("nyquist", fig), ("bode", bode), ("lugar_raizes", locus), ("degrau", degrau)]:
        salvar(pasta, nome, figura)
    csv_gravar(pasta, "validacao.csv",
               ("sistema", "P", "N", "Z", "estado", "w_fase", "MG_dB", "w_ganho", "MF_graus",
                "Routh_RHP", "max_Re_polo_SciPy", "erro_freq_SciPy", "amostras_Nyquist", "dist_min"),
               linhas)


def atraso(pasta):
    def caracteristica(s, tau):
        return s*(s+4) + sqrt(17)*np.exp(-tau*s)
    limites = []
    for tau in (1.0, TAU_CRITICO, 1.35, 1.5):
        m = Malha(D.num, D.den, tau)
        def f(x):
            v = caracteristica(complex(*x), tau)
            return v.real, v.imag
        sol = optimize.root(f, [0, 1], tol=1e-11)
        residual = abs(caracteristica(complex(*sol.x), tau))
        assert residual < 1e-9
        d = nyquist(m)
        limites.append((tau, d["estado"], d["P"], d["N"], d["Z"], *sol.x, residual))
    csv_gravar(pasta, "atraso_exato.csv", ("tau_s", "estado", "P", "N", "Z",
                                         "Re_raiz", "Im_raiz", "residuo"), limites)
    linhas = []
    for n in range(1, 5):
        def modelo(tau):
            return pade(Malha(D.num, D.den, tau), n)
        limite = optimize.brentq(lambda tau: margens(modelo(tau))[1][0][1], 1, 1.8)
        m = modelo(D.atraso)
        polos = planta_polos(m)
        z = routh(m.caracteristico())[1]
        assert z == np.count_nonzero(polos.real > 1e-7)
        linhas.append((n, limite, 100*(limite-TAU_CRITICO)/TAU_CRITICO, z, max(polos.real),
                       margens(m)[1][0][1]))
    csv_gravar(pasta, "pade.csv", ("ordem", "tau_crit_s", "erro_rel_pct", "RHP_Routh",
                                  "max_Re_polo_SciPy", "MF_graus_em_1_35"), linhas)
    fig, ax = plt.subplots(figsize=(7, 4.2), layout="constrained")
    w = np.geomspace(.55, 1.9, 2000)
    for nome, m, cor in [("Exato, τ=1,00 s", Malha(D.num,D.den,1), VERDE),
                          ("Exato, τ=1,35 s", D, VERMELHO),
                          ("Padé [1/1], τ=1,35 s", pade(D,1), AZUL),
                          ("Padé [2/2], τ=1,35 s", pade(D,2), "#7353a4")]:
        z = m.avaliar(1j*w)
        ax.plot(z.real, z.imag, lw=2, label=nome, color=cor)
    ax.plot(-1, 0, "kx", ms=10, mew=2)
    ax.set(xlim=(-1.4, -.6), ylim=(-.31, .30), xlabel="Re L(jω)", ylabel="Im L(jω)")
    ax.legend(loc="best", fontsize=8)
    salvar(pasta, "atraso_pade", fig)


def benchmark(pasta):
    def medir(f):
        for _ in range(4):
            f()
        dados = []
        for _ in range(25):
            inicio = process_time_ns()
            f()
            dados.append((process_time_ns()-inicio)/1000)
        return statistics.median(dados), float(np.percentile(dados, 90))
    linhas = []
    for n in (5, 10, 20):
        alvo = np.array([comb(n, i) for i in range(n+1)], dtype=float)
        den = alvo.copy()
        den[-1] -= .25
        m = Malha((.25,), tuple(den))
        grade = np.geomspace(1e-3, 100, 4096)
        ganhos = np.linspace(0, 2, 301)
        polos = np.roots(alvo)
        r32, z32, _ = routh(alvo, dtype=np.float32)
        r64, z64, _ = routh(alvo)
        dados = {"Routh": lambda: routh(alvo), "Nyquist": lambda: nyquist(m),
                 "Bode+margens": lambda: (m.avaliar(1j*grade), margens(m)),
                 "Lugar das raízes": lambda: lugar_raizes(m, ganhos)}
        erro_polo = max(abs(polos+1))
        residuo = max(abs(np.polyval(alvo, polos)) / np.polyval(abs(alvo), abs(polos)))
        assert z32 == z64 == nyquist(m)["Z"] == 0
        for nome, f in dados.items():
            mediana, p90 = medir(f)
            linhas.append((n, nome, mediana, p90, erro_polo, residuo, np.max(abs(r32-r64))))
    csv_gravar(pasta, "benchmark.csv", ("ordem", "metodo", "mediana_CPU_us", "p90_CPU_us",
                                       "erro_max_polo", "residuo_rel_max", "delta_Routh_32_64"), linhas)
    fig, ax = plt.subplots(figsize=(7, 3.7), layout="constrained")
    for nome in dados:
        pares = [(n, med) for n, metodo, med, *_ in linhas if metodo == nome]
        ax.plot(*zip(*pares), marker="o", label=nome)
    ax.set(yscale="log", xticks=[5,10,20], xlabel="Ordem n", ylabel="Mediana de CPU (µs)")
    ax.legend(fontsize=8)
    salvar(pasta, "benchmark", fig)


def bancada(pasta):
    def ldc(w):
        num, den, _ = discretizacao
        z = np.exp(1j*w*TS)
        inverso = 1 - 1/z
        c = sqrt(17) * (7/16 + TS/(4*inverso) + (9/64)*4*inverso/(inverso+4*TS))
        return c*np.polyval(num[0], z)/np.polyval(den,z)

    discretizacao = signal.cont2discrete(([1.], [1., 2., 1.]), TS, method="zoh")
    wg = optimize.brentq(lambda w: abs(ldc(w))-1, .7, 1.3)
    tau_digital = (pi + np.angle(ldc(wg))) / wg
    linhas = []
    fig, ax = plt.subplots(figsize=(8, 4), layout="constrained")
    with tempfile.TemporaryDirectory() as temp:
        executavel = Path(temp) / "pid_host"
        subprocess.run(["g++", "-std=c++11", "-O2", "-Wall", "-Wextra", "-Werror", "-pedantic",
                        str(PASTA/"hil/pid_host.cpp"), "-o", str(executavel)], check=True)
        with canal(None, executavel) as trocar:
            assert trocar("entrada inválida") == "ERRO"
            for tau, cor in [(1.0, VERDE), (1.32, AZUL), (1.5, VERMELHO)]:
                dados = list(simular(trocar, tau, 60))
                t, y, *_ = np.array(dados).T
                anterior = max(abs(y[(t>=30)&(t<40)]-1))
                recente = max(abs(y[t>=50]-1))
                linhas.append((tau, wg, tau_digital, anterior, recente, recente/anterior))
                csv_gravar(pasta, f"bancada_{tau:.2f}.csv",
                           ("t_s","y","u","u_atrasado","rtt_ms"), dados)
                ax.plot(t, y, color=cor, lw=1.4, label=f"τ={tau:.2f} s")
    ax.axhline(1, color="#64748b", ls="--", lw=.8)
    ax.set(xlabel="Tempo simulado (s)", ylabel="Saída y(t)")
    ax.legend(fontsize=8, loc="upper left")
    salvar(pasta, "bancada", fig)
    csv_gravar(pasta, "bancada_resumo.csv", ("tau_s", "wg_digital", "tau_crit_digital",
                                           "max_erro_30a40", "max_erro_50a60", "razao"), linhas)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--saida", type=Path, default=PASTA / "resultados")
    args = p.parse_args()
    pasta = args.saida
    pasta.mkdir(parents=True, exist_ok=True)
    principal(pasta)
    atraso(pasta)
    benchmark(pasta)
    bancada(pasta)
    cpuinfo = Path("/proc/cpuinfo")
    cpu = next((linha.split(":",1)[1].strip() for linha in cpuinfo.read_text().splitlines()
                if linha.startswith("model name")), "") if cpuinfo.is_file() else platform.processor()
    (pasta/"ambiente.txt").write_text(f"{platform.platform()}\nCPU {cpu}\nPython {platform.python_version()}\n"
                                     f"NumPy {np.__version__}\nSciPy {scipy.__version__}\n"
                                     f"Matplotlib {matplotlib.__version__}\n", encoding="utf-8")
    print(f"Resultados escritos em {pasta}")


if __name__ == "__main__":
    main()

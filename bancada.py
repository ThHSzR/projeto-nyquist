"""Bancada com planta ZOH e atraso FIFO. --porta usa o ESP32 físico."""
import argparse
from collections import deque
from contextlib import contextmanager
import csv
from math import exp, isfinite
from pathlib import Path
import subprocess
import time

TS = 0.02


@contextmanager
def canal(porta, executavel):
    if porta:
        import serial
        with serial.Serial(porta, 115200, timeout=0.25, write_timeout=0.25) as conexao:
            time.sleep(2)
            conexao.reset_input_buffer()
            def trocar(linha):
                conexao.write((linha + "\n").encode("ascii"))
                resposta = conexao.readline()
                if not resposta.endswith(b"\n"):
                    raise TimeoutError("Resposta serial ausente ou incompleta.")
                return resposta.decode("ascii").strip()
            yield trocar
    else:
        with subprocess.Popen([str(Path(executavel).resolve())], stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE, text=True, bufsize=1) as processo:
            def trocar(linha):
                processo.stdin.write(linha + "\n")
                processo.stdin.flush()
                return processo.stdout.readline().strip()
            yield trocar
            processo.stdin.close()


def simular(trocar, atraso, duracao, tempo_real=False):
    if not isfinite(atraso) or not isfinite(duracao) or atraso < 0 or duracao < TS:
        raise ValueError("Informe atraso não negativo e duração de pelo menos uma amostra.")
    if trocar("RESET") != "OK":
        raise RuntimeError("O controlador não confirmou o reset.")
    d = round(atraso / TS)
    fila = deque([0.0] * d)
    a = exp(-TS)
    x1 = x2 = 0.0
    inicio = time.monotonic()
    for k in range(round(duracao / TS) + 1):
        instante = k * TS
        if tempo_real:
            time.sleep(max(0, inicio + instante - time.monotonic()))
        envio = time.monotonic()
        resposta = trocar(f"{k} {1-x2:.15g}").split()
        rtt = time.monotonic() - envio
        if len(resposta) != 2 or resposta[0] != str(k):
            raise RuntimeError(f"Resposta inesperada na amostra {k}: {resposta}")
        u = float(resposta[1])
        if not isfinite(u):
            raise ArithmeticError("Controlador enviou um valor não finito.")
        if tempo_real and time.monotonic() > inicio + instante + TS:
            raise TimeoutError(f"Amostra {k} excedeu o prazo de {TS*1000:.0f} ms.")
        fila.append(u)
        aplicado = fila.popleft()
        yield instante, x2, u, aplicado, rtt * 1000
        # ZOH exato de x1'=-x1+u e x2'=x1-x2, y=x2.
        x1, x2 = a*x1 + (1-a)*aplicado, a*(x2+TS*x1) + (1-a*(1+TS))*aplicado


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--porta", help="Por exemplo, COM4 ou /dev/ttyUSB0.")
    parser.add_argument("--executavel", default="pid_host")
    parser.add_argument("--atraso", type=float, default=1.0)
    parser.add_argument("--duracao", type=float, default=60.0)
    parser.add_argument("--saida", type=Path, default=Path("bancada.csv"))
    args = parser.parse_args()
    with canal(args.porta, args.executavel) as trocar:
        dados = list(simular(trocar, args.atraso, args.duracao, bool(args.porta)))
    with args.saida.open("w", newline="") as arquivo:
        escritor = csv.writer(arquivo)
        escritor.writerow(("t_s", "y", "u", "u_atrasado", "rtt_ms"))
        escritor.writerows(dados)
    print(f"{len(dados)} amostras | atraso efetivo {round(args.atraso/TS)*TS:.2f} s | {args.saida}")


if __name__ == "__main__":
    main()

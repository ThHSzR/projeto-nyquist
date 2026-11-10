"""Métodos numéricos para malhas SISO contínuas, reais e estritamente próprias."""
from dataclasses import dataclass
from math import factorial

import numpy as np


@dataclass(frozen=True)
class Malha:
    num: tuple
    den: tuple
    atraso: float = 0.0

    def __post_init__(self):
        if (len(self.den) <= len(self.num) or not self.num or self.den[0] == 0
                or self.num[0] == 0 or self.atraso < 0
                or not np.all(np.isfinite((*self.num, *self.den, self.atraso)))):
            raise ValueError("Use coeficientes reais finitos, malha estritamente própria e atraso não negativo.")

    def avaliar(self, s):
        return np.polyval(self.num, s) / np.polyval(self.den, s) * np.exp(-self.atraso * s)

    def caracteristico(self, ganho=1.0):
        if self.atraso:
            raise ValueError("A malha com atraso puro não possui polinômio característico.")
        return np.polyadd(self.den, ganho * np.asarray(self.num))


def routh(coef, eps=1e-9, tol=1e-13, dtype=np.float64):
    """Retorna tabela com linhas escaladas positivamente, trocas e eventos."""
    coef = np.asarray(coef, dtype=dtype)
    if coef.ndim != 1 or not len(coef) or not np.all(np.isfinite(coef)) or coef[0] == 0:
        raise ValueError("Polinômio inválido.")
    zeros = 0
    while len(coef) > 1 and coef[-1] == 0:
        coef = coef[:-1]
        zeros += 1
    grau = len(coef) - 1
    tabela = np.zeros((grau + 1, (grau + 2) // 2), dtype=dtype)
    tabela[0, :len(coef[::2])] = coef[::2]
    if grau:
        tabela[1, :len(coef[1::2])] = coef[1::2]
    eventos = [f"{zeros} raiz(es) na origem"] if zeros else []
    for i in range(grau + 1):
        if i >= 2:
            anterior = tabela[i - 1]
            acima = tabela[i - 2]
            tabela[i, :-1] = acima[1:] - acima[0] * anterior[1:] / anterior[0]
        if i and np.max(np.abs(tabela[i])) <= tol:
            potencias = grau - i + 1 - 2 * np.arange(tabela.shape[1])
            tabela[i, :len(potencias)] = tabela[i - 1, :len(potencias)] * np.maximum(potencias, 0)
            eventos.append(f"linha s^{grau-i} nula: derivada do auxiliar")
        escala = np.max(np.abs(tabela[i]))
        if escala:
            tabela[i] /= escala
        if abs(tabela[i, 0]) <= tol:
            tabela[i, 0] = eps
            eventos.append(f"pivô s^{grau-i} substituído por epsilon")
    sinais = np.signbit(tabela[:, 0])
    return tabela, int(np.count_nonzero(sinais[1:] != sinais[:-1])), eventos


def _modulo_quadrado(coef):
    # P(jw)P(-jw), com as potências pares reescritas em x=w².
    produto = np.polymul(coef, np.asarray(coef) * (-1) ** np.arange(len(coef) - 1, -1, -1))
    return produto[::2] * (-1) ** np.arange(len(coef) - 1, -1, -1)


def _raizes_positivas(coef):
    raizes = np.roots(np.trim_zeros(coef, "f"))
    return sorted(float(z.real) for z in raizes if z.real > 0 and abs(z.imag) <= 1e-7 * (1 + z.real))


def cruzamentos_ganho(malha):
    equacao = np.polysub(_modulo_quadrado(malha.den), _modulo_quadrado(malha.num))
    return np.sqrt(_raizes_positivas(equacao))


def _bissecao(f, a, b):
    fa = f(a)
    for _ in range(48):
        meio = (a + b) / 2
        fm = f(meio)
        if np.signbit(fa) == np.signbit(fm):
            a, fa = meio, fm
        else:
            b = meio
    return (a + b) / 2


def margens(malha, wmin=1e-4, wmax=100.0, pontos=4096):
    """Margens assinadas em todos os cruzamentos encontrados na faixa positiva."""
    wg = [w for w in cruzamentos_ganho(malha) if wmin <= w <= wmax]
    if malha.atraso == 0:
        # Im{N(jw)D(-jw)} é um polinômio ímpar em w.
        n, d = np.asarray(malha.num), np.asarray(malha.den)
        n = n * (1j ** np.arange(len(n) - 1, -1, -1))
        d = d * ((-1j) ** np.arange(len(d) - 1, -1, -1))
        candidatos = np.sqrt(_raizes_positivas(np.polymul(n, d).imag[-2::-2][::-1]))
    else:
        w = np.unique(np.r_[np.geomspace(wmin, wmax, pontos),
                            np.linspace(wmin, wmax, max(2, int(wmax * malha.atraso * 16) + 2))])
        im = malha.avaliar(1j * w).imag
        candidatos = [_bissecao(lambda x: malha.avaliar(1j*x).imag, a, b)
                      for a, b in zip(w[:-1][np.signbit(im[:-1]) != np.signbit(im[1:])],
                                      w[1:][np.signbit(im[:-1]) != np.signbit(im[1:])])]
        candidatos.extend(w[im == 0])
    wf = sorted({float(w) for w in candidatos if wmin <= w <= wmax and malha.avaliar(1j*w).real < 0})
    if malha.den[-1] and (malha.num[-1] / malha.den[-1]) < 0:
        wf.insert(0, 0.0)
    mg = [(w, float(-20 * np.log10(abs(malha.avaliar(1j*w))))) for w in wf]
    mf = [(w, float((180 + np.angle(malha.avaliar(1j*w), deg=True) + 180) % 360 - 180))
          for w in wg]
    return mg, mf


def _raio_externo(malha):
    """Limita |L| no arco direito por desigualdades triangulares."""
    d = np.abs(np.asarray(malha.den) / malha.den[0])
    n = np.abs(np.asarray(malha.num) / malha.den[0])
    r = 4.0
    for _ in range(40):
        x = 1 / r
        inferior = 1 - x * np.polynomial.polynomial.polyval(x, d[1:])
        superior = x ** (len(d) - len(n)) * np.polynomial.polynomial.polyval(x, n)
        if inferior > 0 and superior / inferior < 0.05:
            return r
        r *= 2
    raise ArithmeticError("Não foi possível limitar o arco externo.")


def _contorno(raio, polos, desvio, pontos, atraso, cruzamentos):
    curvas = []
    inicio = -raio
    for centro in [*polos, raio + desvio]:
        fim = min(raio, centro - desvio)
        base = np.linspace(inicio, fim, max(2, int(pontos * (fim-inicio) / (2*raio))))
        proximos = [c + sinal * np.geomspace(desvio, 2*raio, max(20, pontos//6))
                    for c in polos for sinal in (-1, 1)]
        criticos = [w + sinal * np.geomspace(1e-8, raio, max(40, pontos//12))
                    for w in cruzamentos for sinal in (-1, 1)]
        valores = np.concatenate([base, *proximos, *criticos])
        curvas.append(1j * np.unique(valores[(valores >= inicio) & (valores <= fim)]))
        if fim < raio:
            theta = np.linspace(-np.pi/2, np.pi/2, max(64, pontos//8))
            curvas.append(1j*centro + desvio*np.exp(1j*theta))
        inicio = centro + desvio
    n_arco = max(64, pontos//4, int(np.pi*raio*atraso*8)+1)
    curvas.append(raio * np.exp(1j*np.linspace(np.pi/2, -np.pi/2, n_arco)))
    return np.concatenate(curvas)


def nyquist(malha, pontos=384, tol=1e-7):
    """Contorno horário com desvios à direita. Retorna N, P, Z e diagnóstico."""
    polos = np.roots(malha.den)
    p = int(np.sum(polos.real > tol))
    eixo = sorted({float(z.imag) for z in polos if abs(z.real) <= tol})
    raio = _raio_externo(malha)
    separacoes = np.diff([-raio, *eixo, raio])
    desvio = min(1e-2, min(separacoes)/4)
    candidatos = list(cruzamentos_ganho(malha))
    if malha.den[-1] != 0:
        candidatos.append(0.0)
    if any(abs(1 + malha.avaliar(1j*w)) <= tol for w in candidatos):
        return {"estado": "fronteira", "N": None, "P": p, "Z": None, "raio": raio}
    ultimo = None
    repeticoes = 0
    base = max(pontos, int(32*raio*malha.atraso)+1)
    for i in range(7):
        s = _contorno(raio, eixo, desvio / 2 ** (i//2), base * 2**i, malha.atraso, [-w for w in candidatos] + candidatos)
        f = 1 + malha.avaliar(s)
        if not np.all(np.isfinite(f)) or np.min(np.abs(f)) <= tol:
            break
        passos = np.angle(np.roll(f, -1) / f)
        bruto = -np.sum(passos) / (2*np.pi)
        n = int(round(bruto))
        if np.max(np.abs(passos)) < np.pi/3 and abs(bruto-n) < 1e-6 and n+p >= 0:
            repeticoes = repeticoes + 1 if n == ultimo else 1
            if repeticoes >= 3:
                return {"estado": "estável" if n+p == 0 else "instável",
                        "N": n, "P": p, "Z": n+p, "raio": raio,
                        "amostras": len(s), "distancia": float(np.min(np.abs(f)))}
            ultimo = n
        else:
            repeticoes = 0
    return {"estado": "inconclusivo", "N": None, "P": p, "Z": None, "raio": raio}


def lugar_raizes(malha, ganhos):
    """Calcula e conecta polos por proximidade, para malha racional."""
    linhas = []
    for k in ganhos:
        raizes = np.roots(malha.caracteristico(k))
        if linhas:
            dist = np.abs(linhas[-1][:, None] - raizes[None, :])
            ordenadas = np.empty_like(raizes)
            for _ in raizes:
                i, j = np.unravel_index(np.argmin(dist), dist.shape)
                ordenadas[i] = raizes[j]
                dist[i, :] = np.inf
                dist[:, j] = np.inf
            raizes = ordenadas
        linhas.append(raizes)
    return np.asarray(linhas)


def pade(malha, ordem):
    """Aproximação diagonal [m/m] de exp(-tau*s)."""
    if ordem < 1 or malha.atraso <= 0:
        raise ValueError("É necessário atraso positivo e ordem pelo menos 1.")
    coef = np.array([factorial(2*ordem-k)*factorial(ordem)*malha.atraso**k /
                     (factorial(2*ordem)*factorial(k)*factorial(ordem-k))
                     for k in range(ordem+1)])
    return Malha(tuple(np.polymul(malha.num, (coef * (-1)**np.arange(ordem+1))[::-1])),
                 tuple(np.polymul(malha.den, coef[::-1])))

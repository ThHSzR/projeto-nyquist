"""Modelos adotados para o projeto e casos de comparação."""
from math import atan, sqrt

from nucleo import Malha

TAU_CRITICO = atan(4.0)
TAU_D = 1.35
A = Malha((12,), (1, 6, 11, 6))
B = Malha((3,), (1, 1, -2))
C = Malha((-1.5, 1.5), (1, 2, 1))
D = Malha((sqrt(17),), (1, 4, 0), TAU_D)
MODELOS = {"A": A, "B": B, "C": C, "D": D}

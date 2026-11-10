import unittest
from math import sqrt

import numpy as np

from modelos import A, B, C, D, TAU_CRITICO
from nucleo import Malha, lugar_raizes, margens, nyquist, pade, routh


class TestNyquist(unittest.TestCase):
    def test_quatro_modelos(self):
        for m, esperado in [(A, (0, 0, 0)), (B, (1, -1, 0)),
                            (C, (0, 0, 0)), (D, (0, 2, 2))]:
            with self.subTest(m=m):
                diag = nyquist(m)
                self.assertEqual(tuple(diag[x] for x in ("P", "N", "Z")), esperado)

    def test_integrador_e_atraso_critico(self):
        self.assertEqual(nyquist(Malha(D.num, D.den, 1))["estado"], "estável")
        diag = nyquist(Malha(D.num, D.den, TAU_CRITICO))
        self.assertEqual(diag["estado"], "fronteira")
        self.assertIsNone(diag["Z"])

    def test_polos_no_eixo_fora_da_origem(self):
        self.assertEqual(nyquist(Malha((1, 1), (1, 2, 1, 2)))["Z"], 0)

    def test_fronteira_em_dc(self):
        self.assertEqual(nyquist(Malha((2,), B.den))["estado"], "fronteira")

    def test_instavel_em_malha_aberta_e_fechada(self):
        self.assertEqual(nyquist(Malha((1,), B.den))["Z"], 1)

    def test_raio_e_estabilidade(self):
        self.assertGreaterEqual(nyquist(A)["raio"], 4)


class TestRouth(unittest.TestCase):
    def test_trocas(self):
        self.assertEqual(routh(A.caracteristico())[1], 0)
        self.assertEqual(routh([1, -2, -3])[1], 1)

    def test_pivo_zero(self):
        _, trocas, eventos = routh([1, 0, 2, 1])
        self.assertEqual(trocas, 2)
        self.assertIn("epsilon", eventos[0])

    def test_linha_nula(self):
        _, trocas, eventos = routh([1, 5, 11, 23, 28, 12])
        self.assertEqual(trocas, 0)
        self.assertIn("auxiliar", eventos[0])

    def test_origem(self):
        self.assertIn("origem", routh([1, 1, 0])[2][0])


class TestComparacao(unittest.TestCase):
    def test_margens(self):
        mg, mf = margens(A)
        self.assertAlmostEqual(mg[0][0], sqrt(11))
        self.assertAlmostEqual(mg[0][1], 20*np.log10(5))
        self.assertAlmostEqual(mf[0][1], 75.636048, places=5)
        self.assertAlmostEqual(margens(B)[0][0][1], 20*np.log10(2/3))
        self.assertAlmostEqual(margens(D)[1][0][0], 1)

    def test_lugar_das_raizes(self):
        r = lugar_raizes(B, [0, 1])
        np.testing.assert_allclose(np.sort_complex(r[0]), [-2, 1])
        self.assertTrue(np.allclose(np.polyval(B.caracteristico(), r[1]), 0))

    def test_pade_troca_de_diagnostico(self):
        self.assertEqual(routh(pade(D, 1).caracteristico())[1], 0)
        self.assertEqual(routh(pade(D, 2).caracteristico())[1], 2)
        self.assertAlmostEqual(abs(pade(D, 4).avaliar(1j)), abs(D.avaliar(1j)))
        with self.assertRaises(ValueError):
            D.caracteristico()

    def test_malha_impropria(self):
        with self.assertRaises(ValueError):
            Malha((1, 2), (1, 3))


if __name__ == "__main__":
    unittest.main()

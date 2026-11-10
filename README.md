# Estabilidade pelo critério de Nyquist

Projeto de pesquisa aplicada da trilha de Engenharia de Computação. O núcleo é implementado em Python com NumPy. SciPy serve como referência para os ensaios, os degraus contínuos e a discretização. O controlador da bancada é escrito em C++ e tem uma versão para Arduino-ESP32 e outra para execução no computador.

## 1. Preparar o ambiente

Use Python 3.12 e um compilador C++ disponível como `g++`. Abra o terminal na raiz do projeto.

```bash
python -m venv .venv
# Linux/macOS
source .venv/bin/activate
# Windows PowerShell: .\.venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
```

Se a política do PowerShell impedir a ativação, use `.\.venv\Scripts\python.exe` no lugar de `python`. A reprodução da bancada requer que `g++ --version` funcione no mesmo terminal.

## 2. Entender as plantas e repetir os cálculos

Adota-se realimentação negativa unitária. Para uma malha racional L(s) = N(s)/D(s), o ganho K multiplica N(s) e o denominador de malha fechada é D(s) + K N(s). A comparação usa K = 1.

| Caso | Malha aberta L(s) | Polinômio fechado em K = 1 | P, N, Z |
|---|---|---|---|
| A, terceira ordem | 12 / [(s+1)(s+2)(s+3)] | s³ + 6s² + 11s + 18 | 0, 0, 0 |
| B, polo aberto instável | 3 / [(s−1)(s+2)] | s² + s + 1 | 1, −1, 0 |
| C, fase não mínima | 1,5(1−s) / (s+1)² | s² + 0,5s + 2,5 | 0, 0, 0 |
| D, atraso puro | √17 exp(−τs) / [s(s+4)] | Equação transcendental, τ = 1,35 s | 0, 2, 2 |

### A. Routh com as contas abertas

1. Multiplique (s+1)(s+2) = s²+3s+2. Em seguida, multiplique por (s+3) para obter s³+6s²+11s+6.
2. Some 12K ao termo constante. As linhas iniciais de Routh são [1, 11] e [6, 6+12K].
3. A linha s tem primeiro termo (6×11 − (6+12K))/6 = 10−2K. A linha constante é 6+12K.
4. Com K positivo, todas as entradas ficam positivas quando 0<K<5. Em K=1, a coluna bruta é [1, 6, 8, 18]. O programa divide cada linha por sua magnitude máxima, mantendo os sinais.
5. Em ω = √11, a parte imaginária do denominador se anula e D(jω) = −60. Portanto, L(jω) = −0,2 e MG = 20 log₁₀(5) = 13,98 dB. A frequência de ganho unitário é 1,22318 rad/s e a MF é 75,64°.

### B. Malha aberta instável, malha fechada estável

1. Expanda (s−1)(s+2) = s²+s−2. Somando 3K, a primeira coluna de Routh fica [1, 1, 3K−2].
2. O termo constante é positivo para K>2/3. Com K=1, as raízes são −0,5 ± j√3/2.
3. Como há um polo aberto à direita, P=1. Nyquist dá uma volta anti-horária, N=−1. Logo, Z=1−1=0.
4. A margem de ganho de −3,52 dB vem de 20 log₁₀(2/3) no limite em corrente contínua. Aqui, perder ganho é que conduz à instabilidade.

### C. Zero à direita e resposta inversa

1. Expanda (s+1)² = s²+2s+1. Some 1,5K(1−s). A coluna de Routh fica [1, 2−1,5K, 1+1,5K].
2. Para K>0, o intervalo estável é 0<K<4/3. Em K=1, o zero fica em +1 e os polos fechados têm parte real −0,25.
3. No degrau, a inclinação inicial é lim s→∞ sT(s) = −1,5. O valor final é T(0)=1,5/2,5=0,6. Por isso, a curva desce antes de subir.

### D. Atraso crítico com as contas abertas

1. Calcule |L(jω)| = √17 / [ω√(ω²+16)]. Em ganho unitário, ω⁴+16ω²−17=0. Fatorando, (ω²−1)(ω²+17)=0 e ωg=1 rad/s.
2. Em ωg, a fase da parte racional é −π/2−arctan(1/4). O atraso acrescenta −ωgτ = −τ rad.
3. Iguale a fase a −π. Resulta τcr = π/2−arctan(1/4) = arctan(4) = **1,3258176637 s**. Em τcr, os polos cruzam o eixo em ±j.
4. Para 1,35 s, MF = (1,3258176637−1,35)×180/π = −1,38555°. Nyquist dá N=2, P=0, Z=2.
5. Em Padé [1/1], a perda de fase é 2 arctan(τ/2). Igualar à margem racional fornece 2 tan(arctan(4)/2) = 1,5615528128 s. O erro relativo é (1,5615528128−1,3258176637)/1,3258176637 × 100 = **17,78036%**.
6. Em 1,35 s, Padé [1/1] ainda dá estável. Padé [2/2] já dá dois polos à direita. Consulte `resultados/pade.csv` para as ordens 1 a 4.

## 3. Rodar os testes e reproduzir as figuras

```bash
python -m unittest discover -s tests -v
python analisar.py --caso B
python analisar.py --caso D --ordem-pade 1
python analisar.py --num 4.123105625617661 --den 1 4 0 --atraso 1.35
python experimento.py --saida repeticao
```

A pasta `resultados/` guarda a execução usada no artigo. A nova pasta `repeticao/` preserva esses dados e permite comparar os tempos no seu computador. Todos os coeficientes informados por CLI devem estar em ordem decrescente e a função precisa ser estritamente própria.

| Arquivo | Leitura |
|---|---|
| `nyquist.png` | Ramos em frequência e ponto −1+j0. O algoritmo fecha o contorno e faz desvios no plano s. |
| `bode.png` | Magnitude, fase e frequência de ganho unitário. A fase de D é ampliada nas frequências de interesse. |
| `lugar_raizes.png` | Caminho dos polos com a variação de K. Em D, usa Padé [2/2]. |
| `degrau.png` | Degraus contínuos. Em D, usa Padé [2/2]. |
| `atraso_pade.png` | Zoom do ramo positivo de Nyquist perto do ponto crítico. |
| `benchmark.png` | Mediana do tempo de CPU versus ordem, com eixo vertical logarítmico. |
| `bancada.png` | Degraus da bancada em software usando a fila de atraso e o PID em C++. |
| `validacao.csv` e `atraso_exato.csv` | Contagens, margens, polos e resíduos. A coluna de Routh em D usa Padé [2/2]. |
| `benchmark.csv` | Mediana, percentil 90, erro direto dos polos, resíduo e diferença float32/float64. |

**Limites do programa.** Os casos próximos do ponto crítico podem retornar fronteira ou inconclusivo, de acordo com a tolerância numérica. O módulo de margens retorna cruzamentos na faixa de 0,0001 a 100 rad/s, mais o limite de fase em DC quando aplicável. Um cruzamento tangencial com atraso pode escapar à varredura por troca de sinal. A conexão de ramos do lugar das raízes por proximidade pode ser ambígua em raízes coincidentes.

## 4. Ensaiar a bancada em software

```bash
g++ -std=c++11 -Wall -Wextra -Werror -pedantic -O2 hil/pid_host.cpp -o pid_host
python bancada.py --executavel ./pid_host --atraso 1.00 --saida software_100.csv
python bancada.py --executavel ./pid_host --atraso 1.32 --saida software_132.csv
python bancada.py --executavel ./pid_host --atraso 1.50 --saida software_150.csv
python comparar_hil.py software_100.csv software_132.csv software_150.csv --saida software.png
```

No Windows, use `-o pid_host.exe` e `--executavel .\pid_host.exe`. O host e o ESP32 compartilham o mesmo cabeçalho PID e o mesmo parser. Os ensaios registrados usaram o processo local, com 3.001 amostras de 20 ms em tempo simulado por cenário. A latência do processo local é distinta da latência serial da placa.

## 5. Realizar a etapa física no ESP32

1. Instale o pacote **esp32 by Espressif Systems** no Arduino IDE e abra `hil/esp32_pid/esp32_pid.ino`.
2. Selecione a placa correspondente ao seu ESP32 e a porta USB. Compile, grave e feche o Monitor Serial.
3. Troque `COM4` pela porta detectada. No Linux, um exemplo é `/dev/ttyUSB0`. Execute uma rodada de 60 s por atraso.

```bash
python bancada.py --porta COM4 --atraso 1.00 --saida esp32_100.csv
python bancada.py --porta COM4 --atraso 1.30 --saida esp32_130.csv
python bancada.py --porta COM4 --atraso 1.32 --saida esp32_132.csv
python bancada.py --porta COM4 --atraso 1.50 --saida esp32_150.csv
python comparar_hil.py esp32_100.csv esp32_130.csv esp32_132.csv esp32_150.csv --saida esp32.png
```

4. Registre modelo da placa, versão do core Arduino-ESP32, porta, RTT p50 e p95 e razão dos envelopes. O atraso discreto previsto é 1,30514 s. 1,30 s e 1,32 s ficam em lados opostos desse limiar, mas o efeito é pequeno perto da fronteira. Use 1,00 s e 1,50 s para um contraste visível.
5. Se houver perda de prazo de amostra, a rodada é abortada. Feche programas que usem a porta, verifique o cabo e tente novamente. Alterar o período exige modificar o firmware, o host e recalcular o limiar.

**Status experimental:** controlador e protocolo verificados em software. A compilação para a placa, gravação e execução física precisam ser feitas com o ESP32 disponível. Os gráficos fornecidos são de simulação.

## 6. Artigo, seminário e publicação

- `artigo/` contém o fonte LaTeX, a classe IEEEtran e o PDF. Em uma instalação TeX Live completa, execute `pdflatex artigo_nyquist.tex` duas vezes a partir de `artigo/`.
- `seminario.pptx` contém as falas nas notas de cada slide. A demonstração mostra a troca de diagnóstico em D e, com a placa disponível, compara dois atrasos em tempo real.
- O pacote inclui um Git bundle. Para recuperar os commits, execute `git clone nyquist.bundle nyquist-versionado` na pasta que contém o bundle. Para publicar na conta GitHub autenticada, entre na pasta clonada e execute os comandos abaixo.

```bash
git log --oneline
gh auth status
gh repo create projeto-nyquist --public --source=. --remote=github --push
```

### Referências

- Åström e Murray, [Feedback Systems, Frequency Domain Analysis](https://fbswiki.org/wiki/index.php/Frequency_Domain_Analysis).
- Ogata, *Engenharia de Controle Moderno*, 5. ed., Pearson, 2010.
- MathWorks, [Time-Delay Approximation](https://www.mathworks.com/help/control/ug/time-delay-approximation.html).
- SciPy, [Signal processing](https://docs.scipy.org/doc/scipy/reference/signal.html) e [companion](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.companion.html).
- NumPy, [roots](https://numpy.org/doc/stable/reference/generated/numpy.roots.html).
- Espressif Systems, [Serial UART para Arduino-ESP32](https://docs.espressif.com/projects/arduino-esp32/en/latest/api/serial.html).
- IEEEtran v1.8b, [classe e condições de licença](https://ctan.org/pkg/ieeetran). A cópia incluída mantém a licença original.

O código original usa a licença MIT. A classe IEEEtran é um componente de terceiros com a licença própria declarada no arquivo `artigo/IEEEtran.cls`.

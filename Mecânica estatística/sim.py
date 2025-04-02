import random
import math

# Constantes físicas básicas
c_na = 6.02214179e23
c_mu = 1.660538921e-27
c_kb = 1.3806488e-23
c_pi = 3.14159265359

# Modelo do Argônio
sigma = 3.41e-10  # sigma
epsilon = 119.8 * c_kb  # epsilon
rc = 3.2 * sigma  # Distância de corte do potencial de Lennard-Jones
m = 39.944 * c_mu  # Massa de um átomo de argônio

# Definindo a região de simulação
n = 500  # Número de partículas na região de simulação
p = [[0.0] * n for _ in range(3)]  # Posicionamento das partículas
teplota = 100.0  # Temperatura
beta = 1.0 / (teplota * c_kb)  # Beta, utilizado no critério de aceitação do MC
h = [0.0, 0.0, 0.0]  # Dimensões da caixa cúbica da simulação

# Parâmetros de Monte Carlo
max_posun = 0.2 * sigma  # Máximo deslocamento por movimento Monte Carlo

# Fracções de tentativas e aceites Monte Carlo
zp_pokusu_p = 0.0
zp_prijato_p = 0.0

# Parâmetros para o perfil de densidade
hp_n = 100
hp_sum = [0.0] * (2 * hp_n + 1)
hp_num = 0

def pocitej_hustotni_profil():
    """Calcula o perfil de densidade."""
    x_cm = sum(p[0]) / n  # Coordenada x do centro de massa
    for c in range(n):
        dx = p[0][c] - x_cm
        if dx > 0.5 * h[0]:
            dx -= h[0]
        if dx < -0.5 * h[0]:
            dx += h[0]
        s = round(dx * hp_n / (0.5 * h[0]))
        if -hp_n <= s <= hp_n:
            hp_sum[s] += 1
    hp_num += 1

def posun_nahodnou_castici():
    """Realiza o movimento de uma partícula aleatória e aplica o critério de aceitação de Monte Carlo."""
    global zp_pokusu_p, zp_prijato_p
    zp_pokusu_p += 1

    # Seleção aleatória da partícula
    c = random.randint(0, n - 1)

    # Geração de um movimento aleatório para a partícula c
    pz = [p[0][c] + (random.random() - 0.5) * max_posun,
          p[1][c] + (random.random() - 0.5) * max_posun,
          p[2][c] + (random.random() - 0.5) * max_posun]

    # Aplicação das condições de contorno periódicas
    for k in range(3):
        if pz[k] > 0.5 * h[k]:
            pz[k] -= h[k]
        if pz[k] < -0.5 * h[k]:
            pz[k] += h[k]

    # Cálculo da energia correspondente à mudança
    e = 0
    ez = 0
    for i in range(n):
        if i == c:
            continue
        e += e_ij(p[0][c], p[0][i], h)
        ez += e_ij(pz, p[0][i], h)

    # Critério de aceitação de Monte Carlo
    if random.random() < math.exp(-beta * (ez - e)):
        p[0][c] = pz[0]
        p[1][c] = pz[1]
        p[2][c] = pz[2]
        zp_prijato_p += 1

def e_pot(p, h):
    """Calcula a energia potencial total."""
    e_pot = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            e_pot += e_ij(p[i], p[j], h)
    return e_pot

def e_ij(p1, p2, h):
    """Calcula a energia entre dois átomos."""
    r = [p2[i] - p1[i] for i in range(3)]
    d2 = sum((0.5 * h[i] - abs(0.5 * h[i] - abs(r[i]))) ** 2 for i in range(3))

    if d2 < rc ** 2:
        return 4 * epsilon * ((sigma ** 2 / d2) ** 6 - (sigma ** 2 / d2) ** 3)
    else:
        return 0


# Inicializando o gerador de números aleatórios
random.seed()

print('Iniciando a simulação do líquido na região cúbica')

# Definindo as dimensões da caixa cúbica, correspondente à densidade do líquido
h = [(n * m / 1400) ** (1.0 / 3)] * 3

# Testando a validade da distância de corte
if rc > 0.5 * h[0]:
    print('A aresta da região de simulação é muito curta ou rc é muito grande')
    exit()

# Gerando a configuração inicial
for i in range(n):
    p[0][i] = random.random() * h[0]
    p[1][i] = random.random() * h[1]
    p[2][i] = random.random() * h[2]

# Zerando os contadores do perfil de densidade
hp_sum = [0.0] * (2 * hp_n + 1)
hp_num = 0

# Simulação Monte Carlo
for i in range(1, 10000000 + 1):
    zp_pokusu_p = 0
    zp_prijato_p = 0

    for j in range(1, 10000 + 1):
        # Função para mover partículas (defina a função posun_nahodnou_castici abaixo)
        posun_nahodnou_castici()

    print(f"i={i}   e_pot={e_pot(p, h) / n * c_na:.3f} zp_p={zp_prijato_p / zp_pokusu_p:.3f}")

    # Expansão da região de simulação (após um equilíbrio razoável)
    if i == 100:
        h[0] = 3 * h[0]
        print('A região de simulação foi expandida e segue a simulação da interface de fase')

    # Cálculo do perfil de densidade
    if i > 100:
        pocitej_hustotni_profil()

    # Gravação dos resultados do perfil de densidade e configuração para animação
    if i % 100 == 0:
        with open('Simulacao6_hp.txt', 'w') as f:
            for j in range(-hp_n, hp_n + 1):
                f.write(f"{j * h[0] / (2 * hp_n)} {m * hp_sum[j] / hp_num / (h[0] / (2 * hp_n) * h[1] * h[2]):.3f}\n")

        # Gravação da configuração da simulação para visualização no VMD
        with open('Simulacao6_konfiguracao.xyz', 'a') as f:
            f.write(f"{n}\n")
            f.write('Simulação da interface de fase do Argônio\n')
            for j in range(n):
                f.write(f"Ar {p[0][j] * 1e10} {p[1][j] * 1e10} {p[2][j] * 1e10}\n")

    # Reinicialização do cálculo do perfil de densidade após atingir o equilíbrio
    if i == 1000:
        hp_sum = [0.0] * (2 * hp_n + 1)
        hp_num = 0

# Funções auxiliares




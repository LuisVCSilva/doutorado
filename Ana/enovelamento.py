from copy import deepcopy
from math import exp
from random import randrange, random
import matplotlib.pyplot as plt

def nova_conformacao_estendida(dna):
    """ Entrada: Sequência de DNA: '0' para resíduo hidrofóbico e '1' para hidrofílico
        Saída: conformação linear estendida
    """
    L = len(dna)
    assert L == dna.count('0') + dna.count('1')
    return [[0, i] for i in range(L)]

def E(s, dna):
    """ Entrada: sequência de DNA 
        Saída: E = número de vizinhos próximos não-covalentes do centro da proteína
    """
    L = len(s)
    n_BB = 0
    for i in range(L - 2):
        for j in range(i + 2, L):
            if dna[i] == '0' and dna[j] == '0':
                dx = s[i][0] - s[j][0]
                dy = s[i][1] - s[j][1]
                if abs(dx) + abs(dy) == 1:
                    n_BB += 1
    return -n_BB

def MMC_mutation_step(s, dna, c_k):
    """ Entrada: conformação e temperatura
        Saída: conformação alterada
    """
    L = len(s)
    s_trial = deepcopy(s)          
    i = randrange(1, L - 1)     
    matrizes_rotacao = [[[0, -1], [1,  0]], [[-1, 0], [0, -1]], [[0, 1], [-1, 0]]]
    rot_mat = matrizes_rotacao[randrange(3)]
    
    # Rotaciona
    for j in range(i):
        dr = [s[j][k] - s[i][k] for k in range(2)]
        for k in range(2):
            s_trial[j][k] = s[i][k]
            for l in range(2):
                s_trial[j][k] += rot_mat[k][l] * dr[l]
    
    # Verifica se é uma conformação válida
    for j in range(i):
        for k in range(i + 1, len(s)):
            if s_trial[j] == s[k]:
                return False
    # Conformação válida, atualiza com algoritmo Metropolis
    dE = E(s_trial, dna) - E(s, dna)
    w = exp(-dE / c_k)
    if random() < w:
        for j in range(i):
            for k in range(2):
                s[j][k] = s_trial[j][k]
    return True

def recozimento_simulado(s, dna, iteracoes):
    with open('output.txt', 'w') as f:
        """ Entrada: conformação s, dna, quantidade de passos para Metropolis/Monte Carlo
            Saída: energia mínima, iterações gastas, conformações válidas
        """
        c_k = 2.0           # temperatura inicial
        alpha = 0.99        # taxa de resfriamento
        valid_confs = 0
        E_min = E(s, dna)
        step_min = 0
        for step in range(1, iteracoes + 1):
            if MMC_mutation_step(s, dna, c_k):
                print(str((s, dna, c_k)), file=f)
                valid_confs += 1
                E_step = E(s, dna)
                if E_step < E_min:
                    E_min = E_step
                    step_min = step
                    print(" E =", E_min, " >> passo = ", step)
            if step % 200000 == 0 and c_k > 0.15:
                c_k *= 0.99
        valid_fraction = valid_confs / float(iteracoes)
        return E_min, step_min, valid_fraction

dna = "11010101100101100101001100101001100101001001"
proteina = nova_conformacao_estendida(dna)
L = len(proteina)
print("0 -> Hidrofóbico")
print("1 -> Hidrofílico")
print("Cadeia de peptídeo = ", dna, "L =", L)
iteracoes = 50 * 10**2
print("Recozimento simulado para ", iteracoes, " iterações")
E_min, step_min, valid_confs = recozimento_simulado(proteina, dna, iteracoes)
print(" E_min =", E_min, " na iteração ", step_min, " fração de configurações válidas =", valid_confs)

x = []
y = []
cores = []  # Lista para armazenar as cores dos vértices

for i, iproteina in enumerate(proteina):
    x.append(iproteina[0])
    y.append(iproteina[1])
    cores.append('blue' if dna[i] == '0' else 'red')  # 'blue' para hidrofóbico, 'red' para hidrofílico

# Plota as conformações com linhas conectando os vértices
for i in range(len(x) - 1):
    plt.plot([x[i], x[i + 1]], [y[i], y[i + 1]], color='black')  # Linha preta para conectar vértices

plt.scatter(x, y, c=cores, marker='o')  # Usa cores diferentes para cada tipo de resíduo
plt.xlabel('x')
plt.ylabel('y')
plt.title('Conformação da Proteína')
plt.show()

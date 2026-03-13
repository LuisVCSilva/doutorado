#modelo HP para protein folding
#R. Unger and J. Moult, J. Mol. Biol. 231, 75-81 (1993)
#https://en.wikipedia.org/wiki/Hydrophobic-polar_protein_folding_model
#dada uma cadeia de peptídeos e uma funcao de energia do sistema, minimiza a funcao computando a configuracao de energia minima/configuracao funcional
from copy import deepcopy
from math import exp
from random import randrange, random
import matplotlib.pyplot as plt

def nova_conformacao_estendida(dna):
    """ entrada:  Sequência de DNA: '0' para residuo hidrofobico e '1'  para hidrofilico
        saida: conformacao linear estendida
    """
    L = len(dna)
    assert L==dna.count('0') + dna.count('1')
    return [[0, i] for i in range(L)]

def E(s, dna):
    """ entrada:  sequencia de dna 
        saida: E = numero de vizinhos proximos nao-covalentes do centro da proteina
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
    """ entrada:  conformaçao e temperatura
        output: conformacao alterada
    """
    L = len(s)
    s_trial = deepcopy(s)          
    i = randrange(1, L - 1)     
    matrizes_rotacao = [[[0, -1],[1,  0]], [ [-1, 0],
                       [0, -1]], [ [ 0, 1],[-1, 0]]]
    rot_mat = matrizes_rotacao[randrange(3)]
    
    #rotaciona
    for j in range(i):
        dr = [ s[j][k] - s[i][k] for k in range(2) ]
        for k in range(2):
            s_trial[j][k] = s[i][k]
            for l in range(2):
                s_trial[j][k] += rot_mat[k][l] * dr[l]
    
    #verifica se eh conformacao valida
    for j in range(i):
        for k in range(i + 1, len(s)):
            if s_trial[j] == s[k]:
                return False
    #conformacao valida, atualiza com algoritmo Metropolis
    dE = E(s_trial, dna) - E(s, dna)
    w = exp(- dE / c_k)
    if random() < w:
        for j in range(i):
            for k in range(2):
                s[j][k] = s_trial[j][k]
    return True

def recozimento_simulado(s, dna, iteracoes):
    with open('output.txt', 'w') as f:
        """ entrada:  conformacao s, dna, qtde de passos para metropolis/monte carlo
            saida: energia minima, iteracoes gastas, conformacoes validas
        """
        c_k = 2.0           # temperatura inicialo
        alpha = 0.99        # taxa de resfriamento
        valid_confs = 0
        E_min = E(s, dna)
        step_min = 0
        for step in range(1, iteracoes + 1):
            if MMC_mutation_step(s, dna, c_k):
                print(str((s,dna,c_k)),file=f)
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


dna = "01011001011010011010"
proteina = nova_conformacao_estendida(dna)
L = len(proteina)
print("0 -> Hidrofobico")
print("1 -> Hidrofobico")
print("Cadeia de peptideo = ", dna, "L =", L)
iteracoes = 50 * 10**6
print("Recozimento simulado para ", iteracoes, " iteracoes")
E_min, step_min, valid_confs = recozimento_simulado(proteina, dna, iteracoes)
print(" E_min =", E_min, " na iteracao ", step_min, "fracao de configuracoes validas =", valid_confs)

x = []
y = []
for iproteina in proteina:
    x.append( iproteina[0] )
    y.append( iproteina[1] )

for i in range(len(x)):
    plt.plot( x,y , marker='o')
    
plt.show()

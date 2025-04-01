import numpy as np
import matplotlib.pyplot as plt

def caminhada_aleatoria(N, p, l=1, tau=1):
    posicao = 0
    trajetoria = [posicao]
    
    for _ in range(N):
        passo = l if np.random.rand() < p else -l
        posicao += passo
        trajetoria.append(posicao)
    
    tempo = np.arange(0, (N + 1) * tau, tau)
    return tempo, trajetoria

# Parâmetros da simulação
N = 1000      # Número de passos
p = 0.5      # Probabilidade de ir para a direita
l = 1        # Tamanho do passo
tau = 1      # Intervalo entre os passos

# Executando a simulação
tempo, trajetoria = caminhada_aleatoria(N, p, l, tau)

# Plotando o resultado
plt.figure(figsize=(10, 6))
plt.plot(tempo, trajetoria, label=f"Caminhada Aleatória (p={p})")
plt.xlabel("Tempo")
plt.ylabel("Posição")
plt.title("Simulação de Caminhada Aleatória Unidimensional")
plt.legend()
plt.grid(True)
plt.show()


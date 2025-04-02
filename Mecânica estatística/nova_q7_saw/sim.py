import math
import matplotlib.pyplot as plt
from tqdm import tqdm
import numpy as np

n_passos = 30
n_caminhadas = 10  
intervalo_plotagem = 5

caminhadas = 0
caminhadas_falhas = 0
valores_r2 = []
valores_theta = []  
valores_r = [] 

caminhos = []  
barra_progresso = tqdm(total=n_caminhadas, desc="Gerando caminhadas")

while caminhadas < n_caminhadas:
    locais = set()  
    x, y = 0, 0  
    locais.add((x, y))
    caminhada_falhou = False
    caminho = [(x, y)]
    direcoes = np.random.rand(n_passos)
    for passo in range(n_passos):
        direcao = direcoes[passo]
        if direcao < 0.25:
            x += 1  # Direita
        elif direcao < 0.50:
            y += 1  # Cima
        elif direcao < 0.75:
            x -= 1  # Esquerda
        else:
            y -= 1  # Baixo
        
        if (x, y) in locais:
            caminhada_falhou = True
            caminhadas_falhas += 1
            break
        
        locais.add((x, y))
        caminho.append((x, y))
        

        r2 = x**2 + y**2
        theta = math.atan2(y, x)  
        r = math.sqrt(r2)  
            
        valores_r2.append(r2)
        valores_theta.append(theta)
        valores_r.append(r)
        '''
        if passo % intervalo_plotagem == 0:     
            plt.figure(figsize=(8, 5))
            x_vals, y_vals = zip(*caminho)
            plt.plot(x_vals, y_vals, marker='o', markersize=3, linestyle='-', alpha=0.7)
            plt.title(f"Caminhada {caminhadas + 1} - Energia (r²) = {r2:.2f}")
            plt.xlabel("Posição x")
            plt.ylabel("Posição y")
            plt.grid()
            plt.savefig(f"caminhada_{caminhadas + 1}_passo_{passo + 1}.png")
            plt.close()
         '''
    if caminhada_falhou:
        continue
    
    caminhadas += 1
    barra_progresso.update(1)

barra_progresso.close()

plt.figure(figsize=(8, 5))
plt.hist(valores_r, bins=30, edgecolor='black', alpha=0.7, density=True)
plt.title("Distribuição da Distância r")
plt.xlabel("r (Distância da origem)")
plt.ylabel("Densidade")
plt.grid()
plt.show()


plt.figure(figsize=(8, 5))
plt.hist(valores_theta, bins=30, edgecolor='black', alpha=0.7, density=True)
plt.title("Distribuição do Ângulo θ")
plt.xlabel("θ (Ângulo com a horizontal)")
plt.ylabel("Densidade")
plt.grid()
plt.show()


r2_medio = sum(valores_r2) / n_caminhadas
coef_difusao = r2_medio / (4 * n_passos)  
print(f"Coeficiente de Difusão D ≈ {coef_difusao:.4f}")

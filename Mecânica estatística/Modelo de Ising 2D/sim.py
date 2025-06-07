import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import json
import csv
import os

# Carregar parâmetros
with open("parametros.json", "r") as f:
    params = json.load(f)

interacao_J = params["interacao_J"]
campo_H = params["campo_H"]
temperatura = params["temperatura"]
tam_x = params["tam_x"]
tam_y = params["tam_y"]
passos_mc = params["passos_mc"]
intervalo_amostra = params["intervalo_amostra"]

n_spins = tam_x * tam_y
rede = np.random.choice([-1, 1], size=(tam_x, tam_y))

# Reservado para futura extensão para spins com três estados
fatores_boltzmann = np.zeros((17, 3))
for delta in range(-8, 9, 4):
    fatores_boltzmann[delta + 8][0] = np.exp(-(delta * interacao_J + 2 * campo_H) / temperatura)
    fatores_boltzmann[delta + 8][2] = np.exp(-(delta * interacao_J - 2 * campo_H) / temperatura)

def passo_metropolis():
    i = np.random.randint(0, tam_x)
    j = np.random.randint(0, tam_y)
    vizinhos = (
        rede[(i + 1) % tam_x, j] +
        rede[(i - 1) % tam_x, j] +
        rede[i, (j + 1) % tam_y] +
        rede[i, (j - 1) % tam_y]
    )
    delta_energia = 2 * rede[i, j] * vizinhos
    prob = fatores_boltzmann[delta_energia + 8][1 + rede[i, j]]
    if np.random.rand() < prob:
        rede[i, j] *= -1

def passo_monte_carlo():
    for _ in range(n_spins):
        passo_metropolis()

def magnetizacao():
    return np.sum(rede) / n_spins

def energia():
    soma_spins = np.sum(rede)
    interacoes = (
        np.sum(rede * np.roll(rede, 1, axis=0)) +
        np.sum(rede * np.roll(rede, 1, axis=1))
    )
    return -(interacao_J * interacoes + campo_H * soma_spins) / n_spins

def autocorrelacao(x, max_lag):
    x = np.array(x)
    media = np.mean(x)
    var = np.var(x)
    n = len(x)
    return [
        np.correlate(x[:n-lag] - media, x[lag:] - media)[0] / ((n - lag) * var)
        for lag in range(max_lag)
    ]

magnetizacoes = []
energias = []

# Visualização ao vivo
figura, (ax_rede, ax_graficos) = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1, 1.5]})
imagem = ax_rede.imshow(rede, cmap='coolwarm', vmin=-1, vmax=1)
ax_rede.set_title(f"Configuração dos Spins (T = {temperatura})")

ax_graficos.set_title("Magnetização e Energia")
ax_graficos.set_xlabel("Iteração")
linha_mag, = ax_graficos.plot([], [], label=r"$\langle M \rangle$", color='blue')
linha_ene, = ax_graficos.plot([], [], label=r"$\langle E \rangle$", color='orange')
ax_graficos.legend()
ax_graficos.set_xlim(0, passos_mc)
ax_graficos.set_ylim(-2, 2)

def inicializar():
    linha_mag.set_data([], [])
    linha_ene.set_data([], [])
    return imagem, linha_mag, linha_ene

def atualizar(frame):
    passo_monte_carlo()
    if frame % intervalo_amostra == 0:
        mag = magnetizacao()
        ene = energia()
        magnetizacoes.append(mag)
        energias.append(ene)

        imagem.set_data(rede)
        linha_mag.set_data(range(len(magnetizacoes)), magnetizacoes)
        linha_ene.set_data(range(len(energias)), energias)
    return imagem, linha_mag, linha_ene

anim = animation.FuncAnimation(
    figura,
    atualizar,
    frames=passos_mc,
    init_func=inicializar,
    interval=50,
    blit=True,
    repeat=False
)

saida_base = f"saida_T{temperatura:.2f}".replace('.', '_')
os.makedirs(saida_base, exist_ok=True)

print("Salvando vídeo...")
writer = animation.FFMpegWriter(fps=20)
anim.save(f"{saida_base}/animacao.mp4", writer=writer)

print("Salvando arquivos CSV...")
# Energia e magnetização
with open(f"{saida_base}/energia.csv", "w", newline="") as f_ene:
    writer = csv.writer(f_ene)
    writer.writerow(["Iteracao", "Energia"])
    writer.writerows(enumerate(energias))

with open(f"{saida_base}/magnetizacao.csv", "w", newline="") as f_mag:
    writer = csv.writer(f_mag)
    writer.writerow(["Iteracao", "Magnetizacao"])
    writer.writerows(enumerate(magnetizacoes))

# Médias
media_magnetizacao = np.mean(magnetizacoes)
media_energia = np.mean(energias)
with open(f"{saida_base}/magnetizacao_media.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["MagnetizacaoMedia"])
    writer.writerow([media_magnetizacao])

# Calor específico
energias_np = np.array(energias)
var_E = np.var(energias_np)
calor_especifico = (n_spins / temperatura**2) * var_E
with open(f"{saida_base}/calor_especifico.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["CalorEspecifico"])
    writer.writerow([calor_especifico])
plt.figure()
plt.plot(energias_np, label="Energia")
plt.title("Energia ao Longo do Tempo")
plt.savefig(f"{saida_base}/energia_plot.png")

# Susceptibilidade magnética
magnetizacoes_np = np.array(magnetizacoes)
var_M = np.var(magnetizacoes_np)
susceptibilidade = (n_spins / temperatura) * var_M
with open(f"{saida_base}/susceptibilidade.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Susceptibilidade"])
    writer.writerow([susceptibilidade])
plt.figure()
plt.plot(magnetizacoes_np, label="Magnetização")
plt.title("Magnetização ao Longo do Tempo")
plt.savefig(f"{saida_base}/magnetizacao_plot.png")

# Autocorrelação
lags = 100
autocorr_mag = autocorrelacao(magnetizacoes_np, lags)
with open(f"{saida_base}/autocorrelacao_magnetizacao.csv", "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["Lag", "Autocorrelacao"])
    writer.writerows(zip(range(lags), autocorr_mag))
plt.figure()
plt.plot(range(lags), autocorr_mag)
plt.title("Autocorrelação da Magnetização")
plt.savefig(f"{saida_base}/autocorrelacao_plot.png")

# Histograma da magnetização
plt.figure()
plt.hist(magnetizacoes_np, bins=40, density=True)
plt.title("Distribuição da Magnetização")
plt.savefig(f"{saida_base}/histograma_magnetizacao.png")

print("Todas as medidas foram salvas.")

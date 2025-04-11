import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


interacao_J = 1.0
campo_H = 0.0
temperatura = 2.2  
tam_x = tam_y = 50
n_spins = tam_x * tam_y
passos_mc = 500
intervalo_amostra = 1


rede = np.random.choice([-1, 1], size=(tam_x, tam_y))


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


magnetizacoes = []
energias = []
configuracoes = []


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
        configuracoes.append(rede.copy())

        imagem.set_data(rede)
        linha_mag.set_data(range(len(magnetizacoes)), magnetizacoes)
        linha_ene.set_data(range(len(energias)), energias)
    return imagem, linha_mag, linha_ene

animacao = animation.FuncAnimation(
    figura,
    atualizar,
    frames=passos_mc,
    init_func=inicializar,
    interval=50,
    blit=True,
    repeat=False
)

plt.tight_layout()
plt.show()


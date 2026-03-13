import numpy as np
import matplotlib.pyplot as plt
import os

# ----------------------------
# Criar pasta de saída
# ----------------------------
output_dir = "lbm_frames"
os.makedirs(output_dir, exist_ok=True)

# ----------------------------
# Parâmetros
# ----------------------------
Nx = 200
Nt = 300
tau = 0.8

c = np.array([0, 1, -1])
w = np.array([2/3, 1/6, 1/6])

# ----------------------------
# Inicialização
# ----------------------------
rho = np.ones(Nx)
rho[Nx//3:Nx//2] = 2.0

u = np.zeros(Nx)

f = np.zeros((3, Nx))

for i in range(3):
    f[i,:] = w[i]*rho

# ----------------------------
# Loop temporal
# ----------------------------
for t in range(Nt):

    # Grandezas macroscópicas
    rho = np.sum(f, axis=0)

    momentum = np.sum(c[:,None]*f, axis=0)
    u = momentum / rho

    e = np.sum(((c[:,None]-u)**2)*f, axis=0) / rho

    # Equilíbrio
    feq = np.zeros_like(f)

    for i in range(3):
        feq[i,:] = w[i]*rho*(1 + 3*c[i]*u)

    # Colisão
    f = f - (1/tau)*(f - feq)

    # Streaming
    f[1,:] = np.roll(f[1,:], 1)
    f[2,:] = np.roll(f[2,:], -1)

    # ----------------------------
    # Heatmap da densidade
    # ----------------------------
    heatmap = np.tile(rho, (20,1))

    plt.figure(figsize=(8,2))
    plt.imshow(heatmap, aspect='auto')
    plt.colorbar(label='densidade')
    plt.title(f"Densidade t = {t}")
    plt.xlabel("x")

    # salvar frame
    plt.savefig(f"{output_dir}/frame_{t:04d}.png")
    plt.close()

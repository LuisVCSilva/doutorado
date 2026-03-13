import numpy as np
import matplotlib.pyplot as plt
import os

# ----------------------------
# Pasta
# ----------------------------
output_dir = "lbm_multiphase"
os.makedirs(output_dir, exist_ok=True)

# ----------------------------
# Parâmetros
# ----------------------------
Nx, Ny = 100, 100
Nt = 400
tau = 1.0
G = -5.0

# ----------------------------
# D2Q9
# ----------------------------
c = np.array([
    [0,0],
    [1,0],[0,1],[-1,0],[0,-1],
    [1,1],[-1,1],[-1,-1],[1,-1]
])

w = np.array([
    4/9,
    1/9,1/9,1/9,1/9,
    1/36,1/36,1/36,1/36
])

# ----------------------------
# Inicialização
# ----------------------------
rho = 1 + 0.01*np.random.randn(Nx,Ny)

f = np.zeros((9,Nx,Ny))

for i in range(9):
    f[i] = w[i]*rho

# ----------------------------
# Loop
# ----------------------------
for t in range(Nt):

    rho = np.sum(f, axis=0)

    ux = np.sum(c[:,0,None,None]*f, axis=0)/rho
    uy = np.sum(c[:,1,None,None]*f, axis=0)/rho

    # ------------------------
    # pseudopotencial
    # ------------------------
    psi = 1 - np.exp(-rho)

    Fx = np.zeros((Nx,Ny))
    Fy = np.zeros((Nx,Ny))

    for i in range(9):
        psi_shift = np.roll(np.roll(psi, c[i,0], axis=0), c[i,1], axis=1)

        Fx += -G*w[i]*psi*psi_shift*c[i,0]
        Fy += -G*w[i]*psi*psi_shift*c[i,1]

    # ------------------------
    # velocidade corrigida
    # ------------------------
    ux = (np.sum(c[:,0,None,None]*f, axis=0) + 0.5*Fx)/rho
    uy = (np.sum(c[:,1,None,None]*f, axis=0) + 0.5*Fy)/rho

    # ------------------------
    # equilíbrio
    # ------------------------
    feq = np.zeros_like(f)

    u2 = ux**2 + uy**2

    for i in range(9):
        cu = 3*(c[i,0]*ux + c[i,1]*uy)

        feq[i] = w[i]*rho*(1 + cu + 0.5*cu**2 - 1.5*u2)

    # ------------------------
    # colisão
    # ------------------------
    f = f - (1/tau)*(f-feq)

    # ------------------------
    # streaming
    # ------------------------
    for i in range(9):
        f[i] = np.roll(f[i], c[i,0], axis=0)
        f[i] = np.roll(f[i], c[i,1], axis=1)

    # ------------------------
    # salvar frames
    # ------------------------
    if t % 5 == 0:
        plt.figure(figsize=(6,6))
        plt.imshow(rho.T, origin='lower')
        plt.colorbar(label='densidade')
        plt.title(f"Shan-Chen t={t}")
        plt.savefig(f"{output_dir}/frame_{t:04d}.png")
        plt.close()

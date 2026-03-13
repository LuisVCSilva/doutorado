import numpy as np
import matplotlib.pyplot as plt
import os

# ----------------------------
# Pasta de saída
# ----------------------------
output_dir = "lbm_d2q9_frames"
os.makedirs(output_dir, exist_ok=True)

# ----------------------------
# Parâmetros
# ----------------------------
Nx, Ny = 100, 100
Nt = 300
tau = 0.8

# ----------------------------
# Velocidades D2Q9
# ----------------------------
c = np.array([
    [0,0],
    [1,0],[0,1],[-1,0],[0,-1],
    [1,1],[-1,1],[-1,-1],[1,-1]
])

# pesos
w = np.array([
    4/9,
    1/9,1/9,1/9,1/9,
    1/36,1/36,1/36,1/36
])

# ----------------------------
# Distribuições
# ----------------------------
f = np.zeros((9, Nx, Ny))

# densidade inicial
rho = np.ones((Nx,Ny))

# perturbação central
rho[Nx//3:2*Nx//3, Ny//3:2*Ny//3] = 2.0

ux = np.zeros((Nx,Ny))
uy = np.zeros((Nx,Ny))

# inicialização equilíbrio
for i in range(9):
    f[i,:,:] = w[i]*rho

# ----------------------------
# Loop temporal
# ----------------------------
for t in range(Nt):

    # ------------------------
    # Grandezas macroscópicas
    # ------------------------
    rho = np.sum(f, axis=0)

    ux = np.sum(c[:,0,None,None]*f, axis=0)/rho
    uy = np.sum(c[:,1,None,None]*f, axis=0)/rho

    # ------------------------
    # Equilíbrio
    # ------------------------
    feq = np.zeros_like(f)

    u2 = ux**2 + uy**2

    for i in range(9):
        cu = 3*(c[i,0]*ux + c[i,1]*uy)

        feq[i,:,:] = w[i]*rho*(
            1 + cu + 0.5*cu**2 - 1.5*u2
        )

    # ------------------------
    # Colisão BGK
    # ------------------------
    f = f - (1/tau)*(f - feq)

    # ------------------------
    # Streaming
    # ------------------------
    for i in range(9):
        f[i,:,:] = np.roll(f[i,:,:], c[i,0], axis=0)
        f[i,:,:] = np.roll(f[i,:,:], c[i,1], axis=1)

    # ------------------------
    # Salvar frame
    # ------------------------
    plt.figure(figsize=(6,6))
    plt.imshow(rho.T, origin='lower')
    plt.colorbar(label='densidade')
    plt.title(f"D2Q9 t={t}")

    plt.savefig(f"{output_dir}/frame_{t:04d}.png")
    plt.close()

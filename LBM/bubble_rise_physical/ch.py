import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

# =============================
# pasta
# =============================
output_dir = "bubble_cahn_hilliard"
os.makedirs(output_dir, exist_ok=True)

# =============================
# domínio
# =============================
Nx, Ny = 180, 300
Nt = 3000

tau = 0.7

# Cahn-Hilliard
M = 0.05
kappa = 0.04

# empuxo
g = 2e-5

# =============================
# D2Q9
# =============================
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

opp = np.array([0,3,4,1,2,7,8,5,6])

# =============================
# campos
# =============================
rho = np.ones((Nx,Ny))
ux = np.zeros((Nx,Ny))
uy = np.zeros((Nx,Ny))

phi = np.ones((Nx,Ny))

# bolha inicial
cx, cy = Nx//2, Ny//6
R = 18

for x in range(Nx):
    for y in range(Ny):
        r = np.sqrt((x-cx)**2 + (y-cy)**2)
        if r < R:
            phi[x,y] = -1.0

# =============================
# distribuição inicial
# =============================
f = np.zeros((9,Nx,Ny))

for i in range(9):
    f[i] = w[i]*rho

# =============================
# laplaciano
# =============================
def laplacian(field):
    return (
        np.roll(field,1,0)+np.roll(field,-1,0)+
        np.roll(field,1,1)+np.roll(field,-1,1)-4*field
    )

# =============================
# loop temporal
# =============================
for t in tqdm(range(Nt), desc="Cahn-Hilliard"):

    # -------------------------
    # LBM macroscópico
    # -------------------------
    rho = np.sum(f, axis=0)

    ux = np.sum(c[:,0,None,None]*f, axis=0)/rho
    uy = np.sum(c[:,1,None,None]*f, axis=0)/rho

    # -------------------------
    # empuxo na bolha
    # -------------------------
    uy += g*(1-phi)

    # -------------------------
    # potencial químico
    # -------------------------
    mu = phi**3 - phi - kappa*laplacian(phi)

    # -------------------------
    # evolução Cahn-Hilliard
    # -------------------------
    dphix = (np.roll(phi,-1,0)-np.roll(phi,1,0))*0.5
    dphiy = (np.roll(phi,-1,1)-np.roll(phi,1,1))*0.5

    adv = ux*dphix + uy*dphiy

    phi += -0.1*adv + M*laplacian(mu)

    phi = np.clip(phi,-1,1)

    # -------------------------
    # equilíbrio LBM
    # -------------------------
    feq = np.zeros_like(f)

    u2 = ux**2 + uy**2

    for i in range(9):
        cu = 3*(c[i,0]*ux + c[i,1]*uy)
        feq[i] = w[i]*rho*(1 + cu + 0.5*cu**2 - 1.5*u2)

    # -------------------------
    # colisão
    # -------------------------
    f = f - (1/tau)*(f-feq)

    # -------------------------
    # streaming
    # -------------------------
    for i in range(9):
        f[i] = np.roll(np.roll(f[i], c[i,0], axis=0), c[i,1], axis=1)

    # -------------------------
    # fundo bounce-back
    # -------------------------
    for i in range(9):
        f[i,:,0] = f[opp[i],:,0]

    # -------------------------
    # topo aberto
    # -------------------------
    for i in range(9):
        f[i,:,-1] = f[i,:,-2]

    # -------------------------
    # salvar
    # -------------------------
    if t % 20 == 0:
        plt.figure(figsize=(6,10))
        plt.imshow(phi.T, origin='lower')
        plt.colorbar(label='phase field')
        plt.title(f"Cahn-Hilliard bubble t={t}")
        plt.savefig(f"{output_dir}/frame_{t:04d}.png")
        plt.close()

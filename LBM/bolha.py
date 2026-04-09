import numpy as np
import matplotlib.pyplot as plt
import os

# =========================
# pasta de saída
# =========================
output_dir = "bubble_rise_physical"
os.makedirs(output_dir, exist_ok=True)

# =========================
# domínio
# =========================
Nx, Ny = 160, 320
Nt = 2500

tau = 0.9
G = -3.8
g = -5e-6

# =========================
# D2Q9
# =========================
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

# =========================
# densidade inicial
# =========================
rho_liq = 1.0
rho_gas = 0.85

rho = rho_liq*np.ones((Nx,Ny))

# bolha no fundo
cx, cy = Nx//2, Ny//6
R = 18

for x in range(Nx):
    for y in range(Ny):
        r = np.sqrt((x-cx)**2 + (y-cy)**2)
        if r < R:
            rho[x,y] = rho_gas

# =========================
# distribuição inicial
# =========================
f = np.zeros((9,Nx,Ny))

for i in range(9):
    f[i] = w[i]*rho

# =========================
# loop temporal
# =========================
for t in range(Nt):

    rho = np.sum(f, axis=0)

    ux = np.sum(c[:,0,None,None]*f, axis=0)/rho
    uy = np.sum(c[:,1,None,None]*f, axis=0)/rho

    # =====================
    # pseudopotencial
    # =====================
    psi = 1 - np.exp(-rho)

    Fx = np.zeros((Nx,Ny))
    Fy = np.zeros((Nx,Ny))

    for i in range(9):
        psi_shift = np.roll(np.roll(psi, c[i,0], axis=0), c[i,1], axis=1)

        Fx += -G*w[i]*psi*psi_shift*c[i,0]
        Fy += -G*w[i]*psi*psi_shift*c[i,1]

    # gravidade suave
    Fy += rho*g

    # =====================
    # velocidade corrigida
    # =====================
    ux = (np.sum(c[:,0,None,None]*f, axis=0)+0.5*Fx)/rho
    uy = (np.sum(c[:,1,None,None]*f, axis=0)+0.5*Fy)/rho

    u2 = ux**2 + uy**2

    # =====================
    # equilíbrio
    # =====================
    feq = np.zeros_like(f)

    for i in range(9):
        cu = 3*(c[i,0]*ux + c[i,1]*uy)
        feq[i] = w[i]*rho*(1 + cu + 0.5*cu**2 - 1.5*u2)

    # =====================
    # colisão
    # =====================
    f = f - (1/tau)*(f-feq)

    # =====================
    # streaming
    # =====================
    f_stream = np.zeros_like(f)

    for i in range(9):
        f_stream[i] = np.roll(np.roll(f[i], c[i,0], axis=0), c[i,1], axis=1)

    f = f_stream

    # =====================
    # paredes inferior bounce-back
    # =====================
    for i in range(9):
        f[i,:,0] = f[opp[i],:,0]

    # =====================
    # topo absorvente
    # =====================
    rho[:, -1] = rho_liq

    for i in range(9):
        f[i,:,-1] = w[i]*rho[:, -1]

    # =====================
    # salvar frames
    # =====================
    if t % 20 == 0:
        plt.figure(figsize=(6,10))
        plt.imshow(rho.T, origin='lower')
        plt.colorbar(label='densidade')
        plt.title(f"Bubble rise t={t}")
        plt.savefig(f"{output_dir}/frame_{t:04d}.png")
        plt.close()

import numpy as np
import matplotlib.pyplot as plt
import os
from tqdm import tqdm

# ======================
# Parâmetros
# ======================
nx, ny = 300, 100
nt = 5000

tau1 = 0.7
tau2 = 0.6

ga = -1.11111111e-7
dp = 1.11111111e-6

G = 1.05
theta = np.pi / 6

Gads1 = -0.05
Gads2 = (Gads1 + G*np.cos(theta)) * 1.03

ux_in = 0.0013
uy_in = 0.0

rho1 = np.full((ny, nx), 0.8)
rho2 = np.zeros((ny, nx))
rho2[:, 0] = 1.0

# D2Q9
ex = np.array([0,1,0,-1,0,1,-1,-1,1])
ey = np.array([0,0,1,0,-1,1,1,-1,-1])
w = np.array([4/9] + [1/9]*4 + [1/36]*4)

# sólidos
solid = np.zeros((ny, nx))
solid[0,:] = 1
solid[-1,:] = 1

# distribuições
f1 = np.zeros((9, ny, nx))
f2 = np.zeros((9, ny, nx))

# ======================
# Funções
# ======================

def psi(rho):
    return 1 - np.exp(-rho)

def feq(rho, ux, uy):
    feq = np.zeros((9, ny, nx))
    u2 = ux**2 + uy**2
    for k in range(9):
        cu = 3*(ex[k]*ux + ey[k]*uy)
        feq[k] = w[k]*rho*(1 + cu + 0.5*cu**2 - 1.5*u2)
    return feq

def macroscopic(f):
    rho = np.sum(f, axis=0)
    px = np.sum(f * ex[:,None,None], axis=0)
    py = np.sum(f * ey[:,None,None], axis=0)

    rho[rho < 1e-12] = 1e-12

    ux = px / rho
    uy = py / rho

    ux[solid == 1] = 0
    uy[solid == 1] = 0

    return rho, ux, uy

def forcing(rho_other, G, Gads):
    Fx = np.zeros_like(rho_other)
    Fy = np.zeros_like(rho_other)

    psi_loc = psi(rho_other)

    for k in range(9):
        shifted = np.roll(np.roll(psi_loc, ey[k], axis=0), ex[k], axis=1)
        Fx += -G * w[k] * shifted * ex[k]
        Fy += -G * w[k] * shifted * ey[k]

    # adsorção (paredes)
    Fx += -Gads * solid
    Fy += -Gads * solid

    return Fx, Fy

def streaming(f):
    f_stream = np.zeros_like(f)
    for k in range(9):
        f_stream[k] = np.roll(f[k], (ey[k], ex[k]), axis=(0,1))
    return f_stream

def collision(f, feq, tau):
    return f - (f - feq)/tau

# ======================
# Inicialização
# ======================
ux1 = np.full((ny, nx), ux_in)
uy1 = np.zeros((ny, nx))

ux2 = np.zeros((ny, nx))
uy2 = np.zeros((ny, nx))

f1 = feq(rho1, ux1, uy1)
f2 = feq(rho2, ux2, uy2)

os.makedirs("frames_multi", exist_ok=True)

# ======================
# Loop principal
# ======================
for t in tqdm(range(nt)):

    # macroscópico
    rho1, ux1, uy1 = macroscopic(f1)
    rho2, ux2, uy2 = macroscopic(f2)

    # velocidade bulk
    rho_tot = rho1 + rho2
    ux_bulk = (rho1*ux1 + rho2*ux2) / rho_tot
    uy_bulk = (rho1*uy1 + rho2*uy2) / rho_tot

    # forças
    Fx1, Fy1 = forcing(rho2, G, Gads1)
    Fx2, Fy2 = forcing(rho1, G, Gads2)

    # velocidades de equilíbrio
    ux1_eq = ux_bulk + tau1*(Fx1 + ga)
    uy1_eq = uy_bulk + tau1*(Fy1)

    ux2_eq = ux_bulk + tau2*(Fx2 + ga)
    uy2_eq = uy_bulk + tau2*(Fy2)

    # gradiente de pressão
    ux1_eq += tau1 * dp
    ux2_eq += tau2 * dp

    # equilíbrio
    feq1 = feq(rho1, ux1_eq, uy1_eq)
    feq2 = feq(rho2, ux2_eq, uy2_eq)

    # colisão
    f1 = collision(f1, feq1, tau1)
    f2 = collision(f2, feq2, tau2)

    # propagação
    f1 = streaming(f1)
    f2 = streaming(f2)

    # saída
    f1[:, :, -1] = f1[:, :, -2]
    f2[:, :, -1] = f2[:, :, -2]

    # ======================
    # OUTPUT
    # ======================
    if t % 50 == 0:

        plt.figure(figsize=(6,3))
        plt.imshow(rho2, origin='lower')
        plt.colorbar(label='densidade água')
        plt.title(f"Step {t}")
        plt.tight_layout()
        plt.savefig(f"frames_multi/frame_{t:05d}.png")
        plt.close()

        print(f"Step {t} | max rho2 = {rho2.max():.4f}")
